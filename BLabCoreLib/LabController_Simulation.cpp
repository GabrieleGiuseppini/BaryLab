/***************************************************************************************
 * Original Author:		Gabriele Giuseppini
 * Created:				2023-07-23
 * Copyright:			Gabriele Giuseppini  (https://github.com/GabrieleGiuseppini)
 ***************************************************************************************/
#include "LabController.h"

#include "BLabMath.h"
#include "Log.h"

#include <limits>

void LabController::InitializeParticleRegime(ElementIndex particleIndex)
{
    std::optional<Particles::StateType::ConstrainedStateType> constrainedState;

    ElementIndex const triangleIndex = FindTriangleContaining(mModel->GetParticles().GetPosition(particleIndex));
    if (triangleIndex != NoneElementIndex)
    {
        vec3f const barycentricCoords = mModel->GetMesh().GetTriangles().ToBarycentricCoordinates(
            mModel->GetParticles().GetPosition(particleIndex),
            triangleIndex,
            mModel->GetMesh().GetVertices());

        constrainedState.emplace(
            triangleIndex,
            barycentricCoords);
    }

    mModel->GetParticles().GetState(particleIndex).ConstrainedState = constrainedState;
}

void LabController::UpdateSimulation(LabParameters const & labParameters)
{
    float const dt = LabParameters::SimulationTimeStepDuration;

    auto & particles = mModel->GetParticles();

    float const particleMass = LabParameters::ParticleMass * labParameters.MassAdjustment;

    for (auto const & p : particles)
    {
        auto & particleState = particles.GetState(p);

        if (!particleState.TargetPosition.has_value())
        {
            //
            // We need a trajectory
            //
            
            vec2f targetPosition;

            if (mCurrentParticleTrajectory.has_value()
                && mCurrentParticleTrajectory->ParticleIndex == p)
            {
                // Use the provided trajectory

                targetPosition = mCurrentParticleTrajectory->TargetPosition;
            }
            else
            {
                // Use physics to calculate trajectory, including momentum

                vec2f const forces = particles.GetWorldForce(p) * labParameters.GravityAdjustment;

                vec2f const deltaPos =
                    particles.GetVelocity(p) * dt
                    + forces / LabParameters::ParticleMass * dt * dt;

                targetPosition = particles.GetPosition(p) + deltaPos;
            }

            // Calculate source position
            vec2f sourcePosition;
            if (particleState.ConstrainedState.has_value())
            {
                // Constrained state: trajectory is from current bary coords in cur triangle up to target position calculated by physics
                sourcePosition = mModel->GetMesh().GetTriangles().FromBarycentricCoordinates(
                    particleState.ConstrainedState->CurrentTriangleBarycentricCoords,
                    particleState.ConstrainedState->CurrentTriangle,
                    mModel->GetMesh().GetVertices());
            }
            else
            {
                // Free state: from current (absolute) position up to target position calculated by physics
                sourcePosition = particles.GetPosition(p);
            }

            LogMessage("TODO: s=", sourcePosition, " t=", targetPosition);

            // Update velocity for trajectory
            particles.SetVelocity(p, (targetPosition - sourcePosition) / dt);

            // Transition state

            particleState.TargetPosition = targetPosition;

            mCurrentParticleTrajectory.emplace(p, targetPosition);
            mCurrentParticleTrajectoryNotification.reset();
        }
        else
        {
            assert(particleState.TargetPosition.has_value());

            bool hasCompleted = UpdateParticleState(p, labParameters);
            if (hasCompleted)
            {
                LogMessage("Particle ", p, " COMPLETED");

                // Destroy state
                particleState.TargetPosition.reset();
                mCurrentParticleTrajectory.reset();
            }
        }
    }
}

bool LabController::UpdateParticleState(
    ElementIndex particleIndex,
    LabParameters const & labParameters)
{
    LogMessage("--------------------------------------");
    LogMessage("P ", particleIndex);

    Particles & particles = mModel->GetParticles();
    Vertices const & vertices = mModel->GetMesh().GetVertices();
    Edges const & edges = mModel->GetMesh().GetEdges();
    Triangles const & triangles = mModel->GetMesh().GetTriangles();

    auto & state = particles.GetState(particleIndex);
    assert(state.TargetPosition.has_value());
    vec2f const targetPosition = *state.TargetPosition;

    //
    // If we are at target, we're done
    //

    if (particles.GetPosition(particleIndex) == targetPosition)
    {
        // Reached destination

        LogMessage("  Reached destination");

        return true;
    }

    //
    // If we're free, we move to target
    //

    if (!state.ConstrainedState)
    {
        LogMessage("  Particle is free, moving to target");

        // ...move to target position directly
        particles.SetPosition(particleIndex, targetPosition);

        return false; //  // We'll end at next iteration
    }

    //
    // We are in constrained state
    //

    LogMessage("  Particle is in constrained state");

    assert(state.ConstrainedState.has_value());    

    ElementIndex const currentTriangle = state.ConstrainedState->CurrentTriangle;

    vec3f const targetBarycentricCoords = triangles.ToBarycentricCoordinates(
        targetPosition,
        currentTriangle,
        vertices);

    //
    // If target is strictly in triangle, we move to target
    //

    if (targetBarycentricCoords.x > 0.0f && targetBarycentricCoords.x < 1.0f
        && targetBarycentricCoords.y > 0.0f && targetBarycentricCoords.y < 1.0f
        && targetBarycentricCoords.z > 0.0f && targetBarycentricCoords.z < 1.0f)
    {
        LogMessage("  Target is in triangle, moving to target");

        // ...move to target position directly
        particles.SetPosition(particleIndex, targetPosition);
        state.ConstrainedState->CurrentTriangleBarycentricCoords = targetBarycentricCoords;

        return false; //  // We'll end at next iteration
    }

    //
    // Check whether we are on an edge
    //

    // Note: if we are exactly at a vertex, we pick an arbitrary edge here
    ElementIndex currentEdge = NoneElementIndex;
    if (state.ConstrainedState->CurrentTriangleBarycentricCoords.x == 0.0f)
    {
        currentEdge = 1;
    }
    else if (state.ConstrainedState->CurrentTriangleBarycentricCoords.y == 0.0f)
    {
        currentEdge = 2;
    }
    else if (state.ConstrainedState->CurrentTriangleBarycentricCoords.z == 0.0f)
    {
        currentEdge = 0;
    }

    if (currentEdge != NoneElementIndex)
    {
        //
        // We are on an edge
        //

        LogMessage("  Particle is on edge ", currentEdge);

        //
        // Check trajectory direction wrt normal to this edge
        //

        // TODOTEST
        {
            vec2f const trajectory = targetPosition - particles.GetPosition(particleIndex);

            // Calculate pseudonormal to edge, considering that we are *inside* the triangle
            // (points outside)
            vec2f const edgePNormal = (
                vertices.GetPosition(triangles.GetVertexIndices(currentTriangle)[(currentEdge + 1) % 3])
                - vertices.GetPosition(triangles.GetVertexIndices(currentTriangle)[currentEdge])
                ).to_perpendicular();

            assert((edgePNormal.dot(trajectory) > 0.0f && targetBarycentricCoords[(currentEdge + 2) % 3] < 0.0f)
                || (edgePNormal.dot(trajectory) < 0.0f && targetBarycentricCoords[(currentEdge + 2) % 3] > 0.0f)
                || (edgePNormal.dot(trajectory) == 0.0f && targetBarycentricCoords[(currentEdge + 2) % 3] == 0.0f));
        }

        if (targetBarycentricCoords[(currentEdge + 2) % 3] < 0.0f)
        {
            //
            // We are on an edge, wanting to go strictly outside
            //

            LogMessage("  Trajectory is strictly outward");

            ElementIndex const currentEdgeElement = triangles.GetSubEdges(currentTriangle).EdgeIndices[currentEdge];

            // Calculate ghost condition - we don't consider a surface floor (and thus allow 
            // particle to ghost through it) if it's in a triangle made fully of floors
            bool const isGhost =
                edges.GetSurfaceType(triangles.GetSubEdgeAIndex(currentTriangle)) == SurfaceType::Floor
                && edges.GetSurfaceType(triangles.GetSubEdgeBIndex(currentTriangle)) == SurfaceType::Floor
                && edges.GetSurfaceType(triangles.GetSubEdgeCIndex(currentTriangle)) == SurfaceType::Floor;

            if (edges.GetSurfaceType(currentEdgeElement) == SurfaceType::Floor
                && !isGhost)
            {
                //
                // Impact
                //

                LogMessage("  Impact");

                vec2f const particleVelocity = particles.GetVelocity(particleIndex);

                // Calculate edge normal (positive pointing into the floor)
                vec2f const edgeNormal = (
                    vertices.GetPosition(triangles.GetVertexIndices(currentTriangle)[currentEdge])
                    - vertices.GetPosition(triangles.GetVertexIndices(currentTriangle)[(currentEdge + 1) % 3])
                    ).to_perpendicular()
                    .normalise();

                // Calculate the component of the particle's velocity along the normal,
                // i.e. towards the interior of the floor...
                float const particleVelocityAlongNormal = particleVelocity.dot(edgeNormal);

                // ...if negative, we have an impact
                if (particleVelocityAlongNormal < 0.0f)
                {
                    // Decompose particle velocity into normal and tangential
                    vec2f const normalVelocity = edgeNormal * particleVelocityAlongNormal;
                    vec2f const tangentialVelocity = particleVelocity - normalVelocity;

                    // Calculate normal reponse: Vn' = -e*Vn (e = elasticity, [0.0 - 1.0])
                    vec2f const normalResponse =
                        -normalVelocity
                        * labParameters.Elasticity;

                    // Calculate tangential response: Vt' = a*Vt (a = (1.0-friction), [0.0 - 1.0])
                    vec2f const tangentialResponse =
                        tangentialVelocity
                        * (1.0f - labParameters.Friction);

                    LogMessage("TODO: n=", normalResponse, " t=", tangentialResponse);

                    // Set velocity to resultant collision velocity
                    particles.SetVelocity(
                        particleIndex,
                        normalResponse + tangentialResponse);
                }
                else
                {
                    // Maintain current velocity
                }

                return true;
            }
            else
            {
                //
                // Climb over edge
                //

                LogMessage("  Climbing over edge");

                // Find opposite triangle
                ElementIndex const oppositeTriangle = edges.GetOppositeTriangle(currentEdgeElement, currentTriangle);
                if (oppositeTriangle == NoneElementIndex)
                {
                    //
                    // Become free
                    //

                    LogMessage("  No opposite triangle found, becoming free");

                    state.ConstrainedState.reset();

                    return false;
                }
                else
                {
                    //
                    // Move to edge of opposite triangle 
                    //

                    int const oppositeTriangleEdgeIndex = triangles.GetSubEdgeIndex(oppositeTriangle, currentEdgeElement);

                    LogMessage("  Moving to edge ", oppositeTriangleEdgeIndex, " of opposite triangle ", oppositeTriangle);

                    state.ConstrainedState->CurrentTriangle = oppositeTriangle;

                    // Calculate new barycentric coords (wrt opposite triangle)
                    vec3f newBarycentricCoords; // In new triangle
                    newBarycentricCoords[(oppositeTriangleEdgeIndex + 2) % 3] = 0.0f;
                    newBarycentricCoords[oppositeTriangleEdgeIndex] = state.ConstrainedState->CurrentTriangleBarycentricCoords[(currentEdge + 1) % 3];
                    newBarycentricCoords[(oppositeTriangleEdgeIndex + 1) % 3] = state.ConstrainedState->CurrentTriangleBarycentricCoords[currentEdge];                    

                    LogMessage("  B-Coords: ", state.ConstrainedState->CurrentTriangleBarycentricCoords, " -> ", newBarycentricCoords);

                    state.ConstrainedState->CurrentTriangleBarycentricCoords = newBarycentricCoords;

                    return false;
                }
            }
        }
        else
        {
            LogMessage("  Trajectory is parallel or towards interior");
        }
    }

    //
    // We're inside triangle, and target is outside triangle,
    // OR we're on edge, with trajectory pointing inside or parallel, and target is outside triangle
    //

    //
    // Find closest intersection with one of the edges, excluding edge we are on
    //
    // Guaranteed to exist and within trajectory: target outside of triangle and 
    // we're inside or on an edge pointing inside
    //

    ElementIndex intersectionVertex = NoneElementIndex;
    float minIntersectionT = std::numeric_limits<float>::max();

    for (ElementIndex vi = 0; vi < 3; ++vi)
    {
        // Skip current edge (if any)
        if ((vi + 1) % 3 != currentEdge)
        {
            float const den = state.ConstrainedState->CurrentTriangleBarycentricCoords[vi] - targetBarycentricCoords[vi];
            float const t = IsAlmostZero(den)
                ? std::numeric_limits<float>::max() // Parallel, meets at infinity
                : state.ConstrainedState->CurrentTriangleBarycentricCoords[vi] / den;

            LogMessage("  t[v", vi, "]=", t);

            // Cull backward intersections
            if (t >= 0.0f)
            {
                if (t < minIntersectionT)
                {
                    intersectionVertex = vi;
                    minIntersectionT = t;
                }
            }
        }
    }

    assert(intersectionVertex != NoneElementIndex); // Guaranteed to exist
    assert(minIntersectionT >= 0.0f && minIntersectionT <= 1.0f); // Guaranteed to exist, within trajectory

    //
    // Move to intersection
    //

    ElementIndex const intersectionEdge = (intersectionVertex + 1) % 3;

    LogMessage("  Moving to intersection with edge ", intersectionEdge);

    // Calculate intersection barycentric coordinates
    vec3f intersectionBarycentricCoords;
    intersectionBarycentricCoords[intersectionVertex] = 0.0f;
    float const lNext = // Barycentric coord of next vertex at intersection
        state.ConstrainedState->CurrentTriangleBarycentricCoords[(intersectionVertex + 1) % 3] * (1.0f - minIntersectionT)
        + targetBarycentricCoords[(intersectionVertex + 1) % 3] * minIntersectionT;
    intersectionBarycentricCoords[(intersectionVertex + 1) % 3] = lNext;
    intersectionBarycentricCoords[(intersectionVertex + 2) % 3] = 1.0f - lNext;

    // Move to intersection

    vec2f const intersectionPosition = triangles.FromBarycentricCoordinates(
        intersectionBarycentricCoords,
        currentTriangle,
        vertices);

    particles.SetPosition(particleIndex, intersectionPosition);
    state.ConstrainedState->CurrentTriangleBarycentricCoords = intersectionBarycentricCoords;
    
    return false;
}