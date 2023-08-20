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
    auto & particles = mModel->GetParticles();

    LogMessage("----------------");
    LogMessage("----------------");

    for (auto const & p : particles)
    {
        LogMessage("----------------");

        auto & particleState = particles.GetState(p);

        if (!particleState.TargetPosition.has_value())
        {
            //
            // We need a trajectory
            //

            LogMessage("Particle ", p, ": trajectory calculation");

            // Calculate target position
            vec2f targetPosition;
            if (mCurrentParticleTrajectory.has_value()
                && mCurrentParticleTrajectory->ParticleIndex == p)
            {
                // We have a user-imposed trajectory for this particle...
                // ...use the provided trajectory

                targetPosition = mCurrentParticleTrajectory->TargetPosition;
            }
            else
            {
                // Use physics to calculate trajectory

                targetPosition = 
                    particles.GetPosition(p) // Physics moves from current pos, before mesh move
                    + CalculatePhysicsDeltaPos(p, labParameters); 
            }

            // Calculate source position
            vec2f sourcePosition;
            if (particleState.ConstrainedState.has_value())
            {
                // Constrained state: trajectory is from current bary coords in cur triangle 
                // (i.e. in mesh's reference frame) up to target position calculated by physics

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

            //
            // Update physics for trajectory
            //

            particles.SetPosition(p, sourcePosition);
            particleState.TargetPosition = targetPosition;
            particles.SetVelocity(p, (targetPosition - sourcePosition) / LabParameters::SimulationTimeStepDuration);

            //
            // Update trajectory notification
            //

            mCurrentParticleTrajectory.emplace(p, targetPosition);
            mCurrentParticleTrajectoryNotification.reset();
        }
        else
        {
            LogMessage("Particle ", p, ": state update");

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

vec2f LabController::CalculatePhysicsDeltaPos(
    ElementIndex particleIndex,
    LabParameters const & labParameters) const
{
    float const dt = LabParameters::SimulationTimeStepDuration;
    float const particleMass = LabParameters::ParticleMass * labParameters.MassAdjustment;

    Particles const & particles = mModel->GetParticles();

    Vertices const & vertices = mModel->GetMesh().GetVertices();
    Edges const & edges = mModel->GetMesh().GetEdges();
    Triangles const & triangles = mModel->GetMesh().GetTriangles();

    vec2f forces = particles.GetWorldForce(particleIndex) * labParameters.GravityAdjustment;

    //
    // Check whether we're on a floor
    //

    auto const & particleState = particles.GetState(particleIndex);
    if (particleState.ConstrainedState.has_value())
    {
        ElementIndex const currentTriangleElementIndex = particleState.ConstrainedState->CurrentTriangle;

        std::int32_t currentEdgeOrdinal = -1;
        if (particleState.ConstrainedState->CurrentTriangleBarycentricCoords.x == 0.0f)
        {
            currentEdgeOrdinal = 1;
        }
        else if (particleState.ConstrainedState->CurrentTriangleBarycentricCoords.y == 0.0f)
        {
            currentEdgeOrdinal = 2;
        }
        else if (particleState.ConstrainedState->CurrentTriangleBarycentricCoords.z == 0.0f)
        {
            currentEdgeOrdinal = 0;
        }

        if (currentEdgeOrdinal >= 0)
        {
            ElementIndex const currentEdgeElementIndex = triangles.GetSubEdges(currentTriangleElementIndex).EdgeIndices[currentEdgeOrdinal];

            // Calculate ghost condition - we don't consider a surface floor (and thus allow 
            // particle to ghost through it) if it's in a triangle made fully of floors
            bool const isGhost =
                edges.GetSurfaceType(triangles.GetSubEdgeAIndex(currentTriangleElementIndex)) == SurfaceType::Floor
                && edges.GetSurfaceType(triangles.GetSubEdgeBIndex(currentTriangleElementIndex)) == SurfaceType::Floor
                && edges.GetSurfaceType(triangles.GetSubEdgeCIndex(currentTriangleElementIndex)) == SurfaceType::Floor;

            if (edges.GetSurfaceType(currentEdgeElementIndex) == SurfaceType::Floor
                && !isGhost)
            {
                //
                // On floor edge
                //                

                //
                // Add friction
                //

                // TODO

                //
                // Integrate, including eventual bounce velocity from a previous impact
                //

                vec2f const deltaPos =
                    particles.GetVelocity(particleIndex) * dt
                    + forces / particleMass * dt * dt;

                //
                // If we're moving against the floor, flatten trajectory
                //

                vec2f const edgeVector = mModel->GetMesh().GetTriangles().GetEdgeVector(
                    currentTriangleElementIndex,
                    currentEdgeOrdinal,
                    vertices);

                if (deltaPos.dot(edgeVector.to_perpendicular()) > 0.0f) // Normal to edge is directed outside of triangle (i.e. towards floor)
                {
                    LogMessage("  Particle is on floor edge, moving against it");

                    //
                    // Flatten trajectory - i.e. take component of deltapos along floor
                    //

                    vec2f const edgeDir = edgeVector.normalise();

                    return edgeDir * deltaPos.dot(edgeDir);
                }
                else
                {
                    LogMessage("  Particle is on floor edge, but not against it");

                    return deltaPos;
                }
            }
        }
    }    

    //
    // Not on a floor edge, use pure force
    //

    LogMessage("  Particle not on floor edge; using pure force");

    //
    // Integrate, including eventual bounce velocity from a previous impact
    //

    return
        particles.GetVelocity(particleIndex) * dt
        + forces / particleMass * dt * dt;

}

bool LabController::UpdateParticleState(
    ElementIndex particleIndex,
    LabParameters const & labParameters)
{
    Particles & particles = mModel->GetParticles();
    Vertices const & vertices = mModel->GetMesh().GetVertices();
    Edges const & edges = mModel->GetMesh().GetEdges();
    Triangles const & triangles = mModel->GetMesh().GetTriangles();

    auto & state = particles.GetState(particleIndex);
    assert(state.TargetPosition.has_value());
    vec2f targetPosition = *state.TargetPosition;
    vec2f trajectory = targetPosition - particles.GetPosition(particleIndex);

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

    vec3f targetBarycentricCoords = triangles.ToBarycentricCoordinates(
        targetPosition,
        currentTriangle,
        vertices);

    //
    // If target is strictly in triangle, we move to target
    //

    if (targetBarycentricCoords.x >= 0.0f && targetBarycentricCoords.x <= 1.0f
        && targetBarycentricCoords.y >= 0.0f && targetBarycentricCoords.y <= 1.0f
        && targetBarycentricCoords.z >= 0.0f && targetBarycentricCoords.z <= 1.0f)
    {
        LogMessage("  Target is on/in triangle, moving to target");

        // ...move to target position directly
        particles.SetPosition(particleIndex, targetPosition);
        state.ConstrainedState->CurrentTriangleBarycentricCoords = targetBarycentricCoords;

        return false; //  // We'll end at next iteration
    }

    // TODOHERE
    (void)labParameters;
    return true;

    /*
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

                //
                // Update velocity with normal response
                //

                vec2f particleVelocity = particles.GetVelocity(particleIndex);

                // Calculate edge direction (from point of view of inside triangle,
                // hence CW)
                vec2f const edgeDir = (
                    vertices.GetPosition(triangles.GetVertexIndices(currentTriangle)[(currentEdge + 1) % 3])
                    - vertices.GetPosition(triangles.GetVertexIndices(currentTriangle)[currentEdge])
                    ).normalise();

                // Calculate edge normal (pointing outside the triangle, into the floor)
                vec2f const edgeNormal = edgeDir.to_perpendicular();

                // Calculate the component of the particle's velocity along the normal,
                // i.e. towards the interior of the floor...
                float const particleVelocityAlongNormal = particleVelocity.dot(edgeNormal);

                // ...if positive, we have an impact
                if (particleVelocityAlongNormal > 0.0f)
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

                //
                // Flatten trajectory
                //

                trajectory = edgeDir * trajectory.dot(edgeDir);

                LogMessage("TODOTEST:  targetPosition: old=", targetPosition, " new=", particles.GetPosition(particleIndex) + trajectory);

                // Update target pos and target barycentric coords
                targetPosition = particles.GetPosition(particleIndex) + trajectory;
                state.TargetPosition = targetPosition;
                targetBarycentricCoords = triangles.ToBarycentricCoordinates(
                    targetPosition,
                    currentTriangle,
                    vertices);

                LogMessage("TODOTEST:  new targetBarycentricCoords=", targetBarycentricCoords);

                // Update velocity
                // TODOHERE

                return false;
            }
            else
            {
                //
                // Not floor, climb over edge
                //

                LogMessage("  Climbing over non-floor edge");

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
    */
}