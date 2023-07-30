/***************************************************************************************
 * Original Author:		Gabriele Giuseppini
 * Created:				2023-07-23
 * Copyright:			Gabriele Giuseppini  (https://github.com/GabrieleGiuseppini)
 ***************************************************************************************/
#include "LabController.h"

#include "Log.h"

#include <limits>

void LabController::InitializeParticleState(ElementIndex particleIndex)
{
    if (mCurrentParticleTrajectory.has_value()
        && mCurrentParticleTrajectory->ParticleIndex == particleIndex)
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

        mModel->GetParticles().GetState(particleIndex).emplace(
            constrainedState,
            mCurrentParticleTrajectory->TargetPosition);
    }
    else
    {
        mModel->GetParticles().GetState(particleIndex).reset();
    }
}

void LabController::UpdateSimulation(LabParameters const & /*labParameters*/)
{
    //////
    ////// Particle physics
    //////

    ////float const dt = LabParameters::SimulationTimeStepDuration;

    ////auto & particles = mModel->GetParticles();

    ////float const particleMass = LabParameters::ParticleMass * labParameters.MassAdjustment;

    ////for (auto const & p : particles)
    ////{
    ////    vec2f const forces = particles.GetWorldForce(p) * labParameters.GravityAdjustment;

    ////    vec2f const deltaPos =
    ////        particles.GetVelocity(p) * dt
    ////        + forces / LabParameters::ParticleMass * dt * dt;

    ////    particles.SetPosition(p, particles.GetPosition(p) + deltaPos);
    ////    particles.SetVelocity(p, deltaPos / dt);
    ////}

    auto & particles = mModel->GetParticles();

    for (auto const & p : particles)
    {
        if (particles.GetState(p).has_value())
        {
            bool hasCompleted = UpdateParticleState(p);
            if (hasCompleted)
            {
                LogMessage("Particle ", p, " COMPLETED");

                // Destroy state
                particles.GetState(p).reset();
                mCurrentParticleTrajectory.reset();
            }
        }
    }
}

bool LabController::UpdateParticleState(ElementIndex particleIndex)
{
    LogMessage("--------------------------------------");
    LogMessage("P ", particleIndex);

    Particles & particles = mModel->GetParticles();
    Vertices const & vertices = mModel->GetMesh().GetVertices();
    Edges const & edges = mModel->GetMesh().GetEdges();
    Triangles const & triangles = mModel->GetMesh().GetTriangles();

    auto & state = particles.GetState(particleIndex);
    assert(state.has_value());

    //
    // If we are at target, we're done
    //

    if (particles.GetPosition(particleIndex) == state->TargetPosition)
    {
        // Reached destination

        LogMessage("  Reached destination");

        return true;
    }

    //
    // If we're free, we move to target
    //

    if (!state->ConstrainedState)
    {
        LogMessage("  Particle is free, moving to target");

        // ...move to target position directly
        particles.SetPosition(particleIndex, state->TargetPosition);

        return false; //  // We'll end at next iteration
    }

    //
    // We are in constrained state
    //

    LogMessage("  Particle is in constrained state");

    assert(state->ConstrainedState.has_value());    

    ElementIndex const currentTriangle = state->ConstrainedState->CurrentTriangle;

    vec3f const targetBarycentricCoords = triangles.ToBarycentricCoordinates(
        state->TargetPosition,
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
        particles.SetPosition(particleIndex, state->TargetPosition);
        state->ConstrainedState->CurrentTriangleBarycentricCoords = targetBarycentricCoords;

        return false; //  // We'll end at next iteration
    }

    //
    // Check whether we are on an edge
    //

    // Note: if we are exactly at a vertex, we pick an arbitrary edge here
    ElementIndex currentEdge = NoneElementIndex;
    if (state->ConstrainedState->CurrentTriangleBarycentricCoords.x == 0.0f)
    {
        currentEdge = 1;
    }
    else if (state->ConstrainedState->CurrentTriangleBarycentricCoords.y == 0.0f)
    {
        currentEdge = 2;
    }
    else if (state->ConstrainedState->CurrentTriangleBarycentricCoords.z == 0.0f)
    {
        currentEdge = 0;
    }

    if (currentEdge != NoneElementIndex)
    {
        LogMessage("  Particle is on edge ", currentEdge);

        //
        // Check trajectory direction wrt normal to this edge
        // TODO: see if can be replaced by checking signs of target b-coords
        //

        vec2f const trajectory = state->TargetPosition - particles.GetPosition(particleIndex);

        // Calculate pseudonormal to edge, considering that we are *inside* the triangle
        // (points outside)
        vec2f const edgePNormal = (
            vertices.GetPosition(triangles.GetVertexIndices(currentTriangle)[(currentEdge + 1) % 3])
            - vertices.GetPosition(triangles.GetVertexIndices(currentTriangle)[currentEdge])
            ).to_perpendicular();

        if (edgePNormal.dot(trajectory) > 0.0f) // Trajectory leaves this edge outbound
        {
            //
            // We are on an edge, wanting to go strictly outside
            //

            LogMessage("  Trajectory is sitrctly outward");

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

                    state->ConstrainedState.reset();

                    return false;
                }
                else
                {
                    //
                    // Move to edge of opposite triangle 
                    //

                    int const oppositeTriangleEdgeIndex = triangles.GetSubEdgeIndex(oppositeTriangle, currentEdgeElement);

                    LogMessage("  Moving to edge ", oppositeTriangleEdgeIndex, " of opposite triangle ", oppositeTriangle);

                    state->ConstrainedState->CurrentTriangle = oppositeTriangle;

                    // Calculate new barycentric coords (wrt opposite triangle)
                    vec3f newBarycentricCoords; // In new triangle
                    newBarycentricCoords[(oppositeTriangleEdgeIndex + 2) % 3] = 0.0f;
                    newBarycentricCoords[oppositeTriangleEdgeIndex] = state->ConstrainedState->CurrentTriangleBarycentricCoords[(currentEdge + 1) % 3];
                    newBarycentricCoords[(oppositeTriangleEdgeIndex + 1) % 3] = state->ConstrainedState->CurrentTriangleBarycentricCoords[currentEdge];                    

                    LogMessage("  B-Coords: ", state->ConstrainedState->CurrentTriangleBarycentricCoords, " -> ", newBarycentricCoords);

                    state->ConstrainedState->CurrentTriangleBarycentricCoords = newBarycentricCoords;

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
            float const den = state->ConstrainedState->CurrentTriangleBarycentricCoords[vi] - targetBarycentricCoords[vi];
            float const t = (den == 0.0f) // TODO: with epsilon (GameMath)
                ? std::numeric_limits<float>::max() // Parallel, meets at infinity
                : state->ConstrainedState->CurrentTriangleBarycentricCoords[vi] / den;

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
        state->ConstrainedState->CurrentTriangleBarycentricCoords[(intersectionVertex + 1) % 3] * (1.0f - minIntersectionT)
        + targetBarycentricCoords[(intersectionVertex + 1) % 3] * minIntersectionT;
    intersectionBarycentricCoords[(intersectionVertex + 1) % 3] = lNext;
    intersectionBarycentricCoords[(intersectionVertex + 2) % 3] = 1.0f - lNext;

    // Move to intersection

    vec2f const newPosition = triangles.FromBarycentricCoordinates(
        intersectionBarycentricCoords,
        currentTriangle,
        vertices);

    particles.SetPosition(particleIndex, newPosition);

    state->ConstrainedState->CurrentTriangleBarycentricCoords = intersectionBarycentricCoords;
    
    return false;
}