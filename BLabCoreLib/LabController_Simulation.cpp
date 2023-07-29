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

    if (!state->ConstrainedState)
    {
        // Free regime...

        // ...move to target position directly
        particles.SetPosition(particleIndex, state->TargetPosition);
    }

    if (particles.GetPosition(particleIndex) == state->TargetPosition)
    {
        // Reached destination

        LogMessage("  Reached destination");

        return true;
    }

    //
    // We must be in constrained regime now
    //

    assert(state->ConstrainedState.has_value());

    ElementIndex const currentTriangle = state->ConstrainedState->CurrentTriangle;

    // Calculate ghost condition - we don't consider a surface floor (and thus allow 
    // particle to ghost through it) if it's in a triangle made fully of floors
    bool const isGhost =
        edges.GetSurfaceType(triangles.GetSubEdgeAIndex(currentTriangle)) == SurfaceType::Floor
        && edges.GetSurfaceType(triangles.GetSubEdgeBIndex(currentTriangle)) == SurfaceType::Floor
        && edges.GetSurfaceType(triangles.GetSubEdgeCIndex(currentTriangle)) == SurfaceType::Floor;

    // Calculate trajectory == Target - CurrentPost
    vec2f const trajectory = state->TargetPosition - particles.GetPosition(particleIndex);

    //
    // If we are on a floor edge and we're moving against it, we're done;
    // Note: we ensure later that we move forcibly to a zero bary-coord when we hit an edge
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

    if (currentEdge != NoneElementIndex
        && edges.GetSurfaceType(triangles.GetSubEdges(currentTriangle).EdgeIndices[currentEdge]) == SurfaceType::Floor
        && !isGhost)
    {
        // Calculate pseudonormal, considering that we are *inside* the triangle
        // (points outside)
        vec2f const edgePNormal = (
            vertices.GetPosition(triangles.GetVertexIndices(currentTriangle)[(currentEdge + 1) % 3])
            - vertices.GetPosition(triangles.GetVertexIndices(currentTriangle)[currentEdge])
            ).to_perpendicular();

        if (edgePNormal.dot(trajectory) >= 0.0f) // Trajectory leaves this edge
        {
            //
            // We're going against the floor - floor impact
            //
            // (actually we've impacted at the previous step, here we just realize it)
            //

            LogMessage("  Impact on edge ", currentEdge);

            return true;
        }            
    }

    //
    // Calculate nearest intersection with triangle's edges
    //
    // The trajectory can be parameterized as (1 − t)*P + t*T, with P and T being the
    // particle and target points, in either coordinate system. By expressing P and T
    // in barycentric coordinates, the touch/cross point of the line with each edge
    // is found by imposing the parameterized trajectory component to be zero, 
    // yielding ti = lpi/(lpi - lti)
    //

    // Calculate barycentric coordinates of target position wrt current triangle
    vec3f const targetBarycentricCoords = triangles.ToBarycentricCoordinates(
        state->TargetPosition,
        currentTriangle,
        vertices);

    LogMessage("  Target pos wrt current triangle ", state->ConstrainedState->CurrentTriangle, ": ", targetBarycentricCoords);

    ElementIndex intersectionVertex = NoneElementIndex;
    float minIntersectionT = std::numeric_limits<float>::max();

    for (ElementIndex vi = 0; vi < 3; ++vi)
    {
        // Skip current edge
        if ((vi + 1) % 3 != currentEdge)
        {
            // TODO: do float_arr
            float const den = state->ConstrainedState->CurrentTriangleBarycentricCoords[vi] - targetBarycentricCoords[vi];
            float const t = (den == 0.0f) // TODO: with epsilon (GameMath)
                ? std::numeric_limits<float>::max() // Parallel, meets at infinity
                : state->ConstrainedState->CurrentTriangleBarycentricCoords[vi] / den;

            LogMessage("    t[v", vi, "]=", t);

            if (t > 0.0f) // If 0.0 we're on the edge itself, want to skip that trivial intersection
            {
                // Meets ahead - in the direction of trajectory
                if (t < minIntersectionT)
                {
                    intersectionVertex = vi;
                    minIntersectionT = t;
                }
            }
        }
    }

    // TODO: see if we find by mistake intersection with verlenging of edge when we are on an edge and traj points outside
    ElementIndex newEdgeIndex;
    vec2f newPosition;
    if (intersectionVertex != NoneElementIndex)
    {
        //
        // We intersect an edge
        //

        assert(minIntersectionT >= 0.0f); // Meets ahead - in the direction of trajectory

        if (minIntersectionT > 1.0f)
        {
            // No intersection before end of trajectory =>
            // trajectory does not touch any edge before end of trajectory =>
            // end-of-trajectory is internal
            assert(
                targetBarycentricCoords.x > 0.0f && targetBarycentricCoords.x < 1.0f
                && targetBarycentricCoords.y > 0.0f && targetBarycentricCoords.y < 1.0f
                && targetBarycentricCoords.z > 0.0f && targetBarycentricCoords.z < 1.0f);

            LogMessage("  Intersection is after end of trajectory");

            // ...move to target and we're done
            particles.SetPosition(particleIndex, state->TargetPosition);
            return true;
        }

        //
        // Trajectory intersects an edge before end-of-trajectory
        //

        ElementIndex const intersectionEdge = (intersectionVertex + 1) % 3;
        LogMessage("  Intersection on edge ", intersectionEdge, " @ t=", minIntersectionT);

        // Calculate intersection's barycentric coordinates
        vec3f intersectionBarycentricCoords;
        intersectionBarycentricCoords[intersectionVertex] = 0.0f;
        float const lNext = // Barycentric coord of next vertex at intersection
            state->ConstrainedState->CurrentTriangleBarycentricCoords[(intersectionVertex + 1) % 3] * (1.0f - minIntersectionT)
            + targetBarycentricCoords[(intersectionVertex + 1) % 3] * minIntersectionT;
        intersectionBarycentricCoords[(intersectionVertex + 1) % 3] = lNext;
        intersectionBarycentricCoords[(intersectionVertex + 2) % 3] = 1.0f - lNext;

        LogMessage("  Intersection b-coords: ", intersectionBarycentricCoords);

        newPosition = triangles.FromBarycentricCoordinates(
            intersectionBarycentricCoords,
            currentTriangle,
            vertices);

        particles.SetPosition(particleIndex, newPosition);
        state->ConstrainedState->CurrentTriangleBarycentricCoords = intersectionBarycentricCoords;

        // Check if edge is floor
        newEdgeIndex = triangles.GetSubEdges(currentTriangle).EdgeIndices[intersectionEdge];
        if (mModel->GetMesh().GetEdges().GetSurfaceType(newEdgeIndex) == SurfaceType::Floor
            && !isGhost)
        {
            //
            // Impact
            //

            LogMessage("  Impact (pre)");

            // Return, we'll then complete since we are on an edge
            return false;
        }
    }
    else
    {
        //
        // Cannot find intersection
        //
        // Must be on an edge
        //

        LogMessage("  No intersection at all, staying put");

        newPosition = particles.GetPosition(particleIndex);
        newEdgeIndex = triangles.GetSubEdges(currentTriangle).EdgeIndices[currentEdge];
    }

    // Find opposite triangle
    ElementIndex const oppositeTriangle = edges.GetOppositeTriangle(newEdgeIndex, state->ConstrainedState->CurrentTriangle);
    if (oppositeTriangle == NoneElementIndex)
    {
        //
        // Become free
        //

        LogMessage("  Becoming free");

        state->ConstrainedState.reset();

        return false;
    }

    // Move to edge of opposite triangle
    state->ConstrainedState->CurrentTriangle = oppositeTriangle;
    // TODO: derive from intersection b-coords
    state->ConstrainedState->CurrentTriangleBarycentricCoords = triangles.ToBarycentricCoordinates(
        newPosition,
        state->ConstrainedState->CurrentTriangle,
        vertices);

    LogMessage("  Moved to pos wrt current triangle ", state->ConstrainedState->CurrentTriangle, ": ", state->ConstrainedState->CurrentTriangleBarycentricCoords);

    return false;
}