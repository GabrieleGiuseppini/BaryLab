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
        return true;
    }

    //
    // We must be in constrained regime now
    //

    assert(state->ConstrainedState.has_value());

    ElementIndex const currentTriangle = state->ConstrainedState->CurrentTriangle;

    // Calculate trajectory == Target - CurrentPost
    vec2f const trajectory = state->TargetPosition - particles.GetPosition(particleIndex);

    //
    // If we are on a floor edge and we're moving against it, we're done;
    // Note: we ensure later that we move forcibly to a zero bary-coord when we hit an edge
    //

    // Note: if we are exactly at a vertex, we pick an arbitrary edge here
    ElementIndex tEdgeIndex = NoneElementIndex;
    if (state->ConstrainedState->CurrentTriangleBarycentricCoords.x == 0.0f)
    {
        tEdgeIndex = 1;
    }
    else if (state->ConstrainedState->CurrentTriangleBarycentricCoords.y == 0.0f)
    {
        tEdgeIndex = 2;
    }
    else if (state->ConstrainedState->CurrentTriangleBarycentricCoords.z == 0.0f)
    {
        tEdgeIndex = 0;
    }

    if (tEdgeIndex != NoneElementIndex 
        && edges.GetSurfaceType(triangles.GetSubEdges(currentTriangle).EdgeIndices[tEdgeIndex]) == SurfaceType::Floor)
    {
        // Allow to move like a ghost through all-floor triangles
        if (edges.GetSurfaceType(triangles.GetSubEdgeAIndex(currentTriangle)) != SurfaceType::Floor
            || edges.GetSurfaceType(triangles.GetSubEdgeBIndex(currentTriangle)) != SurfaceType::Floor
            || edges.GetSurfaceType(triangles.GetSubEdgeCIndex(currentTriangle)) != SurfaceType::Floor)
        {
            // Calculate pseudonormal, considering that we are *inside* the triangle
            // (points outside)
            vec2f const edgePNormal = (
                vertices.GetPosition(triangles.GetVertexIndices(currentTriangle)[(tEdgeIndex + 1) % 3])
                - vertices.GetPosition(triangles.GetVertexIndices(currentTriangle)[tEdgeIndex])
                ).to_perpendicular();

            if (edgePNormal.dot(trajectory) >= 0.0f) // Trajectory leaves this edge
            {
                //
                // We're going against the floor - floor impact
                //
                // (actually we've impacted at the previous step, here we just realize it)
                //

                LogMessage("On edge ", tEdgeIndex);

                return true;
            }            
        }
    }

    //
    // Analyze target position
    //

    // Calculate barycentric coordinates of target position wrt current triangle
    vec3f const targetBarycentricCoords = triangles.ToBarycentricCoordinates(
        state->TargetPosition,
        currentTriangle,
        vertices);

    LogMessage("Target pos wrt current triangle: ", targetBarycentricCoords);

    //
    // Calculate nearest intersection with triangle's edges
    //
    // The trajectory can be parameterized as (1 − t)*P + t*T, with P and T being the
    // particle and target points, in either coordinate system. By expressing P and T
    // in barycentric coordinates, the touch/cross point of the line with each edge
    // is found by imposing the parameterized trajectory component to be zero, 
    // yielding ti = lpi/(lpi - lti)
    //

    ElementIndex intersectionVertex = NoneElementIndex;
    float minIntersectionT = std::numeric_limits<float>::max();

    for (ElementIndex vi = 0; vi < 3; ++vi)
    {
        // TODOTEST: do float_arr

        float const den = state->ConstrainedState->CurrentTriangleBarycentricCoords[vi] - targetBarycentricCoords[vi];
        float const t = (den == 0.0f) // TODO: with epsilon (GameMath)
            ? std::numeric_limits<float>::max() // Parallel, meets at infinity
            : state->ConstrainedState->CurrentTriangleBarycentricCoords[vi] / den;

        if (t > 0.0f)
        {
            // Meets ahead - in the direction of trajectory
            if (t < minIntersectionT)
            {
                intersectionVertex = vi;
                minIntersectionT = t;
            }
        }
        // TODO: if t==0.0f?
    }

    assert(minIntersectionT >= 0.0f); // Meets ahead - in the direction of trajectory

    if (minIntersectionT > 1.0f)
    {
        // No intersection before end of trajectory =>
        // trajectory does not touch any edge before end of trajectory =>
        // end-of-trajectory is internal
        assert(targetBarycentricCoords.x > 0.0f && targetBarycentricCoords.x < 1.0f);
        assert(targetBarycentricCoords.y > 0.0f && targetBarycentricCoords.y < 1.0f);
        assert(targetBarycentricCoords.z > 0.0f && targetBarycentricCoords.z < 1.0f);

        LogMessage("  No intersection before end of trajectory");

        // ...move to target and we're done
        particles.SetPosition(particleIndex, state->TargetPosition);
        return true;
    }

    LogMessage("  Intersection on edge ", (intersectionVertex + 1) % 3, " (", minIntersectionT, ")");

    //
    // Trajectory intersects an edge before end-of-trajectory
    //

    // Calculate intersection's barycentric coordinates
    vec3f intersectionBarycentricCoords;    
    intersectionBarycentricCoords[intersectionVertex] = 0.0f;
    float const lNext = // Barycentric coord of next vertex at intersection
        state->ConstrainedState->CurrentTriangleBarycentricCoords[(intersectionVertex + 1) % 3] * (1.0f - minIntersectionT)
        + targetBarycentricCoords[(intersectionVertex + 1) % 3] * minIntersectionT;    
    intersectionBarycentricCoords[(intersectionVertex + 1) % 3] = lNext;
    intersectionBarycentricCoords[(intersectionVertex + 2) % 3] = 1.0f - lNext;

    LogMessage("  Intersection b-coords: ", intersectionBarycentricCoords);

    ////// TODOOLD
    ////bool const isTargetStrictlyInsideX = (targetBarycentricCoords.x > 0.0f && targetBarycentricCoords.x < 1.0f);
    ////bool const isTargetStrictlyInsideY = (targetBarycentricCoords.y > 0.0f && targetBarycentricCoords.y < 1.0f);
    ////bool const isTargetStrictlyInsideZ = (targetBarycentricCoords.z > 0.0f && targetBarycentricCoords.z < 1.0f);

    ////if (isTargetStrictlyInsideX && isTargetStrictlyInsideY && isTargetStrictlyInsideZ)
    ////{
    ////    // Strictly inside triangle...
    ////    // ...move to target and we're done
    ////    particles.SetPosition(particleIndex, state->TargetPosition);
    ////    return true;
    ////}

    ////// The particle is either on an edge/vertex, or crosses it

    //////
    ////// Find edge and intersection point
    //////
    ////// The trajectory can be parameterized as (1 − t)*P + t*T, with P and T being the
    ////// particle and target points, in either coordinate system. By expressing P and T
    ////// in barycentric coordinates, the touch/cross point of the line with each edge
    ////// is found by imposing the parameterized trajectory component to be zero, 
    ////// yielding ti = lpi/(lpi - lti)
    //////
    
    ////ElementIndex tIntersectionEdgeIndex;
    ////vec3f intersectionBarycentricCoords;

    ////// TODOTEST
    ////tIntersectionEdgeIndex = 0;
    ////intersectionBarycentricCoords = vec3f::zero();

    ////if (targetBarycentricCoords.x <= 0.0f)
    ////{
    ////    if (targetBarycentricCoords.y <= 0.0f)
    ////    {
    ////        // B-C or C-A edge

    ////        assert(targetBarycentricCoords.z >= 0.0f);

    ////        LogMessage("  Touch/cross B-C or C-A");
    ////        
    ////        // TODO
    ////        assert(false);
    ////    }
    ////    else if (targetBarycentricCoords.z <= 0.0f)
    ////    {
    ////        // B-C or A-B edge

    ////        assert(targetBarycentricCoords.y >= 0.0f);

    ////        LogMessage("  Touch/cross B-C or A-B");

    ////        // TODO
    ////        assert(false);
    ////    }
    ////    else
    ////    {
    ////        assert(targetBarycentricCoords.x <= 0.0f);
    ////        assert(targetBarycentricCoords.y > 0.0f);
    ////        assert(targetBarycentricCoords.z > 0.0f);

    ////        // B-C edge

    ////        LogMessage("  Touch/cross B-C");

    ////        tIntersectionEdgeIndex = 1;

    ////        float const den = state->ConstrainedState->CurrentTriangleBarycentricCoords.x - targetBarycentricCoords.x;
    ////        if (den == 0.0f)
    ////        {
    ////            // Parallel to B-C
    ////            // TODO: intersectionT is +inf - good for comparisons, but then we have problems with calculating barycentric coords
    ////            assert(false);
    ////        }
    ////        float const intersectionT = state->ConstrainedState->CurrentTriangleBarycentricCoords.x / den;
    ////        // TODO: use new formula w/trajectory
    ////        intersectionBarycentricCoords = vec3f(
    ////            0.0f,
    ////            state->ConstrainedState->CurrentTriangleBarycentricCoords.y * (1.0f - intersectionT) + targetBarycentricCoords.y * intersectionT,
    ////            state->ConstrainedState->CurrentTriangleBarycentricCoords.z * (1.0f - intersectionT) + targetBarycentricCoords.z * intersectionT);
    ////    }
    ////}
    ////else if (targetBarycentricCoords.y <= 0.0f)
    ////{
    ////    assert(targetBarycentricCoords.x > 0.0f);

    ////    if (targetBarycentricCoords.z <= 0.0f)
    ////    {
    ////        // C-A or A-B edge

    ////        LogMessage("  Touch/cross A-B or C-A");

    ////        // TODO
    ////        assert(false);
    ////    }
    ////    else
    ////    {
    ////        assert(targetBarycentricCoords.x > 0.0f);
    ////        assert(targetBarycentricCoords.y <= 0.0f);
    ////        assert(targetBarycentricCoords.z > 0.0f);

    ////        // C-A edge

    ////        LogMessage("  Touch/cross C-A");

    ////        tIntersectionEdgeIndex = 2;

    ////        float const den = state->ConstrainedState->CurrentTriangleBarycentricCoords.y - targetBarycentricCoords.y;
    ////        if (den == 0.0f)
    ////        {
    ////            // Parallel to C-A
    ////            // TODO
    ////            assert(false);
    ////        }

    ////        float const intersectionT = state->ConstrainedState->CurrentTriangleBarycentricCoords.y / den;
    ////        // TODO: use new formula w/trajectory
    ////        intersectionBarycentricCoords = vec3f(
    ////            state->ConstrainedState->CurrentTriangleBarycentricCoords.x * (1.0f - intersectionT) + targetBarycentricCoords.x * intersectionT,
    ////            0.0f,
    ////            state->ConstrainedState->CurrentTriangleBarycentricCoords.z * (1.0f - intersectionT) + targetBarycentricCoords.z * intersectionT);
    ////    }
    ////}
    ////else
    ////{
    ////    assert(targetBarycentricCoords.x > 0.0f);
    ////    assert(targetBarycentricCoords.y > 0.0f);
    ////    assert(targetBarycentricCoords.z <= 0.0f);

    ////    // A-B edge

    ////    LogMessage("  Touch/cross A-B");

    ////    tIntersectionEdgeIndex = 0;

    ////    float const den = state->ConstrainedState->CurrentTriangleBarycentricCoords.z - targetBarycentricCoords.z;
    ////    if (den == 0.0f)
    ////    {
    ////        // Parallel to A-B
    ////        // TODO
    ////        assert(false);
    ////    }

    ////    float const intersectionT = state->ConstrainedState->CurrentTriangleBarycentricCoords.z / den;
    ////    // TODO: use new formula w/trajectory
    ////    intersectionBarycentricCoords = vec3f(
    ////        state->ConstrainedState->CurrentTriangleBarycentricCoords.x * (1.0f - intersectionT) + targetBarycentricCoords.x * intersectionT,
    ////        state->ConstrainedState->CurrentTriangleBarycentricCoords.y * (1.0f - intersectionT) + targetBarycentricCoords.y * intersectionT,
    ////        0.0f);
    ////}

    ////// Normalize barycentric coords
    ////// TODO: doesn't work if we've just set z to zero
    ////intersectionBarycentricCoords.z = 1.0f - intersectionBarycentricCoords.x - intersectionBarycentricCoords.y;

    //
    // Move to intersection point
    //    

    vec2f const intersectionPosition = triangles.FromBarycentricCoordinates(
        intersectionBarycentricCoords,
        currentTriangle,
        vertices);

    particles.SetPosition(particleIndex, intersectionPosition);
    state->ConstrainedState->CurrentTriangleBarycentricCoords = intersectionBarycentricCoords;

    // Check if edge is floor
    ElementIndex const edgeIndex = triangles.GetSubEdges(currentTriangle).EdgeIndices[(intersectionVertex + 1) % 3];
    if (mModel->GetMesh().GetEdges().GetSurfaceType(edgeIndex) == SurfaceType::Floor)
    {
        //
        // Impact
        //

        LogMessage("  Impact");

        // Return, we'll then complete since we are on an edge
        return false;
    }

    // TODOHERE

    return false;
}