/***************************************************************************************
 * Original Author:		Gabriele Giuseppini
 * Created:				2023-07-23
 * Copyright:			Gabriele Giuseppini  (https://github.com/GabrieleGiuseppini)
 ***************************************************************************************/
#include "LabController.h"

#include "Log.h"

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

    // If we are on a floor spring, we're done;
    // note: we leverage that we move forcibly to a zero b-coord when we hit an edge
    if (state->ConstrainedState->CurrentTriangleBarycentricCoords.x == 0.0f
        || state->ConstrainedState->CurrentTriangleBarycentricCoords.y == 0.0f
        || state->ConstrainedState->CurrentTriangleBarycentricCoords.z == 0.0f)
    {
        return true;
    }

    // Calculate barycentric coordinates of target position wrt current triangle
    vec3f const targetBarycentricCoords = triangles.ToBarycentricCoordinates(
        state->TargetPosition,
        state->ConstrainedState->CurrentTriangle,
        vertices);

    LogMessage("Target pos wrt current triangle: ", targetBarycentricCoords);

    //
    // Analyze target position
    //

    bool const isTargetStrictlyInsideX = (targetBarycentricCoords.x > 0.0f && targetBarycentricCoords.x < 1.0f);
    bool const isTargetStrictlyInsideY = (targetBarycentricCoords.y > 0.0f && targetBarycentricCoords.y < 1.0f);
    bool const isTargetStrictlyInsideZ = (targetBarycentricCoords.z > 0.0f && targetBarycentricCoords.z < 1.0f);

    if (isTargetStrictlyInsideX && isTargetStrictlyInsideY && isTargetStrictlyInsideZ)
    {
        // Strictly inside triangle...
        // ...move to target and we're done
        particles.SetPosition(particleIndex, state->TargetPosition);
        return true;
    }

    // The particle is either on an edge/vertex, or crosses it

    //
    // Find edge and intersection point
    //
    // The trajectory can be parameterized as (1 − t)*P + t*T, with P and T being the
    // particle and target points, in either coordinate system. By expressing P and T
    // in barycentric coordinates, the touch/cross point of the line with each edge
    // is found by imposing the parameterized trajectory component to be zero, 
    // yielding ti = lpi/(lpi - lti)
    //

    ElementIndex tEdgeIndex;
    float intersectionT;

    // TODOTEST
    tEdgeIndex = 0;
    intersectionT = 0.0f;

    if (targetBarycentricCoords.x <= 0.0f)
    {
        if (targetBarycentricCoords.y <= 0.0f)
        {
            // B-C or C-A edge

            assert(targetBarycentricCoords.z >= 0.0f);

            LogMessage("  Touch/cross B-C or C-A");
            
            // TODO
            assert(false);
        }
        else if (targetBarycentricCoords.z <= 0.0f)
        {
            // B-C or A-B edge

            assert(targetBarycentricCoords.y >= 0.0f);

            LogMessage("  Touch/cross B-C or A-B");

            // TODO
            assert(false);
        }
        else
        {
            // B-C edge

            LogMessage("  Touch/cross B-C");

            tEdgeIndex = 1;

            float const den = state->ConstrainedState->CurrentTriangleBarycentricCoords.x - targetBarycentricCoords.x;
            if (den == 0.0f)
            {
                // Parallel to B-C
                // TODO
                assert(false);
            }

            intersectionT = state->ConstrainedState->CurrentTriangleBarycentricCoords.x / den;
        }
    }
    else if (targetBarycentricCoords.y <= 0.0f)
    {
        assert(targetBarycentricCoords.x > 0.0f);

        if (targetBarycentricCoords.z <= 0.0f)
        {
            // C-A or A-B edge

            LogMessage("  Touch/cross A-B or C-A");

            // TODO
            assert(false);
        }
        else
        {
            // C-A edge

            LogMessage("  Touch/cross C-A");

            tEdgeIndex = 2;

            float const den = state->ConstrainedState->CurrentTriangleBarycentricCoords.y - targetBarycentricCoords.y;
            if (den == 0.0f)
            {
                // Parallel to C-A
                // TODO
                assert(false);
            }

            intersectionT = state->ConstrainedState->CurrentTriangleBarycentricCoords.y / den;
        }
    }
    else
    {
        assert(targetBarycentricCoords.x > 0.0f);
        assert(targetBarycentricCoords.y > 0.0f);
        assert(targetBarycentricCoords.z <= 0.0f);

        // A-B edge

        LogMessage("  Touch/cross A-B");

        tEdgeIndex = 0;

        float const den = state->ConstrainedState->CurrentTriangleBarycentricCoords.z - targetBarycentricCoords.z;
        if (den == 0.0f)
        {
            // Parallel to A-B
            // TODO
            assert(false);
        }

        intersectionT = state->ConstrainedState->CurrentTriangleBarycentricCoords.z / den;
    }

    // TODO: since we want to move to a position that is unambiguously over an edge, we need to 
    // construct the intersection barycentric coords from above

    vec3f const intersectionBarycentricCoords =
        state->ConstrainedState->CurrentTriangleBarycentricCoords * (1.0f - intersectionT)
        + targetBarycentricCoords * intersectionT;

    //
    // Move to intersection point
    //    

    vec2f const intersectionPosition = triangles.FromBarycentricCoordinates(
        intersectionBarycentricCoords,
        state->ConstrainedState->CurrentTriangle,
        vertices);

    particles.SetPosition(particleIndex, intersectionPosition);
    state->ConstrainedState->CurrentTriangleBarycentricCoords = intersectionBarycentricCoords;

    // Check if edge is floor
    ElementIndex const edgeIndex = triangles.GetSubEdges(state->ConstrainedState->CurrentTriangle).EdgeIndices[tEdgeIndex];
    if (mModel->GetMesh().GetEdges().GetSurfaceType(edgeIndex) == SurfaceType::Floor)
    {
        //
        // Impact
        //

        // Return, we'll then complete since we are on an edge
        return false;
    }

    // TODOHERE

    return false;
}