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

    auto & state = mModel->GetParticles().GetState(particleIndex);

    assert(state.has_value());

    if (!state->ConstrainedState)
    {
        // Free regime...

        // ...move to target position directly
        mModel->GetParticles().SetPosition(particleIndex, state->TargetPosition);
    }

    if (mModel->GetParticles().GetPosition(particleIndex) == state->TargetPosition)
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
    vec3f const targetBarycentricCoords = mModel->GetMesh().GetTriangles().ToBarycentricCoordinates(
        state->TargetPosition,
        state->ConstrainedState->CurrentTriangle,
        mModel->GetMesh().GetVertices());

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
        mModel->GetParticles().SetPosition(particleIndex, state->TargetPosition);
        return true;
    }

    // TODOHERE
    return false;
}