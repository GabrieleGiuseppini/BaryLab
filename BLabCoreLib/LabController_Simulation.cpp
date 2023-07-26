/***************************************************************************************
 * Original Author:		Gabriele Giuseppini
 * Created:				2023-07-23
 * Copyright:			Gabriele Giuseppini  (https://github.com/GabrieleGiuseppini)
 ***************************************************************************************/
#include "LabController.h"

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
                // Reset state
                particles.GetState(p).reset();
                mCurrentParticleTrajectory.reset();
            }
        }
    }
}

bool LabController::UpdateParticleState(ElementIndex particleIndex)
{
    // TODOHERE
    (void)particleIndex;
    return true;
}