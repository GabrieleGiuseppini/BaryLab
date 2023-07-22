/***************************************************************************************
 * Original Author:		Gabriele Giuseppini
 * Created:				2023-07-23
 * Copyright:			Gabriele Giuseppini  (https://github.com/GabrieleGiuseppini)
 ***************************************************************************************/
#include "LabController.h"

void LabController::UpdateSimulation(LabParameters const & labParameters)
{
    float const dt = LabParameters::SimulationTimeStepDuration;

    //
    // Particles
    //

    auto & particles = mModel->GetParticles();

    float const particleMass = LabParameters::ParticleMass * labParameters.MassAdjustment;

    for (auto const & p : particles)
    {
        vec2f const forces = particles.GetWorldForce(p) * labParameters.GravityAdjustment;

        vec2f const deltaPos =
            particles.GetVelocity(p) * dt
            + forces / LabParameters::ParticleMass * dt * dt;

        particles.SetPosition(p, particles.GetPosition(p) + deltaPos);
        particles.SetVelocity(p, deltaPos / dt);
    }
}
