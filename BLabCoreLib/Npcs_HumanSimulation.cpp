/***************************************************************************************
 * Original Author:		Gabriele Giuseppini
 * Created:				2023-10-20
 * Copyright:			Gabriele Giuseppini  (https://github.com/GabrieleGiuseppini)
 ***************************************************************************************/
#include "Npcs.h"

#include "BLabMath.h"
#include "Log.h"

#include <cmath>

namespace /*anonymous*/ {

	bool IsAtTarget(float currentValue, float targetValue)
	{
		return std::abs(targetValue - currentValue) < 0.01f;
	}

	bool IsOnEdge(vec3f const & barycentricCoords)
	{
		return barycentricCoords.x == 0.0f
			|| barycentricCoords.y == 0.0f
			|| barycentricCoords.z == 0.0f;
	}

	// Head->Feet
	vec2f CalculateHumanVector(ElementIndex primaryParticleIndex, ElementIndex secondaryParticleIndex, NpcParticles const & particles)
	{
		return particles.GetPosition(primaryParticleIndex) - particles.GetPosition(secondaryParticleIndex);
	}

	float CalculateVerticalAlignment(vec2f const & humanVector)
	{
		return humanVector.normalise().dot(LabParameters::GravityDir);
	}

	float CalculateVerticalAlignment(ElementIndex primaryParticleIndex, ElementIndex secondaryParticleIndex, NpcParticles const & particles)
	{
		return CalculateVerticalAlignment(CalculateHumanVector(primaryParticleIndex, secondaryParticleIndex, particles));
	}
}

Npcs::StateType::HumanNpcStateType Npcs::InitializeHuman(
	StateType::NpcParticleStateType const & primaryParticleState,
	StateType::NpcParticleStateType const & secondaryParticleState,
	Mesh const & mesh) const
{
	// TODOHERE
	(void)primaryParticleState;
	(void)secondaryParticleState;
	(void)mesh;

	if (!primaryParticleState.ConstrainedState.has_value()
		&& !secondaryParticleState.ConstrainedState.has_value())
	{
		// Whole NPC is free
		mEventDispatcher.OnHumanNpcBehaviorChanged("Free_KnockedOut");
		return StateType::HumanNpcStateType(StateType::HumanNpcStateType::BehaviorType::Free_KnockedOut, 0.0f, 0.0f);
	}
	else
	{
		// NPC is constrained
		// TODO: check alignment and, if ok, start at Constrained_Equilibrium
		mEventDispatcher.OnHumanNpcBehaviorChanged("Constrained_KnockedOut");
		return StateType::HumanNpcStateType(StateType::HumanNpcStateType::BehaviorType::Constrained_KnockedOut, 0.0f, 0.0f);
	}
}

void Npcs::UpdateHuman(
	StateType::HumanNpcStateType & humanState,
	StateType::NpcParticleStateType const & primaryParticleState,
	StateType::NpcParticleStateType const & secondaryParticleState,
	Mesh const & mesh,
	LabParameters const & labParameters)
{
	float const ToRisingConvergenceRate = 0.05f;
	float constexpr MaxRelativeVelocityForEquilibrium = 5.0f;

	// TODOHERE
	(void)mesh;

	std::optional<std::tuple<std::string, std::string>> publishStateQuantity;

	bool const isFree = !primaryParticleState.ConstrainedState.has_value() && !secondaryParticleState.ConstrainedState.has_value();

	switch (humanState.CurrentBehavior)
	{
		case StateType::HumanNpcStateType::BehaviorType::Constrained_KnockedOut:
		{
			if (isFree)
			{
				// Transition

				humanState.TransitionToState(
					StateType::HumanNpcStateType::BehaviorType::Free_KnockedOut,
					0.0f,
					0.0f);

				mEventDispatcher.OnHumanNpcBehaviorChanged("Free_KnockedOut");

				break;
			}

			// Advance

			humanState.CurrentStateValue += (humanState.TargetStateValue - humanState.CurrentStateValue) * ToRisingConvergenceRate;

			publishStateQuantity = std::make_tuple("CurrentStateValue", std::to_string(humanState.CurrentStateValue));

			if (IsAtTarget(humanState.CurrentStateValue, 1.0f))
			{
				// Transition

				humanState.TransitionToState(
					StateType::HumanNpcStateType::BehaviorType::Constrained_Rising,
					0.0f,
					0.0f);

				mEventDispatcher.OnHumanNpcBehaviorChanged("Constrained_Rising");

				break;
			}

			// Check conditions for this state

			float stateCondition = 0.0f;

			if (primaryParticleState.ConstrainedState.has_value()
				&& IsOnEdge(primaryParticleState.ConstrainedState->CurrentTriangleBarycentricCoords)
				&& primaryParticleState.ConstrainedState->MeshRelativeVelocity.length() < MaxRelativeVelocityForEquilibrium
				&& secondaryParticleState.ConstrainedState.has_value()
				&& IsOnEdge(secondaryParticleState.ConstrainedState->CurrentTriangleBarycentricCoords)
				&& secondaryParticleState.ConstrainedState->MeshRelativeVelocity.length() < MaxRelativeVelocityForEquilibrium)
			{
				stateCondition = 1.0f;
			}

			humanState.TargetStateValue = stateCondition;
			
			break;
		}

		case StateType::HumanNpcStateType::BehaviorType::Constrained_Rising:
		{
			if (isFree)
			{
				// Transition

				humanState.TransitionToState(
					StateType::HumanNpcStateType::BehaviorType::Free_KnockedOut,
					0.0f,
					0.0f);

				mParticles.SetVoluntaryForces(
					secondaryParticleState.ParticleIndex,
					vec2f::zero());

				mEventDispatcher.OnHumanNpcBehaviorChanged("Free_KnockedOut");

				break;
			}

			// Check conditions to stay

			bool stateCondition = false;

			if (primaryParticleState.ConstrainedState.has_value()
				&& IsOnEdge(primaryParticleState.ConstrainedState->CurrentTriangleBarycentricCoords)
				&& primaryParticleState.ConstrainedState->MeshRelativeVelocity.length() < MaxRelativeVelocityForEquilibrium)
			{
				// Make sure secondary's velocity is in direction of becoming erected
				// TODOHERE: secondary: velocity is in direction of becoming erected
				stateCondition = true;
			}

			if (!stateCondition)
			{
				// Transition

				LogMessage("Going to Constrained_KnockedOut; primary's relative velocity: ", primaryParticleState.ConstrainedState->MeshRelativeVelocity.length());

				humanState.TransitionToState(
					StateType::HumanNpcStateType::BehaviorType::Constrained_KnockedOut,
					0.0f,
					0.0f);

				mParticles.SetVoluntaryForces(
					secondaryParticleState.ParticleIndex,
					vec2f::zero());

				mEventDispatcher.OnHumanNpcBehaviorChanged("Constrained_KnockedOut");

				break;
			}

			// Apply torque

			ApplyHumanNpcEquilibriumTorque(
				primaryParticleState.ParticleIndex,
				secondaryParticleState.ParticleIndex,
				mParticles,
				labParameters);

			// Check alignment (note: here so that we may keep torque as we'll be transitioning to Equilibrium)

			float const alignment = CalculateVerticalAlignment(primaryParticleState.ParticleIndex, secondaryParticleState.ParticleIndex, mParticles);

			publishStateQuantity = std::make_tuple("Alignment", std::to_string(alignment));

			if (AreAlmostEqual(alignment, 1.0f, 0.005f))
			{
				// Transition

				humanState.TransitionToState(
					StateType::HumanNpcStateType::BehaviorType::Constrained_Equilibrium,
					0.0f,
					0.0f);

				mEventDispatcher.OnHumanNpcBehaviorChanged("Constrained_Equilibrium");

				break;
			}

			break;
		}

		case StateType::HumanNpcStateType::BehaviorType::Constrained_Equilibrium:
		{
			if (isFree)
			{
				// Transition

				humanState.TransitionToState(
					StateType::HumanNpcStateType::BehaviorType::Free_KnockedOut,
					0.0f,
					0.0f);

				mEventDispatcher.OnHumanNpcBehaviorChanged("Free_KnockedOut");

				break;
			}

			// TODO: When enough time has passed in this state: transition to Walking

			// TODO: knocked off

			// Apply torque

			ApplyHumanNpcEquilibriumTorque(
				primaryParticleState.ParticleIndex,
				secondaryParticleState.ParticleIndex,
				mParticles,
				labParameters);
			
			break;
		}

		case StateType::HumanNpcStateType::BehaviorType::Constrained_Walking:
		{
			// TODOHERE
			break;
		}

		case StateType::HumanNpcStateType::BehaviorType::Free_KnockedOut:
		{
			// TODOHERE
			break;
		}
	}

	mEventDispatcher.OnHumanNpcStateQuantityChanged(publishStateQuantity);
}

void Npcs::ApplyHumanNpcEquilibriumTorque(	
	ElementIndex primaryParticleIndex,
	ElementIndex secondaryParticleIndex,
	NpcParticles & particles,
	LabParameters const & labParameters)
{
	vec2f const headPosition = particles.GetPosition(secondaryParticleIndex);
	vec2f const feetPosition = particles.GetPosition(primaryParticleIndex);

	vec2f const humanVector = headPosition - feetPosition;

	// Torque stiffness

	// Calculate CW angle between head and vertical (pointing up)
	// Positive when human is CW wrt vertical
	float const displacementAngleCW = (-LabParameters::GravityDir).angleCw(humanVector);
	
	// Velocity damping

	// New pos after integration of velocity
	vec2f const positionAfterVelocity =
		headPosition
		+ particles.GetVelocity(secondaryParticleIndex) * LabParameters::SimulationTimeStepDuration;

	// Delta angle of velocity
	// Positive when new position is CW wrt old
	float const velocityAngleCW = -(positionAfterVelocity - feetPosition).angleCw(humanVector);

	// Calculate total torque angle
	float const totalTorqueAngleCW = 
		-displacementAngleCW / 32.0f
		- velocityAngleCW * labParameters.HumanNpcEquilibriumTorqueDampingCoefficient;

	// Calculate displacement now
	vec2f const endPosition = humanVector.rotate(-totalTorqueAngleCW) + feetPosition;
	vec2f const torqueDisplacement = endPosition - headPosition;

	LogMessage("TODOHERE: DisplacementAngleCW:", displacementAngleCW, " VelocityAngleCW:", velocityAngleCW, " Total TorqueAngleCW:", totalTorqueAngleCW, " TorqueDisplacement:", torqueDisplacement);

	vec2f const torqueForce = 
		torqueDisplacement 
		* LabParameters::ParticleMass / (LabParameters::SimulationTimeStepDuration * LabParameters::SimulationTimeStepDuration)
		* labParameters.HumanNpcEquilibriumTorqueStiffnessCoefficient;

	mParticles.SetVoluntaryForces(
		secondaryParticleIndex,
		torqueForce);



	// TODOOLD
	////// TODOTEST
	//////float const torqueStrength = std::min(1.0f - alignment, 0.1f) / 0.1f;
	////float const torqueStrength = 1.0f;

	////// Calculate desired displacement
	////vec2f torqueDisplacement = humanVector.normalise().to_perpendicular() * torqueStrength * 0.04f * labParameters.HumanNpcRisingTorqueFactor;
	////if (humanVector.cross(LabParameters::GravityDir) > 0.0f)
	////	torqueDisplacement *= -1.0f;

	////// Check how much current velocity would already contribute to that
	////vec2f const velocityDisplacement = particles.GetVelocity(secondaryParticleIndex) * labParameters.SimulationTimeStepDuration;

	////// Adjust desired displacement
	////torqueDisplacement = torqueDisplacement - velocityDisplacement;

	////LogMessage("HumanVector(H->F): ", humanVector, " TorqueStrength: ", torqueStrength, " TorqueDisplacement: ", torqueDisplacement);

	////vec2f const torqueForce = torqueDisplacement * LabParameters::ParticleMass / (LabParameters::SimulationTimeStepDuration * LabParameters::SimulationTimeStepDuration);

	////mParticles.SetVoluntaryForces(
	////	secondaryParticleIndex,
	////	torqueForce);
}
