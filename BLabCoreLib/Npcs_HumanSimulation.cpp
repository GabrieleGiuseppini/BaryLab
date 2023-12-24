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
	NpcParticles & particles,
	Mesh const & mesh) const
{
	// TODOHERE
	(void)mesh;

	// Reset voluntary physics

	particles.SetVoluntaryForces(primaryParticleState.ParticleIndex, vec2f::zero());
	particles.SetVoluntaryForces(secondaryParticleState.ParticleIndex, vec2f::zero());

	particles.SetVoluntarySuperimposedDisplacement(primaryParticleState.ParticleIndex, vec2f::zero());
	particles.SetVoluntarySuperimposedDisplacement(secondaryParticleState.ParticleIndex, vec2f::zero());

	particles.SetVoluntaryVelocity(primaryParticleState.ParticleIndex, vec2f::zero());
	particles.SetVoluntaryVelocity(secondaryParticleState.ParticleIndex, vec2f::zero());

	// Return state

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
	// TODOTEST
	(void)mesh;

	float const ToRisingConvergenceRate = 0.067f;
	float const ToWalkingConvergenceRate = 0.09f;
	float constexpr MaxRelativeVelocityForEquilibrium = 2.9f; // So high because we slip a lot while we try to stand up, and thus need to be immune to ourselves

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
				&& (!secondaryParticleState.ConstrainedState.has_value()
					|| secondaryParticleState.ConstrainedState->MeshRelativeVelocity.length() < MaxRelativeVelocityForEquilibrium))
			{
				stateCondition = 1.0f;
			}

			humanState.TargetStateValue = stateCondition;
			
			break;
		}

		case StateType::HumanNpcStateType::BehaviorType::Constrained_Rising:
		case StateType::HumanNpcStateType::BehaviorType::Constrained_Equilibrium:
		case StateType::HumanNpcStateType::BehaviorType::Constrained_Walking:
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

				mParticles.SetVoluntarySuperimposedDisplacement(
					primaryParticleState.ParticleIndex,
					vec2f::zero());

				mParticles.SetVoluntarySuperimposedDisplacement(
					secondaryParticleState.ParticleIndex,
					vec2f::zero());

				mParticles.SetVoluntaryVelocity(
					primaryParticleState.ParticleIndex,
					vec2f::zero());

				mParticles.SetVoluntaryVelocity(
					secondaryParticleState.ParticleIndex,
					vec2f::zero());

				humanState.CurrentEquilibriumTorqueMagnitude = 0.0f;
				humanState.CurrentWalkingMagnitude = 0.0f;

				mEventDispatcher.OnHumanNpcBehaviorChanged("Free_KnockedOut");

				break;
			}

			if (humanState.CurrentBehavior == StateType::HumanNpcStateType::BehaviorType::Constrained_Equilibrium)
			{
				// Advance towards walking

				humanState.CurrentStateValue += (humanState.TargetStateValue - humanState.CurrentStateValue) * ToWalkingConvergenceRate;

				publishStateQuantity = std::make_tuple("CurrentStateValue", std::to_string(humanState.CurrentStateValue));

				if (IsAtTarget(humanState.CurrentStateValue, 1.0f))
				{
					// Transition

					humanState.TransitionToState(
						StateType::HumanNpcStateType::BehaviorType::Constrained_Walking,
						0.0f,
						0.0f);

					humanState.CurrentWalkingMagnitude = 0.1f;

					// Keep torque

					mEventDispatcher.OnHumanNpcBehaviorChanged("Constrained_Walking");

					break;
				}
			}

			// Check conditions to stay

			bool stateCondition = false;

			if (primaryParticleState.ConstrainedState.has_value()
				&& IsOnEdge(primaryParticleState.ConstrainedState->CurrentTriangleBarycentricCoords)
				&& primaryParticleState.ConstrainedState->MeshRelativeVelocity.length() < MaxRelativeVelocityForEquilibrium)
			{
				stateCondition = true;
			}

			if (!stateCondition)
			{
				// Transition

				LogMessage("Going to Constrained_KnockedOut; primary's barycentric coords: ", 
					primaryParticleState.ConstrainedState.has_value() ? primaryParticleState.ConstrainedState->CurrentTriangleBarycentricCoords.toString() : "",
					" primary's relative velocity: ", primaryParticleState.ConstrainedState->MeshRelativeVelocity.length());

				humanState.TransitionToState(
					StateType::HumanNpcStateType::BehaviorType::Constrained_KnockedOut,
					0.0f,
					0.0f);

				mParticles.SetVoluntaryForces(
					secondaryParticleState.ParticleIndex,
					vec2f::zero());

				mParticles.SetVoluntarySuperimposedDisplacement(
					primaryParticleState.ParticleIndex,
					vec2f::zero());

				mParticles.SetVoluntarySuperimposedDisplacement(
					secondaryParticleState.ParticleIndex,
					vec2f::zero());

				mParticles.SetVoluntaryVelocity(
					primaryParticleState.ParticleIndex,
					vec2f::zero());

				mParticles.SetVoluntaryVelocity(
					secondaryParticleState.ParticleIndex,
					vec2f::zero());

				humanState.CurrentEquilibriumTorqueMagnitude = 0.0f;
				humanState.CurrentWalkingMagnitude = 0.0f;

				mEventDispatcher.OnHumanNpcBehaviorChanged("Constrained_KnockedOut");

				break;
			}

			// Maintain equilibrium

			if (!MaintainAndCheckHumanEquilibrium(
				humanState,
				primaryParticleState.ParticleIndex,
				secondaryParticleState.ParticleIndex,
				mParticles,
				labParameters))
			{
				// Transition

				humanState.TransitionToState(
					StateType::HumanNpcStateType::BehaviorType::Constrained_KnockedOut,
					0.0f,
					0.0f);

				mParticles.SetVoluntaryForces(
					secondaryParticleState.ParticleIndex,
					vec2f::zero());

				mParticles.SetVoluntarySuperimposedDisplacement(
					primaryParticleState.ParticleIndex,
					vec2f::zero());

				mParticles.SetVoluntarySuperimposedDisplacement(
					secondaryParticleState.ParticleIndex,
					vec2f::zero());

				mParticles.SetVoluntaryVelocity(
					primaryParticleState.ParticleIndex,
					vec2f::zero());

				mParticles.SetVoluntaryVelocity(
					secondaryParticleState.ParticleIndex,
					vec2f::zero());

				humanState.CurrentEquilibriumTorqueMagnitude = 0.0f;
				humanState.CurrentWalkingMagnitude = 0.0f;

				mEventDispatcher.OnHumanNpcBehaviorChanged("Constrained_KnockedOut");

				break;
			}

			if (humanState.CurrentBehavior == StateType::HumanNpcStateType::BehaviorType::Constrained_Rising)
			{
				// Check if reached alignment (note: here so that we may keep torque as we'll be transitioning to Equilibrium)

				float const alignment = CalculateVerticalAlignment(primaryParticleState.ParticleIndex, secondaryParticleState.ParticleIndex, mParticles);

				publishStateQuantity = std::make_tuple("Alignment", std::to_string(alignment));

				if (AreAlmostEqual(alignment, 1.0f, 0.004f))
				{
					// Transition

					humanState.TransitionToState(
						StateType::HumanNpcStateType::BehaviorType::Constrained_Equilibrium,
						0.0f,
						0.0f);

					mEventDispatcher.OnHumanNpcBehaviorChanged("Constrained_Equilibrium");
				}
			}
			else if (humanState.CurrentBehavior == StateType::HumanNpcStateType::BehaviorType::Constrained_Equilibrium)
			{
				// We are well in a state to advance to walking

				humanState.TargetStateValue = 1.0f;
			}
			else
			{
				assert(humanState.CurrentBehavior == StateType::HumanNpcStateType::BehaviorType::Constrained_Walking);

				// Impart walk displacement & run walking state machine
				RunWalkingHumanStateMachine(
					humanState,
					primaryParticleState,
					secondaryParticleState,
					mesh,
					labParameters);
			}

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

bool Npcs::MaintainAndCheckHumanEquilibrium(
	StateType::HumanNpcStateType & humanState,
	ElementIndex primaryParticleIndex,
	ElementIndex secondaryParticleIndex,
	NpcParticles & particles,
	LabParameters const & labParameters)
{
	//
	// Calculates torque to maintain equilibrium, and makes sure we are not falling out of equilibrium
	//

	vec2f const headPosition = particles.GetPosition(secondaryParticleIndex);
	vec2f const feetPosition = particles.GetPosition(primaryParticleIndex);

	vec2f const humanVector = headPosition - feetPosition;

	// Calculate CW angle between head and vertical (pointing up);
	// positive when human is CW wrt vertical
	float const staticDisplacementAngleCW = (-LabParameters::GravityDir).angleCw(humanVector);

	// Calculate CW angle that would be rotated by (relative to feet) velocity alone;
	// positive when new position is CW wrt old
	vec2f const headPositionAfterVelocity =
		headPosition
		+ (particles.GetVelocity(secondaryParticleIndex) - particles.GetVelocity(primaryParticleIndex)) * LabParameters::SimulationTimeStepDuration;
	float const velocityAngleCW = humanVector.angleCw(headPositionAfterVelocity - feetPosition);

	//
	// Check whether we are still in equulibrium
	//
	// We lose equilibrium if HumanVector is outside of -alpha->alpha sector around vertical, with non-negligible rotation velocity towards outside of sector
	//

	float constexpr MaxStaticAngleForEquilibrium = Pi<float> / 7.0f;

	if (std::abs(staticDisplacementAngleCW) >= MaxStaticAngleForEquilibrium
		&& std::abs(velocityAngleCW) > 0.01f
		&& staticDisplacementAngleCW * velocityAngleCW > 0.0f) // Equal signs
	{
		LogMessage("Losing equilibrium because: StaticDisplacementAngleCW=", staticDisplacementAngleCW, " (Max=", MaxStaticAngleForEquilibrium,
			") VelocityAngleCW=", velocityAngleCW);

		return false;
	}

	// Converge equilibrium torque
	float const ToFullEquilibriumTorqueConvergenceRate = 0.09f;
	humanState.CurrentEquilibriumTorqueMagnitude += (1.0f - humanState.CurrentEquilibriumTorqueMagnitude) * ToFullEquilibriumTorqueConvergenceRate;

	(void)labParameters;

	return true;
}

void Npcs::RunWalkingHumanStateMachine(
	StateType::HumanNpcStateType & humanState,
	StateType::NpcParticleStateType const & primaryParticleState,
	StateType::NpcParticleStateType const & secondaryParticleState,
	Mesh const & mesh,
	LabParameters const & labParameters)
{
	// TODOHERE
	(void)secondaryParticleState;
	(void)mesh;

	// Advance towards 1.0
	// TODO: not needed to be public
	humanState.CurrentWalkingMagnitude = std::min(1.0f, humanState.CurrentWalkingMagnitude + (1.0f - humanState.CurrentWalkingMagnitude) * 0.03f);
	
	LogMessage("        walking: ", humanState.CurrentWalkingMagnitude);

	// TODOTEST
	////// Impart walking velocity to primary particle
	////mParticles.SetVoluntaryVelocity(
	////	primaryParticleState.ParticleIndex,
	////	vec2f(humanState.CurrentFaceDirectionX * labParameters.HumanNpcWalkingSpeed * humanState.CurrentWalkingMagnitude, 0.0f));
	(void)labParameters;
	(void)primaryParticleState;
}
