/***************************************************************************************
 * Original Author:		Gabriele Giuseppini
 * Created:				2023-10-20
 * Copyright:			Gabriele Giuseppini  (https://github.com/GabrieleGiuseppini)
 ***************************************************************************************/
#include "Npcs.h"

#include "BLabMath.h"
#include "Log.h"

#include <cmath>

Npcs::StateType::HumanNpcStateType Npcs::InitializeHuman(
	StateType::NpcParticleStateType const & primaryParticleState,
	StateType::NpcParticleStateType const & secondaryParticleState) const
{
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
	float currentSimulationTime,
	StateType::HumanNpcStateType & humanState,
	StateType::NpcParticleStateType const & primaryParticleState,
	StateType::NpcParticleStateType const & secondaryParticleState,
	Mesh const & mesh,
	LabParameters const & labParameters)
{
	// TODOTEST
	(void)currentSimulationTime;

	float const ToRisingConvergenceRate = 0.067f;
	float const ToWalkingConvergenceRate = 0.09f;
	float constexpr MaxRelativeVelocityForEquilibrium = 3.0f; // So high because we slip a lot while we try to stand up, and thus need to be immune to ourselves

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

				humanState.CurrentEquilibriumSoftTerminationDecision = 0.0f;

				mEventDispatcher.OnHumanNpcBehaviorChanged("Constrained_Rising");

				break;
			}

			// Check conditions for this state

			float stateCondition = 0.0f;

			if (primaryParticleState.ConstrainedState.has_value()
				&& IsOnFloorEdge(*primaryParticleState.ConstrainedState, mesh)
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

				humanState.CurrentWalkMagnitude = 0.0f;

				mEventDispatcher.OnHumanNpcBehaviorChanged("Free_KnockedOut");

				break;
			}

			bool const isOnEdge =
				primaryParticleState.ConstrainedState.has_value()
				&& IsOnFloorEdge(*primaryParticleState.ConstrainedState, mesh);

			if (humanState.CurrentBehavior == StateType::HumanNpcStateType::BehaviorType::Constrained_Equilibrium)
			{
				if (isOnEdge)
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

						humanState.CurrentWalkMagnitude = 0.0f;
						humanState.TargetWalkMagnitude = 1.0f;
						humanState.CurrentWalkFlipDecision = 0.0f;
						humanState.TargetWalkFlipDecision = 0.0f;
						humanState.TotalEdgeTraveledSinceWalkStart = 0.0f;

						// Keep torque

						mEventDispatcher.OnHumanNpcBehaviorChanged("Constrained_Walking");

						break;
					}
				}
			}

			// Check conditions to stay & maintain equilibrium

			bool isStateMaintained;
			if (!isOnEdge)
			{
				// When walking, we want to be a bit more tolerant about "losing the edge"

				// This is a quite important parameter: it's the duration through which we tolerate temporarily losing contact
				// with the ground
				float constexpr ToTerminateEquilibriumConvergenceRate = 0.25f;

				// Advance
				humanState.CurrentEquilibriumSoftTerminationDecision += (1.0f - humanState.CurrentEquilibriumSoftTerminationDecision) * ToTerminateEquilibriumConvergenceRate;

				// Check if enough
				if (IsAtTarget(humanState.CurrentEquilibriumSoftTerminationDecision, 1.0f))
				{
					isStateMaintained = false;
				}
				else
				{
					isStateMaintained = true;
				}
			}
			else
			{
				humanState.CurrentEquilibriumSoftTerminationDecision = 0.0f;

				isStateMaintained = true;
			}

			if (!isStateMaintained
				|| !primaryParticleState.ConstrainedState.has_value()
				|| primaryParticleState.ConstrainedState->MeshRelativeVelocity.length() >= MaxRelativeVelocityForEquilibrium
				|| !CheckAndMaintainHumanEquilibrium(
					primaryParticleState.ParticleIndex,
					secondaryParticleState.ParticleIndex,
					humanState.CurrentBehavior == StateType::HumanNpcStateType::BehaviorType::Constrained_Rising,
					isOnEdge,
					mParticles,
					labParameters))
			{
				// Transition to knocked out

				LogMessage("Going to Constrained_KnockedOut; primary's barycentric coords: ", 
					primaryParticleState.ConstrainedState.has_value() ? primaryParticleState.ConstrainedState->CurrentTriangleBarycentricCoords.toString() : "N/A",
					" primary's relative velocity: ", primaryParticleState.ConstrainedState.has_value() ? std::to_string(primaryParticleState.ConstrainedState->MeshRelativeVelocity.length()) : "N/A", 
					" (max=", MaxRelativeVelocityForEquilibrium, ")");

				humanState.TransitionToState(
					StateType::HumanNpcStateType::BehaviorType::Constrained_KnockedOut,
					0.0f,
					0.0f);

				humanState.CurrentWalkMagnitude = 0.0f;

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

				if (isOnEdge) // Note: no need to silence walk as we don't apply walk displacement in inertial (i.e. not-on-edge) case
				{
					// Impart walk displacement & run walking state machine
					RunWalkingHumanStateMachine(
						humanState,
						primaryParticleState,
						mesh,
						labParameters);
				}

				if (humanState.CurrentWalkFlipDecision != 0.0f)
					publishStateQuantity = std::make_tuple("WalkFlip", std::to_string(humanState.CurrentWalkFlipDecision));
				else
					publishStateQuantity = std::make_tuple("EquilibriumTermination", std::to_string(humanState.CurrentEquilibriumSoftTerminationDecision));
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

bool Npcs::CheckAndMaintainHumanEquilibrium(
	ElementIndex primaryParticleIndex,
	ElementIndex secondaryParticleIndex,
	bool isRisingState,
	bool isOnEdge,
	NpcParticles & particles,
	LabParameters const & /*labParameters*/)
{
	//
	// Make sure we are not falling out of equilibrium
	//

	vec2f const humanVector = particles.GetPosition(secondaryParticleIndex) - particles.GetPosition(primaryParticleIndex);

	// Calculate CW angle between head and vertical (pointing up);
	// positive when human is CW wrt vertical
	//
	// |   H
	// |  /
	// |-/
	// |/
	//
	float const staticDisplacementAngleCW = (-LabParameters::GravityDir).angleCw(humanVector);

	// Calculate CW angle that head would rotate by (relative to feet) due to relative velocity alone;
	// positive when new position is CW wrt old
	//
	// |   H
	// |  /
	// | /\
    // |/__L___H'
	//
	vec2f const relativeVelocityDisplacement = (particles.GetVelocity(secondaryParticleIndex) - particles.GetVelocity(primaryParticleIndex)) * LabParameters::SimulationTimeStepDuration;
	float const relativeVelocityAngleCW = humanVector.angleCw(humanVector + relativeVelocityDisplacement);

	//
	// Check whether we are still in equilibrium
	//
	// We lose equilibrium if HumanVector is outside of sector around vertical, with non-negligible rotation velocity towards outside of sector
	//

	float constexpr MaxStaticAngleForEquilibrium = Pi<float> / 7.0f;

	float const maxRelativeVelocityAngleForEqulibrium = isRisingState
		? 0.01f
		: 0.0f;

	if (std::abs(staticDisplacementAngleCW) >= MaxStaticAngleForEquilibrium
		&& std::abs(relativeVelocityAngleCW) >= maxRelativeVelocityAngleForEqulibrium // Large abs angle == velocity towards divergence
		&& staticDisplacementAngleCW * relativeVelocityAngleCW > 0.0f) // Equal signs
	{
		LogMessage("Losing equilibrium because: StaticDisplacementAngleCW=", staticDisplacementAngleCW, " (Max=+/-", MaxStaticAngleForEquilibrium, ") RelativeVelocityAngleCW=", relativeVelocityAngleCW,
			" (Max=+/-", maxRelativeVelocityAngleForEqulibrium, ")");

		return false;
	}

	//
	// We are in equilibrium
	//

	//
	// Maintain equilibrium
	//

	if (isOnEdge)
	{
		particles.SetEquilibriumTorque(secondaryParticleIndex, 1.0f);
	}

	return true;
}

void Npcs::RunWalkingHumanStateMachine(
	StateType::HumanNpcStateType & humanState,
	StateType::NpcParticleStateType const & primaryParticleState,
	Mesh const & /*mesh*/, // Will come useful when we'll *plan* the walk
	LabParameters const & labParameters)
{
	assert(primaryParticleState.ConstrainedState.has_value());

	// 1. Check condition for potentially flipping: actual (relative) velocity opposite of walking direction,
	// or too small

	float constexpr MinRelativeVelocityAgreementToAcceptWalk = 0.025f;
	float const relativeVelocityAgreement = primaryParticleState.ConstrainedState->MeshRelativeVelocity.dot(vec2f(humanState.CurrentFaceDirectionX * humanState.CurrentWalkMagnitude, 0.0f));
	if (relativeVelocityAgreement < MinRelativeVelocityAgreementToAcceptWalk
		&& labParameters.HumanNpcWalkingAcceleration > 0.0f) // Video
	{ 
		// Flip later
		FlipHumanWalk(humanState, StrongTypedFalse<_DoImmediate>);
	}
	else
	{
		// We're doing good, no flipping at the horizon
		humanState.CurrentWalkFlipDecision = 0.0f;
		humanState.TargetWalkFlipDecision = 0.0f;
	}

	// 2. Advance CurrentFlipDecision towards TargetFlipDecision

	float constexpr ToTargetConvergenceRate = 0.1f;

	humanState.CurrentWalkFlipDecision += (humanState.TargetWalkFlipDecision - humanState.CurrentWalkFlipDecision) * ToTargetConvergenceRate;

	// 3. Check if time to flip
	if (humanState.CurrentWalkFlipDecision >= 0.95f)
	{
		// Flip now
		FlipHumanWalk(humanState, StrongTypedTrue<_DoImmediate>);		
	}

	//
	// 4. Advance walking magnitude towards target
	//

	humanState.CurrentWalkMagnitude += (humanState.TargetWalkMagnitude - humanState.CurrentWalkMagnitude) * labParameters.HumanNpcWalkingAcceleration;
	
	LogMessage("        currentWalkMagnitude: ", humanState.CurrentWalkMagnitude);
}

void Npcs::FlipHumanWalk(
	StateType::HumanNpcStateType & humanState,
	DoImmediate doImmediate) const
{
	if (doImmediate.Value)
	{
		humanState.CurrentFaceDirectionX *= -1.0f;
		humanState.CurrentWalkMagnitude = 0.0f;

		LogMessage("Flipping walk: ", humanState.CurrentFaceDirectionX);

		humanState.TargetWalkFlipDecision = 0.0f;
		humanState.CurrentWalkFlipDecision = 0.0f;
	}
	else
	{
		humanState.TargetWalkFlipDecision = 1.0f;
	}
}