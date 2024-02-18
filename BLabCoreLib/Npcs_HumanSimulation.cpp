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
	StateType::NpcParticleStateType const & secondaryParticleState,
	float currentSimulationTime) const
{
	// Return state

	if (!primaryParticleState.ConstrainedState.has_value()
		&& !secondaryParticleState.ConstrainedState.has_value())
	{
		// Whole NPC is free
		mEventDispatcher.OnHumanNpcBehaviorChanged("Free_Aerial");
		return StateType::HumanNpcStateType(StateType::HumanNpcStateType::BehaviorType::Free_Aerial, currentSimulationTime);
	}
	else
	{
		// NPC is constrained
		mEventDispatcher.OnHumanNpcBehaviorChanged("Constrained_KnockedOut");
		return StateType::HumanNpcStateType(StateType::HumanNpcStateType::BehaviorType::Constrained_KnockedOut, currentSimulationTime);
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
	float const ToRisingConvergenceRate = 0.067f;
	float const ToWalkingConvergenceRate = 0.09f;
	float constexpr MaxRelativeVelocityMagnitudeForEquilibrium = 3.0f; // So high because we slip a lot while we try to stand up, and thus need to be immune to ourselves

	std::optional<std::tuple<std::string, std::string>> publishStateQuantity;

	bool const isFree = !primaryParticleState.ConstrainedState.has_value();

	switch (humanState.CurrentBehavior)
	{
		case StateType::HumanNpcStateType::BehaviorType::Constrained_KnockedOut:
		{
			if (isFree)
			{
				// Transition
				humanState.TransitionToState(StateType::HumanNpcStateType::BehaviorType::Free_Aerial, currentSimulationTime);

				mEventDispatcher.OnHumanNpcBehaviorChanged("Free_Aerial");

				break;
			}

			// Check conditions for rising

			float risingTarget = 0.0f;
			if (primaryParticleState.ConstrainedState.has_value()
				&& IsOnFloorEdge(*primaryParticleState.ConstrainedState, mesh)
				&& primaryParticleState.ConstrainedState->MeshRelativeVelocity.length() < MaxRelativeVelocityMagnitudeForEquilibrium
				&& (!secondaryParticleState.ConstrainedState.has_value()
					|| secondaryParticleState.ConstrainedState->MeshRelativeVelocity.length() < MaxRelativeVelocityMagnitudeForEquilibrium))
			{
				risingTarget = 1.0f;
			}


			// Advance towards rising

			humanState.CurrentBehaviorState.Constrained_KnockedOut.ProgressToRising +=
				(risingTarget - humanState.CurrentBehaviorState.Constrained_KnockedOut.ProgressToRising)
				* ToRisingConvergenceRate;

			publishStateQuantity = std::make_tuple("ProgressToRising", std::to_string(humanState.CurrentBehaviorState.Constrained_KnockedOut.ProgressToRising));

			if (IsAtTarget(humanState.CurrentBehaviorState.Constrained_KnockedOut.ProgressToRising, 1.0f))
			{
				// Transition

				humanState.TransitionToState(StateType::HumanNpcStateType::BehaviorType::Constrained_Rising, currentSimulationTime);

				mEventDispatcher.OnHumanNpcBehaviorChanged("Constrained_Rising");

				break;
			}


			break;
		}

		case StateType::HumanNpcStateType::BehaviorType::Constrained_Rising:
		case StateType::HumanNpcStateType::BehaviorType::Constrained_Equilibrium:
		case StateType::HumanNpcStateType::BehaviorType::Constrained_Walking:
		{
			if (isFree)
			{
				// Transition

				humanState.TransitionToState(StateType::HumanNpcStateType::BehaviorType::Free_Aerial, currentSimulationTime);

				mEventDispatcher.OnHumanNpcBehaviorChanged("Free_Aerial");

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

					humanState.CurrentBehaviorState.Constrained_Equilibrium.ProgressToWalking += (1.0f - humanState.CurrentBehaviorState.Constrained_Equilibrium.ProgressToWalking) * ToWalkingConvergenceRate;

					publishStateQuantity = std::make_tuple("ProgressToWalking", std::to_string(humanState.CurrentBehaviorState.Constrained_Equilibrium.ProgressToWalking));

					if (IsAtTarget(humanState.CurrentBehaviorState.Constrained_Equilibrium.ProgressToWalking, 1.0f))
					{
						// Transition

						humanState.TransitionToState(StateType::HumanNpcStateType::BehaviorType::Constrained_Walking, currentSimulationTime);

						// Keep torque

						mEventDispatcher.OnHumanNpcBehaviorChanged("Constrained_Walking");

						break;
					}
				}
			}

			//
			// Check conditions to stay & maintain equilibrium
			//

			// a. On-edge

			bool isStateMaintained;
			if (!isOnEdge)
			{
				// When walking, we want to be a bit more tolerant about "losing the edge"

				// This is a quite important parameter: it's the duration through which we tolerate temporarily losing contact
				// with the ground
				float const toTerminateEquilibriumConvergenceRate = (humanState.CurrentBehavior != StateType::HumanNpcStateType::BehaviorType::Constrained_Walking)
					? 0.25f
					: 0.1f;

				// Advance
				humanState.CurrentEquilibriumSoftTerminationDecision += (1.0f - humanState.CurrentEquilibriumSoftTerminationDecision) * toTerminateEquilibriumConvergenceRate;

				// Check if enough
				if (IsAtTarget(humanState.CurrentEquilibriumSoftTerminationDecision, 1.0f))
				{
					LogMessage("Been off-edge for too long");

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

			// b. Mesh-relative velocity

			if (humanState.CurrentBehavior != StateType::HumanNpcStateType::BehaviorType::Constrained_Walking)
			{
				// Not walking: we want to be draconian and can't stand a (small) relative velocity
				if (primaryParticleState.ConstrainedState->MeshRelativeVelocity.length() >= MaxRelativeVelocityMagnitudeForEquilibrium)
				{
					isStateMaintained = false;
				}
			}
			else
			{
				// We need to take into account that we _are_ moving because we want it

				assert(humanState.CurrentFaceDirectionX != 0.0f);

				vec2f const idealWalkVelocityDir = vec2f(
					humanState.CurrentFaceDirectionX,
					0.0f);
				float const idealWalkVelocityMagnitude = CalculateActualHumanWalkingAbsoluteSpeed(humanState, labParameters);
				vec2f const idealWalkVelocity = idealWalkVelocityDir * idealWalkVelocityMagnitude;

				float const primaryMeshRelativeVelocityAlongWalkDir = primaryParticleState.ConstrainedState->MeshRelativeVelocity.dot(idealWalkVelocityDir);

				LogMessage("idealWalkVelocity=", idealWalkVelocity, " (mag=", idealWalkVelocityMagnitude, ") ",
					"meshRelativeVelocity=", primaryParticleState.ConstrainedState->MeshRelativeVelocity, " (along dir =", primaryMeshRelativeVelocityAlongWalkDir, ")");

				if (primaryMeshRelativeVelocityAlongWalkDir >= 0.0f)
				{
					// Same direction as walking

					// Stop if it's too much over
					float constexpr MaxAlignedRelativeVelocityMagnitudeForWalking = 5.0f;
					if (primaryMeshRelativeVelocityAlongWalkDir >= MaxAlignedRelativeVelocityMagnitudeForWalking)
					{
						LogMessage("MRV too much in same direction");
						isStateMaintained = false;
					}
				}
				else
				{
					// Opposite direction to walking

					// Note: this check is also at flipping decision in maintaining walking state machine

					// TODOTEST
					////// Stop if we're at least a little bit opposite
					////float constexpr MaxOppositeRelativeVelocityMagnitudeForWalking = 1.0f;
					////if (primaryMeshRelativeVelocityAlongWalkDir <= -MaxOppositeRelativeVelocityMagnitudeForWalking)
					////{
					////	LogMessage("MRV too much opposite");
					////	isStateMaintained = false;
					////}
				}
			}

			// . Check

			assert(primaryParticleState.ConstrainedState.has_value() || !isStateMaintained);

			if (!isStateMaintained
				|| !CheckAndMaintainHumanEquilibrium(
					primaryParticleState.ParticleIndex,
					secondaryParticleState.ParticleIndex,
					humanState.CurrentBehavior == StateType::HumanNpcStateType::BehaviorType::Constrained_Rising,
					isOnEdge, // doMaintainEquilibrium
					mParticles,
					labParameters))
			{
				// Transition to knocked out

				LogMessage("Going to Constrained_KnockedOut; primary's barycentric coords: ",
					primaryParticleState.ConstrainedState.has_value() ? primaryParticleState.ConstrainedState->CurrentTriangleBarycentricCoords.toString() : "N/A",
					" primary's relative velocity mag: ", primaryParticleState.ConstrainedState.has_value() ? std::to_string(primaryParticleState.ConstrainedState->MeshRelativeVelocity.length()) : "N/A",
					" (max=", MaxRelativeVelocityMagnitudeForEquilibrium, ")");

				humanState.TransitionToState(StateType::HumanNpcStateType::BehaviorType::Constrained_KnockedOut, currentSimulationTime);

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

					humanState.TransitionToState(StateType::HumanNpcStateType::BehaviorType::Constrained_Equilibrium, currentSimulationTime);

					mEventDispatcher.OnHumanNpcBehaviorChanged("Constrained_Equilibrium");
				}
			}
			else if (humanState.CurrentBehavior == StateType::HumanNpcStateType::BehaviorType::Constrained_Equilibrium)
			{
				// Nop
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

				if (humanState.CurrentBehaviorState.Constrained_Walking.CurrentFlipDecision != 0.0f)
					publishStateQuantity = std::make_tuple("WalkFlip", std::to_string(humanState.CurrentBehaviorState.Constrained_Walking.CurrentFlipDecision));
				else
					publishStateQuantity = std::make_tuple("EquilibriumTermination", std::to_string(humanState.CurrentEquilibriumSoftTerminationDecision));
			}

			break;
		}

		case StateType::HumanNpcStateType::BehaviorType::Free_Aerial:
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
	bool doMaintainEquilibrium,
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

	if (doMaintainEquilibrium)
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
	assert(humanState.CurrentBehavior == StateType::HumanNpcStateType::BehaviorType::Constrained_Walking);

	auto & walkingState = humanState.CurrentBehaviorState.Constrained_Walking;

	// 1. Check condition for potentially flipping: actual (relative) velocity opposite of walking direction,
	// or too small

	if (walkingState.CurrentWalkMagnitude != 0.0f)
	{
		float constexpr MinRelativeVelocityAgreementToAcceptWalk = 0.025f;
		float const relativeVelocityAgreement = primaryParticleState.ConstrainedState->MeshRelativeVelocity.dot(
			vec2f(
				humanState.CurrentFaceDirectionX * CalculateActualHumanWalkingAbsoluteSpeed(humanState, labParameters),
				0.0f));
		if (relativeVelocityAgreement < MinRelativeVelocityAgreementToAcceptWalk)
		{
			// Flip later
			FlipHumanWalk(humanState, StrongTypedFalse<_DoImmediate>);
		}
		else
		{
			// We're doing good, no flipping at the horizon
			walkingState.CurrentFlipDecision = 0.0f;
			walkingState.TargetFlipDecision = 0.0f;
		}
	}

	// 2. Advance CurrentFlipDecision towards TargetFlipDecision

	float constexpr ToTargetConvergenceRate = 0.1f;

	walkingState.CurrentFlipDecision += (walkingState.TargetFlipDecision - walkingState.CurrentFlipDecision) * ToTargetConvergenceRate;

	// 3. Check if time to flip
	if (walkingState.CurrentFlipDecision >= 0.95f)
	{
		// Flip now
		FlipHumanWalk(humanState, StrongTypedTrue<_DoImmediate>);
	}

	//
	// 4. Advance walking magnitude towards full walk
	//

	float constexpr WalkMagnitudeConvergenceRate = 0.03f;
	walkingState.CurrentWalkMagnitude += (1.0f - walkingState.CurrentWalkMagnitude) * WalkMagnitudeConvergenceRate;

	LogMessage("        currentWalkMagnitude: ", walkingState.CurrentWalkMagnitude);
}

void Npcs::OnHumanImpact(
	vec2f const & /*impactVector*/,
	vec2f const & bounceEdgeNormal,
	StateType & npc,
	bool /*isPrimaryParticle*/) const
{
	switch (npc.HumanNpcState->CurrentBehavior)
	{
		case StateType::HumanNpcStateType::BehaviorType::Constrained_Walking:
		{
			LogMessage("OnHumanImpact: alignment=", bounceEdgeNormal.dot(vec2f(npc.HumanNpcState->CurrentFaceDirectionX, 0.0f)));

			// Check alignment of impact with walking direction; if hit => flip
			float constexpr MaxOppositionSlope = 0.5f;
			if (bounceEdgeNormal.dot(vec2f(npc.HumanNpcState->CurrentFaceDirectionX, 0.0f)) > MaxOppositionSlope
				&& npc.HumanNpcState->CurrentBehaviorState.Constrained_Walking.CurrentWalkMagnitude != 0.0f)
			{
				LogMessage("OnHumanImpact: FLIP!");

				// Flip now
				FlipHumanWalk(*npc.HumanNpcState, StrongTypedTrue<_DoImmediate>);
			}

			break;
		}

		default:
		{
			break;
		}
	}
}

void Npcs::FlipHumanWalk(
	StateType::HumanNpcStateType & humanState,
	DoImmediate doImmediate) const
{
	assert(humanState.CurrentBehavior == StateType::HumanNpcStateType::BehaviorType::Constrained_Walking);
	auto & walkingState = humanState.CurrentBehaviorState.Constrained_Walking;

	if (doImmediate.Value)
	{
		humanState.CurrentFaceDirectionX *= -1.0f;
		walkingState.CurrentWalkMagnitude = 0.0f;

		LogMessage("Flipping walk: ", humanState.CurrentFaceDirectionX);

		walkingState.TargetFlipDecision = 0.0f;
		walkingState.CurrentFlipDecision = 0.0f;
	}
	else
	{
		walkingState.TargetFlipDecision = 1.0f;
	}
}

float Npcs::CalculateActualHumanWalkingAbsoluteSpeed(
	StateType::HumanNpcStateType & humanState,
	LabParameters const & labParameters) const
{
	assert(humanState.CurrentBehavior == StateType::HumanNpcStateType::BehaviorType::Constrained_Walking);

	return std::min(
		labParameters.HumanNpcWalkingSpeed
		* humanState.CurrentBehaviorState.Constrained_Walking.CurrentWalkMagnitude // Note that this is the only one that might be zero
		* (1.0f + humanState.PanicLevel),
		4.0f); // Absolute cap
}
