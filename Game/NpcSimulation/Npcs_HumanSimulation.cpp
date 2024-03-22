/***************************************************************************************
 * Original Author:		Gabriele Giuseppini
 * Created:				2023-10-20
 * Copyright:			Gabriele Giuseppini  (https://github.com/GabrieleGiuseppini)
 ***************************************************************************************/
#include "Physics.h"

#include <GameCore/GameMath.h>
#include <GameCore/GameRandomEngine.h>

#include <cmath>

namespace Physics {

namespace {

	bool IsAtTarget(float currentValue, float targetValue)
	{
		return std::abs(targetValue - currentValue) < 0.01f;
	}

}

Npcs::StateType::KindSpecificStateType::HumanNpcStateType Npcs::CalculateInitialHumanState(
	StateType & npc,
	float currentSimulationTime)
{
	assert(npc.DipoleState.has_value());
	auto const & primaryParticleState = npc.PrimaryParticleState;
	auto const & secondaryParticleState = npc.DipoleState->SecondaryParticleState;

	if (!primaryParticleState.ConstrainedState.has_value()
		&& !secondaryParticleState.ConstrainedState.has_value())
	{
		// Whole NPC is free
		return StateType::KindSpecificStateType::HumanNpcStateType(
			npc.KindSpecificState.HumanNpcState.Kind,
			StateType::KindSpecificStateType::HumanNpcStateType::BehaviorType::Free_Aerial,
			currentSimulationTime);
	}
	else
	{
		// NPC is constrained
		return StateType::KindSpecificStateType::HumanNpcStateType(
			npc.KindSpecificState.HumanNpcState.Kind,
			StateType::KindSpecificStateType::HumanNpcStateType::BehaviorType::Constrained_Aerial,
			currentSimulationTime);
	}
}

void Npcs::UpdateHuman(
	StateType & npc,
	float currentSimulationTime,
	Ship const & shipMesh,
	GameParameters const & gameParameters)
{
	assert(npc.DipoleState.has_value());
	auto const & primaryParticleState = npc.PrimaryParticleState;
	auto const & secondaryParticleState = npc.DipoleState->SecondaryParticleState;

	assert(npc.Kind == NpcKindType::Human);
	auto & humanState = npc.KindSpecificState.HumanNpcState;
	using HumanNpcStateType = StateType::KindSpecificStateType::HumanNpcStateType;

	//
	// Constants
	//

	float constexpr MaxRelativeVelocityMagnitudeForEquilibrium = 3.0f; // So high because we slip a lot while we try to stand up, and thus need to be immune to ourselves

	//
	// Update panic
	//

	humanState.ResultantPanicLevel =
		humanState.ShipOnFirePanicLevel
		+ humanState.GeneralizedPanicLevel;

	// Decay

	humanState.ShipOnFirePanicLevel -= humanState.ShipOnFirePanicLevel * 0.01f;

	//
	// Process human
	//

	std::optional<std::tuple<std::string, std::string>> publishStateQuantity;

	bool const isFree = !primaryParticleState.ConstrainedState.has_value();

	switch (humanState.CurrentBehavior)
	{
		case HumanNpcStateType::BehaviorType::BeingPlaced:
		{
			// Nop
			break;
		}

		case HumanNpcStateType::BehaviorType::Constrained_Aerial:
		{
			if (isFree)
			{
				// Transition
				TransitionHumanToFree(npc, currentSimulationTime);

				break;
			}

			// Check conditions for falling/KO

			bool const isHeadOnFloor = secondaryParticleState.ConstrainedState.has_value() && IsOnFloorEdge(*secondaryParticleState.ConstrainedState, shipMesh);
			bool const areFeetOnFloor = primaryParticleState.ConstrainedState.has_value() && IsOnFloorEdge(*primaryParticleState.ConstrainedState, shipMesh);

			vec2f const floorVector = (primaryParticleState.ConstrainedState.has_value() && primaryParticleState.ConstrainedState->CurrentVirtualEdgeOrdinal >= 0)
				? shipMesh.GetTriangles().GetSubSpringVector(primaryParticleState.ConstrainedState->CurrentTriangle, primaryParticleState.ConstrainedState->CurrentVirtualEdgeOrdinal, shipMesh.GetPoints())
				: vec2f(1.0f, 0.0); // H arbitrarily
			float const headVelocityAlongFloor = secondaryParticleState.GetApplicableVelocity(mParticles).dot(floorVector);
			float const feetVelocityAlongFloor = primaryParticleState.GetApplicableVelocity(mParticles).dot(floorVector);

			float constexpr MinVelocityMagnitudeForFalling = 0.05f;

			float fallingTarget;
			float knockedOutTarget;
			if (isHeadOnFloor || areFeetOnFloor)
			{
				if (std::abs(headVelocityAlongFloor) >= MinVelocityMagnitudeForFalling
					|| std::abs(feetVelocityAlongFloor) >= MinVelocityMagnitudeForFalling)
				{
					// Likely falling
					fallingTarget = 1.0f;

					// Definitely not knocked out
					knockedOutTarget = 0.0f;
				}
				else
				{
					// We're quite still...

					// ...likely knocked out
					knockedOutTarget = 1.0f;

					// Definitely not falling
					fallingTarget = 0.0f;
				}
			}
			else
			{
				// Completely in the air, no transition
				fallingTarget = 0.0f;
				knockedOutTarget = 0.0f;
			}

			// Progress to falling

			float constexpr ToFallingConvergenceRate = 0.75f; // Very high! We do this just to survive micro-instants

			humanState.CurrentBehaviorState.Constrained_Aerial.ProgressToFalling +=
				(fallingTarget - humanState.CurrentBehaviorState.Constrained_Aerial.ProgressToFalling)
				* ToFallingConvergenceRate;

			publishStateQuantity = std::make_tuple("ProgressToFalling", std::to_string(humanState.CurrentBehaviorState.Constrained_Aerial.ProgressToFalling));

			if (IsAtTarget(humanState.CurrentBehaviorState.Constrained_Aerial.ProgressToFalling, 1.0f))
			{
				// Transition

				humanState.TransitionToState(HumanNpcStateType::BehaviorType::Constrained_Falling, currentSimulationTime);

				if (humanState.CurrentFaceOrientation != 0.0f)
				{
					// Face: 0/direction of falling
					humanState.CurrentFaceOrientation = 0.0f;
					humanState.CurrentFaceDirectionX =
						(!secondaryParticleState.ConstrainedState.has_value() && mParticles.GetVelocity(secondaryParticleState.ParticleIndex).x >= 0.0f)
						|| (secondaryParticleState.ConstrainedState.has_value() && secondaryParticleState.ConstrainedState->MeshRelativeVelocity.x >= 0.0f)
						? 1.0f
						: -1.0f;
				}

				if (npc.Id == mCurrentlySelectedNpc)
				{
					mGameEventHandler->OnHumanNpcBehaviorChanged("Constrained_Falling");
				}

				break;
			}

			// Progress to knocked out

			float constexpr ToKnockedOutConvergenceRate = 0.5f;

			humanState.CurrentBehaviorState.Constrained_Aerial.ProgressToKnockedOut +=
				(knockedOutTarget - humanState.CurrentBehaviorState.Constrained_Aerial.ProgressToKnockedOut)
				* ToKnockedOutConvergenceRate;

			if (fallingTarget == 0.0f)
				publishStateQuantity = std::make_tuple("ProgressToKnockedOut", std::to_string(humanState.CurrentBehaviorState.Constrained_Aerial.ProgressToKnockedOut));

			if (IsAtTarget(humanState.CurrentBehaviorState.Constrained_Aerial.ProgressToKnockedOut, 1.0f))
			{
				// Transition

				humanState.TransitionToState(HumanNpcStateType::BehaviorType::Constrained_KnockedOut, currentSimulationTime);

				if (npc.Id == mCurrentlySelectedNpc)
				{
					mGameEventHandler->OnHumanNpcBehaviorChanged("Constrained_KnockedOut");
				}

				break;
			}

			break;
		}

		case HumanNpcStateType::BehaviorType::Constrained_Falling:
		{
			if (isFree)
			{
				// Transition
				TransitionHumanToFree(npc, currentSimulationTime);

				break;
			}

			// Check conditions for knocked out

			bool const areFootOnFloor = primaryParticleState.ConstrainedState.has_value() && IsOnFloorEdge(*primaryParticleState.ConstrainedState, shipMesh);

			float knockedOutTarget = 0.0f;
			float constexpr MaxRelativeVelocityForKnockedOut = 1.2f;
			if (areFootOnFloor
				&& primaryParticleState.ConstrainedState->MeshRelativeVelocity.length() < MaxRelativeVelocityForKnockedOut
				&& secondaryParticleState.GetApplicableVelocity(mParticles).length() < MaxRelativeVelocityForKnockedOut)
			{
				knockedOutTarget = 1.0f;
			}

			// Advance towards knocked out

			float const toKnockedConvergenceRate = 0.14f + std::min(humanState.ResultantPanicLevel, 1.0f) * 0.07f;
			humanState.CurrentBehaviorState.Constrained_Falling.ProgressToKnockedOut +=
				(knockedOutTarget - humanState.CurrentBehaviorState.Constrained_Falling.ProgressToKnockedOut)
				* toKnockedConvergenceRate;

			if (knockedOutTarget > 0.0f)
				publishStateQuantity = std::make_tuple("ProgressToKnockedOut", std::to_string(humanState.CurrentBehaviorState.Constrained_Falling.ProgressToKnockedOut));

			if (IsAtTarget(humanState.CurrentBehaviorState.Constrained_Falling.ProgressToKnockedOut, 1.0f))
			{
				// Transition

				humanState.TransitionToState(HumanNpcStateType::BehaviorType::Constrained_KnockedOut, currentSimulationTime);

				if (npc.Id == mCurrentlySelectedNpc)
				{
					mGameEventHandler->OnHumanNpcBehaviorChanged("Constrained_KnockedOut");
				}

				break;
			}

			// Check conditions for aerial

			bool const isHeadOnFloor = secondaryParticleState.ConstrainedState.has_value() && IsOnFloorEdge(*secondaryParticleState.ConstrainedState, shipMesh);

			if (!areFootOnFloor && !isHeadOnFloor)
			{
				// Advance towards aerial

				float constexpr ToAerialConvergenceRate = 0.35f;

				humanState.CurrentBehaviorState.Constrained_Falling.ProgressToAerial +=
					(1.0f - humanState.CurrentBehaviorState.Constrained_Falling.ProgressToAerial)
					* ToAerialConvergenceRate;

				if (knockedOutTarget == 0.0f)
					publishStateQuantity = std::make_tuple("ProgressToAerial", std::to_string(humanState.CurrentBehaviorState.Constrained_Falling.ProgressToAerial));

				if (IsAtTarget(humanState.CurrentBehaviorState.Constrained_Falling.ProgressToAerial, 1.0f))
				{
					// Transition

					humanState.TransitionToState(HumanNpcStateType::BehaviorType::Constrained_Aerial, currentSimulationTime);

					if (npc.Id == mCurrentlySelectedNpc)
					{
						mGameEventHandler->OnHumanNpcBehaviorChanged("Constrained_Aerial");
					}

					break;
				}
			}
			else
			{
				// Reset progress to aerial
				humanState.CurrentBehaviorState.Constrained_Falling.ProgressToAerial = 0.0f;
				if (knockedOutTarget == 0.0f)
					publishStateQuantity = std::make_tuple("ProgressToAerial", std::to_string(humanState.CurrentBehaviorState.Constrained_Falling.ProgressToAerial));
			}

			break;
		}

		case HumanNpcStateType::BehaviorType::Constrained_KnockedOut:
		{
			if (isFree)
			{
				// Transition
				TransitionHumanToFree(npc, currentSimulationTime);

				break;
			}

			// Check conditions for rising

			bool const areFootOnFloor = primaryParticleState.ConstrainedState.has_value() && IsOnFloorEdge(*primaryParticleState.ConstrainedState, shipMesh);

			float risingTarget = 0.0f;
			if (areFootOnFloor
				&& primaryParticleState.ConstrainedState->MeshRelativeVelocity.length() < MaxRelativeVelocityMagnitudeForEquilibrium
				&& secondaryParticleState.GetApplicableVelocity(mParticles).length() < MaxRelativeVelocityMagnitudeForEquilibrium)
			{
				risingTarget = 1.0f;
			}

			// Advance towards rising

			float const toRisingConvergenceRate = 0.067f + std::min(humanState.ResultantPanicLevel, 1.0f) * 0.07f;
			humanState.CurrentBehaviorState.Constrained_KnockedOut.ProgressToRising +=
				(risingTarget - humanState.CurrentBehaviorState.Constrained_KnockedOut.ProgressToRising)
				* toRisingConvergenceRate;

			publishStateQuantity = std::make_tuple("ProgressToRising", std::to_string(humanState.CurrentBehaviorState.Constrained_KnockedOut.ProgressToRising));

			if (IsAtTarget(humanState.CurrentBehaviorState.Constrained_KnockedOut.ProgressToRising, 1.0f))
			{
				// Transition

				humanState.TransitionToState(HumanNpcStateType::BehaviorType::Constrained_Rising, currentSimulationTime);

				if (npc.Id == mCurrentlySelectedNpc)
				{
					mGameEventHandler->OnHumanNpcBehaviorChanged("Constrained_Rising");
				}

				break;
			}

			// Check conditions for aerial

			bool const isHeadOnFloor = secondaryParticleState.ConstrainedState.has_value() && IsOnFloorEdge(*secondaryParticleState.ConstrainedState, shipMesh);

			if (!areFootOnFloor && !isHeadOnFloor)
			{
				// Advance towards aerial

				float constexpr ToAerialConvergenceRate = 0.15f;

				humanState.CurrentBehaviorState.Constrained_KnockedOut.ProgressToAerial +=
					(1.0f - humanState.CurrentBehaviorState.Constrained_KnockedOut.ProgressToAerial)
					* ToAerialConvergenceRate;

				publishStateQuantity = std::make_tuple("ProgressToAerial", std::to_string(humanState.CurrentBehaviorState.Constrained_KnockedOut.ProgressToAerial));

				if (IsAtTarget(humanState.CurrentBehaviorState.Constrained_KnockedOut.ProgressToAerial, 1.0f))
				{
					// Transition

					humanState.TransitionToState(HumanNpcStateType::BehaviorType::Constrained_Aerial, currentSimulationTime);

					if (npc.Id == mCurrentlySelectedNpc)
					{
						mGameEventHandler->OnHumanNpcBehaviorChanged("Constrained_Aerial");
					}

					break;
				}
			}
			else
			{
				// Reset progress to aerial
				humanState.CurrentBehaviorState.Constrained_KnockedOut.ProgressToAerial = 0.0f;

				if (risingTarget == 0.0f)
					publishStateQuantity = std::make_tuple("ProgressToAerial", std::to_string(humanState.CurrentBehaviorState.Constrained_KnockedOut.ProgressToAerial));
			}



			break;
		}

		case HumanNpcStateType::BehaviorType::Constrained_Rising:
		case HumanNpcStateType::BehaviorType::Constrained_Equilibrium:
		case HumanNpcStateType::BehaviorType::Constrained_Walking:
		{
			if (isFree)
			{
				// Transition
				TransitionHumanToFree(npc, currentSimulationTime);

				break;
			}

			bool const isOnEdge =
				primaryParticleState.ConstrainedState.has_value()
				&& IsOnFloorEdge(*primaryParticleState.ConstrainedState, shipMesh);

			if (humanState.CurrentBehavior == HumanNpcStateType::BehaviorType::Constrained_Equilibrium)
			{
				if (isOnEdge)
				{
					// Advance towards walking

					float const toWalkingConvergenceRate = 0.09f + std::min(humanState.ResultantPanicLevel, 1.0f) * 0.15f;
					humanState.CurrentBehaviorState.Constrained_Equilibrium.ProgressToWalking +=
						(1.0f - humanState.CurrentBehaviorState.Constrained_Equilibrium.ProgressToWalking)
						* toWalkingConvergenceRate;

					publishStateQuantity = std::make_tuple("ProgressToWalking", std::to_string(humanState.CurrentBehaviorState.Constrained_Equilibrium.ProgressToWalking));

					if (IsAtTarget(humanState.CurrentBehaviorState.Constrained_Equilibrium.ProgressToWalking, 1.0f))
					{
						// Transition

						humanState.TransitionToState(HumanNpcStateType::BehaviorType::Constrained_Walking, currentSimulationTime);

						// Face: 0/rnd
						humanState.CurrentFaceOrientation = 0.0f;
						humanState.CurrentFaceDirectionX = GameRandomEngine::GetInstance().GenerateUniformBoolean(0.5f) ? +1.0f : -1.0f;

						// Keep torque

						if (npc.Id == mCurrentlySelectedNpc)
						{
							mGameEventHandler->OnHumanNpcBehaviorChanged("Constrained_Walking");
						}

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
				float const toTerminateEquilibriumConvergenceRate = (humanState.CurrentBehavior != HumanNpcStateType::BehaviorType::Constrained_Walking)
					? 0.25f
					: 0.1f;

				// Advance
				humanState.CurrentEquilibriumSoftTerminationDecision += (1.0f - humanState.CurrentEquilibriumSoftTerminationDecision) * toTerminateEquilibriumConvergenceRate;

				// Check if enough
				if (IsAtTarget(humanState.CurrentEquilibriumSoftTerminationDecision, 1.0f))
				{
					LogNpcDebug("Been off-edge for too long");

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

			if (humanState.CurrentBehavior != HumanNpcStateType::BehaviorType::Constrained_Walking)
			{
				// Not walking: we want to be draconian and can't stand a (small) relative velocity
				if (!primaryParticleState.ConstrainedState.has_value()
					|| primaryParticleState.ConstrainedState->MeshRelativeVelocity.length() >= MaxRelativeVelocityMagnitudeForEquilibrium)
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
				float const idealWalkVelocityMagnitude = CalculateActualHumanWalkingAbsoluteSpeed(humanState, gameParameters);
				vec2f const idealWalkVelocity = idealWalkVelocityDir * idealWalkVelocityMagnitude;

				float const primaryMeshRelativeVelocityAlongWalkDir = primaryParticleState.ConstrainedState->MeshRelativeVelocity.dot(idealWalkVelocityDir);

				LogNpcDebug("idealWalkVelocity=", idealWalkVelocity, " (mag=", idealWalkVelocityMagnitude, ") ",
					"meshRelativeVelocity=", primaryParticleState.ConstrainedState->MeshRelativeVelocity, " (along dir =", primaryMeshRelativeVelocityAlongWalkDir, ")");

				if (primaryMeshRelativeVelocityAlongWalkDir >= 0.0f)
				{
					// Same direction as walking

					// Stop if it's too much over
					float constexpr MaxAlignedRelativeVelocityMagnitudeForWalking = 5.0f;
					if (primaryMeshRelativeVelocityAlongWalkDir >= MaxAlignedRelativeVelocityMagnitudeForWalking)
					{
						LogNpcDebug("MRV too much in same direction");
						isStateMaintained = false;
					}
				}
				else
				{
					// Opposite direction to walking

					// Note: this we check at flipping decision in maintaining walking state machine,
					// so continue
				}
			}

			// . Check

			assert(primaryParticleState.ConstrainedState.has_value() || !isStateMaintained);

			if (!isStateMaintained
				|| !CheckAndMaintainHumanEquilibrium(
					primaryParticleState.ParticleIndex,
					secondaryParticleState.ParticleIndex,
					humanState.CurrentBehavior == HumanNpcStateType::BehaviorType::Constrained_Rising,
					isOnEdge, // doMaintainEquilibrium
					mParticles,
					gameParameters))
			{
				// Transition to aerial/falling, depending on whether we're on an edge

				LogNpcDebug("Going to Constrained_X; primary's barycentric coords: ",
					primaryParticleState.ConstrainedState.has_value() ? primaryParticleState.ConstrainedState->CurrentTriangleBarycentricCoords.toString() : "N/A",
					" primary's relative velocity mag: ", primaryParticleState.ConstrainedState.has_value() ? std::to_string(primaryParticleState.ConstrainedState->MeshRelativeVelocity.length()) : "N/A",
					" (max=", MaxRelativeVelocityMagnitudeForEquilibrium, ")");

				bool const areFootOnEdge = primaryParticleState.ConstrainedState.has_value() && IsOnFloorEdge(*primaryParticleState.ConstrainedState, shipMesh);
				if (areFootOnEdge)
				{
					// Falling

					humanState.TransitionToState(HumanNpcStateType::BehaviorType::Constrained_Falling, currentSimulationTime);

					// Face: 0/direction of falling
					humanState.CurrentFaceOrientation = 0.0f;
					humanState.CurrentFaceDirectionX =
						(!secondaryParticleState.ConstrainedState.has_value() && mParticles.GetVelocity(secondaryParticleState.ParticleIndex).x >= 0.0f)
						|| (secondaryParticleState.ConstrainedState.has_value() && secondaryParticleState.ConstrainedState->MeshRelativeVelocity.x >= 0.0f)
						? 1.0f
						: -1.0f;

					if (npc.Id == mCurrentlySelectedNpc)
					{
						mGameEventHandler->OnHumanNpcBehaviorChanged("Constrained_Falling");
					}
				}
				else
				{
					// Aerial

					humanState.TransitionToState(HumanNpcStateType::BehaviorType::Constrained_Aerial, currentSimulationTime);

					if (npc.Id == mCurrentlySelectedNpc)
					{
						mGameEventHandler->OnHumanNpcBehaviorChanged("Constrained_Aerial");
					}
				}

				break;
			}

			//
			// Update state now
			//

			if (humanState.CurrentBehavior == HumanNpcStateType::BehaviorType::Constrained_Rising)
			{
				// Check if reached alignment (note: here so that we may keep torque as we'll be transitioning to Equilibrium)

				float const alignment = CalculateDipoleVerticalAlignment(primaryParticleState.ParticleIndex, secondaryParticleState.ParticleIndex, mParticles);

				publishStateQuantity = std::make_tuple("Alignment", std::to_string(alignment));

				if (AreAlmostEqual(alignment, 1.0f, 0.004f))
				{
					// Transition

					humanState.TransitionToState(HumanNpcStateType::BehaviorType::Constrained_Equilibrium, currentSimulationTime);

					if (npc.Id == mCurrentlySelectedNpc)
					{
						mGameEventHandler->OnHumanNpcBehaviorChanged("Constrained_Equilibrium");
					}
				}
			}
			else if (humanState.CurrentBehavior == HumanNpcStateType::BehaviorType::Constrained_Equilibrium)
			{
				// Nop
			}
			else
			{
				assert(humanState.CurrentBehavior == HumanNpcStateType::BehaviorType::Constrained_Walking);

				if (isOnEdge) // Note: no need to silence walk as we don't apply walk displacement in inertial (i.e. not-on-edge) case
				{
					// Impart walk displacement & run walking state machine
					RunWalkingHumanStateMachine(
						humanState,
						primaryParticleState,
						shipMesh,
						gameParameters);
				}

				if (humanState.CurrentBehaviorState.Constrained_Walking.CurrentFlipDecision != 0.0f)
					publishStateQuantity = std::make_tuple("WalkFlip", std::to_string(humanState.CurrentBehaviorState.Constrained_Walking.CurrentFlipDecision));
				else
					publishStateQuantity = std::make_tuple("EquilibriumTermination", std::to_string(humanState.CurrentEquilibriumSoftTerminationDecision));
			}

			break;
		}

		case HumanNpcStateType::BehaviorType::Free_Aerial:
		{
			if (!isFree)
			{
				// Transition

				humanState.TransitionToState(HumanNpcStateType::BehaviorType::Constrained_KnockedOut, currentSimulationTime);

				if (npc.Id == mCurrentlySelectedNpc)
				{
					mGameEventHandler->OnHumanNpcBehaviorChanged("Constrained_KnockedOut");
				}

				break;
			}

			// Check if moved to water

			assert(npc.DipoleState.has_value());
			auto const & headPosition = mParticles.GetPosition(npc.DipoleState->SecondaryParticleState.ParticleIndex);
			auto const & feetPosition = mParticles.GetPosition(npc.PrimaryParticleState.ParticleIndex);

			// It's in water if both in water
			if (mParentWorld.GetOceanSurface().GetDepth(headPosition) > 0.0f
				&& mParentWorld.GetOceanSurface().GetDepth(feetPosition) > 0.0f)
			{
				// Transition

				humanState.TransitionToState(HumanNpcStateType::BehaviorType::Free_InWater, currentSimulationTime);

				if (npc.Id == mCurrentlySelectedNpc)
				{
					mGameEventHandler->OnHumanNpcBehaviorChanged("Free_InWater");
				}

				break;
			}

			break;
		}

		case HumanNpcStateType::BehaviorType::Free_InWater:
		case HumanNpcStateType::BehaviorType::Free_Swimming:
		{
			if (!isFree)
			{
				// Transition

				humanState.TransitionToState(HumanNpcStateType::BehaviorType::Constrained_KnockedOut, currentSimulationTime);

				if (npc.Id == mCurrentlySelectedNpc)
				{
					mGameEventHandler->OnHumanNpcBehaviorChanged("Constrained_KnockedOut");
				}

				break;
			}

			// Check if moved to air

			assert(npc.DipoleState.has_value());
			auto const & headPosition = mParticles.GetPosition(npc.DipoleState->SecondaryParticleState.ParticleIndex);
			auto const & feetPosition = mParticles.GetPosition(npc.PrimaryParticleState.ParticleIndex);

			// It's in air if both in air
			if (mParentWorld.GetOceanSurface().GetDepth(headPosition) <= 0.0f
				&& mParentWorld.GetOceanSurface().GetDepth(feetPosition) <= 0.0f)
			{
				// Transition

				humanState.TransitionToState(HumanNpcStateType::BehaviorType::Free_Aerial, currentSimulationTime);

				if (npc.Id == mCurrentlySelectedNpc)
				{
					mGameEventHandler->OnHumanNpcBehaviorChanged("Free_Aerial");
				}

				break;
			}

			// Advance state machine

			if (humanState.CurrentBehavior == HumanNpcStateType::BehaviorType::Free_InWater)
			{
				// Progress to swimming if not rotating and head above feet

				float const rotationMagnitude = (mParticles.GetVelocity(npc.DipoleState->SecondaryParticleState.ParticleIndex) - mParticles.GetVelocity(npc.PrimaryParticleState.ParticleIndex)).length();
				float const targetSwim =
					(1.0f - Step(2.0f, rotationMagnitude))
					* Step(feetPosition.y, headPosition.y);

				float constexpr ToSwimmingConvergenceRate = 0.12f;
				humanState.CurrentBehaviorState.Free_InWater.ProgressToSwimming += (targetSwim - humanState.CurrentBehaviorState.Free_InWater.ProgressToSwimming) * ToSwimmingConvergenceRate;

				publishStateQuantity = std::make_tuple("ProgressToSwimming", std::to_string(humanState.CurrentBehaviorState.Free_InWater.ProgressToSwimming));

				if (IsAtTarget(humanState.CurrentBehaviorState.Free_InWater.ProgressToSwimming, 0.9f)) // We're content with "almost"
				{
					// Transition

					humanState.TransitionToState(HumanNpcStateType::BehaviorType::Free_Swimming, currentSimulationTime);

					// Face: FvB/0
					humanState.CurrentFaceOrientation = 1.0f; // TODO: random: back
					humanState.CurrentFaceDirectionX = 0.0f;

					if (npc.Id == mCurrentlySelectedNpc)
					{
						mGameEventHandler->OnHumanNpcBehaviorChanged("Free_Swimming");
					}

					break;
				}
			}

			break;
		}
	}

	if (npc.Id == mCurrentlySelectedNpc)
	{
		mGameEventHandler->OnHumanNpcStateQuantityChanged(publishStateQuantity);
	}
}

bool Npcs::CheckAndMaintainHumanEquilibrium(
	ElementIndex primaryParticleIndex,
	ElementIndex secondaryParticleIndex,
	bool isRisingState,
	bool doMaintainEquilibrium,
	NpcParticles & particles,
	GameParameters const & /*gameParameters*/)
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
	float const staticDisplacementAngleCW = (-GameParameters::GravityDir).angleCw(humanVector);

	// Calculate CW angle that head would rotate by (relative to feet) due to relative velocity alone;
	// positive when new position is CW wrt old
	//
	// |   H
	// |  /
	// | /\
    // |/__L___H'
	//
	vec2f const relativeVelocityDisplacement = (particles.GetVelocity(secondaryParticleIndex) - particles.GetVelocity(primaryParticleIndex)) * GameParameters::SimulationTimeStepDuration;
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
		LogNpcDebug("Losing equilibrium because: StaticDisplacementAngleCW=", staticDisplacementAngleCW, " (Max=+/-", MaxStaticAngleForEquilibrium, ") RelativeVelocityAngleCW=", relativeVelocityAngleCW,
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
	StateType::KindSpecificStateType::HumanNpcStateType & humanState,
	StateType::NpcParticleStateType const & primaryParticleState,
	Ship const & /*shipMesh*/, // Will come useful when we'll *plan* the walk
	GameParameters const & gameParameters)
{
	assert(primaryParticleState.ConstrainedState.has_value());
	assert(humanState.CurrentBehavior == StateType::KindSpecificStateType::HumanNpcStateType::BehaviorType::Constrained_Walking);

	auto & walkingState = humanState.CurrentBehaviorState.Constrained_Walking;

	// 1. Check condition for potentially flipping: actual (relative) velocity opposite of walking direction,
	// or too small

	if (walkingState.CurrentWalkMagnitude != 0.0f)
	{
		float constexpr MinRelativeVelocityAgreementToAcceptWalk = 0.025f;
		float const relativeVelocityAgreement = primaryParticleState.ConstrainedState->MeshRelativeVelocity.dot(
			vec2f(
				humanState.CurrentFaceDirectionX * CalculateActualHumanWalkingAbsoluteSpeed(humanState, gameParameters),
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

	float const walkMagnitudeConvergenceRate = 0.03f + std::min(humanState.ResultantPanicLevel, 1.0f) * 0.15f;
	walkingState.CurrentWalkMagnitude += (1.0f - walkingState.CurrentWalkMagnitude) * walkMagnitudeConvergenceRate;

	LogNpcDebug("        currentWalkMagnitude: ", walkingState.CurrentWalkMagnitude);
}

void Npcs::OnHumanImpact(
	StateType & npc,
	bool /*isPrimaryParticle*/,
	vec2f const & /*impactVector*/,
	vec2f const & bounceEdgeNormal) const
{
	assert(npc.Kind == NpcKindType::Human);

	switch (npc.KindSpecificState.HumanNpcState.CurrentBehavior)
	{
		case StateType::KindSpecificStateType::HumanNpcStateType::BehaviorType::Constrained_Walking:
		{
			LogNpcDebug("OnHumanImpact: alignment=", bounceEdgeNormal.dot(vec2f(npc.KindSpecificState.HumanNpcState.CurrentFaceDirectionX, 0.0f)));

			// Check alignment of impact with walking direction; if hit => flip
			float constexpr MaxOppositionSlope = 0.5f;
			if (bounceEdgeNormal.dot(vec2f(npc.KindSpecificState.HumanNpcState.CurrentFaceDirectionX, 0.0f)) > MaxOppositionSlope
				&& npc.KindSpecificState.HumanNpcState.CurrentBehaviorState.Constrained_Walking.CurrentWalkMagnitude != 0.0f)
			{
				LogNpcDebug("OnHumanImpact: FLIP!");

				// Flip now
				FlipHumanWalk(npc.KindSpecificState.HumanNpcState, StrongTypedTrue<_DoImmediate>);
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
	StateType::KindSpecificStateType::HumanNpcStateType & humanState,
	DoImmediate doImmediate) const
{
	assert(humanState.CurrentBehavior == StateType::KindSpecificStateType::HumanNpcStateType::BehaviorType::Constrained_Walking);
	auto & walkingState = humanState.CurrentBehaviorState.Constrained_Walking;

	if (doImmediate.Value)
	{
		humanState.CurrentFaceDirectionX *= -1.0f;
		walkingState.CurrentWalkMagnitude = 0.0f;

		LogNpcDebug("Flipping walk: ", humanState.CurrentFaceDirectionX);

		walkingState.TargetFlipDecision = 0.0f;
		walkingState.CurrentFlipDecision = 0.0f;
	}
	else
	{
		walkingState.TargetFlipDecision = 1.0f;
	}
}

void Npcs::TransitionHumanToFree(
	StateType & npc,
	float currentSimulationTime)
{
	assert(npc.Kind == NpcKindType::Human);

	assert(npc.DipoleState.has_value());
	auto const & headPosition = mParticles.GetPosition(npc.DipoleState->SecondaryParticleState.ParticleIndex);
	auto const & feetPosition = mParticles.GetPosition(npc.PrimaryParticleState.ParticleIndex);

	// It's in water if both in water
	if (mParentWorld.GetOceanSurface().GetDepth(headPosition) > 0.0f
		&& mParentWorld.GetOceanSurface().GetDepth(feetPosition) > 0.0f)
	{
		npc.KindSpecificState.HumanNpcState.TransitionToState(StateType::KindSpecificStateType::HumanNpcStateType::BehaviorType::Free_InWater, currentSimulationTime);
		if (npc.Id == mCurrentlySelectedNpc)
		{
			mGameEventHandler->OnHumanNpcBehaviorChanged("Free_InWater");
		}
	}
	else
	{
		npc.KindSpecificState.HumanNpcState.TransitionToState(StateType::KindSpecificStateType::HumanNpcStateType::BehaviorType::Free_Aerial, currentSimulationTime);
		if (npc.Id == mCurrentlySelectedNpc)
		{
			mGameEventHandler->OnHumanNpcBehaviorChanged("Free_Aerial");
		}
	}
}

float Npcs::CalculateActualHumanWalkingAbsoluteSpeed(
	StateType::KindSpecificStateType::HumanNpcStateType & humanState,
	GameParameters const & gameParameters) const
{
	assert(humanState.CurrentBehavior == StateType::KindSpecificStateType::HumanNpcStateType::BehaviorType::Constrained_Walking);

	return std::min(
		gameParameters.HumanNpcWalkingSpeed
		* humanState.CurrentBehaviorState.Constrained_Walking.CurrentWalkMagnitude // Note that this is the only one that might be zero
		* (1.0f + humanState.ResultantPanicLevel * 3.0f),
		4.0f); // Absolute cap
}

}