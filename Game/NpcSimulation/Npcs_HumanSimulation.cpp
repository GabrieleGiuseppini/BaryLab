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

Npcs::StateType::KindSpecificStateType::HumanNpcStateType::BehaviorType Npcs::CalculateHumanBehavior(StateType & npc)
{
	assert(npc.ParticleMesh.Particles.size() == 2);
	auto const & primaryParticleState = npc.ParticleMesh.Particles[0];
	auto const & secondaryParticleState = npc.ParticleMesh.Particles[1];

	if (!primaryParticleState.ConstrainedState.has_value()
		&& !secondaryParticleState.ConstrainedState.has_value())
	{
		// Whole Human is free
		return StateType::KindSpecificStateType::HumanNpcStateType::BehaviorType::Free_Aerial;
	}
	else
	{
		// Human is constrained
		return StateType::KindSpecificStateType::HumanNpcStateType::BehaviorType::Constrained_Aerial;
	}
}

void Npcs::UpdateFurniture(
	StateType & npc,
	float currentSimulationTime,
	Ship & /*homeShip*/ ,
	GameParameters const & /*gameParameters*/)
{
	assert(npc.Kind == NpcKindType::Furniture);
	auto & furnitureState = npc.KindSpecificState.FurnitureNpcState;
	using FurnitureNpcStateType = StateType::KindSpecificStateType::FurnitureNpcStateType;

	//
	// Process furniture
	//

	LogNpcDebug("CurrentBehavior: ", int(furnitureState.CurrentBehavior));

	switch (furnitureState.CurrentBehavior)
	{
		case FurnitureNpcStateType::BehaviorType::Default:
		{
			// Nop
			break;
		}

		case FurnitureNpcStateType::BehaviorType::BeingRemoved:
		{
			// See if we're done
			float const elapsed = currentSimulationTime - furnitureState.CurrentStateTransitionSimulationTimestamp;
			if (elapsed >= FurnitureRemovalDuration)
			{
				// Deferred removal
				assert(std::find(mDeferredRemovalNpcs.cbegin(), mDeferredRemovalNpcs.cend(), npc.Id) == mDeferredRemovalNpcs.cend());
				mDeferredRemovalNpcs.emplace_back(npc.Id);
			}

			break;
		}
	}
}

void Npcs::UpdateHuman(
	StateType & npc,
	float currentSimulationTime,
	Ship & homeShip,
	GameParameters const & gameParameters)
{
	assert(npc.ParticleMesh.Particles.size() == 2);
	auto & primaryParticleState = npc.ParticleMesh.Particles[0];
	auto & secondaryParticleState = npc.ParticleMesh.Particles[1];

	assert(npc.Kind == NpcKindType::Human);
	auto & humanState = npc.KindSpecificState.HumanNpcState;
	using HumanNpcStateType = StateType::KindSpecificStateType::HumanNpcStateType;

	//
	// Constants
	//

	unsigned int constexpr LowFrequencyUpdatePeriod = 4;

	float constexpr MaxRelativeVelocityMagnitudeForEquilibrium = 3.0f; // So high because we slip a lot while we try to stand up, and thus need to be immune to ourselves

	//
	// Reset pulse state variables - variables that we set here and are meant
	// to last for one frame only
	//

	if (primaryParticleState.ConstrainedState.has_value())
	{
		primaryParticleState.ConstrainedState->GhostParticlePulse = false;
	}
	humanState.EquilibriumTorque = 0.0f;

	//
	// Update panic
	//

	humanState.ResultantPanicLevel =
		humanState.OnFirePanicLevel
		+ humanState.BombProximityPanicLevel
		+ humanState.IncomingWaterProximityPanicLevel
		+ humanState.MiscPanicLevel
		+ mGeneralizedPanicLevel;

	// Decay

	humanState.OnFirePanicLevel -= humanState.OnFirePanicLevel * 0.01f;
	humanState.BombProximityPanicLevel -= humanState.BombProximityPanicLevel * 0.0025f;
	humanState.IncomingWaterProximityPanicLevel -= humanState.IncomingWaterProximityPanicLevel * 0.005f;
	humanState.MiscPanicLevel -= humanState.MiscPanicLevel * 0.009f;

	//
	// Process human
	//

#ifdef BARYLAB_PROBING
	std::optional<std::tuple<std::string, std::string>> publishStateQuantity;
#endif

	bool const isFree = !primaryParticleState.ConstrainedState.has_value();

	LogNpcDebug("CurrentBehavior: ", int(humanState.CurrentBehavior));

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
				TransitionHumanBehaviorToFree(npc, currentSimulationTime);

				break;
			}

			// Check conditions for falling/rising

			float fallingTarget;
			float risingTarget;

			bool const isHeadOnFloor = secondaryParticleState.ConstrainedState.has_value() && secondaryParticleState.ConstrainedState->CurrentVirtualFloor.has_value();
			bool const areFeetOnFloor = primaryParticleState.ConstrainedState.has_value() && primaryParticleState.ConstrainedState->CurrentVirtualFloor.has_value();
			if (isHeadOnFloor || areFeetOnFloor)
			{
				vec2f const floorVector = (primaryParticleState.ConstrainedState.has_value() && primaryParticleState.ConstrainedState->CurrentVirtualFloor.has_value())
					? homeShip.GetTriangles().GetSubSpringVector(
						primaryParticleState.ConstrainedState->CurrentVirtualFloor->TriangleElementIndex,
						primaryParticleState.ConstrainedState->CurrentVirtualFloor->EdgeOrdinal,
						homeShip.GetPoints())
					: vec2f(1.0f, 0.0); // H arbitrarily
				float const headVelocityAlongFloor = secondaryParticleState.GetApplicableVelocity(mParticles).dot(floorVector);
				float const feetVelocityAlongFloor = primaryParticleState.GetApplicableVelocity(mParticles).dot(floorVector);

				float constexpr MinVelocityMagnitudeForFalling = 0.05f;

				if (std::abs(headVelocityAlongFloor) >= MinVelocityMagnitudeForFalling
					|| std::abs(feetVelocityAlongFloor) >= MinVelocityMagnitudeForFalling)
				{
					// Likely falling
					fallingTarget = 1.0f;

					// Definitely not rising
					risingTarget = 0.0f;
				}
				else
				{
					// We're quite still...

					// ...likely rising
					risingTarget = 1.0f;

					// Definitely not falling
					fallingTarget = 0.0f;
				}
			}
			else
			{
				// Completely in the air, no transition
				fallingTarget = 0.0f;
				risingTarget = 0.0f;
			}

			// Progress to falling

			float constexpr ToFallingConvergenceRate = 0.75f; // Very high! We do this just to survive micro-instants

			humanState.CurrentBehaviorState.Constrained_Aerial.ProgressToFalling +=
				(fallingTarget - humanState.CurrentBehaviorState.Constrained_Aerial.ProgressToFalling)
				* ToFallingConvergenceRate;

#ifdef BARYLAB_PROBING
			publishStateQuantity = std::make_tuple("ProgressToFalling", std::to_string(humanState.CurrentBehaviorState.Constrained_Aerial.ProgressToFalling));
#endif

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

#ifdef BARYLAB_PROBING
				if (npc.Id == mCurrentlySelectedNpc)
				{
					mGameEventHandler->OnHumanNpcBehaviorChanged("Constrained_Falling");
				}
#endif

				break;
			}

			// Progress to rising

			float constexpr ToRisingConvergenceRate = 0.5f;

			humanState.CurrentBehaviorState.Constrained_Aerial.ProgressToRising +=
				(risingTarget - humanState.CurrentBehaviorState.Constrained_Aerial.ProgressToRising)
				* ToRisingConvergenceRate;

#ifdef BARYLAB_PROBING
			if (fallingTarget == 0.0f)
				publishStateQuantity = std::make_tuple("ProgressToRising", std::to_string(humanState.CurrentBehaviorState.Constrained_Aerial.ProgressToRising));
#endif

			if (IsAtTarget(humanState.CurrentBehaviorState.Constrained_Aerial.ProgressToRising, 1.0f))
			{
				// Transition

				humanState.TransitionToState(HumanNpcStateType::BehaviorType::Constrained_Rising, currentSimulationTime);

#ifdef BARYLAB_PROBING
				if (npc.Id == mCurrentlySelectedNpc)
				{
					mGameEventHandler->OnHumanNpcBehaviorChanged("Constrained_Rising");
				}
#endif

				break;
			}

			// Check if moved to water

			// It's in water if at least one in water
			if (mParticles.GetAnyWaterness(primaryParticleState.ParticleIndex) > 0.5f
				|| mParticles.GetAnyWaterness(secondaryParticleState.ParticleIndex) > 0.5f)
			{
				// Transition

				humanState.TransitionToState(HumanNpcStateType::BehaviorType::Constrained_InWater, currentSimulationTime);

#ifdef BARYLAB_PROBING
				if (npc.Id == mCurrentlySelectedNpc)
				{
					mGameEventHandler->OnHumanNpcBehaviorChanged("Constrained_InWater");
				}
#endif

				break;
			}

			break;
		}

		case HumanNpcStateType::BehaviorType::Constrained_Falling:
		{
			//
			// From here:
			//	- Free_X
			//	- Constrained_PreRising - when "mostly" with feet on floor and almost still
			//	- Constrained_Aerial - when "consistently" with feet and head in air
			//

			if (isFree)
			{
				// Transition
				TransitionHumanBehaviorToFree(npc, currentSimulationTime);

				break;
			}

			assert(primaryParticleState.ConstrainedState.has_value()); // Or else we'd have been free

			bool const areFeetOnFloor = primaryParticleState.ConstrainedState->CurrentVirtualFloor.has_value();
			bool const isHeadOnFloor = secondaryParticleState.ConstrainedState.has_value() && secondaryParticleState.ConstrainedState->CurrentVirtualFloor.has_value();

			// Advance towards pre-rising

			float toPreRisingIncrement;
			float constexpr MaxRelativeVelocityForPreRising = 0.5f;
			if (areFeetOnFloor
				&& primaryParticleState.ConstrainedState->MeshRelativeVelocity.length() < MaxRelativeVelocityForPreRising
				&& secondaryParticleState.GetApplicableVelocity(mParticles).length() < MaxRelativeVelocityForPreRising)
			{
				toPreRisingIncrement = 1.0f;
			}
			else
			{
				toPreRisingIncrement = -1.0f;
			}

			humanState.CurrentBehaviorState.Constrained_Falling.ProgressToPreRising = std::max(
				humanState.CurrentBehaviorState.Constrained_Falling.ProgressToPreRising + toPreRisingIncrement,
				0.0f);

			float const toPreRisingTarget =
				20.0f
				- std::min(humanState.ResultantPanicLevel, 1.0f) * 10.0f;

#ifdef BARYLAB_PROBING
			if (toPreRisingIncrement > 0.0f)
				publishStateQuantity = std::make_tuple("ProgressToPreRising", std::to_string(humanState.CurrentBehaviorState.Constrained_Falling.ProgressToPreRising / toPreRisingTarget));
#endif

			if (humanState.CurrentBehaviorState.Constrained_Falling.ProgressToPreRising >= toPreRisingTarget)
			{
				// Transition

				humanState.TransitionToState(HumanNpcStateType::BehaviorType::Constrained_PreRising, currentSimulationTime);

#ifdef BARYLAB_PROBING
				if (npc.Id == mCurrentlySelectedNpc)
				{
					mGameEventHandler->OnHumanNpcBehaviorChanged("Constrained_PreRising");
				}
#endif

				break;
			}

			// Check conditions for aerial

			if (!areFeetOnFloor && !isHeadOnFloor)
			{
				// Advance towards aerial

				float constexpr ToAerialConvergenceRate = 0.35f;

				humanState.CurrentBehaviorState.Constrained_Falling.ProgressToAerial +=
					(1.0f - humanState.CurrentBehaviorState.Constrained_Falling.ProgressToAerial)
					* ToAerialConvergenceRate;

#ifdef BARYLAB_PROBING
				if (toPreRisingIncrement <= 0.0f)
					publishStateQuantity = std::make_tuple("ProgressToAerial", std::to_string(humanState.CurrentBehaviorState.Constrained_Falling.ProgressToAerial));
#endif

				if (IsAtTarget(humanState.CurrentBehaviorState.Constrained_Falling.ProgressToAerial, 1.0f))
				{
					// Transition

					humanState.TransitionToState(HumanNpcStateType::BehaviorType::Constrained_Aerial, currentSimulationTime);

#ifdef BARYLAB_PROBING
					if (npc.Id == mCurrentlySelectedNpc)
					{
						mGameEventHandler->OnHumanNpcBehaviorChanged("Constrained_Aerial");
					}
#endif

					break;
				}
			}
			else
			{
				// Reset progress to aerial

				humanState.CurrentBehaviorState.Constrained_Falling.ProgressToAerial = 0.0f;

#ifdef BARYLAB_PROBING
				if (toPreRisingIncrement <= 0.0f)
					publishStateQuantity = std::make_tuple("ProgressToAerial", std::to_string(humanState.CurrentBehaviorState.Constrained_Falling.ProgressToAerial));
#endif
			}

			break;
		}

		case HumanNpcStateType::BehaviorType::Constrained_KnockedOut:
		{
			//
			// From here:
			//	- Free_X
			//	- Constrained_PreRising - when "mostly" with feet on floor
			//	- Constrained_Aerial - when "consistently" feet and head in air
			//

			if (isFree)
			{
				// Transition
				TransitionHumanBehaviorToFree(npc, currentSimulationTime);

				break;
			}

			assert(primaryParticleState.ConstrainedState.has_value()); // Or else we'd have been free

			bool const areFeetOnFloor = primaryParticleState.ConstrainedState->CurrentVirtualFloor.has_value();
			bool const isHeadOnFloor = secondaryParticleState.ConstrainedState.has_value() && secondaryParticleState.ConstrainedState->CurrentVirtualFloor.has_value();

			// Advance towards pre-rising

			float toPreRisingIncrement;
			if (areFeetOnFloor
				&& primaryParticleState.ConstrainedState->MeshRelativeVelocity.length() < MaxRelativeVelocityMagnitudeForEquilibrium
				&& secondaryParticleState.GetApplicableVelocity(mParticles).length() < MaxRelativeVelocityMagnitudeForEquilibrium)
			{
				toPreRisingIncrement = 1.0f;
			}
			else
			{
				toPreRisingIncrement = -1.0f;
			}

			humanState.CurrentBehaviorState.Constrained_KnockedOut.ProgressToPreRising = std::max(
				humanState.CurrentBehaviorState.Constrained_KnockedOut.ProgressToPreRising + toPreRisingIncrement,
				0.0f);

			// 30-40-50 (panic 0) -> 15-20-25 (panic +INF)
			float const toPreRisingTarget =
				(40.0f + npc.RandomNormalizedUniformSeed * 10.0f)
				/ (1.0f + std::min(humanState.ResultantPanicLevel, 1.0f));

#ifdef BARYLAB_PROBING
			publishStateQuantity = std::make_tuple("ProgressToPreRising", std::to_string(humanState.CurrentBehaviorState.Constrained_KnockedOut.ProgressToPreRising / toPreRisingTarget));
#endif

			if (humanState.CurrentBehaviorState.Constrained_KnockedOut.ProgressToPreRising >= toPreRisingTarget)
			{
				// Small hack: given that we've established we can rise (and thus we've been static
				// for a while in our current position), see if may be we're hanging by the feet onto
				// a floor, with the head hanging down; if so, free the feet with a ghost pulse

				// Feet to head == head - feet
				vec2f const humanDir = (mParticles.GetPosition(secondaryParticleState.ParticleIndex) - mParticles.GetPosition(primaryParticleState.ParticleIndex)).normalise_approx();

				if (areFeetOnFloor
					&& !isHeadOnFloor
					&& humanDir.y < -0.6f) // ~53deg
				{
					// Free the feet
					assert(primaryParticleState.ConstrainedState.has_value()); // Otherwise we'd have become free
					primaryParticleState.ConstrainedState->GhostParticlePulse = true;

					LogNpcDebug("Pulsed GhostParticle");

					// Since we don't transition out, reset state
					humanState.CurrentBehaviorState.Constrained_KnockedOut.Reset();
				}
				else
				{
					// Transition

					humanState.TransitionToState(HumanNpcStateType::BehaviorType::Constrained_PreRising, currentSimulationTime);

#ifdef BARYLAB_PROBING
					if (npc.Id == mCurrentlySelectedNpc)
					{
						mGameEventHandler->OnHumanNpcBehaviorChanged("Constrained_PreRising");
					}
#endif
				}

				break;
			}

			// Check conditions for aerial

			if (!areFeetOnFloor && !isHeadOnFloor)
			{
				// Advance towards aerial

				float constexpr ToAerialConvergenceRate = 0.2f;

				humanState.CurrentBehaviorState.Constrained_KnockedOut.ProgressToAerial +=
					(1.0f - humanState.CurrentBehaviorState.Constrained_KnockedOut.ProgressToAerial)
					* ToAerialConvergenceRate;

#ifdef BARYLAB_PROBING
				publishStateQuantity = std::make_tuple("ProgressToAerial", std::to_string(humanState.CurrentBehaviorState.Constrained_KnockedOut.ProgressToAerial));
#endif

				if (IsAtTarget(humanState.CurrentBehaviorState.Constrained_KnockedOut.ProgressToAerial, 1.0f))
				{
					// Transition

					humanState.TransitionToState(HumanNpcStateType::BehaviorType::Constrained_Aerial, currentSimulationTime);

#ifdef BARYLAB_PROBING
					if (npc.Id == mCurrentlySelectedNpc)
					{
						mGameEventHandler->OnHumanNpcBehaviorChanged("Constrained_Aerial");
					}
#endif

					break;
				}
			}
			else
			{
				// Reset progress to aerial

				humanState.CurrentBehaviorState.Constrained_KnockedOut.ProgressToAerial = 0.0f;

#ifdef BARYLAB_PROBING
				if (toPreRisingIncrement <= 0.0f)
					publishStateQuantity = std::make_tuple("ProgressToAerial", std::to_string(humanState.CurrentBehaviorState.Constrained_KnockedOut.ProgressToAerial));
#endif
			}

			break;
		}

		case HumanNpcStateType::BehaviorType::Constrained_PreRising:
		{
			//
			// From here:
			//	- Free_X
			//	- Constrained_Rising - when "mostly" with feet on floor
			//	- Constrained_Aerial - when "consistently" feet and head in air
			//

			if (isFree)
			{
				// Transition
				TransitionHumanBehaviorToFree(npc, currentSimulationTime);

				break;
			}

			assert(primaryParticleState.ConstrainedState.has_value()); // Or else we'd have been free

			bool const areFeetOnFloor = primaryParticleState.ConstrainedState->CurrentVirtualFloor.has_value();
			bool const isHeadOnFloor = secondaryParticleState.ConstrainedState.has_value() && secondaryParticleState.ConstrainedState->CurrentVirtualFloor.has_value();

			// Advance towards rising

			float toRisingIncrement;
			if (areFeetOnFloor
				&& primaryParticleState.ConstrainedState->MeshRelativeVelocity.length() < MaxRelativeVelocityMagnitudeForEquilibrium
				&& secondaryParticleState.GetApplicableVelocity(mParticles).length() < MaxRelativeVelocityMagnitudeForEquilibrium)
			{
				toRisingIncrement = 1.0f;
			}
			else
			{
				toRisingIncrement = -1.0f;
			}

			humanState.CurrentBehaviorState.Constrained_PreRising.ProgressToRising = std::max(
				humanState.CurrentBehaviorState.Constrained_PreRising.ProgressToRising + toRisingIncrement,
				0.0f);

			// 10-13-16
			float const toRisingTarget =
				13.0f
				+ npc.RandomNormalizedUniformSeed * 3.0f;

#ifdef BARYLAB_PROBING
			publishStateQuantity = std::make_tuple("ProgressToRising", std::to_string(humanState.CurrentBehaviorState.Constrained_PreRising.ProgressToRising / toRisingTarget));
#endif

			if (humanState.CurrentBehaviorState.Constrained_PreRising.ProgressToRising >= toRisingTarget)
			{
				// Transition

				humanState.TransitionToState(HumanNpcStateType::BehaviorType::Constrained_Rising, currentSimulationTime);

#ifdef BARYLAB_PROBING
				if (npc.Id == mCurrentlySelectedNpc)
				{
					mGameEventHandler->OnHumanNpcBehaviorChanged("Constrained_Rising");
				}
#endif

				break;
			}

			// Check conditions for aerial

			if (!areFeetOnFloor && !isHeadOnFloor)
			{
				// Advance towards aerial

				float constexpr ToAerialConvergenceRate = 0.2f;

				humanState.CurrentBehaviorState.Constrained_PreRising.ProgressToAerial +=
					(1.0f - humanState.CurrentBehaviorState.Constrained_PreRising.ProgressToAerial)
					* ToAerialConvergenceRate;

#ifdef BARYLAB_PROBING
				publishStateQuantity = std::make_tuple("ProgressToAerial", std::to_string(humanState.CurrentBehaviorState.Constrained_PreRising.ProgressToAerial));
#endif

				if (IsAtTarget(humanState.CurrentBehaviorState.Constrained_PreRising.ProgressToAerial, 1.0f))
				{
					// Transition

					humanState.TransitionToState(HumanNpcStateType::BehaviorType::Constrained_Aerial, currentSimulationTime);

#ifdef BARYLAB_PROBING
					if (npc.Id == mCurrentlySelectedNpc)
					{
						mGameEventHandler->OnHumanNpcBehaviorChanged("Constrained_Aerial");
					}
#endif

					break;
				}
			}
			else
			{
				// Reset progress to aerial

				humanState.CurrentBehaviorState.Constrained_PreRising.ProgressToAerial = 0.0f;

#ifdef BARYLAB_PROBING
				if (toRisingIncrement <= 0.0f)
					publishStateQuantity = std::make_tuple("ProgressToAerial", std::to_string(humanState.CurrentBehaviorState.Constrained_PreRising.ProgressToAerial));
#endif
			}

			break;
		}

		case HumanNpcStateType::BehaviorType::Constrained_Rising:
		case HumanNpcStateType::BehaviorType::Constrained_Equilibrium:
		case HumanNpcStateType::BehaviorType::Constrained_Walking:
		case HumanNpcStateType::BehaviorType::Constrained_WalkingUndecided:
		{
			if (isFree)
			{
				// Transition
				TransitionHumanBehaviorToFree(npc, currentSimulationTime);

				break;
			}

			// Check progress to electrified, bomb

			if (mCurrentSimulationSequenceNumber.IsStepOf(npc.Id % LowFrequencyUpdatePeriod, LowFrequencyUpdatePeriod))
			{
				if (humanState.CurrentBehavior == HumanNpcStateType::BehaviorType::Constrained_Equilibrium
					|| humanState.CurrentBehavior == HumanNpcStateType::BehaviorType::Constrained_Walking
					|| humanState.CurrentBehavior == HumanNpcStateType::BehaviorType::Constrained_WalkingUndecided)
				{
					// Check electrification

					if (IsElectrified(npc, homeShip))
					{
						// Transition

						humanState.TransitionToState(HumanNpcStateType::BehaviorType::Constrained_Electrified, currentSimulationTime);

						// Face: rnd/0.0
						humanState.CurrentFaceOrientation = GameRandomEngine::GetInstance().GenerateUniformBoolean(0.5f) ? +1.0f : -1.0f;
						humanState.CurrentFaceDirectionX = 0.0f;

						// Keep torque
						humanState.EquilibriumTorque = 1.0f;

#ifdef BARYLAB_PROBING
						if (npc.Id == mCurrentlySelectedNpc)
						{
							mGameEventHandler->OnHumanNpcBehaviorChanged("Constrained_Electrified");
						}
#endif

						break;
					}

					// Check bomb panic

					if (HasBomb(npc, homeShip))
					{
						if (humanState.BombProximityPanicLevel < 0.6f)
						{
							// Time to flip
							humanState.CurrentFaceDirectionX *= -1.0f;
						}

						// Panic
						humanState.BombProximityPanicLevel = 1.0f;

						// Continue
					}
				}
			}

			if (humanState.CurrentBehavior == HumanNpcStateType::BehaviorType::Constrained_Walking)
			{
				// Check fire panic

				if (npc.CombustionState.has_value())
				{
					if (humanState.OnFirePanicLevel < 0.6f)
					{
						// Time to flip
						// Note: not via helper as we don't want to lose walking magnitude
						humanState.CurrentFaceDirectionX *= -1.0f;
					}

					// Panic
					humanState.OnFirePanicLevel = 1.0f;

					// Continue
				}

				// Check incoming water panic

				if (humanState.BombProximityPanicLevel < 0.6f
					&& humanState.OnFirePanicLevel < 0.6f) // Just make check subordinate to other panics, don't want to flip headlessly
				{
					vec2f const & waterMeshVelocity = mParticles.GetMeshWaterVelocity(primaryParticleState.ParticleIndex);
					if (std::fabsf(waterMeshVelocity.x) > 0.1f)
					{
						// Water velocity passes test

						if (waterMeshVelocity.x * humanState.CurrentFaceDirectionX <= 0.0f
							&& humanState.IncomingWaterProximityPanicLevel < 0.6f)
						{
							// Time to flip
							// Note: not via helper as we don't want to lose walking magnitude
							humanState.CurrentFaceDirectionX *= -1.0f;
						}

						// Panic
						humanState.IncomingWaterProximityPanicLevel = 1.0f;

						// Continue
					}
				}
			}

			// Check progress to walking

			bool doTransitionToWalking = false;

			bool const areFeetOnFloor = primaryParticleState.ConstrainedState.has_value() && primaryParticleState.ConstrainedState->CurrentVirtualFloor.has_value();
			if (humanState.CurrentBehavior == HumanNpcStateType::BehaviorType::Constrained_Equilibrium)
			{
				if (areFeetOnFloor)
				{
					// Advance towards walking

					float const toWalkingConvergenceRate = 0.12f + std::min(humanState.ResultantPanicLevel, 1.0f) * 0.12f;
					humanState.CurrentBehaviorState.Constrained_Equilibrium.ProgressToWalking +=
						(1.0f - humanState.CurrentBehaviorState.Constrained_Equilibrium.ProgressToWalking)
						* toWalkingConvergenceRate;

#ifdef BARYLAB_PROBING
					publishStateQuantity = std::make_tuple("ProgressToWalking", std::to_string(humanState.CurrentBehaviorState.Constrained_Equilibrium.ProgressToWalking));
#endif

					if (IsAtTarget(humanState.CurrentBehaviorState.Constrained_Equilibrium.ProgressToWalking, 1.0f))
					{
						// Transition
						doTransitionToWalking = true;
					}
				}
			}
			else if (humanState.CurrentBehavior == HumanNpcStateType::BehaviorType::Constrained_WalkingUndecided)
			{
				if (areFeetOnFloor)
				{
					float const elapsed = currentSimulationTime - humanState.CurrentStateTransitionSimulationTimestamp;
					if (elapsed >= WalkingUndecidedDuration)
					{
						// Transition
						doTransitionToWalking = true;
					}
				}
			}

			if (doTransitionToWalking)
			{
				// See if we can choose a direction
				auto const direction = DecideWalkingDirection(npc, primaryParticleState, homeShip);
				if (!direction.has_value())
				{
					// Transition to undecided

					humanState.TransitionToState(HumanNpcStateType::BehaviorType::Constrained_WalkingUndecided, currentSimulationTime);

					// Face: 1.0/0.0
					humanState.CurrentFaceOrientation = 1.0f;
					humanState.CurrentFaceDirectionX = 0.0f;

					// Keep torque

#ifdef BARYLAB_PROBING
					if (npc.Id == mCurrentlySelectedNpc)
					{
						mGameEventHandler->OnHumanNpcBehaviorChanged("Constrained_WalkingUndecided");
					}
#endif
				}
				else
				{
					// Transition to walking

					humanState.TransitionToState(HumanNpcStateType::BehaviorType::Constrained_Walking, currentSimulationTime);

					// Face: 0/chosen
					humanState.CurrentFaceOrientation = 0.0f;
					humanState.CurrentFaceDirectionX = *direction;

					// Keep torque

#ifdef BARYLAB_PROBING
					if (npc.Id == mCurrentlySelectedNpc)
					{
						mGameEventHandler->OnHumanNpcBehaviorChanged("Constrained_Walking");
					}
#endif
				}

				break;
			}

			//
			// Check conditions to stay & maintain equilibrium
			//

			// a. Feet on-floor

			bool isStateMaintained;
			if (!areFeetOnFloor)
			{
				float toTerminateEquilibriumConvergenceRate;
				if (humanState.CurrentBehavior == HumanNpcStateType::BehaviorType::Constrained_Walking)
				{
					// When walking, we want to be a bit more tolerant about "losing the edge";
					// this is a quite important parameter: it's the duration through which we tolerate temporarily losing contact
					// with the ground.
					// If we're walking normally we can bear having a short tolerance; we only need a long tolerance
					// when we're walking "fast".
					// Walking speed rel == 1.0 => 0.25
					// Walking speed rel == 1.5 => 0.1
					float const relWalkingSpeed = CalculateHumanWalkingSpeedAdjustment(humanState);
					toTerminateEquilibriumConvergenceRate = Clamp(
						0.25f - (relWalkingSpeed - 1.0f) / (1.5f - 1.0f) * (0.25f - 0.1f),
						0.1f,
						0.25f);
				}
				else
				{
					// When not walking, we lose equilibrium very fast!
					toTerminateEquilibriumConvergenceRate = 0.25f;
				}

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

				float const primaryMeshRelativeVelocityAlongWalkDir = primaryParticleState.ConstrainedState->MeshRelativeVelocity.dot(idealWalkVelocityDir);

				LogNpcDebug("Walk update: mesh-relative velocity check: idealWalkVelocity=", idealWalkVelocityDir * CalculateActualHumanWalkingAbsoluteSpeed(humanState),
					" meshRelativeVelocity=", primaryParticleState.ConstrainedState->MeshRelativeVelocity, " (along dir=", primaryMeshRelativeVelocityAlongWalkDir, ")");

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
					humanState,
					areFeetOnFloor, // doMaintainEquilibrium
					mParticles,
					gameParameters))
			{
				// Transition to aerial/falling, depending on whether we're on an edge

				LogNpcDebug("Going to Constrained_X; primary's barycentric coords: ",
					primaryParticleState.ConstrainedState.has_value() ? primaryParticleState.ConstrainedState->CurrentBCoords.BCoords.toString() : "N/A",
					" primary's relative velocity mag: ", primaryParticleState.ConstrainedState.has_value() ? std::to_string(primaryParticleState.ConstrainedState->MeshRelativeVelocity.length()) : "N/A",
					" (max=", MaxRelativeVelocityMagnitudeForEquilibrium, ")");

				if (areFeetOnFloor)
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

#ifdef BARYLAB_PROBING
					if (npc.Id == mCurrentlySelectedNpc)
					{
						mGameEventHandler->OnHumanNpcBehaviorChanged("Constrained_Falling");
					}
#endif
				}
				else
				{
					// Aerial

					humanState.TransitionToState(HumanNpcStateType::BehaviorType::Constrained_Aerial, currentSimulationTime);

#ifdef BARYLAB_PROBING
					if (npc.Id == mCurrentlySelectedNpc)
					{
						mGameEventHandler->OnHumanNpcBehaviorChanged("Constrained_Aerial");
					}
#endif
				}

				break;
			}

			//
			// Update state now
			//

			if (humanState.CurrentBehavior == HumanNpcStateType::BehaviorType::Constrained_Rising)
			{
				// Check if reached alignment (note: here so that we may keep torque as we'll be transitioning to Equilibrium)

				float const alignment = CalculateSpringVerticalAlignment(primaryParticleState.ParticleIndex, secondaryParticleState.ParticleIndex, mParticles);

#ifdef BARYLAB_PROBING
				publishStateQuantity = std::make_tuple("Alignment", std::to_string(alignment));
#endif

				if (AreAlmostEqual(alignment, 1.0f, 0.004f))
				{
					// Transition

					humanState.TransitionToState(HumanNpcStateType::BehaviorType::Constrained_Equilibrium, currentSimulationTime);

#ifdef BARYLAB_PROBING
					if (npc.Id == mCurrentlySelectedNpc)
					{
						mGameEventHandler->OnHumanNpcBehaviorChanged("Constrained_Equilibrium");
					}
#endif
				}
			}
			else if (humanState.CurrentBehavior == HumanNpcStateType::BehaviorType::Constrained_Equilibrium)
			{
				// Nop
			}
			else if (humanState.CurrentBehavior == HumanNpcStateType::BehaviorType::Constrained_Walking)
			{
				auto & walkingState = humanState.CurrentBehaviorState.Constrained_Walking;

				if (areFeetOnFloor) // Note: no need to silence walk as we don't apply walk displacement in inertial (i.e. not-on-edge) case
				{
					//
					// Good to continue walking!
					//

					// Forward-looking: since it might be expensive, we only do it when we cross the
					// currently-walked triangle edge into its last half

					assert(primaryParticleState.ConstrainedState.has_value());
					auto const & primaryParticleConstrainedState = *primaryParticleState.ConstrainedState;
					int const walkedEdgeOrdinal = primaryParticleConstrainedState.CurrentVirtualFloor->EdgeOrdinal;

					if (primaryParticleConstrainedState.CurrentVirtualFloor.has_value()
						&& primaryParticleConstrainedState.CurrentVirtualFloor->TriangleElementIndex == primaryParticleConstrainedState.CurrentBCoords.TriangleElementIndex) // For safety
					{
						if (!walkingState.LastHalfTriangleEdge.has_value()
							|| primaryParticleConstrainedState.CurrentVirtualFloor->TriangleElementIndex != walkingState.LastHalfTriangleEdge->TriangleElementIndex
							|| walkedEdgeOrdinal != walkingState.LastHalfTriangleEdge->EdgeOrdinal)
						{
							// We have changed edge
							walkingState.LastHalfTriangleEdge = {
								primaryParticleConstrainedState.CurrentVirtualFloor->TriangleElementIndex,
								walkedEdgeOrdinal,
								primaryParticleConstrainedState.CurrentBCoords.BCoords[walkedEdgeOrdinal],
								0 };
						}
						else
						{
							// We were last in this triangle and on this edge

							// Check whether we have crossed the last half
							float const faceDirectionX = npc.KindSpecificState.HumanNpcState.CurrentFaceDirectionX;
							bool const hasCrossed =
								!IsInLastHalfOfEdge(walkingState.LastHalfTriangleEdge->EdgeBCoord, faceDirectionX)
								&& IsInLastHalfOfEdge(primaryParticleConstrainedState.CurrentBCoords.BCoords[walkedEdgeOrdinal], faceDirectionX);

							// Update state
							walkingState.LastHalfTriangleEdge->EdgeBCoord = primaryParticleConstrainedState.CurrentBCoords.BCoords[walkedEdgeOrdinal];

							if (hasCrossed)
							{
								//
								// Just crossed
								//

								LogNpcDebug("Crossed half-edge");

								//
								// Plan ahead
								//

								if (!CanWalkInDirection(npc, *primaryParticleConstrainedState.CurrentVirtualFloor, faceDirectionX, homeShip))
								{
									LogNpcDebug("Cannot continue walking along face direction");

									//
									// Flip
									//
									// We'll figure out whether we're stuck also in other direction after we cross the half-triangle again
									// enough times
									//
									// Also flip when under panic, that will make for funny fall-on-back
									//

									// Check if we've flipped too many times
									if (walkingState.LastHalfTriangleEdge->NumberOfTimesHasFlipped >= 2)
									{
										LogNpcDebug("Crossed half-edge too many times - transitioning to WalkingUndecided")

										// Transition to undecided

										humanState.TransitionToState(HumanNpcStateType::BehaviorType::Constrained_WalkingUndecided, currentSimulationTime);

										// Face: 1.0/0
										humanState.CurrentFaceOrientation = 1.0f;
										humanState.CurrentFaceDirectionX = 0.0f;

										// Keep torque

#ifdef BARYLAB_PROBING
										if (npc.Id == mCurrentlySelectedNpc)
										{
											mGameEventHandler->OnHumanNpcBehaviorChanged("Constrained_WalkingUndecided");
										}
#endif

										break;
									}

									FlipHumanWalk(humanState, StrongTypedTrue<_DoImmediate>);

									//
									// Update number of times we've flipped
									//

									++walkingState.LastHalfTriangleEdge->NumberOfTimesHasFlipped;
								}
							}
						}
					}

					//
					// Impart walk displacement & run walking state machine
					//

					RunWalkingHumanStateMachine(
						humanState,
						primaryParticleState,
						homeShip,
						gameParameters);
				}

#ifdef BARYLAB_PROBING
				if (humanState.CurrentBehaviorState.Constrained_Walking.CurrentFlipDecision != 0.0f)
					publishStateQuantity = std::make_tuple("WalkFlip", std::to_string(humanState.CurrentBehaviorState.Constrained_Walking.CurrentFlipDecision));
				else
					publishStateQuantity = std::make_tuple("EquilibriumTermination", std::to_string(humanState.CurrentEquilibriumSoftTerminationDecision));
#endif
			}
			else
			{
				assert(humanState.CurrentBehavior == HumanNpcStateType::BehaviorType::Constrained_WalkingUndecided);

				// Nop
			}

			break;
		}

		case HumanNpcStateType::BehaviorType::Constrained_InWater:
		case HumanNpcStateType::BehaviorType::Constrained_Swimming_Style1:
		case HumanNpcStateType::BehaviorType::Constrained_Swimming_Style2:
		{
			if (isFree)
			{
				// Transition
				TransitionHumanBehaviorToFree(npc, currentSimulationTime);

				break;
			}

			// Check if moved to air

			// It's in air if both in air
			if (mParticles.GetAnyWaterness(primaryParticleState.ParticleIndex) < 0.25f
				&& mParticles.GetAnyWaterness(secondaryParticleState.ParticleIndex) < 0.25f)
			{
				// Transition

				humanState.TransitionToState(HumanNpcStateType::BehaviorType::Constrained_Aerial, currentSimulationTime);

#ifdef BARYLAB_PROBING
				if (npc.Id == mCurrentlySelectedNpc)
				{
					mGameEventHandler->OnHumanNpcBehaviorChanged("Constrained_Aerial");
				}
#endif

				break;
			}

			// Advance state machine

			if (humanState.CurrentBehavior == HumanNpcStateType::BehaviorType::Constrained_InWater)
			{
				// Progress to swimming after a while here

				float constexpr ToSwimmingConvergenceRate = 0.01f;
				humanState.CurrentBehaviorState.Constrained_InWater.ProgressToSwimming +=
					(1.0f - humanState.CurrentBehaviorState.Constrained_InWater.ProgressToSwimming)
					* ToSwimmingConvergenceRate;

#ifdef BARYLAB_PROBING
				publishStateQuantity = std::make_tuple("ProgressToSwimming", std::to_string(humanState.CurrentBehaviorState.Constrained_InWater.ProgressToSwimming));
#endif

				if (IsAtTarget(humanState.CurrentBehaviorState.Constrained_InWater.ProgressToSwimming, 0.98f))
				{
					// Transition

					StateType::KindSpecificStateType::HumanNpcStateType::BehaviorType const swimStyle = (humanState.CurrentFaceOrientation != 0.0f)
						? HumanNpcStateType::BehaviorType::Constrained_Swimming_Style1
						: HumanNpcStateType::BehaviorType::Constrained_Swimming_Style2;

					humanState.TransitionToState(swimStyle, currentSimulationTime);

#ifdef BARYLAB_PROBING
					if (npc.Id == mCurrentlySelectedNpc)
					{
						mGameEventHandler->OnHumanNpcBehaviorChanged("Constrained_Swimming");
					}
#endif

					break;
				}
			}

			break;
		}

		case HumanNpcStateType::BehaviorType::Constrained_Electrified:
		{
			if (isFree)
			{
				// Transition
				TransitionHumanBehaviorToFree(npc, currentSimulationTime);

				break;
			}

			// Advance towards leaving

			float toLeavingIncrement;
			if (!IsElectrified(npc, homeShip))
			{
				toLeavingIncrement = 1.0f;
			}
			else
			{
				toLeavingIncrement = -1.0f;
			}

			humanState.CurrentBehaviorState.Constrained_Electrified.ProgressToLeaving = std::max(
				humanState.CurrentBehaviorState.Constrained_Electrified.ProgressToLeaving + toLeavingIncrement,
				0.0f);

			float constexpr ToLeavingTarget = 8.0f;

#ifdef BARYLAB_PROBING
			publishStateQuantity = std::make_tuple("ProgressToLeaving", std::to_string(humanState.CurrentBehaviorState.Constrained_Electrified.ProgressToLeaving / ToLeavingTarget));
#endif

			if (humanState.CurrentBehaviorState.Constrained_Electrified.ProgressToLeaving >= ToLeavingTarget)
			{
				// Transition

				humanState.TransitionToState(HumanNpcStateType::BehaviorType::Constrained_KnockedOut, currentSimulationTime);

#ifdef BARYLAB_PROBING
				if (npc.Id == mCurrentlySelectedNpc)
				{
					mGameEventHandler->OnHumanNpcBehaviorChanged("Constrained_KnockedOut");
				}
#endif

				break;
			}

			// Maintain state

			humanState.EquilibriumTorque = 1.0f;

			break;
		}

		case HumanNpcStateType::BehaviorType::Free_Aerial:
		{
			if (!isFree)
			{
				// Transition

				humanState.TransitionToState(HumanNpcStateType::BehaviorType::Constrained_KnockedOut, currentSimulationTime);

#ifdef BARYLAB_PROBING
				if (npc.Id == mCurrentlySelectedNpc)
				{
					mGameEventHandler->OnHumanNpcBehaviorChanged("Constrained_KnockedOut");
				}
#endif

				break;
			}

			// Check if moved to water

			// It's in water if at least one in water
			if (mParticles.GetAnyWaterness(primaryParticleState.ParticleIndex) > 0.0f
				|| mParticles.GetAnyWaterness(secondaryParticleState.ParticleIndex) > 0.0f)
			{
				// Transition

				humanState.TransitionToState(HumanNpcStateType::BehaviorType::Free_InWater, currentSimulationTime);

#ifdef BARYLAB_PROBING
				if (npc.Id == mCurrentlySelectedNpc)
				{
					mGameEventHandler->OnHumanNpcBehaviorChanged("Free_InWater");
				}
#endif

				break;
			}

			// Progress to knocked out (when still)

			float knockedOutTarget;

			auto const velocityMagnitude =
				(primaryParticleState.GetApplicableVelocity(mParticles).length() + secondaryParticleState.GetApplicableVelocity(mParticles).length())
				/ 2.0f;

			if (velocityMagnitude < 0.1f)
			{
				knockedOutTarget = 1.0f;
			}
			else
			{
				knockedOutTarget = 0.0f;
			}

			float constexpr ToKnockedOutConvergenceRate = 0.2f;

			humanState.CurrentBehaviorState.Free_Aerial.ProgressToKnockedOut +=
				(knockedOutTarget - humanState.CurrentBehaviorState.Free_Aerial.ProgressToKnockedOut)
				* ToKnockedOutConvergenceRate;

#ifdef BARYLAB_PROBING
			publishStateQuantity = std::make_tuple("ProgressToKnockedOut", std::to_string(humanState.CurrentBehaviorState.Free_Aerial.ProgressToKnockedOut));
#endif

			if (IsAtTarget(humanState.CurrentBehaviorState.Free_Aerial.ProgressToKnockedOut, 1.0f))
			{
				// Transition

				humanState.TransitionToState(HumanNpcStateType::BehaviorType::Free_KnockedOut, currentSimulationTime);

#ifdef BARYLAB_PROBING
				if (npc.Id == mCurrentlySelectedNpc)
				{
					mGameEventHandler->OnHumanNpcBehaviorChanged("Free_KnockedOut");
				}
#endif

				break;
			}

			break;
		}

		case HumanNpcStateType::BehaviorType::Free_KnockedOut:
		{
			if (!isFree)
			{
				// Transition

				humanState.TransitionToState(HumanNpcStateType::BehaviorType::Constrained_KnockedOut, currentSimulationTime);

#ifdef BARYLAB_PROBING
				if (npc.Id == mCurrentlySelectedNpc)
				{
					mGameEventHandler->OnHumanNpcBehaviorChanged("Constrained_KnockedOut");
				}
#endif

				break;
			}

			// Progress to in-water (when at least one particle in water)

			float inWaterTarget;

			if (mParticles.GetAnyWaterness(primaryParticleState.ParticleIndex) > 0.0f
				|| mParticles.GetAnyWaterness(secondaryParticleState.ParticleIndex) > 0.0f)
			{
				inWaterTarget = 1.0f;
			}
			else
			{
				inWaterTarget = 0.0f;
			}

			float constexpr ToInWaterConvergenceRate = 0.2f;

			humanState.CurrentBehaviorState.Free_KnockedOut.ProgressToInWater +=
				(inWaterTarget - humanState.CurrentBehaviorState.Free_KnockedOut.ProgressToInWater)
				* ToInWaterConvergenceRate;

#ifdef BARYLAB_PROBING
			publishStateQuantity = std::make_tuple("ProgressToInWater", std::to_string(humanState.CurrentBehaviorState.Free_KnockedOut.ProgressToInWater));
#endif

			if (IsAtTarget(humanState.CurrentBehaviorState.Free_KnockedOut.ProgressToInWater, 1.0f))
			{
				// Transition

				humanState.TransitionToState(HumanNpcStateType::BehaviorType::Free_InWater, currentSimulationTime);

#ifdef BARYLAB_PROBING
				if (npc.Id == mCurrentlySelectedNpc)
				{
					mGameEventHandler->OnHumanNpcBehaviorChanged("Free_InWater");
				}
#endif

				break;
			}

			// Progress to aerial (when moving and out-of-water)

			float aerialTarget;

			auto const velocityMagnitude =
				(primaryParticleState.GetApplicableVelocity(mParticles).length() + secondaryParticleState.GetApplicableVelocity(mParticles).length())
				/ 2.0f;

			if (velocityMagnitude > 0.5f)
			{
				aerialTarget = 1.0f;
			}
			else
			{
				aerialTarget = 0.0f;
			}

			float constexpr ToAerialConvergenceRate = 0.2f;

			humanState.CurrentBehaviorState.Free_KnockedOut.ProgressToAerial +=
				(aerialTarget - humanState.CurrentBehaviorState.Free_KnockedOut.ProgressToAerial)
				* ToAerialConvergenceRate;

#ifdef BARYLAB_PROBING
			publishStateQuantity = std::make_tuple("ProgressToAerial", std::to_string(humanState.CurrentBehaviorState.Free_KnockedOut.ProgressToAerial));
#endif

			if (IsAtTarget(humanState.CurrentBehaviorState.Free_KnockedOut.ProgressToAerial, 1.0f))
			{
				// Transition

				humanState.TransitionToState(HumanNpcStateType::BehaviorType::Free_Aerial, currentSimulationTime);

#ifdef BARYLAB_PROBING
				if (npc.Id == mCurrentlySelectedNpc)
				{
					mGameEventHandler->OnHumanNpcBehaviorChanged("Free_Aerial");
				}
#endif

				break;
			}

			break;
		}

		case HumanNpcStateType::BehaviorType::Free_InWater:
		case HumanNpcStateType::BehaviorType::Free_Swimming_Style1:
		case HumanNpcStateType::BehaviorType::Free_Swimming_Style2:
		case HumanNpcStateType::BehaviorType::Free_Swimming_Style3:
		{
			if (!isFree)
			{
				// Transition

				humanState.TransitionToState(HumanNpcStateType::BehaviorType::Constrained_KnockedOut, currentSimulationTime);

#ifdef BARYLAB_PROBING
				if (npc.Id == mCurrentlySelectedNpc)
				{
					mGameEventHandler->OnHumanNpcBehaviorChanged("Constrained_KnockedOut");
				}
#endif

				break;
			}

			// Check if moved to air

			// It's in air if both in air
			if (mParticles.GetAnyWaterness(primaryParticleState.ParticleIndex) == 0.0f
				&& mParticles.GetAnyWaterness(secondaryParticleState.ParticleIndex) == 0.0f)
			{
				// Transition

				humanState.TransitionToState(HumanNpcStateType::BehaviorType::Free_Aerial, currentSimulationTime);

#ifdef BARYLAB_PROBING
				if (npc.Id == mCurrentlySelectedNpc)
				{
					mGameEventHandler->OnHumanNpcBehaviorChanged("Free_Aerial");
				}
#endif

				break;
			}

			// Advance state machine

			vec2f const feetVelocity = mParticles.GetVelocity(primaryParticleState.ParticleIndex);
			vec2f const headVelocity = mParticles.GetVelocity(secondaryParticleState.ParticleIndex);
			float constexpr MaxRelativeVelocityForSwimming = 2.0f;
			float const isRelativeVelocityOk = Step((headVelocity - feetVelocity).length(), MaxRelativeVelocityForSwimming); // Also rotation
			float constexpr MaxVelocityForSwimming = 4.0f;
			float const isAbsoluteVelocityOk = Step(feetVelocity.length(), MaxVelocityForSwimming) * Step(headVelocity.length(), MaxVelocityForSwimming);

			if (humanState.CurrentBehavior == HumanNpcStateType::BehaviorType::Free_InWater)
			{
				// Progress to swimming if velocities ok, and head above feet

				auto const & headPosition = mParticles.GetPosition(secondaryParticleState.ParticleIndex);
				auto const & feetPosition = mParticles.GetPosition(primaryParticleState.ParticleIndex);

				float const targetSwim =
					isRelativeVelocityOk
					* isAbsoluteVelocityOk
					* Step(feetPosition.y, headPosition.y);

				float constexpr ToSwimmingConvergenceRate = 0.12f;
				humanState.CurrentBehaviorState.Free_InWater.ProgressToSwimming += (targetSwim - humanState.CurrentBehaviorState.Free_InWater.ProgressToSwimming) * ToSwimmingConvergenceRate;

#ifdef BARYLAB_PROBING
				publishStateQuantity = std::make_tuple("ProgressToSwimming", std::to_string(humanState.CurrentBehaviorState.Free_InWater.ProgressToSwimming));
#endif

				if (IsAtTarget(humanState.CurrentBehaviorState.Free_InWater.ProgressToSwimming, 0.9f)) // We're content with "almost"
				{
					// Transition

					StateType::KindSpecificStateType::HumanNpcStateType::BehaviorType swimStyle;
					switch (GameRandomEngine::GetInstance().Choose(5))
					{
						case 0:
						{
							swimStyle = HumanNpcStateType::BehaviorType::Free_Swimming_Style2;
							break;
						}

						case 1:
						{
							swimStyle = HumanNpcStateType::BehaviorType::Free_Swimming_Style3;
							break;
						}

						default:
						{
							swimStyle = HumanNpcStateType::BehaviorType::Free_Swimming_Style1;
							break;
						}
					}

					humanState.TransitionToState(swimStyle, currentSimulationTime);

					// Face: FvB/0
					humanState.CurrentFaceOrientation = (GameRandomEngine::GetInstance().Choose(3) == 0) ? -1.0f : 1.0f;
					humanState.CurrentFaceDirectionX = 0.0f;

#ifdef BARYLAB_PROBING
					if (npc.Id == mCurrentlySelectedNpc)
					{
						mGameEventHandler->OnHumanNpcBehaviorChanged("Free_Swimming");
					}
#endif

					break;
				}
			}
			else
			{
				// Leave swimming if velocities not ok

				float const targetLeavingSwimming = 1.0f -
					isRelativeVelocityOk
					* isAbsoluteVelocityOk;

				float constexpr ToLeavingSwimmingConvergenceRate = 0.12f;
				humanState.CurrentBehaviorState.Free_Swimming.ProgressToLeavingSwimming += (targetLeavingSwimming - humanState.CurrentBehaviorState.Free_Swimming.ProgressToLeavingSwimming) * ToLeavingSwimmingConvergenceRate;

#ifdef BARYLAB_PROBING
				publishStateQuantity = std::make_tuple("ProgressToLeavingSwimming", std::to_string(humanState.CurrentBehaviorState.Free_Swimming.ProgressToLeavingSwimming));
#endif

				if (IsAtTarget(humanState.CurrentBehaviorState.Free_Swimming.ProgressToLeavingSwimming, 1.0f))
				{
					// Transition

					humanState.TransitionToState(HumanNpcStateType::BehaviorType::Free_InWater, currentSimulationTime);

#ifdef BARYLAB_PROBING
					if (npc.Id == mCurrentlySelectedNpc)
					{
						mGameEventHandler->OnHumanNpcBehaviorChanged("Free_InWater");
					}
#endif

					break;
				}
			}

			// We are staying at this state...
			// ...update this state then

			// See if time to emit a bubble

			float & nextBubbleSimulationTime = (humanState.CurrentBehavior == HumanNpcStateType::BehaviorType::Free_InWater)
				? humanState.CurrentBehaviorState.Free_InWater.NextBubbleEmissionSimulationTimestamp
				: humanState.CurrentBehaviorState.Free_Swimming.NextBubbleEmissionSimulationTimestamp;

			if (currentSimulationTime > nextBubbleSimulationTime)
			{
				homeShip.SpawnAirBubble(
					mParticles.GetPosition(secondaryParticleState.ParticleIndex),
					GameParameters::NpcAirBubbleFinalScale,
					GameParameters::HumanNpcTemperature,
					currentSimulationTime,
					npc.CurrentPlaneId,
					gameParameters);

				// Calculate next emission timestamp
				nextBubbleSimulationTime =
					currentSimulationTime
					+ 1.0f // Min delay
					+ std::min(GameRandomEngine::GetInstance().GenerateExponentialReal(0.5f), 2.0f);
			}

			break;
		}

		case HumanNpcStateType::BehaviorType::ConstrainedOrFree_Smashed:
		{
			// Advance towards leaving

			humanState.CurrentBehaviorState.ConstrainedOrFree_Smashed.ProgressToLeaving +=
				(1.0f - humanState.CurrentBehaviorState.ConstrainedOrFree_Smashed.ProgressToLeaving)
				* 0.04f;

			if (IsAtTarget(humanState.CurrentBehaviorState.ConstrainedOrFree_Smashed.ProgressToLeaving, 1.0f))
			{
				// Transition

				if (isFree)
				{
					humanState.TransitionToState(HumanNpcStateType::BehaviorType::Free_KnockedOut, currentSimulationTime);

#ifdef BARYLAB_PROBING
					if (npc.Id == mCurrentlySelectedNpc)
					{
						mGameEventHandler->OnHumanNpcBehaviorChanged("Free_KnockedOut");
					}
#endif
				}
				else
				{
					humanState.TransitionToState(HumanNpcStateType::BehaviorType::Constrained_KnockedOut, currentSimulationTime);

#ifdef BARYLAB_PROBING
					if (npc.Id == mCurrentlySelectedNpc)
					{
						mGameEventHandler->OnHumanNpcBehaviorChanged("Constrained_KnockedOut");
					}
#endif
				}

				break;
			}

			break;
		}

		case HumanNpcStateType::BehaviorType::BeingRemoved:
		{
			float const elapsed = currentSimulationTime - humanState.CurrentStateTransitionSimulationTimestamp;

			auto & behaviorState = humanState.CurrentBehaviorState.BeingRemoved;

			switch (behaviorState.CurrentState)
			{
				case HumanNpcStateType::BehaviorStateType::BeingRemovedStateType::StateType::Init:
				{
					auto const & headPosition = mParticles.GetPosition(secondaryParticleState.ParticleIndex);
					auto const & feetPosition = mParticles.GetPosition(primaryParticleState.ParticleIndex);
					vec2f const humanVector = headPosition - feetPosition;
					float const actualAngleCW = humanVector.angleCw(vec2f(0.0f, 1.0f));

					// Initialize angle and duration
					behaviorState.StartUprightAngle = actualAngleCW;
					behaviorState.TotalUprightDuration = std::fabsf(actualAngleCW) / (Pi<float> / 2.0f) * 2.0f; // 2s for PI/2

					// Transition
					behaviorState.CurrentState = HumanNpcStateType::BehaviorStateType::BeingRemovedStateType::StateType::GettingUpright;
					behaviorState.NextTransitionTimestamp = std::max(behaviorState.TotalUprightDuration, HumanRemovalLevitationDuration);

					break;
				}

				case HumanNpcStateType::BehaviorStateType::BeingRemovedStateType::StateType::GettingUpright:
				{
					//
					// Check whether we're done
					//

					if (elapsed >= behaviorState.NextTransitionTimestamp)
					{
						behaviorState.CurrentState = HumanNpcStateType::BehaviorStateType::BeingRemovedStateType::StateType::PreRotation;
						behaviorState.NextTransitionTimestamp += HumanRemovalPreRotationDuration;
						break;
					}

					//
					// Levitation + getting upright
					//

					auto & headPosition = mParticles.GetPosition(secondaryParticleState.ParticleIndex);
					auto & feetPosition = mParticles.GetPosition(primaryParticleState.ParticleIndex);

					//
					// Get upright
					//

					vec2f const humanVector = headPosition - feetPosition;

					// Calculate desired angle
					float const desiredAngle = (1.0f - SmoothStep(0.0f, behaviorState.TotalUprightDuration, elapsed)) * (behaviorState.StartUprightAngle);

					// Get upright
					headPosition = feetPosition + vec2f(0.0f, 1.0f).rotate(desiredAngle) * humanVector.length();

					//
					// Levitate
					//

					float const levitationDepth = std::max(1.0f - (elapsed * elapsed) / (HumanRemovalLevitationDuration * HumanRemovalLevitationDuration), 0.0f);
					vec2f const upTraslation = vec2f(0.0f, GameParameters::SimulationStepTimeDuration<float> / HumanRemovalLevitationDuration * levitationDepth);
					feetPosition += upTraslation;
					headPosition += upTraslation;

					break;
				}

				case HumanNpcStateType::BehaviorStateType::BeingRemovedStateType::StateType::PreRotation:
				{
					//
					// Check whether we're done
					//

					if (elapsed >= behaviorState.NextTransitionTimestamp)
					{
						//
						// Save current limbs
						//

						if (humanState.CurrentFaceOrientation == 0.0f)
						{
							// L/R
							behaviorState.WorkingLimbLRAngles = humanState.AnimationState.LimbAngles;
							behaviorState.WorkingLimbFBAngles = LimbVector{0.0f, 0.0f, 0.0f, 0.0f };
						}
						else
						{
							// F/B
							behaviorState.WorkingLimbFBAngles = humanState.AnimationState.LimbAngles;
							behaviorState.WorkingLimbLRAngles = LimbVector{ 0.0f, 0.0f, 0.0f, 0.0f };
						}

						//
						// Go frontal
						//

						humanState.CurrentFaceDirectionX = 0.0f;
						humanState.CurrentFaceOrientation = 1.0f;

						// Transition
						behaviorState.CurrentState = HumanNpcStateType::BehaviorStateType::BeingRemovedStateType::StateType::Rotating;
						behaviorState.NextTransitionTimestamp += HumanRemovalRotationDuration;

						break;
					}

					break;
				}

				case HumanNpcStateType::BehaviorStateType::BeingRemovedStateType::StateType::Rotating:
				{
					//
					// Check whether we're done
					//

					if (elapsed >= behaviorState.NextTransitionTimestamp)
					{
						// Deferred removal
						assert(std::find(mDeferredRemovalNpcs.cbegin(), mDeferredRemovalNpcs.cend(), npc.Id) == mDeferredRemovalNpcs.cend());
						mDeferredRemovalNpcs.emplace_back(npc.Id);

						break;
					}

					//
					// Rotate
					//

					// TODOHERE

					break;
				}
			}

			// TODOOLD
			////float const elapsed = currentSimulationTime - humanState.CurrentStateTransitionSimulationTimestamp;

			////if (elapsed < HumanRemovalDelay2)
			////{
			////	// Traslate up

			////	float depth = 1.0f - (elapsed * elapsed) / (HumanRemovalDelay2 * HumanRemovalDelay2);
			////	vec2f const upTraslation = vec2f(0.0f, GameParameters::SimulationStepTimeDuration<float> / HumanRemovalDelay2 * depth);
			////	mParticles.SetPosition(
			////		primaryParticleState.ParticleIndex,
			////		mParticles.GetPosition(primaryParticleState.ParticleIndex) + upTraslation);
			////	mParticles.SetPosition(
			////		secondaryParticleState.ParticleIndex,
			////		mParticles.GetPosition(secondaryParticleState.ParticleIndex) + upTraslation);
			////}

			////if (elapsed >= HumanRemovalDelay1 && elapsed < HumanRemovalDelay2)
			////{
			////	// Face forward

			////	humanState.CurrentFaceDirectionX = 0.0f;
			////	humanState.CurrentFaceOrientation = 1.0f;
			////}

			//////if (elapsed >= HumanRemovalDelay1)
			////{
			////	if (elapsed < HumanRemovalDuration)
			////	{
			////		// Get upright

			////		auto & headPosition = mParticles.GetPosition(secondaryParticleState.ParticleIndex);
			////		auto const & feetPosition = mParticles.GetPosition(primaryParticleState.ParticleIndex);
			////		vec2f const humanVector = headPosition - feetPosition;
			////		float const actualAngleCW = humanVector.angleCw(vec2f(0.0f, 1.0f));
			////		headPosition = feetPosition + humanVector.rotate(-actualAngleCW * 0.042f);

			////		// Rotate

			////		if (elapsed > humanState.CurrentBehaviorState.BeingRemoved.NextRotationSimulationTimestamp)
			////		{
			////			if (humanState.CurrentFaceOrientation == 1.0f)
			////			{
			////				humanState.CurrentFaceDirectionX = -1.0f;
			////				humanState.CurrentFaceOrientation = 0.0f;
			////			}
			////			else if (humanState.CurrentFaceDirectionX == -1.0f)
			////			{
			////				humanState.CurrentFaceDirectionX = 0.0f;
			////				humanState.CurrentFaceOrientation = -1.0f;
			////			}
			////			else if (humanState.CurrentFaceOrientation == -1.0f)
			////			{
			////				humanState.CurrentFaceDirectionX = 1.0f;
			////				humanState.CurrentFaceOrientation = 0.0f;
			////			}
			////			else
			////			{
			////				assert(humanState.CurrentFaceDirectionX == 1.0f);
			////				humanState.CurrentFaceDirectionX = 0.0f;
			////				humanState.CurrentFaceOrientation = 1.0f;
			////			}

			////			// Calculate next rotation timestamp
			////			float const progress = Clamp((elapsed - HumanRemovalDelay2) / (HumanRemovalDuration * 0.9f - HumanRemovalDelay2), 0.0f, 1.0f);
			////			humanState.CurrentBehaviorState.BeingRemoved.NextRotationSimulationTimestamp +=
			////				+ (1.0f - std::sqrtf(progress)) * 0.35f;
			////		}

			////	}
			////	else
			////	{
			////		// We're done

			////		// Deferred removal
			////		assert(std::find(mDeferredRemovalNpcs.cbegin(), mDeferredRemovalNpcs.cend(), npc.Id) == mDeferredRemovalNpcs.cend());
			////		mDeferredRemovalNpcs.emplace_back(npc.Id);
			////	}
			////}


			break;
		}
	}

#ifdef BARYLAB_PROBING
	if (npc.Id == mCurrentlySelectedNpc)
	{
		mGameEventHandler->OnHumanNpcStateQuantityChanged(publishStateQuantity);
	}
#endif
}

bool Npcs::CheckAndMaintainHumanEquilibrium(
	ElementIndex primaryParticleIndex,
	ElementIndex secondaryParticleIndex,
	StateType::KindSpecificStateType::HumanNpcStateType & humanState,
	bool doMaintainEquilibrium,
	NpcParticles & particles,
	GameParameters const & /*gameParameters*/)
{
	//
	// Make sure we are not falling out of equilibrium
	//

	vec2f const humanVector = particles.GetPosition(secondaryParticleIndex) - particles.GetPosition(primaryParticleIndex);
	vec2f const humanDir = humanVector.normalise_approx();

	//
	// Static angle necessary condition: human vector outside of sector around vertical
	//
	// We use y component of normalized human vector (i.e. cos(angle with vertical), +1.0 when fully vertical, < 1.0f when less vertical),
	// and we're out if y < cos(MaxAngle)
	//

	float constexpr CosMaxStaticAngleForEquilibrium = 0.84f; //  cos(Pi / 5.5)

	if (humanDir.y < CosMaxStaticAngleForEquilibrium)
	{
		//
		// Rotational velocity necessary condition: non-negligible rotational velocity away from vertical
		//
		// We're out (diverging from vertical) if RelVel component along perpendicular to humanDir (i.e. radialVelocity)
		// is < or > 0 depending on whether head is to the left or to the right of ideal head,
		// i.e. if RelVel dot perp(humanDir) * (IH.x - H.x) > 0
		//

		vec2f const relativeVelocity = (particles.GetVelocity(secondaryParticleIndex) - particles.GetVelocity(primaryParticleIndex));
		float const radialVelocity = relativeVelocity.dot(humanDir.to_perpendicular());

		float const maxRadialVelocityFactor = (humanState.CurrentBehavior == StateType::KindSpecificStateType::HumanNpcStateType::BehaviorType::Constrained_Rising)
			? 1.056f // chord/dt of a human traveling a 0.01 angle
			: 0.0f;

		if (radialVelocity * (-humanVector.x) > maxRadialVelocityFactor)
		{
			LogNpcDebug("Losing equilibrium because: humanDir.y=", humanDir.y, " < ", CosMaxStaticAngleForEquilibrium, " && ",
				"radialVelocity * (-humanVector.x)=", radialVelocity * (-humanVector.x), " > ", maxRadialVelocityFactor);

			return false;
		}
	}

	//
	// We are in equilibrium
	//

	//
	// Maintain equilibrium
	//

	if (doMaintainEquilibrium)
	{
		humanState.EquilibriumTorque = 1.0f;
	}

	return true;
}

void Npcs::RunWalkingHumanStateMachine(
	StateType::KindSpecificStateType::HumanNpcStateType & humanState,
	StateType::NpcParticleStateType const & primaryParticleState,
	Ship const & /*homeShip*/, // Will come useful when we'll *plan* the walk
	GameParameters const & /*gameParameters*/)
{
	assert(primaryParticleState.ConstrainedState.has_value());
	assert(humanState.CurrentBehavior == StateType::KindSpecificStateType::HumanNpcStateType::BehaviorType::Constrained_Walking);

	auto & walkingState = humanState.CurrentBehaviorState.Constrained_Walking;

	// 1. Check condition for growing decision to flip: actual (relative) velocity opposite of walking direction,
	// or too small

	if (walkingState.CurrentWalkMagnitude != 0.0f)
	{
		float const relativeVelocityAgreement = primaryParticleState.ConstrainedState->MeshRelativeVelocity.dot(
			vec2f(
				humanState.CurrentFaceDirectionX * CalculateActualHumanWalkingAbsoluteSpeed(humanState),
				0.0f));

		LogNpcDebug("        relativeVelocityAgreement=", relativeVelocityAgreement);

		float constexpr MinRelativeVelocityAgreementToAcceptWalk = 0.025f;
		if (relativeVelocityAgreement < MinRelativeVelocityAgreementToAcceptWalk)
		{
			// Grow decision to flip

			LogNpcDebug("        Growing impatient because of relative velocity not catching up with walk velocity");

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

		LogNpcDebug("        Reached flip decision");

		FlipHumanWalk(humanState, StrongTypedTrue<_DoImmediate>);
	}

	//
	// 4. Advance walking magnitude towards full walk
	//

	float const walkMagnitudeConvergenceRate = 0.10f + std::min(humanState.ResultantPanicLevel, 1.0f) * 0.08f;
	walkingState.CurrentWalkMagnitude += (1.0f - walkingState.CurrentWalkMagnitude) * walkMagnitudeConvergenceRate;

	LogNpcDebug("        currentWalkMagnitude: ", walkingState.CurrentWalkMagnitude);
}

void Npcs::OnHumanImpact(
	StateType & npc,
	int npcParticleOrdinal,
	vec2f const & normalResponse,
	vec2f const & bounceEdgeNormal,  // Pointing outside of triangle
	float currentSimulationTime) const
{
	assert(npc.Kind == NpcKindType::Human);

	auto & humanState = npc.KindSpecificState.HumanNpcState;

	switch (humanState.CurrentBehavior)
	{
		case StateType::KindSpecificStateType::HumanNpcStateType::BehaviorType::Constrained_Rising:
		{
			if (npcParticleOrdinal == 1 && normalResponse.length() > 0.4f)
			{
				// Hit head while rising

				LogNpcDebug("OnHumanImpact: Hit head while rising - going to KnockedOut");

				humanState.TransitionToState(StateType::KindSpecificStateType::HumanNpcStateType::BehaviorType::Constrained_KnockedOut, currentSimulationTime);
			}

			break;
		}

		case StateType::KindSpecificStateType::HumanNpcStateType::BehaviorType::Constrained_Equilibrium:
		{
			if (npcParticleOrdinal == 1 && normalResponse.length() > 1.5f)
			{
				// Hit head hard while in equilibrium

				LogNpcDebug("OnHumanImpact: Hit head hard while in equilibrium - going to KnockedOut");

				humanState.TransitionToState(StateType::KindSpecificStateType::HumanNpcStateType::BehaviorType::Constrained_KnockedOut, currentSimulationTime);
			}

			break;
		}

		case StateType::KindSpecificStateType::HumanNpcStateType::BehaviorType::Constrained_Walking:
		{
			LogNpcDebug("OnHumanImpact: bounceSlope=", bounceEdgeNormal.dot(vec2f(humanState.CurrentFaceDirectionX, 0.0f)));

			// Check alignment of impact with walking direction; if hit => flip
			// Note: might also want to check *magnitude* of hit
			float const bounceSlope = bounceEdgeNormal.dot(vec2f(humanState.CurrentFaceDirectionX, 0.0f)); // 1.0 when hit wall perpendicular
			if ( ((npcParticleOrdinal == 0 && bounceSlope > 0.85f) || (npcParticleOrdinal == 1 && bounceSlope > 0.50f))
				&& humanState.CurrentBehaviorState.Constrained_Walking.CurrentWalkMagnitude != 0.0f)
			{
				LogNpcDebug("OnHumanImpact: FLIP!");

				// Flip now
				FlipHumanWalk(humanState, StrongTypedTrue<_DoImmediate>);
			}

			break;
		}

		default:
		{
			break;
		}
	}
}

std::optional<float> Npcs::DecideWalkingDirection(
	StateType const & npc,
	StateType::NpcParticleStateType const & primaryParticleState,
	Ship const & homeShip)
{
	FixedSizeVector<float, 2> candidates;

	assert(primaryParticleState.ConstrainedState.has_value());
	if (!primaryParticleState.ConstrainedState->CurrentVirtualFloor.has_value())
	{
		// Not much we can do, just choose randomly
		candidates.emplace_back(1.0f);
		candidates.emplace_back(-1.0f);
	}
	else
	{
		if (CanWalkInDirection(npc, *primaryParticleState.ConstrainedState->CurrentVirtualFloor, 1.0f, homeShip))
		{
			candidates.emplace_back(1.0f);
		}

		if (CanWalkInDirection(npc, *primaryParticleState.ConstrainedState->CurrentVirtualFloor, -1.0f, homeShip))
		{
			candidates.emplace_back(-1.0f);
		}
	}

	if (candidates.empty())
	{
		return std::nullopt;
	}

	return candidates[GameRandomEngine::GetInstance().Choose(candidates.size())];
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

void Npcs::TransitionHumanBehaviorToFree(
	StateType & npc,
	float currentSimulationTime)
{
	assert(npc.Kind == NpcKindType::Human);

	assert(npc.ParticleMesh.Particles.size() == 2);
	auto const & headPosition = mParticles.GetPosition(npc.ParticleMesh.Particles[1].ParticleIndex);
	auto const & feetPosition = mParticles.GetPosition(npc.ParticleMesh.Particles[0].ParticleIndex);

	// It's in water if both in water
	if (mParentWorld.GetOceanSurface().GetDepth(headPosition) > 0.0f
		&& mParentWorld.GetOceanSurface().GetDepth(feetPosition) > 0.0f)
	{
		npc.KindSpecificState.HumanNpcState.TransitionToState(StateType::KindSpecificStateType::HumanNpcStateType::BehaviorType::Free_InWater, currentSimulationTime);

#ifdef BARYLAB_PROBING
		if (npc.Id == mCurrentlySelectedNpc)
		{
			mGameEventHandler->OnHumanNpcBehaviorChanged("Free_InWater");
		}
#endif
	}
	else
	{
		npc.KindSpecificState.HumanNpcState.TransitionToState(StateType::KindSpecificStateType::HumanNpcStateType::BehaviorType::Free_Aerial, currentSimulationTime);

#ifdef BARYLAB_PROBING
		if (npc.Id == mCurrentlySelectedNpc)
		{
			mGameEventHandler->OnHumanNpcBehaviorChanged("Free_Aerial");
		}
#endif
	}
}

}