/***************************************************************************************
 * Original Author:		Gabriele Giuseppini
 * Created:				2023-10-20
 * Copyright:			Gabriele Giuseppini  (https://github.com/GabrieleGiuseppini)
 ***************************************************************************************/
#include "Npcs.h"

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
		mEventDispatcher.OnHumanNpcBehaviorChanged("Constrained_KnockedOut");
		return StateType::HumanNpcStateType(StateType::HumanNpcStateType::BehaviorType::Constrained_KnockedOut, 0.0f, 0.0f);
	}
}

void Npcs::UpdateHuman(
	StateType::HumanNpcStateType & humanState,
	StateType::NpcParticleStateType const & primaryParticleState,
	StateType::NpcParticleStateType const & secondaryParticleState,
	Mesh const & mesh)
{
	// TODOHERE
	(void)mesh;

	std::optional<std::tuple<std::string, std::string>> publishStateQuantity;

	bool const isFree = !primaryParticleState.ConstrainedState.has_value() && !secondaryParticleState.ConstrainedState.has_value();

	switch (humanState.CurrentBehavior)
	{
		case StateType::HumanNpcStateType::BehaviorType::Constrained_KnockedOut:
		{
			float const RisingConvergenceRate = 0.05f;
			float constexpr MaxRelativeVelocity = 0.1f;

			if (isFree)
			{
				// Transition

				humanState.CurrentBehavior = StateType::HumanNpcStateType::BehaviorType::Free_KnockedOut;
				mEventDispatcher.OnHumanNpcBehaviorChanged("Free_KnockedOut");

				break;
			}

			// Advance

			humanState.CurrentStateValue += (humanState.TargetStateValue - humanState.CurrentStateValue) * RisingConvergenceRate;
			publishStateQuantity = std::make_tuple("CurrentStateValue", std::to_string(humanState.CurrentStateValue));

			if (IsAtTarget(humanState.CurrentStateValue, 1.0f))
			{
				// Transition

				humanState.CurrentBehavior = StateType::HumanNpcStateType::BehaviorType::Constrained_Rising;
				mEventDispatcher.OnHumanNpcBehaviorChanged("Constrained_Rising");

				humanState.CurrentStateValue = 0.0f;
				humanState.TargetStateValue = 0.0f;

				break;
			}

			// Check conditions for this state

			float stateCondition = 0.0f;

			if (primaryParticleState.ConstrainedState.has_value()
				&& IsOnEdge(primaryParticleState.ConstrainedState->CurrentTriangleBarycentricCoords)
				&& primaryParticleState.ConstrainedState->MeshRelativeVelocity.length() < MaxRelativeVelocity
				&& secondaryParticleState.ConstrainedState.has_value()
				&& IsOnEdge(secondaryParticleState.ConstrainedState->CurrentTriangleBarycentricCoords)
				&& secondaryParticleState.ConstrainedState->MeshRelativeVelocity.length() < MaxRelativeVelocity)
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

				humanState.CurrentBehavior = StateType::HumanNpcStateType::BehaviorType::Free_KnockedOut;
				mEventDispatcher.OnHumanNpcBehaviorChanged("Free_KnockedOut");

				break;
			}

			// Check conditions to stay
			// TODOHERE

			break;
		}

		case StateType::HumanNpcStateType::BehaviorType::Constrained_Equilibrium:
		{
			// TODOHERE
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