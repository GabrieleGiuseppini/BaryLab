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

	float CalculateVerticalAlignment(ElementIndex primaryParticleIndex, ElementIndex secondaryParticleIndex, NpcParticles const & particles)
	{
		vec2f const humanVerticalVector = particles.GetPosition(primaryParticleIndex) - particles.GetPosition(secondaryParticleIndex);
		return humanVerticalVector.normalise().dot(LabParameters::GravityDir);
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
	Mesh const & mesh)
{
	float const ToRisingConvergenceRate = 0.05f;
	//float constexpr MaxRelativeVelocityForEquilibrium = 0.1f;
	float constexpr MaxRelativeVelocityForEquilibrium = 0.3f;

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

				humanState.TransitionToState(
					StateType::HumanNpcStateType::BehaviorType::Constrained_KnockedOut,
					0.0f,
					0.0f);


				mEventDispatcher.OnHumanNpcBehaviorChanged("Constrained_KnockedOut");

				break;
			}

			// Check alignment

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

			// Apply torque

			vec2f const desiredSecondaryParticlePosition = mParticles.GetPosition(primaryParticleState.ParticleIndex) - LabParameters::GravityDir * LabParameters::HumanNpcLength;
			// TODOHERE			

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