/***************************************************************************************
 * Original Author:		Gabriele Giuseppini
 * Created:				2023-10-20
 * Copyright:			Gabriele Giuseppini  (https://github.com/GabrieleGiuseppini)
 ***************************************************************************************/
#include "Npcs.h"

#include "Log.h"

Npcs::StateType::HumanNpcStateType Npcs::InitializeHuman(
	StateType::NpcParticleStateType const & primaryParticleState,
	StateType::NpcParticleStateType const & secondaryParticleState,
	Mesh const & mesh) const
{
	// TODOHERE
	(void)primaryParticleState;
	(void)secondaryParticleState;
	(void)mesh;
	return StateType::HumanNpcStateType(StateType::HumanNpcStateType::BehaviorType::KnockedOut);
}

void Npcs::UpdateHuman(
	StateType::NpcParticleStateType const & primaryParticleState,
	StateType::NpcParticleStateType const & secondaryParticleState,
	Mesh const & mesh)
{
	// TODOHERE
	(void)primaryParticleState;
	(void)secondaryParticleState;
	(void)mesh;
}