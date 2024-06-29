/***************************************************************************************
* Original Author:		Gabriele Giuseppini
* Created:				2023-07-19
* Copyright:			Gabriele Giuseppini  (https://github.com/GabrieleGiuseppini)
***************************************************************************************/
#include "GameParameters.h"

GameParameters::GameParameters()
	: GravityAdjustment(1.0f)
	, MassAdjustment(1.0f)
	, GlobalDampingAdjustment(1.0f)
	, WaterFrictionDragAdjustment(1.0f)
	, BuoyancyAdjustment(1.0f)
	// NPC
	, NpcMaterialElasticityAdjustment(1.0f)
	, NpcMaterialStaticFrictionAdjustment(1.0f)
	, NpcMaterialKineticFrictionAdjustment(1.0f)
	, NpcSpringReductionFractionAdjustment(1.0f)
	, NpcSpringDampingCoefficientAdjustment(1.0f)
	, NpcSizeAdjustment(1.0f)
	, HumanNpcEquilibriumTorqueStiffnessCoefficient(0.0035f)
	, HumanNpcEquilibriumTorqueDampingCoefficient(0.0012f)
	, HumanNpcWalkingSpeedAdjustment(1.0f)
{
}