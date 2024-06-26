/***************************************************************************************
* Original Author:		Gabriele Giuseppini
* Created:				2023-07-19
* Copyright:			Gabriele Giuseppini  (https://github.com/GabrieleGiuseppini)
***************************************************************************************/
#include "GameParameters.h"

GameParameters::GameParameters()
	: WaterFrictionDragAdjustment(1.0f)
	, BuoyancyAdjustment(1.0f)
	// NPC
	, NpcDamping(0.0078f)
	, NpcMaterialElasticityAdjustment(1.0f)
	, NpcMaterialStaticFrictionAdjustment(1.0f)
	, NpcMaterialKineticFrictionAdjustment(1.0f)
	, NpcSpringReductionFractionAdjustment(1.0f)
	, NpcSpringDampingCoefficientAdjustment(1.0f)
	, HumanNpcEquilibriumTorqueStiffnessCoefficient(0.0035f)
	, HumanNpcEquilibriumTorqueDampingCoefficient(0.0012f)
	, HumanNpcWalkingSpeedAdjustment(1.0f)
	, HumanNpcBodyLengthAdjustment(1.0f)
{
}