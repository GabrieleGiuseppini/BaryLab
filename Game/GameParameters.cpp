/***************************************************************************************
* Original Author:		Gabriele Giuseppini
* Created:				2023-07-19
* Copyright:			Gabriele Giuseppini  (https://github.com/GabrieleGiuseppini)
***************************************************************************************/
#include "GameParameters.h"

GameParameters::GameParameters()
	: GlobalDamping(0.0078f)
	, ElasticityAdjustment(1.0f)
	, StaticFrictionAdjustment(1.0f)
	, KineticFrictionAdjustment(1.0f)
	, WaterFrictionDragCoefficientAdjustment(1.0f)
	, BuoyancyAdjustment(1.0f)
	// Misc
	, ToolSearchRadius(0.3f)
	// NPC
	, NpcSpringReductionFractionAdjustment(1.0f)
	, NpcSpringDampingCoefficientAdjustment(1.0f)
	, HumanNpcEquilibriumTorqueStiffnessCoefficient(0.0035f)
	, HumanNpcEquilibriumTorqueDampingCoefficient(0.0012f)
	, HumanNpcWalkingSpeedAdjustment(1.0f)
	, HumanNpcBodyLengthAdjustment(1.0f)
{
}