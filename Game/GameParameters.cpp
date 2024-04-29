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
	, NpcSpringReductionFraction(0.39f) // At least for humans, larger than this causes heads to pull up feet
	, NpcSpringDampingCoefficient(0.5f)
	, HumanNpcEquilibriumTorqueStiffnessCoefficient(0.004f)
	, HumanNpcEquilibriumTorqueDampingCoefficient(0.0016f)
	, HumanNpcWalkingSpeedAdjustment(1.0f)
	, HumanNpcBodyLengthAdjustment(1.0f)
{
}