/***************************************************************************************
* Original Author:		Gabriele Giuseppini
* Created:				2023-07-19
* Copyright:			Gabriele Giuseppini  (https://github.com/GabrieleGiuseppini)
***************************************************************************************/
#include "GameParameters.h"

GameParameters::GameParameters()
	: GravityAdjustment(1.0f)
	, MassAdjustment(1.0f)
	, AirDensityAdjustment(1.0f)
	, WaterDensityAdjustment(1.0f)
	, GlobalDampingAdjustment(1.0f)
	, ElasticityAdjustment(1.0f)
	, StaticFrictionAdjustment(1.0f)
	, KineticFrictionAdjustment(1.0f)
	, WaterFrictionDragAdjustment(1.0f)
	, BuoyancyAdjustment(1.0f)
	, OceanFloorElasticityCoefficient(0.5f)
	, MoveToolInertia(3.0f)
	, IsUltraViolentMode(false)
	// NPC
	, NpcSpringReductionFractionAdjustment(1.0f)
	, NpcSpringDampingCoefficientAdjustment(1.0f)
	, NpcFrictionAdjustment(1.0f)
	, NpcSizeMultiplier(1.0f)
	, HumanNpcEquilibriumTorqueStiffnessCoefficient(0.0035f)
	, HumanNpcEquilibriumTorqueDampingCoefficient(0.0012f)
	, HumanNpcWalkingSpeedAdjustment(1.0f)
{
}