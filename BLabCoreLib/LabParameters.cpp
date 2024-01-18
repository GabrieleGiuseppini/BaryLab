/***************************************************************************************
* Original Author:		Gabriele Giuseppini
* Created:				2023-07-19
* Copyright:			Gabriele Giuseppini  (https://github.com/GabrieleGiuseppini)
***************************************************************************************/
#include "LabParameters.h"

LabParameters::LabParameters()
	: MassAdjustment(1.0f)
	, GravityAdjustment(1.0f)
	, GlobalDamping(0.0078f)
	, SeaLevel(-2.0f)
	, SpringReductionFraction(0.39f) // At least for humans, larger than this causes heads to pull up feet
	, SpringDampingCoefficient(0.5f)
	, Elasticity(0.6f)
	, StaticFrictionAdjustment(1.0f)
	, KineticFrictionAdjustment(1.0f)
	, WaterFrictionDragCoefficientAdjustment(1.0f)
	, BuoyancyAdjustment(1.0f)
	// NPC
	, HumanNpcEquilibriumTorqueStiffnessCoefficient(0.0032f)
	, HumanNpcEquilibriumTorqueDampingCoefficient(0.055f)
	, HumanNpcWalkingAcceleration(0.027f)
	, HumanNpcWalkingSpeed(1.0f)
{
}