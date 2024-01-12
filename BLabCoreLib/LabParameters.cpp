/***************************************************************************************
* Original Author:		Gabriele Giuseppini
* Created:				2023-07-19
* Copyright:			Gabriele Giuseppini  (https://github.com/GabrieleGiuseppini)
***************************************************************************************/
#include "LabParameters.h"

LabParameters::LabParameters()
	: MassAdjustment(1.0f)
	, GravityAdjustment(1.0f)
	, GlobalDamping(0.00010749653315f)
	, SpringReductionFraction(0.5f)
	, SpringDampingCoefficient(0.5f)
	, Elasticity(0.6f)
	, StaticFrictionAdjustment(1.0f)
	, KineticFrictionAdjustment(1.0f)
	// NPC
	, HumanNpcEquilibriumTorqueStiffnessCoefficient(0.0029f)
	, HumanNpcEquilibriumTorqueDampingCoefficient(0.055f)
	, HumanNpcWalkingSpeed(1.0f)
{
}