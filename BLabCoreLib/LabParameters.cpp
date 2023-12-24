/***************************************************************************************
* Original Author:		Gabriele Giuseppini
* Created:				2023-07-19
* Copyright:			Gabriele Giuseppini  (https://github.com/GabrieleGiuseppini)
***************************************************************************************/
#include "LabParameters.h"

LabParameters::LabParameters()
	: MassAdjustment(1.0f)
	, GravityAdjustment(1.0f)
	, SpringReductionFraction(0.5f)
	, SpringDampingCoefficient(0.5f)
	, Elasticity(0.6f)
	, StaticFriction(0.4f)
	, KineticFriction(0.22f)
	// NPC
	, HumanNpcEquilibriumTorqueStiffnessCoefficient(0.0025f)
	, HumanNpcEquilibriumTorqueDampingCoefficient(0.055f)
	, HumanNpcWalkingSpeed(1.42f)
{
}