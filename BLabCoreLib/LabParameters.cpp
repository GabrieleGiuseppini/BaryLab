/***************************************************************************************
* Original Author:		Gabriele Giuseppini
* Created:				2023-07-19
* Copyright:			Gabriele Giuseppini  (https://github.com/GabrieleGiuseppini)
***************************************************************************************/
#include "LabParameters.h"

LabParameters::LabParameters()
	: MassAdjustment(1.0f)
	, GravityAdjustment(1.0f)
	, SpringReductionFraction(0.1f)
	, SpringDampingCoefficient(0.03f)
	, Elasticity(0.6f)
	, StaticFriction(0.15f)
	, KineticFriction(0.031f)
{
}