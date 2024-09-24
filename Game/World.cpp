/***************************************************************************************
 * Original Author:     Gabriele Giuseppini
 * Created:             2020-05-15
 * Copyright:           Gabriele Giuseppini  (https://github.com/GabrieleGiuseppini)
 ***************************************************************************************/
#include "World.h"

namespace Physics {

// Hiding these from compiler

vec2f World::GetCurrentWindSpeed() const
{
	return vec2f::zero();
}

std::optional<Wind::RadialWindField> World::GetCurrentRadialWindField() const
{
	return std::nullopt;
}

}