/***************************************************************************************
 * Original Author:     Gabriele Giuseppini
 * Created:             2020-05-15
 * Copyright:           Gabriele Giuseppini  (https://github.com/GabrieleGiuseppini)
 ***************************************************************************************/
#include "OceanFloor.h"

namespace Physics {

std::tuple<bool, float, int> OceanFloor::GetHeightIfUnderneathAt(float /*x*/, float /*y*/) const
{
	return { false, 0.0f, 0 };
}

vec2f OceanFloor::GetNormalAt(int) const
{
	return vec2f::zero();
}

}