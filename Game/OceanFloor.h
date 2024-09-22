/***************************************************************************************
 * Original Author:     Gabriele Giuseppini
 * Created:             2020-05-15
 * Copyright:           Gabriele Giuseppini  (https://github.com/GabrieleGiuseppini)
 ***************************************************************************************/
#pragma once

#include <GameCore/Vectors.h>

// Placeholder for real OceanSurface

namespace Physics {

class OceanFloor
{
public:

	OceanFloor()
	{}

	std::tuple<bool, float, int> GetHeightIfUnderneathAt(float x, float y) const;

	vec2f GetNormalAt(int) const;
};

}