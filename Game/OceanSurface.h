/***************************************************************************************
 * Original Author:     Gabriele Giuseppini
 * Created:             2020-05-15
 * Copyright:           Gabriele Giuseppini  (https://github.com/GabrieleGiuseppini)
 ***************************************************************************************/
#pragma once

#include <GameCore/Vectors.h>

// Placeholder for real OceanSurface

namespace Physics {

class OceanSurface
{
public:

	OceanSurface(float depth)
		: mDepth(depth)
	{}

	static float constexpr MinDepth = -100.0f;
	static float constexpr MaxDepth = 100.0f;

	float GetDepth() const
	{
		return mDepth;
	}

	void SetDepth(float value)
	{
		mDepth = value;
	}

	inline float GetDepth(vec2f const & position) const noexcept
	{
		return mDepth - position.y;
	}

private:

	float mDepth;
};

}