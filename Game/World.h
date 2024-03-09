/***************************************************************************************
 * Original Author:     Gabriele Giuseppini
 * Created:             2020-05-15
 * Copyright:           Gabriele Giuseppini  (https://github.com/GabrieleGiuseppini)
 ***************************************************************************************/
#pragma once

#include "OceanSurface.h"

 // Placeholder for real World

namespace Physics {

class World
{
public:

	World()
		: mOceanSurface()
	{}

	OceanSurface const & GetOceanSurface() const
	{
		return mOceanSurface;
	}

	OceanSurface & GetOceanSurface()
	{
		return mOceanSurface;
	}

private:

	OceanSurface mOceanSurface;
};

}