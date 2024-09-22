/***************************************************************************************
* Original Author:		Gabriele Giuseppini
* Created:				2019-10-15
* Copyright:			Gabriele Giuseppini  (https://github.com/GabrieleGiuseppini)
***************************************************************************************/
#pragma once

#include "GameEventDispatcher.h"
#include "GameParameters.h"
#include "RenderContext.h"

#include <GameCore/GameWallClock.h>
#include <GameCore/Vectors.h>

#include <memory>
#include <list>
#include <vector>

namespace Physics
{

class Storm
{
public:

    struct Parameters
    {
        float WindSpeed; // Km/h, absolute (on top of current direction)
        unsigned int NumberOfClouds;
        float CloudsSize; // [0.0f = initial size, 1.0 = full size]
        float CloudDarkening; // [0.0f = full darkness, 1.0 = no darkening]
        float AmbientDarkening; // [0.0f = full darkness, 1.0 = no darkening]
        float RainDensity; // [0.0f = no rain, 1.0f = full rain]
        float RainQuantity; // m/h
        float AirTemperatureDelta; // K

        Parameters()
        {
            Reset();
        }

        void Reset()
        {
            WindSpeed = 0.0f;
            NumberOfClouds = 0;
            CloudsSize = 0.0f;
            CloudDarkening = 1.0f;
            AmbientDarkening = 1.0f;
            RainDensity = 0.0f;
            RainQuantity = 0.0f;
            AirTemperatureDelta = 0.0f;
        }
    };

    Parameters const & GetParameters() const
    {
        return mParameters;
    }

public:

	Storm()
        : mParameters()
    {}

private:

	Parameters mParameters;
};

}
