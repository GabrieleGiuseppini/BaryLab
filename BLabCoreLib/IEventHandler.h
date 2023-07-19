/***************************************************************************************
* Original Author:		Gabriele Giuseppini
* Created:				2020-05-15
* Copyright:			Gabriele Giuseppini  (https://github.com/GabrieleGiuseppini)
***************************************************************************************/
#pragma once

#include "BLabTypes.h"

#include <chrono>
#include <optional>
#include <string>

/*
 * These interfaces define the methods that event handlers must implement.
 *
 * The methods are default-implemented to facilitate implementation of handlers that
 * only care about a subset of the events.
 */

struct IEventHandler
{
    virtual void OnReset(size_t /*numSprings*/)
    {
        // Default-implemented
    }

    virtual void OnCustomProbe(
        std::string const & /*name*/,
        float /*value*/)
    {
        // Default-implemented
    }
};