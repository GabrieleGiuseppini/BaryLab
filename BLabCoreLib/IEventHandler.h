/***************************************************************************************
* Original Author:		Gabriele Giuseppini
* Created:				2020-05-15
* Copyright:			Gabriele Giuseppini  (https://github.com/GabrieleGiuseppini)
***************************************************************************************/
#pragma once

#include "BLabTypes.h"
#include "Vectors.h"

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
    virtual void OnReset()
    {
        // Default-implemented
    }

    virtual void OnSubjectParticleBarycentricCoordinatesWrtOriginTriangleChanged(std::optional<vec3f> const & /*coordinates*/)
    {
        // Default implemented
    }

    virtual void OnSubjectParticleConstrainedRegimeUpdated(std::optional<ConstrainedRegimeParticleProbe> const & /*constrainedRegimeParticleProbe*/)
    {
        // Default implemented
    }

    virtual void OnSubjectParticlePhysicsUpdated(std::optional<PhysicsParticleProbe> const & /*physicsParticleProbe*/)
    {
        // Default implemented
    }

    virtual void OnTrajectoryToggled(bool /*isTrajectorySet*/)
    {
        // Default-implemented
    }

    virtual void OnHumanNpcBehaviorChanged(
        std::optional<std::string> /*behavior*/)
    {
        // Default-implemented
    }

    virtual void OnHumanNpcStateQuantityChanged(
        std::optional<std::tuple<std::string, std::string>> /*nameAndValue*/)
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
