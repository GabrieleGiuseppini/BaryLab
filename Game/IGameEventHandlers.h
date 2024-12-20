/***************************************************************************************
* Original Author:		Gabriele Giuseppini
* Created:				2020-05-15
* Copyright:			Gabriele Giuseppini  (https://github.com/GabrieleGiuseppini)
***************************************************************************************/
#pragma once

#include <GameCore/BarycentricCoords.h>
#include <GameCore/GameTypes.h>
#include <GameCore/Vectors.h>

#include <chrono>
#include <optional>
#include <string>

/*
 * These interfaces define the methods that event handlers must implement.
 *
 * The methods are default-implemented to facilitate implementation of handlers that
 * only care about a subset of the events.
 */

struct IBLabEventHandler
{
    virtual void OnBLabReset()
    {
        // Default-implemented
    }

    virtual void OnSubjectParticleBarycentricCoordinatesWrtOriginTriangleChanged(std::optional<bcoords3f> const & /*coordinates*/)
    {
        // Default implemented
    }

    virtual void OnSubjectParticleConstrainedRegimeUpdated(std::optional<AbsoluteTriangleBCoords> const & /*constrainedRegimeParticleProbe*/)
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

    virtual void OnUpdateTimeMeasured(
        float /*updateDurationMilliseconds*/,
        float /*renderUploadDurationMilliseconds*/)
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

struct IPlaceholderEventHandler
{
    virtual void OnBombExplosion(
        GadgetType /*gadgetType*/,
        bool /*isUnderwater*/,
        unsigned int /*size*/)
    {
        // Default-implemented
    }

    virtual void OnWaterReactionExplosion(
        bool /*isUnderwater*/,
        unsigned int /*size*/)
    {
        // Default-implemented
    }

    virtual void OnPointCombustionBegin()
    {
        // Default-implemented
    }

    virtual void OnPointCombustionEnd()
    {
        // Default-implemented
    }

    virtual void OnCombustionExplosion(
        bool /*isUnderwater*/,
        unsigned int /*size*/)
    {
        // Default-implemented
    }
};

struct INpcGameEventHandler
{
    virtual void OnNpcSelectionChanged(
        std::optional<NpcId> /*selectedNpc*/)
    {
        // Default-implemented
    }

    virtual void OnNpcCountsUpdated(
        size_t /*totalNpcCount*/)
    {
        // Default-implemented
    }

    virtual void OnHumanNpcCountsUpdated(
        size_t /*insideShipCount*/,
        size_t /*outsideShipCount*/)
    {
        // Default-implemented
    }
};
