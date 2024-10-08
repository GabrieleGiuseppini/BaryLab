/***************************************************************************************
* Original Author:		Gabriele Giuseppini
* Created:				2020-05-13
* Copyright:			Gabriele Giuseppini  (https://github.com/GabrieleGiuseppini)
***************************************************************************************/
#pragma once

#include "IGameEventHandlers.h"

#include <vector>

class GameEventDispatcher final : 
    public IBLabEventHandler, 
    public IPlaceholderEventHandler,
    public INpcGameEventHandler
{
public:

    GameEventDispatcher()
        : mBLabSinks()
    {
    }

public:

    void OnBLabReset() override
    {
        for (auto sink : mBLabSinks)
        {
            sink->OnBLabReset();
        }
    }

    void OnSubjectParticleBarycentricCoordinatesWrtOriginTriangleChanged(std::optional<bcoords3f> const & coordinates) override
    {
        for (auto sink : mBLabSinks)
        {
            sink->OnSubjectParticleBarycentricCoordinatesWrtOriginTriangleChanged(coordinates);
        }
    }

    void OnSubjectParticleConstrainedRegimeUpdated(std::optional<AbsoluteTriangleBCoords> const & constrainedRegimeParticleProbe) override
    {
        for (auto sink : mBLabSinks)
        {
            sink->OnSubjectParticleConstrainedRegimeUpdated(constrainedRegimeParticleProbe);
        }
    }

    void OnSubjectParticlePhysicsUpdated(std::optional<PhysicsParticleProbe> const & physicsParticleProbe) override
    {
        for (auto sink : mBLabSinks)
        {
            sink->OnSubjectParticlePhysicsUpdated(physicsParticleProbe);
        }
    }

    void OnTrajectoryToggled(bool isTrajectorySet) override
    {
        for (auto sink : mBLabSinks)
        {
            sink->OnTrajectoryToggled(isTrajectorySet);
        }
    }

    void OnHumanNpcBehaviorChanged(
        std::optional<std::string> behavior) override
    {
        for (auto sink : mBLabSinks)
        {
            sink->OnHumanNpcBehaviorChanged(behavior);
        }
    }

    void OnHumanNpcStateQuantityChanged(
        std::optional<std::tuple<std::string, std::string>> nameAndValue) override
    {
        for (auto sink : mBLabSinks)
        {
            sink->OnHumanNpcStateQuantityChanged(nameAndValue);
        }
    }

    void OnUpdateTimeMeasured(
        float updateDurationMilliseconds,
        float renderUploadDurationMilliseconds) override
    {
        for (auto sink : mBLabSinks)
        {
            sink->OnUpdateTimeMeasured(
                updateDurationMilliseconds,
                renderUploadDurationMilliseconds);
        }
    }

    void OnCustomProbe(
        std::string const & name,
        float value) override
    {
        for (auto sink : mBLabSinks)
        {
            sink->OnCustomProbe(name, value);
        }
    }

    //
    // NPC
    //

    void OnNpcSelectionChanged(
        std::optional<NpcId> selectedNpc) override
    {
        for (auto sink : mNpcSinks)
        {
            sink->OnNpcSelectionChanged(selectedNpc);
        }
    }

    void OnNpcCountsUpdated(
        size_t totalNpcCount) override
    {
        for (auto sink : mNpcSinks)
        {
            sink->OnNpcCountsUpdated(totalNpcCount);
        }
    }

    void OnHumanNpcCountsUpdated(
        size_t insideShipCount,
        size_t outsideShipCount) override
    {
        for (auto sink : mNpcSinks)
        {
            sink->OnHumanNpcCountsUpdated(insideShipCount, outsideShipCount);
        }
    }

public:

    void RegisterBLabEventHandler(IBLabEventHandler * sink)
    {
        mBLabSinks.push_back(sink);
    }

    void RegisterNpcGameEventHandler(INpcGameEventHandler * sink)
    {
        mNpcSinks.push_back(sink);
    }

private:

    // The registered sinks
    std::vector<IBLabEventHandler *> mBLabSinks;
    std::vector<INpcGameEventHandler *> mNpcSinks;
};
