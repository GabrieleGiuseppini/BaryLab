/***************************************************************************************
* Original Author:		Gabriele Giuseppini
* Created:				2020-05-13
* Copyright:			Gabriele Giuseppini  (https://github.com/GabrieleGiuseppini)
***************************************************************************************/
#pragma once

#include "IEventHandler.h"

#include <vector>

class EventDispatcher final : public IEventHandler
{
public:

    EventDispatcher()
        : mSinks()
    {
    }

public:

    void OnReset() override
    {
        for (auto sink : mSinks)
        {
            sink->OnReset();
        }
    }

    void OnSubjectParticleBarycentricCoordinatesWrtOriginTriangleChanged(std::optional<bcoords3f> const & coordinates) override
    {
        for (auto sink : mSinks)
        {
            sink->OnSubjectParticleBarycentricCoordinatesWrtOriginTriangleChanged(coordinates);
        }
    }

    void OnSubjectParticleConstrainedRegimeUpdated(std::optional<ConstrainedRegimeParticleProbe> const & constrainedRegimeParticleProbe) override
    {
        for (auto sink : mSinks)
        {
            sink->OnSubjectParticleConstrainedRegimeUpdated(constrainedRegimeParticleProbe);
        }
    }

    void OnSubjectParticlePhysicsUpdated(std::optional<PhysicsParticleProbe> const & physicsParticleProbe) override
    {
        for (auto sink : mSinks)
        {
            sink->OnSubjectParticlePhysicsUpdated(physicsParticleProbe);
        }
    }

    void OnTrajectoryToggled(bool isTrajectorySet) override
    {
        for (auto sink : mSinks)
        {
            sink->OnTrajectoryToggled(isTrajectorySet);
        }
    }

    void OnHumanNpcBehaviorChanged(
        std::optional<std::string> behavior) override
    {
        for (auto sink : mSinks)
        {
            sink->OnHumanNpcBehaviorChanged(behavior);
        }
    }

    void OnHumanNpcStateQuantityChanged(
        std::optional<std::tuple<std::string, std::string>> nameAndValue) override
    {
        for (auto sink : mSinks)
        {
            sink->OnHumanNpcStateQuantityChanged(nameAndValue);
        }
    }

    void OnCustomProbe(
        std::string const & name,
        float value) override
    {
        for (auto sink : mSinks)
        {
            sink->OnCustomProbe(name, value);
        }
    }

public:

    void RegisterEventHandler(IEventHandler * sink)
    {
        mSinks.push_back(sink);
    }

private:

    // The registered sinks
    std::vector<IEventHandler *> mSinks;
};
