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

    void OnSubjectParticleBarycentricCoordinatesChanged(vec3f const & coordinates) override
    {
        for (auto sink : mSinks)
        {
            sink->OnSubjectParticleBarycentricCoordinatesChanged(coordinates);
        }
    }

    void OnTrajectoryToggled(bool isTrajectorySet) override
    {
        for (auto sink : mSinks)
        {
            sink->OnTrajectoryToggled(isTrajectorySet);
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
