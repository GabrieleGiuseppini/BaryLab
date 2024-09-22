/***************************************************************************************
* Original Author:		Gabriele Giuseppini
* Created:				2019-01-16
* Copyright:			Gabriele Giuseppini  (https://github.com/GabrieleGiuseppini)
***************************************************************************************/
#pragma once

#include <GameCore/Vectors.h>

#include <optional>

namespace Physics
{

/*
 * Wind consists of two components:
 * - A linear (horizontal) wind, whose intensity is modulated by various actors (storm, gusting state machine, etc.)
 * - A radial wind, when triggered interactively
 */
class Wind
{
public:

    struct RadialWindField
    {
        vec2f SourcePos;
        float PreFrontRadius;
        float PreFrontWindForceMagnitude;
        float MainFrontRadius;
        float MainFrontWindForceMagnitude;

        RadialWindField(
            vec2f sourcePos,
            float preFrontRadius,
            float preFrontWindForceMagnitude,
            float mainFrontRadius,
            float mainFrontWindForceMagnitude)
            : SourcePos(sourcePos)
            , PreFrontRadius(preFrontRadius)
            , PreFrontWindForceMagnitude(preFrontWindForceMagnitude)
            , MainFrontRadius(mainFrontRadius)
            , MainFrontWindForceMagnitude(mainFrontWindForceMagnitude)
        {}
    };

public:

    Wind()
        : mCurrentRadialWindField()
    {}

    /*
     * Returns the current radial wind field, if any.
     */
    std::optional<RadialWindField> const & GetCurrentRadialWindField() const
    {
        return mCurrentRadialWindField;
    }

private:

    std::optional<RadialWindField> mCurrentRadialWindField;
};

}
