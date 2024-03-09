/***************************************************************************************
* Original Author:		Gabriele Giuseppini
* Created:				2020-05-15
* Copyright:			Gabriele Giuseppini  (https://github.com/GabrieleGiuseppini)
***************************************************************************************/
#include "GameTypes.h"

#include "GameException.h"
#include "Utils.h"

NpcSurfaceType StrToNpcSurfaceType(std::string const & str)
{
    if (Utils::CaseInsensitiveEquals(str, "Floor"))
        return NpcSurfaceType::Floor;
    else if (Utils::CaseInsensitiveEquals(str, "Open"))
        return NpcSurfaceType::Open;
    else
        throw GameException("Unrecognized NpcSurfaceType \"" + str + "\"");
}