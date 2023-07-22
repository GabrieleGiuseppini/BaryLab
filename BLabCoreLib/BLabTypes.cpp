/***************************************************************************************
* Original Author:		Gabriele Giuseppini
* Created:				2020-05-15
* Copyright:			Gabriele Giuseppini  (https://github.com/GabrieleGiuseppini)
***************************************************************************************/
#include "BLabTypes.h"

#include "BLabException.h"

#include "Utils.h"

SurfaceType StrToSurfaceType(std::string const & str)
{
    if (Utils::CaseInsensitiveEquals(str, "Floor"))
        return SurfaceType::Floor;
    else if (Utils::CaseInsensitiveEquals(str, "Open"))
        return SurfaceType::Open;
    else
        throw BLabException("Unrecognized SurfaceType \"" + str + "\"");
}