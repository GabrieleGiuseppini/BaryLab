/***************************************************************************************
* Original Author:      Gabriele Giuseppini
* Created:              2018-06-20
* Copyright:            Gabriele Giuseppini  (https://github.com/GabrieleGiuseppini)
***************************************************************************************/
#pragma once

#include "BLabMath.h"
#include "BLabTypes.h"
#include "Vectors.h"

#include <cassert>
#include <cmath>

namespace Geometry {

inline bool IsPointInTriangle(
    vec2f const & pPosition,
    vec2f const & aPosition,
    vec2f const & bPosition,
    vec2f const & cPosition)
{
    return (pPosition - aPosition).cross(bPosition - aPosition) >= 0.0f
        && (pPosition - bPosition).cross(cPosition - bPosition) >= 0.0f
        && (pPosition - cPosition).cross(aPosition - cPosition) >= 0.0f;
}

}