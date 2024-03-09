/***************************************************************************************
 * Original Author:      Gabriele Giuseppini
 * Created:              2020-05-16
 * Copyright:            Gabriele Giuseppini  (https://github.com/GabrieleGiuseppini)
 ***************************************************************************************/
#include "Physics.h"

namespace Physics {

void Springs::Add(
    ElementIndex pointAIndex,
    ElementIndex pointBIndex,
    Octant pointAOctant,
    Octant pointBOctant,
    SurfaceType surface,
    TrianglesVector const & triangles)
{
    mEndpointsBuffer.emplace_back(pointAIndex, pointBIndex);
    mEndpointOctantsBuffer.emplace_back(pointAOctant, pointBOctant);
    mTrianglesBuffer.emplace_back(triangles);

    mSurfaceTypeBuffer.emplace_back(surface);

    mRenderColorBuffer.emplace_back(
        surface == SurfaceType::Floor
        ? rgbaColor(0x0a, 0x0a, 0x0a, 0xff)
        : rgbaColor(0xba, 0xba, 0xba, 0xff));
}

}