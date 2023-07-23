/***************************************************************************************
 * Original Author:      Gabriele Giuseppini
 * Created:              2020-05-16
 * Copyright:            Gabriele Giuseppini  (https://github.com/GabrieleGiuseppini)
 ***************************************************************************************/
#include "Edges.h"

void Edges::Add(
    ElementIndex vertexAIndex,
    ElementIndex vertexBIndex,
    Octant pointAOctant,
    Octant pointBOctant,
    SurfaceType surface,
    TrianglesVector const & triangles)
{
    mEndpointsBuffer.emplace_back(vertexAIndex, vertexBIndex);
    mEndpointOctantsBuffer.emplace_back(pointAOctant, pointBOctant);
    mTrianglesBuffer.emplace_back(triangles);

    mSurfaceTypeBuffer.emplace_back(surface);

    mRenderColorBuffer.emplace_back(
        surface == SurfaceType::Floor 
        ? rgbaColor(0x0a, 0x0a, 0x0a, 0xff) 
        : rgbaColor(0xba, 0xba, 0xba, 0xff));
}