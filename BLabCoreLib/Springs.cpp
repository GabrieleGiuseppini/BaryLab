/***************************************************************************************
 * Original Author:      Gabriele Giuseppini
 * Created:              2020-05-16
 * Copyright:            Gabriele Giuseppini  (https://github.com/GabrieleGiuseppini)
 ***************************************************************************************/
#include "Springs.h"

void Springs::Add(
    ElementIndex pointAIndex,
    ElementIndex pointBIndex,
    Octant pointAOctant,
    Octant pointBOctant,
    TrianglesVector const & triangles)
{
    mEndpointsBuffer.emplace_back(pointAIndex, pointBIndex);
    mEndpointOctantsBuffer.emplace_back(pointAOctant, pointBOctant);
    mTrianglesBuffer.emplace_back(triangles);

    mRenderColorBuffer.emplace_back(vec4f::zero());
    mRenderHighlightBuffer.emplace_back(0.0f);
}