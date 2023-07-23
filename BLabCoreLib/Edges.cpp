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
    TrianglesVector const & triangles)
{
    mEndpointsBuffer.emplace_back(vertexAIndex, vertexBIndex);
    mEndpointOctantsBuffer.emplace_back(pointAOctant, pointBOctant);
    mTrianglesBuffer.emplace_back(triangles);

    mRenderColorBuffer.emplace_back(vec4f::zero());
}