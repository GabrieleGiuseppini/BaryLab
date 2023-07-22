/***************************************************************************************
* Original Author:      Gabriele Giuseppini
* Created:              2018-05-13
* Copyright:            Gabriele Giuseppini  (https://github.com/GabrieleGiuseppini)
***************************************************************************************/
#include "Triangles.h"

void Triangles::Add(
    ElementIndex particleAIndex,
    ElementIndex particleBIndex,
    ElementIndex particleCIndex,
    ElementIndex subEdgeAIndex,
    ElementIndex subEdgeBIndex,
    ElementIndex subEdgeCIndex)
{
    mEndpointsBuffer.emplace_back(particleAIndex, particleBIndex, particleCIndex);
    mSubEdgesBuffer.emplace_back(subEdgeAIndex, subEdgeBIndex, subEdgeCIndex);
}

vec3f Triangles::ToBarycentricCoordinates(
    ElementIndex triangleElementIndex,
    ElementIndex vertexIndex,
    Vertices const & vertices)
{
    // TODOHERE
    (void)triangleElementIndex;
    (void)vertexIndex;
    (void)vertices;
    return vec3f::zero();
}
