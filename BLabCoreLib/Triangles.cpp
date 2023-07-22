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
    vec2f const & position,
    ElementIndex triangleElementIndex,
    Vertices const & vertices)
{
    vec2f const & positionA = vertices.GetPosition(mEndpointsBuffer[triangleElementIndex].VertexIndices[0]);
    vec2f const & positionB = vertices.GetPosition(mEndpointsBuffer[triangleElementIndex].VertexIndices[1]);
    vec2f const & positionC = vertices.GetPosition(mEndpointsBuffer[triangleElementIndex].VertexIndices[2]);

    float const denominator =
        (positionB.y - positionC.y) * (positionA.x - positionC.x)
        + (positionC.x - positionB.x) * (positionA.y - positionC.y);

    // TODO: colinearity

    float const l1 =
        (
            (positionB.y - positionC.y) * (position.x - positionC.x)
            + (positionC.x - positionB.x) * (position.y - positionC.y)
        ) / denominator;
    
    float const l2 =
        (
            (positionC.y - positionA.y) * (position.x - positionC.x)
            + (positionA.x - positionC.x) * (position.y - positionC.y)
        ) / denominator;

    float const l3 = 1.0f - l1 - l2;

    return vec3f(
        l1,
        l2,
        l3);
}
