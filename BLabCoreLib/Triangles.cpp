/***************************************************************************************
* Original Author:      Gabriele Giuseppini
* Created:              2018-05-13
* Copyright:            Gabriele Giuseppini  (https://github.com/GabrieleGiuseppini)
***************************************************************************************/
#include "Triangles.h"

#include "BLabMath.h"

#include <cassert>

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

bool Triangles::ContainsPoint(
    vec2f const & position,
    ElementIndex triangleElementIndex,
    Vertices const & vertices) const
{
    vec2f const aPosition = vertices.GetPosition(GetVertexAIndex(triangleElementIndex));
    vec2f const bPosition = vertices.GetPosition(GetVertexBIndex(triangleElementIndex));
    vec2f const cPosition = vertices.GetPosition(GetVertexCIndex(triangleElementIndex));

    return (position - aPosition).cross(bPosition - aPosition) >= 0.0f
        && (position - bPosition).cross(cPosition - bPosition) >= 0.0f
        && (position - cPosition).cross(aPosition - cPosition) >= 0.0f;
}

ElementIndex Triangles::FindContaining(
    vec2f const & position,
    Vertices const & vertices) const
{
    for (auto const t : *this)
    {
        if (ContainsPoint(position, t, vertices))
        {
            return t;
        }
    }

    return NoneElementIndex;
}

vec3f Triangles::ToBarycentricCoordinates(
    vec2f const & position,
    ElementIndex triangleElementIndex,
    Vertices const & vertices) const
{
    vec2f abBaryCoords = InternalToBarycentricCoordinates(
        position,
        triangleElementIndex,
        vertices);

    return vec3f(
        abBaryCoords.x,
        abBaryCoords.y,
        1.0f - abBaryCoords.x - abBaryCoords.y);
}

vec3f Triangles::ToBarycentricCoordinatesFromWithinTriangle(
    vec2f const & position,
    ElementIndex triangleElementIndex,
    Vertices const & vertices) const
{
    assert(ContainsPoint(position, triangleElementIndex, vertices));

    vec2f abBaryCoords = InternalToBarycentricCoordinates(
        position,
        triangleElementIndex,
        vertices).clamp(0.0f, 1.0f, 0.0f, 1.0f);

    return vec3f(
        abBaryCoords.x,
        abBaryCoords.y,
        1.0f - abBaryCoords.x - abBaryCoords.y);
}

vec2f Triangles::FromBarycentricCoordinates(
    vec3f const & barycentricCoordinates,
    ElementIndex triangleElementIndex,
    Vertices const & vertices) const
{
    vec2f const & positionA = vertices.GetPosition(mEndpointsBuffer[triangleElementIndex].VertexIndices[0]);
    vec2f const & positionB = vertices.GetPosition(mEndpointsBuffer[triangleElementIndex].VertexIndices[1]);
    vec2f const & positionC = vertices.GetPosition(mEndpointsBuffer[triangleElementIndex].VertexIndices[2]);

    return
        positionA * barycentricCoordinates.x
        + positionB * barycentricCoordinates.y
        + positionC * barycentricCoordinates.z;
}

//////////////////////////////////////////////////////////

vec2f Triangles::InternalToBarycentricCoordinates(
    vec2f const & position,
    ElementIndex triangleElementIndex,
    Vertices const & vertices) const
{
    vec2f const & positionA = vertices.GetPosition(mEndpointsBuffer[triangleElementIndex].VertexIndices[0]);
    vec2f const & positionB = vertices.GetPosition(mEndpointsBuffer[triangleElementIndex].VertexIndices[1]);
    vec2f const & positionC = vertices.GetPosition(mEndpointsBuffer[triangleElementIndex].VertexIndices[2]);

    float const denominator =
        (positionB.y - positionC.y) * (positionA.x - positionC.x)
        + (positionC.x - positionB.x) * (positionA.y - positionC.y);

    if (IsAlmostZero(denominator))
    {
        // Co-linear, put arbitrarily in center
        float constexpr l = 0.3333333f;
        return vec2f(l, l);
    }
    else
    {
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

        return vec2f(
            l1,
            l2);
    }
}
