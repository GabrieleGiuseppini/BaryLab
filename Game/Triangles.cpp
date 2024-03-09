/***************************************************************************************
* Original Author:      Gabriele Giuseppini
* Created:              2018-05-13
* Copyright:            Gabriele Giuseppini  (https://github.com/GabrieleGiuseppini)
***************************************************************************************/
#include "Physics.h"

#include <GameCore/GameMath.h>

#include <cassert>

namespace Physics {

void Triangles::Add(
    ElementIndex pointAIndex,
    ElementIndex pointBIndex,
    ElementIndex pointCIndex,
    ElementIndex subSpringAIndex,
    ElementIndex subSpringBIndex,
    ElementIndex subSpringCIndex,
    ElementIndex subSpringAOppositeTriangle,
    int subSpringAOppositeTriangleEdgeOrdinal,
    ElementIndex subSpringBOppositeTriangle,
    int subSpringBOppositeTriangleEdgeOrdinal,
    ElementIndex subSpringCOppositeTriangle,
    int subSpringCOppositeTriangleEdgeOrdinal,
    NpcSurfaceType subSpringANpcSurfaceType,
    NpcSurfaceType subSpringBNpcSurfaceType,
    NpcSurfaceType subSpringCNpcSurfaceType)
{
    mEndpointsBuffer.emplace_back(pointAIndex, pointBIndex, pointCIndex);
    mSubSpringsBuffer.emplace_back(subSpringAIndex, subSpringBIndex, subSpringCIndex);
    mOppositeTrianglesBuffer.emplace_back(OppositeTrianglesInfo{
        OppositeTriangleInfo(subSpringAOppositeTriangle, subSpringAOppositeTriangleEdgeOrdinal),
        OppositeTriangleInfo(subSpringBOppositeTriangle, subSpringBOppositeTriangleEdgeOrdinal),
        OppositeTriangleInfo(subSpringCOppositeTriangle, subSpringCOppositeTriangleEdgeOrdinal) });
    mSubSpringNpcSurfaceTypesBuffer.emplace_back(SubSpringNpcSurfaceTypes({ subSpringANpcSurfaceType, subSpringBNpcSurfaceType, subSpringCNpcSurfaceType }));
}

bool Triangles::ContainsPoint(
    vec2f const & position,
    ElementIndex triangleElementIndex,
    Points const & points) const
{
    vec2f const aPosition = points.GetPosition(GetPointAIndex(triangleElementIndex));
    vec2f const bPosition = points.GetPosition(GetPointBIndex(triangleElementIndex));
    vec2f const cPosition = points.GetPosition(GetPointCIndex(triangleElementIndex));

    return (position - aPosition).cross(bPosition - aPosition) >= 0.0f
        && (position - bPosition).cross(cPosition - bPosition) >= 0.0f
        && (position - cPosition).cross(aPosition - cPosition) >= 0.0f;
}

ElementIndex Triangles::FindContaining(
    vec2f const & position,
    Points const & points) const
{
    for (auto const t : *this)
    {
        if (ContainsPoint(position, t, points))
        {
            return t;
        }
    }

    return NoneElementIndex;
}

bcoords3f Triangles::ToBarycentricCoordinates(
    vec2f const & position,
    ElementIndex triangleElementIndex,
    Points const & points) const
{
    vec2f abBaryCoords = InternalToBarycentricCoordinates(
        position,
        triangleElementIndex,
        points);

    return bcoords3f(
        abBaryCoords.x,
        abBaryCoords.y,
        1.0f - abBaryCoords.x - abBaryCoords.y);
}

bcoords3f Triangles::ToBarycentricCoordinates(
    vec2f const & position,
    ElementIndex triangleElementIndex,
    Points const & points,
    float epsilon) const
{
    vec2f abBaryCoords = InternalToBarycentricCoordinates(
        position,
        triangleElementIndex,
        points);

    if (std::abs(abBaryCoords.x) < epsilon)
    {
        abBaryCoords.x = 0.0f;
    }

    if (std::abs(abBaryCoords.y) < epsilon)
    {
        abBaryCoords.y = 0.0f;
    }

    float z = 1.0f - abBaryCoords.x - abBaryCoords.y;
    if (std::abs(z) < epsilon)
    {
        z = 0.0f;
    }

    return bcoords3f(
        abBaryCoords.x,
        abBaryCoords.y,
        z);
}

bcoords3f Triangles::ToBarycentricCoordinatesFromWithinTriangle(
    vec2f const & position,
    ElementIndex triangleElementIndex,
    Points const & points) const
{
    assert(ContainsPoint(position, triangleElementIndex, points));

    vec2f abBaryCoords = InternalToBarycentricCoordinates(
        position,
        triangleElementIndex,
        points).clamp(0.0f, 1.0f, 0.0f, 1.0f);

    return bcoords3f(
        abBaryCoords.x,
        abBaryCoords.y,
        1.0f - abBaryCoords.x - abBaryCoords.y);
}

vec2f Triangles::FromBarycentricCoordinates(
    bcoords3f const & barycentricCoordinates,
    ElementIndex triangleElementIndex,
    Points const & points) const
{
    vec2f const & positionA = points.GetPosition(mEndpointsBuffer[triangleElementIndex].PointIndices[0]);
    vec2f const & positionB = points.GetPosition(mEndpointsBuffer[triangleElementIndex].PointIndices[1]);
    vec2f const & positionC = points.GetPosition(mEndpointsBuffer[triangleElementIndex].PointIndices[2]);

    return
        positionA * barycentricCoordinates[0]
        + positionB * barycentricCoordinates[1]
        + positionC * barycentricCoordinates[2];
}

//////////////////////////////////////////////////////////

vec2f Triangles::InternalToBarycentricCoordinates(
    vec2f const & position,
    ElementIndex triangleElementIndex,
    Points const & points) const
{
    vec2f const & positionA = points.GetPosition(mEndpointsBuffer[triangleElementIndex].PointIndices[0]);
    vec2f const & positionB = points.GetPosition(mEndpointsBuffer[triangleElementIndex].PointIndices[1]);
    vec2f const & positionC = points.GetPosition(mEndpointsBuffer[triangleElementIndex].PointIndices[2]);

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
        // See also: https://gamedev.stackexchange.com/questions/23743/whats-the-most-efficient-way-to-find-barycentric-coordinates

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

        return vec2f(
            l1,
            l2);
    }
}

}