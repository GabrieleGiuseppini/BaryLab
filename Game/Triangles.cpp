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
    NpcFloorKindType subSpringAFloorKind,
    NpcFloorGeometryType subSpringAFloorGeometry,
    NpcFloorKindType subSpringBFloorKind,
    NpcFloorGeometryType subSpringBFloorGeometry,
    NpcFloorKindType subSpringCFloorKind,
    NpcFloorGeometryType subSpringCFloorGeometry)
{
    mEndpointsBuffer.emplace_back(pointAIndex, pointBIndex, pointCIndex);
    mSubSpringsBuffer.emplace_back(subSpringAIndex, subSpringBIndex, subSpringCIndex);
    mOppositeTrianglesBuffer.emplace_back(OppositeTrianglesInfo{
        OppositeTriangleInfo(subSpringAOppositeTriangle, subSpringAOppositeTriangleEdgeOrdinal),
        OppositeTriangleInfo(subSpringBOppositeTriangle, subSpringBOppositeTriangleEdgeOrdinal),
        OppositeTriangleInfo(subSpringCOppositeTriangle, subSpringCOppositeTriangleEdgeOrdinal) });
    mSubSpringNpcFloorKindsBuffer.emplace_back(SubSpringNpcFloorKinds({ subSpringAFloorKind, subSpringBFloorKind, subSpringCFloorKind }));
    mSubSpringNpcFloorGeometriesBuffer.emplace_back(SubSpringNpcFloorGeometries({ subSpringAFloorGeometry, subSpringBFloorGeometry, subSpringCFloorGeometry }));
}

ElementIndex Triangles::FindContaining(
    vec2f const & position,
    Points const & points) const
{
    for (auto const t : *this)
    {
        vec2f const aPosition = points.GetPosition(GetPointAIndex(t));
        vec2f const bPosition = points.GetPosition(GetPointBIndex(t));
        vec2f const cPosition = points.GetPosition(GetPointCIndex(t));

        if (IsPointInTriangle(position, aPosition, bPosition, cPosition))
        {
            return t;
        }
    }

    return NoneElementIndex;
}

}