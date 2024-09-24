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
    std::tuple<ElementIndex, int> subSpringAOppositeTriangleInfo,
    std::tuple<ElementIndex, int> subSpringBOppositeTriangleInfo,
    std::tuple<ElementIndex, int> subSpringCOppositeTriangleInfo,
    std::tuple<NpcFloorKindType, NpcFloorGeometryType> subSpringAFloorInfo,
    std::tuple<NpcFloorKindType, NpcFloorGeometryType> subSpringBFloorInfo,
    std::tuple<NpcFloorKindType, NpcFloorGeometryType> subSpringCFloorInfo,
    std::optional<ElementIndex> /*coveredTraverseSpringIndex*/)
{
    mEndpointsBuffer.emplace_back(pointAIndex, pointBIndex, pointCIndex);
    mSubSpringsBuffer.emplace_back(subSpringAIndex, subSpringBIndex, subSpringCIndex);
    mOppositeTrianglesBuffer.emplace_back(OppositeTrianglesInfo{
        OppositeTriangleInfo(std::get<0>(subSpringAOppositeTriangleInfo), std::get<1>(subSpringAOppositeTriangleInfo)),
        OppositeTriangleInfo(std::get<0>(subSpringBOppositeTriangleInfo), std::get<1>(subSpringBOppositeTriangleInfo)),
        OppositeTriangleInfo(std::get<0>(subSpringCOppositeTriangleInfo), std::get<1>(subSpringCOppositeTriangleInfo)) });
    mSubSpringNpcFloorKindsBuffer.emplace_back(SubSpringNpcFloorKinds({ std::get<0>(subSpringAFloorInfo), std::get<0>(subSpringBFloorInfo), std::get<0>(subSpringCFloorInfo) }));
    mSubSpringNpcFloorGeometriesBuffer.emplace_back(SubSpringNpcFloorGeometries({ std::get<1>(subSpringAFloorInfo), std::get<1>(subSpringBFloorInfo), std::get<1>(subSpringCFloorInfo) }));
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

        if (Geometry::IsPointInTriangle(position, aPosition, bPosition, cPosition))
        {
            return t;
        }
    }

    return NoneElementIndex;
}

}