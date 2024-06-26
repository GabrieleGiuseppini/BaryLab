/***************************************************************************************
* Original Author:      Gabriele Giuseppini
* Created:              2020-05-16
* Copyright:            Gabriele Giuseppini  (https://github.com/GabrieleGiuseppini)
***************************************************************************************/
#pragma once

#include "Materials.h"

#include <GameCore/Colors.h>
#include <GameCore/FixedSizeVector.h>
#include <GameCore/GameTypes.h>
#include <GameCore/Matrix.h>
#include <GameCore/Vectors.h>

#include <cassert>
#include <memory>
#include <unordered_map>
#include <vector>

/*
 * Types describing the intermediate ship structure.
 */

using ShipFactoryPointIndexMatrix = Matrix2<std::optional<ElementIndex>>;

struct ShipFactoryPoint
{
    std::optional<ShipSpaceCoordinates> DefinitionCoordinates; // From any of the layers that provide points
    vec2f Position;
    rgbColor RenderColor;
    StructuralMaterial const & StructuralMtl;

    std::vector<ElementIndex> ConnectedSprings;
    std::vector<ElementIndex> ConnectedTriangles;

    ShipFactoryPoint(
        std::optional<ShipSpaceCoordinates> definitionCoordinates,
        vec2f position,
        rgbColor renderColor,
        StructuralMaterial const & structuralMtl)
        : DefinitionCoordinates(definitionCoordinates)
        , Position(position)
        , RenderColor(renderColor)
        , StructuralMtl(structuralMtl)
        , ConnectedSprings()
        , ConnectedTriangles()
    {
    }

    void AddConnectedSpring(ElementIndex springIndex)
    {
        assert(!ContainsConnectedSpring(springIndex));
        ConnectedSprings.push_back(springIndex);
    }

private:

    inline bool ContainsConnectedSpring(ElementIndex springIndex1) const
    {
        return std::find(
            ConnectedSprings.cbegin(),
            ConnectedSprings.cend(),
            springIndex1)
            != ConnectedSprings.cend();
    }
};

struct ShipFactorySpring
{
    ElementIndex PointAIndex;
    Octant PointAAngle; // Octant from the point of view of PointA

    ElementIndex PointBIndex;
    Octant PointBAngle; // Octant from the point of view of PointB

    FixedSizeVector<ElementIndex, 2> Triangles; // Triangles that have this spring as an edge

    ShipFactorySpring(
        ElementIndex pointAIndex,
        Octant pointAAngle,
        ElementIndex pointBIndex,
        Octant pointBAngle)
        : PointAIndex(pointAIndex)
        , PointAAngle(pointAAngle)
        , PointBIndex(pointBIndex)
        , PointBAngle(pointBAngle)
    {
    }
};

struct ShipFactoryTriangle
{
    // In CW order
    std::array<ElementIndex, 3> PointIndices;

    // In CW order
    FixedSizeVector<ElementIndex, 3> Springs;

    ShipFactoryTriangle(
        std::array<ElementIndex, 3> const & pointIndices)
        : PointIndices(pointIndices)
        , Springs()
    {
    }

    int FindSpringOrdinal(ElementIndex spring) const
    {
        int s = 0;
        for (; s < 3; ++s)
            if (Springs[s] == spring)
                break;

        assert(s < 3);
        return s;
    }
};

struct ShipFactoryPointPair
{
    ElementIndex Endpoint1Index;
    ElementIndex Endpoint2Index;

    ShipFactoryPointPair()
        : Endpoint1Index(NoneElementIndex)
        , Endpoint2Index(NoneElementIndex)
    {}

    ShipFactoryPointPair(
        ElementIndex endpoint1Index,
        ElementIndex endpoint2Index)
        : Endpoint1Index(std::min(endpoint1Index, endpoint2Index))
        , Endpoint2Index(std::max(endpoint1Index, endpoint2Index))
    {}

    bool operator==(ShipFactoryPointPair const & other) const
    {
        return this->Endpoint1Index == other.Endpoint1Index
            && this->Endpoint2Index == other.Endpoint2Index;
    }

    struct Hasher
    {
        size_t operator()(ShipFactoryPointPair const & p) const
        {
            return p.Endpoint1Index * 23
                + p.Endpoint2Index;
        }
    };
};

using ShipFactoryPointPairToIndexMap = std::unordered_map<ShipFactoryPointPair, ElementIndex, ShipFactoryPointPair::Hasher>;

struct ShipFactoryFloorInfo
{
    NpcFloorKindType FloorKind;
    NpcFloorGeometryType FloorGeometry;
    ElementIndex SpringIndex;

    ShipFactoryFloorInfo(
        NpcFloorKindType floorKind,
        NpcFloorGeometryType floorGeometry,
        ElementIndex springIndex)
        : FloorKind(floorKind)
        , FloorGeometry(floorGeometry)
        , SpringIndex(springIndex)
    {}
};

using ShipFactoryFloorPlan = std::unordered_map<ShipFactoryPointPair, ShipFactoryFloorInfo, ShipFactoryPointPair::Hasher>;
