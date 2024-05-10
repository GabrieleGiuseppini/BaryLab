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
    StructuralMaterial const & Material;

    std::vector<ElementIndex> ConnectedSprings;
    std::vector<ElementIndex> ConnectedTriangles;

    ShipFactoryPoint(
        std::optional<ShipSpaceCoordinates> definitionCoordinates,
        vec2f position,
        rgbColor renderColor,
        StructuralMaterial const & material)
        : DefinitionCoordinates(definitionCoordinates)
        , Position(position)
        , RenderColor(renderColor)
        , Material(material)
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
    uint32_t PointAAngle;

    ElementIndex PointBIndex;
    uint32_t PointBAngle;

    FixedSizeVector<ElementIndex, 2> Triangles; // Triangles that have this spring as an edge

    ShipFactorySpring(
        ElementIndex pointAIndex,
        uint32_t pointAAngle,
        ElementIndex pointBIndex,
        uint32_t pointBAngle)
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

    int GetSpringOrdinal(ElementIndex spring) const
    {
        int s = 0;
        for (; s < 3; ++s)
            if (Springs[s] == spring)
                break;

        assert(s < 3);
        return s;
    }
};

struct ShipFactoryFloor
{
    ElementIndex Endpoint1Index;
    ElementIndex Endpoint2Index;
    NpcFloorType FloorType;

    ShipFactoryFloor(
        ElementIndex endpoint1Index,
        ElementIndex endpoint2Index,
        NpcFloorType floorType)
        : Endpoint1Index(endpoint1Index)
        , Endpoint2Index(endpoint2Index)
        , FloorType(floorType)
    {
    }
};

// Utilities for navigating object's structure

struct PointPair
{
    ElementIndex Endpoint1Index;
    ElementIndex Endpoint2Index;

    PointPair()
        : Endpoint1Index(NoneElementIndex)
        , Endpoint2Index(NoneElementIndex)
    {}

    PointPair(
        ElementIndex endpoint1Index,
        ElementIndex endpoint2Index)
        : Endpoint1Index(std::min(endpoint1Index, endpoint2Index))
        , Endpoint2Index(std::max(endpoint1Index, endpoint2Index))
    {}

    bool operator==(PointPair const & other) const
    {
        return this->Endpoint1Index == other.Endpoint1Index
            && this->Endpoint2Index == other.Endpoint2Index;
    }

    struct Hasher
    {
        size_t operator()(PointPair const & p) const
        {
            return p.Endpoint1Index * 23
                + p.Endpoint2Index;
        }
    };
};

using PointPairToIndexMap = std::unordered_map<PointPair, ElementIndex, PointPair::Hasher>;
