/***************************************************************************************
* Original Author:      Gabriele Giuseppini
* Created:              2020-05-16
* Copyright:            Gabriele Giuseppini  (https://github.com/GabrieleGiuseppini)
***************************************************************************************/
#pragma once

#include "BLabTypes.h"
#include "Colors.h"
#include "FixedSizeVector.h"
#include "Matrix.h"
#include "StructuralMaterialDatabase.h"
#include "Vectors.h"

#include <memory>
#include <unordered_map>
#include <vector>

/*
 * Types describing the intermediate ship structure.
 */

using ShipBuildVertexIndexMatrix = Matrix2<std::optional<ElementIndex>>;

struct ShipBuildVertex
{
    vec2f Position;
    rgbColor RenderColor;
    StructuralMaterial const & Material;

    std::vector<ElementIndex> ConnectedEdges;
    std::vector<ElementIndex> ConnectedTriangles;

    ShipBuildVertex(
        vec2f position,
        rgbColor renderColor,
        StructuralMaterial const & material)
        : Position(position)
        , RenderColor(renderColor)
        , Material(material)
        , ConnectedEdges()
        , ConnectedTriangles()
    {
    }

    void AddConnectedEdge(ElementIndex edgeIndex)
    {
        assert(!ContainsConnectedEdge(edgeIndex));
        ConnectedEdges.push_back(edgeIndex);
    }

private:

    inline bool ContainsConnectedEdge(ElementIndex edgeIndex1) const
    {
        return std::find(
            ConnectedEdges.cbegin(),
            ConnectedEdges.cend(),
            edgeIndex1)
            != ConnectedEdges.cend();
    }
};

struct ShipBuildEdge
{
    ElementIndex VertexAIndex;
    uint32_t VertexAAngle;

    ElementIndex VertexBIndex;
    uint32_t VertexBAngle;

    FixedSizeVector<ElementIndex, 2> Triangles; // Triangles that have this spring as an edge

    ShipBuildEdge(
        ElementIndex vertexAIndex,
        uint32_t vertexAAngle,
        ElementIndex vertexBIndex,
        uint32_t vertexBAngle)
        : VertexAIndex(vertexAIndex)
        , VertexAAngle(vertexAAngle)
        , VertexBIndex(vertexBIndex)
        , VertexBAngle(vertexBAngle)
    {
    }
};

struct ShipBuildTriangle
{
    std::array<ElementIndex, 3> VertexIndices;

    FixedSizeVector<ElementIndex, 3> Edges;

    ShipBuildTriangle(
        std::array<ElementIndex, 3> const & vertexIndices)
        : VertexIndices(vertexIndices)
        , Edges()
    {
    }
};

// Utilities for navigating object's structure

struct VertexPair
{
    ElementIndex Endpoint1Index;
    ElementIndex Endpoint2Index;

    VertexPair()
        : Endpoint1Index(NoneElementIndex)
        , Endpoint2Index(NoneElementIndex)
    {}

    VertexPair(
        ElementIndex endpoint1Index,
        ElementIndex endpoint2Index)
        : Endpoint1Index(std::min(endpoint1Index, endpoint2Index))
        , Endpoint2Index(std::max(endpoint1Index, endpoint2Index))
    {}

    bool operator==(VertexPair const & other) const
    {
        return this->Endpoint1Index == other.Endpoint1Index
            && this->Endpoint2Index == other.Endpoint2Index;
    }

    struct Hasher
    {
        size_t operator()(VertexPair const & p) const
        {
            return p.Endpoint1Index * 23
                + p.Endpoint2Index;
        }
    };
};

using VertexPairToIndexMap = std::unordered_map<VertexPair, ElementIndex, VertexPair::Hasher>;
