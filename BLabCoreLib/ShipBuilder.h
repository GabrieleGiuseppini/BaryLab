/***************************************************************************************
* Original Author:      Gabriele Giuseppini
* Created:              2020-05-16
* Copyright:            Gabriele Giuseppini  (https://github.com/GabrieleGiuseppini)
***************************************************************************************/
#pragma once

#include "BLabTypes.h"
#include "Physics.h"
#include "ShipBuilderTypes.h"
#include "ShipDefinition.h"
#include "StructuralMaterialDatabase.h"

#include <memory>

/*
 * This class contains all the logic for building a ship out of a Definition.
 */
class ShipBuilder
{
public:

    static std::unique_ptr<Physics::Ship> BuildShip(
        ShipDefinition && shipDefinition,
        StructuralMaterialDatabase const & structuralMaterialDatabase);

private:

    static void CreateElementInfos(
        ShipBuildVertexIndexMatrix const & vertexIndexMatrix,
        std::vector<ShipBuildVertex> & vertexInfos,
        std::vector<ShipBuildEdge> & edgeInfos,
        std::vector<ShipBuildTriangle> & triangleInfos);

    static void ConnectVerticesToTriangles(
        std::vector<ShipBuildVertex> & vertexInfos,
        std::vector<ShipBuildTriangle> const & triangleInfos);

    static void ConnectEdgesToTriangles(
        std::vector<ShipBuildEdge> & edgeInfos,
        std::vector<ShipBuildTriangle> & triangleInfos);

    static Physics::Vertices CreateVertices(
        std::vector<ShipBuildVertex> const & vertexInfos);

    static Physics::Edges CreateEdges(
        std::vector<ShipBuildEdge> const & edgeInfos,
        std::vector<ShipBuildVertex> & vertexInfos,
        Physics::Vertices & vertices);

    static Physics::Triangles CreateTriangles(
        std::vector<ShipBuildTriangle> const & triangleInfos,
        Physics::Vertices & vertices,
        std::vector<ShipBuildEdge> const & edgeInfos,
        Physics::Edges const & edges);
};
