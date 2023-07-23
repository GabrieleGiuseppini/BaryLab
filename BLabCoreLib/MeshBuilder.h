/***************************************************************************************
* Original Author:      Gabriele Giuseppini
* Created:              2020-05-16
* Copyright:            Gabriele Giuseppini  (https://github.com/GabrieleGiuseppini)
***************************************************************************************/
#pragma once

#include "BLabTypes.h"
#include "Edges.h"
#include "Mesh.h"
#include "MeshBuilderTypes.h"
#include "MeshDefinition.h"
#include "StructuralMaterialDatabase.h"
#include "Triangles.h"
#include "Vertices.h"

#include <memory>

/*
 * This class contains all the logic for building a mesh out of a Definition.
 */
class MeshBuilder
{
public:

    static std::unique_ptr<Mesh> BuildMesh(
        MeshDefinition && meshDefinition,
        StructuralMaterialDatabase const & structuralMaterialDatabase);

private:

    static void CreateElementInfos(
        MeshBuildVertexIndexMatrix const & vertexIndexMatrix,
        std::vector<MeshBuildVertex> & vertexInfos,
        std::vector<MeshBuildEdge> & edgeInfos,
        std::vector<MeshBuildTriangle> & triangleInfos);

    static void ConnectVerticesToTriangles(
        std::vector<MeshBuildVertex> & vertexInfos,
        std::vector<MeshBuildTriangle> const & triangleInfos);

    static void ConnectEdgesToTriangles(
        std::vector<MeshBuildEdge> & edgeInfos,
        std::vector<MeshBuildTriangle> & triangleInfos);

    static Vertices CreateVertices(
        std::vector<MeshBuildVertex> const & vertexInfos);

    static Edges CreateEdges(
        std::vector<MeshBuildEdge> const & edgeInfos,
        std::vector<MeshBuildVertex> & vertexInfos,
        Vertices & vertices);

    static Triangles CreateTriangles(
        std::vector<MeshBuildTriangle> const & triangleInfos,
        Vertices & vertices);
};
