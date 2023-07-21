/***************************************************************************************
* Original Author:		Gabriele Giuseppini
* Created:				2020-05-16
* Copyright:			Gabriele Giuseppini  (https://github.com/GabrieleGiuseppini)
***************************************************************************************/
#include "MeshBuilder.h"

#include "Log.h"

#include <cassert>
#include <utility>


//////////////////////////////////////////////////////////////////////////////

std::unique_ptr<Mesh> MeshBuilder::BuildMesh(std::filesystem::path const & meshDefinitionFilepath)
{
    // TODO
    (void)meshDefinitionFilepath;

    Vertices vertices(3);
    vertices.Add(vec2f(0.0f, 0.0f));
    vertices.Add(vec2f(0.0f, 1.0f));
    vertices.Add(vec2f(1.0f, 0.0f));

    Edges::TrianglesVector etv1;
    etv1.push_back(0);

    Edges edges(3);
    edges.Add(0, 1, Octant(6), Octant(2), etv1);
    edges.Add(1, 2, Octant(1), Octant(5), etv1);
    edges.Add(2, 0, Octant(4), Octant(0), etv1);

    Triangles triangles(1);
    triangles.Add(0, 1, 2, 0, 1, 2);

    vertices.AddConnectedEdge(0, 0, 1);
    vertices.AddConnectedEdge(0, 2, 2);
    vertices.AddConnectedEdge(1, 0, 0);
    vertices.AddConnectedEdge(1, 1, 2);
    vertices.AddConnectedEdge(2, 1, 1);
    vertices.AddConnectedEdge(2, 2, 0);

    vertices.AddConnectedTriangle(0, 0, true);
    vertices.AddConnectedTriangle(1, 0, false);
    vertices.AddConnectedTriangle(2, 0, false);

    auto mesh = std::make_unique<Mesh>(
        std::move(vertices),
        std::move(edges),
        std::move(triangles));

    return mesh;
}
