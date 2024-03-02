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

// This is our local circular order (clockwise, starting from E)
// Note: cardinal directions are labeled according to y growing upwards
static int const TessellationCircularOrderDirections[8][2] = {
    {  1,  0 },  // 0: E
    {  1, -1 },  // 1: SE
    {  0, -1 },  // 2: S
    { -1, -1 },  // 3: SW
    { -1,  0 },  // 4: W
    { -1,  1 },  // 5: NW
    {  0,  1 },  // 6: N
    {  1,  1 }   // 7: NE
};

//////////////////////////////////////////////////////////////////////////////

rgbColor constexpr EmptyMaterialColorKey = rgbColor(255, 255, 255);

std::unique_ptr<Mesh> MeshBuilder::BuildMesh(
    MeshDefinition && meshDefinition,
    StructuralMaterialDatabase const & structuralMaterialDatabase)
{
    //
    // Process structural layer and:
    // - Create MeshBuildVertex's for each vertex
    // - Build a 2D matrix containing indices to the vertices
    //

    int const meshWidth = meshDefinition.StructuralLayerImage.Size.Width;
    float const halfMeshWidth = static_cast<float>(meshWidth / 2); // We want to align on integral world coords

    int const meshHeight = meshDefinition.StructuralLayerImage.Size.Height;
    float const halfMeshHeight = static_cast<float>(meshHeight / 2); // We want to align on integral world coords

    auto const & structuralLayerBuffer = meshDefinition.StructuralLayerImage.Data;

    // Vertices
    std::vector<MeshBuildVertex> vertexInfos;

    // Edges
    std::vector<MeshBuildEdge> edgeInfos;

    // Triangles
    std::vector<MeshBuildTriangle> triangleInfos;

    // Matrix of vertices - we allocate 2 extra dummy rows and cols - around - to avoid checking for boundaries
    MeshBuildVertexIndexMatrix vertexIndexMatrix(meshWidth + 2, meshHeight + 2);

    // Region of actual content
    int minX = meshWidth;
    int maxX = 0;
    int minY = meshHeight;
    int maxY = 0;

    // Visit all columns
    for (int x = 0; x < meshWidth; ++x)
    {
        // From bottom to top
        for (int y = 0; y < meshHeight; ++y)
        {
            StructuralMaterialDatabase::ColorKey const colorKey = structuralLayerBuffer[x + y * meshWidth];
            StructuralMaterial const * structuralMaterial = structuralMaterialDatabase.FindStructuralMaterial(colorKey);
            if (nullptr != structuralMaterial)
            {
                //
                // Make a vertex
                //

                ElementIndex const vertexIndex = static_cast<ElementIndex>(vertexInfos.size());

                vertexIndexMatrix[{x + 1, y + 1}] = static_cast<ElementIndex>(vertexIndex);

                vertexInfos.emplace_back(
                    vec2f(
                        static_cast<float>(x) - halfMeshWidth,
                        static_cast<float>(y) - halfMeshHeight),
                    colorKey,
                    *structuralMaterial);

                //
                // Update min/max coords
                //

                minX = std::min(minX, x);
                maxX = std::max(maxX, x);
                minY = std::min(minY, y);
                maxY = std::max(maxY, y);
            }
            else if (colorKey != EmptyMaterialColorKey)
            {
                throw BLabException("Pixel at coordinate (" + std::to_string(x) + ", " + std::to_string(y) + ") is not a recognized material");
            }
        }
    }


    //
    // Visit vertex matrix and:
    //  - Detect edges and create MeshBuildEdge's for them
    //      - And populate the point pair -> edge index 1 map
    //  - Do tessellation and create MeshBuildTriangle's
    //

    CreateElementInfos(
        vertexIndexMatrix,
        vertexInfos,
        edgeInfos,
        triangleInfos);


    //
    // Connect vertices to triangles
    //

    ConnectVerticesToTriangles(
        vertexInfos,
        triangleInfos);

    //
    // Associate all edges with their triangles
    //

    ConnectEdgesToTriangles(
        edgeInfos,
        triangleInfos);

    //
    // Visit all MeshBuildVertex's and create Vertices, i.e. the entire set of vertices
    //

    Vertices vertices = CreateVertices(
        vertexInfos);

    //
    // Create Edges for all MeshBuildEdge's
    //

    Edges edges = CreateEdges(
        edgeInfos,
        vertexInfos,
        vertices);

    //
    // Create Triangles for all MeshBuildTriangle's
    //

    Triangles triangles = CreateTriangles(
        triangleInfos,
        vertices,
        edges);

    //
    // We're done!
    //

    LogMessage("MeshBuilder: Created mesh: W=", meshWidth, ", H=", meshHeight, ", ",
        vertices.GetBufferElementCount(), "buf vertices, ",
        edges.GetElementCount(), " edges, ",
        triangles.GetElementCount(), " triangles.");

    auto mesh = std::make_unique<Mesh>(
        std::move(vertices),
        std::move(edges),
        std::move(triangles));

    return mesh;
}

void MeshBuilder::CreateElementInfos(
    MeshBuildVertexIndexMatrix const & vertexIndexMatrix,
    std::vector<MeshBuildVertex> & vertexInfos,
    std::vector<MeshBuildEdge> & edgeInfos,
    std::vector<MeshBuildTriangle> & triangleInfos)
{
    //
    // Visit vertex matrix and:
    //  - Detect edges and create MeshBuildEdge's for them
    //  - Do tessellation and create MeshBuildTriangle's
    //

    // From bottom to top - excluding extras at boundaries
    for (int y = 1; y < vertexIndexMatrix.height - 1; ++y)
    {
        // We're starting a new row, so we're not in a mesh now
        bool isInMesh = false;

        // From left to right - excluding extras at boundaries
        for (int x = 1; x < vertexIndexMatrix.width - 1; ++x)
        {
            if (!!vertexIndexMatrix[{x, y}])
            {
                //
                // A vertex exists at these coordinates
                //

                ElementIndex vertexIndex = *vertexIndexMatrix[{x, y}];

                //
                // Check if an edge exists
                //

                // First four directions out of 8: from 0 deg (+x) through to 225 deg (-x -y),
                // i.e. E, SE, S, SW - this covers each pair of vertices in each direction
                for (int i = 0; i < 4; ++i)
                {
                    int adjx1 = x + TessellationCircularOrderDirections[i][0];
                    int adjy1 = y + TessellationCircularOrderDirections[i][1];

                    if (!!vertexIndexMatrix[{adjx1, adjy1}])
                    {
                        // This vertex is adjacent to the first point at one of E, SE, S, SW

                        //
                        // Create MeshBuildEdge
                        //

                        ElementIndex const otherEndpointIndex = *vertexIndexMatrix[{adjx1, adjy1}];

                        // Add edge to edge infos
                        ElementIndex const edgeIndex = static_cast<ElementIndex>(edgeInfos.size());
                        edgeInfos.emplace_back(
                            vertexIndex,
                            i,
                            otherEndpointIndex,
                            (i + 4) % 8);

                        // Add the edge to its endpoints
                        vertexInfos[vertexIndex].AddConnectedEdge(edgeIndex);
                        vertexInfos[otherEndpointIndex].AddConnectedEdge(edgeIndex);


                        //
                        // Check if a triangle exists
                        // - If this is the first vertex that is in a mesh, we check all the way up to W;
                        // - Else, we check only up to S, so to avoid covering areas already covered by the triangulation
                        //   at the previous point
                        //

                        // Check adjacent point in next CW direction
                        int adjx2 = x + TessellationCircularOrderDirections[i + 1][0];
                        int adjy2 = y + TessellationCircularOrderDirections[i + 1][1];
                        if ((!isInMesh || i < 2)
                            && !!vertexIndexMatrix[{adjx2, adjy2}])
                        {
                            // This vertex is adjacent to the first vertex at one of SE, S, SW, W

                            //
                            // Create MeshBuildTriangle
                            //

                            triangleInfos.emplace_back(
                                std::array<ElementIndex, 3>( // Vertices are in CW order
                                {
                                    vertexIndex,
                                    otherEndpointIndex,
                                    * vertexIndexMatrix[{adjx2, adjy2}]
                                }));
                        }

                        // Now, we also want to check whether the single "irregular" triangle from this vertex exists,
                        // i.e. the triangle between this vertex, the vertex at its E, and the vertex at its
                        // S, in case there is no vertex at SE.
                        // We do this so that we can forget the entire W side for inner vertices and yet ensure
                        // full coverage of the area
                        if (i == 0
                            && !vertexIndexMatrix[{x + TessellationCircularOrderDirections[1][0], y + TessellationCircularOrderDirections[1][1]}]
                            && !!vertexIndexMatrix[{x + TessellationCircularOrderDirections[2][0], y + TessellationCircularOrderDirections[2][1]}])
                        {
                            // If we're here, the point at E exists
                            assert(!!vertexIndexMatrix[vec2i(x + TessellationCircularOrderDirections[0][0], y + TessellationCircularOrderDirections[0][1])]);

                            //
                            // Create MeshBuildTriangle
                            //

                            triangleInfos.emplace_back(
                                std::array<ElementIndex, 3>( // Vertices are in CW order
                                {
                                    vertexIndex,
                                    * vertexIndexMatrix[{x + TessellationCircularOrderDirections[0][0], y + TessellationCircularOrderDirections[0][1]}],
                                    * vertexIndexMatrix[{x + TessellationCircularOrderDirections[2][0], y + TessellationCircularOrderDirections[2][1]}]
                                }));
                        }
                    }
                }

                // Remember now that we're in a mesh
                isInMesh = true;
            }
            else
            {
                //
                // No vertex exists at these coordinates
                //

                // From now on we're not in a mesh anymore
                isInMesh = false;
            }
        }
    }
}

void MeshBuilder::ConnectVerticesToTriangles(
    std::vector<MeshBuildVertex> & vertexInfos,
    std::vector<MeshBuildTriangle> const & triangleInfos)
{
    for (ElementIndex t = 0; t < triangleInfos.size(); ++t)
    {
        // Add triangle to its endpoints
        vertexInfos[triangleInfos[t].VertexIndices[0]].ConnectedTriangles.emplace_back(t);
        vertexInfos[triangleInfos[t].VertexIndices[1]].ConnectedTriangles.emplace_back(t);
        vertexInfos[triangleInfos[t].VertexIndices[2]].ConnectedTriangles.emplace_back(t);
    }
}

void MeshBuilder::ConnectEdgesToTriangles(
    std::vector<MeshBuildEdge> & edgeInfos,
    std::vector<MeshBuildTriangle> & triangleInfos)
{
    //
    // 1. Build Vertex Pair -> Edge table
    //

    VertexPairToIndexMap vertexPairToEdgeMap;

    for (ElementIndex s = 0; s < edgeInfos.size(); ++s)
    {
        vertexPairToEdgeMap.emplace(
            std::piecewise_construct,
            std::forward_as_tuple(edgeInfos[s].VertexAIndex, edgeInfos[s].VertexBIndex),
            std::forward_as_tuple(s));
    }

    //
    // 2. Visit all triangles and connect them to their edges
    //

    for (ElementIndex t = 0; t < triangleInfos.size(); ++t)
    {
        for (size_t p = 0; p < triangleInfos[t].VertexIndices.size(); ++p)
        {
            ElementIndex const endpointIndex1 = triangleInfos[t].VertexIndices[p];

            ElementIndex const nextEndpointIndex1 =
                p < triangleInfos[t].VertexIndices.size() - 1
                ? triangleInfos[t].VertexIndices[p + 1]
                : triangleInfos[t].VertexIndices[0];

            // Lookup edge for this pair
            auto const edgeIt = vertexPairToEdgeMap.find({ endpointIndex1, nextEndpointIndex1 });
            assert(edgeIt != vertexPairToEdgeMap.end());

            ElementIndex const edgeIndex = edgeIt->second;

            // Tell this spring that it has this additional super triangle
            edgeInfos[edgeIndex].Triangles.push_back(t);
            assert(edgeInfos[edgeIndex].Triangles.size() <= 2);

            // Tell the triangle about this sub spring
            assert(!triangleInfos[t].Edges.contains(edgeIndex));
            triangleInfos[t].Edges.push_back(edgeIndex);
        }
    }
}

Vertices MeshBuilder::CreateVertices(
    std::vector<MeshBuildVertex> const & vertexInfos)
{
    Vertices vertices(
        static_cast<ElementIndex>(vertexInfos.size()));

    for (size_t v = 0; v < vertexInfos.size(); ++v)
    {
        MeshBuildVertex const & vertexInfo = vertexInfos[v];

        //
        // Create vertex
        //

        vertices.Add(
            vertexInfo.Position);
    }

    return vertices;
}

Edges MeshBuilder::CreateEdges(
    std::vector<MeshBuildEdge> const & edgeInfos,
    std::vector<MeshBuildVertex> & vertexInfos,
    Vertices & vertices)
{
    Edges edges(
        static_cast<ElementIndex>(edgeInfos.size()));

    for (ElementIndex e = 0; e < edgeInfos.size(); ++e)
    {
        // Determine surface type
        SurfaceType const surface =
            (vertexInfos[edgeInfos[e].VertexAIndex].Material.Surface == SurfaceType::Floor && vertexInfos[edgeInfos[e].VertexBIndex].Material.Surface == SurfaceType::Floor)
            ? SurfaceType::Floor : SurfaceType::Open;

        // Create edge
        edges.Add(
            edgeInfos[e].VertexAIndex,
            edgeInfos[e].VertexBIndex,
            edgeInfos[e].VertexAAngle,
            edgeInfos[e].VertexBAngle,
            surface,
            edgeInfos[e].Triangles);

        // Add edge to its endpoints
        vertices.AddConnectedEdge(
            edgeInfos[e].VertexAIndex,
            e,
            edgeInfos[e].VertexBIndex);
        vertices.AddConnectedEdge(
            edgeInfos[e].VertexBIndex,
            e,
            edgeInfos[e].VertexAIndex);
    }

    return edges;
}

Triangles MeshBuilder::CreateTriangles(
    std::vector<MeshBuildTriangle> const & triangleInfos,
    Vertices & vertices,
    Edges const & edges)
{
    Triangles triangles(static_cast<ElementIndex>(triangleInfos.size()));

    for (ElementIndex t = 0; t < triangleInfos.size(); ++t)
    {
        assert(triangleInfos[t].Edges.size() == 3);

        // Create triangle
        triangles.Add(
            triangleInfos[t].VertexIndices[0],
            triangleInfos[t].VertexIndices[1],
            triangleInfos[t].VertexIndices[2],
            triangleInfos[t].Edges[0],
            triangleInfos[t].Edges[1],
            triangleInfos[t].Edges[2],
            edges.GetSurfaceType(triangleInfos[t].Edges[0]),
            edges.GetSurfaceType(triangleInfos[t].Edges[1]),
            edges.GetSurfaceType(triangleInfos[t].Edges[2]));

        // Add triangle to its endpoints
        vertices.AddConnectedTriangle(triangleInfos[t].VertexIndices[0], t, true); // Owner
        vertices.AddConnectedTriangle(triangleInfos[t].VertexIndices[1], t, false); // Not owner
        vertices.AddConnectedTriangle(triangleInfos[t].VertexIndices[2], t, false); // Not owner
    }

    return triangles;
}