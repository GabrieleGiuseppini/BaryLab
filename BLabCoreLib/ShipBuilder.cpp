/***************************************************************************************
* Original Author:		Gabriele Giuseppini
* Created:				2020-05-16
* Copyright:			Gabriele Giuseppini  (https://github.com/GabrieleGiuseppini)
***************************************************************************************/
#include "ShipBuilder.h"

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

std::unique_ptr<Physics::Ship> ShipBuilder::BuildShip(
    ShipDefinition && shipDefinition,
    StructuralMaterialDatabase const & structuralMaterialDatabase)
{
    //
    // Process structural layer and:
    // - Create ShipBuildVertex's for each vertex
    // - Build a 2D matrix containing indices to the vertices
    //

    int const shipWidth = shipDefinition.StructuralLayerImage.Size.Width;
    float const halfShipWidth = static_cast<float>(shipWidth / 2); // We want to align on integral world coords

    int const shipHeight = shipDefinition.StructuralLayerImage.Size.Height;
    float const halfShipHeight = static_cast<float>(shipHeight / 2); // We want to align on integral world coords

    auto const & structuralLayerBuffer = shipDefinition.StructuralLayerImage.Data;

    // Vertices
    std::vector<ShipBuildVertex> vertexInfos;

    // Edges
    std::vector<ShipBuildEdge> edgeInfos;

    // Triangles
    std::vector<ShipBuildTriangle> triangleInfos;

    // Matrix of vertices - we allocate 2 extra dummy rows and cols - around - to avoid checking for boundaries
    ShipBuildVertexIndexMatrix vertexIndexMatrix(shipWidth + 2, shipHeight + 2);

    // Region of actual content
    int minX = shipWidth;
    int maxX = 0;
    int minY = shipHeight;
    int maxY = 0;

    // Visit all columns
    for (int x = 0; x < shipWidth; ++x)
    {
        // From bottom to top
        for (int y = 0; y < shipHeight; ++y)
        {
            StructuralMaterialDatabase::ColorKey const colorKey = structuralLayerBuffer[x + y * shipWidth];
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
                        static_cast<float>(x) - halfShipWidth,
                        static_cast<float>(y) - halfShipHeight),
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
    //  - Detect edges and create ShipBuildEdge's for them
    //      - And populate the point pair -> edge index 1 map
    //  - Do tessellation and create ShipBuildTriangle's
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
    // Visit all ShipBuildVertex's and create Vertices, i.e. the entire set of vertices
    //

    Physics::Vertices vertices = CreateVertices(
        vertexInfos);

    //
    // Create Edges for all ShipBuildEdge's
    //

    Physics::Edges edges = CreateEdges(
        edgeInfos,
        vertexInfos,
        vertices);

    //
    // Create Triangles for all ShipBuildTriangle's
    //

    Physics::Triangles triangles = CreateTriangles(
        triangleInfos,
        vertices,
        edges);

    //
    // We're done!
    //

    LogMessage("ShipBuilder: Created ship: W=", shipWidth, ", H=", shipHeight, ", ",
        vertices.GetBufferElementCount(), "buf vertices, ",
        edges.GetElementCount(), " edges, ",
        triangles.GetElementCount(), " triangles.");

    auto ship = std::make_unique<Physics::Ship>(
        std::move(vertices),
        std::move(edges),
        std::move(triangles));

    return ship;
}

void ShipBuilder::CreateElementInfos(
    ShipBuildVertexIndexMatrix const & vertexIndexMatrix,
    std::vector<ShipBuildVertex> & vertexInfos,
    std::vector<ShipBuildEdge> & edgeInfos,
    std::vector<ShipBuildTriangle> & triangleInfos)
{
    //
    // Visit vertex matrix and:
    //  - Detect edges and create ShipBuildEdge's for them
    //  - Do tessellation and create ShipBuildTriangle's
    //

    // From bottom to top - excluding extras at boundaries
    for (int y = 1; y < vertexIndexMatrix.height - 1; ++y)
    {
        // We're starting a new row, so we're not in a ship now
        bool isInShip = false;

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
                        // Create ShipBuildEdge
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
                        // - If this is the first vertex that is in a ship, we check all the way up to W;
                        // - Else, we check only up to S, so to avoid covering areas already covered by the triangulation
                        //   at the previous point
                        //

                        // Check adjacent point in next CW direction
                        int adjx2 = x + TessellationCircularOrderDirections[i + 1][0];
                        int adjy2 = y + TessellationCircularOrderDirections[i + 1][1];
                        if ((!isInShip || i < 2)
                            && !!vertexIndexMatrix[{adjx2, adjy2}])
                        {
                            // This vertex is adjacent to the first vertex at one of SE, S, SW, W

                            //
                            // Create ShipBuildTriangle
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
                            // Create ShipBuildTriangle
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

                // Remember now that we're in a ship
                isInShip = true;
            }
            else
            {
                //
                // No vertex exists at these coordinates
                //

                // From now on we're not in a ship anymore
                isInShip = false;
            }
        }
    }
}

void ShipBuilder::ConnectVerticesToTriangles(
    std::vector<ShipBuildVertex> & vertexInfos,
    std::vector<ShipBuildTriangle> const & triangleInfos)
{
    for (ElementIndex t = 0; t < triangleInfos.size(); ++t)
    {
        // Add triangle to its endpoints
        vertexInfos[triangleInfos[t].VertexIndices[0]].ConnectedTriangles.emplace_back(t);
        vertexInfos[triangleInfos[t].VertexIndices[1]].ConnectedTriangles.emplace_back(t);
        vertexInfos[triangleInfos[t].VertexIndices[2]].ConnectedTriangles.emplace_back(t);
    }
}

void ShipBuilder::ConnectEdgesToTriangles(
    std::vector<ShipBuildEdge> & edgeInfos,
    std::vector<ShipBuildTriangle> & triangleInfos)
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

Physics::Vertices ShipBuilder::CreateVertices(
    std::vector<ShipBuildVertex> const & vertexInfos)
{
    Physics::Vertices vertices(
        static_cast<ElementIndex>(vertexInfos.size()));

    for (size_t v = 0; v < vertexInfos.size(); ++v)
    {
        ShipBuildVertex const & vertexInfo = vertexInfos[v];

        //
        // Create vertex
        //

        vertices.Add(
            vertexInfo.Position);
    }

    return vertices;
}

Physics::Edges ShipBuilder::CreateEdges(
    std::vector<ShipBuildEdge> const & edgeInfos,
    std::vector<ShipBuildVertex> & vertexInfos,
    Physics::Vertices & vertices)
{
    Physics::Edges edges(static_cast<ElementIndex>(edgeInfos.size()));

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

Physics::Triangles ShipBuilder::CreateTriangles(
    std::vector<ShipBuildTriangle> const & triangleInfos,
    Physics::Vertices & vertices,
    Physics::Edges const & edges)
{
    Physics::Triangles triangles(static_cast<ElementIndex>(triangleInfos.size()));

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