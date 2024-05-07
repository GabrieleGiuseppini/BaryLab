/***************************************************************************************
* Original Author:		Gabriele Giuseppini
* Created:				2020-05-16
* Copyright:			Gabriele Giuseppini  (https://github.com/GabrieleGiuseppini)
***************************************************************************************/
#include "ShipFactory.h"

#include "ShipFloorplanizer.h"

#include <GameCore/Log.h>

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

std::unique_ptr<Physics::Ship> ShipFactory::BuildShip(
    ShipDefinition && shipDefinition,
    MaterialDatabase const & materialDatabase)
{
    //
    // Process structural layer and:
    // - Create ShipFactoryPoint's for each point
    // - Build a 2D matrix containing indices to the points
    //

    int const shipWidth = shipDefinition.StructuralLayerImage.Size.width;
    float const halfShipWidth = static_cast<float>(shipWidth / 2); // We want to align on integral world coords

    int const shipHeight = shipDefinition.StructuralLayerImage.Size.height;
    float const halfShipHeight = static_cast<float>(shipHeight / 2); // We want to align on integral world coords

    auto const & structuralLayerBuffer = shipDefinition.StructuralLayerImage.Data;

    // Points
    std::vector<ShipFactoryPoint> pointInfos;

    // Springs
    std::vector<ShipFactorySpring> springInfos;

    // Triangles
    std::vector<ShipFactoryTriangle> triangleInfos;

    // Matrix of points - we allocate 2 extra dummy rows and cols - around - to avoid checking for boundaries
    ShipFactoryPointIndexMatrix pointIndexMatrix(shipWidth + 2, shipHeight + 2);

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
            MaterialDatabase::ColorKey const colorKey = structuralLayerBuffer[x + y * shipWidth];
            StructuralMaterial const * structuralMaterial = materialDatabase.FindStructuralMaterial(colorKey);
            if (nullptr != structuralMaterial)
            {
                //
                // Make a point
                //

                ShipSpaceCoordinates const coords = ShipSpaceCoordinates(x, y);

                ElementIndex const pointIndex = static_cast<ElementIndex>(pointInfos.size());

                pointIndexMatrix[{x + 1, y + 1}] = static_cast<ElementIndex>(pointIndex);

                pointInfos.emplace_back(
                    coords,
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
                throw GameException("Pixel at coordinate (" + std::to_string(x) + ", " + std::to_string(y) + ") is not a recognized material");
            }
        }
    }


    //
    // Visit point matrix and:
    //  - Detect edges and create ShipFactorySpring's for them
    //      - And populate the point pair -> edge index 1 map
    //  - Do tessellation and create ShipFactoryTriangle's
    //

    CreateElementInfos(
        pointIndexMatrix,
        pointInfos,
        springInfos,
        triangleInfos);


    //
    // Connect points to triangles
    //

    ConnectPointsToTriangles(
        pointInfos,
        triangleInfos);

    //
    // Associate all edges with their triangles
    //

    ConnectSpringsToTriangles(
        springInfos,
        triangleInfos);

    //
    // Create floorplan
    //

    ShipFloorplanizer shipFloorplanizer;

    std::vector<ShipFactoryFloor> floorInfos = shipFloorplanizer.BuildFloorplan(
        pointInfos,
        springInfos);

    //
    // Visit all ShipFactoryPoint's and create Points, i.e. the entire set of points
    //

    Physics::Points points = CreatePoints(
        pointInfos);

    //
    // Create Springs for all ShipFactorySpring's
    //

    Physics::Springs springs = CreateSprings(
        springInfos,
        points);

    //
    // Create Triangles for all ShipFactoryTriangle's
    //

    Physics::Triangles triangles = CreateTriangles(
        triangleInfos,
        points,
        springInfos,
        floorInfos);

    //
    // We're done!
    //

    LogMessage("ShipFactory: Created ship: W=", shipWidth, ", H=", shipHeight, ", ",
        points.GetBufferElementCount(), "buf vertices, ",
        springs.GetElementCount(), " edges, ",
        triangles.GetElementCount(), " triangles.");

    auto ship = std::make_unique<Physics::Ship>(
        0, // First and only ship in BaryLab :-(
        std::move(points),
        std::move(springs),
        std::move(triangles));

    return ship;
}

void ShipFactory::CreateElementInfos(
    ShipFactoryPointIndexMatrix const & pointIndexMatrix,
    std::vector<ShipFactoryPoint> & pointInfos,
    std::vector<ShipFactorySpring> & springInfos,
    std::vector<ShipFactoryTriangle> & triangleInfos)
{
    //
    // Visit point matrix and:
    //  - Detect edges and create ShipFactorySpring's for them
    //  - Do tessellation and create ShipFactoryTriangle's
    //

    // From bottom to top - excluding extras at boundaries
    for (int y = 1; y < pointIndexMatrix.height - 1; ++y)
    {
        // We're starting a new row, so we're not in a ship now
        bool isInShip = false;

        // From left to right - excluding extras at boundaries
        for (int x = 1; x < pointIndexMatrix.width - 1; ++x)
        {
            if (!!pointIndexMatrix[{x, y}])
            {
                //
                // A point exists at these coordinates
                //

                ElementIndex pointIndex = *pointIndexMatrix[{x, y}];

                //
                // Check if an edge exists
                //

                // First four directions out of 8: from 0 deg (+x) through to 225 deg (-x -y),
                // i.e. E, SE, S, SW - this covers each pair of points in each direction
                for (int i = 0; i < 4; ++i)
                {
                    int adjx1 = x + TessellationCircularOrderDirections[i][0];
                    int adjy1 = y + TessellationCircularOrderDirections[i][1];

                    if (!!pointIndexMatrix[{adjx1, adjy1}])
                    {
                        // This point is adjacent to the first point at one of E, SE, S, SW

                        //
                        // Create ShipFactorySpring
                        //

                        ElementIndex const otherEndpointIndex = *pointIndexMatrix[{adjx1, adjy1}];

                        // Add spring to spring infos
                        ElementIndex const springIndex = static_cast<ElementIndex>(springInfos.size());
                        springInfos.emplace_back(
                            pointIndex,
                            i,
                            otherEndpointIndex,
                            (i + 4) % 8);

                        // Add the spring to its endpoints
                        pointInfos[pointIndex].AddConnectedSpring(springIndex);
                        pointInfos[otherEndpointIndex].AddConnectedSpring(springIndex);


                        //
                        // Check if a triangle exists
                        // - If this is the first point that is in a ship, we check all the way up to W;
                        // - Else, we check only up to S, so to avoid covering areas already covered by the triangulation
                        //   at the previous point
                        //

                        // Check adjacent point in next CW direction
                        int adjx2 = x + TessellationCircularOrderDirections[i + 1][0];
                        int adjy2 = y + TessellationCircularOrderDirections[i + 1][1];
                        if ((!isInShip || i < 2)
                            && !!pointIndexMatrix[{adjx2, adjy2}])
                        {
                            // This point is adjacent to the first point at one of SE, S, SW, W

                            //
                            // Create ShipFactoryTriangle
                            //

                            triangleInfos.emplace_back(
                                std::array<ElementIndex, 3>( // Points are in CW order
                                {
                                    pointIndex,
                                    otherEndpointIndex,
                                    * pointIndexMatrix[{adjx2, adjy2}]
                                }));
                        }

                        // Now, we also want to check whether the single "irregular" triangle from this point exists,
                        // i.e. the triangle between this point, the point at its E, and the point at its
                        // S, in case there is no point at SE.
                        // We do this so that we can forget the entire W side for inner points and yet ensure
                        // full coverage of the area
                        if (i == 0
                            && !pointIndexMatrix[{x + TessellationCircularOrderDirections[1][0], y + TessellationCircularOrderDirections[1][1]}]
                            && !!pointIndexMatrix[{x + TessellationCircularOrderDirections[2][0], y + TessellationCircularOrderDirections[2][1]}])
                        {
                            // If we're here, the point at E exists
                            assert(!!pointIndexMatrix[vec2i(x + TessellationCircularOrderDirections[0][0], y + TessellationCircularOrderDirections[0][1])]);

                            //
                            // Create ShipFactoryTriangle
                            //

                            triangleInfos.emplace_back(
                                std::array<ElementIndex, 3>( // Points are in CW order
                                {
                                    pointIndex,
                                    * pointIndexMatrix[{x + TessellationCircularOrderDirections[0][0], y + TessellationCircularOrderDirections[0][1]}],
                                    * pointIndexMatrix[{x + TessellationCircularOrderDirections[2][0], y + TessellationCircularOrderDirections[2][1]}]
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
                // No point exists at these coordinates
                //

                // From now on we're not in a ship anymore
                isInShip = false;
            }
        }
    }
}

void ShipFactory::ConnectPointsToTriangles(
    std::vector<ShipFactoryPoint> & pointInfos,
    std::vector<ShipFactoryTriangle> const & triangleInfos)
{
    for (ElementIndex t = 0; t < triangleInfos.size(); ++t)
    {
        // Add triangle to its endpoints
        pointInfos[triangleInfos[t].PointIndices[0]].ConnectedTriangles.emplace_back(t);
        pointInfos[triangleInfos[t].PointIndices[1]].ConnectedTriangles.emplace_back(t);
        pointInfos[triangleInfos[t].PointIndices[2]].ConnectedTriangles.emplace_back(t);
    }
}

void ShipFactory::ConnectSpringsToTriangles(
    std::vector<ShipFactorySpring> & springInfos,
    std::vector<ShipFactoryTriangle> & triangleInfos)
{
    //
    // 1. Build Point Pair -> Spring table
    //

    PointPairToIndexMap pointPairToSpringMap;

    for (ElementIndex s = 0; s < springInfos.size(); ++s)
    {
        pointPairToSpringMap.emplace(
            std::piecewise_construct,
            std::forward_as_tuple(springInfos[s].PointAIndex, springInfos[s].PointBIndex),
            std::forward_as_tuple(s));
    }

    //
    // 2. Visit all triangles and connect them to their springs
    //

    for (ElementIndex t = 0; t < triangleInfos.size(); ++t)
    {
        for (size_t p = 0; p < triangleInfos[t].PointIndices.size(); ++p)
        {
            ElementIndex const endpointIndex1 = triangleInfos[t].PointIndices[p];

            ElementIndex const nextEndpointIndex1 =
                p < triangleInfos[t].PointIndices.size() - 1
                ? triangleInfos[t].PointIndices[p + 1]
                : triangleInfos[t].PointIndices[0];

            // Lookup spring for this pair
            auto const springIt = pointPairToSpringMap.find({ endpointIndex1, nextEndpointIndex1 });
            assert(springIt != pointPairToSpringMap.end());

            ElementIndex const springIndex = springIt->second;

            // Tell this spring that it has this additional super triangle
            springInfos[springIndex].Triangles.push_back(t);
            assert(springInfos[springIndex].Triangles.size() <= 2);

            // Tell the triangle about this sub spring
            assert(!triangleInfos[t].Springs.contains(springIndex));
            triangleInfos[t].Springs.push_back(springIndex);
        }
    }
}

Physics::Points ShipFactory::CreatePoints(
    std::vector<ShipFactoryPoint> const & pointInfos)
{
    Physics::Points points(
        static_cast<ElementIndex>(pointInfos.size()));

    for (size_t v = 0; v < pointInfos.size(); ++v)
    {
        ShipFactoryPoint const & pointInfo = pointInfos[v];

        //
        // Create point
        //

        points.Add(
            pointInfo.Position);
    }

    return points;
}

Physics::Springs ShipFactory::CreateSprings(
    std::vector<ShipFactorySpring> const & springInfos,
    Physics::Points & points)
{
    Physics::Springs springs(static_cast<ElementIndex>(springInfos.size()));

    for (ElementIndex e = 0; e < springInfos.size(); ++e)
    {
        // Create spring
        springs.Add(
            springInfos[e].PointAIndex,
            springInfos[e].PointBIndex,
            springInfos[e].PointAAngle,
            springInfos[e].PointBAngle,
            springInfos[e].Triangles);

        // Add spring to its endpoints
        points.AddConnectedSpring(
            springInfos[e].PointAIndex,
            e,
            springInfos[e].PointBIndex);
        points.AddConnectedSpring(
            springInfos[e].PointBIndex,
            e,
            springInfos[e].PointAIndex);
    }

    return springs;
}

Physics::Triangles ShipFactory::CreateTriangles(
    std::vector<ShipFactoryTriangle> const & triangleInfos,
    Physics::Points & points,
    std::vector<ShipFactorySpring> const & springInfos,
    std::vector<ShipFactoryFloor> const & floorInfos)
{
    //
    // Prepare floorplan: map of pairs of points that are a floor
    //

    std::unordered_map<PointPair, NpcFloorType, PointPair::Hasher> floorPointPairMap;
    for (auto const & factoryFloor : floorInfos)
    {
        auto const [_, isInserted] = floorPointPairMap.try_emplace(
            { factoryFloor.Endpoint1Index, factoryFloor.Endpoint2Index },
            factoryFloor.FloorType);

        assert(isInserted);
        (void)isInserted;
    }

    //
    // Build triangles
    //

    Physics::Triangles triangles(static_cast<ElementIndex>(triangleInfos.size()));

    for (ElementIndex t = 0; t < triangleInfos.size(); ++t)
    {
        assert(triangleInfos[t].Springs.size() == 3);

        // Derive whether this is a sealed triangle
        bool isSealedTriangle = true;
        for (int iEdge = 0; iEdge < 3; ++iEdge)
        {
            ElementIndex const springElementIndex = triangleInfos[t].Springs[iEdge];

            ElementIndex const edgePointA = springInfos[springElementIndex].PointAIndex;
            ElementIndex const edgePointB = springInfos[springElementIndex].PointBIndex;

            //isSealedTriangle = isSealedTriangle && (floorPointPairMap.find({ edgePointA, edgePointB }) != floorPointPairMap.end());
            isSealedTriangle = isSealedTriangle && (floorPointPairMap.find({ edgePointA, edgePointB }) != floorPointPairMap.end());
        }

        // Calculate opposite triangles and floor types
        std::array<std::pair<ElementIndex, int>, 3> subSpringsOppositeTriangle;
        std::array<NpcFloorType, 3> subSpringsFloorType;
        for (int iEdge = 0; iEdge < 3; ++iEdge)
        {
            ElementIndex const springElementIndex = triangleInfos[t].Springs[iEdge];
            assert(springInfos[springElementIndex].Triangles.size() >= 1 && springInfos[springElementIndex].Triangles.size() <= 2);

            if (springInfos[springElementIndex].Triangles[0] == t)
            {
                if (springInfos[springElementIndex].Triangles.size() >= 2)
                {
                    assert(springInfos[springElementIndex].Triangles.size() == 2);
                    subSpringsOppositeTriangle[iEdge].first = springInfos[springElementIndex].Triangles[1];
                }
                else
                {
                    subSpringsOppositeTriangle[iEdge].first = NoneElementIndex;
                }
            }
            else if (springInfos[springElementIndex].Triangles.size() >= 2)
            {
                assert(springInfos[springElementIndex].Triangles.size() == 2);
                assert(springInfos[springElementIndex].Triangles[1] == t);
                subSpringsOppositeTriangle[iEdge].first = springInfos[springElementIndex].Triangles[0];
            }
            else
            {
                subSpringsOppositeTriangle[iEdge].first = NoneElementIndex;
            }

            if (subSpringsOppositeTriangle[iEdge].first != NoneElementIndex)
            {
                if (triangleInfos[subSpringsOppositeTriangle[iEdge].first].Springs[0] == triangleInfos[t].Springs[iEdge])
                {
                    subSpringsOppositeTriangle[iEdge].second = 0;
                }
                else if (triangleInfos[subSpringsOppositeTriangle[iEdge].first].Springs[1] == triangleInfos[t].Springs[iEdge])
                {
                    subSpringsOppositeTriangle[iEdge].second = 1;
                }
                else
                {
                    assert(triangleInfos[subSpringsOppositeTriangle[iEdge].first].Springs[2] == triangleInfos[t].Springs[iEdge]);
                    subSpringsOppositeTriangle[iEdge].second = 2;
                }
            }

            //
            // Triangle's subedge is floor if:
            //  - Spring is floor, AND
            //  - NOT is sealed, OR (is sealed and) there's no triangle on the other side of this subedge
            //

            ElementIndex const edgePointA = springInfos[springElementIndex].PointAIndex;
            ElementIndex const edgePointB = springInfos[springElementIndex].PointBIndex;

            NpcFloorType floorType;
            if (const auto floorIt = floorPointPairMap.find({ edgePointA, edgePointB });
                floorIt != floorPointPairMap.cend()
                && (!isSealedTriangle || subSpringsOppositeTriangle[iEdge].first == NoneElementIndex))
            {
                floorType = floorIt->second;
            }
            else
            {
                floorType = NpcFloorType::Open;
            }

            subSpringsFloorType[iEdge] = floorType;
        }

        // Create triangle
        triangles.Add(
            triangleInfos[t].PointIndices[0],
            triangleInfos[t].PointIndices[1],
            triangleInfos[t].PointIndices[2],
            triangleInfos[t].Springs[0],
            triangleInfos[t].Springs[1],
            triangleInfos[t].Springs[2],
            subSpringsOppositeTriangle[0].first,
            subSpringsOppositeTriangle[0].second,
            subSpringsOppositeTriangle[1].first,
            subSpringsOppositeTriangle[1].second,
            subSpringsOppositeTriangle[2].first,
            subSpringsOppositeTriangle[2].second,
            subSpringsFloorType[0],
            subSpringsFloorType[1],
            subSpringsFloorType[2]);

        // Add triangle to its endpoints
        points.AddConnectedTriangle(triangleInfos[t].PointIndices[0], t, true); // Owner
        points.AddConnectedTriangle(triangleInfos[t].PointIndices[1], t, false); // Not owner
        points.AddConnectedTriangle(triangleInfos[t].PointIndices[2], t, false); // Not owner
    }

    return triangles;
}