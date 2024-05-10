/***************************************************************************************
 * Original Author:     Gabriele Giuseppini
 * Created:             2024-05-07
 * Copyright:           Gabriele Giuseppini  (https://github.com/GabrieleGiuseppini)
 ***************************************************************************************/
#include "ShipFloorplanizer.h"

#include <GameCore/Log.h>

#include <cassert>

std::vector<ShipFactoryFloor> ShipFloorplanizer::BuildFloorplan(
	ShipFactoryPointIndexMatrix const & pointIndexMatrix,
	std::vector<ShipFactoryPoint> const & pointInfos,
	std::vector<ShipFactorySpring> const & springInfos,
	std::vector<ShipFactoryTriangle> & triangleInfos) const
{
	//
	// Naive implementation: follows spring hullness
	//
	// Note: the floorplan will not necessarily match triangle edges
	// (as we might have suppressed triangles in a quad); this floorplan
	// will anyway be a *superset* of the floor edges
	//

	std::vector<ShipFactoryFloor> floors;
	floors.reserve(springInfos.size());

	//
	// Visit all "hull springs" that derive from structure (e.g. not rope springs)
	//

	for (size_t s = 0; s < springInfos.size(); ++s)
	{
		auto const & springInfo = springInfos[s];

		// Make sure it derives from structure
		if (pointInfos[springInfo.PointAIndex].DefinitionCoordinates.has_value()
			&& pointInfos[springInfo.PointBIndex].DefinitionCoordinates.has_value())
		{
			// Make sure it's a "hull spring" (i.e. that both endpoints are hull)
			if (AreBothEndpointsHull(springInfo, pointInfos))
			{
				// Make sure it's not redundant
				////// TODOTEST
				////(void)triangleInfos;
				////(void)pointIndexMatrix;
				if (!IsRedundantFloorSpring(
					static_cast<ElementIndex>(s),
					pointIndexMatrix,
					pointInfos,
					springInfos,
					triangleInfos))
				{
					NpcFloorType floorType;
					if (pointInfos[springInfo.PointAIndex].DefinitionCoordinates->x == pointInfos[springInfo.PointBIndex].DefinitionCoordinates->x)
					{
						// Horizontal
						assert(std::abs(pointInfos[springInfo.PointAIndex].DefinitionCoordinates->y - pointInfos[springInfo.PointBIndex].DefinitionCoordinates->y) == 1);
						floorType = NpcFloorType::FloorPlane1;
					}
					else if (pointInfos[springInfo.PointAIndex].DefinitionCoordinates->y == pointInfos[springInfo.PointBIndex].DefinitionCoordinates->y)
					{
						// Vertical
						assert(std::abs(pointInfos[springInfo.PointAIndex].DefinitionCoordinates->x - pointInfos[springInfo.PointBIndex].DefinitionCoordinates->x) == 1);
						floorType = NpcFloorType::FloorPlane1;
					}
					else
					{
						// Diagonal
						assert(std::abs(pointInfos[springInfo.PointAIndex].DefinitionCoordinates->x - pointInfos[springInfo.PointBIndex].DefinitionCoordinates->x) == 1);
						assert(std::abs(pointInfos[springInfo.PointAIndex].DefinitionCoordinates->y - pointInfos[springInfo.PointBIndex].DefinitionCoordinates->y) == 1);
						floorType = NpcFloorType::FloorPlane2;
					}

					floors.emplace_back(
						springInfo.PointAIndex,
						springInfo.PointBIndex,
						floorType);
				}
			}
		}
	}

	return floors;
}

//////////////////////////////////////////////////////////////////////////////////////////////////////////

bool ShipFloorplanizer::IsRedundantFloorSpring(
	ElementIndex s,
	ShipFactoryPointIndexMatrix const & pointIndexMatrix,
	std::vector<ShipFactoryPoint> const & pointInfos,
	std::vector<ShipFactorySpring> const & springInfos,
	std::vector<ShipFactoryTriangle> & triangleInfos) const
{
	//
	// S is a redundant floor spring (i.e. a spring that comes out of tessellation
	// but not really meant to be floor) if and only if:
	//	- It's a hull spring
	//	- It's connected to two triangles
	//	- One of the two triangles (T):
	//		- Has all hull sides
	//		- The other two sides (other than S) "continue" (to Psa' and Psb') via hull springs
	//	- The other triangle (T'):
	//		- Has no hull sides (other than S)
	//
	//                  Psb`
	//                   |
	//                   |
	//                   |
	//           Pc-----Psb
	//           | T'   /|
	//           |    /  |
	//           |  /    |
	//	         |/   T  |
	// Psa'-----Psa------Ps
	//

	auto const & springInfo = springInfos[s];

	// It's a hull spring
	assert(AreBothEndpointsHull(springInfo, pointInfos));

	// It's connected to two triangles
	if (springInfo.Triangles.size() != 2)
	{
		return false;
	}

	// The two triangles meet the criteria
	if ((
			IsRedundantFloorSpringTriangle1(triangleInfos[springInfo.Triangles[0]], s, pointIndexMatrix, pointInfos)
			&& IsRedundantFloorSpringTriangle2(triangleInfos[springInfo.Triangles[1]], pointInfos)
		)
		||
		(
			IsRedundantFloorSpringTriangle1(triangleInfos[springInfo.Triangles[1]], s, pointIndexMatrix, pointInfos)
			&& IsRedundantFloorSpringTriangle2(triangleInfos[springInfo.Triangles[0]], pointInfos)
		))
	{
		return true;
	}

	return false;
}

bool ShipFloorplanizer::IsRedundantFloorSpringTriangle1(
	ShipFactoryTriangle const & triangle,
	ElementIndex s,
	ShipFactoryPointIndexMatrix const & pointIndexMatrix,
	std::vector<ShipFactoryPoint> const & pointInfos) const
{
	//	- Has all hull sides

	size_t hullSideCount = CountTriangleHullSides(triangle, pointInfos);
	assert(hullSideCount >= 1 && hullSideCount <= 3);

	// Has all hull sides
	if (hullSideCount != 3)
	{
		return false;
	}

	//	- The other two sides (other than S) "continue" (to Psa' and Psb') via hull springs

	int const sOrdinal = triangle.GetSpringOrdinal(s);
	int const psVertexOrdinal = (sOrdinal + 2) % 3;

	if (!DoesTriangleEdgeExtendViaHullSpring(
		psVertexOrdinal,
		(sOrdinal + 1) % 3, // Psb
		triangle,
		pointIndexMatrix,
		pointInfos))
	{
		return false;
	}

	if (!DoesTriangleEdgeExtendViaHullSpring(
		psVertexOrdinal,
		sOrdinal, // Psa
		triangle,
		pointIndexMatrix,
		pointInfos))
	{
		return false;
	}

	return true;
}

bool ShipFloorplanizer::IsRedundantFloorSpringTriangle2(
	ShipFactoryTriangle const & triangle,
	std::vector<ShipFactoryPoint> const & pointInfos) const
{
	//	- Has no hull sides (other than S)

	size_t hullSideCount = CountTriangleHullSides(triangle, pointInfos);
	assert(hullSideCount >= 1 && hullSideCount <= 3);
	return hullSideCount == 1;
}

bool ShipFloorplanizer::DoesTriangleEdgeExtendViaHullSpring(
	int vertexOrdinal1,
	int vertexOrdinal2,
	ShipFactoryTriangle const & triangle,
	ShipFactoryPointIndexMatrix const & pointIndexMatrix,
	std::vector<ShipFactoryPoint> const & pointInfos) const
{
	ElementIndex p1 = triangle.PointIndices[vertexOrdinal1];
	assert(pointInfos[p1].DefinitionCoordinates.has_value());

	ElementIndex p2 = triangle.PointIndices[vertexOrdinal2];
	assert(pointInfos[p2].DefinitionCoordinates.has_value());

	ShipSpaceCoordinates const targetP =
		*(pointInfos[p2].DefinitionCoordinates)
		+ (*(pointInfos[p2].DefinitionCoordinates) - *(pointInfos[p1].DefinitionCoordinates));

	auto const extensionP = pointIndexMatrix[{targetP.x + 1, targetP.y + 1}];
	return
		extensionP.has_value()
		&& pointInfos[*extensionP].Material.IsHull;
}

size_t ShipFloorplanizer::CountTriangleHullSides(
	ShipFactoryTriangle const & triangle,
	std::vector<ShipFactoryPoint> const & pointInfos) const
{
	size_t count = 0;

	if (pointInfos[triangle.PointIndices[0]].Material.IsHull
		&& pointInfos[triangle.PointIndices[1]].Material.IsHull)
	{
		++count;
	}

	if (pointInfos[triangle.PointIndices[1]].Material.IsHull
		&& pointInfos[triangle.PointIndices[2]].Material.IsHull)
	{
		++count;
	}

	if (pointInfos[triangle.PointIndices[2]].Material.IsHull
		&& pointInfos[triangle.PointIndices[0]].Material.IsHull)
	{
		++count;
	}

	return count;
}

bool ShipFloorplanizer::AreBothEndpointsHull(
	ShipFactorySpring const & spring,
	std::vector<ShipFactoryPoint> const & pointInfos) const
{
	return pointInfos[spring.PointAIndex].Material.IsHull
		&& pointInfos[spring.PointBIndex].Material.IsHull;
}
