/***************************************************************************************
 * Original Author:     Gabriele Giuseppini
 * Created:             2024-05-07
 * Copyright:           Gabriele Giuseppini  (https://github.com/GabrieleGiuseppini)
 ***************************************************************************************/
#include "ShipFloorplanizer.h"

#include <GameCore/Log.h>

#include <cassert>

ShipFactoryFloorPlan ShipFloorplanizer::BuildFloorplan(
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

	ShipFactoryFloorPlan floorPlan;
	floorPlan.reserve(springInfos.size());

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
				if (s == 22)
				{
					LogMessage("HERE");
				}

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

					auto const [_, isInserted] = floorPlan.try_emplace(
						{ springInfo.PointAIndex, springInfo.PointBIndex },
						floorType);

					assert(isInserted);
					(void)isInserted;
				}
			}
		}
	}

	return floorPlan;
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
	//		- From both of the two vertices of S we encounter (in CW or CCW order, depending on vertex) hull springs at "steep" slopes (delta octant neither 4 nor 5)
	//	- The other triangle (T'):
	//		- Has no hull sides (other than S)
	//
	//                  Psb`    /Psb'
	//                   |    /
	//                   |  /
	//                   |/
	//           Pc-----Psb
	//           | T'   /|
	//           |   S/  |
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
			IsRedundantFloorSpringTriangle1(s, springInfo, triangleInfos[springInfo.Triangles[0]], pointIndexMatrix, pointInfos)
			&& IsRedundantFloorSpringTriangle2(triangleInfos[springInfo.Triangles[1]], pointInfos)
		)
		||
		(
			IsRedundantFloorSpringTriangle1(s, springInfo, triangleInfos[springInfo.Triangles[1]], pointIndexMatrix, pointInfos)
			&& IsRedundantFloorSpringTriangle2(triangleInfos[springInfo.Triangles[0]], pointInfos)
		))
	{
		return true;
	}

	return false;
}

bool ShipFloorplanizer::IsRedundantFloorSpringTriangle1(
	ElementIndex s,
	ShipFactorySpring const & spring,
	ShipFactoryTriangle const & triangle,
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

	// - From both of the two vertices of S we encounter(in CW or CCW order, depending on vertex) hull springs at "steep" slopes(delta octant neither 4 nor 5)

	int const sOrdinal = triangle.FindSpringOrdinal(s);

	if (!HasTriangleSteepHullExtensions(
		triangle.PointIndices[sOrdinal], // Psa
		spring,
		7, // CCW
		pointIndexMatrix,
		pointInfos))
	{
		// One smooth path
		return false;
	}

	if (!HasTriangleSteepHullExtensions(
		triangle.PointIndices[(sOrdinal + 1) % 3], // Psb
		spring,
		1, // CW
		pointIndexMatrix,
		pointInfos))
	{
		// One smooth path
		return false;
	}

	// No smooth paths
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

bool ShipFloorplanizer::HasTriangleSteepHullExtensions(
	ElementIndex startPoint,
	ShipFactorySpring const & spring,
	Octant direction,
	ShipFactoryPointIndexMatrix const & pointIndexMatrix,
	std::vector<ShipFactoryPoint> const & pointInfos) const
{
	Octant startOctant;
	if (startPoint == spring.PointAIndex)
	{
		startOctant = spring.PointAAngle;
	}
	else
	{
		assert(startPoint == spring.PointBIndex);
		startOctant = spring.PointBAngle;
	}

	Octant const deltaOctant = FindNextHullSpringOctant(
		startPoint,
		startOctant,
		direction,
		pointIndexMatrix,
		pointInfos);

	return deltaOctant != 4 && deltaOctant != 5;
}

Octant ShipFloorplanizer::FindNextHullSpringOctant(
	ElementIndex startPoint,
	Octant startOctant,
	Octant direction,
	ShipFactoryPointIndexMatrix const & pointIndexMatrix,
	std::vector<ShipFactoryPoint> const & pointInfos) const
{
	assert(pointInfos[startPoint].DefinitionCoordinates.has_value());
	auto const startCoords = *(pointInfos[startPoint].DefinitionCoordinates);
	Octant distanceOctant = 1;
	Octant deltaOctant = direction;
	for (; ; deltaOctant = (deltaOctant + direction) % 8, ++distanceOctant)
	{
		auto const destCoords = vec2i(
			startCoords.x + 1 + TessellationCircularOrderDirections[(startOctant + deltaOctant) % 8][0],
			startCoords.y + 1 + TessellationCircularOrderDirections[(startOctant + deltaOctant) % 8][1]);
		if (pointIndexMatrix[destCoords].has_value()
			&& pointInfos[*pointIndexMatrix[destCoords]].Material.IsHull)
		{
			// Found hull spring
			break;
		}
	}

	return distanceOctant;
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
