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

	//
	// 1. Build first floorplan using all and only hull springs
	//

	ShipFactoryFloorPlan hullSprings;
	hullSprings.reserve(springInfos.size());

	for (size_t s = 0; s < springInfos.size(); ++s)
	{
		auto const & springInfo = springInfos[s];

		// Make sure it derives from structure
		if (pointInfos[springInfo.PointAIndex].DefinitionCoordinates.has_value()
			&& pointInfos[springInfo.PointBIndex].DefinitionCoordinates.has_value())
		{
			// Make sure it's a "hull spring" (i.e. that both endpoints are hull)
			if (pointInfos[springInfo.PointAIndex].Material.IsHull
				&& pointInfos[springInfo.PointBIndex].Material.IsHull)
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

				auto const [_, isInserted] = hullSprings.try_emplace(
					{ springInfo.PointAIndex, springInfo.PointBIndex },
					floorType,
					static_cast<ElementIndex>(s));

				assert(isInserted);
				(void)isInserted;
			}
		}
	}

	//
	// 2. Do a first pass removing redundant springs
	//

	auto intermediateFloorPlan = RemoveRedundantFloors(
		std::move(hullSprings),
		pointIndexMatrix,
		pointInfos,
		springInfos,
		triangleInfos);

	//
	// 3. Do a last pass removing redundant springs
	//

	auto floorPlan = RemoveRedundantFloors(
		std::move(intermediateFloorPlan),
		pointIndexMatrix,
		pointInfos,
		springInfos,
		triangleInfos);

	return floorPlan;
}

//////////////////////////////////////////////////////////////////////////////////////////////////////////

ShipFactoryFloorPlan ShipFloorplanizer::RemoveRedundantFloors(
	ShipFactoryFloorPlan && hullSprings,
	ShipFactoryPointIndexMatrix const & pointIndexMatrix,
	std::vector<ShipFactoryPoint> const & pointInfos,
	std::vector<ShipFactorySpring> const & springInfos,
	std::vector<ShipFactoryTriangle> & triangleInfos) const
{
	ShipFactoryFloorPlan floorPlan;
	floorPlan.reserve(springInfos.size());

	//
	// Visit all "hull springs"
	//

	for (auto const & hullSpringEntry : hullSprings)
	{
		auto const & springInfo = springInfos[hullSpringEntry.second.SpringIndex];

		// TODOTEST
		if (hullSpringEntry.second.SpringIndex == 22)
			LogMessage("TODOTEST");

		// Make sure it's not a redundant spring
		if (!IsRedundantFloorSpring(
			hullSpringEntry.second.SpringIndex,
			hullSprings,
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
				hullSpringEntry.first,
				hullSpringEntry.second);

			assert(isInserted);
			(void)isInserted;
		}
	}

	return floorPlan;
}

bool ShipFloorplanizer::IsRedundantFloorSpring(
	ElementIndex s,
	ShipFactoryFloorPlan const & hullSprings,
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
	assert(hullSprings.count({ springInfo.PointAIndex, springInfo.PointBIndex }) == 1);

	// It's connected to two triangles
	if (springInfo.Triangles.size() != 2)
	{
		return false;
	}

	// The two triangles meet the criteria
	if ((
			IsRedundantFloorSpringTriangle1(s, springInfo, triangleInfos[springInfo.Triangles[0]], hullSprings, pointIndexMatrix, pointInfos)
			&& IsRedundantFloorSpringTriangle2(triangleInfos[springInfo.Triangles[1]], hullSprings)
		)
		||
		(
			IsRedundantFloorSpringTriangle1(s, springInfo, triangleInfos[springInfo.Triangles[1]], hullSprings, pointIndexMatrix, pointInfos)
			&& IsRedundantFloorSpringTriangle2(triangleInfos[springInfo.Triangles[0]], hullSprings)
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
	ShipFactoryFloorPlan const & hullSprings,
	ShipFactoryPointIndexMatrix const & pointIndexMatrix,
	std::vector<ShipFactoryPoint> const & pointInfos) const
{
	//	- Has all hull sides

	size_t hullSideCount = CountTriangleHullSides(triangle, hullSprings);
	assert(hullSideCount >= 1 && hullSideCount <= 3);

	// Has all hull sides
	if (hullSideCount != 3)
	{
		return false;
	}

	// - From both of the two vertices of S we encounter (in CW or CCW order, depending on vertex) hull springs at "steep" slopes(delta octant neither 4 nor 5)

	int const sOrdinal = triangle.FindSpringOrdinal(s);

	if (!HasTriangleSteepHullExtensions(
		triangle.PointIndices[sOrdinal], // Psa
		spring,
		7, // CCW
		hullSprings,
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
		hullSprings,
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
	ShipFactoryFloorPlan const & hullSprings) const
{
	//	- Has no hull sides (other than S)

	size_t hullSideCount = CountTriangleHullSides(triangle, hullSprings);
	assert(hullSideCount >= 1 && hullSideCount <= 3);
	return hullSideCount == 1;
}

bool ShipFloorplanizer::HasTriangleSteepHullExtensions(
	ElementIndex startPoint,
	ShipFactorySpring const & spring,
	Octant direction,
	ShipFactoryFloorPlan const & hullSprings,
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
		hullSprings,
		pointIndexMatrix,
		pointInfos);

	return deltaOctant != 4 && deltaOctant != 5;
}

Octant ShipFloorplanizer::FindNextHullSpringOctant(
	ElementIndex startPoint,
	Octant startOctant,
	Octant direction,
	ShipFactoryFloorPlan const & hullSprings,
	ShipFactoryPointIndexMatrix const & pointIndexMatrix,
	std::vector<ShipFactoryPoint> const & pointInfos) const
{
	// TODOHERE: cleanup
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
			&& hullSprings.count({startPoint, *pointIndexMatrix[destCoords] }) > 0)
		{
			// Found hull spring
			break;
		}
	}

	return distanceOctant;
}

size_t ShipFloorplanizer::CountTriangleHullSides(
	ShipFactoryTriangle const & triangle,
	ShipFactoryFloorPlan const & hullSprings) const
{
	size_t const count =
		hullSprings.count({ triangle.PointIndices[0], triangle.PointIndices[1] })
		+ hullSprings.count({ triangle.PointIndices[1], triangle.PointIndices[2] })
		+ hullSprings.count({ triangle.PointIndices[2], triangle.PointIndices[0] });

	return count;
}
