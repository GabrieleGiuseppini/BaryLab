/***************************************************************************************
 * Original Author:     Gabriele Giuseppini
 * Created:             2024-05-07
 * Copyright:           Gabriele Giuseppini  (https://github.com/GabrieleGiuseppini)
 ***************************************************************************************/
#include "ShipFloorplanizer.h"

#include "GameParameters.h"

#include <GameCore/Log.h>

#include <cassert>

ShipFactoryFloorPlan ShipFloorplanizer::BuildFloorplan(
	ShipFactoryPointIndexMatrix const & pointIndexMatrix,
	std::vector<ShipFactoryPoint> const & pointInfos,
	std::vector<ShipFactorySpring> const & springInfos,
	std::vector<ShipFactoryTriangle> & triangleInfos) const
{
	//
	// 1. Build list of springs that we do not want to use as floors;
	//    we do so by detecting specific vertex patterns in 3x3 blocks
	//

	SpringExclusionSet springExclusionSet;

	// Process all 3x3 blocks - including the 1-wide "borders"
	VertexBlock vertexBlock(IntegralRectSize(3, 3));
	for (int y = 0; y < pointIndexMatrix.height - 2; ++y)
	{
		for (int x = 0; x < pointIndexMatrix.width - 2; ++x)
		{
			// TODOTEST
			if (x == 5 && y == 9)
				LogMessage("TODOTEST");

			// Build block
			for (int yb = 0; yb < 3; ++yb)
			{
				for (int xb = 0; xb < 3; ++xb)
				{
					if (pointIndexMatrix[{x + xb, y + yb}].has_value()
						&& pointInfos[*pointIndexMatrix[{x + xb, y + yb}]].Material.IsHull)
					{
						vertexBlock[{xb, yb}] = *pointIndexMatrix[{x + xb, y + yb}];
					}
					else
					{
						vertexBlock[{xb, yb}] = NoneElementIndex;
					}
				}
			}

			// Process block
			ProcessVertexBlock(
				vertexBlock,
				springExclusionSet);
		}
	}

	//
	// 2. Build floorplan with All and ONLY "hull" springs which:
	//	- Are directly derived from structure, and
	//	- Are on the side of a triangle, and
	//	- Are not in the exclusione set
	//

	ShipFactoryFloorPlan floorPlan;
	floorPlan.reserve(springInfos.size());

	for (size_t s = 0; s < springInfos.size(); ++s)
	{
		auto const & springInfo = springInfos[s];

		////// TODOTEST
		////if (s == 23)
		////	LogMessage("TODOTEST");

		// Make sure it's viable as a floor and, if it's a non-external edge, it's not in the exclusion list
		if (IsSpringViableForFloor(springInfo, pointInfos)
			&& (springInfo.Triangles.size() == 1 || springExclusionSet.count({springInfo.PointAIndex, springInfo.PointBIndex}) == 0))
		{
			//
			// Take this spring
			//

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
				floorType,
				static_cast<ElementIndex>(s));

			assert(isInserted);
			(void)isInserted;
		}
	}

	// TODOTEST
	(void)triangleInfos;

	return floorPlan;
}

bool ShipFloorplanizer::IsSpringViableForFloor(
	ShipFactorySpring const & springInfo,
	std::vector<ShipFactoryPoint> const & pointInfos) const
{
	return
		pointInfos[springInfo.PointAIndex].DefinitionCoordinates.has_value() // Is point derived directly from structure?
		&& pointInfos[springInfo.PointAIndex].Material.IsHull // Is point hull?
		&& pointInfos[springInfo.PointBIndex].DefinitionCoordinates.has_value() // Is point derived directly from structure?
		&& pointInfos[springInfo.PointBIndex].Material.IsHull // Is point hull?
		&& springInfo.Triangles.size() > 0 // Is it an edge of a triangle?
		;
}

void ShipFloorplanizer::ProcessVertexBlock(
	VertexBlock & vertexBlock,
	/*out*/ SpringExclusionSet & springExclusionSet) const
{
	// 1. All rotations of symmetry 1

	for (int i = 0; i < 4; ++i)
	{
		ProcessVertexBlockPatterns(
			vertexBlock,
			springExclusionSet);

		vertexBlock.Rotate90(RotationDirectionType::Clockwise);

	}

	// 2. All rotations of symmetry 2

	vertexBlock.Flip(DirectionType::Vertical);

	for (int i = 0; i < 4; ++i)
	{
		ProcessVertexBlockPatterns(
			vertexBlock,
			springExclusionSet);

		vertexBlock.Rotate90(RotationDirectionType::Clockwise);

	}
}

void ShipFloorplanizer::ProcessVertexBlockPatterns(
	VertexBlock const & vertexBlock,
	/*out*/ SpringExclusionSet & springExclusionSet) const
{
	//
	// Check for a set of specific patterns; once one is found,
	// exclude specific springs (which might not even exist)
	//

	//
	// Pattern 1: "under a stair" (_\): take care of redundant /
	//
	//   *?o
	//   o*?
	//   ***
	//

	if (vertexBlock[{0, 0}] != NoneElementIndex && vertexBlock[{1, 0}] != NoneElementIndex && vertexBlock[{2, 0}] != NoneElementIndex
		&& vertexBlock[{0, 1}] == NoneElementIndex && vertexBlock[{1, 1}] != NoneElementIndex
		&& vertexBlock[{0, 2}] != NoneElementIndex && vertexBlock[{2, 2}] == NoneElementIndex)
	{
		springExclusionSet.insert({ vertexBlock[{0, 0}] , vertexBlock[{1, 1}] });
	}

	//
	// Pattern 2: "under a stair" (_\): take care of redundant |
	//
	//   *oo
	//   o*?
	//   ***
	//

	if (vertexBlock[{0, 0}] != NoneElementIndex && vertexBlock[{1, 0}] != NoneElementIndex && vertexBlock[{2, 0}] != NoneElementIndex
		&& vertexBlock[{0, 1}] == NoneElementIndex && vertexBlock[{1, 1}] != NoneElementIndex
		&& vertexBlock[{0, 2}] != NoneElementIndex && vertexBlock[{1, 2}] == NoneElementIndex && vertexBlock[{2, 2}] == NoneElementIndex)
	{
		springExclusionSet.insert({ vertexBlock[{1, 0}] , vertexBlock[{1, 1}] });
	}

	//
	// Pattern 3: "wall-on-floor" (|-): take care of redundant \ and /
	//
	//  o?o
	//  o*o
	//  ***
	//

	if (vertexBlock[{0, 0}] != NoneElementIndex && vertexBlock[{1, 0}] != NoneElementIndex && vertexBlock[{2, 0}] != NoneElementIndex
		&& vertexBlock[{0, 1}] == NoneElementIndex && vertexBlock[{1, 1}] != NoneElementIndex && vertexBlock[{2, 1}] == NoneElementIndex
		&& vertexBlock[{0, 2}] == NoneElementIndex && vertexBlock[{2, 2}] == NoneElementIndex)
	{
		springExclusionSet.insert({ vertexBlock[{0, 0}] , vertexBlock[{1, 1}] });
		springExclusionSet.insert({ vertexBlock[{2, 0}] , vertexBlock[{1, 1}] });
	}

	//
	// Pattern 4: "corner" (|_): take care of redundant \
	//
	//  *o?
	//  *o?
	//  ***
	//

	if (vertexBlock[{0, 0}] != NoneElementIndex && vertexBlock[{1, 0}] != NoneElementIndex && vertexBlock[{2, 0}] != NoneElementIndex
		&& vertexBlock[{0, 1}] != NoneElementIndex && vertexBlock[{1, 1}] == NoneElementIndex
		&& vertexBlock[{0, 2}] != NoneElementIndex && vertexBlock[{1, 2}] == NoneElementIndex)
	{
		springExclusionSet.insert({ vertexBlock[{0, 1}] , vertexBlock[{1, 0}] });
	}

	//
	// Pattern 5: "TODO" (_||): take care of redundant /
	//
	//  o**
	//  o**
	//  ***
	//

	if (vertexBlock[{0, 0}] != NoneElementIndex && vertexBlock[{1, 0}] != NoneElementIndex && vertexBlock[{2, 0}] != NoneElementIndex
		&& vertexBlock[{0, 1}] == NoneElementIndex && vertexBlock[{1, 1}] != NoneElementIndex && vertexBlock[{2, 1}] != NoneElementIndex
		&& vertexBlock[{0, 2}] == NoneElementIndex && vertexBlock[{1, 2}] != NoneElementIndex && vertexBlock[{2, 2}] != NoneElementIndex)
	{
		springExclusionSet.insert({ vertexBlock[{0, 0}] , vertexBlock[{1, 1}] });
		springExclusionSet.insert({ vertexBlock[{2, 0}] , vertexBlock[{1, 1}] });
	}

}

// TODOOLD

////ShipFactoryFloorPlan ShipFloorplanizer::BuildFloorplan(
////	ShipFactoryPointIndexMatrix const & pointIndexMatrix,
////	std::vector<ShipFactoryPoint> const & pointInfos,
////	std::vector<ShipFactorySpring> const & springInfos,
////	std::vector<ShipFactoryTriangle> & triangleInfos) const
////{
////	//
////	// Floorplan is:
////	//	- All and only "hull" springs - directly derived from structure - which are on the side of a triangle
////	//	AND is NOT:
////	//	- One single spring long (i.e. there is nothing before or after them along their exact direction), and
////	//	- Incident on "long" floors (i.e. with two aligned floors on both sides of both endpoints)
////	//		- So we keep traverse one-spring in U patterns
////	//
////
////	ShipFactoryFloorPlan floorPlan;
////	floorPlan.reserve(springInfos.size());
////
////	for (size_t s = 0; s < springInfos.size(); ++s)
////	{
////		auto const & springInfo = springInfos[s];
////
////		// TODOTEST
////		if (s == 23)
////			LogMessage("TODOTEST");
////
////		// Make sure it's viable as a floor
////		if (IsSpringViableForFloor(springInfo, pointInfos, triangleInfos))
////		{
////			// Make sure it's not a one-spring floor, and if it is,
////			// it's not incident on long floors
////			if (!IsIsolatedFloor(springInfo, pointInfos, springInfos, triangleInfos)
////				|| !IsFloorTwiceIncidentOnLongFloors(springInfo, static_cast<ElementIndex>(s), pointInfos, springInfos, triangleInfos))
////			{
////				//
////				// Take this hull spring
////				//
////
////				NpcFloorType floorType;
////				if (pointInfos[springInfo.PointAIndex].DefinitionCoordinates->x == pointInfos[springInfo.PointBIndex].DefinitionCoordinates->x)
////				{
////					// Horizontal
////					assert(std::abs(pointInfos[springInfo.PointAIndex].DefinitionCoordinates->y - pointInfos[springInfo.PointBIndex].DefinitionCoordinates->y) == 1);
////					floorType = NpcFloorType::FloorPlane1;
////				}
////				else if (pointInfos[springInfo.PointAIndex].DefinitionCoordinates->y == pointInfos[springInfo.PointBIndex].DefinitionCoordinates->y)
////				{
////					// Vertical
////					assert(std::abs(pointInfos[springInfo.PointAIndex].DefinitionCoordinates->x - pointInfos[springInfo.PointBIndex].DefinitionCoordinates->x) == 1);
////					floorType = NpcFloorType::FloorPlane1;
////				}
////				else
////				{
////					// Diagonal
////					assert(std::abs(pointInfos[springInfo.PointAIndex].DefinitionCoordinates->x - pointInfos[springInfo.PointBIndex].DefinitionCoordinates->x) == 1);
////					assert(std::abs(pointInfos[springInfo.PointAIndex].DefinitionCoordinates->y - pointInfos[springInfo.PointBIndex].DefinitionCoordinates->y) == 1);
////					floorType = NpcFloorType::FloorPlane2;
////				}
////
////				auto const [_, isInserted] = floorPlan.try_emplace(
////					{ springInfo.PointAIndex, springInfo.PointBIndex },
////					floorType,
////					static_cast<ElementIndex>(s));
////
////				assert(isInserted);
////				(void)isInserted;
////			}
////		}
////	}
////
////	// TODOTEST
////	(void)pointIndexMatrix;
////
////	return floorPlan;
////}

////////////////////////////////

////bool ShipFloorplanizer::IsSealedTriangle(
////	ElementIndex t,
////	std::vector<ShipFactoryPoint> const & pointInfos,
////	std::vector<ShipFactoryTriangle> & triangleInfos) const
////{
////	for (int v = 0; v < 3; ++v)
////	{
////		ElementIndex const vertexIndex = triangleInfos[t].PointIndices[v];
////		if (!pointInfos[vertexIndex].Material.IsHull)
////		{
////			return false;
////		}
////	}
////
////	return true;
////}
////
////bool ShipFloorplanizer::IsIsolatedFloor(
////	ShipFactorySpring const & spring,
////	std::vector<ShipFactoryPoint> const & pointInfos,
////	std::vector<ShipFactorySpring> const & springInfos,
////	std::vector<ShipFactoryTriangle> & triangleInfos) const
////{
////	//
////	// Check whether it's a one-long floor
////	//
////
////	// Before A
////	if (DoesFloorContinue(
////		spring.PointAIndex,
////		spring.PointBAngle,
////		pointInfos,
////		springInfos,
////		triangleInfos))
////	{
////		return false;
////	}
////
////	// After B
////	if (DoesFloorContinue(
////		spring.PointBIndex,
////		spring.PointAAngle,
////		pointInfos,
////		springInfos,
////		triangleInfos))
////	{
////		return false;
////	}
////
////	return true;
////}
////
////bool ShipFloorplanizer::DoesFloorContinue(
////	ElementIndex pointIndex,
////	Octant direction,
////	std::vector<ShipFactoryPoint> const & pointInfos,
////	std::vector<ShipFactorySpring> const & springInfos,
////	std::vector<ShipFactoryTriangle> & triangleInfos) const
////{
////	//
////	// Check whether the spring continues from the specified point in the specified direction
////	//
////
////	assert(pointInfos[pointIndex].DefinitionCoordinates.has_value());
////	assert(pointInfos[pointIndex].Material.IsHull);
////
////	// Check whether there is a spring in the specified direction
////	auto const floorContinuationSpring = FindSpringAtPointAndDirection(
////		pointIndex,
////		direction,
////		pointInfos,
////		springInfos);
////
////	if (!floorContinuationSpring.has_value())
////	{
////		return false;
////	}
////
////	// Check whether this spring we've found is a viable spring
////	if (!IsSpringViableForFloor(
////			springInfos[*floorContinuationSpring],
////			pointInfos,
////			triangleInfos))
////	{
////		return false;
////	}
////
////	return true;
////}
////
////bool ShipFloorplanizer::IsFloorTwiceIncidentOnLongFloors(
////	ShipFactorySpring const & springInfo,
////	ElementIndex springIndex,
////	std::vector<ShipFactoryPoint> const & pointInfos,
////	std::vector<ShipFactorySpring> const & springInfos,
////	std::vector<ShipFactoryTriangle> & triangleInfos) const
////{
////	if (!IsFloorIncidentOnLongFloors(
////		springInfo.PointAIndex,
////		springIndex,
////		pointInfos,
////		springInfos,
////		triangleInfos))
////	{
////		return false;
////	}
////
////	if (!IsFloorIncidentOnLongFloors(
////		springInfo.PointBIndex,
////		springIndex,
////		pointInfos,
////		springInfos,
////		triangleInfos))
////	{
////		return false;
////	}
////
////	return true;
////}
////
////bool ShipFloorplanizer::IsFloorIncidentOnLongFloors(
////	ElementIndex endpointIndex,
////	ElementIndex floorSpringIndex,
////	std::vector<ShipFactoryPoint> const & pointInfos,
////	std::vector<ShipFactorySpring> const & springInfos,
////	std::vector<ShipFactoryTriangle> & triangleInfos) const
////{
////	// TODOTEST
////	(void)triangleInfos;
////
////	assert(pointInfos[endpointIndex].DefinitionCoordinates.has_value());
////	assert(pointInfos[endpointIndex].Material.IsHull);
////
////	// Find all other floors incident to this point
////	ElementIndex otherFloorSprings[GameParameters::MaxSpringsPerPoint];
////	int floorSpringCount = 0;
////	for (auto const s : pointInfos[endpointIndex].ConnectedSprings)
////	{
////		if (s != floorSpringIndex
////			&& IsSpringViableForFloor(springInfos[s], pointInfos, triangleInfos))
////		{
////			otherFloorSprings[floorSpringCount] = s;
////			++floorSpringCount;
////		}
////	}
////
////	// A bit fuzzy here: at this moment we consider floors that are incident *ONLY* to
////	// a long floor (thus with only 2 hull springs connected to the endpoint)
////	if (floorSpringCount != 2)
////	{
////		return false;
////	}
////
////	// Check if the second spring is in the same direction as the first spring - but
////	// on the other side of the endpoint
////
////	auto const & firstSpringInfo = springInfos[otherFloorSprings[0]];
////	Octant checkDirection;
////	if (firstSpringInfo.PointAIndex == endpointIndex)
////	{
////		checkDirection = firstSpringInfo.PointBAngle;
////	}
////	else
////	{
////		assert(firstSpringInfo.PointBIndex == endpointIndex);
////		checkDirection = firstSpringInfo.PointAAngle;
////	}
////
////	auto const & secondSpringInfo = springInfos[otherFloorSprings[1]];
////	Octant secondSpringDirection; // From the point of view of endpointIndex
////	if (secondSpringInfo.PointAIndex == endpointIndex)
////	{
////		secondSpringDirection = firstSpringInfo.PointAAngle;
////	}
////	else
////	{
////		assert(secondSpringInfo.PointBIndex == endpointIndex);
////		secondSpringDirection = firstSpringInfo.PointBAngle;
////	}
////
////	if (checkDirection != secondSpringDirection)
////	{
////		// Not aligned
////		return false;
////	}
////
////	return true;
////
////	////// TODOOLD
////	////// TODOHERE: hull springs only
////	////if (pointInfos[endpointIndex].ConnectedSprings.size() != 3
////	////	|| !IsSpringViableForFloor(springInfos[pointInfos[endpointIndex].ConnectedSprings[0]], pointInfos)
////	////	|| !IsSpringViableForFloor(springInfos[pointInfos[endpointIndex].ConnectedSprings[0]], pointInfos)
////	////{
////	////	return false;
////	////}
////
////	////// Find indices of springs that are not self
////	////int otherFloorSpringIndices[2] = { -1, -1 };
////	////for (int search = 0, target = 0; search < 3; ++search)
////	////{
////	////	if (pointInfos[endpointIndex].ConnectedSprings[search] != floorSpringIndex)
////	////	{
////	////		otherFloorSpringIndices[target] = pointInfos[endpointIndex].ConnectedSprings[search];
////	////		++target;
////	////	}
////	////}
////	////assert(otherFloorSpringIndices[0] >= 0 && otherFloorSpringIndices[1] >= 0);
////
////	////ElementIndex firstSpringOtherEndpoint;
////	////Octant secondSpringCheckDirection; // From endpointIndex
////	////ShipFactorySpring const & firstSpring = springInfos[otherFloorSpringIndices[0]];
////	////if (endpointIndex == firstSpring.PointAIndex)
////	////{
////	////	firstSpringOtherEndpoint = firstSpring.PointBIndex;
////	////	secondSpringCheckDirection = firstSpring.PointBAngle;
////	////}
////	////else
////	////{
////	////	assert(endpointIndex == firstSpring.PointBIndex);
////	////	firstSpringOtherEndpoint = firstSpring.PointAIndex;
////	////	secondSpringCheckDirection = firstSpring.PointAAngle;
////	////}
////
////	////// Make sure first spring is a floor
////	////if (!IsViableSpringEndpointForFloor(firstSpringOtherEndpoint, firstSpring, pointIndexMatrix, pointInfos))
////	////{
////	////	return false;
////	////}
////
////	////// Check if the second spring is in the same direction as the first spring - but
////	////// on the other side of the endpoint
////	////// TODOHERE: use FindSpringAtPointAndDirection
////
////	////// Check if the second spring is viable
////	////// TODOHERE: IsViable
////
////	////(void)floorSpring;
////	////return false;
////}
////
////std::optional<ElementIndex> ShipFloorplanizer::FindSpringAtPointAndDirection(
////	ElementIndex pointIndex,
////	Octant direction,
////	std::vector<ShipFactoryPoint> const & pointInfos,
////	std::vector<ShipFactorySpring> const & springInfos) const
////{
////	// Visit all springs at this point
////	for (ElementIndex s : pointInfos[pointIndex].ConnectedSprings)
////	{
////		auto const & springInfo = springInfos[s];
////		if ((springInfo.PointAIndex == pointIndex && springInfo.PointAAngle == direction)
////			|| (springInfo.PointBIndex == pointIndex && springInfo.PointBAngle == direction))
////		{
////			return s;
////		}
////	}
////
////	return std::nullopt;
////}
////
////////bool ShipFloorplanizer::IsSpringViableForFloor(
////////	ShipFactorySpring const & springInfo,
////////	std::vector<ShipFactoryPoint> const & pointInfos,
////////	std::vector<ShipFactoryTriangle> & triangleInfos) const
////////{
////////	return
////////		pointInfos[springInfo.PointAIndex].DefinitionCoordinates.has_value() // Is point derived directly from structure?
////////		&& pointInfos[springInfo.PointAIndex].Material.IsHull // Is point hull?
////////		&& pointInfos[springInfo.PointBIndex].DefinitionCoordinates.has_value() // Is point derived directly from structure?
////////		&& pointInfos[springInfo.PointBIndex].Material.IsHull // Is point hull?
////////		&& (springInfo.Triangles.size() == 1 // Is it a non-internal edge of a triangle?
////////			|| (springInfo.Triangles.size() == 2 &&
////////				(!IsSealedTriangle(springInfo.Triangles[0], pointInfos, triangleInfos) || !IsSealedTriangle(springInfo.Triangles[1], pointInfos, triangleInfos))))
////////		;
////////}
////
////bool ShipFloorplanizer::IsSpringViableForFloorAlsoInternal(
////	ShipFactorySpring const & springInfo,
////	std::vector<ShipFactoryPoint> const & pointInfos) const
////{
////	return
////		pointInfos[springInfo.PointAIndex].DefinitionCoordinates.has_value() // Is point derived directly from structure?
////		&& pointInfos[springInfo.PointAIndex].Material.IsHull // Is point hull?
////		&& pointInfos[springInfo.PointBIndex].DefinitionCoordinates.has_value() // Is point derived directly from structure?
////		&& pointInfos[springInfo.PointBIndex].Material.IsHull // Is point hull?
////		&& springInfo.Triangles.size() > 0 // Is it an edge of a triangle?
////		;
////}
////
////bool ShipFloorplanizer::IsSpringEndpointViableForFloor(
////	ElementIndex pointIndex,
////	ShipFactorySpring const & spring,
////	std::vector<ShipFactoryPoint> const & pointInfos) const
////{
////	return
////		pointInfos[pointIndex].DefinitionCoordinates.has_value() // Is point derived directly from structure?
////		&& pointInfos[pointIndex].Material.IsHull // Is point hull?
////		&& spring.Triangles.size() > 0 // Is it an edge of a triangle?
////		;
////}
////
////// TODOOLD
////
////ShipFactoryFloorPlan ShipFloorplanizer::BuildFloorplan_Old(
////	ShipFactoryPointIndexMatrix const & pointIndexMatrix,
////	std::vector<ShipFactoryPoint> const & pointInfos,
////	std::vector<ShipFactorySpring> const & springInfos,
////	std::vector<ShipFactoryTriangle> & triangleInfos) const
////{
////	//
////	// Naive implementation: follows spring hullness
////	//
////	// Note: the floorplan will not necessarily match triangle edges
////	// (as we might have suppressed triangles in a quad); this floorplan
////	// will anyway be a *superset* of the floor edges
////	//
////
////	//
////	// 1. Build first floorplan using all and only hull springs
////	//
////
////	ShipFactoryFloorPlan hullSprings;
////	hullSprings.reserve(springInfos.size());
////
////	for (size_t s = 0; s < springInfos.size(); ++s)
////	{
////		auto const & springInfo = springInfos[s];
////
////		// Make sure it derives from structure
////		if (pointInfos[springInfo.PointAIndex].DefinitionCoordinates.has_value()
////			&& pointInfos[springInfo.PointBIndex].DefinitionCoordinates.has_value())
////		{
////			// Make sure it's a "hull spring" (i.e. that both endpoints are hull)
////			if (pointInfos[springInfo.PointAIndex].Material.IsHull
////				&& pointInfos[springInfo.PointBIndex].Material.IsHull)
////			{
////				NpcFloorType floorType;
////				if (pointInfos[springInfo.PointAIndex].DefinitionCoordinates->x == pointInfos[springInfo.PointBIndex].DefinitionCoordinates->x)
////				{
////					// Horizontal
////					assert(std::abs(pointInfos[springInfo.PointAIndex].DefinitionCoordinates->y - pointInfos[springInfo.PointBIndex].DefinitionCoordinates->y) == 1);
////					floorType = NpcFloorType::FloorPlane1;
////				}
////				else if (pointInfos[springInfo.PointAIndex].DefinitionCoordinates->y == pointInfos[springInfo.PointBIndex].DefinitionCoordinates->y)
////				{
////					// Vertical
////					assert(std::abs(pointInfos[springInfo.PointAIndex].DefinitionCoordinates->x - pointInfos[springInfo.PointBIndex].DefinitionCoordinates->x) == 1);
////					floorType = NpcFloorType::FloorPlane1;
////				}
////				else
////				{
////					// Diagonal
////					assert(std::abs(pointInfos[springInfo.PointAIndex].DefinitionCoordinates->x - pointInfos[springInfo.PointBIndex].DefinitionCoordinates->x) == 1);
////					assert(std::abs(pointInfos[springInfo.PointAIndex].DefinitionCoordinates->y - pointInfos[springInfo.PointBIndex].DefinitionCoordinates->y) == 1);
////					floorType = NpcFloorType::FloorPlane2;
////				}
////
////				auto const [_, isInserted] = hullSprings.try_emplace(
////					{ springInfo.PointAIndex, springInfo.PointBIndex },
////					floorType,
////					static_cast<ElementIndex>(s));
////
////				assert(isInserted);
////				(void)isInserted;
////			}
////		}
////	}
////
////	////// TODOTEST
////	////(void)pointIndexMatrix;
////	////(void)triangleInfos;
////	////return hullSprings;
////
////	//
////	// 2. Do a first pass removing redundant springs
////	//
////
////	auto intermediateFloorPlan = RemoveRedundantFloors(
////		std::move(hullSprings),
////		pointIndexMatrix,
////		pointInfos,
////		springInfos,
////		triangleInfos);
////
////	//
////	// 3. Do a last pass removing redundant springs
////	//
////
////	// TODOTEST
////	return intermediateFloorPlan;
////
////	////auto floorPlan = RemoveRedundantFloors(
////	////	std::move(intermediateFloorPlan),
////	////	pointIndexMatrix,
////	////	pointInfos,
////	////	springInfos,
////	////	triangleInfos);
////
////	////return floorPlan;
////}
////
//////////////////////////////////////////////////////////////////////////////////////////////////////////////
////
////ShipFactoryFloorPlan ShipFloorplanizer::RemoveRedundantFloors(
////	ShipFactoryFloorPlan && hullSprings,
////	ShipFactoryPointIndexMatrix const & pointIndexMatrix,
////	std::vector<ShipFactoryPoint> const & pointInfos,
////	std::vector<ShipFactorySpring> const & springInfos,
////	std::vector<ShipFactoryTriangle> & triangleInfos) const
////{
////	ShipFactoryFloorPlan floorPlan;
////	floorPlan.reserve(springInfos.size());
////
////	//
////	// Visit all "hull springs"
////	//
////
////	for (auto const & hullSpringEntry : hullSprings)
////	{
////		auto const & springInfo = springInfos[hullSpringEntry.second.SpringIndex];
////
////		// TODOTEST
////		if (hullSpringEntry.second.SpringIndex == 22)
////			LogMessage("TODOTEST");
////
////		// Make sure it's not a redundant spring
////		if (!IsRedundantFloorSpring(
////			hullSpringEntry.second.SpringIndex,
////			hullSprings,
////			pointIndexMatrix,
////			pointInfos,
////			springInfos,
////			triangleInfos))
////		{
////			NpcFloorType floorType;
////			if (pointInfos[springInfo.PointAIndex].DefinitionCoordinates->x == pointInfos[springInfo.PointBIndex].DefinitionCoordinates->x)
////			{
////				// Horizontal
////				assert(std::abs(pointInfos[springInfo.PointAIndex].DefinitionCoordinates->y - pointInfos[springInfo.PointBIndex].DefinitionCoordinates->y) == 1);
////				floorType = NpcFloorType::FloorPlane1;
////			}
////			else if (pointInfos[springInfo.PointAIndex].DefinitionCoordinates->y == pointInfos[springInfo.PointBIndex].DefinitionCoordinates->y)
////			{
////				// Vertical
////				assert(std::abs(pointInfos[springInfo.PointAIndex].DefinitionCoordinates->x - pointInfos[springInfo.PointBIndex].DefinitionCoordinates->x) == 1);
////				floorType = NpcFloorType::FloorPlane1;
////			}
////			else
////			{
////				// Diagonal
////				assert(std::abs(pointInfos[springInfo.PointAIndex].DefinitionCoordinates->x - pointInfos[springInfo.PointBIndex].DefinitionCoordinates->x) == 1);
////				assert(std::abs(pointInfos[springInfo.PointAIndex].DefinitionCoordinates->y - pointInfos[springInfo.PointBIndex].DefinitionCoordinates->y) == 1);
////				floorType = NpcFloorType::FloorPlane2;
////			}
////
////			auto const [_, isInserted] = floorPlan.try_emplace(
////				hullSpringEntry.first,
////				hullSpringEntry.second);
////
////			assert(isInserted);
////			(void)isInserted;
////		}
////	}
////
////	return floorPlan;
////}
////
////bool ShipFloorplanizer::IsRedundantFloorSpring(
////	ElementIndex s,
////	ShipFactoryFloorPlan const & hullSprings,
////	ShipFactoryPointIndexMatrix const & pointIndexMatrix,
////	std::vector<ShipFactoryPoint> const & pointInfos,
////	std::vector<ShipFactorySpring> const & springInfos,
////	std::vector<ShipFactoryTriangle> & triangleInfos) const
////{
////	//
////	// S is a redundant floor spring (i.e. a spring that comes out of tessellation
////	// but not really meant to be floor) if and only if:
////	//	- It's a hull spring
////	//	- It's connected to two triangles
////	//	- One of the two triangles (T):
////	//		- Has all hull sides
////	//		- From both of the two vertices of S we encounter (in CW or CCW order, depending on vertex) hull springs at "steep" slopes (delta octant neither 4 nor 5)
////	//	- The other triangle (T'):
////	//		- Has at least one non-hull side (at least one is open) (obviously not S)
////	//
////	//                  Psb`    /Psb'
////	//                   |    /
////	//                   |  /
////	//                   |/
////	//           Pc-----Psb
////	//           | T'   /|
////	//           |   S/  |
////	//           |  /    |
////	//	         |/   T  |
////	// Psa'-----Psa------Ps
////	//
////
////	auto const & springInfo = springInfos[s];
////
////	// It's a hull spring
////	assert(hullSprings.count({ springInfo.PointAIndex, springInfo.PointBIndex }) == 1);
////
////	// It's connected to two triangles
////	if (springInfo.Triangles.size() != 2)
////	{
////		return false;
////	}
////
////	// The two triangles meet the criteria
////	if ((
////			IsRedundantFloorSpringTriangle1(s, springInfo, triangleInfos[springInfo.Triangles[0]], hullSprings, pointIndexMatrix, pointInfos)
////			&& IsRedundantFloorSpringTriangle2(triangleInfos[springInfo.Triangles[1]], hullSprings)
////		)
////		||
////		(
////			IsRedundantFloorSpringTriangle1(s, springInfo, triangleInfos[springInfo.Triangles[1]], hullSprings, pointIndexMatrix, pointInfos)
////			&& IsRedundantFloorSpringTriangle2(triangleInfos[springInfo.Triangles[0]], hullSprings)
////		))
////	{
////		return true;
////	}
////
////	return false;
////}
////
////bool ShipFloorplanizer::IsRedundantFloorSpringTriangle1(
////	ElementIndex s,
////	ShipFactorySpring const & spring,
////	ShipFactoryTriangle const & triangle,
////	ShipFactoryFloorPlan const & hullSprings,
////	ShipFactoryPointIndexMatrix const & pointIndexMatrix,
////	std::vector<ShipFactoryPoint> const & pointInfos) const
////{
////	//	- Has all hull sides
////
////	size_t hullSideCount = CountTriangleHullSides(triangle, hullSprings);
////	assert(hullSideCount >= 1 && hullSideCount <= 3);
////
////	// Has all hull sides
////	if (hullSideCount != 3)
////	{
////		return false;
////	}
////
////	// - From both of the two vertices of S we encounter (in CW or CCW order, depending on vertex) hull springs at "steep" slopes(delta octant neither 4 nor 5)
////
////	int const sOrdinal = triangle.FindSpringOrdinal(s);
////
////	if (!HasTriangleSteepHullExtensions(
////		triangle.PointIndices[sOrdinal], // Psa
////		spring,
////		7, // CCW
////		hullSprings,
////		pointIndexMatrix,
////		pointInfos))
////	{
////		// One smooth path
////		return false;
////	}
////
////	if (!HasTriangleSteepHullExtensions(
////		triangle.PointIndices[(sOrdinal + 1) % 3], // Psb
////		spring,
////		1, // CW
////		hullSprings,
////		pointIndexMatrix,
////		pointInfos))
////	{
////		// One smooth path
////		return false;
////	}
////
////	// No smooth paths
////	return true;
////}
////
////bool ShipFloorplanizer::IsRedundantFloorSpringTriangle2(
////	ShipFactoryTriangle const & triangle,
////	ShipFactoryFloorPlan const & hullSprings) const
////{
////	//	- Has no hull sides (other than S)
////
////	size_t hullSideCount = CountTriangleHullSides(triangle, hullSprings);
////	assert(hullSideCount >= 1 && hullSideCount <= 3);
////	return hullSideCount != 3;
////}
////
////bool ShipFloorplanizer::HasTriangleSteepHullExtensions(
////	ElementIndex startPoint,
////	ShipFactorySpring const & spring,
////	Octant direction,
////	ShipFactoryFloorPlan const & hullSprings,
////	ShipFactoryPointIndexMatrix const & pointIndexMatrix,
////	std::vector<ShipFactoryPoint> const & pointInfos) const
////{
////	Octant startOctant;
////	if (startPoint == spring.PointAIndex)
////	{
////		startOctant = spring.PointAAngle;
////	}
////	else
////	{
////		assert(startPoint == spring.PointBIndex);
////		startOctant = spring.PointBAngle;
////	}
////
////	Octant const deltaOctant = FindNextHullSpringOctant(
////		startPoint,
////		startOctant,
////		direction,
////		hullSprings,
////		pointIndexMatrix,
////		pointInfos);
////
////	return deltaOctant != 4 && deltaOctant != 5;
////}
////
////Octant ShipFloorplanizer::FindNextHullSpringOctant(
////	ElementIndex startPoint,
////	Octant startOctant,
////	Octant direction,
////	ShipFactoryFloorPlan const & hullSprings,
////	ShipFactoryPointIndexMatrix const & pointIndexMatrix,
////	std::vector<ShipFactoryPoint> const & pointInfos) const
////{
////	assert(pointInfos[startPoint].DefinitionCoordinates.has_value());
////	auto const startCoords = *(pointInfos[startPoint].DefinitionCoordinates);
////	Octant distanceOctant = 1;
////	for (Octant octant = (startOctant + direction) % 8; ; octant = (octant + direction) % 8, ++distanceOctant)
////	{
////		auto const destCoords = vec2i(
////			startCoords.x + 1 + TessellationCircularOrderDirections[octant][0],
////			startCoords.y + 1 + TessellationCircularOrderDirections[octant][1]);
////		if (pointIndexMatrix[destCoords].has_value()
////			&& hullSprings.count({ startPoint, *pointIndexMatrix[destCoords] }) > 0)
////		{
////			// Found hull spring
////			break;
////		}
////	}
////
////	return distanceOctant;
////}
////
////size_t ShipFloorplanizer::CountTriangleHullSides(
////	ShipFactoryTriangle const & triangle,
////	ShipFactoryFloorPlan const & hullSprings) const
////{
////	size_t const count =
////		hullSprings.count({ triangle.PointIndices[0], triangle.PointIndices[1] })
////		+ hullSprings.count({ triangle.PointIndices[1], triangle.PointIndices[2] })
////		+ hullSprings.count({ triangle.PointIndices[2], triangle.PointIndices[0] });
////
////	return count;
////}
