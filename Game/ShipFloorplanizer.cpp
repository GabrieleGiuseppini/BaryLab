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
	std::vector<ShipFactorySpring> const & springInfos) const
{
	//
	// 1. Build list of springs that we do not want to use as floors;
	//    we do so by detecting specific vertex patterns in 3x3 blocks
	//

	SpringExclusionSet springExclusionSet;

	// Process all 3x3 blocks - including the 1-wide "borders"
	VertexBlock vertexBlock;
	for (int y = 0; y < pointIndexMatrix.height - 2; ++y)
	{
		for (int x = 0; x < pointIndexMatrix.width - 2; ++x)
		{
			// Build block
			for (int yb = 0; yb < 3; ++yb)
			{
				for (int xb = 0; xb < 3; ++xb)
				{
					if (pointIndexMatrix[{x + xb, y + yb}].has_value()
						&& pointInfos[*pointIndexMatrix[{x + xb, y + yb}]].Material.IsHull)
					{
						vertexBlock[xb][yb] = *pointIndexMatrix[{x + xb, y + yb}];
					}
					else
					{
						vertexBlock[xb][yb] = NoneElementIndex;
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
				// Vertical
				assert(std::abs(pointInfos[springInfo.PointAIndex].DefinitionCoordinates->y - pointInfos[springInfo.PointBIndex].DefinitionCoordinates->y) == 1);
				floorType = NpcFloorType::FloorPlane1V;
			}
			else if (pointInfos[springInfo.PointAIndex].DefinitionCoordinates->y == pointInfos[springInfo.PointBIndex].DefinitionCoordinates->y)
			{
				// Horizontal
				assert(std::abs(pointInfos[springInfo.PointAIndex].DefinitionCoordinates->x - pointInfos[springInfo.PointBIndex].DefinitionCoordinates->x) == 1);
				floorType = NpcFloorType::FloorPlane1H;
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

		Rotate90CW(vertexBlock);

	}

	// 2. All rotations of symmetry 2

	FlipV(vertexBlock);

	for (int i = 0; ; ++i)
	{
		ProcessVertexBlockPatterns(
			vertexBlock,
			springExclusionSet);

		// Save one rotation
		if (i == 3)
		{
			break;
		}

		Rotate90CW(vertexBlock);
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

	if (vertexBlock[0][0] != NoneElementIndex && vertexBlock[1][0] != NoneElementIndex && vertexBlock[2][0] != NoneElementIndex
		&& vertexBlock[0][1] == NoneElementIndex && vertexBlock[1][1] != NoneElementIndex
		&& vertexBlock[0][2] != NoneElementIndex && vertexBlock[2][2] == NoneElementIndex)
	{
		springExclusionSet.insert({ vertexBlock[0][0] , vertexBlock[1][1] });
	}

	//
	// Pattern 2: "under a stair" (_\): take care of redundant |
	//
	//   *oo
	//   o*?
	//   ***
	//

	if (vertexBlock[0][0] != NoneElementIndex && vertexBlock[1][0] != NoneElementIndex && vertexBlock[2][0] != NoneElementIndex
		&& vertexBlock[0][1] == NoneElementIndex && vertexBlock[1][1] != NoneElementIndex
		&& vertexBlock[0][2] != NoneElementIndex && vertexBlock[1][2] == NoneElementIndex && vertexBlock[2][2] == NoneElementIndex)
	{
		springExclusionSet.insert({ vertexBlock[1][0] , vertexBlock[1][1] });
	}

	//
	// Pattern 3: "wall-on-floor" (|-): take care of redundant \ and /
	//
	//  o?o
	//  o*o
	//  ***
	//

	if (vertexBlock[0][0] != NoneElementIndex && vertexBlock[1][0] != NoneElementIndex && vertexBlock[2][0] != NoneElementIndex
		&& vertexBlock[0][1] == NoneElementIndex && vertexBlock[1][1] != NoneElementIndex && vertexBlock[2][1] == NoneElementIndex
		&& vertexBlock[0][2] == NoneElementIndex && vertexBlock[2][2] == NoneElementIndex)
	{
		springExclusionSet.insert({ vertexBlock[0][0] , vertexBlock[1][1] });
		springExclusionSet.insert({ vertexBlock[2][0] , vertexBlock[1][1] });
	}

	//
	// Pattern 4: "corner" (|_): take care of redundant \
	//
	//  *o?
	//  *o?
	//  ***
	//

	if (vertexBlock[0][0] != NoneElementIndex && vertexBlock[1][0] != NoneElementIndex && vertexBlock[2][0] != NoneElementIndex
		&& vertexBlock[0][1] != NoneElementIndex && vertexBlock[1][1] == NoneElementIndex
		&& vertexBlock[0][2] != NoneElementIndex && vertexBlock[1][2] == NoneElementIndex)
	{
		springExclusionSet.insert({ vertexBlock[0][1] , vertexBlock[1][0] });
	}

	//
	// Pattern 5: "floor-meets-high-wall" (_||): take care of redundant /
	//
	//  o**
	//  o**
	//  ***
	//

	if (vertexBlock[0][0] != NoneElementIndex && vertexBlock[1][0] != NoneElementIndex && vertexBlock[2][0] != NoneElementIndex
		&& vertexBlock[0][1] == NoneElementIndex && vertexBlock[1][1] != NoneElementIndex && vertexBlock[2][1] != NoneElementIndex
		&& vertexBlock[0][2] == NoneElementIndex && vertexBlock[1][2] != NoneElementIndex && vertexBlock[2][2] != NoneElementIndex)
	{
		springExclusionSet.insert({ vertexBlock[0][0] , vertexBlock[1][1] });
	}
}

void ShipFloorplanizer::Rotate90CW(VertexBlock & vertexBlock) const
{
	auto const tmp1 = vertexBlock[0][0];
	vertexBlock[0][0] = vertexBlock[2][0];
	vertexBlock[2][0] = vertexBlock[2][2];
	vertexBlock[2][2] = vertexBlock[0][2];
	vertexBlock[0][2] = tmp1;

	auto const tmp2 = vertexBlock[1][0];
	vertexBlock[1][0] = vertexBlock[2][1];
	vertexBlock[2][1] = vertexBlock[1][2];
	vertexBlock[1][2] = vertexBlock[0][1];
	vertexBlock[0][1] = tmp2;
}

void ShipFloorplanizer::FlipV(VertexBlock & vertexBlock) const
{
	std::swap(vertexBlock[0][0], vertexBlock[0][2]);
	std::swap(vertexBlock[1][0], vertexBlock[1][2]);
	std::swap(vertexBlock[2][0], vertexBlock[2][2]);
}
