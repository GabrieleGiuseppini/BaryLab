/***************************************************************************************
 * Original Author:     Gabriele Giuseppini
 * Created:             2024-05-07
 * Copyright:           Gabriele Giuseppini  (https://github.com/GabrieleGiuseppini)
 ***************************************************************************************/
#include "ShipFloorplanizer.h"

std::vector<ShipFactoryFloor> ShipFloorplanizer::BuildFloorplan(
	std::vector<ShipFactoryPoint> & pointInfos,
	std::vector<ShipFactorySpring> & springInfos)
{
	//
	// Naive implementation: follows spring hullness
	//

	std::vector<ShipFactoryFloor> floors;
	floors.reserve(springInfos.size());

	for (auto const & springInfo : springInfos)
	{
		bool const isPointAHull = pointInfos[springInfo.PointAIndex].Material.IsHull;
		bool const isPointBHull = pointInfos[springInfo.PointBIndex].Material.IsHull;
		if (isPointAHull && isPointBHull)
		{
			NpcFloorType floorType = NpcFloorType::Open;
			if (pointInfos[springInfo.PointAIndex].DefinitionCoordinates.has_value()
				&& pointInfos[springInfo.PointBIndex].DefinitionCoordinates.has_value())
			{
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
			}

			floors.emplace_back(
				springInfo.PointAIndex,
				springInfo.PointBIndex,
				floorType);
		}
	}

	return floors;
}