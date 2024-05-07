/***************************************************************************************
 * Original Author:     Gabriele Giuseppini
 * Created:             2024-05-07
 * Copyright:           Gabriele Giuseppini  (https://github.com/GabrieleGiuseppini)
 ***************************************************************************************/
#pragma once

#include "ShipFactoryTypes.h"

#include <GameCore/GameTypes.h>

#include <vector>

class ShipFloorplanizer
{
public:

	ShipFloorplanizer() = default;

	std::vector<ShipFactoryFloor> BuildFloorplan(
		std::vector<ShipFactoryPoint> & pointInfos,
		std::vector<ShipFactorySpring> & springInfos);

private:

};
