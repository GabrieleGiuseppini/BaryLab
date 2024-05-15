/***************************************************************************************
 * Original Author:     Gabriele Giuseppini
 * Created:             2024-05-07
 * Copyright:           Gabriele Giuseppini  (https://github.com/GabrieleGiuseppini)
 ***************************************************************************************/
#pragma once

#include "ShipFactoryTypes.h"

#include <GameCore/Buffer2D.h>
#include <GameCore/GameTypes.h>

#include <unordered_set>

class ShipFloorplanizer
{
public:

	ShipFloorplanizer() = default;

	ShipFactoryFloorPlan BuildFloorplan(
		ShipFactoryPointIndexMatrix const & pointIndexMatrix,
		std::vector<ShipFactoryPoint> const & pointInfos,
		std::vector<ShipFactorySpring> const & springInfos) const;

private:

	bool IsSpringViableForFloor(
		ShipFactorySpring const & springInfo,
		std::vector<ShipFactoryPoint> const & pointInfos) const;

	using VertexBlock = Buffer2D<ElementIndex, struct IntegralTag>;
	using SpringExclusionSet = std::unordered_set<ShipFactoryPointPair, ShipFactoryPointPair::Hasher>;

	void ProcessVertexBlock(
		VertexBlock & vertexBlock,
		/*out*/ SpringExclusionSet & springExclusionSet) const;

	void ProcessVertexBlockPatterns(
		VertexBlock const & vertexBlock,
		/*out*/ SpringExclusionSet & springExclusionSet) const;
};
