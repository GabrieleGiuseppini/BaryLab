/***************************************************************************************
 * Original Author:     Gabriele Giuseppini
 * Created:             2024-05-07
 * Copyright:           Gabriele Giuseppini  (https://github.com/GabrieleGiuseppini)
 ***************************************************************************************/
#pragma once

#include "ShipFactoryTypes.h"

#include <GameCore/GameTypes.h>

#include <array>
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

	using VertexBlock = std::array<std::array<ElementIndex, 3>, 3>;
	using SpringExclusionSet = std::unordered_set<ShipFactoryPointPair, ShipFactoryPointPair::Hasher>;

	void ProcessVertexBlock(
		VertexBlock & vertexBlock,
		/*out*/ SpringExclusionSet & springExclusionSet) const;

	void ProcessVertexBlockPatterns(
		VertexBlock const & vertexBlock,
		/*out*/ SpringExclusionSet & springExclusionSet) const;

	void Rotate90CW(VertexBlock & vertexBlock) const;
	void FlipV(VertexBlock & vertexBlock) const;
};
