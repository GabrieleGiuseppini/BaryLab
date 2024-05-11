/***************************************************************************************
 * Original Author:     Gabriele Giuseppini
 * Created:             2024-05-07
 * Copyright:           Gabriele Giuseppini  (https://github.com/GabrieleGiuseppini)
 ***************************************************************************************/
#pragma once

#include "ShipFactoryTypes.h"

#include <GameCore/GameTypes.h>

class ShipFloorplanizer
{
public:

	ShipFloorplanizer() = default;

	ShipFactoryFloorPlan BuildFloorplan(
		ShipFactoryPointIndexMatrix const & pointIndexMatrix,
		std::vector<ShipFactoryPoint> const & pointInfos,
		std::vector<ShipFactorySpring> const & springInfos,
		std::vector<ShipFactoryTriangle> & triangleInfos) const;

private:

	bool IsRedundantFloorSpring(
		ElementIndex s,
		ShipFactoryPointIndexMatrix const & pointIndexMatrix,
		std::vector<ShipFactoryPoint> const & pointInfos,
		std::vector<ShipFactorySpring> const & springInfos,
		std::vector<ShipFactoryTriangle> & triangleInfos) const;

	bool IsRedundantFloorSpringTriangle1(
		ElementIndex s,
		ShipFactorySpring const & spring,
		ShipFactoryTriangle const & triangle,
		ShipFactoryPointIndexMatrix const & pointIndexMatrix,
		std::vector<ShipFactoryPoint> const & pointInfos) const;

	bool IsRedundantFloorSpringTriangle2(
		ShipFactoryTriangle const & triangle,
		std::vector<ShipFactoryPoint> const & pointInfos) const;

	bool HasTriangleSteepHullExtensions(
		ElementIndex startPoint,
		ShipFactorySpring const & spring,
		Octant direction,
		ShipFactoryPointIndexMatrix const & pointIndexMatrix,
		std::vector<ShipFactoryPoint> const & pointInfos) const;

	Octant FindNextHullSpringOctant(
		ElementIndex startPoint,
		Octant startOctant,
		Octant direction,
		ShipFactoryPointIndexMatrix const & pointIndexMatrix,
		std::vector<ShipFactoryPoint> const & pointInfos) const;

	size_t CountTriangleHullSides(
		ShipFactoryTriangle const & triangle,
		std::vector<ShipFactoryPoint> const & pointInfos) const;

	bool AreBothEndpointsHull(
		ShipFactorySpring const & spring,
		std::vector<ShipFactoryPoint> const & pointInfos) const;
};
