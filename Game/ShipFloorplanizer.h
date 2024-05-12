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

	ShipFactoryFloorPlan BuildFloorplan_Old(
		ShipFactoryPointIndexMatrix const & pointIndexMatrix,
		std::vector<ShipFactoryPoint> const & pointInfos,
		std::vector<ShipFactorySpring> const & springInfos,
		std::vector<ShipFactoryTriangle> & triangleInfos) const;

private:

	// TODO:needed?
	bool IsSealedTriangle(
		ElementIndex t,
		std::vector<ShipFactoryPoint> const & pointInfos,
		std::vector<ShipFactoryTriangle> & triangleInfos) const;

	bool IsIsolatedFloor(
		ShipFactorySpring const & spring,
		ShipFactoryPointIndexMatrix const & pointIndexMatrix,
		std::vector<ShipFactoryPoint> const & pointInfos) const;

	bool DoesFloorContinue(
		ElementIndex pointIndex,
		Octant direction,
		ShipFactoryPointIndexMatrix const & pointIndexMatrix,
		std::vector<ShipFactoryPoint> const & pointInfos) const;

	// TODOOLD

	ShipFactoryFloorPlan RemoveRedundantFloors(
		ShipFactoryFloorPlan && hullSprings,
		ShipFactoryPointIndexMatrix const & pointIndexMatrix,
		std::vector<ShipFactoryPoint> const & pointInfos,
		std::vector<ShipFactorySpring> const & springInfos,
		std::vector<ShipFactoryTriangle> & triangleInfos) const;

	bool IsRedundantFloorSpring(
		ElementIndex s,
		ShipFactoryFloorPlan const & hullSprings,
		ShipFactoryPointIndexMatrix const & pointIndexMatrix,
		std::vector<ShipFactoryPoint> const & pointInfos,
		std::vector<ShipFactorySpring> const & springInfos,
		std::vector<ShipFactoryTriangle> & triangleInfos) const;

	bool IsRedundantFloorSpringTriangle1(
		ElementIndex s,
		ShipFactorySpring const & spring,
		ShipFactoryTriangle const & triangle,
		ShipFactoryFloorPlan const & hullSprings,
		ShipFactoryPointIndexMatrix const & pointIndexMatrix,
		std::vector<ShipFactoryPoint> const & pointInfos) const;

	bool IsRedundantFloorSpringTriangle2(
		ShipFactoryTriangle const & triangle,
		ShipFactoryFloorPlan const & hullSprings) const;

	bool HasTriangleSteepHullExtensions(
		ElementIndex startPoint,
		ShipFactorySpring const & spring,
		Octant direction,
		ShipFactoryFloorPlan const & hullSprings,
		ShipFactoryPointIndexMatrix const & pointIndexMatrix,
		std::vector<ShipFactoryPoint> const & pointInfos) const;

	Octant FindNextHullSpringOctant(
		ElementIndex startPoint,
		Octant startOctant,
		Octant direction,
		ShipFactoryFloorPlan const & hullSprings,
		ShipFactoryPointIndexMatrix const & pointIndexMatrix,
		std::vector<ShipFactoryPoint> const & pointInfos) const;

	size_t CountTriangleHullSides(
		ShipFactoryTriangle const & triangle,
		ShipFactoryFloorPlan const & hullSprings) const;
};
