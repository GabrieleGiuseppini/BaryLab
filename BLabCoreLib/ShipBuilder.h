/***************************************************************************************
* Original Author:      Gabriele Giuseppini
* Created:              2020-05-16
* Copyright:            Gabriele Giuseppini  (https://github.com/GabrieleGiuseppini)
***************************************************************************************/
#pragma once

#include "GameTypes.h"
#include "MaterialDatabase.h"
#include "Physics.h"
#include "ShipBuilderTypes.h"
#include "ShipDefinition.h"

#include <memory>

/*
 * This class contains all the logic for building a ship out of a Definition.
 */
class ShipBuilder
{
public:

    static std::unique_ptr<Physics::Ship> BuildShip(
        ShipDefinition && shipDefinition,
        MaterialDatabase const & materialDatabase);

private:

    static void CreateElementInfos(
        ShipBuildPointIndexMatrix const & pointIndexMatrix,
        std::vector<ShipBuildPoint> & pointInfos,
        std::vector<ShipBuildSpring> & springInfos,
        std::vector<ShipBuildTriangle> & triangleInfos);

    static void ConnectPointsToTriangles(
        std::vector<ShipBuildPoint> & pointInfos,
        std::vector<ShipBuildTriangle> const & triangleInfos);

    static void ConnectSpringsToTriangles(
        std::vector<ShipBuildSpring> & springInfos,
        std::vector<ShipBuildTriangle> & triangleInfos);

    static Physics::Points CreatePoints(
        std::vector<ShipBuildPoint> const & pointInfos);

    static Physics::Springs CreateSprings(
        std::vector<ShipBuildSpring> const & springInfos,
        std::vector<ShipBuildPoint> & pointInfos,
        Physics::Points & points);

    static Physics::Triangles CreateTriangles(
        std::vector<ShipBuildTriangle> const & triangleInfos,
        Physics::Points & points,
        std::vector<ShipBuildSpring> const & springInfos,
        Physics::Springs const & springs);
};
