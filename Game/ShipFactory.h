/***************************************************************************************
* Original Author:      Gabriele Giuseppini
* Created:              2020-05-16
* Copyright:            Gabriele Giuseppini  (https://github.com/GabrieleGiuseppini)
***************************************************************************************/
#pragma once

#include "MaterialDatabase.h"
#include "Physics.h"
#include "ShipDefinition.h"
#include "ShipFactoryTypes.h"

#include <GameCore/GameTypes.h>

#include <memory>

/*
 * This class contains all the logic for building a ship out of a Definition.
 */
class ShipFactory
{
public:

    static std::unique_ptr<Physics::Ship> BuildShip(
        ShipDefinition && shipDefinition,
        MaterialDatabase const & materialDatabase);

private:

    static void CreateElementInfos(
        ShipFactoryPointIndexMatrix const & pointIndexMatrix,
        std::vector<ShipFactoryPoint> & pointInfos,
        std::vector<ShipFactorySpring> & springInfos,
        std::vector<ShipFactoryTriangle> & triangleInfos);

    static void ConnectPointsToTriangles(
        std::vector<ShipFactoryPoint> & pointInfos,
        std::vector<ShipFactoryTriangle> const & triangleInfos);

    static void ConnectSpringsToTriangles(
        std::vector<ShipFactorySpring> & springInfos,
        std::vector<ShipFactoryTriangle> & triangleInfos);

    static Physics::Points CreatePoints(
        std::vector<ShipFactoryPoint> const & pointInfos);

    static Physics::Springs CreateSprings(
        std::vector<ShipFactorySpring> const & springInfos,
        Physics::Points & points);

    static Physics::Triangles CreateTriangles(
        std::vector<ShipFactoryTriangle> const & triangleInfos,
        Physics::Points & points,
        std::vector<ShipFactorySpring> const & springInfos,
        std::vector<ShipFactoryFloor> const & floorInfos);
};
