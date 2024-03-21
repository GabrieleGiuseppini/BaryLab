/***************************************************************************************
 * Original Author:		Gabriele Giuseppini
 * Created:				2023-07-19
 * Copyright:			Gabriele Giuseppini  (https://github.com/GabrieleGiuseppini)
 ***************************************************************************************/
#pragma once

#include "Physics.h"

 // Placeholder for real Ship

namespace Physics {

class Ship
{
public:

    Ship(
        ShipId id,
        Points && points,
        Springs && springs,
        Triangles && triangles)
        : mId(id)
        , mPoints(std::move(points))
        , mSprings(std::move(springs))
        , mTriangles(std::move(triangles))
    {}

    ShipId GetId() const
    {
        return mId;
    }

    PlaneId GetMaxPlaneId() const
    {
        // Placeholder
        return 0;
    }

    Points const & GetPoints() const
    {
        return mPoints;
    }

    Points & GetPoints()
    {
        return mPoints;
    }

    Springs const & GetSprings() const
    {
        return mSprings;
    }

    Springs & GetSprings()
    {
        return mSprings;
    }

    Triangles const & GetTriangles() const
    {
        return mTriangles;
    }

    Triangles & GetTriangles()
    {
        return mTriangles;
    }

private:

    ShipId const mId;
    Points mPoints;
    Springs mSprings;
    Triangles mTriangles;
};

}