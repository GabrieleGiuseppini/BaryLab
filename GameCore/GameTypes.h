/***************************************************************************************
* Original Author:		Gabriele Giuseppini
* Created:				2020-05-15
* Copyright:			Gabriele Giuseppini  (https://github.com/GabrieleGiuseppini)
***************************************************************************************/
#pragma once

#include "BarycentricCoords.h"
#include "Vectors.h"

#include <picojson.h>

#include <cassert>
#include <cstdint>
#include <limits>
#include <sstream>
#include <string>

////////////////////////////////////////////////////////////////////////////////////////////////
// Data Structures
////////////////////////////////////////////////////////////////////////////////////////////////

/*
 * These types define the cardinality of elements in the ElementContainer.
 *
 * Indices are equivalent to pointers in OO terms. Given that we don't believe
 * we'll ever have more than 4 billion elements, a 32-bit integer suffices.
 *
 * This also implies that where we used to store one pointer, we can now store two indices,
 * resulting in even better data locality.
 */
using ElementCount = std::uint32_t;
using ElementIndex = std::uint32_t;
static constexpr ElementIndex NoneElementIndex = std::numeric_limits<ElementIndex>::max();

////////////////////////////////////////////////////////////////////////////////////////////////
// Geometry
////////////////////////////////////////////////////////////////////////////////////////////////

/*
 * Octants, i.e. the direction of a spring connecting two neighbors.
 *
 * Octant 0 is E, octant 1 is SE, ..., Octant 7 is NE.
 */
using Octant = std::int32_t;

struct Quadf
{
    vec2f TopRight;
    vec2f TopLeft;
    vec2f BottomRight;
    vec2f BottomLeft;

    Quadf() = default;

    Quadf(
        vec2f topRight,
        vec2f topLeft,
        vec2f bottomRight,
        vec2f bottomLeft)
        : TopRight(topRight)
        , TopLeft(topLeft)
        , BottomRight(bottomRight)
        , BottomLeft(bottomLeft)
    {}
};

////////////////////////////////////////////////////////////////////////////////////////////////
// Misc
////////////////////////////////////////////////////////////////////////////////////////////////

enum class SimulationControlStateType
{
    Paused = 0,
    Play
};

struct ConstrainedRegimeParticleProbe
{
    ElementIndex CurrentTriangle;
    bcoords3f CurrentTriangleBarycentricCoords;

    ConstrainedRegimeParticleProbe(
        ElementIndex currentTriangle,
        bcoords3f const & currentTriangleBarycentricCoords)
        : CurrentTriangle(currentTriangle)
        , CurrentTriangleBarycentricCoords(currentTriangleBarycentricCoords)
    {}
};

struct PhysicsParticleProbe
{
    vec2f Velocity;

    PhysicsParticleProbe(
        vec2f const & velocity)
        : Velocity(velocity)
    {}
};

struct ParticleTrajectory
{
    ElementIndex ParticleIndex;
    vec2f TargetPosition;

    ParticleTrajectory(
        ElementIndex particleIndex,
        vec2f const & targetPosition)
        : ParticleIndex(particleIndex)
        , TargetPosition(targetPosition)
    {}
};

enum class SurfaceType
{
    Floor,
    Open
};

SurfaceType StrToSurfaceType(std::string const & str);

enum class NpcRenderMode
{
    Physical,
    Limbs
};
