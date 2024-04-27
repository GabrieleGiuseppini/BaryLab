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

/*
 * Ship identifiers.
 *
 * Comparable and ordered. Start from 0.
 */
using ShipId = std::uint32_t;
static ShipId constexpr NoneShip = std::numeric_limits<ShipId>::max();

/*
 * Plane (depth) identifiers.
 *
 * Comparable and ordered. Start from 0.
 */
using PlaneId = std::uint32_t;
static PlaneId constexpr NonePlaneId = std::numeric_limits<PlaneId>::max();

/*
 * NPC (global) identifiers.
 *
 * Comparable and ordered. Start from 0.
 */
using NpcId = std::uint32_t;
static NpcId constexpr NoneNpcId = std::numeric_limits<NpcId>::max();

/*
 * Object ID's, identifying objects of ships across ships.
 *
 * An ObjectId is unique only in the context in which it's used; for example,
 * a gadget might have the same object ID as a switch. That's where the type tag
 * comes from.
 *
 * Not comparable, not ordered.
 */
template<typename TLocalObjectId, typename TTypeTag>
struct ObjectId
{
    using LocalObjectId = TLocalObjectId;

    ObjectId(
        ShipId shipId,
        LocalObjectId localObjectId)
        : mShipId(shipId)
        , mLocalObjectId(localObjectId)
    {}

    inline ShipId GetShipId() const noexcept
    {
        return mShipId;
    };

    inline LocalObjectId GetLocalObjectId() const noexcept
    {
        return mLocalObjectId;
    }

    ObjectId & operator=(ObjectId const & other) = default;

    inline bool operator==(ObjectId const & other) const
    {
        return this->mShipId == other.mShipId
            && this->mLocalObjectId == other.mLocalObjectId;
    }

    inline bool operator<(ObjectId const & other) const
    {
        return this->mShipId < other.mShipId
            || (this->mShipId == other.mShipId && this->mLocalObjectId < other.mLocalObjectId);
    }

    std::string ToString() const
    {
        std::stringstream ss;

        ss << static_cast<int>(mShipId) << ":" << static_cast<int>(mLocalObjectId);

        return ss.str();
    }

private:

    ShipId mShipId;
    LocalObjectId mLocalObjectId;
};

// Generic ID for generic elements (points, springs, etc.)
using ElementId = ObjectId<ElementIndex, struct _ElementIdTag>;

/*
 * Return type of picking an object.
 */
template<typename TObjectId>
struct PickedObjectId
{
    TObjectId ObjectId;
    vec2f WorldOffset;

    PickedObjectId(
        TObjectId objectId,
        vec2f const & worldOffset)
        : ObjectId(objectId)
        , WorldOffset(worldOffset)
    {}
};

/*
 * Identifies the edge of a triangle among all edges on a ship.
 */
struct TriangleAndEdge
{
    ElementIndex TriangleElementIndex;
    int EdgeOrdinal;

    TriangleAndEdge(
        ElementIndex triangleElementIndex,
        int edgeOrdinal)
        : TriangleElementIndex(triangleElementIndex)
        , EdgeOrdinal(edgeOrdinal)
    {
        assert(triangleElementIndex != NoneElementIndex);
        assert(edgeOrdinal >= 0 && edgeOrdinal < 3);
    }
};

////////////////////////////////////////////////////////////////////////////////////////////////
// Geometry
////////////////////////////////////////////////////////////////////////////////////////////////

/*
 * Octants, i.e. the direction of a spring connecting two neighbors.
 *
 * Octant 0 is E, octant 1 is SE, ..., Octant 7 is NE.
 */
using Octant = std::int32_t;

////////////////////////////////////////////////////////////////////////////////////////////////
// Rendering
////////////////////////////////////////////////////////////////////////////////////////////////

struct TextureQuad
{
    vec2f TopLeftPosition;
    vec2f TopLeftTexture;
    vec2f TopRightPosition;
    vec2f TopRightTexture;
    vec2f BottomLeftPosition;
    vec2f BottomLeftTexture;
    vec2f BottomRightPosition;
    vec2f BottomRightTexture;

    TextureQuad() = default;

    TextureQuad(
        vec2f topLeftPosition,
        vec2f topLeftTexture,
        vec2f topRightPosition,
        vec2f topRightTexture,
        vec2f bottomLeftPosition,
        vec2f bottomLeftTexture,
        vec2f bottomRightPosition,
        vec2f bottomRightTexture)
        : TopLeftPosition(topLeftPosition)
        , TopLeftTexture(topLeftTexture)
        , TopRightPosition(topRightPosition)
        , TopRightTexture(topRightTexture)
        , BottomLeftPosition(bottomLeftPosition)
        , BottomLeftTexture(bottomLeftTexture)
        , BottomRightPosition(bottomRightPosition)
        , BottomRightTexture(bottomRightTexture)
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

////////////////////////////////////////////////////////////////////////////////////////////////
// NPC
////////////////////////////////////////////////////////////////////////////////////////////////

/*
 * Top level of NPC type hierarchy.
 */
enum class NpcKindType
{
    Furniture,
    Human
};

/*
 * Furniture NPC types (second level).
 */
enum class FurnitureNpcKindType
{
    Particle
};

/*
 * Human NPC types (second level).
 */
enum class HumanNpcKindType
{
    Passenger,
    Programmer
};

// Futurework: FurnitureNpcKind

enum class NpcSurfaceType
{
    Floor,
    Open
};

NpcSurfaceType StrToNpcSurfaceType(std::string const & str);

/*
 * Types of hightlight for NPCs.
 */

enum class NpcHighlightType
{
    None,
    Candidate
};

enum class NpcRenderModeType
{
    Physical, // Only in Barylab
    Limbs
};

