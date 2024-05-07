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

    TriangleAndEdge() = default;

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

/*
 * Integral system
 */

#pragma pack(push, 1)

template<typename TIntegralTag>
struct _IntegralSize
{
    using integral_type = int;

    integral_type width;
    integral_type height;

    constexpr _IntegralSize(
        integral_type _width,
        integral_type _height)
        : width(_width)
        , height(_height)
    {}

    static _IntegralSize<TIntegralTag> FromFloatRound(vec2f const & vec)
    {
        return _IntegralSize<TIntegralTag>(
            static_cast<integral_type>(std::round(vec.x)),
            static_cast<integral_type>(std::round(vec.y)));
    }

    static _IntegralSize<TIntegralTag> FromFloatFloor(vec2f const & vec)
    {
        return _IntegralSize<TIntegralTag>(
            static_cast<integral_type>(std::floor(vec.x)),
            static_cast<integral_type>(std::floor(vec.y)));
    }

    inline bool operator==(_IntegralSize<TIntegralTag> const & other) const
    {
        return this->width == other.width
            && this->height == other.height;
    }

    inline bool operator!=(_IntegralSize<TIntegralTag> const & other) const
    {
        return !(*this == other);
    }

    inline _IntegralSize<TIntegralTag> operator+(_IntegralSize<TIntegralTag> const & sz) const
    {
        return _IntegralSize<TIntegralTag>(
            this->width + sz.width,
            this->height + sz.height);
    }

    inline void operator+=(_IntegralSize<TIntegralTag> const & sz)
    {
        this->width += sz.width;
        this->height += sz.height;
    }

    inline _IntegralSize<TIntegralTag> operator*(integral_type factor) const
    {
        return _IntegralSize<TIntegralTag>(
            this->width * factor,
            this->height * factor);
    }

    inline size_t GetLinearSize() const
    {
        return this->width * this->height;
    }

    inline void Rotate90()
    {
        std::swap(width, height);
    }

    inline _IntegralSize<TIntegralTag> Union(_IntegralSize<TIntegralTag> const & other) const
    {
        return _IntegralSize<TIntegralTag>(
            std::max(this->width, other.width),
            std::max(this->height, other.height));
    }

    inline _IntegralSize<TIntegralTag> Intersection(_IntegralSize<TIntegralTag> const & other) const
    {
        return _IntegralSize<TIntegralTag>(
            std::min(this->width, other.width),
            std::min(this->height, other.height));
    }

    vec2f ToFloat() const
    {
        return vec2f(
            static_cast<float>(width),
            static_cast<float>(height));
    }

    template<typename TCoordsRatio>
    vec2f ToFractionalCoords(TCoordsRatio const & coordsRatio) const
    {
        assert(coordsRatio.inputUnits != 0.0f);

        return vec2f(
            static_cast<float>(width) / coordsRatio.inputUnits * coordsRatio.outputUnits,
            static_cast<float>(height) / coordsRatio.inputUnits * coordsRatio.outputUnits);
    }

    std::string ToString() const
    {
        std::stringstream ss;
        ss << "(" << width << " x " << height << ")";
        return ss.str();
    }
};

#pragma pack(pop)

template<typename TTag>
inline std::basic_ostream<char> & operator<<(std::basic_ostream<char> & os, _IntegralSize<TTag> const & is)
{
    os << is.ToString();
    return os;
}

using ImageSize = _IntegralSize<struct ImageTag>;

#pragma pack(push, 1)

template<typename TIntegralTag>
struct _IntegralCoordinates
{
    using integral_type = int;

    integral_type x;
    integral_type y;

    constexpr _IntegralCoordinates(
        integral_type _x,
        integral_type _y)
        : x(_x)
        , y(_y)
    {}

    static _IntegralCoordinates<TIntegralTag> FromFloatRound(vec2f const & vec)
    {
        return _IntegralCoordinates<TIntegralTag>(
            static_cast<integral_type>(std::round(vec.x)),
            static_cast<integral_type>(std::round(vec.y)));
    }

    static _IntegralCoordinates<TIntegralTag> FromFloatFloor(vec2f const & vec)
    {
        return _IntegralCoordinates<TIntegralTag>(
            static_cast<integral_type>(std::floor(vec.x)),
            static_cast<integral_type>(std::floor(vec.y)));
    }

    inline bool operator==(_IntegralCoordinates<TIntegralTag> const & other) const
    {
        return this->x == other.x
            && this->y == other.y;
    }

    inline bool operator!=(_IntegralCoordinates<TIntegralTag> const & other) const
    {
        return !(*this == other);
    }

    inline _IntegralCoordinates<TIntegralTag> operator+(_IntegralSize<TIntegralTag> const & sz) const
    {
        return _IntegralCoordinates<TIntegralTag>(
            this->x + sz.width,
            this->y + sz.height);
    }

    inline void operator+=(_IntegralSize<TIntegralTag> const & sz)
    {
        this->x += sz.width;
        this->y += sz.height;
    }

    inline _IntegralCoordinates<TIntegralTag> operator-() const
    {
        return _IntegralCoordinates<TIntegralTag>(
            -this->x,
            -this->y);
    }

    inline _IntegralSize<TIntegralTag> operator-(_IntegralCoordinates<TIntegralTag> const & other) const
    {
        return _IntegralSize<TIntegralTag>(
            this->x - other.x,
            this->y - other.y);
    }

    inline _IntegralCoordinates<TIntegralTag> operator-(_IntegralSize<TIntegralTag> const & offset) const
    {
        return _IntegralCoordinates<TIntegralTag>(
            this->x - offset.width,
            this->y - offset.height);
    }

    inline _IntegralCoordinates<TIntegralTag> scale(_IntegralCoordinates<TIntegralTag> const & multiplier) const
    {
        return _IntegralCoordinates<TIntegralTag>(
            this->x * multiplier.x,
            this->y * multiplier.y);
    }

    template<typename TSize>
    bool IsInSize(TSize const & size) const
    {
        return x >= 0 && x < size.width && y >= 0 && y < size.height;
    }

    template<typename TRect>
    bool IsInRect(TRect const & rect) const
    {
        return x >= rect.origin.x && x < rect.origin.x + rect.size.width
            && y >= rect.origin.y && y < rect.origin.y + rect.size.height;
    }

    vec2f ToFloat() const
    {
        return vec2f(
            static_cast<float>(x),
            static_cast<float>(y));
    }

    template<typename TCoordsRatio>
    vec2f ToFractionalCoords(TCoordsRatio const & coordsRatio) const
    {
        assert(coordsRatio.inputUnits != 0.0f);

        return vec2f(
            static_cast<float>(x) / coordsRatio.inputUnits * coordsRatio.outputUnits,
            static_cast<float>(y) / coordsRatio.inputUnits * coordsRatio.outputUnits);
    }

    std::string ToString() const
    {
        std::stringstream ss;
        ss << "(" << x << ", " << y << ")";
        return ss.str();
    }
};

#pragma pack(pop)

template<typename TTag>
inline bool operator<(_IntegralCoordinates<TTag> const & l, _IntegralCoordinates<TTag> const & r)
{
    return l.x < r.x
        || (l.x == r.x && l.y < r.y);
}

template<typename TTag>
inline std::basic_ostream<char> & operator<<(std::basic_ostream<char> & os, _IntegralCoordinates<TTag> const & p)
{
    os << p.ToString();
    return os;
}

using IntegralCoordinates = _IntegralCoordinates<struct IntegralTag>; // Generic integer
using ShipSpaceCoordinates = _IntegralCoordinates<struct ShipSpaceTag>; // Y=0 at bottom

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

enum class NpcFloorType
{
    Open,
    FloorPlane1, // Planes N: a category of Horiz|Vert or Diag
    FloorPlane2
};

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

