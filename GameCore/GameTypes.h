/***************************************************************************************
* Original Author:		Gabriele Giuseppini
* Created:				2020-05-15
* Copyright:			Gabriele Giuseppini  (https://github.com/GabrieleGiuseppini)
***************************************************************************************/
#pragma once

#include "BarycentricCoords.h"
#include "Colors.h"
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
 * Our local circular order (clockwise, starting from E), indexes by Octant.
 * Note: cardinal directions are labeled according to x growing to the right and y growing upwards
 */
extern int const TessellationCircularOrderDirections[8][2];

/*
 * Generic directions.
 */
enum class DirectionType
{
    Horizontal = 1,
    Vertical = 2
};

//template <> struct is_flag<DirectionType> : std::true_type {};

/*
 * Generic rotation directions.
 */
enum class RotationDirectionType
{
    Clockwise,
    CounterClockwise
};

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

using IntegralRectSize = _IntegralSize<struct IntegralTag>;
using ImageSize = _IntegralSize<struct ImageTag>;
using ShipSpaceSize = _IntegralSize<struct ShipSpaceTag>;

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

    _IntegralCoordinates<TIntegralTag> FlipX(integral_type width) const
    {
        assert(width > x);
        return _IntegralCoordinates<TIntegralTag>(width - 1 - x, y);
    }

    _IntegralCoordinates<TIntegralTag> FlipY(integral_type height) const
    {
        assert(height > y);
        return _IntegralCoordinates<TIntegralTag>(x, height - 1 - y);
    }

    // Returns coords of this point after being rotated (and assuming
    // the size will also get rotated).
    template<RotationDirectionType TDirection>
    _IntegralCoordinates<TIntegralTag> Rotate90(_IntegralSize<TIntegralTag> const & sz) const
    {
        if constexpr (TDirection == RotationDirectionType::Clockwise)
        {
            return _IntegralCoordinates<TIntegralTag>(
                this->y,
                sz.width - 1 - this->x);
        }
        else
        {
            static_assert(TDirection == RotationDirectionType::CounterClockwise);

            return _IntegralCoordinates<TIntegralTag>(
                sz.height - 1 - this->y,
                this->x);
        }
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

#pragma pack(push)

template<typename TIntegralTag>
struct _IntegralRect
{
    _IntegralCoordinates<TIntegralTag> origin;
    _IntegralSize<TIntegralTag> size;

    constexpr _IntegralRect()
        : origin(0, 0)
        , size(0, 0)
    {}

    constexpr _IntegralRect(
        _IntegralCoordinates<TIntegralTag> const & _origin,
        _IntegralSize<TIntegralTag> const & _size)
        : origin(_origin)
        , size(_size)
    {}

    explicit constexpr _IntegralRect(_IntegralCoordinates<TIntegralTag> const & _origin)
        : origin(_origin)
        , size(1, 1)
    {}

    constexpr _IntegralRect(
        _IntegralCoordinates<TIntegralTag> const & _origin,
        _IntegralCoordinates<TIntegralTag> const & _oppositeCorner)
        : origin(
            std::min(_origin.x, _oppositeCorner.x),
            std::min(_origin.y, _oppositeCorner.y))
        , size(
            std::abs(_oppositeCorner.x - _origin.x),
            std::abs(_oppositeCorner.y - _origin.y))
    {
    }

    /*
     * Makes a rectangle from {0, 0} of the specified size.
     */
    explicit constexpr _IntegralRect(_IntegralSize<TIntegralTag> const & _size)
        : origin(0, 0)
        , size(_size)
    {}

    _IntegralCoordinates<TIntegralTag> MinMin() const
    {
        return origin;
    }

    _IntegralCoordinates<TIntegralTag> MaxMin() const
    {
        return _IntegralCoordinates<TIntegralTag>(
            origin.x + size.width,
            origin.y);
    }

    _IntegralCoordinates<TIntegralTag> MaxMax() const
    {
        return _IntegralCoordinates<TIntegralTag>(
            origin.x + size.width,
            origin.y + size.height);
    }

    _IntegralCoordinates<TIntegralTag> MinMax() const
    {
        return _IntegralCoordinates<TIntegralTag>(
            origin.x,
            origin.y + size.height);
    }

    _IntegralCoordinates<TIntegralTag> Center() const
    {
        return _IntegralCoordinates<TIntegralTag>(
            origin.x + size.width / 2,
            origin.y + size.height / 2);
    }

    inline bool operator==(_IntegralRect<TIntegralTag> const & other) const
    {
        return origin == other.origin
            && size == other.size;
    }

    inline bool operator!=(_IntegralRect<TIntegralTag> const & other) const
    {
        return !(*this == other);
    }

    bool IsEmpty() const
    {
        return size.width == 0 || size.height == 0;
    }

    bool IsContainedInRect(_IntegralRect<TIntegralTag> const & container) const
    {
        return origin.x >= container.origin.x
            && origin.y >= container.origin.y
            && origin.x + size.width <= container.origin.x + container.size.width
            && origin.y + size.height <= container.origin.y + container.size.height;
    }

    void UnionWith(_IntegralCoordinates<TIntegralTag> const & other)
    {
        auto const newOrigin = _IntegralCoordinates<TIntegralTag>(
            std::min(origin.x, other.x),
            std::min(origin.y, other.y));

        auto const newSize = _IntegralSize<TIntegralTag>(
            std::max(origin.x + size.width, other.x + 1) - newOrigin.x,
            std::max(origin.y + size.height, other.y + 1) - newOrigin.y);

        assert(newSize.width >= 0 && newSize.height >= 0);

        origin = newOrigin;
        size = newSize;
    }

    void UnionWith(_IntegralRect<TIntegralTag> const & other)
    {
        auto const newOrigin = _IntegralCoordinates<TIntegralTag>(
            std::min(origin.x, other.origin.x),
            std::min(origin.y, other.origin.y));

        auto const newSize = _IntegralSize<TIntegralTag>(
            std::max(origin.x + size.width, other.origin.x + other.size.width) - newOrigin.x,
            std::max(origin.y + size.height, other.origin.y + other.size.height) - newOrigin.y);

        assert(newSize.width >= 0 && newSize.height >= 0);

        origin = newOrigin;
        size = newSize;
    }

    std::optional<_IntegralRect<TIntegralTag>> MakeIntersectionWith(_IntegralRect<TIntegralTag> const & other) const
    {
        auto const newOrigin = _IntegralCoordinates<TIntegralTag>(
            std::max(origin.x, other.origin.x),
            std::max(origin.y, other.origin.y));

        auto const newSize = _IntegralSize<TIntegralTag>(
            std::min(size.width - (newOrigin.x - origin.x), other.size.width - (newOrigin.x - other.origin.x)),
            std::min(size.height - (newOrigin.y - origin.y), other.size.height - (newOrigin.y - other.origin.y)));

        if (newSize.width <= 0 || newSize.height <= 0)
        {
            return std::nullopt;
        }
        else
        {
            return _IntegralRect<TIntegralTag>(
                newOrigin,
                newSize);
        }
    }

    std::string ToString() const
    {
        std::stringstream ss;
        ss << "(" << origin.x << ", " << origin.y << " -> " << size.width << " x " << size.height << ")";
        return ss.str();
    }
};

#pragma pack(pop)

template<typename TTag>
inline std::basic_ostream<char> & operator<<(std::basic_ostream<char> & os, _IntegralRect<TTag> const & p)
{
    os << p.origin.ToString() << "x" << p.size.ToString();
    return os;
}


using IntegralRect = _IntegralRect<struct IntegralTag>;
using ImageRect = _IntegralRect<struct ImageTag>;
using ShipSpaceRect = _IntegralRect<struct ShipSpaceTag>;


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

/*
 * Barycentric coordinates in a specific triangle.
 */
struct AbsoluteTriangleBCoords
{
    ElementIndex TriangleElementIndex;
    bcoords3f BCoords;

    AbsoluteTriangleBCoords() = default;

    AbsoluteTriangleBCoords(
        ElementIndex triangleElementIndex,
        bcoords3f bCoords)
        : TriangleElementIndex(triangleElementIndex)
        , BCoords(bCoords)
    {
        assert(triangleElementIndex != NoneElementIndex);
    }

    bool operator==(AbsoluteTriangleBCoords const & other) const
    {
        return this->TriangleElementIndex == other.TriangleElementIndex
            && this->BCoords == other.BCoords;
    }

    std::string ToString() const
    {
        std::stringstream ss;
        ss << TriangleElementIndex << ":" << BCoords;
        return ss.str();
    }
};

inline std::basic_ostream<char> & operator<<(std::basic_ostream<char> & os, AbsoluteTriangleBCoords const & is)
{
    os << is.ToString();
    return os;
}

////////////////////////////////////////////////////////////////////////////////////////////////
// Rendering
////////////////////////////////////////////////////////////////////////////////////////////////

////////////////////////////////////////////////////////////////////////////////////////////////
// Misc
////////////////////////////////////////////////////////////////////////////////////////////////

/*
 * The color key of materials.
 */
using MaterialColorKey = rgbColor;

enum class SimulationControlStateType
{
    Paused = 0,
    Play
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
 * Second level of NPC type hierarchy; domain is open
 * as it may be expanded after compile time, via NPC packs.
 * The unique identifier of an NPC kind is the whole
 * <NpcKindType,NpcSubKindIdType> tuple; so, for example,
 * NpcSubKindIdType=X means one thing for Humans and another
 * thing for Furniture.
 */
using NpcSubKindIdType = std::uint32_t;

/*
 * Roles for humans.
 */
enum class NpcHumanRoleType
{
    Passenger
};

NpcHumanRoleType StrToNpcHumanRoleType(std::string const & str);

enum class NpcFloorKindType
{
    NotAFloor,
    DefaultFloor // Futurework: areas, etc.
};

enum class NpcFloorGeometryDepthType
{
    NotAFloor,
    Depth1, // Main depth: H-V
    Depth2 // Staircases: S-S
};

enum class NpcFloorGeometryType
{
    NotAFloor,
    // Depth 1: main depth
    Depth1H,
    Depth1V,
    // Depth 2: staircases
    Depth2S1,
    Depth2S2
};

inline NpcFloorGeometryDepthType NpcFloorGeometryDepth(NpcFloorGeometryType geometry)
{
    switch (geometry)
    {
        case NpcFloorGeometryType::NotAFloor:
        {
            return NpcFloorGeometryDepthType::NotAFloor;
        }

        case NpcFloorGeometryType::Depth1H:
        case NpcFloorGeometryType::Depth1V:
        {
            return NpcFloorGeometryDepthType::Depth1;
        }

        case NpcFloorGeometryType::Depth2S1:
        case NpcFloorGeometryType::Depth2S2:
        {
            return NpcFloorGeometryDepthType::Depth2;
        }
    }

    assert(false);
    return NpcFloorGeometryDepthType::NotAFloor;
}

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

