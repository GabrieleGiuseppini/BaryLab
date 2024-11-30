/***************************************************************************************
* Original Author:      Gabriele Giuseppini
* Created:              2023-07-20
* Copyright:            Gabriele Giuseppini  (https://github.com/GabrieleGiuseppini)
***************************************************************************************/
#pragma once

#include "GameParameters.h"
#include "Materials.h"
#include "Physics.h"

#include <GameCore/AABB.h>
#include <GameCore/Buffer.h>
#include <GameCore/ElementContainer.h>
#include <GameCore/ElementIndexRangeIterator.h>
#include <GameCore/FixedSizeVector.h>
#include <GameCore/GameTypes.h>
#include <GameCore/Log.h>
#include <GameCore/Vectors.h>

#include <algorithm>
#include <cassert>
#include <optional>

namespace Physics {

class Points : public ElementContainer
{
public:

    /*
     * The metadata of a single edge connected to a vertex.
     */
    struct ConnectedSpring
    {
        ElementIndex SpringIndex;
        ElementIndex OtherEndpointIndex;

        ConnectedSpring()
            : SpringIndex(NoneElementIndex)
            , OtherEndpointIndex(NoneElementIndex)
        {}

        ConnectedSpring(
            ElementIndex springIndex,
            ElementIndex otherEndpointIndex)
            : SpringIndex(springIndex)
            , OtherEndpointIndex(otherEndpointIndex)
        {}
    };

    using ConnectedSpringsVector = FixedSizeVector<ConnectedSpring, GameParameters::MaxSpringsPerPoint>;

    /*
     * The metadata of all the triangles connected to a vertex.
     */
    struct ConnectedTrianglesVector
    {
        FixedSizeVector<ElementIndex, GameParameters::MaxTrianglesPerPoint> ConnectedTriangles;
        size_t OwnedConnectedTrianglesCount;

        ConnectedTrianglesVector()
            : ConnectedTriangles()
            , OwnedConnectedTrianglesCount(0u)
        {}

        inline void ConnectTriangle(
            ElementIndex triangleElementIndex,
            bool isAtOwner)
        {
            // Add so that all triangles owned by this vertex come first
            if (isAtOwner)
            {
                ConnectedTriangles.emplace_front(triangleElementIndex);
                ++OwnedConnectedTrianglesCount;
            }
            else
            {
                ConnectedTriangles.emplace_back(triangleElementIndex);
            }
        }

        inline void DisconnectTriangle(
            ElementIndex triangleElementIndex,
            bool isAtOwner)
        {
            bool found = ConnectedTriangles.erase_first(
                [triangleElementIndex](ElementIndex c)
                {
                    return c == triangleElementIndex;
                });

            assert(found);
            (void)found;

            // Update count of owned triangles, if this triangle is owned
            if (isAtOwner)
            {
                assert(OwnedConnectedTrianglesCount > 0);
                --OwnedConnectedTrianglesCount;
            }
        }
    };

public:

    Points(ElementCount pointCount)
        : ElementContainer(pointCount)
        //////////////////////////////////
        // Buffers
        //////////////////////////////////
        , mPositionBuffer(mBufferElementCount, pointCount, vec2f::zero())
        , mMaterialBuffer(mBufferElementCount, pointCount, nullptr)
        , mConnectedSpringsBuffer(mBufferElementCount, pointCount, ConnectedSpringsVector())
        , mConnectedTrianglesBuffer(mBufferElementCount, pointCount, ConnectedTrianglesVector())
    {
    }

    Points(Points && other) = default;

    void Add(
        vec2f const & position,
        StructuralMaterial const * material);

    Geometry::AABB GetAABB() const
    {
        Geometry::AABB box;
        for (ElementIndex pointIndex : *this)
        {
            box.ExtendTo(GetPosition(pointIndex));
        }

        return box;
    }

    bool QueryAt(vec2f const & worldCoordinates) const;

    void Query(ElementIndex pointElementIndex) const;

public:

    vec2f const & GetPosition(ElementIndex pointElementIndex) const noexcept
    {
        return mPositionBuffer[pointElementIndex];
    }

    vec2f const * GetPositionBuffer() const noexcept
    {
        return mPositionBuffer.data();
    }

    vec2f * GetPositionBuffer() noexcept
    {
        return mPositionBuffer.data();
    }

    void SetPosition(
        ElementIndex pointElementIndex,
        vec2f const & value) noexcept
    {
        mPositionBuffer[pointElementIndex] = value;
    }

    vec2f GetVelocity(ElementIndex pointElementIndex) const noexcept
    {
        (void)pointElementIndex;
        return vec2f::zero();
    }

    void AddStaticForce(
        ElementIndex /*pointElementIndex*/,
        vec2f const & force) noexcept
    {
        LogDebug("Static force: ", force, " (", force / GameParameters::GravityMagnitude, " kg)");
        (void)force;
    }

    StructuralMaterial const & GetStructuralMaterial(ElementIndex pointElementIndex) const noexcept
    {
        return *(mMaterialBuffer[pointElementIndex]);
    }

    PlaneId GetPlaneId(ElementIndex /*pointElementIndex*/) const
    {
        // Placeholder
        return 0;
    }

    ConnectedComponentId GetConnectedComponentId(ElementIndex /*pointElementIndex*/) const
    {
        // Placeholder
        return 0;
    }

    auto const & GetConnectedSprings(ElementIndex pointElementIndex) const
    {
        return mConnectedSpringsBuffer[pointElementIndex];
    }

    void AddConnectedSpring(
        ElementIndex pointElementIndex,
        ElementIndex springElementIndex,
        ElementIndex otherEndpointElementIndex)
    {
        assert(!mConnectedSpringsBuffer[pointElementIndex].contains(
            [springElementIndex](auto const & cs)
            {
                return cs.SpringIndex == springElementIndex;
            }));

        mConnectedSpringsBuffer[pointElementIndex].emplace_back(
            springElementIndex,
            otherEndpointElementIndex);
    }

    auto const & GetConnectedTriangles(ElementIndex pointElementIndex) const
    {
        return mConnectedTrianglesBuffer[pointElementIndex];
    }

    void AddConnectedTriangle(
        ElementIndex pointElementIndex,
        ElementIndex triangleElementIndex,
        bool isAtOwner)
    {
        assert(!mConnectedTrianglesBuffer[pointElementIndex].ConnectedTriangles.contains(
            [triangleElementIndex](auto const & ct)
            {
                return ct == triangleElementIndex;
            }));

        mConnectedTrianglesBuffer[pointElementIndex].ConnectTriangle(
            triangleElementIndex,
            isAtOwner);
    }

    // Placeholders

    float GetWater(ElementIndex pointElementIndex) const;

    vec2f GetWaterVelocity(ElementIndex pointElementIndex) const;

    void AddTransientAdditionalMass(
        ElementIndex pointElementIndex,
        float mass);

    float GetTemperature(ElementIndex pointElementIndex) const;

    float GetLight(ElementIndex pointElementIndex) const;

    bool GetIsHull(ElementIndex pointElementIndex) const;

    bool IsBurning(ElementIndex pointElementIndex) const;

    bool GetIsElectrified(ElementIndex pointElementIndex) const;

private:

    Buffer<vec2f> mPositionBuffer;
    Buffer<StructuralMaterial const *> mMaterialBuffer;
    Buffer<ConnectedSpringsVector> mConnectedSpringsBuffer;
    Buffer<ConnectedTrianglesVector> mConnectedTrianglesBuffer;
};

}