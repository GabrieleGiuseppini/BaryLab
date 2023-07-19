/***************************************************************************************
* Original Author:      Gabriele Giuseppini
* Created:              2020-05-16
* Copyright:            Gabriele Giuseppini  (https://github.com/GabrieleGiuseppini)
***************************************************************************************/
#pragma once

#include "AABB.h"
#include "BLabTypes.h"
#include "Buffer.h"
#include "ElementContainer.h"
#include "ElementIndexRangeIterator.h"
#include "FixedSizeVector.h"
#include "LabParameters.h"
#include "Vectors.h"

#include <algorithm>
#include <cassert>
#include <optional>

class Points : public ElementContainer
{
public:

    /*
     * The metadata of a single spring connected to a point.
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

    using ConnectedSpringsVector = FixedSizeVector<ConnectedSpring, LabParameters::MaxSpringsPerPoint>;

    /*
     * The metadata of all the triangles connected to a point.
     */
    struct ConnectedTrianglesVector
    {
        FixedSizeVector<ElementIndex, LabParameters::MaxTrianglesPerPoint> ConnectedTriangles;
        size_t OwnedConnectedTrianglesCount;

        ConnectedTrianglesVector()
            : ConnectedTriangles()
            , OwnedConnectedTrianglesCount(0u)
        {}

        inline void ConnectTriangle(
            ElementIndex triangleElementIndex,
            bool isAtOwner)
        {
            // Add so that all triangles owned by this point come first
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
        // Physics
        , mPositionBuffer(mBufferElementCount, pointCount, vec2f::zero())
        , mVelocityBuffer(mBufferElementCount, pointCount, vec2f::zero())
        , mAssignedForceBuffer(mBufferElementCount, pointCount, vec2f::zero())
        // Structure        
        , mConnectedSpringsBuffer(mBufferElementCount, pointCount, ConnectedSpringsVector())
        , mConnectedTrianglesBuffer(mBufferElementCount, pointCount, ConnectedTrianglesVector())
        // Render
        , mRenderColorBuffer(mBufferElementCount, pointCount, vec4f::zero())
    {
    }

    Points(Points && other) = default;

    void Add(
        vec2f const & position,
        vec3f const & color);

    AABB GetAABB() const
    {
        AABB box;
        for (ElementIndex pointIndex : *this)
        {
            box.ExtendTo(GetPosition(pointIndex));
        }

        return box;
    }

    void Query(ElementIndex pointElementIndex) const;

public:

    //
    // Physics
    //

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

    vec2f const & GetVelocity(ElementIndex pointElementIndex) const noexcept
    {
        return mVelocityBuffer[pointElementIndex];
    }

    vec2f const * GetVelocityBuffer() const noexcept
    {
        return mVelocityBuffer.data();
    }

    vec2f * GetVelocityBuffer() noexcept
    {
        return mVelocityBuffer.data();
    }

    void SetVelocity(
        ElementIndex pointElementIndex,
        vec2f const & value) noexcept
    {
        mVelocityBuffer[pointElementIndex] = value;
    }

    vec2f const & GetAssignedForce(ElementIndex pointElementIndex) const noexcept
    {
        return mAssignedForceBuffer[pointElementIndex];
    }

    vec2f * GetAssignedForceBuffer() noexcept
    {
        return mAssignedForceBuffer.data();
    }

    void SetAssignedForce(
        ElementIndex pointElementIndex,
        vec2f const & value) noexcept
    {
        mAssignedForceBuffer[pointElementIndex] = value;
    }

    //
    // Structure
    //

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

    void ConnectTriangle(
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

    //
    // Render
    //

    vec4f const & GetRenderColor(ElementIndex pointElementIndex) const
    {
        return mRenderColorBuffer[pointElementIndex];
    }

    void SetRenderColor(
        ElementIndex pointElementIndex,
        vec4f const & color)
    {
        mRenderColorBuffer[pointElementIndex] = color;
    }

    vec4f const * GetRenderColorBuffer() const
    {
        return mRenderColorBuffer.data();
    }

private:

    //////////////////////////////////////////////////////////
    // Buffers
    //////////////////////////////////////////////////////////

    //
    // Physics
    //

    Buffer<vec2f> mPositionBuffer;
    Buffer<vec2f> mVelocityBuffer;
    Buffer<vec2f> mAssignedForceBuffer;
    
    //
    // Structure
    //

    Buffer<ConnectedSpringsVector> mConnectedSpringsBuffer;
    Buffer<ConnectedTrianglesVector> mConnectedTrianglesBuffer;

    //
    // Render
    //

    Buffer<vec4f> mRenderColorBuffer;
};
