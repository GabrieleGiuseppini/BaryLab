/***************************************************************************************
* Original Author:      Gabriele Giuseppini
* Created:              2023-07-20
* Copyright:            Gabriele Giuseppini  (https://github.com/GabrieleGiuseppini)
***************************************************************************************/
#pragma once

#include "AABB.h"
#include "Buffer.h"
#include "ElementContainer.h"
#include "ElementIndexRangeIterator.h"
#include "FixedSizeVector.h"
#include "GameTypes.h"
#include "LabParameters.h"
#include "Physics.h"
#include "Vectors.h"

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

    using ConnectedSpringsVector = FixedSizeVector<ConnectedSpring, LabParameters::MaxSpringsPerPoint>;

    /*
     * The metadata of all the triangles connected to a vertex.
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
        , mConnectedSpringsBuffer(mBufferElementCount, pointCount, ConnectedSpringsVector())
        , mConnectedTrianglesBuffer(mBufferElementCount, pointCount, ConnectedTrianglesVector())
    {
    }

    Points(Points && other) = default;

    void Add(vec2f const & position);

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

private:

    Buffer<vec2f> mPositionBuffer;
    Buffer<ConnectedSpringsVector> mConnectedSpringsBuffer;
    Buffer<ConnectedTrianglesVector> mConnectedTrianglesBuffer;
};

}