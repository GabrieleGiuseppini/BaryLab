/***************************************************************************************
* Original Author:      Gabriele Giuseppini
* Created:              2023-07-20
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
#include "Physics.h"
#include "Vectors.h"

#include <algorithm>
#include <cassert>
#include <optional>

namespace Physics {

class Vertices : public ElementContainer
{
public:

    /*
     * The metadata of a single edge connected to a vertex.
     */
    struct ConnectedEdge
    {
        ElementIndex EdgeIndex;
        ElementIndex OtherEndpointIndex;

        ConnectedEdge()
            : EdgeIndex(NoneElementIndex)
            , OtherEndpointIndex(NoneElementIndex)
        {}

        ConnectedEdge(
            ElementIndex edgeIndex,
            ElementIndex otherEndpointIndex)
            : EdgeIndex(edgeIndex)
            , OtherEndpointIndex(otherEndpointIndex)
        {}
    };

    using ConnectedEdgesVector = FixedSizeVector<ConnectedEdge, LabParameters::MaxEdgesPerVertex>;

    /*
     * The metadata of all the triangles connected to a vertex.
     */
    struct ConnectedTrianglesVector
    {
        FixedSizeVector<ElementIndex, LabParameters::MaxTrianglesPerVertex> ConnectedTriangles;
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

    Vertices(ElementCount vertexCount)
        : ElementContainer(vertexCount)
        //////////////////////////////////
        // Buffers
        //////////////////////////////////
        , mPositionBuffer(mBufferElementCount, vertexCount, vec2f::zero())
        , mConnectedEdgesBuffer(mBufferElementCount, vertexCount, ConnectedEdgesVector())
        , mConnectedTrianglesBuffer(mBufferElementCount, vertexCount, ConnectedTrianglesVector())
    {
    }

    Vertices(Vertices && other) = default;

    void Add(vec2f const & position);

    AABB GetAABB() const
    {
        AABB box;
        for (ElementIndex vertexIndex : *this)
        {
            box.ExtendTo(GetPosition(vertexIndex));
        }

        return box;
    }

    void Query(ElementIndex vertexElementIndex) const;

public:

    vec2f const & GetPosition(ElementIndex vertexElementIndex) const noexcept
    {
        return mPositionBuffer[vertexElementIndex];
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
        ElementIndex vertexElementIndex,
        vec2f const & value) noexcept
    {
        mPositionBuffer[vertexElementIndex] = value;
    }

    auto const & GetConnectedEdges(ElementIndex vertexElementIndex) const
    {
        return mConnectedEdgesBuffer[vertexElementIndex];
    }

    void AddConnectedEdge(
        ElementIndex vertexElementIndex,
        ElementIndex edgeElementIndex,
        ElementIndex otherEndpointElementIndex)
    {
        assert(!mConnectedEdgesBuffer[vertexElementIndex].contains(
            [edgeElementIndex](auto const & cs)
            {
                return cs.EdgeIndex == edgeElementIndex;
            }));

        mConnectedEdgesBuffer[vertexElementIndex].emplace_back(
            edgeElementIndex,
            otherEndpointElementIndex);
    }

    auto const & GetConnectedTriangles(ElementIndex vertexElementIndex) const
    {
        return mConnectedTrianglesBuffer[vertexElementIndex];
    }

    void AddConnectedTriangle(
        ElementIndex vertexElementIndex,
        ElementIndex triangleElementIndex,
        bool isAtOwner)
    {
        assert(!mConnectedTrianglesBuffer[vertexElementIndex].ConnectedTriangles.contains(
            [triangleElementIndex](auto const & ct)
            {
                return ct == triangleElementIndex;
            }));

        mConnectedTrianglesBuffer[vertexElementIndex].ConnectTriangle(
            triangleElementIndex,
            isAtOwner);
    }

private:

    Buffer<vec2f> mPositionBuffer;
    Buffer<ConnectedEdgesVector> mConnectedEdgesBuffer;
    Buffer<ConnectedTrianglesVector> mConnectedTrianglesBuffer;
};

}