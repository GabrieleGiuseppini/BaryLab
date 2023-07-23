/***************************************************************************************
* Original Author:      Gabriele Giuseppini
* Created:              2020-05-16
* Copyright:            Gabriele Giuseppini  (https://github.com/GabrieleGiuseppini)
***************************************************************************************/
#pragma once

#include "Buffer.h"
#include "Colors.h"
#include "ElementContainer.h"
#include "FixedSizeVector.h"
#include "Vertices.h"

#include <cassert>

class Edges : public ElementContainer
{
public:

    /*
     * The endpoints of a edge.
     */
    struct Endpoints
    {
        ElementIndex VertexAIndex;
        ElementIndex VertexBIndex;

        Endpoints(
            ElementIndex vertexAIndex,
            ElementIndex vertexBIndex)
            : VertexAIndex(vertexAIndex)
            , VertexBIndex(vertexBIndex)
        {}
    };

    /*
     * The factory angle of the edge from the point of view
     * of each endpoint.
     *
     * Angle 0 is E, angle 1 is SE, ..., angle 7 is NE.
     */
    struct EndpointOctants
    {
        Octant VertexAOctant;
        Octant VertexBOctant;

        EndpointOctants(
            Octant vertexAOctant,
            Octant vertexBOctant)
            : VertexAOctant(vertexAOctant)
            , VertexBOctant(vertexBOctant)
        {}
    };

    /*
     * The triangles that have an edge along this edge.
     */
    using TrianglesVector = FixedSizeVector<ElementIndex, 2>;

public:

    Edges(ElementCount elementCount)
        : ElementContainer(elementCount)
        //////////////////////////////////
        // Buffers
        //////////////////////////////////
        // Structure
        , mEndpointsBuffer(mBufferElementCount, mElementCount, Endpoints(NoneElementIndex, NoneElementIndex))
        , mEndpointOctantsBuffer(mBufferElementCount, mElementCount, EndpointOctants(0, 4))
        , mTrianglesBuffer(mBufferElementCount, mElementCount, TrianglesVector())
        // Render
        , mRenderColorBuffer(mBufferElementCount, mElementCount, vec4f::zero())
    {
    }

    Edges(Edges && other) = default;

    void Add(
        ElementIndex vertexAIndex,
        ElementIndex vertexBIndex,
        Octant vertexAOctant,
        Octant vertexBOctant,
        TrianglesVector const & triangles);

public:

    //
    // Structure
    //

    ElementIndex GetEndpointAIndex(ElementIndex edgeElementIndex) const noexcept
    {
        return mEndpointsBuffer[edgeElementIndex].VertexAIndex;
    }

    ElementIndex GetEndpointBIndex(ElementIndex edgeElementIndex) const noexcept
    {
        return mEndpointsBuffer[edgeElementIndex].VertexBIndex;
    }

    ElementIndex GetOtherEndpointIndex(
        ElementIndex edgeElementIndex,
        ElementIndex vertexElementIndex) const
    {
        if (vertexElementIndex == mEndpointsBuffer[edgeElementIndex].VertexAIndex)
            return mEndpointsBuffer[edgeElementIndex].VertexBIndex;
        else
        {
            assert(vertexElementIndex == mEndpointsBuffer[edgeElementIndex].VertexBIndex);
            return mEndpointsBuffer[edgeElementIndex].VertexAIndex;
        }
    }

    Endpoints const * restrict GetEndpointsBuffer() const noexcept
    {
        return mEndpointsBuffer.data();
    }

    // Returns +1.0 if the edge is directed outward from the specified vertex;
    // otherwise, -1.0.
    float GetEdgeDirectionFrom(
        ElementIndex edgeElementIndex,
        ElementIndex vertexIndex) const
    {
        if (vertexIndex == mEndpointsBuffer[edgeElementIndex].VertexAIndex)
            return 1.0f;
        else
            return -1.0f;
    }

    vec2f const & GetEndpointAPosition(
        ElementIndex edgeElementIndex,
        Vertices const & vertices) const
    {
        return vertices.GetPosition(mEndpointsBuffer[edgeElementIndex].VertexAIndex);
    }

    vec2f const & GetEndpointBPosition(
        ElementIndex edgeElementIndex,
        Vertices const & vertices) const
    {
        return vertices.GetPosition(mEndpointsBuffer[edgeElementIndex].VertexBIndex);
    }

    //
    // Endpoint octants
    //

    Octant GetEndpointAOctant(ElementIndex edgeElementIndex) const
    {
        return mEndpointOctantsBuffer[edgeElementIndex].VertexAOctant;
    }

    Octant GetEndpointBOctant(ElementIndex edgeElementIndex) const
    {
        return mEndpointOctantsBuffer[edgeElementIndex].VertexBOctant;
    }

    Octant GetEndpointOctant(
        ElementIndex edgeElementIndex,
        ElementIndex vertexElementIndex) const
    {
        if (vertexElementIndex == GetEndpointAIndex(edgeElementIndex))
            return GetEndpointAOctant(edgeElementIndex);
        else
        {
            assert(vertexElementIndex == GetEndpointBIndex(edgeElementIndex));
            return GetEndpointBOctant(edgeElementIndex);
        }
    }

    Octant GetOtherEndpointOctant(
        ElementIndex edgeElementIndex,
        ElementIndex vertexElementIndex) const
    {
        if (vertexElementIndex == GetEndpointAIndex(edgeElementIndex))
            return GetEndpointBOctant(edgeElementIndex);
        else
        {
            assert(vertexElementIndex == GetEndpointBIndex(edgeElementIndex));
            return GetEndpointAOctant(edgeElementIndex);
        }
    }

    //
    // Triangles
    //

    inline auto const & GetTriangles(ElementIndex edgeElementIndex) const
    {
        return mTrianglesBuffer[edgeElementIndex];
    }

    inline void AddTriangle(
        ElementIndex edgeElementIndex,
        ElementIndex triangleElementIndex)
    {
        assert(mTrianglesBuffer[edgeElementIndex].contains(
            [triangleElementIndex](auto st)
            {
                return st == triangleElementIndex;
            }));

        mTrianglesBuffer[edgeElementIndex].push_back(triangleElementIndex);
    }

    inline void RemoveTriangle(
        ElementIndex edgeElementIndex,
        ElementIndex triangleElementIndex)
    {
        bool found = mTrianglesBuffer[edgeElementIndex].erase_first(triangleElementIndex);

        assert(found);
        (void)found;
    }

    //
    // Render
    //

    vec4f const & GetRenderColor(ElementIndex edgeElementIndex) const
    {
        return mRenderColorBuffer[edgeElementIndex];
    }

private:

    //////////////////////////////////////////////////////////
    // Buffers
    //////////////////////////////////////////////////////////

    //
    // Structure
    //

    Buffer<Endpoints> mEndpointsBuffer;
    Buffer<EndpointOctants> mEndpointOctantsBuffer;
    Buffer<TrianglesVector> mTrianglesBuffer;

    //
    // Render
    //

    Buffer<vec4f> mRenderColorBuffer;
};
