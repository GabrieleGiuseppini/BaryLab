/***************************************************************************************
* Original Author:      Gabriele Giuseppini
* Created:              2018-05-13
* Copyright:            Gabriele Giuseppini  (https://github.com/GabrieleGiuseppini)
***************************************************************************************/
#pragma once

#include "BarycentricCoords.h"
#include "Buffer.h"
#include "Colors.h"
#include "ElementContainer.h"
#include "FixedSizeVector.h"
#include "Vertices.h"

#include <algorithm>
#include <array>
#include <cassert>
#include <functional>
#include <optional>

class Triangles : public ElementContainer
{
private:

    /*
     * The endpoints of a triangle, in CW order.
     */
    struct Endpoints
    {
        // A, B, C
        std::array<ElementIndex, 3u> VertexIndices;

        Endpoints(
            ElementIndex vertexAIndex,
            ElementIndex vertexBIndex,
            ElementIndex vertexCIndex)
            : VertexIndices({ vertexAIndex, vertexBIndex, vertexCIndex })
        {}
    };

    /*
     * The edges of a triangle, in CW order.
     */
    struct SubEdges
    {
        // A, B, C
        std::array<ElementIndex, 3u> EdgeIndices;

        SubEdges(
            ElementIndex subEdgeAIndex,
            ElementIndex subEdgeBIndex,
            ElementIndex subEdgeCIndex)
            : EdgeIndices({ subEdgeAIndex, subEdgeBIndex, subEdgeCIndex })
        {}
    };

public:

    Triangles(ElementCount elementCount)
        : ElementContainer(elementCount)
        //////////////////////////////////
        // Buffers
        //////////////////////////////////
        // Endpoints
        , mEndpointsBuffer(mBufferElementCount, mElementCount, Endpoints(NoneElementIndex, NoneElementIndex, NoneElementIndex))
        // Sub wdges
        , mSubEdgesBuffer(mBufferElementCount, mElementCount, SubEdges(NoneElementIndex, NoneElementIndex, NoneElementIndex))
    {
    }

    Triangles(Triangles && other) = default;

    void Add(
        ElementIndex vertexAIndex,
        ElementIndex vertexBIndex,
        ElementIndex vertexCIndex,
        ElementIndex subEdgeAIndex,
        ElementIndex subEdgeBIndex,
        ElementIndex subEdgeCIndex);

public:

    inline bool IsDeleted(ElementIndex triangleElementIndex) const
    {
        (void)triangleElementIndex;
        return false;
    }

    //
    // Endpoints
    //

    inline auto const & GetVertexIndices(ElementIndex triangleElementIndex) const
    {
        return mEndpointsBuffer[triangleElementIndex].VertexIndices;
    }

    inline ElementIndex GetVertexAIndex(ElementIndex triangleElementIndex) const
    {
        return mEndpointsBuffer[triangleElementIndex].VertexIndices[0];
    }

    inline ElementIndex GetVertexBIndex(ElementIndex triangleElementIndex) const
    {
        return mEndpointsBuffer[triangleElementIndex].VertexIndices[1];
    }

    inline ElementIndex GetVertexCIndex(ElementIndex triangleElementIndex) const
    {
        return mEndpointsBuffer[triangleElementIndex].VertexIndices[2];
    }

    inline bool AreVerticesInCwOrder(
        ElementIndex triangleElementIndex,
        ElementIndex vertex1Index,
        ElementIndex vertex2Index) const
    {
        return (GetVertexAIndex(triangleElementIndex) == vertex1Index && GetVertexBIndex(triangleElementIndex) == vertex2Index)
            || (GetVertexBIndex(triangleElementIndex) == vertex1Index && GetVertexCIndex(triangleElementIndex) == vertex2Index)
            || (GetVertexCIndex(triangleElementIndex) == vertex1Index && GetVertexAIndex(triangleElementIndex) == vertex2Index);
    }

    bool ContainsPoint(
        vec2f const & position,
        ElementIndex triangleElementIndex,
        Vertices const & vertices) const;

    ElementIndex FindContaining(
        vec2f const & position,
        Vertices const & vertices) const;

    bcoords3f ToBarycentricCoordinates(
        vec2f const & position,
        ElementIndex triangleElementIndex,
        Vertices const & vertices) const;

    bcoords3f ToBarycentricCoordinates(
        vec2f const & position,
        ElementIndex triangleElementIndex,
        Vertices const & vertices,
        float epsilon) const;

    bcoords3f ToBarycentricCoordinatesFromWithinTriangle(
        vec2f const & position,
        ElementIndex triangleElementIndex,
        Vertices const & vertices) const;

    vec2f FromBarycentricCoordinates(
        bcoords3f const & barycentricCoordinates,
        ElementIndex triangleElementIndex,
        Vertices const & vertices) const;

    //
    // Sub edges
    //

    auto const & GetSubEdges(ElementIndex triangleElementIndex) const
    {
        return mSubEdgesBuffer[triangleElementIndex];
    }

    inline ElementIndex GetSubEdgeAIndex(ElementIndex triangleElementIndex) const
    {
        return mSubEdgesBuffer[triangleElementIndex].EdgeIndices[0];
    }

    inline ElementIndex GetSubEdgeBIndex(ElementIndex triangleElementIndex) const
    {
        return mSubEdgesBuffer[triangleElementIndex].EdgeIndices[1];
    }

    inline ElementIndex GetSubEdgeCIndex(ElementIndex triangleElementIndex) const
    {
        return mSubEdgesBuffer[triangleElementIndex].EdgeIndices[2];
    }

    inline int GetSubEdgeOrdinal(
        ElementIndex triangleElementIndex,
        ElementIndex edgeElementIndex) const
    {
        if (mSubEdgesBuffer[triangleElementIndex].EdgeIndices[0] == edgeElementIndex)
            return 0;
        else if (mSubEdgesBuffer[triangleElementIndex].EdgeIndices[1] == edgeElementIndex)
            return 1;
        else
        {
            assert(mSubEdgesBuffer[triangleElementIndex].EdgeIndices[2] == edgeElementIndex);
            return 2;
        }
    }

    /*
     * Returns the vector representing the specified edge (ordinal), oriented
     * according to the triangle's point of view (thus CW).
     */
    inline vec2f GetSubEdgeVector(
        ElementIndex triangleElementIndex,
        int edgeOrdinal,
        Vertices const & vertices) const
    {
        ElementIndex const v2 = mEndpointsBuffer[triangleElementIndex].VertexIndices[(edgeOrdinal + 1) % 3];
        ElementIndex const v1 = mEndpointsBuffer[triangleElementIndex].VertexIndices[edgeOrdinal];

        return vertices.GetPosition(v2) - vertices.GetPosition(v1);
    }

private:

    inline vec2f InternalToBarycentricCoordinates(
        vec2f const & position,
        ElementIndex triangleElementIndex,
        Vertices const & vertices) const;

private:

    //////////////////////////////////////////////////////////
    // Buffers
    //////////////////////////////////////////////////////////

    // Endpoints
    Buffer<Endpoints> mEndpointsBuffer;

    // Sub edges
    Buffer<SubEdges> mSubEdgesBuffer;
};
