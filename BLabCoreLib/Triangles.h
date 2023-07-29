/***************************************************************************************
* Original Author:      Gabriele Giuseppini
* Created:              2018-05-13
* Copyright:            Gabriele Giuseppini  (https://github.com/GabrieleGiuseppini)
***************************************************************************************/
#pragma once

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

    vec3f ToBarycentricCoordinates(
        vec2f const & position,
        ElementIndex triangleElementIndex,
        Vertices const & vertices) const;

    vec2f FromBarycentricCoordinates(
        vec3f const & barycentricCoordinates,
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

    inline ElementIndex GetSubEdgeIndex(
        ElementIndex triangleElementIndex,
        ElementIndex edgeIndex) const
    {
        if (mSubEdgesBuffer[triangleElementIndex].EdgeIndices[0] == edgeIndex)
            return 0;
        else if (mSubEdgesBuffer[triangleElementIndex].EdgeIndices[1] == edgeIndex)
            return 1;
        else
        {
            assert(mSubEdgesBuffer[triangleElementIndex].EdgeIndices[2] == edgeIndex);
            return 2;
        }
    }

private:

    //////////////////////////////////////////////////////////
    // Buffers
    //////////////////////////////////////////////////////////

    // Endpoints
    Buffer<Endpoints> mEndpointsBuffer;

    // Sub edges
    Buffer<SubEdges> mSubEdgesBuffer;
};
