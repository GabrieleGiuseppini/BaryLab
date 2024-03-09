/***************************************************************************************
* Original Author:      Gabriele Giuseppini
* Created:              2018-05-13
* Copyright:            Gabriele Giuseppini  (https://github.com/GabrieleGiuseppini)
***************************************************************************************/
#pragma once

#include "Physics.h"

#include <GameCore/BarycentricCoords.h>
#include <GameCore/Buffer.h>
#include <GameCore/Colors.h>
#include <GameCore/ElementContainer.h>
#include <GameCore/FixedSizeVector.h>
#include <GameCore/GameTypes.h>

#include <algorithm>
#include <array>
#include <cassert>
#include <functional>

namespace Physics {

class Triangles : public ElementContainer
{
private:

    /*
     * The endpoints of a triangle, in CW order.
     */
    struct Endpoints
    {
        // A, B, C
        std::array<ElementIndex, 3u> PointIndices;

        Endpoints(
            ElementIndex pointAIndex,
            ElementIndex pointBIndex,
            ElementIndex pointCIndex)
            : PointIndices({ pointAIndex, pointBIndex, pointCIndex })
        {}
    };

    /*
     * The edges of a triangle, in CW order.
     */
    struct SubSprings
    {
        // A, B, C
        std::array<ElementIndex, 3u> SpringIndices;

        SubSprings(
            ElementIndex subSpringAIndex,
            ElementIndex subSpringBIndex,
            ElementIndex subSpringCIndex)
            : SpringIndices({ subSpringAIndex, subSpringBIndex, subSpringCIndex })
        {}
    };

    using SubSpringSurfaceTypes = std::array<SurfaceType, 3>;

    /*
     * The opposite triangles of an edge, by edge ordinal.
     */

    struct OppositeTriangleInfo
    {
        ElementIndex TriangleElementIndex;
        int EdgeOrdinal;

        OppositeTriangleInfo(
            ElementIndex triangleElementIndex,
            int edgeOrdinal)
            : TriangleElementIndex(triangleElementIndex)
            , EdgeOrdinal(edgeOrdinal)
        {}
    };

    using OppositeTrianglesInfo = std::array<OppositeTriangleInfo, 3>;

public:

    Triangles(ElementCount elementCount)
        : ElementContainer(elementCount)
        //////////////////////////////////
        // Buffers
        //////////////////////////////////
        // Container
        , mEndpointsBuffer(mBufferElementCount, mElementCount, Endpoints(NoneElementIndex, NoneElementIndex, NoneElementIndex))
        , mSubSpringsBuffer(mBufferElementCount, mElementCount, SubSprings(NoneElementIndex, NoneElementIndex, NoneElementIndex))
        , mOppositeTrianglesBuffer(mBufferElementCount, mElementCount, { OppositeTriangleInfo(NoneElementIndex, -1), OppositeTriangleInfo(NoneElementIndex, -1), OppositeTriangleInfo(NoneElementIndex, -1) })
        , mSubSpringSurfaceTypesBuffer(mBufferElementCount, mElementCount, {SurfaceType::Open, SurfaceType::Open, SurfaceType::Open})
    {
    }

    Triangles(Triangles && other) = default;

    void Add(
        ElementIndex pointAIndex,
        ElementIndex pointBIndex,
        ElementIndex pointCIndex,
        ElementIndex subSpringAIndex,
        ElementIndex subSpringBIndex,
        ElementIndex subSpringCIndex,
        ElementIndex subSpringAOppositeTriangle,
        int subSpringAOppositeTriangleEdgeOrdinal,
        ElementIndex subSpringBOppositeTriangle,
        int subSpringBOppositeTriangleEdgeOrdinal,
        ElementIndex subSpringCOppositeTriangle,
        int subSpringCOppositeTriangleEdgeOrdinal,
        SurfaceType subSpringASurfaceType,
        SurfaceType subSpringBSurfaceType,
        SurfaceType subSpringCSurfaceType);

public:

    inline bool IsDeleted(ElementIndex triangleElementIndex) const
    {
        (void)triangleElementIndex;
        return false;
    }

    //
    // Endpoints
    //

    inline auto const & GetPointIndices(ElementIndex triangleElementIndex) const
    {
        return mEndpointsBuffer[triangleElementIndex].PointIndices;
    }

    inline ElementIndex GetPointAIndex(ElementIndex triangleElementIndex) const
    {
        return mEndpointsBuffer[triangleElementIndex].PointIndices[0];
    }

    inline ElementIndex GetPointBIndex(ElementIndex triangleElementIndex) const
    {
        return mEndpointsBuffer[triangleElementIndex].PointIndices[1];
    }

    inline ElementIndex GetPointCIndex(ElementIndex triangleElementIndex) const
    {
        return mEndpointsBuffer[triangleElementIndex].PointIndices[2];
    }

    inline bool ArePointsInCwOrder(
        ElementIndex triangleElementIndex,
        ElementIndex point1Index,
        ElementIndex point2Index) const
    {
        return (GetPointAIndex(triangleElementIndex) == point1Index && GetPointBIndex(triangleElementIndex) == point2Index)
            || (GetPointBIndex(triangleElementIndex) == point1Index && GetPointCIndex(triangleElementIndex) == point2Index)
            || (GetPointCIndex(triangleElementIndex) == point1Index && GetPointAIndex(triangleElementIndex) == point2Index);
    }

    bool ContainsPoint(
        vec2f const & position,
        ElementIndex triangleElementIndex,
        Points const & points) const;

    ElementIndex FindContaining(
        vec2f const & position,
        Points const & points) const;

    bcoords3f ToBarycentricCoordinates(
        vec2f const & position,
        ElementIndex triangleElementIndex,
        Points const & points) const;

    bcoords3f ToBarycentricCoordinates(
        vec2f const & position,
        ElementIndex triangleElementIndex,
        Points const & points,
        float epsilon) const;

    bcoords3f ToBarycentricCoordinatesFromWithinTriangle(
        vec2f const & position,
        ElementIndex triangleElementIndex,
        Points const & points) const;

    vec2f FromBarycentricCoordinates(
        bcoords3f const & barycentricCoordinates,
        ElementIndex triangleElementIndex,
        Points const & points) const;

    //
    // Sub springs
    //

    auto const & GetSubSprings(ElementIndex triangleElementIndex) const
    {
        return mSubSpringsBuffer[triangleElementIndex];
    }

    inline ElementIndex GetSubSpringAIndex(ElementIndex triangleElementIndex) const
    {
        return mSubSpringsBuffer[triangleElementIndex].SpringIndices[0];
    }

    inline ElementIndex GetSubSpringBIndex(ElementIndex triangleElementIndex) const
    {
        return mSubSpringsBuffer[triangleElementIndex].SpringIndices[1];
    }

    inline ElementIndex GetSubSpringCIndex(ElementIndex triangleElementIndex) const
    {
        return mSubSpringsBuffer[triangleElementIndex].SpringIndices[2];
    }

    inline int GetSubSpringOrdinal(
        ElementIndex triangleElementIndex,
        ElementIndex springElementIndex) const
    {
        if (mSubSpringsBuffer[triangleElementIndex].SpringIndices[0] == springElementIndex)
            return 0;
        else if (mSubSpringsBuffer[triangleElementIndex].SpringIndices[1] == springElementIndex)
            return 1;
        else
        {
            assert(mSubSpringsBuffer[triangleElementIndex].SpringIndices[2] == springElementIndex);
            return 2;
        }
    }

    /*
     * Returns the vector representing the specified edge (ordinal), oriented
     * according to the triangle's point of view (thus CW).
     */
    inline vec2f GetSubSpringVector(
        ElementIndex triangleElementIndex,
        int springOrdinal,
        Points const & points) const
    {
        assert(springOrdinal >= 0 && springOrdinal < 3);

        ElementIndex const v2 = mEndpointsBuffer[triangleElementIndex].PointIndices[(springOrdinal + 1) % 3];
        ElementIndex const v1 = mEndpointsBuffer[triangleElementIndex].PointIndices[springOrdinal];

        return points.GetPosition(v2) - points.GetPosition(v1);
    }

    // Opposite triangles

    OppositeTrianglesInfo const & GetOppositeTriangles(ElementIndex triangleElementIndex) const
    {
        return mOppositeTrianglesBuffer[triangleElementIndex];
    }

    OppositeTriangleInfo const & GetOppositeTriangle(
        ElementIndex triangleElementIndex,
        int springOrdinal) const
    {
        assert(springOrdinal >= 0 && springOrdinal < 3);
        return mOppositeTrianglesBuffer[triangleElementIndex][springOrdinal];
    }

    // Surface types

    SurfaceType GetSubSpringSurfaceType(
        ElementIndex triangleElementIndex,
        int springOrdinal) const
    {
        assert(springOrdinal >= 0 && springOrdinal < 3);
        return mSubSpringSurfaceTypesBuffer[triangleElementIndex][springOrdinal];
    }

private:

    inline vec2f InternalToBarycentricCoordinates(
        vec2f const & position,
        ElementIndex triangleElementIndex,
        Points const & points) const;

private:

    //////////////////////////////////////////////////////////
    // Buffers
    //////////////////////////////////////////////////////////

    // Endpoints
    Buffer<Endpoints> mEndpointsBuffer;

    // Sub springs
    Buffer<SubSprings> mSubSpringsBuffer;

    // Opposite triangles
    Buffer<OppositeTrianglesInfo> mOppositeTrianglesBuffer;

    // Surface types
    Buffer<SubSpringSurfaceTypes> mSubSpringSurfaceTypesBuffer;
};

}