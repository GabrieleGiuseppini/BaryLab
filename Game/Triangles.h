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
#include <GameCore/GameGeometry.h>
#include <GameCore/GameTypes.h>

#include <algorithm>
#include <array>
#include <cassert>
#include <functional>
#include <optional>

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

    using SubSpringNpcFloorKinds = std::array<NpcFloorKindType, 3>;
    using SubSpringNpcFloorGeometries = std::array<NpcFloorGeometryType, 3>;

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
        , mSubSpringNpcFloorKindsBuffer(mBufferElementCount, mElementCount, { NpcFloorKindType::NotAFloor, NpcFloorKindType::NotAFloor, NpcFloorKindType::NotAFloor })
        , mSubSpringNpcFloorGeometriesBuffer(mBufferElementCount, mElementCount, { NpcFloorGeometryType::NotAFloor, NpcFloorGeometryType::NotAFloor, NpcFloorGeometryType::NotAFloor })
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
        std::tuple<ElementIndex, int> subSpringAOppositeTriangleInfo,
        std::tuple<ElementIndex, int> subSpringBOppositeTriangleInfo,
        std::tuple<ElementIndex, int> subSpringCOppositeTriangleInfo,
        std::tuple<NpcFloorKindType, NpcFloorGeometryType> subSpringAFloorInfo,
        std::tuple<NpcFloorKindType, NpcFloorGeometryType> subSpringBFloorInfo,
        std::tuple<NpcFloorKindType, NpcFloorGeometryType> subSpringCFloorInfo,
        std::optional<ElementIndex> coveredTraverseSpringIndex);

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

    inline bool AreVerticesInCwOrder(
        ElementIndex triangleElementIndex,
        Points const & points) const
    {
        auto const pa = points.GetPosition(GetPointAIndex(triangleElementIndex));
        auto const pb = points.GetPosition(GetPointBIndex(triangleElementIndex));
        auto const pc = points.GetPosition(GetPointCIndex(triangleElementIndex));

        return (pb.x - pa.x) * (pc.y - pa.y) - (pc.x - pa.x) * (pb.y - pa.y) < 0;
    }

    ElementIndex FindContaining(
        vec2f const & position,
        Points const & points) const;

    //
    // Barycentric coordinates
    //

    bcoords3f ToBarycentricCoordinates(
        vec2f const & position,
        ElementIndex triangleElementIndex,
        Points const & points) const noexcept
    {
        vec2f const abBaryCoords = InternalToBarycentricCoordinates<2, 0, 1>(
            position,
            triangleElementIndex,
            points);

        return bcoords3f(
            abBaryCoords.x,
            abBaryCoords.y,
            1.0f - abBaryCoords.x - abBaryCoords.y);
    }

    bcoords3f ToBarycentricCoordinatesInsideEdge(
        vec2f const & position,
        ElementIndex triangleElementIndex,
        Points const & points,
        int insideEdge) const noexcept
    {
        //
        // Calculate bary coords enforcing that the coord wrt the specified edge
        // is not negative; to be used when we know that we're entering this
        // triangle from that edge. Avoids the infamous "around edge" oscillations
        // that happen when we cross an edge and we re-cross it again ad infinitum
        // because of numerical slack issues.
        //
        // - Calculate coords using any of the not-that-edge vertices as anchors (=> that-edge is one of the two coords calc'd)
        // - Then clamp that-edge and calc 3rd coord via 1-...
        //

        if (insideEdge == 0)
        {
            // Vertex is 2

            vec2f v1v2BaryCoords = InternalToBarycentricCoordinates<0, 1, 2>(
                position,
                triangleElementIndex,
                points);

            v1v2BaryCoords.y = std::max(v1v2BaryCoords.y, 0.0f);

            return bcoords3f(
                1.0f - v1v2BaryCoords.x - v1v2BaryCoords.y,
                v1v2BaryCoords.x,
                v1v2BaryCoords.y);
        }
        else if (insideEdge == 1)
        {
            // Vertex is 0

            vec2f v2v0BaryCoords = InternalToBarycentricCoordinates<1, 2, 0>(
                position,
                triangleElementIndex,
                points);

            v2v0BaryCoords.y = std::max(v2v0BaryCoords.y, 0.0f);

            return bcoords3f(
                v2v0BaryCoords.y,
                1.0f - v2v0BaryCoords.x - v2v0BaryCoords.y,
                v2v0BaryCoords.x);
        }
        else
        {
            assert(insideEdge == 2);

            // Vertex is 1

            vec2f v0v1BaryCoords = InternalToBarycentricCoordinates<2, 0, 1>(
                position,
                triangleElementIndex,
                points);

            v0v1BaryCoords.y = std::max(v0v1BaryCoords.y, 0.0f);

            return bcoords3f(
                v0v1BaryCoords.x,
                v0v1BaryCoords.y,
                1.0f - v0v1BaryCoords.x - v0v1BaryCoords.y);
        }
    }

    bcoords3f ToBarycentricCoordinatesFromWithinTriangle(
        vec2f const & position,
        ElementIndex triangleElementIndex,
        Points const & points) const noexcept
    {
        assert(Geometry::IsPointInTriangle(
            position,
            points.GetPosition(GetPointAIndex(triangleElementIndex)),
            points.GetPosition(GetPointBIndex(triangleElementIndex)),
            points.GetPosition(GetPointCIndex(triangleElementIndex))));

        vec2f abBaryCoords = InternalToBarycentricCoordinates<2, 0, 1>(
            position,
            triangleElementIndex,
            points).clamp(0.0f, 1.0f, 0.0f, 1.0f);

        return bcoords3f(
            abBaryCoords.x,
            abBaryCoords.y,
            1.0f - abBaryCoords.x - abBaryCoords.y);
    }

    vec2f FromBarycentricCoordinates(
        bcoords3f const & barycentricCoordinates,
        ElementIndex triangleElementIndex,
        Points const & points) const noexcept
    {
        vec2f const & positionA = points.GetPosition(mEndpointsBuffer[triangleElementIndex].PointIndices[0]);
        vec2f const & positionB = points.GetPosition(mEndpointsBuffer[triangleElementIndex].PointIndices[1]);
        vec2f const & positionC = points.GetPosition(mEndpointsBuffer[triangleElementIndex].PointIndices[2]);

        return
            positionA * barycentricCoordinates[0]
            + positionB * barycentricCoordinates[1]
            + positionC * barycentricCoordinates[2];
    }

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

    // Floor types

    NpcFloorKindType GetSubSpringNpcFloorKind(
        ElementIndex triangleElementIndex,
        int springOrdinal) const
    {
        assert(springOrdinal >= 0 && springOrdinal < 3);
        return mSubSpringNpcFloorKindsBuffer[triangleElementIndex][springOrdinal];
    }

    NpcFloorGeometryType GetSubSpringNpcFloorGeometry(
        ElementIndex triangleElementIndex,
        int springOrdinal) const
    {
        assert(springOrdinal >= 0 && springOrdinal < 3);
        return mSubSpringNpcFloorGeometriesBuffer[triangleElementIndex][springOrdinal];
    }

private:

    template<unsigned int anchorVertex, unsigned int vertex1, unsigned int vertex2>
    inline vec2f InternalToBarycentricCoordinates(
        vec2f const & position,
        ElementIndex triangleElementIndex,
        Points const & points) const noexcept
    {
        vec2f const & positionAnchor = points.GetPosition(mEndpointsBuffer[triangleElementIndex].PointIndices[anchorVertex]);

        vec2f const v1 = points.GetPosition(mEndpointsBuffer[triangleElementIndex].PointIndices[vertex1]) - positionAnchor;
        vec2f const v2 = points.GetPosition(mEndpointsBuffer[triangleElementIndex].PointIndices[vertex2]) - positionAnchor;
        vec2f const vp = position - positionAnchor;

        float const denominator = v2.y * v1.x - v2.x * v1.y;

        if (IsAlmostZero(denominator))
        {
            // Co-linear, put arbitrarily in center
            float constexpr l = 0.3333333f;
            return vec2f(l, l);
        }
        else
        {
            // See also: https://gamedev.stackexchange.com/questions/23743/whats-the-most-efficient-way-to-find-barycentric-coordinates

            float const l1 =
                ((v2.y * vp.x) - (v2.x * vp.y))
                / denominator;

            float const l2 =
                ((v1.x * vp.y) - (v1.y * vp.x))
                / denominator;

            return vec2f(
                l1,
                l2);
        }
    }

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

    // NPC Floor types
    Buffer<SubSpringNpcFloorKinds> mSubSpringNpcFloorKindsBuffer;
    Buffer<SubSpringNpcFloorGeometries> mSubSpringNpcFloorGeometriesBuffer;
};

}