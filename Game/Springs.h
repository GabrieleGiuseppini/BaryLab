/***************************************************************************************
* Original Author:      Gabriele Giuseppini
* Created:              2020-05-16
* Copyright:            Gabriele Giuseppini  (https://github.com/GabrieleGiuseppini)
***************************************************************************************/
#pragma once

#include "Physics.h"

#include <GameCore/Buffer.h>
#include <GameCore/Colors.h>
#include <GameCore/ElementContainer.h>
#include <GameCore/FixedSizeVector.h>

#include <cassert>

namespace Physics {

class Springs : public ElementContainer
{
public:

    /*
     * The endpoints of a edge.
     */
    struct Endpoints
    {
        ElementIndex PointAIndex;
        ElementIndex PointBIndex;

        Endpoints(
            ElementIndex pointAIndex,
            ElementIndex pointBIndex)
            : PointAIndex(pointAIndex)
            , PointBIndex(pointBIndex)
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
        Octant PointAOctant;
        Octant PointBOctant;

        EndpointOctants(
            Octant pointAOctant,
            Octant pointBOctant)
            : PointAOctant(pointAOctant)
            , PointBOctant(pointBOctant)
        {}
    };

    /*
     * The triangles that have an edge along this edge.
     */
    using TrianglesVector = FixedSizeVector<ElementIndex, 2>;

public:

    Springs(ElementCount elementCount)
        : ElementContainer(elementCount)
        //////////////////////////////////
        // Buffers
        //////////////////////////////////
        // Structure
        , mEndpointsBuffer(mBufferElementCount, mElementCount, Endpoints(NoneElementIndex, NoneElementIndex))
        , mEndpointOctantsBuffer(mBufferElementCount, mElementCount, EndpointOctants(0, 4))
        , mTrianglesBuffer(mBufferElementCount, mElementCount, TrianglesVector())
    {
    }

    Springs(Springs && other) = default;

    void Add(
        ElementIndex pointAIndex,
        ElementIndex pointBIndex,
        Octant pointAOctant,
        Octant pointBOctant,
        TrianglesVector const & triangles);

public:

    //
    // Structure
    //

    ElementIndex GetEndpointAIndex(ElementIndex springElementIndex) const noexcept
    {
        return mEndpointsBuffer[springElementIndex].PointAIndex;
    }

    ElementIndex GetEndpointBIndex(ElementIndex springElementIndex) const noexcept
    {
        return mEndpointsBuffer[springElementIndex].PointBIndex;
    }

    ElementIndex GetOtherEndpointIndex(
        ElementIndex springElementIndex,
        ElementIndex pointElementIndex) const
    {
        if (pointElementIndex == mEndpointsBuffer[springElementIndex].PointAIndex)
            return mEndpointsBuffer[springElementIndex].PointBIndex;
        else
        {
            assert(pointElementIndex == mEndpointsBuffer[springElementIndex].PointBIndex);
            return mEndpointsBuffer[springElementIndex].PointAIndex;
        }
    }

    Endpoints const * restrict GetEndpointsBuffer() const noexcept
    {
        return mEndpointsBuffer.data();
    }

    // Returns +1.0 if the edge is directed outward from the specified vertex;
    // otherwise, -1.0.
    float GetSpringDirectionFrom(
        ElementIndex springElementIndex,
        ElementIndex pointIndex) const
    {
        if (pointIndex == mEndpointsBuffer[springElementIndex].PointAIndex)
            return 1.0f;
        else
            return -1.0f;
    }

    vec2f const & GetEndpointAPosition(
        ElementIndex springElementIndex,
        Points const & points) const
    {
        return points.GetPosition(mEndpointsBuffer[springElementIndex].PointAIndex);
    }

    vec2f const & GetEndpointBPosition(
        ElementIndex springElementIndex,
        Points const & points) const
    {
        return points.GetPosition(mEndpointsBuffer[springElementIndex].PointBIndex);
    }

    //
    // Endpoint octants
    //

    Octant GetEndpointAOctant(ElementIndex springElementIndex) const
    {
        return mEndpointOctantsBuffer[springElementIndex].PointAOctant;
    }

    Octant GetEndpointBOctant(ElementIndex springElementIndex) const
    {
        return mEndpointOctantsBuffer[springElementIndex].PointBOctant;
    }

    Octant GetEndpointOctant(
        ElementIndex springElementIndex,
        ElementIndex pointElementIndex) const
    {
        if (pointElementIndex == GetEndpointAIndex(springElementIndex))
            return GetEndpointAOctant(springElementIndex);
        else
        {
            assert(pointElementIndex == GetEndpointBIndex(springElementIndex));
            return GetEndpointBOctant(springElementIndex);
        }
    }

    Octant GetOtherEndpointOctant(
        ElementIndex springElementIndex,
        ElementIndex pointElementIndex) const
    {
        if (pointElementIndex == GetEndpointAIndex(springElementIndex))
            return GetEndpointBOctant(springElementIndex);
        else
        {
            assert(pointElementIndex == GetEndpointBIndex(springElementIndex));
            return GetEndpointAOctant(springElementIndex);
        }
    }

    //
    // Triangles
    //

    inline auto const & GetTriangles(ElementIndex springElementIndex) const
    {
        return mTrianglesBuffer[springElementIndex];
    }

    inline void AddTriangle(
        ElementIndex springElementIndex,
        ElementIndex triangleElementIndex)
    {
        assert(mTrianglesBuffer[springElementIndex].contains(
            [triangleElementIndex](auto st)
            {
                return st == triangleElementIndex;
            }));

        mTrianglesBuffer[springElementIndex].push_back(triangleElementIndex);
    }

    inline void RemoveTriangle(
        ElementIndex springElementIndex,
        ElementIndex triangleElementIndex)
    {
        bool found = mTrianglesBuffer[springElementIndex].erase_first(triangleElementIndex);

        assert(found);
        (void)found;
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
};

}