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
#include "Points.h"

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
        std::array<ElementIndex, 3u> PointIndices;

        Endpoints(
            ElementIndex pointAIndex,
            ElementIndex pointBIndex,
            ElementIndex pointCIndex)
            : PointIndices({ pointAIndex, pointBIndex, pointCIndex })
        {}
    };

    /*
     * The springs along the edges of a triangle, in CW order.
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

public:

    Triangles(ElementCount elementCount)
        : ElementContainer(elementCount)
        //////////////////////////////////
        // Buffers
        //////////////////////////////////
        // Endpoints
        , mEndpointsBuffer(mBufferElementCount, mElementCount, Endpoints(NoneElementIndex, NoneElementIndex, NoneElementIndex))
        // Sub springs
        , mSubSpringsBuffer(mBufferElementCount, mElementCount, SubSprings(NoneElementIndex, NoneElementIndex, NoneElementIndex))
    {
    }

    Triangles(Triangles && other) = default;

    void Add(
        ElementIndex pointAIndex,
        ElementIndex pointBIndex,
        ElementIndex pointCIndex,
        ElementIndex subSpringAIndex,
        ElementIndex subSpringBIndex,
        ElementIndex subSpringCIndex);

public:

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

private:

    //////////////////////////////////////////////////////////
    // Buffers
    //////////////////////////////////////////////////////////

    // Endpoints
    Buffer<Endpoints> mEndpointsBuffer;

    // Sub springs - the springs that have this triangle among their super triangles
    Buffer<SubSprings> mSubSpringsBuffer;
};
