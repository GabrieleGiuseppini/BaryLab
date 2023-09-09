/***************************************************************************************
* Original Author:      Gabriele Giuseppini
* Created:              2020-05-16
* Copyright:            Gabriele Giuseppini  (https://github.com/GabrieleGiuseppini)
***************************************************************************************/
#pragma once

#include "BLabTypes.h"
#include "Buffer.h"
#include "Colors.h"
#include "ElementContainer.h"
#include "ElementIndexRangeIterator.h"
#include "LabParameters.h"
#include "Vectors.h"

#include <algorithm>
#include <cassert>
#include <optional>

class Particles : public ElementContainer
{
public:

    struct StateType
    {
        struct ConstrainedStateType
        {
            ElementIndex CurrentTriangle;
            vec3f CurrentTriangleBarycentricCoords;
            vec2f MeshRelativeVelocity; // Velocity of particle (as in velocity buffer), but relative to mesh at the moment velocity was calculated

            ConstrainedStateType(
                ElementIndex currentTriangle,
                vec3f const & currentTriangleBarycentricCoords)
                : CurrentTriangle(currentTriangle)
                , CurrentTriangleBarycentricCoords(currentTriangleBarycentricCoords)
                , MeshRelativeVelocity(vec2f::zero())
            {}
        };

        std::optional<ConstrainedStateType> ConstrainedState;

        struct TrajectoryStateType
        {
            vec2f SourcePosition;
            vec2f TargetPosition;

            struct ConstrainedStateType
            {
                vec3f TargetPositionCurrentTriangleBarycentricCoords;
                vec2f MeshDisplacement;

                ConstrainedStateType(
                    vec3f const & targetPositionCurrentTriangleBarycentricCoords,
                    vec2f const & meshDisplacement)
                    : TargetPositionCurrentTriangleBarycentricCoords(targetPositionCurrentTriangleBarycentricCoords)
                    , MeshDisplacement(meshDisplacement)
                {}
            };

            std::optional<ConstrainedStateType> ConstrainedState; // Always set when in constrained state; updated when current triangle changes

            vec2f CurrentPosition;

            TrajectoryStateType(
                vec2f const & sourcePosition,
                vec2f const & targetPosition,
                std::optional<ConstrainedStateType> constrainedState)
                : SourcePosition(sourcePosition)
                , TargetPosition(targetPosition)
                , ConstrainedState(std::move(constrainedState))
                , CurrentPosition(sourcePosition)
            {}
        };

        std::optional<TrajectoryStateType> TrajectoryState; // When set, we have a trajectory target; when not set, we have to calculate a target

        StateType()
            : ConstrainedState()
            , TrajectoryState()
        {}
    };

public:

    Particles(ElementCount particleCount)
        : ElementContainer(particleCount)
        //////////////////////////////////
        // Buffers
        //////////////////////////////////
        // Physics
        , mPositionBuffer(mBufferElementCount, particleCount, vec2f::zero())
        , mVelocityBuffer(mBufferElementCount, particleCount, vec2f::zero())
        , mWorldForceBuffer(mBufferElementCount, particleCount, vec2f::zero())
        // Render
        , mRenderColorBuffer(mBufferElementCount, particleCount, rgbaColor::zero())
        // State
        , mStateBuffer(mBufferElementCount, particleCount, StateType())
    {
    }

    Particles(Particles && other) = default;

    void Add(
        vec2f const & position,
        rgbaColor const & color);

    void Query(ElementIndex particleElementIndex) const;

public:

    //
    // Physics
    //

    vec2f const & GetPosition(ElementIndex particleElementIndex) const noexcept
    {
        return mPositionBuffer[particleElementIndex];
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
        ElementIndex particleElementIndex,
        vec2f const & value) noexcept
    {
        mPositionBuffer[particleElementIndex] = value;
    }

    vec2f const & GetVelocity(ElementIndex particleElementIndex) const noexcept
    {
        return mVelocityBuffer[particleElementIndex];
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
        ElementIndex particleElementIndex,
        vec2f const & value) noexcept
    {
        mVelocityBuffer[particleElementIndex] = value;
    }

    vec2f const & GetWorldForce(ElementIndex particleElementIndex) const noexcept
    {
        return mWorldForceBuffer[particleElementIndex];
    }

    vec2f * GetWorldForceBuffer() noexcept
    {
        return mWorldForceBuffer.data();
    }

    void SetWorldForce(
        ElementIndex particleElementIndex,
        vec2f const & value) noexcept
    {
        mWorldForceBuffer[particleElementIndex] = value;
    }

    //
    // Render
    //

    rgbaColor const & GetRenderColor(ElementIndex particleElementIndex) const
    {
        return mRenderColorBuffer[particleElementIndex];
    }

    void SetRenderColor(
        ElementIndex particleElementIndex,
        rgbaColor const & color)
    {
        mRenderColorBuffer[particleElementIndex] = color;
    }

    rgbaColor const * GetRenderColorBuffer() const
    {
        return mRenderColorBuffer.data();
    }

    //
    // State
    //

    StateType const & GetState(ElementIndex particleElementIndex) const
    {
        return mStateBuffer[particleElementIndex];
    }

    StateType & GetState(ElementIndex particleElementIndex)
    {
        return mStateBuffer[particleElementIndex];
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
    Buffer<vec2f> mWorldForceBuffer;

    //
    // Render
    //

    Buffer<rgbaColor> mRenderColorBuffer;

    //
    // State
    //

    Buffer<StateType> mStateBuffer;
};
