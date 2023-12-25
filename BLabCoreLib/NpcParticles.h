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

class NpcParticles final : public ElementContainer
{
public:

    struct PhysicalProperties
    {
        float Mass;
        float StaticFriction;
        float KineticFriction;

        PhysicalProperties(
            float mass,
            float staticFriction,
            float kineticFriction)
            : Mass(mass)
            , StaticFriction(staticFriction)
            , KineticFriction(kineticFriction)
        {}
    };

    NpcParticles(ElementCount particleCount)
        : ElementContainer(particleCount)
        , mParticleCount(0)
        //////////////////////////////////
        // Buffers
        //////////////////////////////////
        // Physics
        , mPhysicalPropertiesBuffer(mBufferElementCount, particleCount, PhysicalProperties(0.0f, 0.0f, 0.0f))
        , mPositionBuffer(mBufferElementCount, particleCount, vec2f::zero())
        , mVelocityBuffer(mBufferElementCount, particleCount, vec2f::zero())
        , mExternalForcesBuffer(mBufferElementCount, particleCount, vec2f::zero())
        , mSpringForcesBuffer(mBufferElementCount, particleCount, vec2f::zero())
        , mVoluntaryForcesBuffer(mBufferElementCount, particleCount, vec2f::zero())
        , mVoluntarySuperimposedDisplacementBuffer(mBufferElementCount, particleCount, vec2f::zero())
        , mVoluntaryVelocityBuffer(mBufferElementCount, particleCount, vec2f::zero())
        // Render
        , mRenderColorBuffer(mBufferElementCount, particleCount, rgbaColor::zero())
    {
    }

    NpcParticles(NpcParticles && other) = default;

    ElementCount GetParticleCount() const
    {
        return mParticleCount;
    }

    void Add(
        float mass,
        float staticFriction,
        float kineticFriction,
        vec2f const & position,
        rgbaColor const & color);

    void Query(ElementIndex particleElementIndex) const;

public:

    //
    // Physics
    //

    PhysicalProperties const & GetPhysicalProperties(ElementIndex particleElementIndex) const noexcept
    {
        return mPhysicalPropertiesBuffer[particleElementIndex];
    }

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

    vec2f const & GetExternalForces(ElementIndex particleElementIndex) const noexcept
    {
        return mExternalForcesBuffer[particleElementIndex];
    }

    vec2f * GetExternalForcesBuffer() noexcept
    {
        return mExternalForcesBuffer.data();
    }

    void SetExternalForces(
        ElementIndex particleElementIndex,
        vec2f const & value) noexcept
    {
        mExternalForcesBuffer[particleElementIndex] = value;
    }

    vec2f const & GetSpringForces(ElementIndex particleElementIndex) const noexcept
    {
        return mSpringForcesBuffer[particleElementIndex];
    }

    vec2f * GetSpringForcesBuffer() noexcept
    {
        return mSpringForcesBuffer.data();
    }

    void SetSpringForces(
        ElementIndex particleElementIndex,
        vec2f const & value) noexcept
    {
        mSpringForcesBuffer[particleElementIndex] = value;
    }

    vec2f const & GetVoluntaryForces(ElementIndex particleElementIndex) const noexcept
    {
        return mVoluntaryForcesBuffer[particleElementIndex];
    }

    vec2f * GetVoluntaryForcesBuffer() noexcept
    {
        return mVoluntaryForcesBuffer.data();
    }

    void SetVoluntaryForces(
        ElementIndex particleElementIndex,
        vec2f const & value) noexcept
    {
        mVoluntaryForcesBuffer[particleElementIndex] = value;
    }

    vec2f const & GetVoluntarySuperimposedDisplacement(ElementIndex particleElementIndex) const noexcept
    {
        return mVoluntarySuperimposedDisplacementBuffer[particleElementIndex];
    }

    vec2f * GetVoluntarySuperimposedDisplacementBuffer() noexcept
    {
        return mVoluntarySuperimposedDisplacementBuffer.data();
    }

    void SetVoluntarySuperimposedDisplacement(
        ElementIndex particleElementIndex,
        vec2f const & value) noexcept
    {
        mVoluntarySuperimposedDisplacementBuffer[particleElementIndex] = value;
    }

    vec2f const & GetVoluntaryVelocity(ElementIndex particleElementIndex) const noexcept
    {
        return mVoluntaryVelocityBuffer[particleElementIndex];
    }

    vec2f * GetVoluntaryVelocityBuffer() noexcept
    {
        return mVoluntaryVelocityBuffer.data();
    }

    void SetVoluntaryVelocity(
        ElementIndex particleElementIndex,
        vec2f const & value) noexcept
    {
        mVoluntaryVelocityBuffer[particleElementIndex] = value;
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

private:

    ElementCount mParticleCount;

    //////////////////////////////////////////////////////////
    // Buffers
    //////////////////////////////////////////////////////////

    //
    // Physics
    //

    Buffer<PhysicalProperties> mPhysicalPropertiesBuffer;
    Buffer<vec2f> mPositionBuffer;
    Buffer<vec2f> mVelocityBuffer;
    Buffer<vec2f> mExternalForcesBuffer;
    Buffer<vec2f> mSpringForcesBuffer; // Hooke+Damp
    Buffer<vec2f> mVoluntaryForcesBuffer;
    Buffer<vec2f> mVoluntarySuperimposedDisplacementBuffer;
    Buffer<vec2f> mVoluntaryVelocityBuffer;

    //
    // Render
    //

    Buffer<rgbaColor> mRenderColorBuffer;
};
