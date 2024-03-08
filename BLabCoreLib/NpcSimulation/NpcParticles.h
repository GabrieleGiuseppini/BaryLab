/***************************************************************************************
* Original Author:      Gabriele Giuseppini
* Created:              2020-05-16
* Copyright:            Gabriele Giuseppini  (https://github.com/GabrieleGiuseppini)
***************************************************************************************/
#pragma once

#include "Buffer.h"
#include "Colors.h"
#include "ElementContainer.h"
#include "ElementIndexRangeIterator.h"
#include "GameParameters.h"
#include "GameTypes.h"
#include "Physics.h"
#include "Vectors.h"

#include <algorithm>
#include <cassert>
#include <optional>

namespace Physics {

class NpcParticles final : public ElementContainer
{
public:

    struct PhysicalProperties
    {
        float Mass;
        float StaticFriction;
        float KineticFriction;
        float Elasticity;
        float BuoyancyVolumeFill;

        PhysicalProperties(
            float mass,
            float staticFriction,
            float kineticFriction,
            float elasticity,
            float buoyancyVolumeFill)
            : Mass(mass)
            , StaticFriction(staticFriction)
            , KineticFriction(kineticFriction)
            , Elasticity(elasticity)
            , BuoyancyVolumeFill(buoyancyVolumeFill)
        {}
    };

    NpcParticles(ElementCount particleCount)
        : ElementContainer(particleCount)
        , mParticleCount(0)
        //////////////////////////////////
        // Buffers
        //////////////////////////////////
        // Physics
        , mPhysicalPropertiesBuffer(mBufferElementCount, particleCount, PhysicalProperties(0.0f, 0.0f, 0.0f, 0.0f, 0.0f))
        , mPositionBuffer(mBufferElementCount, particleCount, vec2f::zero())
        , mVelocityBuffer(mBufferElementCount, particleCount, vec2f::zero())
        , mEquilibriumTorqueBuffer(mBufferElementCount, particleCount, 0.0f)
        , mPreliminaryForcesBuffer(mBufferElementCount, particleCount, vec2f::zero())
        , mExternalForcesBuffer(mBufferElementCount, particleCount, vec2f::zero())
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
        float elasticity,
        float buoyancyVolumeFill,
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

    float const & GetEquilibriumTorque(ElementIndex particleElementIndex) const noexcept
    {
        return mEquilibriumTorqueBuffer[particleElementIndex];
    }

    void SetEquilibriumTorque(
        ElementIndex particleElementIndex,
        float const & value) noexcept
    {
        mEquilibriumTorqueBuffer[particleElementIndex] = value;
    }

    void ResetEquilibriumTorque()
    {
        mEquilibriumTorqueBuffer.fill(0.0f);
    }

    vec2f const & GetPreliminaryForces(ElementIndex particleElementIndex) const noexcept
    {
        return mPreliminaryForcesBuffer[particleElementIndex];
    }

    vec2f * GetPreliminaryForces() noexcept
    {
        return mPreliminaryForcesBuffer.data();
    }

    void SetPreliminaryForces(
        ElementIndex particleElementIndex,
        vec2f const & value) noexcept
    {
        mPreliminaryForcesBuffer[particleElementIndex] = value;
    }

    vec2f const & GetExternalForces(ElementIndex particleElementIndex) const noexcept
    {
        return mExternalForcesBuffer[particleElementIndex];
    }

    vec2f * GetExternalForces() noexcept
    {
        return mExternalForcesBuffer.data();
    }

    void SetExternalForces(
        ElementIndex particleElementIndex,
        vec2f const & value) noexcept
    {
        mExternalForcesBuffer[particleElementIndex] = value;
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
    Buffer<float> mEquilibriumTorqueBuffer;
    Buffer<vec2f> mPreliminaryForcesBuffer;
    Buffer<vec2f> mExternalForcesBuffer;

    //
    // Render
    //

    Buffer<rgbaColor> mRenderColorBuffer;
};

}