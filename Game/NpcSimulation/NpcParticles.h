/***************************************************************************************
* Original Author:      Gabriele Giuseppini
* Created:              2020-05-16
* Copyright:            Gabriele Giuseppini  (https://github.com/GabrieleGiuseppini)
***************************************************************************************/
#pragma once

#include "GameParameters.h"
#include "Materials.h"

#include <GameCore/Buffer.h>
#include <GameCore/Colors.h>
#include <GameCore/ElementContainer.h>
#include <GameCore/ElementIndexRangeIterator.h>
#include <GameCore/GameTypes.h>
#include <GameCore/SysSpecifics.h>
#include <GameCore/Vectors.h>

#include <algorithm>
#include <cassert>
#include <optional>

namespace Physics {

class NpcParticles final : public ElementContainer
{
public:

    NpcParticles(ElementCount maxParticleCount)
        : ElementContainer(make_aligned_float_element_count(maxParticleCount))
        , mMaxParticleCount(maxParticleCount)
        //////////////////////////////////
        // Buffers
        //////////////////////////////////
        , mIsInUseBuffer(maxParticleCount, false)
        // Physics
        , mMaterialBuffer(maxParticleCount, nullptr)
        , mMassBuffer(maxParticleCount, 0.0f)
        , mBuoyancyFactorBuffer(maxParticleCount, 0.0f)
        , mPositionBuffer(maxParticleCount, vec2f::zero())
        , mVelocityBuffer(maxParticleCount, vec2f::zero())
        , mPreliminaryForcesBuffer(maxParticleCount, vec2f::zero())
        , mExternalForcesBuffer(maxParticleCount, vec2f::zero())
        // Render
        , mRenderColorBuffer(maxParticleCount, rgbaColor::zero())
        //////////////////////////////////
        // Container
        //////////////////////////////////
        , mParticleInUseCount(0)
        , mFreeParticleSearchStartIndex(0)
    {
    }

    NpcParticles(NpcParticles && other) = default;

    ElementCount GetRemainingParticlesCount() const
    {
        assert(mParticleInUseCount <= mMaxParticleCount);
        return mMaxParticleCount - mParticleInUseCount;
    }

    ElementIndex Add(
        float mass,
        float buoyancyFactor,
        StructuralMaterial const * material,
        vec2f const & position,
        rgbaColor const & color);

    void Remove(
        ElementIndex particleIndex);

    void Query(ElementIndex particleElementIndex) const;

public:

    //
    // Physics
    //

    StructuralMaterial const & GetMaterial(ElementIndex particleElementIndex) const noexcept
    {
        assert(mMaterialBuffer[particleElementIndex] != nullptr);
        return *(mMaterialBuffer[particleElementIndex]);
    }

    float const GetMass(ElementIndex particleElementIndex) const noexcept
    {
        return mMassBuffer[particleElementIndex];
    }

    void SetMass(
        ElementIndex particleElementIndex,
        float value) noexcept
    {
        mMassBuffer[particleElementIndex] = value;
    }

    float const GetBuoyancyFactor(ElementIndex particleElementIndex) const noexcept
    {
        return mBuoyancyFactorBuffer[particleElementIndex];
    }

    void SetBuoyancyFactor(
        ElementIndex particleElementIndex,
        float value) noexcept
    {
        mBuoyancyFactorBuffer[particleElementIndex] = value;
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

    ElementIndex FindFreeParticleIndex();

private:

    ElementCount const mMaxParticleCount;

    //////////////////////////////////////////////////////////
    // Buffers
    //////////////////////////////////////////////////////////

    // In use: true when the particle is occupied
    Buffer<bool> mIsInUseBuffer;

    //
    // Physics
    //

    Buffer<StructuralMaterial const *> mMaterialBuffer;
    Buffer<float> mMassBuffer; // Adjusted
    Buffer<float> mBuoyancyFactorBuffer; // Adjusted
    Buffer<vec2f> mPositionBuffer;
    Buffer<vec2f> mVelocityBuffer;
    Buffer<vec2f> mPreliminaryForcesBuffer;
    Buffer<vec2f> mExternalForcesBuffer;

    //
    // Render
    //

    Buffer<rgbaColor> mRenderColorBuffer;

    //////////////////////////////////////////////////////////
    // Container
    //////////////////////////////////////////////////////////

    // Convenience counter
    ElementCount mParticleInUseCount;

    // The index at which to start searching for free particles
    // (just an optimization over restarting from zero each time)
    ElementIndex mFreeParticleSearchStartIndex;
};

}