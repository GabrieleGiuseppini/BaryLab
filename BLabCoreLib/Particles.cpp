/***************************************************************************************
* Original Author:      Gabriele Giuseppini
* Created:              2020-05-15
* Copyright:            Gabriele Giuseppini  (https://github.com/GabrieleGiuseppini)
***************************************************************************************/
#include "Particles.h"

#include "Log.h"

void Particles::Add(
    vec2f const & position,
    rgbaColor const & color)
{
    mPositionBuffer.emplace_back(position);
    mVelocityBuffer.emplace_back(vec2f::zero());
    mAssignedForceBuffer.emplace_back(vec2f::zero());

    mRenderColorBuffer.emplace_back(color);
}

void Particles::Query(ElementIndex particleElementIndex) const
{
    LogMessage("ParticleIndex: ", particleElementIndex);
    LogMessage("P=", mPositionBuffer[particleElementIndex].toString(), " V=", mVelocityBuffer[particleElementIndex].toString());
}