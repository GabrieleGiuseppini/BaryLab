/***************************************************************************************
* Original Author:      Gabriele Giuseppini
* Created:              2020-05-15
* Copyright:            Gabriele Giuseppini  (https://github.com/GabrieleGiuseppini)
***************************************************************************************/
#include "NpcParticles.h"

#include "Log.h"

void NpcParticles::Add(
    float mass,
    float staticFriction,
    float kineticFriction,
    vec2f const & position,
    rgbaColor const & color)
{
    mPhysicalPropertiesBuffer.emplace_back(mass, staticFriction, kineticFriction);
    mPositionBuffer.emplace_back(position);
    mVelocityBuffer.emplace_back(vec2f::zero());
    mExternalForcesBuffer.emplace_back(vec2f::zero());
    mSpringForcesBuffer.emplace_back(vec2f::zero());
    mVoluntaryForcesBuffer.emplace_back(vec2f::zero());
    mVoluntarySuperimposedDisplacementBuffer.emplace_back(vec2f::zero());
    mVoluntaryVelocityBuffer.emplace_back(vec2f::zero());

    mRenderColorBuffer.emplace_back(color);

    ++mParticleCount;
}

void NpcParticles::Query(ElementIndex particleElementIndex) const
{
    LogMessage("ParticleIndex: ", particleElementIndex);
    LogMessage("P=", mPositionBuffer[particleElementIndex].toString(), " V=", mVelocityBuffer[particleElementIndex].toString());
}