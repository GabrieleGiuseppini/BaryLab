/***************************************************************************************
* Original Author:      Gabriele Giuseppini
* Created:              2020-05-15
* Copyright:            Gabriele Giuseppini  (https://github.com/GabrieleGiuseppini)
***************************************************************************************/
#include "Points.h"

#include "Log.h"

void Points::Add(
    vec2f const & position,
    vec3f const & color)
{
    mPositionBuffer.emplace_back(position);
    mVelocityBuffer.emplace_back(vec2f::zero());
    mAssignedForceBuffer.emplace_back(vec2f::zero());

    mConnectedSpringsBuffer.emplace_back();
    mConnectedTrianglesBuffer.emplace_back();

    mRenderColorBuffer.emplace_back(vec4f(color, 1.0f));
}

void Points::Query(ElementIndex pointElementIndex) const
{
    LogMessage("PointIndex: ", pointElementIndex);
    LogMessage("P=", mPositionBuffer[pointElementIndex].toString(), " V=", mVelocityBuffer[pointElementIndex].toString());
    LogMessage("Springs: ", mConnectedSpringsBuffer[pointElementIndex].size());
    LogMessage("Triangles: ", mConnectedTrianglesBuffer[pointElementIndex].ConnectedTriangles.size());
}