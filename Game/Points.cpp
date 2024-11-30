/***************************************************************************************
* Original Author:      Gabriele Giuseppini
* Created:              2023-07-20
* Copyright:            Gabriele Giuseppini  (https://github.com/GabrieleGiuseppini)
***************************************************************************************/
#include "Physics.h"

#include <GameCore/Log.h>

#include <limits>

namespace Physics {

void Points::Add(
    vec2f const & position,
    StructuralMaterial const * material)
{
    mPositionBuffer.emplace_back(position);
    mMaterialBuffer.emplace_back(material);
    mConnectedSpringsBuffer.emplace_back();
    mConnectedTrianglesBuffer.emplace_back();
}

bool Points::QueryAt(vec2f const & worldCoordinates) const
{
    ElementIndex bestPointIndex = NoneElementIndex;
    float bestDistance = std::numeric_limits<float>::max();
    for (auto p : *this)
    {
        float const distance = (GetPosition(p) - worldCoordinates).length();
        if (distance < 5.0f
            && distance < bestDistance)
        {
            bestPointIndex = p;
            bestDistance = distance;
        }
    }

    if (bestPointIndex != NoneElementIndex)
    {
        Query(bestPointIndex);
        return true;
    }

    return false;
}

void Points::Query(ElementIndex vertexElementIndex) const
{
    LogMessage("VertexIndex: ", vertexElementIndex);
    LogMessage("P=", mPositionBuffer[vertexElementIndex].toString());
    LogMessage("Springs: ", mConnectedSpringsBuffer[vertexElementIndex].size());
    for (auto s : mConnectedSpringsBuffer[vertexElementIndex])
        LogMessage("  ", s.SpringIndex);
    LogMessage("Triangles: ", mConnectedTrianglesBuffer[vertexElementIndex].ConnectedTriangles.size());
}

float Points::GetWater(ElementIndex /*pointElementIndex*/) const
{
    return 0.0f;
}

vec2f Points::GetWaterVelocity(ElementIndex /*pointElementIndex*/) const
{
    return vec2f::zero();
}

void Points::AddTransientAdditionalMass(
    ElementIndex /*pointElementIndex*/,
    float /*mass*/)
{
}

float Points::GetTemperature(ElementIndex /*pointElementIndex*/) const
{
    return GameParameters::Temperature0;
}

float Points::GetLight(ElementIndex /*pointElementIndex*/) const
{
    return 0.0f;
}

bool Points::GetIsHull(ElementIndex /*pointElementIndex*/) const
{
    return false;
}

bool Points::IsBurning(ElementIndex /*pointElementIndex*/) const
{
    return false;
}

bool Points::GetIsElectrified(ElementIndex /*pointElementIndex*/) const
{
    return false;
}

}