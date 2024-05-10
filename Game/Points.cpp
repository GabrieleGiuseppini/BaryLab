/***************************************************************************************
* Original Author:      Gabriele Giuseppini
* Created:              2023-07-20
* Copyright:            Gabriele Giuseppini  (https://github.com/GabrieleGiuseppini)
***************************************************************************************/
#include "Physics.h"

#include <GameCore/Log.h>

#include <limits>

namespace Physics {

void Points::Add(vec2f const & position)
{
    mPositionBuffer.emplace_back(position);
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

}