/***************************************************************************************
* Original Author:      Gabriele Giuseppini
* Created:              2023-07-20
* Copyright:            Gabriele Giuseppini  (https://github.com/GabrieleGiuseppini)
***************************************************************************************/
#include "Physics.h"

#include "Log.h"

namespace Physics {

void Vertices::Add(vec2f const & position)
{
    mPositionBuffer.emplace_back(position);
    mConnectedEdgesBuffer.emplace_back();
    mConnectedTrianglesBuffer.emplace_back();
}

void Vertices::Query(ElementIndex vertexElementIndex) const
{
    LogMessage("VertexIndex: ", vertexElementIndex);
    LogMessage("P=", mPositionBuffer[vertexElementIndex].toString());
    LogMessage("Edges: ", mConnectedEdgesBuffer[vertexElementIndex].size());
    LogMessage("Triangles: ", mConnectedTrianglesBuffer[vertexElementIndex].ConnectedTriangles.size());
}

}