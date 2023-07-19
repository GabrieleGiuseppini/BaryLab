/***************************************************************************************
* Original Author:      Gabriele Giuseppini
* Created:              2018-05-13
* Copyright:            Gabriele Giuseppini  (https://github.com/GabrieleGiuseppini)
***************************************************************************************/
#include "Triangles.h"

void Triangles::Add(
    ElementIndex pointAIndex,
    ElementIndex pointBIndex,
    ElementIndex pointCIndex,
    ElementIndex subSpringAIndex,
    ElementIndex subSpringBIndex,
    ElementIndex subSpringCIndex)
{
    mEndpointsBuffer.emplace_back(pointAIndex, pointBIndex, pointCIndex);
    mSubSpringsBuffer.emplace_back(subSpringAIndex, subSpringBIndex, subSpringCIndex);
}
