/***************************************************************************************
* Original Author:      Gabriele Giuseppini
* Created:              2023-07-19
* Copyright:            Gabriele Giuseppini  (https://github.com/GabrieleGiuseppini)
***************************************************************************************/
#pragma once

#include "Vectors.h"

/*
 * Parameters that affect the lab.
 */
struct LabParameters
{
    LabParameters();

    //
    // Physical Constants
    //

    // Gravity
    static constexpr vec2f Gravity = vec2f(0.0f, -9.80f);
    static constexpr vec2f GravityNormalized = vec2f(0.0f, -1.0f);
    static float constexpr GravityMagnitude = 9.80f; // m/s

    //
    // Structural constants
    //

    static size_t constexpr MaxSpringsPerPoint = 8u + 1u; // 8 neighbours and 1 rope spring, when this is a rope endpoint
    static size_t constexpr MaxTrianglesPerPoint = 8u;
};
