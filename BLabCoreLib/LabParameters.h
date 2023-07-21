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

    static size_t constexpr MaxEdgesPerVertex = 8u;
    static size_t constexpr MaxTrianglesPerVertex = 8u;

    //
    // Simulation
    //

    static float constexpr SimulationTimeStepDuration = 1.0f / 64.0f;

    static float constexpr ParticleMass = 68.95f; // Kg

    float MassAdjustment;
    static float constexpr MinMassAdjustment = 0.0001f;
    static float constexpr MaxMassAdjustment = 1000.0f;

    float GravityAdjustment;
    static float constexpr MinGravityAdjustment = 0.0f;
    static float constexpr MaxGravityAdjustment = 1000.0f;

    static float constexpr VertexRadius = 0.05f;
    static float constexpr EdgeThickness = 0.03f;
    static float constexpr ParticleRadius = 0.15f;
};
