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
    static constexpr vec2f GravityDir = vec2f(0.0f, -1.0f);
    static float constexpr GravityMagnitude = 9.80f; // m/s
    static constexpr vec2f Gravity = GravityDir * GravityMagnitude;

    //
    // Structural constants
    //

    static size_t constexpr MaxEdgesPerVertex = 8u;
    static size_t constexpr MaxTrianglesPerVertex = 8u;
    static size_t constexpr MaxNpcs = 32u;
    static size_t constexpr MaxParticlesPerNpc = 2u;

    //
    // Simulation
    //

    static float constexpr SimulationTimeStepDuration = 1.0f / 64.0f;

    float MassAdjustment;
    static float constexpr MinMassAdjustment = 0.0001f;
    static float constexpr MaxMassAdjustment = 1000.0f;

    float GravityAdjustment;
    static float constexpr MinGravityAdjustment = 0.0f;
    static float constexpr MaxGravityAdjustment = 1000.0f;

    float GlobalDamping;
    static float constexpr MinGlobalDamping = 0.0f;
    static float constexpr MaxGlobalDamping = 1.0f;

    float SpringReductionFraction;
    static float constexpr MinSpringReductionFraction = 0.0f;
    static float constexpr MaxSpringReductionFraction = 2.0f;

    float SpringDampingCoefficient;
    static float constexpr MinSpringDampingCoefficient = 0.0f;
    static float constexpr MaxSpringDampingCoefficient = 1.0f;

    float Elasticity;
    static float constexpr MinElasticity = 0.0f;
    static float constexpr MaxElasticity = 1.0f;

    float StaticFrictionAdjustment;
    static float constexpr MinStaticFrictionAdjustment = 0.0f;
    static float constexpr MaxStaticFrictionAdjustment = 40.0f;

    float KineticFrictionAdjustment;
    static float constexpr MinKineticFrictionAdjustment = 0.0f;
    static float constexpr MaxKineticFrictionAdjustment = 40.0f;

    // Misc

    static float constexpr VertexRadius = 0.05f;
    static float constexpr EdgeThickness = 0.03f;
    static float constexpr SpringThickness = 0.06f;
    static float constexpr ParticleRadius = 0.15f;
    static float constexpr ParticleTrajectoryThickness = 0.04f;

    //
    // NPCs
    //

    static float constexpr HumanNpcLength = 1.65f;

    float HumanNpcEquilibriumTorqueStiffnessCoefficient;
    static float constexpr MinHumanNpcEquilibriumTorqueStiffnessCoefficient = 0.0f;
    static float constexpr MaxHumanNpcEquilibriumTorqueStiffnessCoefficient = 0.01f;

    float HumanNpcEquilibriumTorqueDampingCoefficient;
    static float constexpr MinHumanNpcEquilibriumTorqueDampingCoefficient = 0.0f;
    static float constexpr MaxHumanNpcEquilibriumTorqueDampingCoefficient = 0.1f;
    
    float HumanNpcWalkingSpeed;
    static float constexpr MinHumanNpcWalkingSpeed = 0.0f;
    static float constexpr MaxHumanNpcWalkingSpeed = 2.5f;
};
