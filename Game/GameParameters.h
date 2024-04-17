/***************************************************************************************
* Original Author:      Gabriele Giuseppini
* Created:              2023-07-19
* Copyright:            Gabriele Giuseppini  (https://github.com/GabrieleGiuseppini)
***************************************************************************************/
#pragma once

#include <GameCore/Vectors.h>

/*
 * Parameters that affect the lab.
 */
struct GameParameters
{
    GameParameters();

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

    static size_t constexpr MaxSpringsPerPoint = 8u;
    static size_t constexpr MaxTrianglesPerPoint = 8u;
    static size_t constexpr MaxNpcs = 8192u;
    static size_t constexpr MaxParticlesPerNpc = 2u;

    //
    // Physics
    //

    static float constexpr SimulationTimeStepDuration = 1.0f / 64.0f;

    float GlobalDamping;
    static float constexpr MinGlobalDamping = 0.0f;
    static float constexpr MaxGlobalDamping = 1.0f;

    // TODO: the following 3 should probably be renamed as "Npc_" for the port

    float ElasticityAdjustment;
    static float constexpr MinElasticityAdjustment = 0.0f;
    static float constexpr MaxElasticityAdjustment = 4.0f;

    float StaticFrictionAdjustment;
    static float constexpr MinStaticFrictionAdjustment = 0.0f;
    static float constexpr MaxStaticFrictionAdjustment = 40.0f;

    float KineticFrictionAdjustment;
    static float constexpr MinKineticFrictionAdjustment = 0.0f;
    static float constexpr MaxKineticFrictionAdjustment = 40.0f;

    static float constexpr WaterFrictionDragCoefficient = 6.0f;

    float WaterFrictionDragCoefficientAdjustment;
    static float constexpr MinWaterFrictionDragCoefficientAdjustment = 0.0f;
    static float constexpr MaxWaterFrictionDragCoefficientAdjustment = 4.0f;

    float BuoyancyAdjustment;
    static float constexpr MinBuoyancyAdjustment = 0.0f;
    static float constexpr MaxBuoyancyAdjustment = 4.0f;

    // Misc

    float ToolSearchRadius;

    static float constexpr VertexRadius = 0.05f;
    static float constexpr EdgeThickness = 0.03f;
    static float constexpr SpringThickness = 0.06f;
    static float constexpr ParticleRadius = 0.15f;
    static float constexpr ParticleTrajectoryThickness = 0.04f;

    //
    // NPCs
    //

    float NpcSpringReductionFraction;
    static float constexpr MinNpcSpringReductionFraction = 0.0f;
    static float constexpr MaxNpcSpringReductionFraction = 2.0f;

    float NpcSpringDampingCoefficient;
    static float constexpr MinNpcSpringDampingCoefficient = 0.0f;
    static float constexpr MaxNpcSpringDampingCoefficient = 1.0f;

    float HumanNpcEquilibriumTorqueStiffnessCoefficient;
    static float constexpr MinHumanNpcEquilibriumTorqueStiffnessCoefficient = 0.0f;
    static float constexpr MaxHumanNpcEquilibriumTorqueStiffnessCoefficient = 0.01f;

    float HumanNpcEquilibriumTorqueDampingCoefficient;
    static float constexpr MinHumanNpcEquilibriumTorqueDampingCoefficient = 0.0f;
    static float constexpr MaxHumanNpcEquilibriumTorqueDampingCoefficient = 0.01f;

    static float constexpr HumanNpcWalkingSpeedBase = 1.0f; // m/s
    float HumanNpcWalkingSpeedAdjustment;
    static float constexpr MinHumanNpcWalkingSpeedAdjustment = 0.5f;
    static float constexpr MaxHumanNpcWalkingSpeedAdjustment = 2.5f;

    static float constexpr MaxHumanNpcTotalWalkingSpeedAdjustment = 4.0f; // Absolute cap

    struct HumanNpcGeometry
    {
        static float constexpr BodyLengthMean = 1.65f;
        static float constexpr BodyLengthStdDev = 0.065f;
        static float constexpr BodyWidthNarrowMultiplierStdDev = 0.05f;
        static float constexpr BodyWidthWideMultiplierStdDev = 0.15f;

        // All fractions below ar relative to BodyLength

        static float constexpr HeadLengthFraction = 1.0f / 8.0f;
        static float constexpr HeadWidthFraction = 1.0f / 10.0f;
        static float constexpr HeadDepthFraction = HeadWidthFraction;

        static float constexpr TorsoLengthFraction = 1.0f / 2.0f - HeadLengthFraction;
        static float constexpr TorsoWidthFraction = 1.0f / 7.0f;
        static float constexpr TorsoDepthFraction = 1.0f / 6.0f;

        static float constexpr ArmLengthFraction = 3.0f / 8.0f;
        static float constexpr ArmWidthFraction = 1.0f / 10.0f;
        static float constexpr ArmDepthFraction = ArmWidthFraction;

        static float constexpr LegLengthFraction = 1.0f / 2.0f;
        static float constexpr LegWidthFraction = 1.0f / 10.0f;
        static float constexpr LegDepthFraction = LegWidthFraction;

        static_assert(LegLengthFraction + TorsoLengthFraction + HeadLengthFraction == 1.0f);

        static float constexpr StepLengthFraction = 0.43f; // From foot to foot at longest separation
    };

    float HumanNpcBodyLengthAdjustment;
    static float constexpr MinHumanNpcBodyLengthAdjustment = 0.2f;
    static float constexpr MaxHumanNpcBodyLengthAdjustment = 10.0f;
};
