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
    static size_t constexpr MaxParticlesPerNpc = 4u;
    static size_t constexpr MaxSpringsPerNpc = 6u;
    static size_t constexpr MaxSpringsPerNpcParticle = 3u;

    //
    // Physics
    //

    template <typename T>
    static T constexpr SimulationStepTimeDuration = 1.0f / 64.0f;

    float GravityAdjustment;
    static float constexpr MinGravityAdjustment = 0.0f;
    static float constexpr MaxGravityAdjustment = 100.0f;

    float MassAdjustment;
    static float constexpr MinMassAdjustment = 0.0001f;
    static float constexpr MaxMassAdjustment = 100.0f;

    float GlobalDampingAdjustment;
    static float constexpr MinGlobalDampingAdjustment = 0.0f;
    static float constexpr MaxGlobalDampingAdjustment = 10.0f;

    float ElasticityAdjustment;
    static float constexpr MinElasticityAdjustment = 0.0f;
    static float constexpr MaxElasticityAdjustment = 10.0f;

    float StaticFrictionAdjustment;
    static float constexpr MinStaticFrictionAdjustment = 0.0f;
    static float constexpr MaxStaticFrictionAdjustment = 10.0f;

    float KineticFrictionAdjustment;
    static float constexpr MinKineticFrictionAdjustment = 0.0f;
    static float constexpr MaxKineticFrictionAdjustment = 10.0f;

    static float constexpr WaterFrictionDragCoefficient = 6.0f;

    float WaterFrictionDragAdjustment;
    static float constexpr MinWaterFrictionDragAdjustment = 0.0f;
    static float constexpr MaxWaterFrictionDragAdjustment = 4.0f;

    float BuoyancyAdjustment;
    static float constexpr MinBuoyancyAdjustment = 0.0f;
    static float constexpr MaxBuoyancyAdjustment = 4.0f;

    // Misc

    static float constexpr VertexRadius = 0.05f;
    static float constexpr EdgeThickness = 0.03f;
    static float constexpr SpringThickness = 0.06f;
    static float constexpr ParticleRadius = 0.15f;
    static float constexpr ParticleTrajectoryThickness = 0.04f;

    float MoveToolInertia;
    bool IsUltraViolentMode;

    //
    // NPCs
    //

    static float constexpr NpcDamping = 0.0078f;

    float NpcSpringReductionFractionAdjustment;
    static float constexpr MinNpcSpringReductionFractionAdjustment = 0.0f;
    static float constexpr MaxNpcSpringReductionFractionAdjustment = 5.0f;

    float NpcSpringDampingCoefficientAdjustment;
    static float constexpr MinNpcSpringDampingCoefficientAdjustment = 0.0f;
    static float constexpr MaxNpcSpringDampingCoefficientAdjustment = 2.0f;

    float NpcSizeMultiplier;
    static float constexpr MinNpcSizeMultiplier = 0.2f;
    static float constexpr MaxNpcSizeMultiplier = 10.0f;

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

    static float constexpr MaxHumanNpcTotalWalkingSpeedAdjustment = 3.5f; // Absolute cap

    static float constexpr MaxHumanNpcWalkSinSlope = 0.87f; // Max sin of slope we're willing to climb up: ___*\___<W---  (60.5 degrees)

    struct HumanNpcGeometry
    {
        static float constexpr BodyLengthMean = 1.65f;
        static float constexpr BodyLengthStdDev = 0.065f;
        static float constexpr BodyWidthNarrowMultiplierStdDev = 0.045f;
        static float constexpr BodyWidthWideMultiplierStdDev = 0.15f;

        // All fractions below are relative to BodyLength
        static float constexpr HeadWidthFraction = 1.0f / 8.0f; // Our DB has head as a square
        static float constexpr TorsoLengthFraction = 2.0f / 5.0f; // Width then depends on texture frame
        static float constexpr ArmLengthFraction = 2.0f / 5.0f; // Width then depends on texture frame
        static float constexpr LegLengthFraction = 19.0f / 40.0f; // Width then depends on texture frame

        static_assert(LegLengthFraction + TorsoLengthFraction + HeadWidthFraction  == 1.0f);

        static float constexpr StepLengthFraction = 0.43f; // From foot to foot at longest separation
    };
};
