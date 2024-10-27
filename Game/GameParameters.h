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

    // Air
    static float constexpr AirMass = 1.2754f; // Kg
    static float constexpr AirPressureAtSeaLevel = 101325.0f; // No immediate relation to mass of air, due to compressibility of air; Pa == 1atm

    // Water
    static float constexpr WaterMass = 1000.0f; // Kg

    // Temperature at which all the constants are taken at
    static float constexpr Temperature0 = 298.15f; // 25C

    //
    // Structural constants
    //

    static float constexpr MaxWorldWidth = 10000.0f;
    static float constexpr HalfMaxWorldWidth = MaxWorldWidth / 2.0f;

    static float constexpr MaxWorldHeight = 22000.0f;
    static float constexpr HalfMaxWorldHeight = MaxWorldHeight / 2.0f;

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

    float AirDensityAdjustment;
    static float constexpr MinAirDensityAdjustment = 0.001f;
    static float constexpr MaxAirDensityAdjustment = 1000.0f;

    float WaterDensityAdjustment;
    static float constexpr MinWaterDensityAdjustment = 0.001f;
    static float constexpr MaxWaterDensityAdjustment = 100.0f;

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

    float OceanFloorElasticityCoefficient;

    // Heat and combustion

    float AirTemperature = 298.15f; // 25C, Kelvin

    static float constexpr AirThermalExpansionCoefficient = 0.0034f; // 1/K

    float WaterTemperature = 288.15f; // 15C, Kelvin

    static float constexpr WaterThermalExpansionCoefficient = 0.000207f; // 1/K

    float IgnitionTemperatureAdjustment = 1.0f;

    // Misc

    static float constexpr NpcAirBubbleFinalScale = 0.05f;

    static float constexpr VertexRadius = 0.05f;
    static float constexpr EdgeThickness = 0.03f;
    static float constexpr SpringThickness = 0.06f;
    static float constexpr ParticleRadius = 0.15f;
    static float constexpr ParticleTrajectoryThickness = 0.04f;

    float MoveToolInertia;
    bool IsUltraViolentMode;

    // Placeholders

    float AntiMatterBombImplosionStrength = 3.0f;

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

    float NpcFrictionAdjustment; // Both frictions; for fun
    static float constexpr MinNpcFrictionAdjustment = 0.0f;
    static float constexpr MaxNpcFrictionAdjustment = 4.0f;

    float NpcWindReceptivityAdjustment = 1.0f;

    float NpcSizeMultiplier;
    static float constexpr MinNpcSizeMultiplier = 0.2f;
    static float constexpr MaxNpcSizeMultiplier = 10.0f;

    static float constexpr HumanNpcTemperature = 310.15f; // 37 Celsius

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

    static float constexpr MaxNpcToolMoveVelocityMagnitude = 10.0f; // We cap velocity gained with move tool to prevent crazy spinning

    static float constexpr MaxHumanNpcWalkSinSlope = 0.87f; // Max sin of slope we're willing to climb up: ___*\___<W---  (60.5 degrees)

    struct HumanNpcGeometry
    {
        static float constexpr BodyLengthMean = 1.65f;
        static float constexpr BodyLengthStdDev = 0.065f;
        static float constexpr BodyWidthNarrowMultiplierStdDev = 0.045f;
        static float constexpr BodyWidthWideMultiplierStdDev = 0.15f;

        // All fractions below are relative to BodyLength
        //
        // Lengths are from Leonardo's Vitruvian man; subkinds may override
        static float constexpr HeadLengthFraction = 1.0f / 7.0f;
        static float constexpr QuadModeHeadWidthFraction = 1.0f / 8.0f;
        static float constexpr TorsoLengthFraction = 1.0f / 2.0f - HeadLengthFraction;
        static float constexpr QuadModeTorsoWidthFraction = 1.0f / 7.0f;
        static float constexpr ArmLengthFraction = 3.0f / 8.0f;
        static float constexpr QuadModeArmWidthFraction = 1.0f / 10.0f;
        static float constexpr LegLengthFraction = 1.0f / 2.0f;
        static float constexpr QuadModeLegWidthFraction = 1.0f / 10.0f;

        static_assert(LegLengthFraction + TorsoLengthFraction + HeadLengthFraction == 1.0f);

        static float constexpr StepLengthFraction = 0.43f; // From foot to foot at longest separation
    };

    static size_t constexpr NpcsPerGroup = 32; // The number of NPCs we add when we "add NPC group"
};
