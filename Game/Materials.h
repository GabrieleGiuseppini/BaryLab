/***************************************************************************************
 * Original Author:		Gabriele Giuseppini
 * Created:				2020-05-16
 * Copyright:			Gabriele Giuseppini  (https://github.com/GabrieleGiuseppini)
 ***************************************************************************************/
#pragma once

#include <GameCore/Colors.h>
#include <GameCore/GameTypes.h>

#include <picojson.h>

#include <string>

struct StructuralMaterial
{

public:

    MaterialColorKey ColorKey;

    std::string Name;
    rgbaColor RenderColor;
    float Strength;
    float NominalMass;
    float Density;
    float BuoyancyVolumeFill;
    float Stiffness;
    float StrainThresholdFraction;
    float ElasticityCoefficient;
    float KineticFrictionCoefficient;
    float StaticFrictionCoefficient;

    float Opacity;

    // Water
    bool IsHull;
    float WaterIntake;
    float WaterDiffusionSpeed;
    float WaterRetention;
    float RustReceptivity;

    // Heat
    float IgnitionTemperature; // K
    float MeltingTemperature; // K
    float ThermalConductivity; // W/(m*K)
    float ThermalExpansionCoefficient; // 1/K
    float SpecificHeat; // J/(Kg*K)
    float ExplosiveCombustionRadius; // m
    float ExplosiveCombustionStrength; // adimensional

    // Misc
    float WindReceptivity;
    float WaterReactivity;
    bool IsLegacyElectrical;

    // NPC-specific
    float NpcSpringReductionFraction;
    float NpcSpringDampingCoefficient;

public:

    static StructuralMaterial Create(
        MaterialColorKey const & colorKey,
        picojson::object const & structuralMaterialJson);

    /*
     * Returns the mass of this particle, calculated assuming that the particle is a cubic meter
     * full of a quantity of material equal to the density; for example, an iron truss has a lower
     * density than solid iron.
     */
    float GetMass() const
    {
        return NominalMass * Density;
    }

    /*
     * Returns the heat capacity of the material, in J/K.
     */
    float GetHeatCapacity() const
    {
        return SpecificHeat * GetMass();
    }

    StructuralMaterial(
        MaterialColorKey const & colorKey,
        std::string name,
        rgbaColor const & renderColor,
        float strength,
        float nominalMass,
        float density,
        float buoyancyVolumeFill,
        float stiffness,
        float strainThresholdFraction,
        float elasticityCoefficient,
        float kineticFrictionCoefficient,
        float staticFrictionCoefficient,
        float opacity,
        // Water
        bool isHull,
        float waterIntake,
        float waterDiffusionSpeed,
        float waterRetention,
        float rustReceptivity,
        // Heat
        float ignitionTemperature,
        float meltingTemperature,
        float thermalConductivity,
        float thermalExpansionCoefficient,
        float specificHeat,
        float explosiveCombustionRadius,
        float explosiveCombustionStrength,
        // Misc
        float windReceptivity,
        float waterReactivity,
        bool isLegacyElectrical,
        // NPC-specific
        float npcSpringReductionFraction,
        float npcSpringDampingCoefficient)
        : ColorKey(colorKey)
        , Name(name)
        , RenderColor(renderColor)
        , Strength(strength)
        , NominalMass(nominalMass)
        , Density(density)
        , BuoyancyVolumeFill(buoyancyVolumeFill)
        , Stiffness(stiffness)
        , StrainThresholdFraction(strainThresholdFraction)
        , ElasticityCoefficient(elasticityCoefficient)
        , KineticFrictionCoefficient(kineticFrictionCoefficient)
        , StaticFrictionCoefficient(staticFrictionCoefficient)
        , Opacity(opacity)
        , IsHull(isHull)
        , WaterIntake(waterIntake)
        , WaterDiffusionSpeed(waterDiffusionSpeed)
        , WaterRetention(waterRetention)
        , RustReceptivity(rustReceptivity)
        , IgnitionTemperature(ignitionTemperature)
        , MeltingTemperature(meltingTemperature)
        , ThermalConductivity(thermalConductivity)
        , ThermalExpansionCoefficient(thermalExpansionCoefficient)
        , SpecificHeat(specificHeat)
        , ExplosiveCombustionRadius(explosiveCombustionRadius)
        , ExplosiveCombustionStrength(explosiveCombustionStrength)
        , WindReceptivity(windReceptivity)
        , WaterReactivity(waterReactivity)
        , IsLegacyElectrical(isLegacyElectrical)
        , NpcSpringReductionFraction(npcSpringReductionFraction)
        , NpcSpringDampingCoefficient(npcSpringDampingCoefficient)
    {}
};
