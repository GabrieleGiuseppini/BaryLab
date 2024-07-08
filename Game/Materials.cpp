/***************************************************************************************
 * Original Author:		Gabriele Giuseppini
 * Created:				2020-05-16
 * Copyright:			Gabriele Giuseppini  (https://github.com/GabrieleGiuseppini)
 ***************************************************************************************/
#include "Materials.h"

#include <GameCore/GameException.h>
#include <GameCore/Utils.h>

StructuralMaterial StructuralMaterial::Create(
    MaterialColorKey const & colorKey,
    picojson::object const & structuralMaterialJson)
{
    std::string const name = Utils::GetMandatoryJsonMember<std::string>(structuralMaterialJson, "name");

    try
    {
        float const strength = Utils::GetMandatoryJsonMember<float>(structuralMaterialJson, "strength");

        picojson::object massJson = Utils::GetMandatoryJsonObject(structuralMaterialJson, "mass");
        float const nominalMass = Utils::GetMandatoryJsonMember<float>(massJson, "nominal_mass");
        float const density = Utils::GetMandatoryJsonMember<float>(massJson, "density");
        float const buoyancyVolumeFill = Utils::GetMandatoryJsonMember<float>(structuralMaterialJson, "buoyancy_volume_fill");
        float const stiffness = Utils::GetOptionalJsonMember<float>(structuralMaterialJson, "stiffness", 1.0);
        float const strainThresholdFraction = Utils::GetOptionalJsonMember<float>(structuralMaterialJson, "strain_threshold_fraction", 0.5f);
        float const elasticityCoefficient = Utils::GetOptionalJsonMember<float>(structuralMaterialJson, "elasticity_coefficient", 0.5f);
        float const kineticFrictionCoefficient = Utils::GetOptionalJsonMember<float>(structuralMaterialJson, "friction_kinetic_coefficient", 0.25f);
        float const staticFrictionCoefficient = Utils::GetOptionalJsonMember<float>(structuralMaterialJson, "friction_static_coefficient", 0.25f);

        float const opacity = Utils::GetOptionalJsonMember<float>(structuralMaterialJson, "opacity", 1.0f);

        // Water

        bool const isHull = Utils::GetMandatoryJsonMember<bool>(structuralMaterialJson, "is_hull");
        float const waterIntake = Utils::GetOptionalJsonMember<float>(structuralMaterialJson, "water_intake", 1.0);
        float const waterDiffusionSpeed = Utils::GetOptionalJsonMember<float>(structuralMaterialJson, "water_diffusion_speed", 0.5f);
        float const waterRetention = Utils::GetOptionalJsonMember<float>(structuralMaterialJson, "water_retention", 0.05f);
        float const rustReceptivity = Utils::GetOptionalJsonMember<float>(structuralMaterialJson, "rust_receptivity", 1.0);

        // Heat

        float const ignitionTemperature = Utils::GetMandatoryJsonMember<float>(structuralMaterialJson, "ignition_temperature");
        float const meltingTemperature = Utils::GetMandatoryJsonMember<float>(structuralMaterialJson, "melting_temperature");
        float const thermalConductivity = Utils::GetOptionalJsonMember<float>(structuralMaterialJson, "thermal_conductivity", 50.0f);
        float const thermalExpansionCoefficient = Utils::GetOptionalJsonMember<float>(structuralMaterialJson, "thermal_expansion_coefficient", 0.0);
        float const specificHeat = Utils::GetOptionalJsonMember<float>(structuralMaterialJson, "specific_heat", 100.0f);
        float const explosiveCombustionRadius = Utils::GetOptionalJsonMember<float>(structuralMaterialJson, "explosive_combustion_radius", 0.0);
        float const explosiveCombustionStrength = Utils::GetOptionalJsonMember<float>(structuralMaterialJson, "explosive_combustion_strength", 1.0);

        // Misc

        float const windReceptivity = Utils::GetMandatoryJsonMember<float>(structuralMaterialJson, "wind_receptivity");
        float const waterReactivity = Utils::GetOptionalJsonMember<float>(structuralMaterialJson, "water_reactivity", 0.0f);
        bool isLegacyElectrical = Utils::GetOptionalJsonMember<bool>(structuralMaterialJson, "is_legacy_electrical", false);
        bool isExemptFromPalette = Utils::GetOptionalJsonMember<bool>(structuralMaterialJson, "is_exempt_from_palette", false);

        // NPC-specific

        float const npcSpringReductionFraction = Utils::GetOptionalJsonMember<float>(structuralMaterialJson, "npc_spring_reduction_fraction", 0.97f);
        float const npcSpringDampingCoefficient = Utils::GetOptionalJsonMember<float>(structuralMaterialJson, "npc_spring_damping_coefficient", 0.5f);

        return StructuralMaterial(
            colorKey,
            name,
            rgbaColor(colorKey, static_cast<rgbaColor::data_type>(255.0f * opacity)), // renderColor
            strength,
            nominalMass,
            density,
            buoyancyVolumeFill,
            stiffness,
            strainThresholdFraction,
            elasticityCoefficient,
            kineticFrictionCoefficient,
            staticFrictionCoefficient,
            opacity,
            isHull,
            waterIntake,
            waterDiffusionSpeed,
            waterRetention,
            rustReceptivity,
            // Heat
            ignitionTemperature,
            meltingTemperature,
            thermalConductivity,
            thermalExpansionCoefficient,
            specificHeat,
            explosiveCombustionRadius,
            explosiveCombustionStrength,
            // Misc
            windReceptivity,
            waterReactivity,
            isLegacyElectrical,
            // NPC-specific
            npcSpringReductionFraction,
            npcSpringDampingCoefficient);
    }
    catch (GameException const & ex)
    {
        throw GameException(std::string("Error parsing structural material \"") + name + "\": " + ex.what());
    }
}
