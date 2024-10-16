/***************************************************************************************
* Original Author:		Gabriele Giuseppini
* Created:				2024-07-13
* Copyright:			Gabriele Giuseppini  (https://github.com/GabrieleGiuseppini)
***************************************************************************************/
#include "NpcDatabase.h"

#include "ResourceLocator.h"

#include <GameCore/GameException.h>
#include <GameCore/Utils.h>

#include <cassert>

NpcDatabase NpcDatabase::Load(MaterialDatabase const & materialDatabase)
{
    ResourceLocator resourceLocator;
    picojson::value const root = Utils::ParseJSONFile(resourceLocator.GetNpcDatabaseFilePath());
    if (!root.is<picojson::object>())
    {
        throw GameException("NPC database is not a JSON object");
    }

    picojson::object const & rootObject = root.get<picojson::object>();

    std::map<NpcSubKindIdType, HumanSubKind> humanSubKinds;
    std::map<NpcSubKindIdType, FurnitureSubKind> furnitureSubKinds;

    //
    // Humans
    //

    {
        auto const humansObject = Utils::GetMandatoryJsonObject(rootObject, "humans");

        auto const humansGlobalObject = Utils::GetMandatoryJsonObject(humansObject, "global");

        StructuralMaterial const & headMaterial = materialDatabase.GetStructuralMaterial(
            Utils::GetMandatoryJsonMember<std::string>(humansGlobalObject, "head_material"));

        StructuralMaterial const & feetMaterial = materialDatabase.GetStructuralMaterial(
            Utils::GetMandatoryJsonMember<std::string>(humansGlobalObject, "feet_material"));

        ParticleAttributesType globalHeadParticleAttributes = MakeParticleAttributes(humansGlobalObject, "head_particle_attributes_overrides", MakeDefaultParticleAttributes(headMaterial));

        ParticleAttributesType globalFeetParticleAttributes = MakeParticleAttributes(humansGlobalObject, "feet_particle_attributes_overrides", MakeDefaultParticleAttributes(feetMaterial));

        NpcSubKindIdType nextSubKindId = 0;
        auto const humanSubKindsArray = Utils::GetMandatoryJsonArray(humansObject, "sub_kinds");
        for (auto const & humanSubKindArrayElement : humanSubKindsArray)
        {
            if (!humanSubKindArrayElement.is<picojson::object>())
            {
                throw GameException("Human NPC kind array element is not a JSON object");
            }

            HumanSubKind subKind = ParseHumanSubKind(
                humanSubKindArrayElement.get<picojson::object>(),
                headMaterial,
                feetMaterial,
                globalHeadParticleAttributes,
                globalFeetParticleAttributes);

            humanSubKinds.try_emplace(nextSubKindId, std::move(subKind));
            ++nextSubKindId;
        }
    }

    //
    // Furniture
    //

    {
        auto const furnitureObject = Utils::GetMandatoryJsonObject(rootObject, "furniture");

        NpcSubKindIdType nextSubKindId = 0;
        auto const furnitureSubKindsArray = Utils::GetMandatoryJsonArray(furnitureObject, "sub_kinds");
        for (auto const & furnitureSubKindArrayElement : furnitureSubKindsArray)
        {
            if (!furnitureSubKindArrayElement.is<picojson::object>())
            {
                throw GameException("Furniture NPC kind array element is not a JSON object");
            }

            FurnitureSubKind subKind = ParseFurnitureSubKind(
                furnitureSubKindArrayElement.get<picojson::object>(),
                materialDatabase);

            furnitureSubKinds.try_emplace(nextSubKindId, std::move(subKind));
            ++nextSubKindId;
        }
    }

    //
    // Wrap it up
    //

    return NpcDatabase(
        std::move(humanSubKinds),
        std::move(furnitureSubKinds));
}

///////////////////////////////////////////////////////////////////////////////////////////////////////

NpcDatabase::NpcDatabase(
    std::map<NpcSubKindIdType, HumanSubKind> && humanSubKinds,
    std::map<NpcSubKindIdType, FurnitureSubKind> && furnitureSubKinds)
    : mHumanSubKinds(std::move(humanSubKinds))
    , mFurnitureSubKinds(std::move(furnitureSubKinds))
{
    // Build lookup tables

    mHumanSubKindIdsByRole.resize(static_cast<size_t>(NpcHumanRoleType::_Last) + 1);
    for (auto const & entry : mHumanSubKinds)
    {
        mHumanSubKindIdsByRole[static_cast<uint32_t>(entry.second.Role)].push_back(entry.first);
    }

    mFurnitureSubKindIdsByRole.resize(static_cast<size_t>(NpcFurnitureRoleType::_Last) + 1);
    for (auto const & entry : mFurnitureSubKinds)
    {
        mFurnitureSubKindIdsByRole[static_cast<uint32_t>(entry.second.Role)].push_back(entry.first);
    }
}

NpcDatabase::HumanSubKind NpcDatabase::ParseHumanSubKind(
    picojson::object const & subKindObject,
    StructuralMaterial const & headMaterial,
    StructuralMaterial const & feetMaterial,
    ParticleAttributesType const & globalHeadParticleAttributes,
    ParticleAttributesType const & globalFeetParticleAttributes)
{
    std::string const name = Utils::GetMandatoryJsonMember<std::string>(subKindObject, "name");
    NpcHumanRoleType const role = StrToNpcHumanRoleType(Utils::GetMandatoryJsonMember<std::string>(subKindObject, "role"));
    rgbColor const renderColor = Utils::Hex2RgbColor(Utils::GetMandatoryJsonMember<std::string>(subKindObject, "render_color"));

    ParticleAttributesType headParticleAttributes = MakeParticleAttributes(subKindObject, "head_particle_attributes_overrides", globalHeadParticleAttributes);
    ParticleAttributesType feetParticleAttributes = MakeParticleAttributes(subKindObject, "feet_particle_attributes_overrides", globalFeetParticleAttributes);

    float const sizeMultiplier = Utils::GetOptionalJsonMember<float>(subKindObject, "size_multiplier", 1.0f);
    float const bodyWidthRandomizationSensitivity = Utils::GetOptionalJsonMember<float>(subKindObject, "body_width_randomization_sensitivity", 1.0f);

    HumanTextureFramesType textureFrames({
        {-1.0f, 1.0f, -1.0f, 1.0f},
        {-1.0f, 1.0f, -1.0f, 1.0f},
        {-1.0f, 1.0f, -1.0f, 1.0f},
        {-1.0f, 1.0f, -1.0f, 1.0f},
        {-1.0f, 1.0f, -1.0f, 1.0f},
        {-1.0f, 1.0f, -1.0f, 1.0f},
        {-1.0f, 1.0f, -1.0f, 1.0f},
        {-1.0f, 1.0f, -1.0f, 1.0f},
        {-1.0f, 1.0f, -1.0f, 1.0f},
        {-1.0f, 1.0f, -1.0f, 1.0f},
        {-1.0f, 1.0f, -1.0f, 1.0f},
        {-1.0f, 1.0f, -1.0f, 1.0f}
        });

    auto const textureGeometry = HumanTextureGeometryType{
        // These are not used in BaryLab
        0.0f,
        0.0f,
        0.0f,
        0.0f,
        0.0f,
        0.0f,
        0.0f,
        0.0f
    };

    return HumanSubKind({
        std::move(name),
        role,
        renderColor,
        headMaterial,
        feetMaterial,
        {feetParticleAttributes, headParticleAttributes},
        sizeMultiplier,
        bodyWidthRandomizationSensitivity,
        textureFrames,
        textureGeometry });
}

NpcDatabase::FurnitureSubKind NpcDatabase::ParseFurnitureSubKind(
    picojson::object const & subKindObject,
    MaterialDatabase const & materialDatabase)
{
    std::string const name = Utils::GetMandatoryJsonMember<std::string>(subKindObject, "name");
    NpcFurnitureRoleType const role = StrToNpcFurnitureRoleType(Utils::GetMandatoryJsonMember<std::string>(subKindObject, "role"));
    rgbColor const renderColor = Utils::Hex2RgbColor(Utils::GetMandatoryJsonMember<std::string>(subKindObject, "render_color"));

    StructuralMaterial const & material = materialDatabase.GetStructuralMaterial(
        Utils::GetMandatoryJsonMember<std::string>(subKindObject, "material"));

    auto const & particleMeshObject = Utils::GetMandatoryJsonObject(subKindObject, "particle_mesh");
    ParticleMeshKindType const particleMeshKind = StrToParticleMeshKindType(
        Utils::GetMandatoryJsonMember<std::string>(particleMeshObject, "kind"));

    int particleCount = 0;
    FurnitureGeometryType geometry{ 0.0f, 0.0f };
    switch (particleMeshKind)
    {
        case ParticleMeshKindType::Dipole:
        {
            particleCount = 2;

            geometry = FurnitureGeometryType({
                0,
                Utils::GetMandatoryJsonMember<float>(particleMeshObject, "height") });

            break;
        }

        case ParticleMeshKindType::Particle:
        {
            particleCount = 1;

            geometry = FurnitureGeometryType({
                0,
                0 });

            break;
        }

        case ParticleMeshKindType::Quad:
        {
            particleCount = 4;

            float const height = Utils::GetMandatoryJsonMember<float>(particleMeshObject, "height");
            float const width = height;

            geometry = FurnitureGeometryType({
                width,
                height });

            break;
        }
    }

    ParticleAttributesType const defaultParticleAttributes = MakeDefaultParticleAttributes(material);
    std::vector<ParticleAttributesType> particleAttributes;
    auto const jsonParticleAttributesOverrides = Utils::GetOptionalJsonArray(particleMeshObject, "particle_attributes_overrides");
    if (jsonParticleAttributesOverrides.has_value())
    {
        if (jsonParticleAttributesOverrides->size() == particleCount)
        {
            // As-is
            for (auto const & bvfV : *jsonParticleAttributesOverrides)
            {
                particleAttributes.push_back(
                    MakeParticleAttributes(
                        Utils::GetJsonValueAs<picojson::object>(bvfV, "particle_attributes_overrides"),
                        defaultParticleAttributes));
            }
        }
        else if (jsonParticleAttributesOverrides->size() == 1)
        {
            // Repeat
            particleAttributes = std::vector<ParticleAttributesType>(
                particleCount,
                MakeParticleAttributes(
                    Utils::GetJsonValueAs<picojson::object>((*jsonParticleAttributesOverrides)[0], "particle_attributes_overrides"),
                    defaultParticleAttributes));
        }
        else
        {
            throw GameException("Invalid size of particle_attributes_overrides for furniture NPC \"" + name + "\"");
        }
    }
    else
    {
        // Use material's and defaults for all particles
        particleAttributes = std::vector<ParticleAttributesType>(
            particleCount,
            defaultParticleAttributes);
    }

    Render::TextureCoordinatesQuad textureCoordinatesQuad = Render::TextureCoordinatesQuad({ -1.0f, 1.0f, -1.0f, 1.0f });

    return FurnitureSubKind({
        std::move(name),
        role,
        renderColor,
        material,
        std::move(particleAttributes),
        particleMeshKind,
        geometry,
        std::move(textureCoordinatesQuad) });
}

NpcDatabase::ParticleAttributesType NpcDatabase::MakeParticleAttributes(
    picojson::object const & containerObject,
    std::string const & particleAttributesOverrideMemberName,
    ParticleAttributesType const & defaultParticleAttributes)
{
    auto const overridesJson = Utils::GetOptionalJsonObject(containerObject, particleAttributesOverrideMemberName);
    if (overridesJson.has_value())
    {
        return MakeParticleAttributes(*overridesJson, defaultParticleAttributes);
    }
    else
    {
        return defaultParticleAttributes;
    }
}

NpcDatabase::ParticleAttributesType NpcDatabase::MakeParticleAttributes(
    picojson::object const & particleAttributesOverrideJsonObject,
    ParticleAttributesType const & defaultParticleAttributes)
{
    float const buoyancyVolumeFill = Utils::GetOptionalJsonMember<float>(particleAttributesOverrideJsonObject, "buoyancy_volume_fill", defaultParticleAttributes.BuoyancyVolumeFill);
    float const springReductionFraction = Utils::GetOptionalJsonMember<float>(particleAttributesOverrideJsonObject, "spring_reduction_fraction", defaultParticleAttributes.SpringReductionFraction);
    float const springDampingCoefficient = Utils::GetOptionalJsonMember<float>(particleAttributesOverrideJsonObject, "spring_damping_coefficient", defaultParticleAttributes.SpringDampingCoefficient);

    return ParticleAttributesType{
        buoyancyVolumeFill,
        springReductionFraction,
        springDampingCoefficient
    };
}

NpcDatabase::ParticleAttributesType NpcDatabase::MakeDefaultParticleAttributes(StructuralMaterial const & baseMaterial)
{
    float constexpr DefaultSpringReductionFraction = 0.97f;
    float constexpr DefaultSpringDampingCoefficient = 0.5f * 0.906f;

    return ParticleAttributesType{
        baseMaterial.BuoyancyVolumeFill,
        DefaultSpringReductionFraction,
        DefaultSpringDampingCoefficient
    };
}

template<typename TNpcSubKindContainer>
std::vector<std::tuple<NpcSubKindIdType, std::string>> NpcDatabase::GetSubKinds(
    TNpcSubKindContainer const & container,
    std::string const & language)
{
    std::vector<std::tuple<NpcSubKindIdType, std::string>> kinds;

    for (auto const & it : container)
    {
        kinds.emplace_back(it.first, it.second.Name.Get(language));
    }

    return kinds;
}

NpcDatabase::ParticleMeshKindType NpcDatabase::StrToParticleMeshKindType(std::string const & str)
{
    if (Utils::CaseInsensitiveEquals(str, "Dipole"))
        return ParticleMeshKindType::Dipole;
    else if (Utils::CaseInsensitiveEquals(str, "Particle"))
        return ParticleMeshKindType::Particle;
    else if (Utils::CaseInsensitiveEquals(str, "Quad"))
        return ParticleMeshKindType::Quad;
    else
        throw GameException("Unrecognized ParticleMeshKindType \"" + str + "\"");
}