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

    std::map<NpcSubKindIdType, HumanKind> humanKinds;
    std::map<NpcSubKindIdType, FurnitureKind> furnitureKinds;

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

        NpcSubKindIdType nextKindId = 0;
        auto const humanKindsArray = Utils::GetMandatoryJsonArray(humansObject, "kinds");
        for (auto const & humanKindArrayElement : humanKindsArray)
        {
            if (!humanKindArrayElement.is<picojson::object>())
            {
                throw GameException("Human NPC kind array element is not a JSON object");
            }

            HumanKind kind = ParseHumanKind(
                humanKindArrayElement.get<picojson::object>(),
                headMaterial,
                feetMaterial,
                globalHeadParticleAttributes,
                globalFeetParticleAttributes);

            humanKinds.try_emplace(nextKindId, std::move(kind));
            ++nextKindId;
        }
    }

    //
    // Furniture
    //

    {
        auto const furnitureObject = Utils::GetMandatoryJsonObject(rootObject, "furniture");

        NpcSubKindIdType nextKindId = 0;
        auto const furnitureKindsArray = Utils::GetMandatoryJsonArray(furnitureObject, "kinds");
        for (auto const & furnitureKindArrayElement : furnitureKindsArray)
        {
            if (!furnitureKindArrayElement.is<picojson::object>())
            {
                throw GameException("Furniture NPC kind array element is not a JSON object");
            }

            FurnitureKind kind = ParseFurnitureKind(
                furnitureKindArrayElement.get<picojson::object>(),
                materialDatabase);

            furnitureKinds.try_emplace(nextKindId, std::move(kind));
            ++nextKindId;
        }
    }

    //
    // Wrap it up
    //

    return NpcDatabase(
        std::move(humanKinds),
        std::move(furnitureKinds));
}

std::vector<std::tuple<NpcSubKindIdType, std::string>> NpcDatabase::GetHumanSubKinds(std::string const & language) const
{
    return GetSubKinds(mHumanKinds, language);
}

std::vector<std::tuple<NpcSubKindIdType, std::string>> NpcDatabase::GetFurnitureSubKinds(std::string const & language) const
{
    return GetSubKinds(mFurnitureKinds, language);
}

///////////////////////////////////////////////////////////////////////////////////////////////////////

NpcDatabase::HumanKind NpcDatabase::ParseHumanKind(
    picojson::object const & kindObject,
    StructuralMaterial const & headMaterial,
    StructuralMaterial const & feetMaterial,
    ParticleAttributesType const & globalHeadParticleAttributes,
    ParticleAttributesType const & globalFeetParticleAttributes)
{
    MultiLingualText name = ParseMultilingualText(kindObject, "name");
    NpcHumanRoleType const role = StrToNpcHumanRoleType(Utils::GetMandatoryJsonMember<std::string>(kindObject, "role"));
    ParticleAttributesType headParticleAttributes = MakeParticleAttributes(kindObject, "head_particle_attributes_overrides", globalHeadParticleAttributes);
    ParticleAttributesType feetParticleAttributes = MakeParticleAttributes(kindObject, "feet_particle_attributes_overrides", globalFeetParticleAttributes);
    float const sizeMultiplier = Utils::GetOptionalJsonMember<float>(kindObject, "size_multiplier", 1.0f);

    auto const dimensions = HumanDimensionsType({
        1.0f,
        // More or less from Vitruvian Man
        8.0f / 8.0f,
        8.0f / 21.0f,
        8.0f / 30.0f,
        2.0f / 10.0f });

    HumanTextureFramesType humanTextureFrames({
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

    float const bodyWidthRandomizationSensitivity = Utils::GetOptionalJsonMember<float>(kindObject, "body_width_randomization_sensitivity", 1.0f);

    return HumanKind({
        std::move(name),
        role,
        headMaterial,
        feetMaterial,
        {feetParticleAttributes, headParticleAttributes},
        sizeMultiplier,
        dimensions,
        bodyWidthRandomizationSensitivity,
        humanTextureFrames });
}


NpcDatabase::FurnitureKind NpcDatabase::ParseFurnitureKind(
    picojson::object const & kindObject,
    MaterialDatabase const & materialDatabase)
{
    MultiLingualText name = ParseMultilingualText(kindObject, "name");

    StructuralMaterial const & material = materialDatabase.GetStructuralMaterial(
        Utils::GetMandatoryJsonMember<std::string>(kindObject, "material"));

    auto const & particleMeshObject = Utils::GetMandatoryJsonObject(kindObject, "particle_mesh");
    ParticleMeshKindType const particleMeshKind = StrToParticleMeshKindType(
        Utils::GetMandatoryJsonMember<std::string>(particleMeshObject, "kind"));

    int particleCount = 0;
    FurnitureDimensionsType dimensions{ 0.0f, 0.0f };
    switch (particleMeshKind)
    {
    case ParticleMeshKindType::Dipole:
    {
        particleCount = 2;

        dimensions = FurnitureDimensionsType({
            0,
            Utils::GetMandatoryJsonMember<float>(particleMeshObject, "height") });

        break;
    }

    case ParticleMeshKindType::Particle:
    {
        particleCount = 1;

        dimensions = FurnitureDimensionsType({
            0,
            0 });

        break;
    }

    case ParticleMeshKindType::Quad:
    {
        particleCount = 4;

        float const height = Utils::GetMandatoryJsonMember<float>(particleMeshObject, "height");
        float const width = height;

        dimensions = FurnitureDimensionsType({
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
            throw GameException("Invalid size of particle_attributes_overrides for furniture NPC \"" + name.Get("") + "\"");
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

    return FurnitureKind({
        std::move(name),
        material,
        std::move(particleAttributes),
        particleMeshKind,
        dimensions,
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

NpcDatabase::MultiLingualText NpcDatabase::ParseMultilingualText(
    picojson::object const & containerObject,
    std::string const & textName)
{
    std::map<std::string, std::string> valuesByLanguageMap;

    std::string defaultValue;

    auto const entriesArray = Utils::GetMandatoryJsonArray(containerObject, textName);
    for (auto const & entryElement : entriesArray)
    {
        if (!entryElement.is<picojson::object>())
        {
            throw GameException("Multi-lingual text element for text \"" + textName + "\" in NPC database is not a JSON object");
        }

        auto const & entry = entryElement.get<picojson::object>();

        std::string const & language = Utils::GetMandatoryJsonMember<std::string>(entry, "language");
        std::string const & value = Utils::GetMandatoryJsonMember<std::string>(entry, "value");

        if (language.empty())
        {
            throw GameException("Language \"" + language + "\" is invalid for text \"" + textName + "\" in NPC database");
        }

        if (value.empty())
        {
            throw GameException("Value \"" + value + "\" is invalid for text \"" + textName + "\" in NPC database");
        }

        auto const [_, isInserted] = valuesByLanguageMap.try_emplace(Utils::ToLower(language), value);
        if (!isInserted)
        {
            throw GameException("Language \"" + language + "\" is duplicated for text \"" + textName + "\" in NPC database");
        }

        // Check if this could be a default
        if (Utils::ToLower(language) == "en")
            defaultValue = value;
    }

    if (valuesByLanguageMap.size() == 0)
    {
        throw GameException("Multi-lingual text element for text \"" + textName + "\" in NPC database is empty");
    }

    // Ensure we have a default

    if (defaultValue.empty())
    {
        // Take first one arbitrarily
        defaultValue = valuesByLanguageMap.begin()->second;
    }

    auto const [_, isInserted] = valuesByLanguageMap.try_emplace("", defaultValue);
    assert(isInserted);

    // Wrap it up

    return MultiLingualText(std::move(valuesByLanguageMap));
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