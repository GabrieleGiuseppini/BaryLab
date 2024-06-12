/***************************************************************************************
 * Original Author:		Gabriele Giuseppini
 * Created:				2020-05-16
 * Copyright:			Gabriele Giuseppini  (https://github.com/GabrieleGiuseppini)
 ***************************************************************************************/
#include "Materials.h"

#include <GameCore/GameException.h>
#include <GameCore/Utils.h>

StructuralMaterial StructuralMaterial::Create(
    rgbaColor const & renderColor,
    picojson::object const & structuralMaterialJson)
{
    std::string const name = Utils::GetMandatoryJsonMember<std::string>(structuralMaterialJson, "name");

    try
    {
        bool const isHull = Utils::GetMandatoryJsonMember<bool>(structuralMaterialJson, "is_hull");

        return StructuralMaterial(
            name,
            renderColor,
            isHull);
    }
    catch (GameException const & ex)
    {
        throw GameException(std::string("Error parsing structural material \"") + name + "\": " + ex.what());
    }
}

NpcMaterial NpcMaterial::Create(picojson::object const & npcMaterialJson)
{
    std::string const strKind = Utils::GetMandatoryJsonMember<std::string>(npcMaterialJson, "kind");

    try
    {
        KindType const kind = StrToKindType(strKind);
        rgbColor const renderColor = Utils::Hex2RgbColor(Utils::GetMandatoryJsonMember<std::string>(npcMaterialJson, "render_color"));

        float const mass = Utils::GetMandatoryJsonMember<float>(npcMaterialJson, "mass");
        float const springReductionFraction = Utils::GetMandatoryJsonMember<float>(npcMaterialJson, "spring_reduction_fraction");
        float const springDampingCoefficient = Utils::GetMandatoryJsonMember<float>(npcMaterialJson, "spring_damping_coefficient");
        float const staticFriction = Utils::GetMandatoryJsonMember<float>(npcMaterialJson, "static_friction");
        float const kineticFriction = Utils::GetMandatoryJsonMember<float>(npcMaterialJson, "kinetic_friction");
        float const elasticity = Utils::GetMandatoryJsonMember<float>(npcMaterialJson, "elasticity");
        float const buoyancyVolumeFill = Utils::GetMandatoryJsonMember<float>(npcMaterialJson, "buoyancy_volume_fill");

        return NpcMaterial(
            strKind,
            kind,
            rgbaColor(renderColor, rgbaColor::data_type_max),
            mass,
            springReductionFraction,
            springDampingCoefficient,
            staticFriction,
            kineticFriction,
            elasticity,
            buoyancyVolumeFill);
    }
    catch (GameException const & ex)
    {
        throw GameException(std::string("Error parsing NPC material \"") + strKind + "\": " + ex.what());
    }
}

NpcMaterial::KindType NpcMaterial::StrToKindType(std::string const & strKind)
{
    if (strKind == "Furniture")
        return KindType::Furniture;
    else if (strKind == "HumanHead")
        return KindType::HumanHead;
    else if (strKind == "HumanFeet")
        return KindType::HumanFeet;
    else
        throw GameException("Unrecognized NPC kind \"" + strKind + "\"");
}