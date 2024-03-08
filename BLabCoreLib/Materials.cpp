/***************************************************************************************
 * Original Author:		Gabriele Giuseppini
 * Created:				2020-05-16
 * Copyright:			Gabriele Giuseppini  (https://github.com/GabrieleGiuseppini)
 ***************************************************************************************/
#include "Materials.h"

#include "GameException.h"
#include "Utils.h"

StructuralMaterial StructuralMaterial::Create(
    rgbaColor const & renderColor,
    picojson::object const & structuralMaterialJson)
{
    std::string const name = Utils::GetMandatoryJsonMember<std::string>(structuralMaterialJson, "name");

    try
    {
        std::string const surfaceStr = Utils::GetMandatoryJsonMember<std::string>(structuralMaterialJson, "surface");
        SurfaceType const surface = StrToSurfaceType(surfaceStr);

        return StructuralMaterial(
            name,
            renderColor,
            surface);
    }
    catch (GameException const & ex)
    {
        throw GameException(std::string("Error parsing structural material \"") + name + "\": " + ex.what());
    }
}

NpcMaterial NpcMaterial::Create(picojson::object const & npcMaterialJson)
{
    std::string const name = Utils::GetMandatoryJsonMember<std::string>(npcMaterialJson, "name");

    try
    {
        KindType const kind = StrToKindType(Utils::GetMandatoryJsonMember<std::string>(npcMaterialJson, "kind"));
        rgbColor const renderColor = Utils::Hex2RgbColor(Utils::GetMandatoryJsonMember<std::string>(npcMaterialJson, "render_color"));

        float const mass = Utils::GetOptionalJsonMember<float>(npcMaterialJson, "mass", 1.0f);
        float const staticFriction = Utils::GetOptionalJsonMember<float>(npcMaterialJson, "static_friction", 1.0f);
        float const kineticFriction = Utils::GetOptionalJsonMember<float>(npcMaterialJson, "kinetic_friction", 1.0f);
        float const elasticity = Utils::GetOptionalJsonMember<float>(npcMaterialJson, "elasticity", 0.0f);
        float const buoyancyVolumeFill = Utils::GetOptionalJsonMember<float>(npcMaterialJson, "buoyancy_volume_fill", 1.0f);

        std::string const surfaceStr = Utils::GetMandatoryJsonMember<std::string>(npcMaterialJson, "surface");
        SurfaceType const surface = StrToSurfaceType(surfaceStr);

        return NpcMaterial(
            name,
            kind,
            rgbaColor(renderColor, rgbaColor::data_type_max),
            mass,
            staticFriction,
            kineticFriction,
            elasticity,
            buoyancyVolumeFill,
            surface);
    }
    catch (GameException const & ex)
    {
        throw GameException(std::string("Error parsing NPC material \"") + name + "\": " + ex.what());
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