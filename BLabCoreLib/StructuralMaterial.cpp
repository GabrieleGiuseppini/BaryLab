/***************************************************************************************
 * Original Author:		Gabriele Giuseppini
 * Created:				2020-05-16
 * Copyright:			Gabriele Giuseppini  (https://github.com/GabrieleGiuseppini)
 ***************************************************************************************/
#include "StructuralMaterial.h"

#include "BLabException.h"

#include "Utils.h"

StructuralMaterial StructuralMaterial::Create(
    rgbaColor const & renderColor,
    picojson::object const & structuralMaterialJson)
{
    std::string const name = Utils::GetMandatoryJsonMember<std::string>(structuralMaterialJson, "name");

    try
    {
        float const mass = Utils::GetOptionalJsonMember<float>(structuralMaterialJson, "mass", 1.0f);
        float const staticFriction = Utils::GetOptionalJsonMember<float>(structuralMaterialJson, "static_friction", 1.0f);
        float const kineticFriction = Utils::GetOptionalJsonMember<float>(structuralMaterialJson, "kinetic_friction", 1.0f);

        std::string const surfaceStr = Utils::GetMandatoryJsonMember<std::string>(structuralMaterialJson, "surface");
        SurfaceType const surface = StrToSurfaceType(surfaceStr);

        return StructuralMaterial(
            name,
            renderColor,
            mass,
            staticFriction,
            kineticFriction,
            surface);
    }
    catch (BLabException const & ex)
    {
        throw BLabException(std::string("Error parsing structural material \"") + name + "\": " + ex.what());
    }
}