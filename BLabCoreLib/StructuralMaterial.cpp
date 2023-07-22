/***************************************************************************************
 * Original Author:		Gabriele Giuseppini
 * Created:				2020-05-16
 * Copyright:			Gabriele Giuseppini  (https://github.com/GabrieleGiuseppini)
 ***************************************************************************************/
#include "StructuralMaterial.h"

#include "BLabException.h"

#include "Utils.h"

StructuralMaterial StructuralMaterial::Create(picojson::object const & structuralMaterialJson)
{
    std::string const name = Utils::GetMandatoryJsonMember<std::string>(structuralMaterialJson, "name");

    try
    {
        std::string const surfaceStr = Utils::GetMandatoryJsonMember<std::string>(structuralMaterialJson, "surface");
        SurfaceType const surface = StrToSurfaceType(surfaceStr);

        return StructuralMaterial(
            name,
            surface);
    }
    catch (BLabException const & ex)
    {
        throw BLabException(std::string("Error parsing structural material \"") + name + "\": " + ex.what());
    }
}