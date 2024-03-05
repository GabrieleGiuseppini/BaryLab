/***************************************************************************************
* Original Author:		Gabriele Giuseppini
* Created:				2020-05-16
* Copyright:			Gabriele Giuseppini  (https://github.com/GabrieleGiuseppini)
***************************************************************************************/
#include "ShipDefinition.h"

#include "ImageFileTools.h"

#include <cassert>

ShipDefinition ShipDefinition::Load(std::filesystem::path const & filepath)
{
    ImageData structuralImage = ImageFileTools::LoadImageRgb(filepath);

    return ShipDefinition(
        std::move(structuralImage),
        filepath.stem().string());
}