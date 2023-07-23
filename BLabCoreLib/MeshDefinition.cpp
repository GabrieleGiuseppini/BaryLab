/***************************************************************************************
* Original Author:		Gabriele Giuseppini
* Created:				2020-05-16
* Copyright:			Gabriele Giuseppini  (https://github.com/GabrieleGiuseppini)
***************************************************************************************/
#include "MeshDefinition.h"

#include "ImageFileTools.h"

#include <cassert>

MeshDefinition MeshDefinition::Load(std::filesystem::path const & filepath)
{
    ImageData structuralImage = ImageFileTools::LoadImageRgb(filepath);

    return MeshDefinition(
        std::move(structuralImage),
        filepath.stem().string());
}