/***************************************************************************************
* Original Author:		Gabriele Giuseppini
* Created:				2020-05-16
* Copyright:			Gabriele Giuseppini  (https://github.com/GabrieleGiuseppini)
***************************************************************************************/
#pragma once

#include "ImageData.h"

#include <filesystem>
#include <string>

/*
* The complete definition of a mesh.
*/
struct MeshDefinition
{
public:

    RgbImageData StructuralLayerImage;
    std::string const MeshName;

    static MeshDefinition Load(std::filesystem::path const & filepath);

private:

    MeshDefinition(
        RgbImageData structuralLayerImage,
        std::string const & meshName)
        : StructuralLayerImage(std::move(structuralLayerImage))
        , MeshName(meshName)
    {
    }
};
