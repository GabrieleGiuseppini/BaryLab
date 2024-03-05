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
* The complete definition of a mesh (ship).
*/
struct ShipDefinition
{
public:

    RgbImageData StructuralLayerImage;
    std::string const ShipName;

    static ShipDefinition Load(std::filesystem::path const & filepath);

private:

    ShipDefinition(
        RgbImageData structuralLayerImage,
        std::string const & shipName)
        : StructuralLayerImage(std::move(structuralLayerImage))
        , ShipName(shipName)
    {
    }
};
