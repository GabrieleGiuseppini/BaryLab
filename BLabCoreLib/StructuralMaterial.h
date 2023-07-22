/***************************************************************************************
 * Original Author:		Gabriele Giuseppini
 * Created:				2020-05-16
 * Copyright:			Gabriele Giuseppini  (https://github.com/GabrieleGiuseppini)
 ***************************************************************************************/
#pragma once

#include "BLabTypes.h"

#include <picojson.h>

#include <string>

struct StructuralMaterial
{
public:

    std::string Name;

    SurfaceType Surface;

public:

    static StructuralMaterial Create(picojson::object const & structuralMaterialJson);

    StructuralMaterial(
        std::string name,
        SurfaceType surface)
        : Name(name)
        , Surface(surface)
    {}
};
