/***************************************************************************************
 * Original Author:		Gabriele Giuseppini
 * Created:				2020-05-16
 * Copyright:			Gabriele Giuseppini  (https://github.com/GabrieleGiuseppini)
 ***************************************************************************************/
#pragma once

#include "BLabTypes.h"
#include "Colors.h"

#include <picojson.h>

#include <string>

struct StructuralMaterial
{
public:

    std::string Name;
    rgbaColor RenderColor;

    float Mass;
    float StaticFriction;
    float KineticFriction;

    SurfaceType Surface;

public:

    static StructuralMaterial Create(
        rgbaColor const & renderColor,
        picojson::object const & structuralMaterialJson);

    StructuralMaterial(
        std::string name,
        rgbaColor const & renderColor,
        float mass,
        float staticFriction,
        float kineticFriction,
        SurfaceType surface)
        : Name(name)
        , RenderColor(renderColor)
        , Mass(mass)
        , StaticFriction(staticFriction)
        , KineticFriction(kineticFriction)
        , Surface(surface)
    {}
};
