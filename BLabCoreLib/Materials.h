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

    SurfaceType Surface;

public:

    static StructuralMaterial Create(
        rgbaColor const & renderColor,
        picojson::object const & structuralMaterialJson);

    StructuralMaterial(
        std::string name,
        rgbaColor const & renderColor,
        SurfaceType surface)
        : Name(name)
        , RenderColor(renderColor)
        , Surface(surface)
    {}
};

struct NpcMaterial
{
public:

    enum class KindType
    {
        Furniture,
        HumanHead,
        HumanFeet
    };

public:

    std::string Name;
    KindType Kind;
    rgbaColor RenderColor;

    float Mass;
    float StaticFriction;
    float KineticFriction;
    float Elasticity;
    float BuoyancyVolumeFill;

    SurfaceType Surface;

public:

    static NpcMaterial Create(picojson::object const & npcMaterialJson);

    NpcMaterial(
        std::string name,
        KindType kind,
        rgbaColor const & renderColor,
        float mass,
        float staticFriction,
        float kineticFriction,
        float elasticity,
        float buoyancyVolumeFill,
        SurfaceType surface)
        : Name(name)
        , Kind(kind)
        , RenderColor(renderColor)
        , Mass(mass)
        , StaticFriction(staticFriction)
        , KineticFriction(kineticFriction)
        , Elasticity(elasticity)
        , BuoyancyVolumeFill(buoyancyVolumeFill)
        , Surface(surface)
    {}

private:

    static KindType StrToKindType(std::string const & strKind);
};
