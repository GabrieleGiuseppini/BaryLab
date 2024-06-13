/***************************************************************************************
 * Original Author:		Gabriele Giuseppini
 * Created:				2020-05-16
 * Copyright:			Gabriele Giuseppini  (https://github.com/GabrieleGiuseppini)
 ***************************************************************************************/
#pragma once

#include <GameCore/Colors.h>
#include <GameCore/GameTypes.h>

#include <picojson.h>

#include <string>

struct StructuralMaterial
{
public:

    std::string Name;
    rgbaColor RenderColor;

    bool IsHull;

public:

    static StructuralMaterial Create(
        rgbaColor const & renderColor,
        picojson::object const & structuralMaterialJson);

    StructuralMaterial(
        std::string name,
        rgbaColor const & renderColor,
        bool isHull)
        : Name(name)
        , RenderColor(renderColor)
        , IsHull(isHull)
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
    float SpringReductionFraction;
    float SpringDampingCoefficient;
    float StaticFriction;
    float KineticFriction;
    float Elasticity;
    float BuoyancyVolumeFill;

public:

    static NpcMaterial Create(picojson::object const & npcMaterialJson);

    NpcMaterial(
        std::string name,
        KindType kind,
        rgbaColor const & renderColor,
        float mass,
        float springReductionFraction,
        float springDampingCoefficient,
        float staticFriction,
        float kineticFriction,
        float elasticity,
        float buoyancyVolumeFill)
        : Name(name)
        , Kind(kind)
        , RenderColor(renderColor)
        , Mass(mass)
        , SpringReductionFraction(springReductionFraction)
        , SpringDampingCoefficient(springDampingCoefficient)
        , StaticFriction(staticFriction)
        , KineticFriction(kineticFriction)
        , Elasticity(elasticity)
        , BuoyancyVolumeFill(buoyancyVolumeFill)
    {}

private:

    static KindType StrToKindType(std::string const & strKind);
};
