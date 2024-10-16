/***************************************************************************************
* Original Author:		Gabriele Giuseppini
* Created:				2024-07-13
* Copyright:			Gabriele Giuseppini  (https://github.com/GabrieleGiuseppini)
***************************************************************************************/
#pragma once

#include "MaterialDatabase.h"
#include "RenderTypes.h"

#include <map>
#include <string>
#include <vector>

/*
 * Information over the different sub-kinds of NPCs.
 */
class NpcDatabase
{
public:

    enum class ParticleMeshKindType
    {
        Particle,
        Dipole,
        Quad
    };

    struct ParticleAttributesType
    {
        float BuoyancyVolumeFill;
        float SpringReductionFraction;
        float SpringDampingCoefficient;
    };

    struct HumanTextureGeometryType
    {
        float HeadLengthFraction; // Fraction of dipole length; can overshoot 1.0 - leg+torso for e.g. hats
        float HeadWHRatio; // To recover head quad width from height
        float TorsoLengthFraction; // Fraction of dipole length
        float TorsoWHRatio; // To recover torso quad width from height
        float ArmLengthFraction; // Fraction of dipole length
        float ArmWHRatio; // To recover arm quad width from height
        float LegLengthFraction; // Fraction of dipole length
        float LegWHRatio; // To recover leg quad width from height
    };

    struct FurnitureGeometryType
    {
        float Width;
        float Height;
    };

    struct HumanTextureFramesType
    {
        Render::TextureCoordinatesQuad HeadFront;
        Render::TextureCoordinatesQuad HeadBack;
        Render::TextureCoordinatesQuad HeadSide;

        Render::TextureCoordinatesQuad TorsoFront;
        Render::TextureCoordinatesQuad TorsoBack;
        Render::TextureCoordinatesQuad TorsoSide;

        Render::TextureCoordinatesQuad ArmFront;
        Render::TextureCoordinatesQuad ArmBack;
        Render::TextureCoordinatesQuad ArmSide;

        Render::TextureCoordinatesQuad LegFront;
        Render::TextureCoordinatesQuad LegBack;
        Render::TextureCoordinatesQuad LegSide;
    };

public:

    // Only movable
    NpcDatabase(NpcDatabase const & other) = delete;
    NpcDatabase(NpcDatabase && other) = default;
    NpcDatabase & operator=(NpcDatabase const & other) = delete;
    NpcDatabase & operator=(NpcDatabase && other) = default;

    static NpcDatabase Load(MaterialDatabase const & materialDatabase);

    // Humans

    std::vector<std::vector<NpcSubKindIdType>> const & GetHumanSubKindIdsByRole() const
    {
        return mHumanSubKindIdsByRole;
    }

    NpcHumanRoleType GetHumanRole(NpcSubKindIdType subKindId) const
    {
        return mHumanSubKinds.at(subKindId).Role;
    }

    rgbColor GetHumanRenderColor(NpcSubKindIdType subKindId) const
    {
        return mHumanSubKinds.at(subKindId).RenderColor;
    }

    StructuralMaterial const & GetHumanHeadMaterial(NpcSubKindIdType subKindId) const
    {
        return mHumanSubKinds.at(subKindId).HeadMaterial;
    }

    StructuralMaterial const & GetHumanFeetMaterial(NpcSubKindIdType subKindId) const
    {
        return mHumanSubKinds.at(subKindId).FeetMaterial;
    }

    ParticleAttributesType const & GetHumanHeadParticleAttributes(NpcSubKindIdType subKindId) const
    {
        assert(mHumanSubKinds.at(subKindId).ParticleAttributes.size() == 2);
        return mHumanSubKinds.at(subKindId).ParticleAttributes[1];
    }

    ParticleAttributesType const & GetHumanFeetParticleAttributes(NpcSubKindIdType subKindId) const
    {
        assert(mHumanSubKinds.at(subKindId).ParticleAttributes.size() == 2);
        return mHumanSubKinds.at(subKindId).ParticleAttributes[0];
    }

    float GetHumanSizeMultiplier(NpcSubKindIdType subKindId) const
    {
        return mHumanSubKinds.at(subKindId).SizeMultiplier;
    }

    HumanTextureGeometryType const & GetHumanTextureGeometry(NpcSubKindIdType subKindId) const
    {
        return mHumanSubKinds.at(subKindId).TextureGeometry;
    }

    float GetHumanBodyWidthRandomizationSensitivity(NpcSubKindIdType subKindId) const
    {
        return mHumanSubKinds.at(subKindId).BodyWidthRandomizationSensitivity;
    }

    HumanTextureFramesType const & GetHumanTextureCoordinatesQuads(NpcSubKindIdType subKindId) const
    {
        return mHumanSubKinds.at(subKindId).TextureCoordinatesQuads;
    }

    // Furniture

    std::vector<std::vector<NpcSubKindIdType>> const & GetFurnitureSubKindIdsByRole() const
    {
        return mFurnitureSubKindIdsByRole;
    }

    NpcFurnitureRoleType GetFurnitureRole(NpcSubKindIdType subKindId) const
    {
        return mFurnitureSubKinds.at(subKindId).Role;
    }

    rgbColor GetFurnitureRenderColor(NpcSubKindIdType subKindId) const
    {
        return mFurnitureSubKinds.at(subKindId).RenderColor;
    }

    StructuralMaterial const & GetFurnitureMaterial(NpcSubKindIdType subKindId) const
    {
        return mFurnitureSubKinds.at(subKindId).Material;
    }

    ParticleAttributesType const & GetFurnitureParticleAttributes(NpcSubKindIdType subKindId, int particleOrdinal) const
    {
        assert(particleOrdinal < mFurnitureSubKinds.at(subKindId).ParticleAttributes.size());
        return mFurnitureSubKinds.at(subKindId).ParticleAttributes[particleOrdinal];
    }

    ParticleMeshKindType const & GetFurnitureParticleMeshKindType(NpcSubKindIdType subKindId) const
    {
        return mFurnitureSubKinds.at(subKindId).ParticleMeshKind;
    }

    FurnitureGeometryType const & GetFurnitureGeometry(NpcSubKindIdType subKindId) const
    {
        return mFurnitureSubKinds.at(subKindId).Geometry;
    }

    Render::TextureCoordinatesQuad const & GetFurnitureTextureCoordinatesQuad(NpcSubKindIdType subKindId) const
    {
        return mFurnitureSubKinds.at(subKindId).TextureCoordinatesQuad;
    }

private:

    struct HumanSubKind
    {
        std::string Name;
        NpcHumanRoleType Role;
        rgbColor RenderColor;

        StructuralMaterial const & HeadMaterial;
        StructuralMaterial const & FeetMaterial;

        std::array<ParticleAttributesType, 2> ParticleAttributes;

        float SizeMultiplier;
        float BodyWidthRandomizationSensitivity;

        HumanTextureFramesType const TextureCoordinatesQuads;
        HumanTextureGeometryType TextureGeometry;
    };

    struct FurnitureSubKind
    {
        std::string Name;
        NpcFurnitureRoleType Role;
        rgbColor RenderColor;

        StructuralMaterial const & Material;

        std::vector<ParticleAttributesType> ParticleAttributes;

        ParticleMeshKindType ParticleMeshKind;

        FurnitureGeometryType Geometry;

        Render::TextureCoordinatesQuad TextureCoordinatesQuad;
    };

private:

    NpcDatabase(
        std::map<NpcSubKindIdType, HumanSubKind> && humanSubKinds,
        std::map<NpcSubKindIdType, FurnitureSubKind> && furnitureSubKinds);

    static HumanSubKind ParseHumanSubKind(
        picojson::object const & subKindObject,
        StructuralMaterial const & headMaterial,
        StructuralMaterial const & feetMaterial,
        ParticleAttributesType const & globalHeadParticleAttributes,
        ParticleAttributesType const & globalFeetParticleAttributes);

    static FurnitureSubKind ParseFurnitureSubKind(
        picojson::object const & subKindObject,
        MaterialDatabase const & materialDatabase);

    static ParticleAttributesType MakeParticleAttributes(
        picojson::object const & containerObject,
        std::string const & particleAttributesOverrideMemberName,
        ParticleAttributesType const & defaultParticleAttributes);

    static ParticleAttributesType MakeParticleAttributes(
        picojson::object const & particleAttributesOverrideJsonObject,
        ParticleAttributesType const & defaultParticleAttributes);

    static ParticleAttributesType MakeDefaultParticleAttributes(StructuralMaterial const & baseMaterial);

    template<typename TNpcSubKindContainer>
    static std::vector<std::tuple<NpcSubKindIdType, std::string>> GetSubKinds(
        TNpcSubKindContainer const & container,
        std::string const & language);

    static ParticleMeshKindType StrToParticleMeshKindType(std::string const & str);

private:

    std::map<NpcSubKindIdType, HumanSubKind> mHumanSubKinds;
    std::map<NpcSubKindIdType, FurnitureSubKind> mFurnitureSubKinds;

    std::vector<std::vector<NpcSubKindIdType>> mHumanSubKindIdsByRole;
    std::vector<std::vector<NpcSubKindIdType>> mFurnitureSubKindIdsByRole;
};
