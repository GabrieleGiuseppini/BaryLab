/***************************************************************************************
* Original Author:		Gabriele Giuseppini
* Created:				2020-05-16
* Copyright:			Gabriele Giuseppini  (https://github.com/GabrieleGiuseppini)
***************************************************************************************/
#pragma once

#include "Materials.h"
#include "ResourceLocator.h"

#include <GameCore/Colors.h>
#include <GameCore/GameException.h>
#include <GameCore/Utils.h>

#include <picojson.h>

#include <cassert>
#include <cstdint>
#include <map>

class MaterialDatabase
{
public:

    using ColorKey = rgbColor;

public:

    static MaterialDatabase Load()
    {
        //
        // Structural
        //

        picojson::value structuralMaterialsRoot = Utils::ParseJSONFile(ResourceLocator::GetStructuralMaterialDatabaseFilePath());

        if (!structuralMaterialsRoot.is<picojson::array>())
        {
            throw GameException("Structural materials definition is not a JSON array");
        }

        // Read into map

        std::map<ColorKey, StructuralMaterial> structuralMaterialsMap;

        picojson::array const & structuralMaterialsRootArray = structuralMaterialsRoot.get<picojson::array>();
        for (auto const & materialElem : structuralMaterialsRootArray)
        {
            if (!materialElem.is<picojson::object>())
            {
                throw GameException("Found a non-object in structural materials definition");
            }

            picojson::object const & materialObject = materialElem.get<picojson::object>();

            ColorKey colorKey = Utils::Hex2RgbColor(
                Utils::GetMandatoryJsonMember<std::string>(materialObject, "color_key"));

            StructuralMaterial material = StructuralMaterial::Create(
                rgbaColor(colorKey, rgbaColor::data_type_max),
                materialObject);

            // Make sure there are no dupes
            if (structuralMaterialsMap.count(colorKey) != 0)
            {
                throw GameException("Structural material \"" + material.Name + "\" has a duplicate color key");
            }

            // Store
            auto const storedEntry = structuralMaterialsMap.emplace(
                std::make_pair(
                    colorKey,
                    material));

            assert(storedEntry.second);
        }

        //
        // NPC
        //

        picojson::value npcMaterialsRoot = Utils::ParseJSONFile(ResourceLocator::GetNpcMaterialDatabaseFilePath());

        if (!npcMaterialsRoot.is<picojson::array>())
        {
            throw GameException("NPC materials definition is not a JSON array");
        }

        // Read into map

        std::map<NpcMaterial::KindType, NpcMaterial> npcMaterialsMap;

        picojson::array const & npcMaterialsRootArray = npcMaterialsRoot.get<picojson::array>();
        for (auto const & materialElem : npcMaterialsRootArray)
        {
            if (!materialElem.is<picojson::object>())
            {
                throw GameException("Found a non-object in NPC materials definition");
            }

            picojson::object const & materialObject = materialElem.get<picojson::object>();

            NpcMaterial material = NpcMaterial::Create(materialObject);

            // Store
            auto const storedEntry = npcMaterialsMap.emplace(
                material.Kind,
                material);

            if (!storedEntry.second)
            {
                throw GameException("NPC material \"" + material.Name + "\" has a duplicate kind");
            }
        }

        return MaterialDatabase(
            std::move(structuralMaterialsMap),
            std::move(npcMaterialsMap));
    }

    StructuralMaterial const * FindStructuralMaterial(ColorKey const & colorKey) const
    {
        if (auto srchIt = mStructuralMaterialMap.find(colorKey);
            srchIt != mStructuralMaterialMap.end())
        {
            // Found color key verbatim!
            return &(srchIt->second);
        }

        // No luck
        return nullptr;
    }

    NpcMaterial const & GetNpcMaterial(NpcMaterial::KindType kind) const
    {
        return mNpcMaterialMap.at(kind);
    }

private:

    MaterialDatabase(
        std::map<ColorKey, StructuralMaterial> && structuralMaterialMap,
        std::map<NpcMaterial::KindType, NpcMaterial> && npcMaterialMap)
        : mStructuralMaterialMap(std::move(structuralMaterialMap))
        , mNpcMaterialMap(std::move(npcMaterialMap))
    {
    }

    // Structural
    std::map<ColorKey, StructuralMaterial> const mStructuralMaterialMap;

    // NPC
    std::map<NpcMaterial::KindType, NpcMaterial> const mNpcMaterialMap;
};
