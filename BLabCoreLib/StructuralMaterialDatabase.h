/***************************************************************************************
* Original Author:		Gabriele Giuseppini
* Created:				2020-05-16
* Copyright:			Gabriele Giuseppini  (https://github.com/GabrieleGiuseppini)
***************************************************************************************/
#pragma once

#include "ResourceLocator.h"
#include "StructuralMaterial.h"

#include "BLabException.h"
#include "Colors.h"
#include "Utils.h"

#include <picojson.h>

#include <cassert>
#include <cstdint>
#include <map>

class StructuralMaterialDatabase
{
public:

    using ColorKey = rgbColor;

    enum class UniqueMaterialKeyType
    {
        Furniture,
        HumanHead,
        HumanFeet
    };

public:

    static StructuralMaterialDatabase Load()
    {
        picojson::value structuralMaterialsRoot = Utils::ParseJSONFile(ResourceLocator::GetStructuralMaterialDatabaseFilePath());

        if (!structuralMaterialsRoot.is<picojson::array>())
        {
            throw BLabException("Structural materials definition is not a JSON array");
        }

        //
        // Read into map
        //

        std::map<ColorKey, StructuralMaterial> structuralMaterialsMap;
        std::map<UniqueMaterialKeyType, ColorKey> structuralMaterialByUniqueKeyMap;

        picojson::array const & structuralMaterialsRootArray = structuralMaterialsRoot.get<picojson::array>();
        for (auto const & materialElem : structuralMaterialsRootArray)
        {
            if (!materialElem.is<picojson::object>())
            {
                throw BLabException("Found a non-object in structural materials definition");
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
                throw BLabException("Structural material \"" + material.Name + "\" has a duplicate color key");
            }

            // Store
            auto const storedEntry = structuralMaterialsMap.emplace(
                std::make_pair(
                    colorKey,
                    material));

            assert(storedEntry.second);

            // Store by unique key, if any
            auto const uniqueKeyStr = Utils::GetOptionalJsonMember<std::string>(materialObject, "unique_material");
            if (uniqueKeyStr.has_value())
            {
                UniqueMaterialKeyType const uniqueMaterialKey = StrToUniqueMaterialKey(*uniqueKeyStr);
                auto [_, isUniqueMaterialInserted] = structuralMaterialByUniqueKeyMap.emplace(
                    std::make_pair(
                        uniqueMaterialKey,
                        colorKey));
                if (!isUniqueMaterialInserted)
                {
                    throw BLabException("Found more than one material with unique key \"" + *uniqueKeyStr + "\"");
                }
            }
        }

        return StructuralMaterialDatabase(
            std::move(structuralMaterialsMap),
            std::move(structuralMaterialByUniqueKeyMap));
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

    StructuralMaterial const & GetStructuralMaterial(UniqueMaterialKeyType const & uniqueKey) const
    {
        return mStructuralMaterialMap.at(mStructuralMaterialByUniqueKeyMap.at(uniqueKey));
    }

private:

    StructuralMaterialDatabase(
        std::map<ColorKey, StructuralMaterial> && structuralMaterialMap,
        std::map<UniqueMaterialKeyType, ColorKey> && structuralMaterialByUniqueKeyMap)
        : mStructuralMaterialMap(std::move(structuralMaterialMap))
        , mStructuralMaterialByUniqueKeyMap(std::move(structuralMaterialByUniqueKeyMap))
    {
    }

    static UniqueMaterialKeyType StrToUniqueMaterialKey(std::string const & str)
    {
        if (Utils::CaseInsensitiveEquals(str, "Furniture"))
            return UniqueMaterialKeyType::Furniture;
        else if (Utils::CaseInsensitiveEquals(str, "HumanFeet"))
            return UniqueMaterialKeyType::HumanFeet;
        else if (Utils::CaseInsensitiveEquals(str, "HumanHead"))
            return UniqueMaterialKeyType::HumanHead;
        else
            throw BLabException("Unrecognized unique material key \"" + str + "\"");
    }

    std::map<ColorKey, StructuralMaterial> const mStructuralMaterialMap;
    std::map<UniqueMaterialKeyType, ColorKey> const mStructuralMaterialByUniqueKeyMap;
};
