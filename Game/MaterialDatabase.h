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

    template<typename TMaterial>
    using MaterialColorMap = std::map<MaterialColorKey, TMaterial>;

    template<typename TMaterial>
    using MaterialNameMap = std::map<std::string, MaterialColorKey>;

public:

    static MaterialDatabase Load()
    {
        //
        // Structural
        //

        picojson::value structuralMaterialsRoot = Utils::ParseJSONFile(ResourceLocator::GetStructuralMaterialDatabaseFilePath());

        if (!structuralMaterialsRoot.is<picojson::object>())
        {
            throw GameException("Structural materials definition is not a JSON object");
        }

        picojson::object const & structuralMaterialsRootObj = structuralMaterialsRoot.get<picojson::object>();

        // Read into map

        MaterialColorMap<StructuralMaterial> structuralMaterialColorMap;
        MaterialNameMap<StructuralMaterial> structuralMaterialNameMap;

        picojson::array const & structuralMaterialsRootArray = Utils::GetMandatoryJsonArray(structuralMaterialsRootObj, "materials");
        for (auto const & materialElem : structuralMaterialsRootArray)
        {
            if (!materialElem.is<picojson::object>())
            {
                throw GameException("Found a non-object in structural materials definition");
            }

            picojson::object const & materialObject = materialElem.get<picojson::object>();

            // Parse color key

            auto const & memberIt = materialObject.find("color_key");
            if (materialObject.end() == memberIt)
            {
                throw GameException("Error parsing JSON: cannot find member \"color_key\"");
            }

            MaterialColorKey colorKey;
            if (memberIt->second.is<std::string>())
            {
                colorKey = Utils::Hex2RgbColor(memberIt->second.get<std::string>());
            }
            else if (memberIt->second.is<picojson::array>())
            {
                auto const & memberArray = memberIt->second.get<picojson::array>();

                // Take first
                assert(memberArray.size() > 0);
                if (!memberArray[0].is<std::string>())
                {
                    throw GameException("Error parsing JSON: an element of the material's \"color_key\" array is not a 'string'");
                }

                // Take first
                colorKey = Utils::Hex2RgbColor(memberArray[0].get<std::string>());
            }
            else
            {
                throw GameException("Error parsing JSON: material's \"color_key\" member is neither a 'string' nor an 'array'");
            }

            // Parse material

            StructuralMaterial material = StructuralMaterial::Create(
                colorKey,
                materialObject);

            // Make sure there are no dupes
            if (structuralMaterialColorMap.count(colorKey) != 0)
            {
                throw GameException("Structural material \"" + material.Name + "\" has a duplicate color key");
            }

            // Store
            auto const storedEntry = structuralMaterialColorMap.emplace(
                std::make_pair(
                    colorKey,
                    material));

            assert(storedEntry.second);

            // Store by name - making sure there are no dupes
            auto const [_, isNameInserted] = structuralMaterialNameMap.emplace(
                std::make_pair(
                    material.Name,
                    colorKey));
            if (!isNameInserted)
            {
                throw GameException("Material name \"" + material.Name + "\" already belongs to another material");
            }
        }

        return MaterialDatabase(
            std::move(structuralMaterialColorMap),
            std::move(structuralMaterialNameMap));
    }

    StructuralMaterial const * FindStructuralMaterial(MaterialColorKey const & colorKey) const
    {
        if (auto srchIt = mStructuralMaterialColorMap.find(colorKey);
            srchIt != mStructuralMaterialColorMap.end())
        {
            // Found color key verbatim!
            return &(srchIt->second);
        }

        // No luck
        return nullptr;
    }

    StructuralMaterial const & GetStructuralMaterial(std::string const & name) const
    {
        if (auto const srchIt = mStructuralMaterialNameMap.find(name);
            srchIt != mStructuralMaterialNameMap.cend())
        {
            // Found
            return mStructuralMaterialColorMap.at(srchIt->second);
        }

        throw GameException("Cannot find material \"" + name + "\"");
    }

private:

    MaterialDatabase(
        MaterialColorMap<StructuralMaterial> && structuralMaterialColorMap,
        MaterialNameMap<StructuralMaterial> && structuralMaterialNameMap)
        : mStructuralMaterialColorMap(std::move(structuralMaterialColorMap))
        , mStructuralMaterialNameMap(std::move(structuralMaterialNameMap))
    {
    }

    // Structural
    MaterialColorMap<StructuralMaterial> mStructuralMaterialColorMap;
    MaterialNameMap<StructuralMaterial> mStructuralMaterialNameMap;
};
