/***************************************************************************************
* Original Author:		Gabriele Giuseppini
* Created:				2020-05-16
* Copyright:			Gabriele Giuseppini  (https://github.com/GabrieleGiuseppini)
***************************************************************************************/
#pragma once

#include <filesystem>

class ResourceLocator
{
public:

    static std::filesystem::path GetInstalledMeshesFolderPath()
    {
        return std::filesystem::absolute(std::filesystem::path("Meshes"));
    }

    static std::filesystem::path GetDefaultMeshDefinitionFilePath()
    {
        return GetInstalledMeshesFolderPath() / "default_mesh.png";
    }

    static std::filesystem::path GetStructuralMaterialDatabaseFilePath()
    {
        return std::filesystem::absolute(std::filesystem::path("Data") / "materials_structural.json");
    }
    
    static std::filesystem::path GetResourcesFolderPath()
    {
        return std::filesystem::absolute(std::filesystem::path("Data") / "Resources");
    }

    static std::filesystem::path GetShadersRootFolderPath()
    {
        return std::filesystem::absolute(std::filesystem::path("Data") / "Shaders");
    }
};
