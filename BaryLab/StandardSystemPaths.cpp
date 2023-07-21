/***************************************************************************************
* Original Author:		Gabriele Giuseppini
* Created:				2020-05-15
* Copyright:			Gabriele Giuseppini  (https://github.com/GabrieleGiuseppini)
***************************************************************************************/
#include "StandardSystemPaths.h"

#include <BLabCoreLib/Version.h>

#include <wx/stdpaths.h>

StandardSystemPaths * StandardSystemPaths::mSingleInstance = nullptr;

std::filesystem::path StandardSystemPaths::GetUserSettingsRootFolderPath() const
{
    auto userFolder = wxStandardPaths::Get().GetUserConfigDir();

    return std::filesystem::path(userFolder.ToStdString())
        / ApplicationName // Without version - we want this to be sticky across upgrades
        / "Settings";
}