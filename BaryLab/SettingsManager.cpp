/***************************************************************************************
* Original Author:		Gabriele Giuseppini
* Created:				2020-05-23
* Copyright:			Gabriele Giuseppini  (https://github.com/GabrieleGiuseppini)
***************************************************************************************/
#include "SettingsManager.h"

#include <cctype>
#include <sstream>

std::string MangleSettingName(std::string && settingName);

#define ADD_SETTING(type, name)                         \
    factory.AddSetting<type>(                           \
        SLabSettings::name,                             \
        MangleSettingName(#name),                       \
        [labController]() -> type { return labController->Get##name(); }, \
        [labController](auto const & v) { labController->Set##name(v); }, \
		[labController](auto const & v) { labController->Set##name(v); });

BaseSettingsManager<SLabSettings>::BaseSettingsManagerFactory SettingsManager::MakeSettingsFactory(
    std::shared_ptr<LabController> labController)
{
    BaseSettingsManagerFactory factory;

    ADD_SETTING(float, ElasticityAdjustment);
    ADD_SETTING(float, StaticFrictionAdjustment);
    ADD_SETTING(float, KineticFrictionAdjustment);
    ADD_SETTING(float, MassAdjustment);
    ADD_SETTING(float, GravityAdjustment);
    ADD_SETTING(float, GlobalDampingAdjustment);
    ADD_SETTING(float, SeaLevel);
    ADD_SETTING(float, SpringReductionFractionAdjustment);
    ADD_SETTING(float, SpringDampingCoefficientAdjustment);
    ADD_SETTING(float, WaterFrictionDragAdjustment);
    ADD_SETTING(float, BuoyancyAdjustment);
    ADD_SETTING(float, NpcSizeMultiplier);
    ADD_SETTING(float, HumanNpcEquilibriumTorqueStiffnessCoefficient);
    ADD_SETTING(float, HumanNpcEquilibriumTorqueDampingCoefficient);
    ADD_SETTING(float, HumanNpcWalkingSpeedAdjustment);

    return factory;
}

SettingsManager::SettingsManager(
    std::shared_ptr<LabController> labController,
    std::filesystem::path const & rootUserSettingsDirectoryPath)
    : BaseSettingsManager<SLabSettings>(
        MakeSettingsFactory(labController),
        rootUserSettingsDirectoryPath,
        rootUserSettingsDirectoryPath)
{}

std::string MangleSettingName(std::string && settingName)
{
    std::stringstream ss;

    bool isFirst = true;
    for (char ch : settingName)
    {
        if (std::isupper(ch))
        {
            if (!isFirst)
            {
                ss << '_';
            }

            ss << static_cast<char>(std::tolower(ch));
        }
        else
        {
            ss << ch;
        }

        isFirst = false;
    }

    return ss.str();
}