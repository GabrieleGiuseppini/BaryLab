/***************************************************************************************
* Original Author:		Gabriele Giuseppini
* Created:				2020-05-23
* Copyright:			Gabriele Giuseppini  (https://github.com/GabrieleGiuseppini)
***************************************************************************************/
#pragma once

#include <BLabCoreLib/LabController.h>
#include <BLabCoreLib/Settings.h>

enum class SLabSettings : size_t
{
    Elasticity = 0,
    StaticFriction,
    KineticFriction,
    MassAdjustment,
    GravityAdjustment,
    SpringReductionFraction,
    SpringDampingCoefficient,
    HumanNpcRisingTorqueFactor,

    _Last = HumanNpcRisingTorqueFactor
};

class SettingsManager final : public BaseSettingsManager<SLabSettings>
{
public:

    SettingsManager(
        std::shared_ptr<LabController> labController,
        std::filesystem::path const & rootUserSettingsDirectoryPath);

private:

    static BaseSettingsManagerFactory MakeSettingsFactory(
        std::shared_ptr<LabController> labController);
};