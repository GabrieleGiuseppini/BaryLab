/***************************************************************************************
* Original Author:		Gabriele Giuseppini
* Created:				2020-05-23
* Copyright:			Gabriele Giuseppini  (https://github.com/GabrieleGiuseppini)
***************************************************************************************/
#pragma once

#include <Game/LabController.h>

#include <GameCore/Settings.h>

enum class SLabSettings : size_t
{
    ElasticityAdjustment = 0,
    StaticFrictionAdjustment,
    KineticFrictionAdjustment,
    MassAdjustment,
    GravityAdjustment,
    GlobalDampingAdjustment,
    SeaLevel,
    SpringReductionFractionAdjustment,
    SpringDampingCoefficientAdjustment,
    WaterFrictionDragAdjustment,
    BuoyancyAdjustment,
    NpcSizeMultiplier,
    HumanNpcEquilibriumTorqueStiffnessCoefficient,
    HumanNpcEquilibriumTorqueDampingCoefficient,
    HumanNpcWalkingSpeedAdjustment,

    _Last = HumanNpcWalkingSpeedAdjustment
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