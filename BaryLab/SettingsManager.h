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
    NpcMaterialElasticityAdjustment = 0,
    NpcMaterialStaticFrictionAdjustment,
    NpcMaterialKineticFrictionAdjustment,
    MassAdjustment,
    GravityAdjustment,
    NpcDamping,
    SeaLevel,
    SpringReductionFractionAdjustment,
    SpringDampingCoefficientAdjustment,
    WaterFrictionDragAdjustment,
    BuoyancyAdjustment,
    HumanNpcEquilibriumTorqueStiffnessCoefficient,
    HumanNpcEquilibriumTorqueDampingCoefficient,
    HumanNpcWalkingSpeedAdjustment,
    HumanNpcBodyLengthAdjustment,

    _Last = HumanNpcBodyLengthAdjustment
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