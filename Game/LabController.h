/***************************************************************************************
* Original Author:		Gabriele Giuseppini
* Created:				2020-05-16
* Copyright:			Gabriele Giuseppini  (https://github.com/GabrieleGiuseppini)
***************************************************************************************/
#pragma once

#include "GameEventDispatcher.h"
#include "GameParameters.h"
#include "MaterialDatabase.h"
#include "PerfStats.h"
#include "Physics.h"
#include "RenderContext.h"
#include "ShipDefinition.h"

#include <GameCore/GameChronometer.h>
#include <GameCore/GameTypes.h>
#include <GameCore/Vectors.h>

#include <cassert>
#include <filesystem>
#include <memory>
#include <optional>

/*
 * This class is responsible for managing the lab - both its lifetime and the user
 * interactions.
 */
class LabController
{
public:

    static std::unique_ptr<LabController> Create(
        int initialCanvasWidth,
        int initialCanvasHeight);

public:

    void RegisterBLabEventHandler(IBLabEventHandler * handler)
    {
        mGameEventHandler->RegisterBLabEventHandler(handler);
    }

    void RegisterNpcGameEventHandler(INpcGameEventHandler * handler)
    {
        mGameEventHandler->RegisterNpcGameEventHandler(handler);
    }

    //
    // Lab
    //

    SimulationControlStateType GetSimulationControlState() const;

    void SetSimulationControlState(SimulationControlStateType state);

    void SetSimulationControlPulse();

    void LoadShip(std::filesystem::path const & shipDefinitionFilepath);

    void Update();

    void Render();

    void Reset();

    //
    // Interactions
    //

    std::optional<ElementIndex> TryPickVertex(vec2f const & screenCoordinates) const;

    void MoveVertexBy(ElementIndex vertexIndex, vec2f const & screenOffset);

    void RotateShipBy(vec2f const & centerScreenCoordinates, float screenStride);

    void RotateShipBy(ElementIndex particleIndex, float screenStride);

    bool TrySelectOriginTriangle(vec2f const & screenCoordinates);

    bool TrySelectParticle(vec2f const & screenCoordinates);

    std::optional<ElementIndex> TryPickNpcParticle(vec2f const & screenCoordinates) const;

    void MoveNpcParticleBy(ElementIndex npcParticleIndex, vec2f const & screenOffset, vec2f const & inertialStride);

    void NotifyNpcParticleTrajectory(ElementIndex npcParticleIndex, vec2f const & targetScreenCoordinates);

    void SetNpcParticleTrajectory(ElementIndex npcParticleIndex, vec2f const & targetScreenCoordinates);

    void QueryNearestNpcParticleAt(vec2f const & screenCoordinates) const;

    vec2f const & GetShipVelocity() const;

    void SetShipVelocity(vec2f const & velocity);

    void SetWavesAmplitude(float wavesAmplitude);

    void SetWavesSpeed(float wavesSpeed);

    void QueryPointAt(vec2f const & screenCoordinates) const;

    void FlipCurrentlySelectedHuman();

    void DoStepForVideo();

    //
    // NPC
    //

    std::optional<PickedObjectId<NpcId>> BeginPlaceNewFurnitureNpc(FurnitureNpcKindType humanKind, vec2f const & screenCoordinates);
    std::optional<PickedObjectId<NpcId>> BeginPlaceNewHumanNpc(HumanNpcKindType humanKind, vec2f const & screenCoordinates);
    std::optional<PickedObjectId<NpcId>> ProbeNpcAt(vec2f const & screenCoordinates) const;
    void BeginMoveNpc(NpcId id);
    void MoveNpcTo(NpcId id, vec2f const & screenCoordinates, vec2f const & worldOffset);
    void EndMoveNpc(NpcId id);
    void CompleteNewNpc(NpcId id);
    void RemoveNpc(NpcId id);
    void AbortNewNpc(NpcId id);
    void HighlightNpc(NpcId id, NpcHighlightType highlight);

    void SetNpcPanicLevelForAllHumans(float panicLevel);

    //
    // Render controls
    //

    void SetCanvasSize(int width, int height)
    {
        assert(!!mRenderContext);
        mRenderContext->SetCanvasSize(width, height);
    }

    void Pan(vec2f const & screenOffset)
    {
        assert(!!mRenderContext);
        vec2f const worldOffset = mRenderContext->ScreenOffsetToWorldOffset(screenOffset);
        mRenderContext->SetCameraWorldPosition(mRenderContext->GetCameraWorldPosition() + worldOffset);
    }

    void ResetPan()
    {
        assert(!!mRenderContext);
        mRenderContext->SetCameraWorldPosition(vec2f::zero());
    }

    void AdjustZoom(float amount)
    {
        assert(!!mRenderContext);
        mRenderContext->SetZoom(mRenderContext->GetZoom() * amount);
    }

    void ResetZoom()
    {
        assert(!!mRenderContext);
        mRenderContext->SetZoom(1.0f);
    }

    vec2f ScreenToWorld(vec2f const & screenCoordinates) const
    {
        assert(!!mRenderContext);
        return mRenderContext->ScreenToWorld(screenCoordinates);
    }

    vec2f ScreenOffsetToWorldOffset(vec2f const & screenOffset) const
    {
        assert(!!mRenderContext);
        return mRenderContext->ScreenOffsetToWorldOffset(screenOffset);
    }

    float ScreenOffsetToWorldOffset(float screenOffset) const
    {
        assert(!!mRenderContext);
        return mRenderContext->ScreenOffsetToWorldOffset(screenOffset);
    }

    vec2f WorldToScreen(vec2f const & worldCoordinates) const
    {
        assert(!!mRenderContext);
        return mRenderContext->WorldToScreen(worldCoordinates);
    }

    bool IsViewGridEnabled() const
    {
        assert(!!mRenderContext);
        return mRenderContext->IsGridEnabled();
    }

    void SetViewGridEnabled(bool value)
    {
        assert(!!mRenderContext);
        mRenderContext->SetGridEnabled(value);
    }

    NpcRenderModeType GetNpcRenderMode() const
    {
        assert(!!mRenderContext);
        return mRenderContext->GetNpcRenderMode();
    }

    void SetNpcRenderMode(NpcRenderModeType value)
    {
        assert(!!mRenderContext);
        mRenderContext->SetNpcRenderMode(value);
    }

    //
    // Simmulation parameters
    //

    float GetGlobalDampingAdjustment() const { return mGameParameters.GlobalDampingAdjustment; }
    void SetGlobalDampingAdjustment(float value) { mGameParameters.GlobalDampingAdjustment = value; }
    float GetMinGlobalDampingAdjustment() const { return GameParameters::MinGlobalDampingAdjustment; }
    float GetMaxGlobalDampingAdjustment() const { return GameParameters::MaxGlobalDampingAdjustment; }

    float GetNpcMaterialElasticityAdjustment() const { return mGameParameters.NpcMaterialElasticityAdjustment; }
    void SetNpcMaterialElasticityAdjustment(float value) { mGameParameters.NpcMaterialElasticityAdjustment = value; }
    float GetMinNpcMaterialElasticityAdjustment() const { return GameParameters::MinNpcMaterialElasticityAdjustment; }
    float GetMaxNpcMaterialElasticityAdjustment() const { return GameParameters::MaxNpcMaterialElasticityAdjustment; }

    float GetNpcMaterialStaticFrictionAdjustment() const { return mGameParameters.NpcMaterialStaticFrictionAdjustment; }
    void SetNpcMaterialStaticFrictionAdjustment(float value) { mGameParameters.NpcMaterialStaticFrictionAdjustment = value; }
    float GetMinNpcMaterialStaticFrictionAdjustment() const { return GameParameters::MinNpcMaterialStaticFrictionAdjustment; }
    float GetMaxNpcMaterialStaticFrictionAdjustment() const { return GameParameters::MaxNpcMaterialStaticFrictionAdjustment; }

    float GetNpcMaterialKineticFrictionAdjustment() const { return mGameParameters.NpcMaterialKineticFrictionAdjustment; }
    void SetNpcMaterialKineticFrictionAdjustment(float value) { mGameParameters.NpcMaterialKineticFrictionAdjustment = value; }
    float GetMinNpcMaterialKineticFrictionAdjustment() const { return GameParameters::MinNpcMaterialKineticFrictionAdjustment; }
    float GetMaxNpcMaterialKineticFrictionAdjustment() const { return GameParameters::MaxNpcMaterialKineticFrictionAdjustment; }

    float GetGravityAdjustment() const { return mGameParameters.GravityAdjustment; }
    void SetGravityAdjustment(float value) { mGameParameters.GravityAdjustment = value; }
    float GetMinGravityAdjustment() const { return GameParameters::MinGravityAdjustment; }
    float GetMaxGravityAdjustment() const { return GameParameters::MaxGravityAdjustment; }

    float GetMassAdjustment() const { return mGameParameters.MassAdjustment; }
    void SetMassAdjustment(float value) { mGameParameters.MassAdjustment = value; }
    float GetMinMassAdjustment() const { return GameParameters::MinMassAdjustment; }
    float GetMaxMassAdjustment() const { return GameParameters::MaxMassAdjustment;  }

    float GetSeaLevel() const { return mOceanDepth; }
    void SetSeaLevel(float value) { mOceanDepth = value; if (mWorld) mWorld->GetOceanSurface().SetDepth(value); }
    float GetMinSeaLevel() const { return Physics::OceanSurface::MinDepth; }
    float GetMaxSeaLevel() const { return Physics::OceanSurface::MaxDepth; }

    float GetSpringReductionFractionAdjustment() const { return mGameParameters.NpcSpringReductionFractionAdjustment; }
    void SetSpringReductionFractionAdjustment(float value) { mGameParameters.NpcSpringReductionFractionAdjustment = value; }
    float GetMinSpringReductionFractionAdjustment() const { return GameParameters::MinNpcSpringReductionFractionAdjustment; }
    float GetMaxSpringReductionFractionAdjustment() const { return GameParameters::MaxNpcSpringReductionFractionAdjustment; }

    float GetSpringDampingCoefficientAdjustment() const { return mGameParameters.NpcSpringDampingCoefficientAdjustment; }
    void SetSpringDampingCoefficientAdjustment(float value) { mGameParameters.NpcSpringDampingCoefficientAdjustment = value; }
    float GetMinSpringDampingCoefficientAdjustment() const { return GameParameters::MinNpcSpringDampingCoefficientAdjustment; }
    float GetMaxSpringDampingCoefficientAdjustment() const { return GameParameters::MaxNpcSpringDampingCoefficientAdjustment; }

    float GetWaterFrictionDragAdjustment() const { return mGameParameters.WaterFrictionDragAdjustment; }
    void SetWaterFrictionDragAdjustment(float value) { mGameParameters.WaterFrictionDragAdjustment = value; }
    float GetMinWaterFrictionDragAdjustment() const { return GameParameters::MinWaterFrictionDragAdjustment; }
    float GetMaxWaterFrictionDragAdjustment() const { return GameParameters::MaxWaterFrictionDragAdjustment; }

    float GetBuoyancyAdjustment() const { return mGameParameters.BuoyancyAdjustment; }
    void SetBuoyancyAdjustment(float value) { mGameParameters.BuoyancyAdjustment = value; }
    float GetMinBuoyancyAdjustment() const { return GameParameters::MinBuoyancyAdjustment; }
    float GetMaxBuoyancyAdjustment() const { return GameParameters::MaxBuoyancyAdjustment; }

    float GetNpcSizeAdjustment() const { return mGameParameters.NpcSizeAdjustment; }
    void SetNpcSizeAdjustment(float value) { mGameParameters.NpcSizeAdjustment = value; }
    float GetMinNpcSizeAdjustment() const { return GameParameters::MinNpcSizeAdjustment; }
    float GetMaxNpcSizeAdjustment() const { return GameParameters::MaxNpcSizeAdjustment; }

    float GetHumanNpcEquilibriumTorqueStiffnessCoefficient() const { return mGameParameters.HumanNpcEquilibriumTorqueStiffnessCoefficient; }
    void SetHumanNpcEquilibriumTorqueStiffnessCoefficient(float value) { mGameParameters.HumanNpcEquilibriumTorqueStiffnessCoefficient = value; }
    float GetMinHumanNpcEquilibriumTorqueStiffnessCoefficient() const { return GameParameters::MinHumanNpcEquilibriumTorqueStiffnessCoefficient; }
    float GetMaxHumanNpcEquilibriumTorqueStiffnessCoefficient() const { return GameParameters::MaxHumanNpcEquilibriumTorqueStiffnessCoefficient; }

    float GetHumanNpcEquilibriumTorqueDampingCoefficient() const { return mGameParameters.HumanNpcEquilibriumTorqueDampingCoefficient; }
    void SetHumanNpcEquilibriumTorqueDampingCoefficient(float value) { mGameParameters.HumanNpcEquilibriumTorqueDampingCoefficient = value; }
    float GetMinHumanNpcEquilibriumTorqueDampingCoefficient() const { return GameParameters::MinHumanNpcEquilibriumTorqueDampingCoefficient; }
    float GetMaxHumanNpcEquilibriumTorqueDampingCoefficient() const { return GameParameters::MaxHumanNpcEquilibriumTorqueDampingCoefficient; }

    float GetHumanNpcWalkingSpeedAdjustment() const { return mGameParameters.HumanNpcWalkingSpeedAdjustment; }
    void SetHumanNpcWalkingSpeedAdjustment(float value) { mGameParameters.HumanNpcWalkingSpeedAdjustment = value; }
    float GetMinHumanNpcWalkingSpeedAdjustment() const { return GameParameters::MinHumanNpcWalkingSpeedAdjustment; }
    float GetMaxHumanNpcWalkingSpeedAdjustment() const { return GameParameters::MaxHumanNpcWalkingSpeedAdjustment; }

private:

    LabController(
        MaterialDatabase && materialDatabase,
        std::unique_ptr<Render::RenderContext> renderContext);

    void Reset(
        std::unique_ptr<Physics::Ship> ship,
        std::filesystem::path const & shipDefinitionFilepath);

    void UpdateShipTransformations();

private:

    MaterialDatabase mMaterialDatabase;
    std::unique_ptr<Render::RenderContext> mRenderContext;
    std::shared_ptr<GameEventDispatcher> mGameEventHandler;

    //
    // Current state
    //

    GameParameters mGameParameters;
    std::unique_ptr<Physics::World> mWorld;
    std::optional<std::filesystem::path> mCurrentShipFilePath;

    float mCurrentSimulationTime;

    PerfStats mPerfStats;
    PerfStats mLastPublishedPerfStats;
    GameChronometer::time_point mLastPerfPublishTimestamp;

    //
    // Simulation parameters owned by us
    //

    float mOceanDepth;

    vec2f mCurrentShipTranslationVelocity;
    float mCurrentShipTranslationAccelerationIndicator;

    float mCurrentWavesAmplitude;
    float mTargetWavesAmplitude;
    float mCurrentWavesSpeed;
    float mTargetWavesSpeed;
    float mLastWaveRotationAngle;
    float mLastWaveTimeArg;

    int mCurrentVideoStep;

    //
    // Simulation control
    //

    SimulationControlStateType mSimulationControlState;
    bool mSimulationControlImpulse;
};
