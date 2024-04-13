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

    float GetElasticityAdjustment() const { return mGameParameters.ElasticityAdjustment; }
    void SetElasticityAdjustment(float value) { mGameParameters.ElasticityAdjustment = value; }
    float GetMinElasticityAdjustment() const { return GameParameters::MinElasticityAdjustment; }
    float GetMaxElasticityAdjustment() const { return GameParameters::MaxElasticityAdjustment; }

    float GetStaticFrictionAdjustment() const { return mGameParameters.StaticFrictionAdjustment; }
    void SetStaticFrictionAdjustment(float value) { mGameParameters.StaticFrictionAdjustment = value; }
    float GetMinStaticFrictionAdjustment() const { return GameParameters::MinStaticFrictionAdjustment; }
    float GetMaxStaticFrictionAdjustment() const { return GameParameters::MaxStaticFrictionAdjustment; }

    float GetKineticFrictionAdjustment() const { return mGameParameters.KineticFrictionAdjustment; }
    void SetKineticFrictionAdjustment(float value) { mGameParameters.KineticFrictionAdjustment = value; }
    float GetMinKineticFrictionAdjustment() const { return GameParameters::MinKineticFrictionAdjustment; }
    float GetMaxKineticFrictionAdjustment() const { return GameParameters::MaxKineticFrictionAdjustment; }

    float GetMassAdjustment() const { return mMassAdjustment; }
    void SetMassAdjustment(float value) { mMassAdjustment = value; if (mWorld) mWorld->GetNpcs().OnMassAdjustmentChanged(mMassAdjustment); }
    float GetMinMassAdjustment() const { return 0.0001f; }
    float GetMaxMassAdjustment() const { return 100.0f; }

    float GetGravityAdjustment() const { return mGravityAdjustment; }
    void SetGravityAdjustment(float value) { mGravityAdjustment = value; if (mWorld) mWorld->GetNpcs().OnGravityAdjustmentChanged(mGravityAdjustment); }
    float GetMinGravityAdjustment() const { return 0.0f; }
    float GetMaxGravityAdjustment() const { return 100.0f; }

    float GetGlobalDamping() const { return mGameParameters.GlobalDamping; }
    void SetGlobalDamping(float value) { mGameParameters.GlobalDamping = value; }
    float GetMinGlobalDamping() const { return GameParameters::MinGlobalDamping; }
    float GetMaxGlobalDamping() const { return GameParameters::MaxGlobalDamping; }

    float GetSeaLevel() const { return mOceanDepth; }
    void SetSeaLevel(float value) { mOceanDepth = value; if (mWorld) mWorld->GetOceanSurface().SetDepth(value); }
    float GetMinSeaLevel() const { return Physics::OceanSurface::MinDepth; }
    float GetMaxSeaLevel() const { return Physics::OceanSurface::MaxDepth; }

    float GetSpringReductionFraction() const { return mGameParameters.NpcSpringReductionFraction; }
    void SetSpringReductionFraction(float value) { mGameParameters.NpcSpringReductionFraction = value; }
    float GetMinSpringReductionFraction() const { return GameParameters::MinNpcSpringReductionFraction; }
    float GetMaxSpringReductionFraction() const { return GameParameters::MaxNpcSpringReductionFraction; }

    float GetSpringDampingCoefficient() const { return mGameParameters.NpcSpringDampingCoefficient; }
    void SetSpringDampingCoefficient(float value) { mGameParameters.NpcSpringDampingCoefficient = value; }
    float GetMinSpringDampingCoefficient() const { return GameParameters::MinNpcSpringDampingCoefficient; }
    float GetMaxSpringDampingCoefficient() const { return GameParameters::MaxNpcSpringDampingCoefficient; }

    float GetWaterFrictionDragCoefficientAdjustment() const { return mGameParameters.WaterFrictionDragCoefficientAdjustment; }
    void SetWaterFrictionDragCoefficientAdjustment(float value) { mGameParameters.WaterFrictionDragCoefficientAdjustment = value; }
    float GetMinWaterFrictionDragCoefficientAdjustment() const { return GameParameters::MinWaterFrictionDragCoefficientAdjustment; }
    float GetMaxWaterFrictionDragCoefficientAdjustment() const { return GameParameters::MaxWaterFrictionDragCoefficientAdjustment; }

    float GetBuoyancyAdjustment() const { return mGameParameters.BuoyancyAdjustment; }
    void SetBuoyancyAdjustment(float value) { mGameParameters.BuoyancyAdjustment = value; }
    float GetMinBuoyancyAdjustment() const { return GameParameters::MinBuoyancyAdjustment; }
    float GetMaxBuoyancyAdjustment() const { return GameParameters::MaxBuoyancyAdjustment; }

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

    float GetHumanNpcBodyLengthAdjustment() const { return mGameParameters.HumanNpcBodyLengthAdjustment; }
    void SetHumanNpcBodyLengthAdjustment(float value) { mGameParameters.HumanNpcBodyLengthAdjustment = value; }
    float GetMinHumanNpcBodyLengthAdjustment() const { return GameParameters::MinHumanNpcBodyLengthAdjustment; }
    float GetMaxHumanNpcBodyLengthAdjustment() const { return GameParameters::MaxHumanNpcBodyLengthAdjustment; }

private:

    LabController(
        MaterialDatabase && materialDatabase,
        std::unique_ptr<RenderContext> renderContext);

    void Reset(
        std::unique_ptr<Physics::Ship> ship,
        std::filesystem::path const & shipDefinitionFilepath);

    void UpdateShipTransformations();

private:

    MaterialDatabase mMaterialDatabase;
    std::unique_ptr<RenderContext> mRenderContext;
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

    float mMassAdjustment;
    float mGravityAdjustment;
    float mOceanDepth;

    vec2f mCurrentShipTranslationVelocity;
    float mCurrentShipTranslationAccelerationIndicator;

    float mCurrentWavesAmplitude;
    float mTargetWavesAmplitude;
    float mCurrentWavesSpeed;
    float mTargetWavesSpeed;
    float mLastWaveRotationAngle;
    float mCurrentWaveTimeArg;

    int mCurrentVideoStep;

    //
    // Simulation control
    //

    SimulationControlStateType mSimulationControlState;
    bool mSimulationControlImpulse;
};
