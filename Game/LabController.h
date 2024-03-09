/***************************************************************************************
* Original Author:		Gabriele Giuseppini
* Created:				2020-05-16
* Copyright:			Gabriele Giuseppini  (https://github.com/GabrieleGiuseppini)
***************************************************************************************/
#pragma once

#include "GameEventDispatcher.h"
#include "GameParameters.h"
#include "MaterialDatabase.h"
#include "Model.h"
#include "Physics.h"
#include "RenderContext.h"
#include "ShipDefinition.h"

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
        mGameEventDispatcher.RegisterBLabEventHandler(handler);
    }

    //
    // Lab
    //

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

    bool IsGravityEnabled() const;

    void SetGravityEnabled(bool isEnabled);

    vec2f const & GetShipVelocity() const;

    void SetShipVelocity(vec2f const & velocity);

    void SetNpcPanicLevelForAllHumans(float panicLevel);

    void DoStepForVideo();

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

    NpcRenderMode GetNpcRenderMode() const
    {
        assert(mModel);
        return mModel->GetNpcs().GetNpcRenderMode();
    }

    void SetNpcRenderMode(NpcRenderMode value)
    {
        assert(mModel);
        mModel->GetNpcs().SetNpcRenderMode(value);
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

    float GetMassAdjustment() const { return mGameParameters.MassAdjustment; }
    void SetMassAdjustment(float value) { mGameParameters.MassAdjustment = value; }
    float GetMinMassAdjustment() const { return GameParameters::MinMassAdjustment; }
    float GetMaxMassAdjustment() const { return GameParameters::MaxMassAdjustment; }

    float GetGravityAdjustment() const { return mGameParameters.GravityAdjustment; }
    void SetGravityAdjustment(float value) { mGameParameters.GravityAdjustment = value; }
    float GetMinGravityAdjustment() const { return GameParameters::MinGravityAdjustment; }
    float GetMaxGravityAdjustment() const { return GameParameters::MaxGravityAdjustment; }

    float GetGlobalDamping() const { return mGameParameters.GlobalDamping; }
    void SetGlobalDamping(float value) { mGameParameters.GlobalDamping = value; }
    float GetMinGlobalDamping() const { return GameParameters::MinGlobalDamping; }
    float GetMaxGlobalDamping() const { return GameParameters::MaxGlobalDamping; }

    float GetSeaLevel() const { return mWorld.GetOceanSurface().GetDepth(); }
    void SetSeaLevel(float value) { mWorld.GetOceanSurface().SetDepth(value); }
    float GetMinSeaLevel() const { return Physics::OceanSurface::MinDepth; }
    float GetMaxSeaLevel() const { return Physics::OceanSurface::MaxDepth; }

    float GetSpringReductionFraction() const { return mGameParameters.SpringReductionFraction; }
    void SetSpringReductionFraction(float value) { mGameParameters.SpringReductionFraction = value; }
    float GetMinSpringReductionFraction() const { return GameParameters::MinSpringReductionFraction; }
    float GetMaxSpringReductionFraction() const { return GameParameters::MaxSpringReductionFraction; }

    float GetSpringDampingCoefficient() const { return mGameParameters.SpringDampingCoefficient; }
    void SetSpringDampingCoefficient(float value) { mGameParameters.SpringDampingCoefficient = value; }
    float GetMinSpringDampingCoefficient() const { return GameParameters::MinSpringDampingCoefficient; }
    float GetMaxSpringDampingCoefficient() const { return GameParameters::MaxSpringDampingCoefficient; }

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

    float GetHumanNpcWalkingSpeed() const { return mGameParameters.HumanNpcWalkingSpeed; }
    void SetHumanNpcWalkingSpeed(float value) { mGameParameters.HumanNpcWalkingSpeed = value; }
    float GetMinHumanNpcWalkingSpeed() const { return GameParameters::MinHumanNpcWalkingSpeed; }
    float GetMaxHumanNpcWalkingSpeed() const { return GameParameters::MaxHumanNpcWalkingSpeed; }

private:

    LabController(
        MaterialDatabase && materialDatabase,
        std::unique_ptr<RenderContext> renderContext);

    void LoadShip(
        std::filesystem::path const & shipDefinitionFilepath,
        bool addExperimentalNpc);

    void Reset(
        std::unique_ptr<Model> newModel,
        std::filesystem::path const & shipDefinitionFilepath);

    void UpdateShipTransformations();

private:

    MaterialDatabase mMaterialDatabase;
    std::unique_ptr<RenderContext> mRenderContext;
    GameEventDispatcher mGameEventDispatcher;

    //
    // Current state
    //

    GameParameters mGameParameters;
    std::unique_ptr<Model> mModel;
    Physics::World mWorld; // Dummy placeholder
    std::optional<std::filesystem::path> mCurrentShipFilePath;

    float mCurrentSimulationTime;

    bool mIsGravityEnabled;
    vec2f mCurrentShipTranslationVelocity;
    float mCurrentShipTranslationAccelerationIndicator;

    int mCurrentVideoStep;

    //
    // Simulation control
    //

    SimulationControlStateType mSimulationControlState;
    bool mSimulationControlImpulse;
};
