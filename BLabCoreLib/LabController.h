/***************************************************************************************
* Original Author:		Gabriele Giuseppini
* Created:				2020-05-16
* Copyright:			Gabriele Giuseppini  (https://github.com/GabrieleGiuseppini)
***************************************************************************************/
#pragma once

#include "BLabTypes.h"
#include "EventDispatcher.h"
#include "LabParameters.h"
#include "MeshDefinition.h"
#include "Model.h"
#include "RenderContext.h"
#include "StructuralMaterialDatabase.h"
#include "Vectors.h"

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

    void RegisterEventHandler(IEventHandler * handler)
    {
        mEventDispatcher.RegisterEventHandler(handler);
    }

    //
    // Lab
    //

    void SetSimulationControlState(SimulationControlStateType state);

    void SetSimulationControlPulse();

    void LoadMesh(std::filesystem::path const & meshDefinitionFilepath);

    void Update();

    void Render();

    void Reset();

    //
    // Interactions
    //

    std::optional<ElementIndex> TryPickVertex(vec2f const & screenCoordinates) const;

    void MoveVertexBy(ElementIndex vertexIndex, vec2f const & screenOffset);

    void RotateMeshBy(vec2f const & centerScreenCoordinates, float screenAngle);

    void RotateMeshBy(ElementIndex particleIndex, float screenAngle);

    bool TrySelectOriginTriangle(vec2f const & screenCoordinates);

    bool TrySelectParticle(vec2f const & screenCoordinates);

    std::optional<ElementIndex> TryPickNpcParticle(vec2f const & screenCoordinates) const;

    void MoveNpcParticleBy(ElementIndex npcParticleIndex, vec2f const & screenOffset, vec2f const & inertialStride);

    void NotifyNpcParticleTrajectory(ElementIndex npcParticleIndex, vec2f const & targetScreenCoordinates);

    void SetNpcParticleTrajectory(ElementIndex npcParticleIndex, vec2f const & targetScreenCoordinates);

    void QueryNearestNpcParticleAt(vec2f const & screenCoordinates) const;

    bool IsGravityEnabled() const;

    void SetGravityEnabled(bool isEnabled);

    vec2f const & GetMeshVelocity() const;

    void SetMeshVelocity(vec2f const & velocity);

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

    //
    // Simmulation parameters
    //

    float GetElasticity() const { return mLabParameters.Elasticity; }
    void SetElasticity(float value) { mLabParameters.Elasticity = value; }
    float GetMinElasticity() const { return LabParameters::MinElasticity; }
    float GetMaxElasticity() const { return LabParameters::MaxElasticity; }

    float GetStaticFrictionAdjustment() const { return mLabParameters.StaticFrictionAdjustment; }
    void SetStaticFrictionAdjustment(float value) { mLabParameters.StaticFrictionAdjustment = value; }
    float GetMinStaticFrictionAdjustment() const { return LabParameters::MinStaticFrictionAdjustment; }
    float GetMaxStaticFrictionAdjustment() const { return LabParameters::MaxStaticFrictionAdjustment; }

    float GetKineticFrictionAdjustment() const { return mLabParameters.KineticFrictionAdjustment; }
    void SetKineticFrictionAdjustment(float value) { mLabParameters.KineticFrictionAdjustment = value; }
    float GetMinKineticFrictionAdjustment() const { return LabParameters::MinKineticFrictionAdjustment; }
    float GetMaxKineticFrictionAdjustment() const { return LabParameters::MaxKineticFrictionAdjustment; }

    float GetMassAdjustment() const { return mLabParameters.MassAdjustment; }
    void SetMassAdjustment(float value) { mLabParameters.MassAdjustment = value; }
    float GetMinMassAdjustment() const { return LabParameters::MinMassAdjustment; }
    float GetMaxMassAdjustment() const { return LabParameters::MaxMassAdjustment; }

    float GetGravityAdjustment() const { return mLabParameters.GravityAdjustment; }
    void SetGravityAdjustment(float value) { mLabParameters.GravityAdjustment = value; }
    float GetMinGravityAdjustment() const { return LabParameters::MinGravityAdjustment; }
    float GetMaxGravityAdjustment() const { return LabParameters::MaxGravityAdjustment; }

    float GetGlobalDamping() const { return mLabParameters.GlobalDamping; }
    void SetGlobalDamping(float value) { mLabParameters.GlobalDamping = value; }
    float GetMinGlobalDamping() const { return LabParameters::MinGlobalDamping; }
    float GetMaxGlobalDamping() const { return LabParameters::MaxGlobalDamping; }

    float GetSeaLevel() const { return mLabParameters.SeaLevel; }
    void SetSeaLevel(float value) { mLabParameters.SeaLevel = value; }
    float GetMinSeaLevel() const { return LabParameters::MinSeaLevel; }
    float GetMaxSeaLevel() const { return LabParameters::MaxSeaLevel; }

    float GetSpringReductionFraction() const { return mLabParameters.SpringReductionFraction; }
    void SetSpringReductionFraction(float value) { mLabParameters.SpringReductionFraction = value; }
    float GetMinSpringReductionFraction() const { return LabParameters::MinSpringReductionFraction; }
    float GetMaxSpringReductionFraction() const { return LabParameters::MaxSpringReductionFraction; }

    float GetSpringDampingCoefficient() const { return mLabParameters.SpringDampingCoefficient; }
    void SetSpringDampingCoefficient(float value) { mLabParameters.SpringDampingCoefficient = value; }
    float GetMinSpringDampingCoefficient() const { return LabParameters::MinSpringDampingCoefficient; }
    float GetMaxSpringDampingCoefficient() const { return LabParameters::MaxSpringDampingCoefficient; }

    float GetWaterFrictionDragCoefficientAdjustment() const { return mLabParameters.WaterFrictionDragCoefficientAdjustment; }
    void SetWaterFrictionDragCoefficientAdjustment(float value) { mLabParameters.WaterFrictionDragCoefficientAdjustment = value; }
    float GetMinWaterFrictionDragCoefficientAdjustment() const { return LabParameters::MinWaterFrictionDragCoefficientAdjustment; }
    float GetMaxWaterFrictionDragCoefficientAdjustment() const { return LabParameters::MaxWaterFrictionDragCoefficientAdjustment; }

    float GetBuoyancyAdjustment() const { return mLabParameters.BuoyancyAdjustment; }
    void SetBuoyancyAdjustment(float value) { mLabParameters.BuoyancyAdjustment = value; }
    float GetMinBuoyancyAdjustment() const { return LabParameters::MinBuoyancyAdjustment; }
    float GetMaxBuoyancyAdjustment() const { return LabParameters::MaxBuoyancyAdjustment; }

    float GetHumanNpcEquilibriumTorqueStiffnessCoefficient() const { return mLabParameters.HumanNpcEquilibriumTorqueStiffnessCoefficient; }
    void SetHumanNpcEquilibriumTorqueStiffnessCoefficient(float value) { mLabParameters.HumanNpcEquilibriumTorqueStiffnessCoefficient = value; }
    float GetMinHumanNpcEquilibriumTorqueStiffnessCoefficient() const { return LabParameters::MinHumanNpcEquilibriumTorqueStiffnessCoefficient; }
    float GetMaxHumanNpcEquilibriumTorqueStiffnessCoefficient() const { return LabParameters::MaxHumanNpcEquilibriumTorqueStiffnessCoefficient; }

    float GetHumanNpcEquilibriumTorqueDampingCoefficient() const { return mLabParameters.HumanNpcEquilibriumTorqueDampingCoefficient; }
    void SetHumanNpcEquilibriumTorqueDampingCoefficient(float value) { mLabParameters.HumanNpcEquilibriumTorqueDampingCoefficient = value; }
    float GetMinHumanNpcEquilibriumTorqueDampingCoefficient() const { return LabParameters::MinHumanNpcEquilibriumTorqueDampingCoefficient; }
    float GetMaxHumanNpcEquilibriumTorqueDampingCoefficient() const { return LabParameters::MaxHumanNpcEquilibriumTorqueDampingCoefficient; }

    float GetHumanNpcWalkingSpeed() const { return mLabParameters.HumanNpcWalkingSpeed; }
    void SetHumanNpcWalkingSpeed(float value) { mLabParameters.HumanNpcWalkingSpeed = value; }
    float GetMinHumanNpcWalkingSpeed() const { return LabParameters::MinHumanNpcWalkingSpeed; }
    float GetMaxHumanNpcWalkingSpeed() const { return LabParameters::MaxHumanNpcWalkingSpeed; }

private:

    explicit LabController(
        StructuralMaterialDatabase && structuralMaterialDatabase,
        std::unique_ptr<RenderContext> renderContext);

    void Reset(
        std::unique_ptr<Model> newModel,
        std::filesystem::path const & meshDefinitionFilepath);

    void UpdateMeshTransformations();

private:
    
    StructuralMaterialDatabase mStructuralMaterialDatabase;
    std::unique_ptr<RenderContext> mRenderContext;
    EventDispatcher mEventDispatcher;

    //
    // Current state
    //

    LabParameters mLabParameters;
    std::unique_ptr<Model> mModel;
    std::optional<std::filesystem::path> mCurrentMeshFilePath;

    bool mIsGravityEnabled;
    vec2f mCurrentMeshTranslationVelocity;
    float mCurrentMeshTranslationAccelerationIndicator;

    //
    // Simulation control
    //

    SimulationControlStateType mSimulationControlState;
    bool mSimulationControlImpulse;
};
