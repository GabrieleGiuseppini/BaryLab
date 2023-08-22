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

    bool TrySelectOriginTriangle(vec2f const & screenCoordinates);

    std::optional<ElementIndex> TryPickParticle(vec2f const & screenCoordinates) const;

    void MoveParticleBy(ElementIndex particleIndex, vec2f const & screenOffset, vec2f const & inertialStride);

    void NotifyParticleTrajectory(ElementIndex particleIndex, vec2f const & targetScreenCoordinates);

    void SetParticleTrajectory(ElementIndex particleIndex, vec2f const & targetScreenCoordinates);

    void QueryNearestParticleAt(vec2f const & screenCoordinates) const;

    bool IsParticleGravityEnabled() const;

    void SetParticleGravityEnabled(bool isEnabled);

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

    vec2f WorldToScreen(vec2f const & worldCoordinates) const
    {
        assert(!!mRenderContext);
        return mRenderContext->WorldToScreen(worldCoordinates);
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

    float GetFriction() const { return mLabParameters.Friction; }
    void SetFriction(float value) { mLabParameters.Friction = value; }
    float GetMinFriction() const { return LabParameters::MinFriction; }
    float GetMaxFriction() const { return LabParameters::MaxFriction; }

    float GetMassAdjustment() const { return mLabParameters.MassAdjustment; }
    void SetMassAdjustment(float value) { mLabParameters.MassAdjustment = value; }
    float GetMinMassAdjustment() const { return LabParameters::MinMassAdjustment; }
    float GetMaxMassAdjustment() const { return LabParameters::MaxMassAdjustment; }

    float GetGravityAdjustment() const { return mLabParameters.GravityAdjustment; }
    void SetGravityAdjustment(float value) { mLabParameters.GravityAdjustment = value; }
    float GetMinGravityAdjustment() const { return LabParameters::MinGravityAdjustment; }
    float GetMaxGravityAdjustment() const { return LabParameters::MaxGravityAdjustment; }

private:

    explicit LabController(
        StructuralMaterialDatabase structuralMaterialDatabase,
        std::unique_ptr<RenderContext> renderContext);

    void Reset(
        std::unique_ptr<Model> newModel,
        std::filesystem::path const & meshDefinitionFilepath);

    ElementIndex FindTriangleContaining(vec2f const & position) const;

    //
    // Simulation
    //

    void InitializeParticleRegime(ElementIndex particleIndex);

    void UpdateSimulation(LabParameters const & labParameters);

    vec2f CalculatePhysicsDeltaPos(
        ElementIndex particleIndex,
        LabParameters const & labParameters) const;

    struct FinalPerticleState
    {
        vec2f Position;
        vec2f Velocity;

        FinalPerticleState(
            vec2f const & position,
            vec2f const & velocity)
            : Position(position)
            , Velocity(velocity)
        {}
    };

    std::optional<FinalPerticleState> UpdateParticleTrajectoryState(
        Particles::StateType & particleState,
        LabParameters const & labParameters);

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

    std::optional<ElementIndex> mCurrentlySelectedParticleProbe;
    std::optional<ElementIndex> mCurrentOriginTriangle;

    std::optional<ParticleTrajectory> mCurrentParticleTrajectory;
    std::optional<ParticleTrajectory> mCurrentParticleTrajectoryNotification;

    //
    // Simulation control
    //

    SimulationControlStateType mSimulationControlState;
    bool mSimulationControlImpulse;

    //
    // Our parameters
    //

    bool mIsGravityEnabled;
};
