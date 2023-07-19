/***************************************************************************************
* Original Author:      Gabriele Giuseppini
* Created:              2020-05-16
* Copyright:            Gabriele Giuseppini  (https://github.com/GabrieleGiuseppini)
***************************************************************************************/
#include "LabController.h"

#include "AABB.h"
#include "MeshBuilder.h"
#include "ResourceLocator.h"

#include <sstream>

std::unique_ptr<LabController> LabController::Create(
    int initialCanvasWidth,
    int initialCanvasHeight)
{
    LogMessage("InitialCanvasSize: ", initialCanvasWidth, "x", initialCanvasHeight);

    // Create render context
    std::unique_ptr<RenderContext> renderContext = std::make_unique<RenderContext>(
        initialCanvasWidth,
        initialCanvasHeight);

    //
    // Create controller
    //

    return std::unique_ptr<LabController>(
        new LabController(std::move(renderContext)));
}

LabController::LabController(std::unique_ptr<RenderContext> renderContext)
    : mEventDispatcher()
    , mRenderContext(std::move(renderContext))
    // Simulation state
    , mLabParameters()
    , mModel()
    , mCurrentMeshFilePath()
{    
}

void LabController::LoadMesh(std::filesystem::path const & meshDefinitionFilepath)
{
    // Load mesh
    std::unique_ptr<Mesh> mesh = MeshBuilder::BuildMesh(meshDefinitionFilepath);

    // Create a new model
    std::unique_ptr<Model> newModel = std::make_unique<Model>(
        std::move(mesh),
        mEventDispatcher);

    //
    // No errors, so we may continue
    //

    Reset(
        std::move(newModel),
        meshDefinitionFilepath);
}

void LabController::Update()
{
}

void LabController::Render()
{
    assert(!!mRenderContext);

    mRenderContext->RenderStart();

    if (mModel)
    {
        // TODOHERE

        //////
        ////// Triangles
        //////

        ////mRenderContext->UploadSpringsStart(mObject->GetSprings().GetElementCount());

        ////for (auto s : mObject->GetSprings())
        ////{
        ////    mRenderContext->UploadSpring(
        ////        mObject->GetPoints().GetPosition(mObject->GetSprings().GetEndpointAIndex(s)),
        ////        mObject->GetPoints().GetPosition(mObject->GetSprings().GetEndpointBIndex(s)),
        ////        mObject->GetSprings().GetRenderColor(s),
        ////        mObject->GetSprings().GetRenderNormThickness(s),
        ////        mObject->GetSprings().GetRenderHighlight(s));
        ////}

        ////mRenderContext->UploadSpringsEnd();

        //////
        ////// Particles
        //////

        ////mRenderContext->UploadPoints(
        ////    mObject->GetPoints().GetElementCount(),
        ////    mObject->GetPoints().GetPositionBuffer(),
        ////    mObject->GetPoints().GetRenderColorBuffer(),
        ////    mObject->GetPoints().GetRenderNormRadiusBuffer(),
        ////    mObject->GetPoints().GetRenderHighlightBuffer(),
        ////    mObject->GetPoints().GetFrozenCoefficientBuffer());
    }

    mRenderContext->RenderEnd();
}

void LabController::Reset()
{
    assert(mCurrentMeshFilePath);
    LoadMesh(*mCurrentMeshFilePath);
}

bool LabController::TrySelectOriginTriangle(vec2f const & screenCoordinates)
{
    // TODOHERE
}

std::optional<ElementIndex> LabController::TryPickParticle(vec2f const & screenCoordinates)
{
    // TODOHERE
}

void LabController::MoveParticleTo(ElementIndex particleIndex, vec2f const & targetScreenCoordinates)
{
    // TODOHERE
}

void LabController::SetParticleTrajectory(vec2f const & targetScreenCoordinates)
{
    // TODOHERE
}

void LabController::RunTrajectoryStep()
{
    // TODOHERE
}

////////////////////////////////////////////////

void LabController::Reset(
    std::unique_ptr<Model> newModel,
    std::filesystem::path const & meshDefinitionFilepath)
{
    //
    // Take object in
    //

    mModel.reset();
    mModel = std::move(newModel);
    mCurrentMeshFilePath = meshDefinitionFilepath;

    //
    // Auto-zoom & center
    //

    {
        AABB const objectAABB = mModel->GetMesh().GetPoints().GetAABB();

        vec2f const objectSize = objectAABB.GetSize();

        // Zoom to fit width and height (plus a nicely-looking margin)
        float const newZoom = std::min(
            mRenderContext->CalculateZoomForWorldWidth(objectSize.x + 5.0f),
            mRenderContext->CalculateZoomForWorldHeight(objectSize.y + 3.0f));
        mRenderContext->SetZoom(newZoom);

        // Center
        vec2f const objectCenter(
            (objectAABB.BottomLeft.x + objectAABB.TopRight.x) / 2.0f,
            (objectAABB.BottomLeft.y + objectAABB.TopRight.y) / 2.0f);
        mRenderContext->SetCameraWorldPosition(objectCenter);

    }

    //
    // Publish reset
    //

    mEventDispatcher.OnReset();
}
