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

    // Load materials
    StructuralMaterialDatabase structuralMaterialDatabase = StructuralMaterialDatabase::Load();

    // Create render context
    std::unique_ptr<RenderContext> renderContext = std::make_unique<RenderContext>(
        initialCanvasWidth,
        initialCanvasHeight);

    //
    // Create controller
    //

    return std::unique_ptr<LabController>(
        new LabController(
        std::move(structuralMaterialDatabase),
            std::move(renderContext)));
}

LabController::LabController(
    StructuralMaterialDatabase structuralMaterialDatabase,
    std::unique_ptr<RenderContext> renderContext)
    : mStructuralMaterialDatabase(std::move(structuralMaterialDatabase))
    , mRenderContext(std::move(renderContext))
    , mEventDispatcher()
    // Simulation state
    , mLabParameters()
    , mModel()
    , mCurrentMeshFilePath()
    , mCurrentOriginTriangle()
    // Simulation control
    , mSimulationControlState(SimulationControlStateType::Paused)
    , mSimulationControlImpulse(false)
{    
}

void LabController::SetSimulationControlState(SimulationControlStateType state)
{
    mSimulationControlState = state;
}

void LabController::SetSimulationControlPulse()
{
    mSimulationControlImpulse = true;
}

void LabController::LoadMesh(std::filesystem::path const & meshDefinitionFilepath)
{
    // Load mesh definition
    auto meshDefinition = MeshDefinition::Load(meshDefinitionFilepath);

    // Make mesh
    std::unique_ptr<Mesh> mesh = MeshBuilder::BuildMesh(
        std::move(meshDefinition),
        mStructuralMaterialDatabase);

    // Create particles
    std::unique_ptr<Particles> particles = std::make_unique<Particles>(1);
    particles->Add(vec2f(0.25f, 0.25f), rgbaColor(0x60, 0x60, 0x60, 0xff));

    // Create a new model
    std::unique_ptr<Model> newModel = std::make_unique<Model>(
        std::move(mesh),
        std::move(particles),
        mEventDispatcher);

    // Select first triangle as origin triangle
    mCurrentOriginTriangle = 0;

    //
    // No errors, so we may continue
    //

    Reset(
        std::move(newModel),
        meshDefinitionFilepath);
}

void LabController::Update()
{
    assert(mModel);

    if (mSimulationControlState == SimulationControlStateType::Play
        || mSimulationControlImpulse)
    {
        UpdateSimulation(mLabParameters);

        // Update state
        mSimulationControlImpulse = false;
    }

    std::optional<vec3f> barycentricCoordinates;
    if (mCurrentOriginTriangle)
    {
        barycentricCoordinates = mModel->GetMesh().GetTriangles().ToBarycentricCoordinates(
            mModel->GetParticles().GetPosition(0),
            *mCurrentOriginTriangle, 
            mModel->GetMesh().GetVertices());
    }

    mEventDispatcher.OnSubjectParticleBarycentricCoordinatesChanged(barycentricCoordinates);
}

void LabController::Render()
{
    assert(!!mRenderContext);

    mRenderContext->RenderStart();

    if (mModel)
    {
        //
        // Vertices
        //

        mRenderContext->UploadVertices(
            mModel->GetMesh().GetVertices().GetElementCount(),
            mModel->GetMesh().GetVertices().GetPositionBuffer());

        //
        // Edges
        //

        mRenderContext->UploadEdgesStart(mModel->GetMesh().GetEdges().GetElementCount());

        for (auto e : mModel->GetMesh().GetEdges())
        {
            mRenderContext->UploadEdge(
                mModel->GetMesh().GetEdges().GetEndpointAPosition(e, mModel->GetMesh().GetVertices()),
                mModel->GetMesh().GetEdges().GetEndpointBPosition(e, mModel->GetMesh().GetVertices()));
        }

        mRenderContext->UploadEdgesEnd();

        //
        // Particles
        //

        mRenderContext->UploadParticlesStart(mModel->GetParticles().GetElementCount());

        for (auto p : mModel->GetParticles())
        {
            mRenderContext->UploadParticle(
                mModel->GetParticles().GetPosition(p),
                mModel->GetParticles().GetRenderColor(p));
        }

        mRenderContext->UploadParticlesEnd();
    }

    mRenderContext->RenderEnd();
}

void LabController::Reset()
{
    assert(mCurrentMeshFilePath);
    LoadMesh(*mCurrentMeshFilePath);
}

std::optional<ElementIndex> LabController::TryPickVertex(vec2f const & screenCoordinates) const
{
    assert(!!mModel);

    //
    // Find closest vertex within the radius
    //

    vec2f const worldCoordinates = ScreenToWorld(screenCoordinates);

    float constexpr SquareSearchRadius = LabParameters::VertexRadius * LabParameters::VertexRadius;

    float bestSquareDistance = std::numeric_limits<float>::max();
    ElementIndex bestVertex = NoneElementIndex;

    auto const & vertices = mModel->GetMesh().GetVertices();
    for (auto v : vertices)
    {
        float const squareDistance = (vertices.GetPosition(v) - worldCoordinates).squareLength();
        if (squareDistance < SquareSearchRadius
            && squareDistance < bestSquareDistance)
        {
            bestSquareDistance = squareDistance;
            bestVertex = v;
        }
    }

    if (bestVertex != NoneElementIndex)
        return bestVertex;
    else
        return std::nullopt;
}

void LabController::MoveVertexBy(
    ElementIndex vertexIndex, 
    vec2f const & screenOffset)
{
    assert(!!mModel);

    vec2f const worldOffset = ScreenOffsetToWorldOffset(screenOffset);

    mModel->GetMesh().GetVertices().SetPosition(
        vertexIndex,
        mModel->GetMesh().GetVertices().GetPosition(vertexIndex) + worldOffset);
}

bool LabController::TrySelectOriginTriangle(vec2f const & screenCoordinates)
{
    // TODOHERE
    (void)screenCoordinates;
    return false;
}

std::optional<ElementIndex> LabController::TryPickParticle(vec2f const & screenCoordinates) const
{
    assert(!!mModel);

    //
    // Find closest particle within the radius
    //

    vec2f const worldCoordinates = ScreenToWorld(screenCoordinates);

    float constexpr SquareSearchRadius = LabParameters::ParticleRadius * LabParameters::ParticleRadius;

    float bestSquareDistance = std::numeric_limits<float>::max();
    ElementIndex bestParticle = NoneElementIndex;

    auto const & particles = mModel->GetParticles();
    for (auto p : particles)
    {
        float const squareDistance = (particles.GetPosition(p) - worldCoordinates).squareLength();
        if (squareDistance < SquareSearchRadius
            && squareDistance < bestSquareDistance)
        {
            bestSquareDistance = squareDistance;
            bestParticle = p;
        }
    }

    if (bestParticle != NoneElementIndex)
        return bestParticle;
    else
        return std::nullopt;
}

void LabController::MoveParticleBy(
    ElementIndex particleIndex, 
    vec2f const & screenOffset,
    vec2f const & inertialStride)
{
    assert(!!mModel);

    vec2f const worldOffset = ScreenOffsetToWorldOffset(screenOffset);
    vec2f const worldStride = ScreenOffsetToWorldOffset(inertialStride);

    mModel->GetParticles().SetPosition(
        particleIndex, 
        mModel->GetParticles().GetPosition(particleIndex) + worldOffset);

    mModel->GetParticles().SetVelocity(
        particleIndex, 
        worldStride / LabParameters::SimulationTimeStepDuration * 0.5f); // Magic adjustment
}

void LabController::SetParticleTrajectory(
    ElementIndex particleIndex,
    vec2f const & targetScreenCoordinates)
{
    // TODOHERE
    (void)particleIndex;
    (void)targetScreenCoordinates;
}

void LabController::QueryNearestParticleAt(vec2f const & screenCoordinates) const
{
    assert(!!mModel);

    auto const nearestParticle = TryPickParticle(screenCoordinates);
    if (nearestParticle.has_value())
    {
        mModel->GetParticles().Query(*nearestParticle);
    }
}

void LabController::SetParticleGravityEnabled(bool isEnabled)
{
    vec2f const gravity = isEnabled ? LabParameters::Gravity : vec2f::zero();

    auto & particles = mModel->GetParticles();
    for (auto p : particles)
    {
        particles.SetWorldForce(p, gravity);
    }
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
        AABB const objectAABB = mModel->GetMesh().GetVertices().GetAABB();

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
