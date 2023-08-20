/***************************************************************************************
* Original Author:      Gabriele Giuseppini
* Created:              2020-05-16
* Copyright:            Gabriele Giuseppini  (https://github.com/GabrieleGiuseppini)
***************************************************************************************/
#include "LabController.h"

#include "AABB.h"
#include "Geometry.h"
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
    , mCurrentParticleTrajectory()
    , mCurrentParticleTrajectoryNotification()
    // Simulation control
    , mSimulationControlState(SimulationControlStateType::Paused)
    , mSimulationControlImpulse(false)
    // Our own parameters
    , mIsGravityEnabled(false)
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

    if (mesh->GetTriangles().GetElementCount() < 1)
    {
        throw BLabException("Mesh must contain at least one triangle");
    }

    // Create particle
    std::unique_ptr<Particles> particles = std::make_unique<Particles>(1);
    {
        // TODOTEST
        ////vec2f const center = (
        ////    mesh->GetVertices().GetPosition(mesh->GetTriangles().GetVertexAIndex(0))
        ////    + mesh->GetVertices().GetPosition(mesh->GetTriangles().GetVertexBIndex(0))
        ////    + mesh->GetVertices().GetPosition(mesh->GetTriangles().GetVertexCIndex(0))) / 3.0f;
        vec2f const center =
            mesh->GetVertices().GetPosition(mesh->GetTriangles().GetVertexBIndex(0))
            - vec2f(0.5f, 0.5f);

        particles->Add(center, rgbaColor(0x60, 0x60, 0x60, 0xff));
    }    

    // Create a new model
    std::unique_ptr<Model> newModel = std::make_unique<Model>(
        std::move(mesh),
        std::move(particles),
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
    assert(mModel);

    // Reconcile gravity

    vec2f const gravity = mIsGravityEnabled ? LabParameters::Gravity : vec2f::zero();
    auto & particles = mModel->GetParticles();
    for (auto p : particles)
    {
        particles.SetWorldForce(p, gravity);
    }

    // Update simulation

    if (mSimulationControlState == SimulationControlStateType::Play
        || mSimulationControlImpulse)
    {
        UpdateSimulation(mLabParameters);

        // Update state
        mSimulationControlImpulse = false;
    }

    // Publish particle probe, if constrained

    std::optional<ParticleProbe> particleProbe;
    if (mCurrentlySelectedParticleProbe.has_value()
        && mModel->GetParticles().GetState(*mCurrentlySelectedParticleProbe).ConstrainedState.has_value())
    {
        particleProbe.emplace(
            mModel->GetParticles().GetState(*mCurrentlySelectedParticleProbe).ConstrainedState->CurrentTriangle,
            mModel->GetParticles().GetState(*mCurrentlySelectedParticleProbe).ConstrainedState->CurrentTriangleBarycentricCoords);
    }

    mEventDispatcher.OnSubjectParticleUpdated(particleProbe);

    // Publish barycentric coords wrt origin

    std::optional<vec3f> barycentricCoordinates;
    if (mCurrentOriginTriangle)
    {
        barycentricCoordinates = mModel->GetMesh().GetTriangles().ToBarycentricCoordinates(
            mModel->GetParticles().GetPosition(0),
            *mCurrentOriginTriangle, 
            mModel->GetMesh().GetVertices());
    }

    mEventDispatcher.OnSubjectParticleBarycentricCoordinatesWrtOriginTriangleChanged(barycentricCoordinates);
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
        // Edges (of triangles)
        //
        // Note: for lazyness we load lots of repetitions
        //

        mRenderContext->UploadEdgesStart();

        for (auto t : mModel->GetMesh().GetTriangles())
        {
            mRenderContext->UploadEdge(
                mModel->GetMesh().GetVertices().GetPosition(mModel->GetMesh().GetTriangles().GetVertexAIndex(t)),
                mModel->GetMesh().GetVertices().GetPosition(mModel->GetMesh().GetTriangles().GetVertexBIndex(t)),
                mModel->GetMesh().GetEdges().GetRenderColor(mModel->GetMesh().GetTriangles().GetSubEdgeAIndex(t)));

            mRenderContext->UploadEdge(
                mModel->GetMesh().GetVertices().GetPosition(mModel->GetMesh().GetTriangles().GetVertexBIndex(t)),
                mModel->GetMesh().GetVertices().GetPosition(mModel->GetMesh().GetTriangles().GetVertexCIndex(t)),
                mModel->GetMesh().GetEdges().GetRenderColor(mModel->GetMesh().GetTriangles().GetSubEdgeBIndex(t)));

            mRenderContext->UploadEdge(
                mModel->GetMesh().GetVertices().GetPosition(mModel->GetMesh().GetTriangles().GetVertexCIndex(t)),
                mModel->GetMesh().GetVertices().GetPosition(mModel->GetMesh().GetTriangles().GetVertexAIndex(t)),
                mModel->GetMesh().GetEdges().GetRenderColor(mModel->GetMesh().GetTriangles().GetSubEdgeCIndex(t)));
        }

        mRenderContext->UploadEdgesEnd();

        //
        // Particles
        //

        mRenderContext->UploadParticlesStart();

        for (auto p : mModel->GetParticles())
        {
            mRenderContext->UploadParticle(
                mModel->GetParticles().GetPosition(p),
                mModel->GetParticles().GetRenderColor(p));
        }

        mRenderContext->UploadParticlesEnd();

        //
        // Particle trajectory
        //

        mRenderContext->UploadParticleTrajectoriesStart();

        if (mCurrentParticleTrajectoryNotification)
        {
            mRenderContext->UploadParticleTrajectory(
                mModel->GetParticles().GetPosition(mCurrentParticleTrajectoryNotification->ParticleIndex),
                mCurrentParticleTrajectoryNotification->TargetPosition,
                rgbaColor(0xc0, 0xc0, 0xc0, 0xff));
        }

        if (mCurrentParticleTrajectory)
        {
            mRenderContext->UploadParticleTrajectory(
                mModel->GetParticles().GetPosition(mCurrentParticleTrajectory->ParticleIndex),
                mCurrentParticleTrajectory->TargetPosition,
                rgbaColor(0x99, 0x99, 0x99, 0xff));
        }

        mRenderContext->UploadParticleTrajectoriesEnd();

        //
        // Selected triangles
        //

        mRenderContext->UploadSelectedTrianglesStart();

        if (mCurrentOriginTriangle)
        {
            mRenderContext->UploadSelectedTriangle(
                mModel->GetMesh().GetVertices().GetPosition(mModel->GetMesh().GetTriangles().GetVertexAIndex(*mCurrentOriginTriangle)),
                mModel->GetMesh().GetVertices().GetPosition(mModel->GetMesh().GetTriangles().GetVertexBIndex(*mCurrentOriginTriangle)),
                mModel->GetMesh().GetVertices().GetPosition(mModel->GetMesh().GetTriangles().GetVertexCIndex(*mCurrentOriginTriangle)),
                rgbaColor(227, 107, 107, 77));
        }

        if (mCurrentlySelectedParticleProbe.has_value()
            && mModel->GetParticles().GetState(*mCurrentlySelectedParticleProbe).ConstrainedState.has_value())
        {
            auto const t = mModel->GetParticles().GetState(*mCurrentlySelectedParticleProbe).ConstrainedState->CurrentTriangle;

            mRenderContext->UploadSelectedTriangle(
                mModel->GetMesh().GetVertices().GetPosition(mModel->GetMesh().GetTriangles().GetVertexAIndex(t)),
                mModel->GetMesh().GetVertices().GetPosition(mModel->GetMesh().GetTriangles().GetVertexBIndex(t)),
                mModel->GetMesh().GetVertices().GetPosition(mModel->GetMesh().GetTriangles().GetVertexCIndex(t)),
                rgbaColor(107, 227, 107, 77));
        }

        mRenderContext->UploadSelectedTrianglesEnd();        
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

    assert(mModel->GetParticles().GetElementCount() >= 1);
    InitializeParticleRegime(0);
}

bool LabController::TrySelectOriginTriangle(vec2f const & screenCoordinates)
{
    assert(!!mModel);

    vec2f const worldCoordinates = ScreenToWorld(screenCoordinates);

    ElementIndex const t = FindTriangleContaining(worldCoordinates);
    if (t != NoneElementIndex)
    {
        mCurrentOriginTriangle = t;
        return true;
    }
    else
    {
        mCurrentOriginTriangle.reset();
        return false;
    }
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
        vec2f::zero()); // Zero-out velocity

    InitializeParticleRegime(particleIndex);

    // Select particle
    mCurrentlySelectedParticleProbe.emplace(particleIndex);

    // Re-calcualate trajectory
    mModel->GetParticles().GetState(particleIndex).TargetPosition.reset();
    mCurrentParticleTrajectory.reset();
    mCurrentParticleTrajectory.reset();
}

void LabController::NotifyParticleTrajectory(
    ElementIndex particleIndex,
    vec2f const & targetScreenCoordinates)
{
    mCurrentParticleTrajectoryNotification.emplace(
        particleIndex,
        ScreenToWorld(targetScreenCoordinates));

    mCurrentParticleTrajectory.reset();
}

void LabController::SetParticleTrajectory(
    ElementIndex particleIndex,
    vec2f const & targetScreenCoordinates)
{
    mCurrentParticleTrajectory.emplace(
        particleIndex,
        ScreenToWorld(targetScreenCoordinates));
    
    mCurrentParticleTrajectoryNotification.reset();

    // Reset state to needing to calculate a trajectory
    mModel->GetParticles().GetState(particleIndex).TargetPosition.reset();
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

bool LabController::IsParticleGravityEnabled() const
{
    return mIsGravityEnabled;
}

void LabController::SetParticleGravityEnabled(bool isEnabled)
{
    mIsGravityEnabled = isEnabled;
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

    // Reset state
    assert(mModel->GetParticles().GetElementCount() == 1);
    InitializeParticleRegime(0);
    mCurrentlySelectedParticleProbe.emplace(0);
    mCurrentOriginTriangle.reset();
    mCurrentParticleTrajectory.reset();
    mCurrentParticleTrajectoryNotification.reset();    

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

ElementIndex LabController::FindTriangleContaining(vec2f const & position) const
{
    for (auto const t : mModel->GetMesh().GetTriangles())
    {
        if (Geometry::IsPointInTriangle(
            position,
            mModel->GetMesh().GetVertices().GetPosition(mModel->GetMesh().GetTriangles().GetVertexAIndex(t)),
            mModel->GetMesh().GetVertices().GetPosition(mModel->GetMesh().GetTriangles().GetVertexBIndex(t)),
            mModel->GetMesh().GetVertices().GetPosition(mModel->GetMesh().GetTriangles().GetVertexCIndex(t))))
        {
            return t;
        }
    }

    return NoneElementIndex;
}