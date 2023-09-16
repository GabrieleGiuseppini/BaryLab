/***************************************************************************************
* Original Author:      Gabriele Giuseppini
* Created:              2020-05-16
* Copyright:            Gabriele Giuseppini  (https://github.com/GabrieleGiuseppini)
***************************************************************************************/
#include "LabController.h"

#include "AABB.h"
#include "MeshBuilder.h"
#include "ResourceLocator.h"

#include <cmath>
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
    , mIsGravityEnabled(false)
    , mCurrentMeshTranslationVelocity(vec2f::zero())
    , mCurrentMeshTranslationAccelerationIndicator(0.0f)
    , mRenderSimulationSteps(true)
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
        ////vec2f const center =
        ////    mesh->GetVertices().GetPosition(mesh->GetTriangles().GetVertexBIndex(0))
        ////    - vec2f(0.5f, 0.5f);
        vec2f const center = vec2f(0.5f, -2.0f);

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
        particles.SetWorldForce(p, gravity * mLabParameters.ParticleMass * mLabParameters.MassAdjustment);
    }

    // Update simulation

    if (mSimulationControlState == SimulationControlStateType::Play
        || mSimulationControlImpulse)
    {
        UpdateSimulation(mLabParameters);

        mCurrentMeshTranslationAccelerationIndicator *= 0.98f;

        // Update state
        mSimulationControlImpulse = false;
    }

    // Publish particle probe, if constrained

    std::optional<ConstrainedRegimeParticleProbe> constrainedRegimeParticleProbe;
    if (mCurrentlySelectedParticleProbe.has_value()
        && mModel->GetParticles().GetState(*mCurrentlySelectedParticleProbe).ConstrainedState.has_value())
    {
        constrainedRegimeParticleProbe.emplace(
            mModel->GetParticles().GetState(*mCurrentlySelectedParticleProbe).ConstrainedState->CurrentTriangle,
            mModel->GetParticles().GetState(*mCurrentlySelectedParticleProbe).ConstrainedState->CurrentTriangleBarycentricCoords);
    }

    mEventDispatcher.OnSubjectParticleConstrainedRegimeUpdated(constrainedRegimeParticleProbe);

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

    // Publish particle physics

    std::optional<PhysicsParticleProbe> physicsParticleProbe;
    if (mCurrentlySelectedParticleProbe.has_value())
    {
        physicsParticleProbe.emplace(mModel->GetParticles().GetVelocity(0));
    }

    mEventDispatcher.OnSubjectParticlePhysicsUpdated(physicsParticleProbe);
}

void LabController::Render()
{
    assert(!!mRenderContext);

    mRenderContext->RenderStart();

    if (mModel)
    {
        auto const & vertices = mModel->GetMesh().GetVertices();
        auto const & edges = mModel->GetMesh().GetEdges();
        auto const & triangles = mModel->GetMesh().GetTriangles();

        //
        // Vertices
        //

        mRenderContext->UploadVertices(
            vertices.GetElementCount(),
            vertices.GetPositionBuffer());

        //
        // Edges (of triangles)
        //
        // Note: for lazyness we load lots of repetitions
        //

        mRenderContext->UploadEdgesStart();

        for (auto t : triangles)
        {
            mRenderContext->UploadEdge(
                vertices.GetPosition(triangles.GetVertexAIndex(t)),
                vertices.GetPosition(triangles.GetVertexBIndex(t)),
                edges.GetRenderColor(triangles.GetSubEdgeAIndex(t)));

            mRenderContext->UploadEdge(
                vertices.GetPosition(triangles.GetVertexBIndex(t)),
                vertices.GetPosition(triangles.GetVertexCIndex(t)),
                edges.GetRenderColor(triangles.GetSubEdgeBIndex(t)));

            mRenderContext->UploadEdge(
                vertices.GetPosition(triangles.GetVertexCIndex(t)),
                vertices.GetPosition(triangles.GetVertexAIndex(t)),
                edges.GetRenderColor(triangles.GetSubEdgeCIndex(t)));
        }

        mRenderContext->UploadEdgesEnd();

        //
        // Triangles
        //

        mRenderContext->UploadTrianglesStart();

        for (auto t : triangles)
        {
            std::optional<rgbaColor> color;

            if (mCurrentlySelectedParticleProbe.has_value()
                && mModel->GetParticles().GetState(*mCurrentlySelectedParticleProbe).ConstrainedState.has_value()
                && t == mModel->GetParticles().GetState(*mCurrentlySelectedParticleProbe).ConstrainedState->CurrentTriangle)
            {
                color = rgbaColor(107, 227, 107, 77);
            }
            else if (t == mCurrentOriginTriangle)
            {
                color = rgbaColor(227, 107, 107, 77);
            }
            else if (
                edges.GetSurfaceType(triangles.GetSubEdgeAIndex(t)) == SurfaceType::Floor
                && edges.GetSurfaceType(triangles.GetSubEdgeBIndex(t)) == SurfaceType::Floor
                && edges.GetSurfaceType(triangles.GetSubEdgeCIndex(t)) == SurfaceType::Floor)
            {
                color = rgbaColor(35, 35, 35, 77);
            }

            if (color)
            {
                mRenderContext->UploadTriangle(
                    vertices.GetPosition(triangles.GetVertexAIndex(t)),
                    vertices.GetPosition(triangles.GetVertexBIndex(t)),
                    vertices.GetPosition(triangles.GetVertexCIndex(t)),
                    *color);
            }
        }

        mRenderContext->UploadTrianglesEnd();

        //
        // Particles
        //

        mRenderContext->UploadParticlesStart();

        for (auto p : mModel->GetParticles())
        {
            mRenderContext->UploadParticle(
                mModel->GetParticles().GetPosition(p),
                mModel->GetParticles().GetRenderColor(p),
                1.0f);

            if (mModel->GetParticles().GetState(p).TrajectoryState.has_value())
            {
                mRenderContext->UploadParticle(
                    mModel->GetParticles().GetState(p).TrajectoryState->CurrentPosition,
                    mModel->GetParticles().GetRenderColor(p),
                    0.5f);
            }
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
            vec2f sourcePosition;
            if (mModel->GetParticles().GetState(mCurrentParticleTrajectory->ParticleIndex).TrajectoryState.has_value())
            {
                sourcePosition = mModel->GetParticles().GetState(mCurrentParticleTrajectory->ParticleIndex).TrajectoryState->CurrentPosition;
            }
            else
            {
                sourcePosition = mModel->GetParticles().GetPosition(mCurrentParticleTrajectory->ParticleIndex);
            }

            mRenderContext->UploadParticleTrajectory(
                sourcePosition,
                mCurrentParticleTrajectory->TargetPosition,
                rgbaColor(0x99, 0x99, 0x99, 0xff));
        }

        mRenderContext->UploadParticleTrajectoriesEnd();

        //
        // Mesh velocity
        //

        mRenderContext->UploadMeshVelocity(
            mCurrentMeshTranslationVelocity,
            mCurrentMeshTranslationAccelerationIndicator);
    }

    mRenderContext->RenderEnd();    
}

void LabController::Reset()
{
    assert(mCurrentMeshFilePath);
    LoadMesh(*mCurrentMeshFilePath);
}

void LabController::UpdateMeshTransformations()
{
    LogMessage("UpdateMeshTransformations()");

    vec2f const translation = mCurrentMeshTranslationVelocity * LabParameters::SimulationTimeStepDuration;

    // Update mesh
    auto & vertices = mModel->GetMesh().GetVertices();
    for (auto v : vertices)
    {
        vertices.SetPosition(v, vertices.GetPosition(v) + translation);
    }

    // Update pan
    mRenderContext->SetCameraWorldPosition(mRenderContext->GetCameraWorldPosition() + translation);
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

void LabController::RotateMeshBy(
    vec2f const & centerScreenCoordinates, 
    float screenAngle)
{
    assert(!!mModel);

    vec2f const worldCenter = ScreenToWorld(centerScreenCoordinates);
    float const worldAngle = ScreenOffsetToWorldOffset(screenAngle);

    float const cosAngle = std::cos(worldAngle);
    float const sinAngle = std::sin(worldAngle);

    auto & vertices = mModel->GetMesh().GetVertices();
    for (auto v : vertices)
    {
        vec2f const centeredPos = vertices.GetPosition(v) - worldCenter;
        vec2f const rotatedPos = vec2f(
            centeredPos.x * cosAngle - centeredPos.y * sinAngle,
            centeredPos.x * sinAngle + centeredPos.y * cosAngle);

        vertices.SetPosition(v, rotatedPos + worldCenter);
    }

    assert(mModel->GetParticles().GetElementCount() >= 1);

    // Reset regime
    InitializeParticleRegime(0);

    // Reset trajectory state
    mModel->GetParticles().GetState(0).TrajectoryState.reset();
}

void LabController::RotateMeshBy(
    ElementIndex particleIndex, 
    float screenAngle)
{
    assert(!!mModel);

    assert(particleIndex < mModel->GetParticles().GetElementCount());

    //
    // Rotate mesh
    //

    vec2f const worldCenter = mModel->GetParticles().GetPosition(particleIndex);
    float const worldAngle = ScreenOffsetToWorldOffset(screenAngle);

    float const cosAngle = std::cos(worldAngle);
    float const sinAngle = std::sin(worldAngle);

    auto & vertices = mModel->GetMesh().GetVertices();
    for (auto v : vertices)
    {
        vec2f const centeredPos = vertices.GetPosition(v) - worldCenter;
        vec2f const rotatedPos = vec2f(
            centeredPos.x * cosAngle - centeredPos.y * sinAngle,
            centeredPos.x * sinAngle + centeredPos.y * cosAngle);

        vertices.SetPosition(v, rotatedPos + worldCenter);
    }

    assert(mModel->GetParticles().GetElementCount() >= 1);

    //
    // Make sure that on-edgeness of particle, if any, is maintained
    //

    if (mModel->GetParticles().GetState(particleIndex).ConstrainedState.has_value())
    {
        ElementIndex const triangleIndex = mModel->GetParticles().GetState(particleIndex).ConstrainedState->CurrentTriangle;

        vec3f const oldBarycentricCoords = mModel->GetParticles().GetState(particleIndex).ConstrainedState->CurrentTriangleBarycentricCoords;
        int edgeOrdinal = -1;
        if (oldBarycentricCoords[0] == 0.0f)
        {
            edgeOrdinal = 0;
        }
        else if (oldBarycentricCoords[1] == 0.0f)
        {
            edgeOrdinal = 1;
        }
        else if (oldBarycentricCoords[2] == 0.0f)
        {
            edgeOrdinal = 2;
        }

        if (edgeOrdinal >= 0)
        {
            vec3f newBarycentricCoords = mModel->GetMesh().GetTriangles().ToBarycentricCoordinates(
                mModel->GetParticles().GetPosition(particleIndex),
                triangleIndex,
                mModel->GetMesh().GetVertices());

            newBarycentricCoords[edgeOrdinal] = 0.0f;
            newBarycentricCoords[(edgeOrdinal + 1) % 3] = 1.0f - newBarycentricCoords[(edgeOrdinal + 2) % 3];

            mModel->GetParticles().GetState(particleIndex).ConstrainedState->CurrentTriangleBarycentricCoords = newBarycentricCoords;
        }
    }
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

    // Reset trajectory state
    mModel->GetParticles().GetState(particleIndex).TrajectoryState.reset();
    mCurrentParticleTrajectory.reset();
    mCurrentParticleTrajectoryNotification.reset();
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
    mModel->GetParticles().GetState(particleIndex).TrajectoryState.reset();
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

vec2f const & LabController::GetMeshVelocity() const
{
    return mCurrentMeshTranslationVelocity;
}

void LabController::SetMeshVelocity(vec2f const & velocity)
{
    mCurrentMeshTranslationVelocity = velocity;
    mCurrentMeshTranslationAccelerationIndicator = 1.0f;
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
        if (mModel->GetMesh().GetTriangles().ContainsPoint(position, t, mModel->GetMesh().GetVertices()))
        {
            return t;
        }
    }

    return NoneElementIndex;
}