/***************************************************************************************
* Original Author:      Gabriele Giuseppini
* Created:              2020-05-16
* Copyright:            Gabriele Giuseppini  (https://github.com/GabrieleGiuseppini)
***************************************************************************************/
#include "LabController.h"

#include "AABB.h"
#include "GameRandomEngine.h"
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
    StructuralMaterialDatabase && structuralMaterialDatabase,
    std::unique_ptr<RenderContext> renderContext)
    : mStructuralMaterialDatabase(std::move(structuralMaterialDatabase))
    , mRenderContext(std::move(renderContext))
    , mEventDispatcher()
    // Simulation state
    , mLabParameters()
    , mModel()
    , mWorld()
    , mCurrentMeshFilePath()
    , mCurrentSimulationTime(0.0f)
    , mIsGravityEnabled(true)
    , mCurrentMeshTranslationVelocity(vec2f::zero())
    , mCurrentMeshTranslationAccelerationIndicator(0.0f)
    , mCurrentVideoStep(0)
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
    LoadMesh(meshDefinitionFilepath, true);
}

void LabController::Update()
{
    assert(mModel);

    if (mSimulationControlState == SimulationControlStateType::Play
        || mSimulationControlImpulse)
    {
        UpdateMeshTransformations();            

        mModel->GetNpcs().Update(
            mCurrentSimulationTime,
            mModel->GetMesh(),
            mLabParameters);        

        // Update state
        mSimulationControlImpulse = false;

        // Update rendering
        mCurrentMeshTranslationAccelerationIndicator *= 0.98f;
    }

    mCurrentSimulationTime += LabParameters::SimulationTimeStepDuration;
}

void LabController::Render()
{
    assert(!!mRenderContext);

    mRenderContext->RenderStart();

    mRenderContext->UploadSeaLevel(mWorld.GetOceanSurface().GetDepth());

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
            rgbaColor color = edges.GetRenderColor(triangles.GetSubEdgeAIndex(t));
            if (mModel->GetNpcs().IsEdgeHostingCurrentlySelectedParticle(triangles.GetSubEdgeAIndex(t)))
            {
                color.r = static_cast<rgbaColor::data_type>(std::min(color.r + 0x90, static_cast<int>(rgbaColor::data_type_max)));
            }
            mRenderContext->UploadEdge(
                vertices.GetPosition(triangles.GetVertexAIndex(t)),
                vertices.GetPosition(triangles.GetVertexBIndex(t)),
                color);

            color = edges.GetRenderColor(triangles.GetSubEdgeBIndex(t));
            if (mModel->GetNpcs().IsEdgeHostingCurrentlySelectedParticle(triangles.GetSubEdgeBIndex(t)))
            {
                color.r = static_cast<rgbaColor::data_type>(std::min(color.r + 0x90, static_cast<int>(rgbaColor::data_type_max)));
            }
            mRenderContext->UploadEdge(
                vertices.GetPosition(triangles.GetVertexBIndex(t)),
                vertices.GetPosition(triangles.GetVertexCIndex(t)),
                color);

            color = edges.GetRenderColor(triangles.GetSubEdgeCIndex(t));
            if (mModel->GetNpcs().IsEdgeHostingCurrentlySelectedParticle(triangles.GetSubEdgeCIndex(t)))
            {
                color.r = static_cast<rgbaColor::data_type>(std::min(color.r + 0x90, static_cast<int>(rgbaColor::data_type_max)));
            }
            mRenderContext->UploadEdge(
                vertices.GetPosition(triangles.GetVertexCIndex(t)),
                vertices.GetPosition(triangles.GetVertexAIndex(t)),
                color);
        }

        ////for (auto e : edges)
        ////{
        ////    rgbaColor color = edges.GetRenderColor(e);
        ////    if (mModel->GetNpcs().IsEdgeHostingCurrentlySelectedParticle(e))
        ////    {
        ////        color.r = static_cast<rgbaColor::data_type>(std::min(color.r + 0x90, static_cast<int>(rgbaColor::data_type_max)));
        ////    }

        ////    mRenderContext->UploadEdge(
        ////        vertices.GetPosition(edges.GetEndpointAIndex(e)),
        ////        vertices.GetPosition(edges.GetEndpointBIndex(e)),
        ////        color);
        ////}

        mRenderContext->UploadEdgesEnd();

        //
        // Triangles
        //

        mRenderContext->UploadTrianglesStart();

        for (auto t : triangles)
        {
            std::optional<rgbaColor> color;

            if (mModel->GetNpcs().IsTriangleConstrainingCurrentlySelectedParticle(t))
            {
                color = rgbaColor(107, 227, 107, 77);
            }
            else if (t == mModel->GetNpcs().GetCurrentOriginTriangle())
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
        // Npcs
        //

        mModel->GetNpcs().Render(*mRenderContext);

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
    vec2f const translation = mCurrentMeshTranslationVelocity * LabParameters::SimulationTimeStepDuration;

    // Update mesh
    auto & vertices = mModel->GetMesh().GetVertices();
    for (auto v : vertices)
    {
        vertices.SetPosition(v, vertices.GetPosition(v) + translation);
    }

    // Update camera pan
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

    mModel->GetNpcs().OnVertexMoved(mCurrentSimulationTime, mModel->GetMesh());
}

void LabController::RotateMeshBy(
    vec2f const & centerScreenCoordinates, 
    float screenStride)
{
    assert(!!mModel);

    vec2f const worldCenter = ScreenToWorld(centerScreenCoordinates);
    float const worldAngle = ScreenOffsetToWorldOffset(screenStride) * 0.05f;

    float const cosAngle = std::cos(worldAngle);
    float const sinAngle = std::sin(worldAngle);

    //
    // Rotate mesh
    //

    auto & vertices = mModel->GetMesh().GetVertices();
    for (auto v : vertices)
    {
        vec2f const centeredPos = vertices.GetPosition(v) - worldCenter;
        vec2f const rotatedPos = vec2f(
            centeredPos.x * cosAngle - centeredPos.y * sinAngle,
            centeredPos.x * sinAngle + centeredPos.y * cosAngle);

        vertices.SetPosition(v, rotatedPos + worldCenter);
    }

    //
    // Rotate particles
    //

    mModel->GetNpcs().RotateParticlesWithMesh(
        worldCenter,
        cosAngle,
        sinAngle,
        mModel->GetMesh());
}

void LabController::RotateMeshBy(
    ElementIndex particleIndex, 
    float screenStride)
{
    assert(!!mModel);

    assert(particleIndex < mModel->GetNpcs().GetParticles().GetParticleCount());

    //
    // Rotate mesh
    //

    vec2f const worldCenter = mModel->GetNpcs().GetParticles().GetPosition(particleIndex);
    float const worldAngle = ScreenOffsetToWorldOffset(screenStride) * 0.05f;

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

    //
    // Rotate particles
    //

    mModel->GetNpcs().RotateParticlesWithMesh(
        worldCenter,
        cosAngle,
        sinAngle,
        mModel->GetMesh());
}

bool LabController::TrySelectOriginTriangle(vec2f const & screenCoordinates)
{
    assert(!!mModel);

    vec2f const worldCoordinates = ScreenToWorld(screenCoordinates);

    ElementIndex const t = mModel->GetMesh().GetTriangles().FindContaining(worldCoordinates, mModel->GetMesh().GetVertices());
    if (t != NoneElementIndex)
    {
        mModel->GetNpcs().SelectOriginTriangle(t);
        return true;
    }
    else
    {
        mModel->GetNpcs().ResetOriginTriangle();
        return false;
    }
}

bool LabController::TrySelectParticle(vec2f const & screenCoordinates)
{
    assert(!!mModel);

    //
    // Find closest particle within the radius
    //

    vec2f const worldCoordinates = ScreenToWorld(screenCoordinates);

    float constexpr SquareSearchRadius = LabParameters::ParticleRadius * LabParameters::ParticleRadius;

    float bestSquareDistance = std::numeric_limits<float>::max();
    ElementIndex bestParticle = NoneElementIndex;

    auto const & particles = mModel->GetNpcs().GetParticles();
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
    {
        mModel->GetNpcs().SelectParticle(bestParticle, mModel->GetMesh());
        return true;
    }
    else
    {
        return false;
    }
}

std::optional<ElementIndex> LabController::TryPickNpcParticle(vec2f const & screenCoordinates) const
{
    assert(!!mModel);

    //
    // Find closest particle within the radius
    //

    vec2f const worldCoordinates = ScreenToWorld(screenCoordinates);

    float constexpr SquareSearchRadius = LabParameters::ParticleRadius * LabParameters::ParticleRadius;

    float bestSquareDistance = std::numeric_limits<float>::max();
    ElementIndex bestParticle = NoneElementIndex;

    auto const & particles = mModel->GetNpcs().GetParticles();
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

void LabController::MoveNpcParticleBy(
    ElementIndex npcParticleIndex, 
    vec2f const & screenOffset,
    vec2f const & inertialStride)
{
    (void)inertialStride;

    assert(!!mModel);

    vec2f const worldOffset = ScreenOffsetToWorldOffset(screenOffset);

    mModel->GetNpcs().MoveParticleBy(
        npcParticleIndex,
        worldOffset,
        mCurrentSimulationTime,
        mModel->GetMesh());
}

void LabController::NotifyNpcParticleTrajectory(
    ElementIndex npcParticleIndex,
    vec2f const & targetScreenCoordinates)
{
    assert(!!mModel);

    mModel->GetNpcs().NotifyParticleTrajectory(
        npcParticleIndex,
        ScreenToWorld(targetScreenCoordinates));
}

void LabController::SetNpcParticleTrajectory(
    ElementIndex npcParticleIndex,
    vec2f const & targetScreenCoordinates)
{
    mModel->GetNpcs().SetParticleTrajectory(
        npcParticleIndex,
        ScreenToWorld(targetScreenCoordinates));
}

void LabController::QueryNearestNpcParticleAt(vec2f const & screenCoordinates) const
{
    assert(!!mModel);

    auto const nearestNpcParticle = TryPickNpcParticle(screenCoordinates);
    if (nearestNpcParticle.has_value())
    {
        mModel->GetNpcs().GetParticles().Query(*nearestNpcParticle);
    }
}

bool LabController::IsGravityEnabled() const
{
    return mIsGravityEnabled;
}

void LabController::SetGravityEnabled(bool isEnabled)
{
    mIsGravityEnabled = isEnabled;

    if (mModel)
    {
        mModel->GetNpcs().SetGravityEnabled(isEnabled);
    }
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

void LabController::DoStepForVideo()
{
    assert(mModel);

    // TODOTEST
    RotateMeshBy(0, -0.5f);
    return;

    ////++mCurrentVideoStep;

    ////int matchIndex = 0;

    ////////++matchIndex;
    ////////if (mCurrentVideoStep == matchIndex)
    ////////{
    ////////    //
    ////////    // Load mesh and stay clean
    ////////    //

    ////////    LoadMesh(std::filesystem::absolute("Meshes\\video_mesh.png"), false);

    ////////    // Enable gravity
    ////////    SetGravityEnabled(true);

    ////////    // Enable auto-play            
    ////////    SetSimulationControlState(SimulationControlStateType::Play);

    ////////    // Other settings
    ////////    SetSeaLevel(-7.0f);        

    ////////    return;
    ////////}

    ////////++matchIndex;
    ////////if (mCurrentVideoStep == matchIndex)
    ////////{
    ////////    //
    ////////    // Create one ball
    ////////    //

    ////////    mModel->GetNpcs().Add(
    ////////        Npcs::NpcType::Furniture,
    ////////        vec2f(0.0f, 0.0f),
    ////////        std::nullopt,
    ////////        mStructuralMaterialDatabase,
    ////////        mModel->GetMesh());

    ////////    return;
    ////////}

    ////////++matchIndex;
    ////////if (mCurrentVideoStep == matchIndex)
    ////////{
    ////////    //
    ////////    // Reset
    ////////    //

    ////////    LoadMesh(std::filesystem::absolute("Meshes\\video_mesh.png"), false);

    ////////    return;
    ////////}

    ////////++matchIndex;
    ////////if (mCurrentVideoStep == matchIndex)
    ////////{
    ////////    //
    ////////    // Create a few balls
    ////////    //

    ////////    for (int i = 0; i < 9; ++i)
    ////////    {
    ////////        vec2f const position = vec2f(
    ////////            GameRandomEngine::GetInstance().GenerateUniformReal(-9.0f, 8.0f),
    ////////            GameRandomEngine::GetInstance().GenerateUniformReal(-5.0f, 5.0f));

    ////////        mModel->GetNpcs().Add(
    ////////            Npcs::NpcType::Furniture,
    ////////            position,
    ////////            std::nullopt,
    ////////            mStructuralMaterialDatabase,
    ////////            mModel->GetMesh());
    ////////    }

    ////////    return;
    ////////}

    ////////++matchIndex;
    ////////if (mCurrentVideoStep == matchIndex)
    ////////{
    ////////    //
    ////////    // Reset
    ////////    //

    ////////    LoadMesh(std::filesystem::absolute("Meshes\\video_mesh.png"), false);

    ////////    return;
    ////////}

    ////////++matchIndex;
    ////////if (mCurrentVideoStep == matchIndex)
    ////////{
    ////////    //
    ////////    // Create many balls
    ////////    //

    ////////    for (int i = 0; i < 25; ++i)
    ////////    {
    ////////        vec2f const position = vec2f(
    ////////            GameRandomEngine::GetInstance().GenerateUniformReal(-9.0f, 8.0f),
    ////////            GameRandomEngine::GetInstance().GenerateUniformReal(-5.0f, 5.0f));

    ////////        mModel->GetNpcs().Add(
    ////////            Npcs::NpcType::Furniture,
    ////////            position,
    ////////            std::nullopt,
    ////////            mStructuralMaterialDatabase,
    ////////            mModel->GetMesh());
    ////////    }

    ////////    return;
    ////////}

    /////////////////////////////////////////////////////////////
    ////// Humans
    /////////////////////////////////////////////////////////////

    ////++matchIndex;
    ////if (mCurrentVideoStep == matchIndex)
    ////{
    ////    //
    ////    // Load mesh and stay clean
    ////    //

    ////    LoadMesh(std::filesystem::absolute("Meshes\\video_mesh.png"), false);

    ////    // Enable auto-play            
    ////    SetSimulationControlState(SimulationControlStateType::Play);

    ////    // Disable gravity
    ////    SetGravityEnabled(false);

    ////    // Other settings
    ////    SetHumanNpcEquilibriumTorqueStiffnessCoefficient(0.0f);
    ////    SetSeaLevel(-7.0f);        

    ////    return;
    ////}

    ////++matchIndex;
    ////if (mCurrentVideoStep == matchIndex)
    ////{
    ////    //
    ////    // Add one horizontal human
    ////    //

    ////    vec2f primaryPosition = vec2f(1.5f, -2.0f);
    ////    vec2f secondaryPosition = primaryPosition + vec2f(1.0f, 0.0f) * LabParameters::HumanNpcGeometry::BodyLength * mLabParameters.HumanNpcBodyLengthAdjustment;
    ////    mModel->GetNpcs().Add(
    ////        Npcs::NpcType::Human,
    ////        primaryPosition,
    ////        secondaryPosition,
    ////        mCurrentSimulationTime,
    ////        mStructuralMaterialDatabase,
    ////        mModel->GetMesh(),
    ////        mLabParameters);

    ////    return;
    ////}

    ////++matchIndex;
    ////if (mCurrentVideoStep == matchIndex)
    ////{
    ////    //
    ////    // Enable gravity
    ////    //


    ////    SetGravityEnabled(true);

    ////    return;
    ////}

    ////++matchIndex;
    ////if (mCurrentVideoStep == matchIndex)
    ////{
    ////    //
    ////    // Low-rising 1
    ////    //

    ////    SetHumanNpcEquilibriumTorqueStiffnessCoefficient(0.0010f);

    ////    return;
    ////}

    ////++matchIndex;
    ////if (mCurrentVideoStep == matchIndex)
    ////{
    ////    //
    ////    // Stop rising
    ////    //

    ////    SetHumanNpcEquilibriumTorqueStiffnessCoefficient(0.0f);

    ////    return;
    ////}

    ////++matchIndex;
    ////if (mCurrentVideoStep == matchIndex)
    ////{
    ////    //
    ////    // Low-rising 2
    ////    //

    ////    SetHumanNpcEquilibriumTorqueStiffnessCoefficient(0.0013f);

    ////    return;
    ////}

    ////++matchIndex;
    ////if (mCurrentVideoStep == matchIndex)
    ////{
    ////    //
    ////    // Stop rising
    ////    //

    ////    SetHumanNpcEquilibriumTorqueStiffnessCoefficient(0.0f);

    ////    return;
    ////}

    ////++matchIndex;
    ////if (mCurrentVideoStep == matchIndex)
    ////{
    ////    //
    ////    // Low-rising 3
    ////    //

    ////    SetHumanNpcEquilibriumTorqueStiffnessCoefficient(0.0018f);

    ////    return;
    ////}

    ////++matchIndex;
    ////if (mCurrentVideoStep == matchIndex)
    ////{
    ////    //
    ////    // Stop rising
    ////    //

    ////    SetHumanNpcEquilibriumTorqueStiffnessCoefficient(0.0f);

    ////    return;
    ////}

    ////++matchIndex;
    ////if (mCurrentVideoStep == matchIndex)
    ////{
    ////    //
    ////    // Definitive rising  - but no walking
    ////    //

    ////    SetHumanNpcEquilibriumTorqueStiffnessCoefficient(0.0032f);

    ////    // But no walking
    ////    mLabParameters.HumanNpcWalkingAcceleration = 0.0f;

    ////    return;
    ////}

    ////++matchIndex;
    ////if (mCurrentVideoStep == matchIndex)
    ////{
    ////    //
    ////    // Walking
    ////    //

    ////    mLabParameters.HumanNpcWalkingAcceleration = 0.027f;

    ////    return;
    ////}

    ////++matchIndex;
    ////if (mCurrentVideoStep == matchIndex)
    ////{
    ////    return;
    ////}

    //////
    ////// Wrap around
    //////

    ////mCurrentVideoStep = 0;
}

////////////////////////////////////////////////

void LabController::LoadMesh(
    std::filesystem::path const & meshDefinitionFilepath,
    bool addExperimentalNpc)
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

    // Create NPCs

    std::unique_ptr<Npcs> npcs = std::make_unique<Npcs>(
        mWorld, 
        mEventDispatcher, 
        mLabParameters,
        mIsGravityEnabled);

    if (addExperimentalNpc)
    {
        // TODOTEST
        ////vec2f const position = (
        ////    mesh->GetVertices().GetPosition(mesh->GetTriangles().GetVertexAIndex(0))
        ////    + mesh->GetVertices().GetPosition(mesh->GetTriangles().GetVertexBIndex(0))
        ////    + mesh->GetVertices().GetPosition(mesh->GetTriangles().GetVertexCIndex(0))) / 3.0f;
        //vec2f const position =
        //    mesh->GetVertices().GetPosition(mesh->GetTriangles().GetVertexBIndex(0))
        //    - vec2f(0.5f, 0.5f);

        // TODOTEST
        //vec2f const position = vec2f(0.5f, -2.0f);
        // TODO: for small mesh, in the middle
        //vec2f const position = vec2f(0.5f, 0.0f);
        // TODO: for small mesh, on the floor, left triangle
        //vec2f const position = vec2f(-0.5f, -2.0f);
        // TODO: for small mesh, on the floor, right triangle
        //vec2f const position = vec2f(0.5f, -2.0f);
        // TODO: for large mesh, on floor
        //vec2f const position = vec2f(5.5f, -6.0f);

        // TODO: for repro
        //vec2f const position = vec2f(-0.5f, -2.0f);
        //vec2f const position = vec2f(-0.6395f, -2.0f);
        vec2f const position = vec2f(-0.634f, -2.0f);

        npcs->Add(
            // TODOTEST
            Npcs::NpcType::Furniture,
            //Npcs::NpcType::Human,
            position,
            std::nullopt,
            mCurrentSimulationTime,
            mStructuralMaterialDatabase,
            *mesh,
            mLabParameters);

        // Select particle
        assert(npcs->GetParticles().GetParticleCount() > 0);
        npcs->SelectParticle(0, *mesh);
    }

    // Create a new model
    std::unique_ptr<Model> newModel = std::make_unique<Model>(
        std::move(mesh),
        std::move(npcs),
        mEventDispatcher);

    //
    // No errors, so we may continue
    //

    Reset(
        std::move(newModel),
        meshDefinitionFilepath);
}

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
    // Reset our state
    //

    mCurrentSimulationTime = 0.0f;

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

