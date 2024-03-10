/***************************************************************************************
* Original Author:      Gabriele Giuseppini
* Created:              2020-05-16
* Copyright:            Gabriele Giuseppini  (https://github.com/GabrieleGiuseppini)
***************************************************************************************/
#include "LabController.h"

#include "ShipBuilder.h"

#include <GameCore/AABB.h>
#include <GameCore/GameRandomEngine.h>
#include <GameCore/Log.h>

#include <cmath>
#include <sstream>

std::unique_ptr<LabController> LabController::Create(
    int initialCanvasWidth,
    int initialCanvasHeight)
{
    LogMessage("InitialCanvasSize: ", initialCanvasWidth, "x", initialCanvasHeight);

    // Load materials
    MaterialDatabase structuralMaterialDatabase = MaterialDatabase::Load();

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
    MaterialDatabase && materialDatabase,
    std::unique_ptr<RenderContext> renderContext)
    : mMaterialDatabase(std::move(materialDatabase))
    , mRenderContext(std::move(renderContext))
    , mGameEventDispatcher()
    // Simulation state
    , mGameParameters()
    , mModel()
    , mWorld()
    , mCurrentShipFilePath()
    , mCurrentSimulationTime(0.0f)
    , mIsGravityEnabled(true)
    , mCurrentShipTranslationVelocity(vec2f::zero())
    , mCurrentShipTranslationAccelerationIndicator(0.0f)
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

void LabController::LoadShip(std::filesystem::path const & shipDefinitionFilepath)
{
    LoadShip(shipDefinitionFilepath, true);
}

void LabController::Update()
{
    assert(mModel);

    if (mSimulationControlState == SimulationControlStateType::Play
        || mSimulationControlImpulse)
    {
        UpdateShipTransformations();

        mModel->GetNpcs().Update(
            mCurrentSimulationTime,
            mModel->GetShip(),
            mGameParameters);

        // Update state
        mSimulationControlImpulse = false;

        // Update rendering
        mCurrentShipTranslationAccelerationIndicator *= 0.98f;
    }

    mCurrentSimulationTime += GameParameters::SimulationTimeStepDuration;
}

void LabController::Render()
{
    assert(!!mRenderContext);

    mRenderContext->RenderStart();

    mRenderContext->UploadSeaLevel(mWorld.GetOceanSurface().GetDepth());

    if (mModel)
    {
        auto const & points = mModel->GetShip().GetPoints();
        auto const & springs = mModel->GetShip().GetSprings();
        auto const & triangles = mModel->GetShip().GetTriangles();

        //
        // Vertices
        //

        mRenderContext->UploadVertices(
            points.GetElementCount(),
            points.GetPositionBuffer());

        //
        // Edges (of triangles)
        //
        // Note: for lazyness we load lots of repetitions
        //

        mRenderContext->UploadEdgesStart();

        for (auto t : triangles)
        {
            rgbaColor color = springs.GetRenderColor(triangles.GetSubSpringAIndex(t));
            if (mModel->GetNpcs().IsSpringHostingCurrentlySelectedParticle(triangles.GetSubSpringAIndex(t), mModel->GetShip()))
            {
                color.r = static_cast<rgbaColor::data_type>(std::min(color.r + 0x90, static_cast<int>(rgbaColor::data_type_max)));
            }
            mRenderContext->UploadEdge(
                points.GetPosition(triangles.GetPointAIndex(t)),
                points.GetPosition(triangles.GetPointBIndex(t)),
                color);

            color = springs.GetRenderColor(triangles.GetSubSpringBIndex(t));
            if (mModel->GetNpcs().IsSpringHostingCurrentlySelectedParticle(triangles.GetSubSpringBIndex(t), mModel->GetShip()))
            {
                color.r = static_cast<rgbaColor::data_type>(std::min(color.r + 0x90, static_cast<int>(rgbaColor::data_type_max)));
            }
            mRenderContext->UploadEdge(
                points.GetPosition(triangles.GetPointBIndex(t)),
                points.GetPosition(triangles.GetPointCIndex(t)),
                color);

            color = springs.GetRenderColor(triangles.GetSubSpringCIndex(t));
            if (mModel->GetNpcs().IsSpringHostingCurrentlySelectedParticle(triangles.GetSubSpringCIndex(t), mModel->GetShip()))
            {
                color.r = static_cast<rgbaColor::data_type>(std::min(color.r + 0x90, static_cast<int>(rgbaColor::data_type_max)));
            }
            mRenderContext->UploadEdge(
                points.GetPosition(triangles.GetPointCIndex(t)),
                points.GetPosition(triangles.GetPointAIndex(t)),
                color);
        }

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
                springs.GetNpcSurfaceType(triangles.GetSubSpringAIndex(t)) == NpcSurfaceType::Floor
                && springs.GetNpcSurfaceType(triangles.GetSubSpringBIndex(t)) == NpcSurfaceType::Floor
                && springs.GetNpcSurfaceType(triangles.GetSubSpringCIndex(t)) == NpcSurfaceType::Floor)
            {
                color = rgbaColor(35, 35, 35, 77);
            }

            if (color)
            {
                mRenderContext->UploadTriangle(
                    points.GetPosition(triangles.GetPointAIndex(t)),
                    points.GetPosition(triangles.GetPointBIndex(t)),
                    points.GetPosition(triangles.GetPointCIndex(t)),
                    *color);
            }
        }

        mRenderContext->UploadTrianglesEnd();

        //
        // Npcs
        //

        mModel->GetNpcs().Render(*mRenderContext);

        //
        // Ship velocity
        //

        mRenderContext->UploadShipVelocity(
            mCurrentShipTranslationVelocity,
            mCurrentShipTranslationAccelerationIndicator);
    }

    mRenderContext->RenderEnd();
}

void LabController::Reset()
{
    assert(mCurrentShipFilePath);
    LoadShip(*mCurrentShipFilePath);
}

void LabController::UpdateShipTransformations()
{
    vec2f const translation = mCurrentShipTranslationVelocity * GameParameters::SimulationTimeStepDuration;

    // Update ship
    auto & points = mModel->GetShip().GetPoints();
    for (auto p : points)
    {
        points.SetPosition(p, points.GetPosition(p) + translation);
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

    float constexpr SquareSearchRadius = GameParameters::VertexRadius * GameParameters::VertexRadius;

    float bestSquareDistance = std::numeric_limits<float>::max();
    ElementIndex bestPoint = NoneElementIndex;

    auto const & points = mModel->GetShip().GetPoints();
    for (auto p : points)
    {
        float const squareDistance = (points.GetPosition(p) - worldCoordinates).squareLength();
        if (squareDistance < SquareSearchRadius
            && squareDistance < bestSquareDistance)
        {
            bestSquareDistance = squareDistance;
            bestPoint = p;
        }
    }

    if (bestPoint != NoneElementIndex)
        return bestPoint;
    else
        return std::nullopt;
}

void LabController::MoveVertexBy(
    ElementIndex pointIndex,
    vec2f const & screenOffset)
{
    assert(!!mModel);

    vec2f const worldOffset = ScreenOffsetToWorldOffset(screenOffset);

    mModel->GetShip().GetPoints().SetPosition(
        pointIndex,
        mModel->GetShip().GetPoints().GetPosition(pointIndex) + worldOffset);

    mModel->GetNpcs().OnPointMoved(mCurrentSimulationTime, mModel->GetShip());
}

void LabController::RotateShipBy(
    vec2f const & centerScreenCoordinates,
    float screenStride)
{
    assert(!!mModel);

    vec2f const worldCenter = ScreenToWorld(centerScreenCoordinates);
    float const worldAngle = ScreenOffsetToWorldOffset(screenStride) * 0.05f;

    float const cosAngle = std::cos(worldAngle);
    float const sinAngle = std::sin(worldAngle);

    //
    // Rotate ship
    //

    auto & points = mModel->GetShip().GetPoints();
    for (auto p : points)
    {
        vec2f const centeredPos = points.GetPosition(p) - worldCenter;
        vec2f const rotatedPos = vec2f(
            centeredPos.x * cosAngle - centeredPos.y * sinAngle,
            centeredPos.x * sinAngle + centeredPos.y * cosAngle);

        points.SetPosition(p, rotatedPos + worldCenter);
    }

    // TODOTEST
    //////
    ////// Rotate particles
    //////

    ////mModel->GetNpcs().RotateParticlesWithShip(
    ////    worldCenter,
    ////    cosAngle,
    ////    sinAngle,
    ////    mModel->GetShip());
}

void LabController::RotateShipBy(
    ElementIndex particleIndex,
    float screenStride)
{
    assert(!!mModel);

    //
    // Rotate ship
    //

    vec2f const worldCenter = mModel->GetNpcs().GetParticles().GetPosition(particleIndex);
    float const worldAngle = ScreenOffsetToWorldOffset(screenStride) * 0.05f;

    float const cosAngle = std::cos(worldAngle);
    float const sinAngle = std::sin(worldAngle);

    auto & points = mModel->GetShip().GetPoints();
    for (auto p : points)
    {
        vec2f const centeredPos = points.GetPosition(p) - worldCenter;
        vec2f const rotatedPos = vec2f(
            centeredPos.x * cosAngle - centeredPos.y * sinAngle,
            centeredPos.x * sinAngle + centeredPos.y * cosAngle);

        points.SetPosition(p, rotatedPos + worldCenter);
    }

    // TODOTEST
    //////
    ////// Rotate particles
    //////

    ////mModel->GetNpcs().RotateParticlesWithShip(
    ////    worldCenter,
    ////    cosAngle,
    ////    sinAngle,
    ////    mModel->GetShip());
}

bool LabController::TrySelectOriginTriangle(vec2f const & screenCoordinates)
{
    assert(!!mModel);

    vec2f const worldCoordinates = ScreenToWorld(screenCoordinates);

    ElementIndex const t = mModel->GetShip().GetTriangles().FindContaining(worldCoordinates, mModel->GetShip().GetPoints());
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

    float constexpr SquareSearchRadius = GameParameters::ParticleRadius * GameParameters::ParticleRadius;

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
        mModel->GetNpcs().SelectParticle(bestParticle, mModel->GetShip());
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

    float constexpr SquareSearchRadius = GameParameters::ParticleRadius * GameParameters::ParticleRadius;

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
        mModel->GetShip());
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

vec2f const & LabController::GetShipVelocity() const
{
    return mCurrentShipTranslationVelocity;
}

void LabController::SetShipVelocity(vec2f const & velocity)
{
    mCurrentShipTranslationVelocity = velocity;
    mCurrentShipTranslationAccelerationIndicator = 1.0f;
}

void LabController::DoStepForVideo()
{
    assert(mModel);

    mModel->GetNpcs().FlipHumanWalk(0);
    //mModel->GetNpcs().FlipHumanFrontBack(0);

    ////// TODOTEST
    //////RotateShipBy(0, 3.1f * 3.0f); // Stuck against wall
    ////RotateShipBy(0, -3.1f * 12.0f); // Sfarfallio
    ////return;

    ////++mCurrentVideoStep;

    ////int matchIndex = 0;

    ////////++matchIndex;
    ////////if (mCurrentVideoStep == matchIndex)
    ////////{
    ////////    //
    ////////    // Load ship and stay clean
    ////////    //

    ////////    LoadShip(std::filesystem::absolute("Meshes\\video_mesh.png"), false);

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
    ////////        mModel->GetShip());

    ////////    return;
    ////////}

    ////////++matchIndex;
    ////////if (mCurrentVideoStep == matchIndex)
    ////////{
    ////////    //
    ////////    // Reset
    ////////    //

    ////////    LoadShip(std::filesystem::absolute("Meshes\\video_mesh.png"), false);

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
    ////////            mModel->GetShip());
    ////////    }

    ////////    return;
    ////////}

    ////////++matchIndex;
    ////////if (mCurrentVideoStep == matchIndex)
    ////////{
    ////////    //
    ////////    // Reset
    ////////    //

    ////////    LoadShip(std::filesystem::absolute("Meshes\\video_mesh.png"), false);

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
    ////////            mModel->GetShip());
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
    ////    // Load ship and stay clean
    ////    //

    ////    LoadShip(std::filesystem::absolute("Meshes\\video_mesh.png"), false);

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
    ////    vec2f secondaryPosition = primaryPosition + vec2f(1.0f, 0.0f) * GameParameters::HumanNpcGeometry::BodyLength * mLabParameters.HumanNpcBodyLengthAdjustment;
    ////    mModel->GetNpcs().Add(
    ////        Npcs::NpcType::Human,
    ////        primaryPosition,
    ////        secondaryPosition,
    ////        mCurrentSimulationTime,
    ////        mStructuralMaterialDatabase,
    ////        mModel->GetShip(),
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

std::optional<PickedObjectId<NpcId>> LabController::BeginPlaceNewHumanNpc(HumanNpcKindType humanKind, vec2f const & screenCoordinates)
{
    vec2f const worldCoordinates = ScreenToWorld(screenCoordinates);

    // TODOHERE
    (void)humanKind;
    return std::nullopt;
}

std::optional<PickedObjectId<NpcId>> LabController::ProbeNpc(vec2f const & screenCoordinates) const
{
    vec2f const worldCoordinates = ScreenToWorld(screenCoordinates);

    // TODOHERE
    return std::nullopt;
}

void LabController::MoveNpcTo(NpcId id, vec2f const & screenCoordinates, vec2f const & worldOffset)
{
    vec2f const worldCoordinates = ScreenToWorld(screenCoordinates);

    // TODOHERE
    (void)id;
    (void)worldOffset;
}

void LabController::EndMoveNpc(NpcId id)
{
    // TODOHERE
    (void)id;
}

void LabController::CompleteNewNpc(NpcId id)
{
    // TODOHERE
    (void)id;
}

void LabController::AbortNewNpc(NpcId id)
{
    // TODOHERE
    (void)id;
}

void LabController::RemoveNpc(NpcId id)
{
    // TODOHERE
    (void)id;
}

void LabController::HighlightNpc(NpcId id, NpcHighlightType highlight)
{
    // TODOHERE
    (void)id;
    (void)highlight;
}

void LabController::SetNpcPanicLevelForAllHumans(float panicLevel)
{
    mModel->GetNpcs().SetPanicLevelForAllHumans(panicLevel);
}

////////////////////////////////////////////////

void LabController::LoadShip(
    std::filesystem::path const & shipDefinitionFilepath,
    bool addExperimentalNpc)
{
    // Load ship definition
    auto shipDefinition = ShipDefinition::Load(shipDefinitionFilepath);

    // Make ship

    std::unique_ptr<Physics::Ship> ship = ShipBuilder::BuildShip(
        std::move(shipDefinition),
        mMaterialDatabase);

    if (ship->GetTriangles().GetElementCount() < 1)
    {
        throw GameException("Ship must contain at least one triangle");
    }

    // Create NPCs

    std::unique_ptr<Physics::Npcs> npcs = std::make_unique<Physics::Npcs>(
        mWorld,
        mMaterialDatabase,
        mGameEventDispatcher,
        mGameParameters,
        mIsGravityEnabled);

    if (addExperimentalNpc)
    {
        // TODO: for small ship, in the middle
        //vec2f const position = vec2f(0.5f, 0.0f);
        // TODO: for small ship, on the floor, left triangle
        vec2f const position = vec2f(-0.5f, -2.0f);
        // TODO: for small ship, on the floor, right triangle
        //vec2f const position = vec2f(0.5f, -2.0f);
        // TODO: for large ship, on floor
        //vec2f const position = vec2f(5.5f, -6.0f);

        // TODO: for repro of traj acceleration w/human
        //vec2f const position = vec2f(-0.634f, -2.0f);

        ////// TODO: for repro of traj acceleration w/ball
        ////ElementIndex const triangleIndex = 46;
        ////float const TODO = 0.99435f;
        ////bcoords3f const baryCoords = bcoords3f(TODO, 0.0f, 1.0f - TODO);
        ////vec2f const position = ship->GetTriangles().FromBarycentricCoordinates(baryCoords, triangleIndex, ship->GetPoints());

        npcs->Add(
            // TODOTEST
            //NpcKindType::Furniture,
            NpcKindType::Human,
            position,
            std::nullopt, // Secondary position
            mCurrentSimulationTime,
            *ship,
            mGameParameters);

        ////// TODOTEST: multiple balls
        ////for (int i = 0; i < 40; ++i)
        ////{
        ////    vec2f const p = vec2f(
        ////        GameRandomEngine::GetInstance().GenerateUniformReal(-9.0f, 8.0f),
        ////        GameRandomEngine::GetInstance().GenerateUniformReal(-5.0f, 5.0f));

        ////    npcs->Add(
        ////        Npcs::NpcType::Furniture,
        ////        p,
        ////        std::nullopt, // Secondary position
        ////        mCurrentSimulationTime,
        ////        mStructuralMaterialDatabase,
        ////        *ship,
        ////        mLabParameters);
        ////}

        ////// TODO: for repro w/ball, part II
        ////assert(npcs->GetState(0).PrimaryParticleState.ConstrainedState.has_value());
        ////assert(npcs->GetState(0).PrimaryParticleState.ConstrainedState->CurrentTriangle == triangleIndex);
        ////npcs->GetState(0).PrimaryParticleState.ConstrainedState->CurrentTriangleBarycentricCoords = baryCoords;
        ////npcs->GetParticles().SetVelocity(npcs->GetState(0).PrimaryParticleState.ParticleIndex, vec2f(-1.0f, 0.0f));
        ////// TODOTEST: trying now for inertial (no G and no friction)
        //////npcs->GetParticles().SetVelocity(npcs->GetState(0).PrimaryParticleState.ParticleIndex, vec2f(-1.0f, (GameParameters::GravityMagnitude + 16.25f) * GameParameters::SimulationTimeStepDuration));

        // Select particle
        npcs->SelectParticle(0, *ship);
    }

    // Create a new model
    std::unique_ptr<Model> newModel = std::make_unique<Model>(
        std::move(ship),
        std::move(npcs),
        mGameEventDispatcher);

    //
    // No errors, so we may continue
    //

    Reset(
        std::move(newModel),
        shipDefinitionFilepath);
}

void LabController::Reset(
    std::unique_ptr<Model> newModel,
    std::filesystem::path const & shipDefinitionFilepath)
{
    //
    // Take object in
    //

    mModel.reset();
    mModel = std::move(newModel);
    mCurrentShipFilePath = shipDefinitionFilepath;

    //
    // Reset our state
    //

    mCurrentSimulationTime = 0.0f;

    //
    // Auto-zoom & center
    //

    {
        AABB const objectAABB = mModel->GetShip().GetPoints().GetAABB();

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

    mGameEventDispatcher.OnBLabReset();
}

