/***************************************************************************************
* Original Author:      Gabriele Giuseppini
* Created:              2020-05-16
* Copyright:            Gabriele Giuseppini  (https://github.com/GabrieleGiuseppini)
***************************************************************************************/
#include "LabController.h"

#include "ShipFactory.h"

#include <GameCore/AABB.h>
#include <GameCore/GameRandomEngine.h>
#include <GameCore/Log.h>

#include <chrono>
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
    std::unique_ptr<Render::RenderContext> renderContext = std::make_unique<Render::RenderContext>(
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
    std::unique_ptr<Render::RenderContext> renderContext)
    : mMaterialDatabase(std::move(materialDatabase))
    , mRenderContext(std::move(renderContext))
    , mGameEventHandler(std::make_shared<GameEventDispatcher>())
    // Simulation state
    , mGameParameters()
    , mWorld()
    , mCurrentShipFilePath()
    , mCurrentSimulationTime(0.0f)
    , mPerfStats()
    , mLastPublishedPerfStats()
    , mLastPerfPublishTimestamp(GameChronometer::now())
    , mMassAdjustment(1.0f)
    , mGravityAdjustment(1.0f)
    , mOceanDepth(-7.0f)
    , mCurrentShipTranslationVelocity(vec2f::zero())
    , mCurrentShipTranslationAccelerationIndicator(0.0f)
    , mCurrentWavesAmplitude(0.0f)
    , mTargetWavesAmplitude(0.0f)
    , mCurrentWavesSpeed(0.0f)
    , mTargetWavesSpeed(0.0f)
    , mLastWaveRotationAngle(0.0f)
    , mLastWaveTimeArg(0.0f)
    , mCurrentVideoStep(0)
    // Simulation control
    , mSimulationControlState(SimulationControlStateType::Play)
    , mSimulationControlImpulse(false)
{
}

SimulationControlStateType LabController::GetSimulationControlState() const
{
    return mSimulationControlState;
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
    // Load ship definition
    auto shipDefinition = ShipDefinition::Load(shipDefinitionFilepath);

    // Make ship

    std::unique_ptr<Physics::Ship> ship = ShipFactory::BuildShip(
        std::move(shipDefinition),
        mMaterialDatabase);

    if (ship->GetTriangles().GetElementCount() < 1)
    {
        throw GameException("Ship must contain at least one triangle");
    }

    //
    // No errors, so we may continue
    //

    Reset(
        std::move(ship),
        shipDefinitionFilepath);
}

void LabController::Update()
{
    assert(mWorld);

    if (mSimulationControlState == SimulationControlStateType::Play
        || mSimulationControlImpulse)
    {
        UpdateShipTransformations();

        mWorld->Update(
            mCurrentSimulationTime,
            mGameParameters,
            mPerfStats);

#ifdef _DEBUG
        mWorld->GetNpcs().Publish();
#endif

        // Update state
        mSimulationControlImpulse = false;

        // Update rendering
        mCurrentShipTranslationAccelerationIndicator *= 0.98f;

        mCurrentSimulationTime += GameParameters::SimulationStepTimeDuration<float>;

        auto const now = GameChronometer::now();
        if (now > mLastPerfPublishTimestamp + std::chrono::milliseconds(500))
        {
            // Publish perf stats
            auto const deltaStats = mPerfStats - mLastPublishedPerfStats;
            mGameEventHandler->OnUpdateTimeMeasured(
                deltaStats.TotalNpcUpdateDuration.ToRatio<std::chrono::milliseconds>(),
                deltaStats.TotalNpcRenderUploadDuration.ToRatio<std::chrono::milliseconds>());

            mLastPublishedPerfStats = mPerfStats;
            mLastPerfPublishTimestamp = now;
        }
    }
}

void LabController::Render()
{
    assert(!!mRenderContext);

    mRenderContext->RenderStart();

    if (mWorld)
    {
        mRenderContext->UploadSeaLevel(mWorld->GetOceanSurface().GetDepth());

        auto const & points = mWorld->GetShip().GetPoints();
        auto const & springs = mWorld->GetShip().GetSprings();
        auto const & triangles = mWorld->GetShip().GetTriangles();

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

        auto const colorChooser = [&](ElementIndex t, int e) -> std::tuple<rgbaColor, float>
            {
                float constexpr FloorThicknessAdjustment = 3.0f;

                NpcFloorGeometryType floorGeometry = triangles.GetSubSpringNpcFloorGeometry(t, e);
                switch (floorGeometry)
                {
                    case NpcFloorGeometryType::NotAFloor:
                    {
                        if (points.GetMaterial(springs.GetEndpointAIndex(triangles.GetSubSprings(t).SpringIndices[e]))->IsHull
                            && points.GetMaterial(springs.GetEndpointBIndex(triangles.GetSubSprings(t).SpringIndices[e]))->IsHull)
                        {
                            return { { 0x6a, 0x6a, 0x6a, 0x80 }, 1.5f };
                        }
                        else
                        {
                            return { { 0xca, 0xca, 0xca, 0xc0 }, 1.0f };
                        }
                    }

                    case NpcFloorGeometryType::Depth1H:
                    case NpcFloorGeometryType::Depth1V:
                    {
                        return { { 0x12, 0x1b, 0x54, 0xff }, FloorThicknessAdjustment };
                    }

                    case NpcFloorGeometryType::Depth2S1:
                    case NpcFloorGeometryType::Depth2S2:
                    {
                        return { { 0x12, 0x54, 0x1b, 0xff }, FloorThicknessAdjustment };
                    }
                }

                assert(false);
                return {{ 0x00, 0x00, 0x00, 0xff }, 0.0f};
            };

        for (auto t : triangles)
        {
            {
                auto [color, thicknessAdjustment] = colorChooser(t, 0);
                if (mWorld->GetNpcs().IsSpringHostingCurrentlySelectedParticle(triangles.GetSubSpringAIndex(t)))
                {
                    color.r = static_cast<rgbaColor::data_type>(std::min(color.r + 0x90, static_cast<int>(rgbaColor::data_type_max)));
                }
                mRenderContext->UploadEdge(
                    points.GetPosition(triangles.GetPointAIndex(t)),
                    points.GetPosition(triangles.GetPointBIndex(t)),
                    color,
                    thicknessAdjustment);
            }

            {
                auto [color, thicknessAdjustment] = colorChooser(t, 1);
                if (mWorld->GetNpcs().IsSpringHostingCurrentlySelectedParticle(triangles.GetSubSpringBIndex(t)))
                {
                    color.r = static_cast<rgbaColor::data_type>(std::min(color.r + 0x90, static_cast<int>(rgbaColor::data_type_max)));
                }
                mRenderContext->UploadEdge(
                    points.GetPosition(triangles.GetPointBIndex(t)),
                    points.GetPosition(triangles.GetPointCIndex(t)),
                    color,
                    thicknessAdjustment);
            }

            {
                auto [color, thicknessAdjustment] = colorChooser(t, 2);
                if (mWorld->GetNpcs().IsSpringHostingCurrentlySelectedParticle(triangles.GetSubSpringCIndex(t)))
                {
                    color.r = static_cast<rgbaColor::data_type>(std::min(color.r + 0x90, static_cast<int>(rgbaColor::data_type_max)));
                }
                mRenderContext->UploadEdge(
                    points.GetPosition(triangles.GetPointCIndex(t)),
                    points.GetPosition(triangles.GetPointAIndex(t)),
                    color,
                    thicknessAdjustment);
            }
        }

        mRenderContext->UploadEdgesEnd();

        //
        // Triangles
        //

        mRenderContext->UploadTrianglesStart();

        for (auto t : triangles)
        {
            std::optional<rgbaColor> color;

            if (mWorld->GetNpcs().IsTriangleConstrainingCurrentlySelectedParticle(t))
            {
                color = rgbaColor(107, 227, 107, 77);
            }
            else if (t == mWorld->GetNpcs().GetCurrentOriginTriangle())
            {
                color = rgbaColor(227, 107, 107, 77);
            }
            else if (
                triangles.GetSubSpringNpcFloorKind(t, 0) != NpcFloorKindType::NotAFloor
                && triangles.GetSubSpringNpcFloorKind(t, 1) != NpcFloorKindType::NotAFloor
                && triangles.GetSubSpringNpcFloorKind(t, 2) != NpcFloorKindType::NotAFloor)
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

        mWorld->GetNpcs().Upload(*mRenderContext);

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
    auto & points = mWorld->GetShip().GetPoints();

    //
    // Translation
    //

    vec2f const translation = mCurrentShipTranslationVelocity * GameParameters::SimulationStepTimeDuration<float>;

    // Update ship
    for (auto p : points)
    {
        points.SetPosition(p, points.GetPosition(p) + translation);
    }

    // Update camera pan
    mRenderContext->SetCameraWorldPosition(mRenderContext->GetCameraWorldPosition() + translation);

    //
    // Waves
    //

    mCurrentWavesAmplitude += (mTargetWavesAmplitude - mCurrentWavesAmplitude) * 0.05f;
    mCurrentWavesSpeed += (mTargetWavesSpeed - mCurrentWavesSpeed) * 0.05f;
    if (mCurrentWavesAmplitude != 0.0f && mCurrentWavesSpeed != 0.0f)
    {
        float const newArg = mLastWaveTimeArg + mCurrentWavesSpeed * GameParameters::SimulationStepTimeDuration<float>;
        float const rotationAngle = std::sinf(newArg) * mCurrentWavesAmplitude;
        float const deltaRotationAngle = rotationAngle - mLastWaveRotationAngle;
        for (auto p : points)
        {
            points.SetPosition(p, points.GetPosition(p).rotate(deltaRotationAngle));
        }

        mLastWaveRotationAngle = rotationAngle;
        mLastWaveTimeArg = newArg;
    }
}

std::optional<ElementIndex> LabController::TryPickVertex(vec2f const & screenCoordinates) const
{
    assert(!!mWorld);

    //
    // Find closest vertex within the radius
    //

    vec2f const worldCoordinates = ScreenToWorld(screenCoordinates);

    float constexpr SquareSearchRadius = GameParameters::VertexRadius * GameParameters::VertexRadius;

    float bestSquareDistance = std::numeric_limits<float>::max();
    ElementIndex bestPoint = NoneElementIndex;

    auto const & points = mWorld->GetShip().GetPoints();
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
    assert(!!mWorld);

    vec2f const worldOffset = ScreenOffsetToWorldOffset(screenOffset);

    mWorld->GetShip().GetPoints().SetPosition(
        pointIndex,
        mWorld->GetShip().GetPoints().GetPosition(pointIndex) + worldOffset);

    mWorld->GetNpcs().OnPointMoved(mCurrentSimulationTime);
}

void LabController::RotateShipBy(
    vec2f const & centerScreenCoordinates,
    float screenStride)
{
    assert(!!mWorld);

    vec2f const worldCenter = ScreenToWorld(centerScreenCoordinates);
    float const worldAngle = ScreenOffsetToWorldOffset(screenStride) * 0.05f;

    float const cosAngle = std::cos(worldAngle);
    float const sinAngle = std::sin(worldAngle);

    //
    // Rotate ship
    //

    auto & points = mWorld->GetShip().GetPoints();
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
    ////    sinAngle);
}

void LabController::RotateShipBy(
    ElementIndex particleIndex,
    float screenStride)
{
    assert(!!mWorld);

    //
    // Rotate ship
    //

    vec2f const worldCenter = mWorld->GetNpcs().GetParticles().GetPosition(particleIndex);
    float const worldAngle = ScreenOffsetToWorldOffset(screenStride) * 0.05f;

    float const cosAngle = std::cos(worldAngle);
    float const sinAngle = std::sin(worldAngle);

    auto & points = mWorld->GetShip().GetPoints();
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
    ////    sinAngle);
}

bool LabController::TrySelectOriginTriangle(vec2f const & screenCoordinates)
{
    assert(!!mWorld);

    vec2f const worldCoordinates = ScreenToWorld(screenCoordinates);

    ElementIndex const t = mWorld->GetShip().GetTriangles().FindContaining(worldCoordinates, mWorld->GetShip().GetPoints());
    if (t != NoneElementIndex)
    {
        mWorld->GetNpcs().SelectOriginTriangle(t);
        return true;
    }
    else
    {
        mWorld->GetNpcs().ResetOriginTriangle();
        return false;
    }
}

bool LabController::TrySelectParticle(vec2f const & screenCoordinates)
{
    assert(!!mWorld);

    //
    // Find closest particle within the radius
    //

    vec2f const worldCoordinates = ScreenToWorld(screenCoordinates);

    float constexpr SquareSearchRadius = GameParameters::ParticleRadius * GameParameters::ParticleRadius;

    float bestSquareDistance = std::numeric_limits<float>::max();
    ElementIndex bestParticle = NoneElementIndex;

    auto const & particles = mWorld->GetNpcs().GetParticles();
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
        mWorld->GetNpcs().SelectParticle(bestParticle);
        return true;
    }
    else
    {
        return false;
    }
}

std::optional<ElementIndex> LabController::TryPickNpcParticle(vec2f const & screenCoordinates) const
{
    assert(!!mWorld);

    //
    // Find closest particle within the radius
    //

    vec2f const worldCoordinates = ScreenToWorld(screenCoordinates);

    float constexpr SquareSearchRadius = GameParameters::ParticleRadius * GameParameters::ParticleRadius;

    float bestSquareDistance = std::numeric_limits<float>::max();
    ElementIndex bestParticle = NoneElementIndex;

    auto const & particles = mWorld->GetNpcs().GetParticles();
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

    assert(!!mWorld);

    vec2f const worldOffset = ScreenOffsetToWorldOffset(screenOffset);

    mWorld->GetNpcs().MoveParticleBy(
        npcParticleIndex,
        worldOffset,
        mCurrentSimulationTime);
}

void LabController::NotifyNpcParticleTrajectory(
    ElementIndex npcParticleIndex,
    vec2f const & targetScreenCoordinates)
{
    assert(!!mWorld);

    mWorld->GetNpcs().NotifyParticleTrajectory(
        npcParticleIndex,
        ScreenToWorld(targetScreenCoordinates));
}

void LabController::SetNpcParticleTrajectory(
    ElementIndex npcParticleIndex,
    vec2f const & targetScreenCoordinates)
{
    mWorld->GetNpcs().SetParticleTrajectory(
        npcParticleIndex,
        ScreenToWorld(targetScreenCoordinates));
}

void LabController::QueryNearestNpcParticleAt(vec2f const & screenCoordinates) const
{
    assert(!!mWorld);

    auto const nearestNpcParticle = TryPickNpcParticle(screenCoordinates);
    if (nearestNpcParticle.has_value())
    {
        mWorld->GetNpcs().GetParticles().Query(*nearestNpcParticle);
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

void LabController::SetWavesAmplitude(float wavesAmplitude)
{
    mTargetWavesAmplitude = wavesAmplitude;
}

void LabController::SetWavesSpeed(float wavesSpeed)
{
    mTargetWavesSpeed = wavesSpeed;
}

void LabController::QueryPointAt(vec2f const & screenCoordinates) const
{
    assert(!!mWorld);

    vec2f const worldCoordinates = ScreenToWorld(screenCoordinates);
    mWorld->QueryPointAt(worldCoordinates);
}

void LabController::FlipCurrentlySelectedHuman()
{
    assert(mWorld);

    if (mWorld->GetNpcs().GetCurrentlySelectedNpc().has_value())
    {
        mWorld->GetNpcs().FlipHumanWalk(*mWorld->GetNpcs().GetCurrentlySelectedNpc());
        //mModel->GetNpcs().FlipHumanFrontBack(0);
    }
}

void LabController::DoStepForVideo()
{
    assert(mWorld);

    ++mCurrentVideoStep;

    ////// Segment 1

    ////int matchIndex = 1;
    ////if (mCurrentVideoStep == matchIndex)
    ////{
    ////    // Setup camera
    ////    mRenderContext->SetCameraWorldPosition(vec2f(-4.42721558f + 8.0f, -5.58747005f));
    ////}

    ////++matchIndex;
    ////if (mCurrentVideoStep == matchIndex)
    ////{
    ////    // Place NPC
    ////    mWorld->GetNpcs().AddHumanNpc(
    ////        HumanNpcKindType::Passenger,
    ////        vec2f(-1.0f, -6.0f + 1.69f),
    ////        mCurrentSimulationTime);

    ////    // Select NPC
    ////    mWorld->GetNpcs().SelectNpc(0);
    ////}

    ////++matchIndex;
    ////if (mCurrentVideoStep == matchIndex)
    ////{
    ////    // Pause
    ////    SetSimulationControlState(SimulationControlStateType::Paused);
    ////}

    ////++matchIndex;
    ////if (mCurrentVideoStep == matchIndex)
    ////{
    ////    // Play
    ////    SetSimulationControlState(SimulationControlStateType::Play);
    ////}

    ////++matchIndex;
    ////if (mCurrentVideoStep == matchIndex)
    ////{
    ////    // Zoom
    ////    mRenderContext->SetZoom(1.5f);
    ////}

    ////// Segment 2

    ////int matchIndex = 1;
    ////if (mCurrentVideoStep == matchIndex)
    ////{
    ////    // Place NPC1
    ////    mWorld->GetNpcs().AddHumanNpc(
    ////        HumanNpcKindType::Passenger,
    ////        vec2f(-3.0f, -6.0f + 1.69f),
    ////        mCurrentSimulationTime);

    ////    // Place NPC2
    ////    mWorld->GetNpcs().AddHumanNpc(
    ////        HumanNpcKindType::Passenger,
    ////        vec2f(1.0f, -6.0f + 1.69f),
    ////        mCurrentSimulationTime);
    ////}

    ////++matchIndex;
    ////if (mCurrentVideoStep == matchIndex)
    ////{
    ////    // Pause
    ////    SetSimulationControlState(SimulationControlStateType::Paused);
    ////}

    ////++matchIndex;
    ////if (mCurrentVideoStep == matchIndex)
    ////{
    ////    // Play
    ////    SetSimulationControlState(SimulationControlStateType::Play);
    ////}

    ////++matchIndex;
    ////if (mCurrentVideoStep == matchIndex)
    ////{
    ////    AABB const shipAABB = mWorld->GetShip().GetPoints().GetAABB();

    ////    // Many more NPCs
    ////    for (int i = 0; i < 10; ++i)
    ////    {
    ////        float const posX = GameRandomEngine::GetInstance().GenerateUniformReal(shipAABB.BottomLeft.x, shipAABB.TopRight.x);
    ////        float const posY = GameRandomEngine::GetInstance().GenerateUniformReal(shipAABB.BottomLeft.y + GameParameters::HumanNpcGeometry::BodyLengthMean * 1.5f * mGameParameters.HumanNpcBodyLengthAdjustment, shipAABB.TopRight.y);

    ////        mWorld->GetNpcs().AddHumanNpc(
    ////            HumanNpcKindType::Passenger,
    ////            vec2f(posX, posY),
    ////            mCurrentSimulationTime);
    ////    }
    ////}

    ////// Segment 3

    ////int matchIndex = 1;
    ////if (mCurrentVideoStep == matchIndex)
    ////{
    ////    AABB const shipAABB = mWorld->GetShip().GetPoints().GetAABB();

    ////    for (size_t i = 0; i < GameParameters::MaxNpcs; ++i)
    ////    {
    ////        float const posX = GameRandomEngine::GetInstance().GenerateUniformReal(shipAABB.BottomLeft.x, shipAABB.TopRight.x);
    ////        float const posY = GameRandomEngine::GetInstance().GenerateUniformReal(shipAABB.BottomLeft.y + GameParameters::HumanNpcGeometry::BodyLengthMean * 1.5f * mGameParameters.HumanNpcBodyLengthAdjustment, shipAABB.TopRight.y);

    ////        bool const hasBeenAdded = mWorld->GetNpcs().AddHumanNpc(
    ////            HumanNpcKindType::Passenger,
    ////            vec2f(posX, posY),
    ////            mCurrentSimulationTime);

    ////        if (!hasBeenAdded)
    ////        {
    ////            throw GameException("Cannot add NPC!");
    ////        }
    ////    }
    ////}
    ////else
    ////{
    ////    AdjustZoom(0.995f);
    ////}
}

////////////////////////////////////////////////

std::optional<PickedObjectId<NpcId>> LabController::BeginPlaceNewFurnitureNpc(
    FurnitureNpcKindType furnitureKind,
    vec2f const & screenCoordinates)
{
    assert(!!mWorld);

    vec2f const worldCoordinates = ScreenToWorld(screenCoordinates);

    auto const pickedObjectId = mWorld->GetNpcs().BeginPlaceNewFurnitureNpc(
        furnitureKind,
        worldCoordinates,
        mCurrentSimulationTime);

    if (pickedObjectId.has_value())
    {
        mWorld->GetNpcs().SelectNpc(pickedObjectId->ObjectId);
    }

    return pickedObjectId;
}

std::optional<PickedObjectId<NpcId>> LabController::BeginPlaceNewHumanNpc(
    HumanNpcKindType humanKind,
    vec2f const & screenCoordinates)
{
    assert(!!mWorld);

    vec2f const worldCoordinates = ScreenToWorld(screenCoordinates);

    auto const pickedObjectId = mWorld->GetNpcs().BeginPlaceNewHumanNpc(
        humanKind,
        worldCoordinates,
        mCurrentSimulationTime);

    if (pickedObjectId.has_value())
    {
        mWorld->GetNpcs().SelectNpc(pickedObjectId->ObjectId);
    }

    return pickedObjectId;
}

std::optional<PickedObjectId<NpcId>> LabController::ProbeNpcAt(vec2f const & screenCoordinates) const
{
    assert(!!mWorld);

    vec2f const worldCoordinates = ScreenToWorld(screenCoordinates);

    return mWorld->GetNpcs().ProbeNpcAt(
        worldCoordinates,
        0.3f);
}

void LabController::BeginMoveNpc(NpcId id)
{
    assert(!!mWorld);

    mWorld->GetNpcs().BeginMoveNpc(
        id,
        mCurrentSimulationTime);
}

void LabController::MoveNpcTo(
    NpcId id,
    vec2f const & screenCoordinates,
    vec2f const & worldOffset)
{
    assert(!!mWorld);

    vec2f const worldCoordinates = ScreenToWorld(screenCoordinates);

    mWorld->GetNpcs().MoveNpcTo(
        id,
        worldCoordinates,
        worldOffset);
}

void LabController::EndMoveNpc(NpcId id)
{
    assert(!!mWorld);

    mWorld->GetNpcs().EndMoveNpc(id, mCurrentSimulationTime);
}

void LabController::CompleteNewNpc(NpcId id)
{
    assert(!!mWorld);

    mWorld->GetNpcs().CompleteNewNpc(id, mCurrentSimulationTime);
}

void LabController::RemoveNpc(NpcId id)
{
    assert(!!mWorld);

    mWorld->GetNpcs().RemoveNpc(id);

    if (id == mWorld->GetNpcs().GetCurrentlySelectedNpc())
    {
        mWorld->GetNpcs().DeselectNpc();
    }
}

void LabController::AbortNewNpc(NpcId id)
{
    assert(!!mWorld);

    mWorld->GetNpcs().AbortNewNpc(id);

    if (id == mWorld->GetNpcs().GetCurrentlySelectedNpc())
    {
        mWorld->GetNpcs().DeselectNpc();
    }
}

void LabController::HighlightNpc(
    NpcId id,
    NpcHighlightType highlight)
{
    assert(!!mWorld);

    mWorld->GetNpcs().HighlightNpc(
        id,
        highlight);
}

void LabController::SetNpcPanicLevelForAllHumans(float panicLevel)
{
    assert(!!mWorld);
    mWorld->GetNpcs().SetGeneralizedPanicLevelForAllHumans(panicLevel);
}

////////////////////////////////////////////////

void LabController::Reset(
    std::unique_ptr<Physics::Ship> ship,
    std::filesystem::path const & shipDefinitionFilepath)
{
    AABB const shipAABB = ship->GetPoints().GetAABB();

    //
    // Make new world
    //

    mWorld.reset();

    mWorld = std::make_unique<Physics::World>(
        mMaterialDatabase,
        mGameEventHandler,
        mOceanDepth);

    // Because we need to use a separate interface
    mWorld->GetNpcs().OnMassAdjustmentChanged(mMassAdjustment);
    mWorld->GetNpcs().OnGravityAdjustmentChanged(mGravityAdjustment);

    mWorld->AddShip(std::move(ship));

    mCurrentShipFilePath = shipDefinitionFilepath;

    //
    // Reset our state
    //

    mCurrentSimulationTime = 0.0f;

    mPerfStats.Reset();

    //
    // Auto-zoom & center
    //

    {
        vec2f const shipSize = shipAABB.GetSize();

        // Zoom to fit width and height (plus a nicely-looking margin)
        float const newZoom = std::min(
            mRenderContext->CalculateZoomForWorldWidth(shipSize.x + 5.0f),
            mRenderContext->CalculateZoomForWorldHeight(shipSize.y + 3.0f));
        mRenderContext->SetZoom(newZoom);

        // Center
        vec2f const objectCenter(
            (shipAABB.BottomLeft.x + shipAABB.TopRight.x) / 2.0f,
            (shipAABB.BottomLeft.y + shipAABB.TopRight.y) / 2.0f);
        mRenderContext->SetCameraWorldPosition(objectCenter);
    }

#ifndef _DEBUG

    //
    // Add initial NPCs
    //

    for (size_t i = 0; i < GameParameters::MaxNpcs; ++i)
    {
        float const posX = GameRandomEngine::GetInstance().GenerateUniformReal(shipAABB.BottomLeft.x, shipAABB.TopRight.x);
        float const posY = GameRandomEngine::GetInstance().GenerateUniformReal(shipAABB.BottomLeft.y + GameParameters::HumanNpcGeometry::BodyLengthMean * 1.5f * mGameParameters.HumanNpcBodyLengthAdjustment, shipAABB.TopRight.y);

        bool const hasBeenAdded = mWorld->GetNpcs().AddHumanNpc(
            HumanNpcKindType::Passenger,
            vec2f(posX, posY),
            mCurrentSimulationTime);

        if (!hasBeenAdded)
        {
            throw GameException("Cannot add NPC!");
        }
    }

#endif

    //
    // Publish reset
    //

    mGameEventHandler->OnBLabReset();
}

