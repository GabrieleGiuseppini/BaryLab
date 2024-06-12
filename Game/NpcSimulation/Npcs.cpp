/***************************************************************************************
* Original Author:      Gabriele Giuseppini
* Created:              2023-10-06
* Copyright:            Gabriele Giuseppini  (https://github.com/GabrieleGiuseppini)
***************************************************************************************/
#include "Physics.h"

#include <GameCore/Colors.h>
#include <GameCore/GameGeometry.h>

#include <algorithm>
#include <cassert>

namespace Physics {

/*
Main principles:
	- Global damping: when constrained, we only apply it to velocity *relative* to the mesh ("air moves with the ship")
*/

void Npcs::Update(
	float currentSimulationTime,
	GameParameters const & gameParameters)
{
	//
	// Update parameters
	//

	if (gameParameters.HumanNpcBodyLengthAdjustment != mCurrentHumanNpcBodyLengthAdjustment)
	{
		mCurrentHumanNpcBodyLengthAdjustment = gameParameters.HumanNpcBodyLengthAdjustment;

		RecalculateHumanNpcDipoleLengths();
	}

	if (gameParameters.HumanNpcWalkingSpeedAdjustment != mCurrentHumanNpcWalkingSpeedAdjustment)
	{
		mCurrentHumanNpcWalkingSpeedAdjustment = gameParameters.HumanNpcWalkingSpeedAdjustment;
	}

	if (gameParameters.NpcSpringReductionFractionAdjustment != mCurrentSpringReductionFractionAdjustment
		|| gameParameters.NpcSpringDampingCoefficientAdjustment != mCurrentSpringDampingCoefficientAdjustment)
	{
		mCurrentSpringReductionFractionAdjustment = gameParameters.NpcSpringReductionFractionAdjustment;
		mCurrentSpringDampingCoefficientAdjustment = gameParameters.NpcSpringDampingCoefficientAdjustment;

		RecalculateSpringForceParameters();
	}

	//
	// Update NPCs' state
	//

	UpdateNpcs(currentSimulationTime, gameParameters);
}

void Npcs::RenderUpload(
	RenderContext & renderContext,
	PerfStats & perfStats)
{
#ifdef IN_BARYLAB
	if (renderContext.GetNpcRenderMode() == NpcRenderModeType::Physical)
	{
		renderContext.UploadNpcParticlesStart();
		renderContext.UploadNpcSpringsStart();

		for (ShipId shipId = 0; shipId < mShips.size(); ++shipId)
		{
			if (mShips[shipId].has_value())
			{
				ShipRenderContext & shipRenderContext = renderContext.GetShipRenderContext(shipId);

				for (NpcId const npcId : mShips[shipId]->Npcs)
				{
					assert(mStateBuffer[npcId].has_value());
					auto const & state = *mStateBuffer[npcId];

					auto const planeId = state.CurrentPlaneId.has_value()
						? *(state.CurrentPlaneId)
						: mShips[shipId]->ShipMesh.GetMaxPlaneId();

					// Particles
					for (auto const & particle: state.ParticleMesh.Particles)
					{
						shipRenderContext.UploadNpcParticle(
							planeId,
							mParticles.GetPosition(particle.ParticleIndex),
							mParticles.GetRenderColor(particle.ParticleIndex),
							1.0f,
							state.Highlight);
					}

					// Springs
					for (auto const & spring : state.ParticleMesh.Springs)
					{
						renderContext.UploadNpcSpring(
							planeId,
							mParticles.GetPosition(spring.EndpointAIndex),
							mParticles.GetPosition(spring.EndpointBIndex),
							rgbaColor(0x4a, 0x4a, 0x4a, 0xff));
					}
				}
			}
		}

		renderContext.UploadNpcSpringsEnd();
		renderContext.UploadNpcParticlesEnd();

		return;
	}
#endif

	assert(renderContext.GetNpcRenderMode() == NpcRenderModeType::Limbs);

	// Futurework
	renderContext.UploadNpcQuadsStart(mFurnitureNpcCount * 4 + mHumanNpcCount * (6 + 2));

#ifdef IN_BARYLAB
	// For furniture
	renderContext.UploadNpcParticlesStart();
#endif

	auto const startTime = GameChronometer::now();

	for (ShipId shipId = 0; shipId < mShips.size(); ++shipId)
	{
		if (mShips[shipId].has_value())
		{
			ShipRenderContext & shipRenderContext = renderContext.GetShipRenderContext(shipId);

			for (NpcId const npcId : mShips[shipId]->Npcs)
			{
				assert(mStateBuffer[npcId].has_value());
				auto const & state = *mStateBuffer[npcId];

				RenderNpc(state, shipRenderContext);
			}
		}
	}

	perfStats.TotalNpcRenderUploadDuration.Update(GameChronometer::now() - startTime);

#ifdef IN_BARYLAB
	// For furniture
	renderContext.UploadNpcParticlesEnd();
#endif

	renderContext.UploadNpcQuadsEnd();

#ifdef IN_BARYLAB

	//
	// Particle trajectories
	//

	renderContext.UploadParticleTrajectoriesStart();

	if (mCurrentParticleTrajectoryNotification)
	{
		renderContext.UploadParticleTrajectory(
			mParticles.GetPosition(mCurrentParticleTrajectoryNotification->ParticleIndex),
			mCurrentParticleTrajectoryNotification->TargetPosition,
			rgbaColor(0xc0, 0xc0, 0xc0, 0xff));
	}

	if (mCurrentParticleTrajectory)
	{
		renderContext.UploadParticleTrajectory(
			mParticles.GetPosition(mCurrentParticleTrajectory->ParticleIndex),
			mCurrentParticleTrajectory->TargetPosition,
			rgbaColor(0x99, 0x99, 0x99, 0xff));
	}

	renderContext.UploadParticleTrajectoriesEnd();

#endif
}

///////////////////////////////

void Npcs::OnShipAdded(Ship const & ship)
{
	size_t const s = static_cast<size_t>(ship.GetId());

	// Make room for ship
	if (s >= mShips.size())
	{
		mShips.resize(s + 1);
	}

	// We do not know about this ship yet
	assert(!mShips[s].has_value());

	// Initialize NPC Ship
	mShips[s].emplace(ship);
}

void Npcs::OnShipRemoved(ShipId shipId)
{
	size_t const s = static_cast<size_t>(shipId);

	// We know about this ship
	assert(s < mShips.size());
	assert(mShips[s].has_value());

	//
	// Handle destruction of all NPCs of this NPC ship
	//

	bool doPublishHumanNpcStats = false;

	for (auto const npcId : mShips[s]->Npcs)
	{
		assert(mStateBuffer[npcId].has_value());

		if (mStateBuffer[npcId]->Kind == NpcKindType::Human)
		{
			if (mStateBuffer[npcId]->CurrentRegime == StateType::RegimeType::Constrained)
			{
				assert(mConstrainedRegimeHumanNpcCount > 0);
				--mConstrainedRegimeHumanNpcCount;
				doPublishHumanNpcStats = true;
			}
			else if (mStateBuffer[npcId]->CurrentRegime == StateType::RegimeType::Free)
			{
				assert(mFreeRegimeHumanNpcCount > 0);
				--mFreeRegimeHumanNpcCount;
				doPublishHumanNpcStats = true;
			}

			--mHumanNpcCount;
		}
		else
		{
			assert(mStateBuffer[npcId]->Kind == NpcKindType::Furniture);

			--mFurnitureNpcCount;
		}

		mStateBuffer[npcId].reset();
	}

	if (doPublishHumanNpcStats)
	{
		PublishHumanNpcStats();
	}

	//
	// Destroy NPC ship
	//

	mShips[s].reset();
}

std::optional<PickedObjectId<NpcId>> Npcs::BeginPlaceNewFurnitureNpc(
	FurnitureNpcKindType furnitureKind,
	vec2f const & worldCoordinates,
	float /*currentSimulationTime*/)
{
	//
	// Check if there are enough NPCs
	//

	if ((mHumanNpcCount + mFurnitureNpcCount) >= GameParameters::MaxNpcs)
	{
		return std::nullopt;
	}

	//
	// Create NPC
	//

	StateType::ParticleMeshType particleMesh;

	switch (furnitureKind)
	{
		case FurnitureNpcKindType::Quad:
		{
			// Check if there are enough particles

			if (mParticles.GetRemainingParticlesCount() < 4)
			{
				return std::nullopt;
			}

			// Particles

			float constexpr SideLength = 1.5f;
			float constexpr DiagonalSideLength = SideLength * 1.41421356f; // sqrt(2)

			auto const & furnitureMaterial = mMaterialDatabase.GetNpcMaterial(NpcMaterial::KindType::Furniture);

			for (int p = 0; p < 4; ++p)
			{
				// CW order
				vec2f particlePosition = worldCoordinates;
				if (p == 1 || p == 2)
				{
					particlePosition.x += SideLength;
				}
				if (p == 2 || p == 3)
				{
					particlePosition.y -= SideLength;
				}

				auto const particleIndex = mParticles.Add(
					furnitureMaterial.Mass,
					furnitureMaterial.StaticFriction,
					furnitureMaterial.KineticFriction,
					furnitureMaterial.Elasticity,
					furnitureMaterial.BuoyancyVolumeFill,
					particlePosition,
					furnitureMaterial.RenderColor);

				particleMesh.Particles.emplace_back(particleIndex, std::nullopt);
			}

			// Springs

			float const massFactor =
				(furnitureMaterial.Mass * furnitureMaterial.Mass)
				/ (furnitureMaterial.Mass + furnitureMaterial.Mass);

			StateType::NpcSpringStateType baseSpring(
				NoneElementIndex,
				NoneElementIndex,
				0,
				furnitureMaterial.SpringReductionFraction,
				furnitureMaterial.SpringDampingCoefficient,
				massFactor);

			// 0 - 1
			{
				baseSpring.EndpointAIndex = particleMesh.Particles[0].ParticleIndex;
				baseSpring.EndpointBIndex = particleMesh.Particles[1].ParticleIndex;
				baseSpring.DipoleLength = SideLength;
				auto & spring = particleMesh.Springs.emplace_back(baseSpring);
				RecalculateSpringForceParameters(spring);
			}

			// 0 | 3
			{
				baseSpring.EndpointAIndex = particleMesh.Particles[0].ParticleIndex;
				baseSpring.EndpointBIndex = particleMesh.Particles[3].ParticleIndex;
				baseSpring.DipoleLength = SideLength;
				auto & spring = particleMesh.Springs.emplace_back(baseSpring);
				RecalculateSpringForceParameters(spring);
			}

			// 0 \ 2
			{
				baseSpring.EndpointAIndex = particleMesh.Particles[0].ParticleIndex;
				baseSpring.EndpointBIndex = particleMesh.Particles[2].ParticleIndex;
				baseSpring.DipoleLength = DiagonalSideLength;
				auto & spring = particleMesh.Springs.emplace_back(baseSpring);
				RecalculateSpringForceParameters(spring);
			}

			// 1 | 2
			{
				baseSpring.EndpointAIndex = particleMesh.Particles[1].ParticleIndex;
				baseSpring.EndpointBIndex = particleMesh.Particles[2].ParticleIndex;
				baseSpring.DipoleLength = SideLength;
				auto & spring = particleMesh.Springs.emplace_back(baseSpring);
				RecalculateSpringForceParameters(spring);
			}

			// 2 - 3
			{
				baseSpring.EndpointAIndex = particleMesh.Particles[2].ParticleIndex;
				baseSpring.EndpointBIndex = particleMesh.Particles[3].ParticleIndex;
				baseSpring.DipoleLength = SideLength;
				auto & spring = particleMesh.Springs.emplace_back(baseSpring);
				RecalculateSpringForceParameters(spring);
			}

			// 1 / 3
			{
				baseSpring.EndpointAIndex = particleMesh.Particles[1].ParticleIndex;
				baseSpring.EndpointBIndex = particleMesh.Particles[3].ParticleIndex;
				baseSpring.DipoleLength = DiagonalSideLength;
				auto & spring = particleMesh.Springs.emplace_back(baseSpring);
				RecalculateSpringForceParameters(spring);
			}

			break;
		}

		case FurnitureNpcKindType::Particle:
		{
			// Check if there are enough particles

			if (mParticles.GetRemainingParticlesCount() < 1)
			{
				return std::nullopt;
			}

			// Primary

			auto const & furnitureMaterial = mMaterialDatabase.GetNpcMaterial(NpcMaterial::KindType::Furniture);

			auto const primaryParticleIndex = mParticles.Add(
				furnitureMaterial.Mass,
				furnitureMaterial.StaticFriction,
				furnitureMaterial.KineticFriction,
				furnitureMaterial.Elasticity,
				furnitureMaterial.BuoyancyVolumeFill,
				worldCoordinates,
				furnitureMaterial.RenderColor);

			particleMesh.Particles.emplace_back(primaryParticleIndex, std::nullopt);

			break;
		}
	}

	// Furniture

	StateType::KindSpecificStateType::FurnitureNpcStateType furnitureState = StateType::KindSpecificStateType::FurnitureNpcStateType(
		furnitureKind);

	//
	// Store NPC
	//

	NpcId const npcId = GetNewNpcId();

	// This NPC begins its journey on the topmost ship, just
	// to make sure it's at the nearest Z
	ShipId const shipId = GetTopmostShipId();

	mStateBuffer[npcId].emplace(
		npcId,
		NpcKindType::Furniture,
		shipId, // Topmost ship ID
		std::nullopt, // Topmost plane ID
		StateType::RegimeType::BeingPlaced,
		std::move(particleMesh),
		StateType::KindSpecificStateType(std::move(furnitureState)));

	assert(mShips[shipId].has_value());
	mShips[shipId]->Npcs.push_back(npcId);

	//
	// Update stats
	//

	++mFurnitureNpcCount;

	return PickedObjectId<NpcId>(npcId, vec2f::zero());
}

std::optional<PickedObjectId<NpcId>> Npcs::BeginPlaceNewHumanNpc(
	HumanNpcKindType humanKind,
	vec2f const & worldCoordinates,
	float currentSimulationTime)
{
	//
	// Check if there are enough NPCs and particles
	//

	if ((mHumanNpcCount + mFurnitureNpcCount) >= GameParameters::MaxNpcs || mParticles.GetRemainingParticlesCount() < 2)
	{
		return std::nullopt;
	}

	//
	// Create NPC
	//

	StateType::ParticleMeshType particleMesh;

	// Decide height

	float const height = GameRandomEngine::GetInstance().GenerateNormalReal(
		GameParameters::HumanNpcGeometry::BodyLengthMean,
		GameParameters::HumanNpcGeometry::BodyLengthStdDev);

	// Feet (primary)

	auto const & feetMaterial = mMaterialDatabase.GetNpcMaterial(NpcMaterial::KindType::HumanFeet);

	auto const primaryParticleIndex = mParticles.Add(
		feetMaterial.Mass,
		feetMaterial.StaticFriction,
		feetMaterial.KineticFriction,
		feetMaterial.Elasticity,
		feetMaterial.BuoyancyVolumeFill,
		worldCoordinates - vec2f(0.0f, 1.0f) * height * mCurrentHumanNpcBodyLengthAdjustment,
		feetMaterial.RenderColor);

	auto & primaryParticle = particleMesh.Particles.emplace_back(primaryParticleIndex, std::nullopt);

	// Head (secondary)

	auto const & headMaterial = mMaterialDatabase.GetNpcMaterial(NpcMaterial::KindType::HumanHead);

	auto const secondaryParticleIndex = mParticles.Add(
		headMaterial.Mass,
		headMaterial.StaticFriction,
		headMaterial.KineticFriction,
		headMaterial.Elasticity,
		headMaterial.BuoyancyVolumeFill,
		worldCoordinates,
		headMaterial.RenderColor);

	auto & secondaryParticle = particleMesh.Particles.emplace_back(secondaryParticleIndex, std::nullopt);

	// Dipole spring

	float const massFactor =
		(feetMaterial.Mass * headMaterial.Mass)
		/ (feetMaterial.Mass + headMaterial.Mass);

	auto & dipoleSpring = particleMesh.Springs.emplace_back(
		primaryParticleIndex,
		secondaryParticleIndex,
		CalculateHumanNpcDipoleLength(height),
		headMaterial.SpringReductionFraction,
		headMaterial.SpringDampingCoefficient,
		massFactor);

	RecalculateSpringForceParameters(dipoleSpring);

	// Human

	float widthMultiplier;
	if (GameRandomEngine::GetInstance().Choose(2) == 0)
	{
		// Narrow
		widthMultiplier = 1.0f - std::min(
			std::abs(GameRandomEngine::GetInstance().GenerateNormalReal(0.0f, GameParameters::HumanNpcGeometry::BodyWidthNarrowMultiplierStdDev)),
			2.0f * GameParameters::HumanNpcGeometry::BodyWidthNarrowMultiplierStdDev);
	}
	else
	{
		// Wide
		widthMultiplier = 1.0f + std::min(
			std::abs(GameRandomEngine::GetInstance().GenerateNormalReal(0.0f, GameParameters::HumanNpcGeometry::BodyWidthWideMultiplierStdDev)),
			3.0f * GameParameters::HumanNpcGeometry::BodyWidthWideMultiplierStdDev);
	}

	float const walkingSpeedBase =
		1.0f
		* height / 1.65f; // Just comes from 1m/s looking good when human is 1.65

	StateType::KindSpecificStateType::HumanNpcStateType humanState = StateType::KindSpecificStateType::HumanNpcStateType(
		humanKind,
		height,
		widthMultiplier,
		walkingSpeedBase,
		StateType::KindSpecificStateType::HumanNpcStateType::BehaviorType::BeingPlaced,
		currentSimulationTime);

	// Frontal
	humanState.CurrentFaceOrientation = 1.0f;
	humanState.CurrentFaceDirectionX = 0.0f;

	//
	// Store NPC
	//

	NpcId const npcId = GetNewNpcId();

	// This NPC begins its journey on the topmost ship, just
	// to make sure it's at the nearest Z
	ShipId const shipId = GetTopmostShipId();

	mStateBuffer[npcId].emplace(
		npcId,
		NpcKindType::Human,
		shipId, // Topmost ship ID
		std::nullopt, // Topmost plane ID
		StateType::RegimeType::BeingPlaced,
		std::move(particleMesh),
		StateType::KindSpecificStateType(std::move(humanState)));

	assert(mShips[shipId].has_value());
	mShips[shipId]->Npcs.push_back(npcId);

	//
	// Update stats
	//

	++mHumanNpcCount;

	return PickedObjectId<NpcId>(npcId, vec2f::zero());
}

std::optional<PickedObjectId<NpcId>> Npcs::ProbeNpcAt(
	vec2f const & position,
	GameParameters const & gameParameters) const
{
	struct NearestNpcType
	{
		NpcId Id{ NoneNpcId };
		float SquareDistance{ std::numeric_limits<float>::max() };
	};

	NearestNpcType nearestOnPlaneNpc;
	NearestNpcType nearestOffPlaneNpc;

	//
	// Determine ship and plane of this position - if any
	//

	std::pair<ShipId, PlaneId> probeDepth;

	// Find topmost triangle containing this position
	auto const topmostTriangle = FindTopmostTriangleContaining(position);
	if (topmostTriangle)
	{
		assert(topmostTriangle->GetShipId() < mShips.size());
		assert(mShips[topmostTriangle->GetShipId()].has_value());
		auto const & ship = *mShips[topmostTriangle->GetShipId()];

		ElementIndex const trianglePointIndex = ship.ShipMesh.GetTriangles().GetPointAIndex(topmostTriangle->GetLocalObjectId());
		PlaneId const planeId = ship.ShipMesh.GetPoints().GetPlaneId(trianglePointIndex);

		probeDepth = { ship.ShipMesh.GetId(), planeId };
	}
	else
	{
		probeDepth = { 0, 0 }; // Bottommost
	}

	//
	// Visit all NPCs and find winner. if any
	//

	auto particleVisitor = [&](ElementIndex candidateParticleIndex, float squareSearchRadius, StateType const & npc, ShipId const & shipId)
		{
			vec2f const candidateNpcPosition = mParticles.GetPosition(candidateParticleIndex);
			float const squareDistance = (candidateNpcPosition - position).squareLength();
			if (squareDistance < squareSearchRadius)
			{
				auto const candidatePlaneId = npc.CurrentPlaneId.value_or(std::numeric_limits<PlaneId>::max());
				if (std::make_pair(shipId, candidatePlaneId) >= probeDepth)
				{
					// It's on-plane
					if (squareDistance < nearestOnPlaneNpc.SquareDistance)
					{
						nearestOnPlaneNpc = { npc.Id, squareDistance };
					}
				}
				else
				{
					// It's off-plane
					if (squareDistance < nearestOffPlaneNpc.SquareDistance)
					{
						nearestOffPlaneNpc = { npc.Id, squareDistance };
					}
				}
			}
		};

	for (auto const & state : mStateBuffer)
	{
		if (state.has_value())
		{
			auto const & ship = mShips[state->CurrentShipId];

			switch (state->Kind)
			{
				case NpcKindType::Furniture:
				{
					// All particles

					float const squareSearchRadius = std::max(
						gameParameters.ToolSearchRadius * gameParameters.ToolSearchRadius,
						1.2f * 1.2f); // Allow large furniture to be found when picked inside

					for (auto const & particle : state->ParticleMesh.Particles)
					{
						particleVisitor(
							particle.ParticleIndex,
							squareSearchRadius,
							*state,
							state->CurrentShipId);
					}

					break;
				}

				case NpcKindType::Human:
				{
					// Head

					float const squareSearchRadius = gameParameters.ToolSearchRadius * gameParameters.ToolSearchRadius;

					assert(state->ParticleMesh.Particles.size() == 2);

					particleVisitor(
						state->ParticleMesh.Particles[1].ParticleIndex,
						squareSearchRadius,
						*state,
						state->CurrentShipId);

					break;
				}
			}
		}
	}

	//
	// Pick a winner - on-plane has higher prio than off-place
	//

	NpcId foundId = NoneNpcId;
	if (nearestOnPlaneNpc.Id != NoneNpcId)
	{
		foundId = nearestOnPlaneNpc.Id;
	}
	else if (nearestOffPlaneNpc.Id != NoneNpcId)
	{
		foundId = nearestOffPlaneNpc.Id;
	}

	if (foundId != NoneNpcId)
	{
		assert(mStateBuffer[foundId].has_value());

		int referenceParticleOrdinal = mStateBuffer[foundId]->Kind == NpcKindType::Furniture
			? 0
			: 1;

		ElementIndex referenceParticleIndex = mStateBuffer[foundId]->ParticleMesh.Particles[referenceParticleOrdinal].ParticleIndex;

		return PickedObjectId<NpcId>(
			foundId,
			position - mParticles.GetPosition(referenceParticleIndex));
	}
	else
	{
		return std::nullopt;
	}
}

void Npcs::BeginMoveNpc(
	NpcId id,
	float currentSimulationTime)
{
	assert(mStateBuffer[id].has_value());
	auto & npc = *mStateBuffer[id];

	//
	// Move NPC to topmost ship and its topmost plane
	//

	TransferNpcToShip(npc, GetTopmostShipId());
	npc.CurrentPlaneId = std::nullopt;

	//
	// Move NPC to BeingPlaced
	//

	// All particles become free
	for (auto & particle : npc.ParticleMesh.Particles)
	{
		particle.ConstrainedState.reset();
	}

	// Change regime
	auto const oldRegime = npc.CurrentRegime;
	npc.CurrentRegime = StateType::RegimeType::BeingPlaced;

	if (npc.Kind == NpcKindType::Human)
	{
		// Change behavior
		npc.KindSpecificState.HumanNpcState.TransitionToState(
			StateType::KindSpecificStateType::HumanNpcStateType::BehaviorType::BeingPlaced,
			currentSimulationTime);
	}

	//
	// Maintain stats
	//

	if (npc.Kind == NpcKindType::Human)
	{
		if (oldRegime == StateType::RegimeType::Constrained)
		{
			assert(mConstrainedRegimeHumanNpcCount > 0);
			--mConstrainedRegimeHumanNpcCount;
			PublishHumanNpcStats();
		}
		else if (oldRegime == StateType::RegimeType::Free)
		{
			assert(mFreeRegimeHumanNpcCount > 0);
			--mFreeRegimeHumanNpcCount;
			PublishHumanNpcStats();
		}
	}
}

void Npcs::MoveNpcTo(
	NpcId id,
	vec2f const & position,
	vec2f const & offset)
{
	assert(mStateBuffer[id].has_value());
	assert(mStateBuffer[id]->CurrentRegime == StateType::RegimeType::BeingPlaced);

	vec2f const newTargetPosition = position - offset;

	float constexpr InertialVelocityFactor = 0.5f; // Magic number for how much velocity we impart

	switch (mStateBuffer[id]->Kind)
	{
		case NpcKindType::Furniture:
		{
			// Move all particles

			// Get position of primary as reference
			assert(mStateBuffer[id]->ParticleMesh.Particles.size() > 0);
			vec2f const oldPrimaryPosition = mParticles.GetPosition(mStateBuffer[id]->ParticleMesh.Particles[0].ParticleIndex);

			// Calculate new absolute velocity of primary
			vec2f const newPrimaryAbsoluteVelocity = (newTargetPosition - oldPrimaryPosition) / GameParameters::SimulationTimeStepDuration * InertialVelocityFactor;

			// Move all particles
			for (auto & particle : mStateBuffer[id]->ParticleMesh.Particles)
			{
				auto const particleIndex = particle.ParticleIndex;
				vec2f const newPosition = newTargetPosition + (mParticles.GetPosition(particleIndex) - oldPrimaryPosition); // Set to new position + offset from reference
				mParticles.SetPosition(particleIndex, newPosition);
				mParticles.SetVelocity(particleIndex, newPrimaryAbsoluteVelocity);

				if (particle.ConstrainedState.has_value())
				{
					// We can only assume here, and we assume the ship is still and since the user doesn't move with the ship,
					// all this velocity is also relative to mesh
					particle.ConstrainedState->MeshRelativeVelocity = newPrimaryAbsoluteVelocity;
				}
			}

			break;
		}

		case NpcKindType::Human:
		{
			// Move secondary particle

			assert(mStateBuffer[id]->ParticleMesh.Particles.size() == 2);

			auto const particleIndex = mStateBuffer[id]->ParticleMesh.Particles[1].ParticleIndex;
			vec2f const oldPosition = mParticles.GetPosition(particleIndex);
			mParticles.SetPosition(particleIndex, newTargetPosition);
			vec2f const absoluteVelocity = (newTargetPosition - oldPosition) / GameParameters::SimulationTimeStepDuration * InertialVelocityFactor;
			mParticles.SetVelocity(particleIndex, absoluteVelocity);

			if (mStateBuffer[id]->ParticleMesh.Particles[1].ConstrainedState.has_value())
			{
				// We can only assume here, and we assume the ship is still and since the user doesn't move with the ship,
				// all this velocity is also relative to mesh
				mStateBuffer[id]->ParticleMesh.Particles[1].ConstrainedState->MeshRelativeVelocity = absoluteVelocity;
			}

			break;
		}
	}
}

void Npcs::EndMoveNpc(
	NpcId id,
	float currentSimulationTime)
{
	assert(mStateBuffer[id].has_value());

	auto & npc = *mStateBuffer[id];

	assert(npc.CurrentRegime == StateType::RegimeType::BeingPlaced);

	ResetNpcStateToWorld(
		npc,
		currentSimulationTime);

	OnMayBeNpcRegimeChanged(
		StateType::RegimeType::BeingPlaced,
		npc);

#ifdef IN_BARYLAB
	// Select NPC's primary particle
	SelectParticle(npc.ParticleMesh.Particles[0].ParticleIndex);
#endif
}

void Npcs::CompleteNewNpc(
	NpcId id,
	float currentSimulationTime)
{
	EndMoveNpc(id, currentSimulationTime);
}

void Npcs::RemoveNpc(NpcId id)
{
	assert(mStateBuffer[id].has_value());

	//
	// Maintain stats
	//

	if (mStateBuffer[id]->Kind == NpcKindType::Human)
	{
		if (mStateBuffer[id]->CurrentRegime == StateType::RegimeType::Constrained)
		{
			assert(mConstrainedRegimeHumanNpcCount > 0);
			--mConstrainedRegimeHumanNpcCount;
			PublishHumanNpcStats();
		}
		else if (mStateBuffer[id]->CurrentRegime == StateType::RegimeType::Free)
		{
			assert(mFreeRegimeHumanNpcCount > 0);
			--mFreeRegimeHumanNpcCount;
			PublishHumanNpcStats();
		}

		--mHumanNpcCount;
	}
	else
	{
		assert(mStateBuffer[id]->Kind == NpcKindType::Furniture);

		--mFurnitureNpcCount;
	}

	//
	// Update ship indices
	//

	assert(mShips[mStateBuffer[id]->CurrentShipId].has_value());

	auto it = std::find(
		mShips[mStateBuffer[id]->CurrentShipId]->Npcs.begin(),
		mShips[mStateBuffer[id]->CurrentShipId]->Npcs.end(),
		id);

	assert(it != mShips[mStateBuffer[id]->CurrentShipId]->Npcs.end());

	mShips[mStateBuffer[id]->CurrentShipId]->Npcs.erase(it);

	//
	// Get rid of NPC
	//

	mStateBuffer[id].reset();
}

void Npcs::AbortNewNpc(NpcId id)
{
	RemoveNpc(id);
}

void Npcs::HighlightNpc(
	NpcId id,
	NpcHighlightType highlight)
{
	assert(mStateBuffer[id].has_value());
	mStateBuffer[id]->Highlight = highlight;
}

void Npcs::SetGeneralizedPanicLevelForAllHumans(float panicLevel)
{
	for (auto & npc : mStateBuffer)
	{
		if (npc.has_value())
		{
			if (npc->Kind == NpcKindType::Human)
			{
				npc->KindSpecificState.HumanNpcState.GeneralizedPanicLevel = panicLevel;
			}
		}
	}
}

/////////////////////////////// Barylab-specific

#ifdef IN_BARYLAB

bool Npcs::AddHumanNpc(
	HumanNpcKindType humanKind,
	vec2f const & worldCoordinates,
	float currentSimulationTime)
{
	auto const npcId = BeginPlaceNewHumanNpc(
		humanKind,
		worldCoordinates,
		currentSimulationTime);

	if (npcId.has_value())
	{
		CompleteNewNpc(
			npcId->ObjectId,
			currentSimulationTime);

		return true;
	}
	else
	{
		return false;
	}
}

void Npcs::FlipHumanWalk(int npcIndex)
{
	if (npcIndex < mStateBuffer.size()
		&& mStateBuffer[npcIndex].has_value()
		&& mStateBuffer[npcIndex]->Kind == NpcKindType::Human
		&& mStateBuffer[npcIndex]->KindSpecificState.HumanNpcState.CurrentBehavior == StateType::KindSpecificStateType::HumanNpcStateType::BehaviorType::Constrained_Walking)
	{
		FlipHumanWalk(mStateBuffer[npcIndex]->KindSpecificState.HumanNpcState, StrongTypedTrue<_DoImmediate>);
	}
}

void Npcs::FlipHumanFrontBack(int npcIndex)
{
	if (npcIndex < mStateBuffer.size()
		&& mStateBuffer[npcIndex].has_value()
		&& mStateBuffer[npcIndex]->Kind == NpcKindType::Human
		&& mStateBuffer[npcIndex]->KindSpecificState.HumanNpcState.CurrentBehavior == StateType::KindSpecificStateType::HumanNpcStateType::BehaviorType::Constrained_Walking)
	{
		auto & humanState = mStateBuffer[npcIndex]->KindSpecificState.HumanNpcState;

		if (humanState.CurrentFaceOrientation != 0.0f)
		{
			humanState.CurrentFaceOrientation *= -1.0f;
		}
	}
}

void Npcs::MoveParticleBy(
	ElementIndex particleIndex,
	vec2f const & offset,
	float currentSimulationTime)
{
	//
	// Move particle
	//

	mParticles.SetPosition(
		particleIndex,
		mParticles.GetPosition(particleIndex) + offset);

	mParticles.SetVelocity(
		particleIndex,
		vec2f::zero()); // Zero-out velocity

	//
	// Re-initialize state of NPC that contains this particle
	//

	for (auto & state : mStateBuffer)
	{
		if (state.has_value())
		{
			if (std::any_of(
				state->ParticleMesh.Particles.cbegin(),
				state->ParticleMesh.Particles.cend(),
				[particleIndex](auto const & p) { return p.ParticleIndex == particleIndex;}))
			{
				auto const oldRegime = state->CurrentRegime;

				ResetNpcStateToWorld(*state, currentSimulationTime);

				OnMayBeNpcRegimeChanged(oldRegime, *state);

				//
				// Select particle
				//

				SelectParticle(particleIndex);

				//
				// Reset trajectories
				//

				mCurrentParticleTrajectory.reset();
				mCurrentParticleTrajectoryNotification.reset();

				break;
			}
		}
	}
}

void Npcs::RotateParticlesWithShip(
	vec2f const & centerPos,
	float cosAngle,
	float sinAngle)
{
	//
	// Rotate particles
	//

	for (auto const & ship : mShips)
	{
		if (ship.has_value())
		{
			for (auto const n : ship->Npcs)
			{
				auto & state = mStateBuffer[n];
				assert(state.has_value());

				for (auto & particle : state->ParticleMesh.Particles)
				{
					RotateParticleWithShip(
						particle,
						centerPos,
						cosAngle,
						sinAngle,
						ship->ShipMesh);
				}
			}
		}
	}
}

void Npcs::RotateParticleWithShip(
	StateType::NpcParticleStateType const & npcParticleState,
	vec2f const & centerPos,
	float cosAngle,
	float sinAngle,
	Ship const & shipMesh)
{
	vec2f newPosition;

	if (npcParticleState.ConstrainedState.has_value())
	{
		// Simply set position from current bary coords

		newPosition = shipMesh.GetTriangles().FromBarycentricCoordinates(
			npcParticleState.ConstrainedState->CurrentBCoords.BCoords,
			npcParticleState.ConstrainedState->CurrentBCoords.TriangleElementIndex,
			shipMesh.GetPoints());
	}
	else
	{
		// Rotate particle

		vec2f const centeredPos = mParticles.GetPosition(npcParticleState.ParticleIndex) - centerPos;
		vec2f const rotatedPos = vec2f(
			centeredPos.x * cosAngle - centeredPos.y * sinAngle,
			centeredPos.x * sinAngle + centeredPos.y * cosAngle);

		newPosition = rotatedPos + centerPos;
	}

	mParticles.SetPosition(npcParticleState.ParticleIndex, newPosition);
}

void Npcs::OnPointMoved(float currentSimulationTime)
{
	//
	// Recalculate state of all NPCs
	//

	for (auto & state : mStateBuffer)
	{
		if (state.has_value())
		{
			auto const oldRegime = state->CurrentRegime;

			ResetNpcStateToWorld(*state, currentSimulationTime);

			OnMayBeNpcRegimeChanged(oldRegime, *state);
		}
	}
}

bool Npcs::IsTriangleConstrainingCurrentlySelectedParticle(ElementIndex triangleIndex) const
{
#ifdef _DEBUG
	if (mCurrentlySelectedParticle.has_value())
	{
		for (auto const & state : mStateBuffer)
		{
			if (state.has_value())
			{
				if (std::any_of(
					state->ParticleMesh.Particles.cbegin(),
					state->ParticleMesh.Particles.cend(),
					[triangleIndex, selectedParticleIndex = *mCurrentlySelectedParticle](auto const & p)
					{
						return p.ParticleIndex == selectedParticleIndex
							&& p.ConstrainedState.has_value()
							&& p.ConstrainedState->CurrentBCoords.TriangleElementIndex == triangleIndex;
					}))
				{
					return true;
				}
			}
		}
	}
#else
	(void)triangleIndex;
#endif

	return false;
}

bool Npcs::IsSpringHostingCurrentlySelectedParticle(ElementIndex springIndex) const
{
#ifdef _DEBUG
	if (mCurrentlySelectedParticle.has_value())
	{
		for (auto const & ship : mShips)
		{
			if (ship.has_value())
			{
				for (auto const n : ship->Npcs)
				{
					auto & state = mStateBuffer[n];
					assert(state.has_value());

					for (auto const & particle : state->ParticleMesh.Particles)
					{
						if (particle.ParticleIndex == *mCurrentlySelectedParticle
							&& particle.ConstrainedState.has_value()
							&& particle.ConstrainedState->CurrentVirtualFloor.has_value()
							&& ship->ShipMesh.GetTriangles()
							.GetSubSprings(particle.ConstrainedState->CurrentVirtualFloor->TriangleElementIndex).SpringIndices[particle.ConstrainedState->CurrentVirtualFloor->EdgeOrdinal] == springIndex)
						{
							return true;
						}
					}
				}
			}
		}
	}
#else
	(void)springIndex;
#endif

	return false;
}

void Npcs::Publish() const
{
#ifdef _DEBUG
	std::optional<AbsoluteTriangleBCoords> constrainedRegimeParticleProbe;
	std::optional<int> constrainedRegimeLastEnteredFloorDepth;
	std::optional<bcoords3f> subjectParticleBarycentricCoordinatesWrtOriginTriangleChanged;
	std::optional<PhysicsParticleProbe> physicsParticleProbe;

	if (mCurrentlySelectedParticle.has_value())
	{
		for (size_t n = 0; n < mStateBuffer.size(); ++n)
		{
			if (mStateBuffer[n].has_value())
			{
				auto const & state = *mStateBuffer[n];

				assert(mShips[state.CurrentShipId].has_value());
				auto const & shipMesh = mShips[state.CurrentShipId]->ShipMesh;

				for (auto const & particle : state.ParticleMesh.Particles)
				{
					if (particle.ParticleIndex == *mCurrentlySelectedParticle)
					{
						if (particle.ConstrainedState.has_value())
						{
							constrainedRegimeParticleProbe.emplace(particle.ConstrainedState->CurrentBCoords);

							if (mCurrentOriginTriangle.has_value())
							{
								subjectParticleBarycentricCoordinatesWrtOriginTriangleChanged = shipMesh.GetTriangles().ToBarycentricCoordinates(
									mParticles.GetPosition(particle.ParticleIndex),
									*mCurrentOriginTriangle,
									shipMesh.GetPoints());
							}
						}

						physicsParticleProbe.emplace(mParticles.GetVelocity(particle.ParticleIndex));
					}
				}
			}
		}
	}

	mGameEventHandler->OnSubjectParticleConstrainedRegimeUpdated(constrainedRegimeParticleProbe);
	mGameEventHandler->OnSubjectParticleBarycentricCoordinatesWrtOriginTriangleChanged(subjectParticleBarycentricCoordinatesWrtOriginTriangleChanged);
	mGameEventHandler->OnSubjectParticlePhysicsUpdated(physicsParticleProbe);

	if (mCurrentlySelectedNpc.has_value()
		&& mStateBuffer[*mCurrentlySelectedNpc].has_value())
	{
		if (mStateBuffer[*mCurrentlySelectedNpc]->Kind == NpcKindType::Human)
		{
			switch (mStateBuffer[*mCurrentlySelectedNpc]->KindSpecificState.HumanNpcState.CurrentBehavior)
			{
				case StateType::KindSpecificStateType::HumanNpcStateType::BehaviorType::BeingPlaced:
				{
					mGameEventHandler->OnHumanNpcBehaviorChanged("BeingPlaced");
					break;
				}

				case StateType::KindSpecificStateType::HumanNpcStateType::BehaviorType::Constrained_Aerial:
				{
					mGameEventHandler->OnHumanNpcBehaviorChanged("Constrained_Aerial");
					break;
				}

				case StateType::KindSpecificStateType::HumanNpcStateType::BehaviorType::Constrained_Equilibrium:
				{
					mGameEventHandler->OnHumanNpcBehaviorChanged("Constrained_Equilibrium");
					break;
				}

				case StateType::KindSpecificStateType::HumanNpcStateType::BehaviorType::Constrained_Falling:
				{
					mGameEventHandler->OnHumanNpcBehaviorChanged("Constrained_Falling");
					break;
				}

				case StateType::KindSpecificStateType::HumanNpcStateType::BehaviorType::Constrained_KnockedOut:
				{
					mGameEventHandler->OnHumanNpcBehaviorChanged("Constrained_KnockedOut");
					break;
				}

				case StateType::KindSpecificStateType::HumanNpcStateType::BehaviorType::Constrained_PreRising:
				{
					mGameEventHandler->OnHumanNpcBehaviorChanged("Constrained_PreRising");
					break;
				}

				case StateType::KindSpecificStateType::HumanNpcStateType::BehaviorType::Constrained_Rising:
				{
					mGameEventHandler->OnHumanNpcBehaviorChanged("Constrained_Rising");
					break;
				}

				case StateType::KindSpecificStateType::HumanNpcStateType::BehaviorType::Constrained_Walking:
				{
					mGameEventHandler->OnHumanNpcBehaviorChanged("Constrained_Walking");
					break;
				}

				case StateType::KindSpecificStateType::HumanNpcStateType::BehaviorType::Free_Aerial:
				{
					mGameEventHandler->OnHumanNpcBehaviorChanged("Free_Aerial");
					break;
				}

				case StateType::KindSpecificStateType::HumanNpcStateType::BehaviorType::Free_InWater:
				{
					mGameEventHandler->OnHumanNpcBehaviorChanged("Free_InWater");
					break;
				}

				case StateType::KindSpecificStateType::HumanNpcStateType::BehaviorType::Free_Swimming_Style1:
				case StateType::KindSpecificStateType::HumanNpcStateType::BehaviorType::Free_Swimming_Style2:
				case StateType::KindSpecificStateType::HumanNpcStateType::BehaviorType::Free_Swimming_Style3:
				{
					mGameEventHandler->OnHumanNpcBehaviorChanged("Free_Swimming");
					break;
				}
			}
		}
	}
#endif
}

#endif

///////////////////////////////

NpcId Npcs::GetNewNpcId()
{
	// See if we can find a hole, so we stay compact
	for (size_t n = 0; n < mStateBuffer.size(); ++n)
	{
		if (!mStateBuffer[n].has_value())
		{
			return static_cast<NpcId>(n);
		}
	}

	// No luck, add new entry
	NpcId const newNpcId = static_cast<NpcId>(mStateBuffer.size());
	mStateBuffer.emplace_back(std::nullopt);
	return newNpcId;
}

ShipId Npcs::GetTopmostShipId() const
{
	assert(mShips.size() > 0);

	for (size_t s = mShips.size() - 1; ;)
	{
		if (mShips[s].has_value())
		{
			return static_cast<ShipId>(s);
		}

		if (s == 0)
		{
			break;
		}

		--s;
	}

	assert(false);
	return 0;
}

std::optional<ElementId> Npcs::FindTopmostTriangleContaining(vec2f const & position) const
{
	// Visit all ships in reverse ship ID order (i.e. from topmost to bottommost)
	assert(mShips.size() > 0);
	for (size_t s = mShips.size() - 1; ;)
	{
		if (mShips[s].has_value())
		{
			// Find the triangle in this ship containing this position and having the highest plane ID

			auto const & shipMesh = mShips[s]->ShipMesh;

			// TODO: this might be optimized

			std::optional<ElementIndex> bestTriangleIndex;
			PlaneId bestPlaneId = std::numeric_limits<PlaneId>::lowest();
			for (auto const triangleIndex : shipMesh.GetTriangles())
			{
				vec2f const aPosition = shipMesh.GetPoints().GetPosition(shipMesh.GetTriangles().GetPointAIndex(triangleIndex));
				vec2f const bPosition = shipMesh.GetPoints().GetPosition(shipMesh.GetTriangles().GetPointBIndex(triangleIndex));
				vec2f const cPosition = shipMesh.GetPoints().GetPosition(shipMesh.GetTriangles().GetPointCIndex(triangleIndex));

				if (IsPointInTriangle(position, aPosition, bPosition, cPosition)
					&& (!bestTriangleIndex || shipMesh.GetPoints().GetPlaneId(shipMesh.GetTriangles().GetPointAIndex(triangleIndex)) > bestPlaneId))
				{
					bestTriangleIndex = triangleIndex;
					bestPlaneId = shipMesh.GetPoints().GetPlaneId(shipMesh.GetTriangles().GetPointAIndex(triangleIndex));
				}
			}

			if (bestTriangleIndex)
			{
				// Found a triangle on this ship
				return ElementId(static_cast<ShipId>(s), *bestTriangleIndex);
			}
		}

		if (s == 0)
			break;
		--s;
	}

	// No triangle found
	return std::nullopt;
}

ElementIndex Npcs::FindTriangleContaining(
	vec2f const & position,
	Ship const & shipMesh)
{
	for (auto const triangleIndex : shipMesh.GetTriangles())
	{
		vec2f const aPosition = shipMesh.GetPoints().GetPosition(shipMesh.GetTriangles().GetPointAIndex(triangleIndex));
		vec2f const bPosition = shipMesh.GetPoints().GetPosition(shipMesh.GetTriangles().GetPointBIndex(triangleIndex));
		vec2f const cPosition = shipMesh.GetPoints().GetPosition(shipMesh.GetTriangles().GetPointCIndex(triangleIndex));

		if (IsPointInTriangle(position, aPosition, bPosition, cPosition))
		{
			return triangleIndex;
		}
	}

	return NoneElementIndex;
}

void Npcs::TransferNpcToShip(
	StateType & npc,
	ShipId newShip)
{
	if (npc.CurrentShipId == newShip)
	{
		return;
	}

	//
	// Maintain indices
	//

	assert(mShips[npc.CurrentShipId].has_value());

	auto it = std::find(
		mShips[npc.CurrentShipId]->Npcs.begin(),
		mShips[npc.CurrentShipId]->Npcs.end(),
		npc.Id);

	assert(it != mShips[npc.CurrentShipId]->Npcs.end());

	mShips[npc.CurrentShipId]->Npcs.erase(it);

	assert(mShips[newShip].has_value());

	mShips[newShip]->Npcs.push_back(newShip);

	//
	// Set ShipId in npc
	//

	npc.CurrentShipId = newShip;
}

void Npcs::RenderNpc(
	StateType const & npc,
	ShipRenderContext & shipRenderContext)
{
	assert(mShips[npc.CurrentShipId].has_value());

	auto const planeId = npc.CurrentPlaneId.has_value()
		? *(npc.CurrentPlaneId)
		: mShips[npc.CurrentShipId]->ShipMesh.GetMaxPlaneId();

	switch(npc.Kind)
	{
		case NpcKindType::Human:
		{
			assert(npc.ParticleMesh.Particles.size() == 2);
			auto const & humanNpcState = npc.KindSpecificState.HumanNpcState;
			auto const & animationState = humanNpcState.AnimationState;

			// Note:
			// - head, neck, shoulder, crotch, feet: based on current dipole length; anchor point is feet
			// - arms, legs : based on ideal (incl. adjustment)
			//

			vec2f const feetPosition = mParticles.GetPosition(npc.ParticleMesh.Particles[0].ParticleIndex);

			vec2f const actualBodyVector = mParticles.GetPosition(npc.ParticleMesh.Particles[1].ParticleIndex) - feetPosition; // From feet to head
			vec2f const actualBodyVDir = -actualBodyVector.normalise_approx(); // From head to feet - facilitates arm and length angle-making
			vec2f const actualBodyHDir = actualBodyVDir.to_perpendicular(); // Points R (of the screen)

			vec2f const crotchPosition = feetPosition + actualBodyVector * (GameParameters::HumanNpcGeometry::LegLengthFraction * animationState.LowerExtremityLengthMultiplier);
			vec2f const headPosition = crotchPosition + actualBodyVector * (GameParameters::HumanNpcGeometry::HeadLengthFraction + GameParameters::HumanNpcGeometry::TorsoLengthFraction);
			vec2f const neckPosition = headPosition - actualBodyVector * GameParameters::HumanNpcGeometry::HeadLengthFraction;
			vec2f const shoulderPosition = neckPosition - actualBodyVector * GameParameters::HumanNpcGeometry::ArmDepthFraction / 2.0f;

			// Arm and Leg lengths are relative to ideal
			float const adjustedIdealHumanHeight = humanNpcState.Height * mCurrentHumanNpcBodyLengthAdjustment;
			float const leftArmLength = adjustedIdealHumanHeight * GameParameters::HumanNpcGeometry::ArmLengthFraction * animationState.LimbLengthMultipliers.LeftArm;
			float const rightArmLength = adjustedIdealHumanHeight * GameParameters::HumanNpcGeometry::ArmLengthFraction * animationState.LimbLengthMultipliers.RightArm;
			float const leftLegLength = adjustedIdealHumanHeight * GameParameters::HumanNpcGeometry::LegLengthFraction * animationState.LimbLengthMultipliers.LeftLeg;
			float const rightLegLength = adjustedIdealHumanHeight * GameParameters::HumanNpcGeometry::LegLengthFraction * animationState.LimbLengthMultipliers.RightLeg;

			if (humanNpcState.CurrentFaceOrientation != 0.0f)
			{
				//
				// Front-back
				//

				float const halfHeadW = (adjustedIdealHumanHeight * GameParameters::HumanNpcGeometry::HeadWidthFraction * humanNpcState.WidthMultipier) / 2.0f;
				float const halfTorsoW = (adjustedIdealHumanHeight * GameParameters::HumanNpcGeometry::TorsoWidthFraction * humanNpcState.WidthMultipier) / 2.0f;
				float const halfArmW = (adjustedIdealHumanHeight * GameParameters::HumanNpcGeometry::ArmWidthFraction * humanNpcState.WidthMultipier) / 2.0f;
				float const halfLegW = (adjustedIdealHumanHeight * GameParameters::HumanNpcGeometry::LegWidthFraction * humanNpcState.WidthMultipier) / 2.0f;

				// Head

				shipRenderContext.UploadNpcQuad(
					planeId,
					TextureQuad(
						headPosition - actualBodyHDir * halfHeadW,
						vec2f(-1.0f, +1.0f),
						headPosition + actualBodyHDir * halfHeadW,
						vec2f(+1.0f, +1.0f),
						neckPosition - actualBodyHDir * halfHeadW,
						vec2f(-1.0f, -1.0f),
						neckPosition + actualBodyHDir * halfHeadW,
						vec2f(+1.0f, -1.0f)),
					humanNpcState.CurrentFaceOrientation,
					npc.Highlight);

				// Arms and legs

				vec2f const leftArmJointPosition = shoulderPosition - actualBodyHDir * halfTorsoW;
				vec2f const rightArmJointPosition = shoulderPosition + actualBodyHDir * halfTorsoW;

				vec2f const leftLegJointPosition = crotchPosition - actualBodyHDir * (halfTorsoW - halfLegW);
				vec2f const rightLegJointPosition = crotchPosition + actualBodyHDir * (halfTorsoW - halfLegW);

				if (humanNpcState.CurrentFaceOrientation > 0.0f)
				{
					// Front

					// Left arm (on left side of the screen)
					vec2f const leftArmDir = actualBodyVDir.rotate(animationState.LimbAnglesCos.LeftArm, animationState.LimbAnglesSin.LeftArm);
					vec2f const leftArmVector = leftArmDir * leftArmLength;
					vec2f const leftArmTraverseVector = leftArmDir.to_perpendicular() * halfArmW;
					shipRenderContext.UploadNpcQuad(
						planeId,
						TextureQuad(
							leftArmJointPosition - leftArmTraverseVector,
							vec2f(-1.0f, +1.0f),
							leftArmJointPosition + leftArmTraverseVector,
							vec2f(+1.0f, +1.0f),
							leftArmJointPosition + leftArmVector - leftArmTraverseVector,
							vec2f(-1.0f, -1.0f),
							leftArmJointPosition + leftArmVector + leftArmTraverseVector,
							vec2f(+1.0f, -1.0f)),
						humanNpcState.CurrentFaceOrientation,
						npc.Highlight);

					// Right arm (on right side of the screen)
					vec2f const rightArmDir = actualBodyVDir.rotate(animationState.LimbAnglesCos.RightArm, animationState.LimbAnglesSin.RightArm);
					vec2f const rightArmVector = rightArmDir * rightArmLength;
					vec2f const rightArmTraverseVector = rightArmDir.to_perpendicular() * halfArmW;
					shipRenderContext.UploadNpcQuad(
						planeId,
						TextureQuad(
							rightArmJointPosition - rightArmTraverseVector,
							vec2f(-1.0f, +1.0f),
							rightArmJointPosition + rightArmTraverseVector,
							vec2f(+1.0f, +1.0f),
							rightArmJointPosition + rightArmVector - rightArmTraverseVector,
							vec2f(-1.0f, -1.0f),
							rightArmJointPosition + rightArmVector + rightArmTraverseVector,
							vec2f(+1.0f, -1.0f)),
						humanNpcState.CurrentFaceOrientation,
						npc.Highlight);

					// Left leg (on left side of the screen)
					vec2f const leftLegDir = actualBodyVDir.rotate(animationState.LimbAnglesCos.LeftLeg, animationState.LimbAnglesSin.LeftLeg);
					vec2f const leftLegVector = leftLegDir * leftLegLength;
					vec2f const leftLegTraverseVector = leftLegDir.to_perpendicular() * halfLegW;
					shipRenderContext.UploadNpcQuad(
						planeId,
						TextureQuad(
							leftLegJointPosition - leftLegTraverseVector,
							vec2f(-1.0f, +1.0f),
							leftLegJointPosition + leftLegTraverseVector,
							vec2f(+1.0f, +1.0f),
							leftLegJointPosition + leftLegVector - leftLegTraverseVector,
							vec2f(-1.0f, -1.0f),
							leftLegJointPosition + leftLegVector + leftLegTraverseVector,
							vec2f(+1.0f, -1.0f)),
						humanNpcState.CurrentFaceOrientation,
						npc.Highlight);

					// Right leg (on right side of the screen)
					vec2f const rightLegDir = actualBodyVDir.rotate(animationState.LimbAnglesCos.RightLeg, animationState.LimbAnglesSin.RightLeg);
					vec2f const rightLegVector = rightLegDir * rightLegLength;
					vec2f const rightLegTraverseVector = rightLegDir.to_perpendicular() * halfLegW;
					shipRenderContext.UploadNpcQuad(
						planeId,
						TextureQuad(
							rightLegJointPosition - rightLegTraverseVector,
							vec2f(-1.0f, +1.0f),
							rightLegJointPosition + rightLegTraverseVector,
							vec2f(+1.0f, +1.0f),
							rightLegJointPosition + rightLegVector - rightLegTraverseVector,
							vec2f(-1.0f, -1.0f),
							rightLegJointPosition + rightLegVector + rightLegTraverseVector,
							vec2f(+1.0f, -1.0f)),
						humanNpcState.CurrentFaceOrientation,
						npc.Highlight);
				}
				else
				{
					// Back

					// Left arm (on right side of screen)
					vec2f const leftArmDir = actualBodyVDir.rotate(animationState.LimbAnglesCos.LeftArm, -animationState.LimbAnglesSin.LeftArm);
					vec2f const leftArmVector = leftArmDir * leftArmLength;
					vec2f const leftArmTraverseVector = leftArmDir.to_perpendicular() * halfArmW;
					shipRenderContext.UploadNpcQuad(
						planeId,
						TextureQuad(
							rightArmJointPosition - leftArmTraverseVector,
							vec2f(-1.0f, +1.0f),
							rightArmJointPosition + leftArmTraverseVector,
							vec2f(+1.0f, +1.0f),
							rightArmJointPosition + leftArmVector - leftArmTraverseVector,
							vec2f(-1.0f, -1.0f),
							rightArmJointPosition + leftArmVector + leftArmTraverseVector,
							vec2f(+1.0f, -1.0f)),
						humanNpcState.CurrentFaceOrientation,
						npc.Highlight);

					// Right arm (on left side of the screen)
					vec2f const rightArmDir = actualBodyVDir.rotate(animationState.LimbAnglesCos.RightArm, -animationState.LimbAnglesSin.RightArm);
					vec2f const rightArmVector = rightArmDir * rightArmLength;
					vec2f const rightArmTraverseVector = rightArmDir.to_perpendicular() * halfArmW;
					shipRenderContext.UploadNpcQuad(
						planeId,
						TextureQuad(
							leftArmJointPosition - rightArmTraverseVector,
							vec2f(-1.0f, +1.0f),
							leftArmJointPosition + rightArmTraverseVector,
							vec2f(+1.0f, +1.0f),
							leftArmJointPosition + rightArmVector - rightArmTraverseVector,
							vec2f(-1.0f, -1.0f),
							leftArmJointPosition + rightArmVector + rightArmTraverseVector,
							vec2f(+1.0f, -1.0f)),
						humanNpcState.CurrentFaceOrientation,
						npc.Highlight);

					// Left leg (on right side of the screen)
					vec2f const leftLegDir = actualBodyVDir.rotate(animationState.LimbAnglesCos.LeftLeg, -animationState.LimbAnglesSin.LeftLeg);
					vec2f const leftLegVector = leftLegDir * leftLegLength;
					vec2f const leftLegTraverseVector = leftLegDir.to_perpendicular() * halfLegW;
					shipRenderContext.UploadNpcQuad(
						planeId,
						TextureQuad(
							rightLegJointPosition - leftLegTraverseVector,
							vec2f(-1.0f, +1.0f),
							rightLegJointPosition + leftLegTraverseVector,
							vec2f(+1.0f, +1.0f),
							rightLegJointPosition + leftLegVector - leftLegTraverseVector,
							vec2f(-1.0f, -1.0f),
							rightLegJointPosition + leftLegVector + leftLegTraverseVector,
							vec2f(+1.0f, -1.0f)),
						humanNpcState.CurrentFaceOrientation,
						npc.Highlight);

					// Right leg (on left side of the screen)
					vec2f const rightLegDir = actualBodyVDir.rotate(animationState.LimbAnglesCos.RightLeg, -animationState.LimbAnglesSin.RightLeg);
					vec2f const rightLegVector = rightLegDir * rightLegLength;
					vec2f const rightLegTraverseVector = rightLegDir.to_perpendicular() * halfLegW;
					shipRenderContext.UploadNpcQuad(
						planeId,
						TextureQuad(
							leftLegJointPosition - rightLegTraverseVector,
							vec2f(-1.0f, +1.0f),
							leftLegJointPosition + rightLegTraverseVector,
							vec2f(+1.0f, +1.0f),
							leftLegJointPosition + rightLegVector - rightLegTraverseVector,
							vec2f(-1.0f, -1.0f),
							leftLegJointPosition + rightLegVector + rightLegTraverseVector,
							vec2f(+1.0f, -1.0f)),
						humanNpcState.CurrentFaceOrientation,
						npc.Highlight);
				}

				// Torso

				shipRenderContext.UploadNpcQuad(
					planeId,
					TextureQuad(
						neckPosition - actualBodyHDir * halfTorsoW,
						vec2f(-1.0f, +1.0f),
						neckPosition + actualBodyHDir * halfTorsoW,
						vec2f(+1.0f, +1.0f),
						crotchPosition - actualBodyHDir * halfTorsoW,
						vec2f(-1.0f, -1.0f),
						crotchPosition + actualBodyHDir * halfTorsoW,
						vec2f(+1.0f, -1.0f)),
					humanNpcState.CurrentFaceOrientation,
					npc.Highlight);
			}
			else
			{
				//
				// Left-Right
				//

				float const halfHeadD = (adjustedIdealHumanHeight * GameParameters::HumanNpcGeometry::HeadDepthFraction * humanNpcState.WidthMultipier) / 2.0f;
				float const halfTorsoD = (adjustedIdealHumanHeight * GameParameters::HumanNpcGeometry::TorsoDepthFraction * humanNpcState.WidthMultipier) / 2.0f;
				float const halfArmD = (adjustedIdealHumanHeight * GameParameters::HumanNpcGeometry::ArmDepthFraction * humanNpcState.WidthMultipier) / 2.0f;
				float const halfLegD = (adjustedIdealHumanHeight * GameParameters::HumanNpcGeometry::LegDepthFraction * humanNpcState.WidthMultipier) / 2.0f;

				vec2f topLeftTexture(-humanNpcState.CurrentFaceDirectionX, +1.0f);
				vec2f topRightTexture(humanNpcState.CurrentFaceDirectionX, +1.0f);
				vec2f bottomLeftTexture(-humanNpcState.CurrentFaceDirectionX, -1.0f);
				vec2f bottomRightTexture(humanNpcState.CurrentFaceDirectionX, -1.0f);

				// Note: angles are with body vertical, regardless of L/R

				vec2f const leftArmDir = actualBodyVDir.rotate(animationState.LimbAnglesCos.LeftArm, animationState.LimbAnglesSin.LeftArm);
				vec2f const leftArmVector = leftArmDir * leftArmLength;
				vec2f const leftArmTraverseVector = leftArmDir.to_perpendicular() * halfArmD;
				TextureQuad leftArmQuad(
					shoulderPosition - leftArmTraverseVector,
					topLeftTexture,
					shoulderPosition + leftArmTraverseVector,
					topRightTexture,
					shoulderPosition + leftArmVector - leftArmTraverseVector,
					bottomLeftTexture,
					shoulderPosition + leftArmVector + leftArmTraverseVector,
					bottomRightTexture);

				vec2f const rightArmDir = actualBodyVDir.rotate(animationState.LimbAnglesCos.RightArm, animationState.LimbAnglesSin.RightArm);
				vec2f const rightArmVector = rightArmDir * rightArmLength;
				vec2f const rightArmTraverseVector = rightArmDir.to_perpendicular() * halfArmD;
				TextureQuad rightArmQuad(
					shoulderPosition - rightArmTraverseVector,
					topLeftTexture,
					shoulderPosition + rightArmTraverseVector,
					topRightTexture,
					shoulderPosition + rightArmVector - rightArmTraverseVector,
					bottomLeftTexture,
					shoulderPosition + rightArmVector + rightArmTraverseVector,
					bottomRightTexture);

				vec2f const leftUpperLegDir = actualBodyVDir.rotate(animationState.LimbAnglesCos.LeftLeg, animationState.LimbAnglesSin.LeftLeg);
				vec2f const leftUpperLegVector = leftUpperLegDir * leftLegLength * animationState.UpperLegLengthFraction;
				vec2f const leftUpperLegTraverseDir = leftUpperLegDir.to_perpendicular();
				vec2f const leftKneeOrFootPosition = crotchPosition + leftUpperLegVector; // When UpperLegLengthFraction is 1.0 (whole leg), this is the (virtual) foot
				TextureQuad leftUpperLegQuad;
				std::optional<TextureQuad> leftLowerLegQuad;

				vec2f const rightUpperLegDir = actualBodyVDir.rotate(animationState.LimbAnglesCos.RightLeg, animationState.LimbAnglesSin.RightLeg);
				vec2f const rightUpperLegVector = rightUpperLegDir * rightLegLength * animationState.UpperLegLengthFraction;
				vec2f const rightUpperLegTraverseDir = rightUpperLegDir.to_perpendicular();
				vec2f const rightKneeOrFootPosition = crotchPosition + rightUpperLegVector; // When UpperLegLengthFraction is 1.0 (whole leg), this is the (virtual) foot
				TextureQuad rightUpperLegQuad;
				std::optional<TextureQuad> rightLowerLegQuad;

				float const lowerLegLengthFraction = 1.0f - animationState.UpperLegLengthFraction;
				if (lowerLegLengthFraction != 0.0f)
				{
					//
					// Both upper and lower legs
					//

					float const kneeTextureY = 1.0f - 2.0f * animationState.UpperLegLengthFraction;

					// We extrude the corners to make them join nicely to the previous
					// and next segments. The calculation of the extrusion (J) between two
					// segments is based on these observations:
					//  * The direction of the extrusion is along the resultant of the normals
					//    to the two segments
					//  * The magnitude of the extrusion is (W/2) / cos(alpha), where alpha is
					//    the angle between a normal and the direction of the extrusion

					float constexpr MinJ = 0.9f;

					vec2f const leftLowerLegDir = (feetPosition - leftKneeOrFootPosition).normalise_approx();
					vec2f const leftLowerLegVector = leftLowerLegDir * leftLegLength * lowerLegLengthFraction;
					vec2f const leftLowerLegTraverseDir = leftLowerLegDir.to_perpendicular();
					vec2f const leftLegResultantNormal = leftUpperLegTraverseDir + leftLowerLegTraverseDir;
					vec2f const leftLegJ = leftLegResultantNormal / std::max(MinJ, leftUpperLegTraverseDir.dot(leftLegResultantNormal)) * halfLegD;

					leftUpperLegQuad = TextureQuad(
						crotchPosition - leftUpperLegTraverseDir * halfLegD,
						topLeftTexture,
						crotchPosition + leftUpperLegTraverseDir * halfLegD,
						topRightTexture,
						leftKneeOrFootPosition - leftLegJ,
						vec2f(bottomLeftTexture.x, kneeTextureY),
						leftKneeOrFootPosition + leftLegJ,
						vec2f(bottomRightTexture.x, kneeTextureY));

					leftLowerLegQuad = TextureQuad(
						leftKneeOrFootPosition - leftLegJ,
						vec2f(bottomLeftTexture.x, kneeTextureY),
						leftKneeOrFootPosition + leftLegJ,
						vec2f(bottomRightTexture.x, kneeTextureY),
						leftKneeOrFootPosition + leftLowerLegVector - leftLowerLegTraverseDir * halfLegD,
						bottomLeftTexture,
						leftKneeOrFootPosition + leftLowerLegVector + leftLowerLegTraverseDir * halfLegD,
						bottomRightTexture);

					vec2f const rightLowerLegDir = (feetPosition - rightKneeOrFootPosition).normalise_approx();
					vec2f const rightLowerLegVector = rightLowerLegDir * rightLegLength * lowerLegLengthFraction;
					vec2f const rightLowerLegTraverseDir = rightLowerLegDir.to_perpendicular();
					vec2f const rightLegResultantNormal = rightUpperLegTraverseDir + rightLowerLegTraverseDir;
					vec2f const rightLegJ = rightLegResultantNormal / std::max(MinJ, rightUpperLegTraverseDir.dot(rightLegResultantNormal)) * halfLegD;

					rightUpperLegQuad = TextureQuad(
						crotchPosition - rightUpperLegTraverseDir * halfLegD,
						topLeftTexture,
						crotchPosition + rightUpperLegTraverseDir * halfLegD,
						topRightTexture,
						rightKneeOrFootPosition - rightLegJ,
						vec2f(bottomLeftTexture.x, kneeTextureY),
						rightKneeOrFootPosition + rightLegJ,
						vec2f(bottomRightTexture.x, kneeTextureY));

					rightLowerLegQuad = TextureQuad(
						rightKneeOrFootPosition - rightLegJ,
						vec2f(bottomLeftTexture.x, kneeTextureY),
						rightKneeOrFootPosition + rightLegJ,
						vec2f(bottomRightTexture.x, kneeTextureY),
						rightKneeOrFootPosition + rightLowerLegVector - rightLowerLegTraverseDir * halfLegD,
						bottomLeftTexture,
						rightKneeOrFootPosition + rightLowerLegVector + rightLowerLegTraverseDir * halfLegD,
						bottomRightTexture);
				}
				else
				{
					// Just upper leg

					leftUpperLegQuad = TextureQuad(
						crotchPosition - leftUpperLegTraverseDir * halfLegD,
						topLeftTexture,
						crotchPosition + leftUpperLegTraverseDir * halfLegD,
						topRightTexture,
						leftKneeOrFootPosition - leftUpperLegTraverseDir * halfLegD,
						bottomLeftTexture,
						leftKneeOrFootPosition + leftUpperLegTraverseDir * halfLegD,
						bottomRightTexture);

					rightUpperLegQuad = TextureQuad(
						crotchPosition - rightUpperLegTraverseDir * halfLegD,
						topLeftTexture,
						crotchPosition + rightUpperLegTraverseDir * halfLegD,
						topRightTexture,
						rightKneeOrFootPosition - rightUpperLegTraverseDir * halfLegD,
						bottomLeftTexture,
						rightKneeOrFootPosition + rightUpperLegTraverseDir * halfLegD,
						bottomRightTexture);
				}

				// Arm and leg far

				if (humanNpcState.CurrentFaceDirectionX > 0.0f)
				{
					// Left arm
					shipRenderContext.UploadNpcQuad(
						planeId,
						leftArmQuad,
						humanNpcState.CurrentFaceOrientation,
						npc.Highlight);

					// Left leg
					shipRenderContext.UploadNpcQuad(
						planeId,
						leftUpperLegQuad,
						humanNpcState.CurrentFaceOrientation,
						npc.Highlight);
					if (leftLowerLegQuad.has_value())
					{
						shipRenderContext.UploadNpcQuad(
							planeId,
							*leftLowerLegQuad,
							humanNpcState.CurrentFaceOrientation,
							npc.Highlight);
					}
				}
				else
				{
					// Right arm
					shipRenderContext.UploadNpcQuad(
						planeId,
						rightArmQuad,
						humanNpcState.CurrentFaceOrientation,
						npc.Highlight);

					// Right leg
					shipRenderContext.UploadNpcQuad(
						planeId,
						rightUpperLegQuad,
						humanNpcState.CurrentFaceOrientation,
						npc.Highlight);
					if (rightLowerLegQuad.has_value())
					{
						shipRenderContext.UploadNpcQuad(
							planeId,
							*rightLowerLegQuad,
							humanNpcState.CurrentFaceOrientation,
							npc.Highlight);
					}
				}

				// Head

				shipRenderContext.UploadNpcQuad(
					planeId,
					TextureQuad(
						headPosition - actualBodyHDir * halfHeadD,
						topLeftTexture,
						headPosition + actualBodyHDir * halfHeadD,
						topRightTexture,
						neckPosition - actualBodyHDir * halfHeadD,
						bottomLeftTexture,
						neckPosition + actualBodyHDir * halfHeadD,
						bottomRightTexture),
					humanNpcState.CurrentFaceOrientation,
					npc.Highlight);

				// Torso

				shipRenderContext.UploadNpcQuad(
					planeId,
					TextureQuad(
						neckPosition - actualBodyHDir * halfTorsoD,
						topLeftTexture,
						neckPosition + actualBodyHDir * halfTorsoD,
						topRightTexture,
						crotchPosition - actualBodyHDir * halfTorsoD,
						bottomLeftTexture,
						crotchPosition + actualBodyHDir * halfTorsoD,
						bottomRightTexture),
					humanNpcState.CurrentFaceOrientation,
					npc.Highlight);

				// Arm and leg near

				if (humanNpcState.CurrentFaceDirectionX > 0.0f)
				{
					// Right arm
					shipRenderContext.UploadNpcQuad(
						planeId,
						rightArmQuad,
						humanNpcState.CurrentFaceOrientation,
						npc.Highlight);

					// Right leg
					shipRenderContext.UploadNpcQuad(
						planeId,
						rightUpperLegQuad,
						humanNpcState.CurrentFaceOrientation,
						npc.Highlight);
					if (rightLowerLegQuad.has_value())
					{
						shipRenderContext.UploadNpcQuad(
							planeId,
							*rightLowerLegQuad,
							humanNpcState.CurrentFaceOrientation,
							npc.Highlight);
					}
				}
				else
				{
					// Left arm
					shipRenderContext.UploadNpcQuad(
						planeId,
						leftArmQuad,
						humanNpcState.CurrentFaceOrientation,
						npc.Highlight);

					// Left leg
					shipRenderContext.UploadNpcQuad(
						planeId,
						leftUpperLegQuad,
						humanNpcState.CurrentFaceOrientation,
						npc.Highlight);
					if (leftLowerLegQuad.has_value())
					{
						shipRenderContext.UploadNpcQuad(
							planeId,
							*leftLowerLegQuad,
							humanNpcState.CurrentFaceOrientation,
							npc.Highlight);
					}
				}
			}

			break;
		}

		case NpcKindType::Furniture:
		{
			// Futurework
			if (npc.ParticleMesh.Particles.size() == 4)
			{
				// Quad
				shipRenderContext.UploadNpcQuad(
					planeId,
					TextureQuad(
						mParticles.GetPosition(npc.ParticleMesh.Particles[0].ParticleIndex),
						vec2f{-1.0f, 1.0f},
						mParticles.GetPosition(npc.ParticleMesh.Particles[1].ParticleIndex),
						vec2f{ 1.0f, 1.0f },
						mParticles.GetPosition(npc.ParticleMesh.Particles[3].ParticleIndex),
						vec2f{ -1.0f, -1.0f },
						mParticles.GetPosition(npc.ParticleMesh.Particles[2].ParticleIndex),
						vec2f{ 1.0f, -1.0f }
					),
					1.0f,
					npc.Highlight);
			}
			else
			{
				// Bunch-of-particles

				for (auto const & particle : npc.ParticleMesh.Particles)
				{
					shipRenderContext.UploadNpcParticle(
						planeId,
						mParticles.GetPosition(particle.ParticleIndex),
						mParticles.GetRenderColor(particle.ParticleIndex),
						1.0f,
						npc.Highlight);
				}
			}

			break;
		}
	}
}

void Npcs::PublishHumanNpcStats()
{
	// TODOTEST
	////mGameEventHandler->OnHumanNpcCountsUpdated(
	////	mConstrainedRegimeHumanNpcCount,
	////	mFreeRegimeHumanNpcCount);
}

}