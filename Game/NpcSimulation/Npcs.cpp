/***************************************************************************************
* Original Author:      Gabriele Giuseppini
* Created:              2023-10-06
* Copyright:            Gabriele Giuseppini  (https://github.com/GabrieleGiuseppini)
***************************************************************************************/
#include "Physics.h"

#include <GameCore/Colors.h>

#include <cassert>

namespace Physics {

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
	}

	//
	// Update NPCs' state
	//

	UpdateNpcs(currentSimulationTime, gameParameters);

	//
	// Publish
	//

	Publish();
}

void Npcs::Upload(RenderContext & renderContext)
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

					shipRenderContext.UploadNpcParticle(
						mParticles.GetPosition(state.PrimaryParticleState.ParticleIndex),
						mParticles.GetRenderColor(state.PrimaryParticleState.ParticleIndex),
						1.0f);

					if (state.DipoleState.has_value())
					{
						shipRenderContext.UploadNpcParticle(
							mParticles.GetPosition(state.DipoleState->SecondaryParticleState.ParticleIndex),
							mParticles.GetRenderColor(state.DipoleState->SecondaryParticleState.ParticleIndex),
							1.0f);

						renderContext.UploadNpcSpring(
							mParticles.GetPosition(state.PrimaryParticleState.ParticleIndex),
							mParticles.GetPosition(state.DipoleState->SecondaryParticleState.ParticleIndex),
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

	renderContext.UploadNpcQuadsStart();

#ifdef IN_BARYLAB
	// For furniture
	renderContext.UploadNpcParticlesStart();
#endif

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
	mShips[s].emplace(
		ShipMeshType(
			ship.GetPoints(),
			ship.GetSprings(),
			ship.GetTriangles()));
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

	for (auto const npcId : mShips[s]->Npcs)
	{
		OnNpcDestroyed(npcId);
		mStateBuffer[npcId].reset();
	}

	mGameEventHandler->OnNpcCountsUpdated(
		mNpcCount,
		mConstrainedRegimeHumanNpcCount,
		mFreeRegimeHumanNpcCount,
		GameParameters::MaxNpcs - mNpcCount);

	//
	// Destroy NPC ship
	//

	mShips[s].reset();
}

std::optional<PickedObjectId<NpcId>> Npcs::BeginPlaceNewHumanNpc(
	HumanNpcKindType humanKind,
	vec2f const & worldCoordinates,
	float currentSimulationTime)
{
	//
	// Check if there are enough NPCs and particles
	//

	if (mNpcCount >= GameParameters::MaxNpcs || mParticles.GetRemainingParticlesCount() < 2)
	{
		return std::nullopt;
	}

	//
	// Create NPC
	//

	// Feet (primary)

	auto const & feetMaterial = mMaterialDatabase.GetNpcMaterial(NpcMaterial::KindType::HumanFeet);

	auto const primaryParticleIndex = mParticles.Add(
		feetMaterial.Mass,
		feetMaterial.StaticFriction,
		feetMaterial.KineticFriction,
		feetMaterial.Elasticity,
		feetMaterial.BuoyancyVolumeFill,
		worldCoordinates - vec2f(0.0f, 1.0f) * GameParameters::HumanNpcGeometry::BodyLength * mCurrentHumanNpcBodyLengthAdjustment,
		feetMaterial.RenderColor);

	StateType::NpcParticleStateType primaryParticleState = StateType::NpcParticleStateType(
		primaryParticleIndex,
		std::nullopt);

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

	StateType::NpcParticleStateType secondaryParticleState = StateType::NpcParticleStateType(
		secondaryParticleIndex,
		std::nullopt);

	// Dipole

	float const massFactor =
		(feetMaterial.Mass * headMaterial.Mass)
		/ (feetMaterial.Mass + headMaterial.Mass);

	StateType::DipoleStateType dipoleState = StateType::DipoleStateType(
		std::move(secondaryParticleState),
		StateType::DipolePropertiesType(
			GameParameters::HumanNpcGeometry::BodyLength,
			massFactor,
			1.0f));

	// Human

	StateType::KindSpecificStateType::HumanNpcStateType humanState = StateType::KindSpecificStateType::HumanNpcStateType(
		humanKind,
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
	ShipId const shipId = FindTopmostShipId();

	mStateBuffer[npcId].emplace(
		npcId,
		NpcKindType::Human,
		shipId,
		StateType::RegimeType::BeingPlaced,
		std::move(primaryParticleState),
		std::move(dipoleState),
		StateType::KindSpecificStateType(std::move(humanState)));

	assert(mShips[shipId].has_value());
	mShips[shipId]->Npcs.push_back(npcId);

	OnNpcCreated(npcId);

	return PickedObjectId<NpcId>(npcId, vec2f::zero());
}

void Npcs::MoveNpcTo(
	NpcId id,
	vec2f const & position,
	vec2f const & offset)
{
	assert(mStateBuffer[id].has_value());
	assert(mStateBuffer[id]->CurrentRegime == StateType::RegimeType::BeingPlaced);

	vec2f const newPosition = position - offset;

	switch (mStateBuffer[id]->Kind)
	{
		case NpcKindType::Furniture:
		{
			// Act on primary particle

			auto const particleIndex = mStateBuffer[id]->PrimaryParticleState.ParticleIndex;
			vec2f const oldPosition = mParticles.GetPosition(particleIndex);
			mParticles.SetPosition(particleIndex, newPosition);
			mParticles.SetVelocity(particleIndex, (newPosition - oldPosition) / GameParameters::SimulationTimeStepDuration);

			if (mStateBuffer[id]->PrimaryParticleState.ConstrainedState.has_value())
			{
				// TODO
				mStateBuffer[id]->PrimaryParticleState.ConstrainedState->MeshRelativeVelocity = vec2f::zero();
			}

			break;
		}

		case NpcKindType::Human:
		{
			// Act on secondary particle

			assert(mStateBuffer[id]->DipoleState.has_value());

			auto const particleIndex = mStateBuffer[id]->DipoleState->SecondaryParticleState.ParticleIndex;
			vec2f const oldPosition = mParticles.GetPosition(particleIndex);
			mParticles.SetPosition(particleIndex, newPosition);
			mParticles.SetVelocity(particleIndex, (newPosition - oldPosition) / GameParameters::SimulationTimeStepDuration);

			if (mStateBuffer[id]->DipoleState->SecondaryParticleState.ConstrainedState.has_value())
			{
				// TODO
				mStateBuffer[id]->DipoleState->SecondaryParticleState.ConstrainedState->MeshRelativeVelocity = vec2f::zero();
			}

			mStateBuffer[id]->KindSpecificState.HumanNpcState.TotalDistanceTraveledOffEdgeSinceStateTransition += (newPosition - oldPosition).length();

			break;
		}
	}
}

void Npcs::SetPanicLevelForAllHumans(float panicLevel)
{
	for (auto & npc : mStateBuffer)
	{
		if (npc.has_value())
		{
			if (npc->Kind == NpcKindType::Human)
			{
				npc->KindSpecificState.HumanNpcState.PanicLevel = panicLevel;
			}
		}
	}
}

/////////////////////////////// Barylab-specific

#ifdef IN_BARYLAB

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

	for (auto const & ship : mShips)
	{
		if (ship.has_value())
		{
			for (auto const n : ship->Npcs)
			{
				auto & state = mStateBuffer[n];
				assert(state.has_value());

				if (state->PrimaryParticleState.ParticleIndex == particleIndex
					|| (state->DipoleState.has_value() && state->DipoleState->SecondaryParticleState.ParticleIndex == particleIndex))
				{
					ResetNpcStateToWorld(*state, currentSimulationTime, ship->ShipMesh);

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

				RotateParticleWithShip(
					state->PrimaryParticleState,
					centerPos,
					cosAngle,
					sinAngle,
					ship->ShipMesh);

				if (state->DipoleState.has_value())
				{
					RotateParticleWithShip(
						state->DipoleState->SecondaryParticleState,
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
	ShipMeshType const & shipMesh)
{
	vec2f newPosition;

	if (npcParticleState.ConstrainedState.has_value())
	{
		// Simply set position from current bary coords

		newPosition = shipMesh.ShipTriangles.FromBarycentricCoordinates(
			npcParticleState.ConstrainedState->CurrentTriangleBarycentricCoords,
			npcParticleState.ConstrainedState->CurrentTriangle,
			shipMesh.ShipPoints);
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

	for (auto const & ship : mShips)
	{
		if (ship.has_value())
		{
			for (auto const n : ship->Npcs)
			{
				auto & state = mStateBuffer[n];
				assert(state.has_value());

				ResetNpcStateToWorld(*state, currentSimulationTime, ship->ShipMesh);
			}
		}
	}
}

bool Npcs::IsTriangleConstrainingCurrentlySelectedParticle(ElementIndex triangleIndex) const
{
	if (mCurrentlySelectedParticle.has_value())
	{
		for (auto const & state : mStateBuffer)
		{
			if (state.has_value())
			{
				if (state->PrimaryParticleState.ParticleIndex == *mCurrentlySelectedParticle
					&& state->PrimaryParticleState.ConstrainedState.has_value()
					&& triangleIndex == state->PrimaryParticleState.ConstrainedState->CurrentTriangle)
				{
					return true;
				}

				if (state->DipoleState.has_value()
					&& state->DipoleState->SecondaryParticleState.ParticleIndex == *mCurrentlySelectedParticle
					&& state->DipoleState->SecondaryParticleState.ConstrainedState.has_value()
					&& triangleIndex == state->DipoleState->SecondaryParticleState.ConstrainedState->CurrentTriangle)
				{
					return true;
				}
			}
		}
	}

	return false;
}

bool Npcs::IsSpringHostingCurrentlySelectedParticle(ElementIndex springIndex) const
{
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

					if (state->PrimaryParticleState.ParticleIndex == *mCurrentlySelectedParticle
						|| (state->DipoleState.has_value() && state->DipoleState->SecondaryParticleState.ParticleIndex == *mCurrentlySelectedParticle))
					{
						if (state->PrimaryParticleState.ConstrainedState.has_value()
							&& state->PrimaryParticleState.ConstrainedState->CurrentVirtualEdgeOrdinal >= 0
							&& springIndex == ship->ShipMesh.ShipTriangles.GetSubSprings(state->PrimaryParticleState.ConstrainedState->CurrentTriangle).SpringIndices[state->PrimaryParticleState.ConstrainedState->CurrentVirtualEdgeOrdinal])
						{
							return true;
						}

						if (state->DipoleState.has_value()
							&& state->DipoleState->SecondaryParticleState.ConstrainedState.has_value()
							&& state->DipoleState->SecondaryParticleState.ConstrainedState->CurrentVirtualEdgeOrdinal >= 0
							&& springIndex == ship->ShipMesh.ShipTriangles.GetSubSprings(state->DipoleState->SecondaryParticleState.ConstrainedState->CurrentTriangle).SpringIndices[state->DipoleState->SecondaryParticleState.ConstrainedState->CurrentVirtualEdgeOrdinal])
						{
							return true;
						}
					}
				}
			}
		}
	}

	return false;
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

ShipId Npcs::FindTopmostShipId() const
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

void Npcs::OnNpcCreated(NpcId id)
{
	assert(mStateBuffer[id].has_value());
	auto const & state = *mStateBuffer[id];

	//
	// Update stats
	//

	if (state.CurrentRegime == StateType::RegimeType::Constrained)
	{
		++mConstrainedRegimeHumanNpcCount;
	}
	else if (state.CurrentRegime == StateType::RegimeType::Free)
	{
		++mFreeRegimeHumanNpcCount;
	}

	++mNpcCount;

	mGameEventHandler->OnNpcCountsUpdated(
		mNpcCount,
		mConstrainedRegimeHumanNpcCount,
		mFreeRegimeHumanNpcCount,
		GameParameters::MaxNpcs - mNpcCount);
}

void Npcs::OnNpcDestroyed(NpcId id)
{
	// We are invoked _before_ actual destruction happens
	assert(mStateBuffer[id].has_value());
	auto const & state = *mStateBuffer[id];

	//
	// Update stats
	//

	if (state.CurrentRegime == StateType::RegimeType::Constrained)
	{
		assert(mConstrainedRegimeHumanNpcCount > 0);
		--mConstrainedRegimeHumanNpcCount;
	}
	else if (state.CurrentRegime == StateType::RegimeType::Free)
	{
		assert(mFreeRegimeHumanNpcCount > 0);
		--mFreeRegimeHumanNpcCount;
	}

	--mNpcCount;

	mGameEventHandler->OnNpcCountsUpdated(
		mNpcCount,
		mConstrainedRegimeHumanNpcCount,
		mFreeRegimeHumanNpcCount,
		GameParameters::MaxNpcs - mNpcCount);
}

void Npcs::RenderNpc(
	StateType const & npc,
	ShipRenderContext & shipRenderContext)
{
	switch(npc.Kind)
	{
		case NpcKindType::Human:
		{
			assert(npc.DipoleState.has_value());
			auto const & humanNpcState = npc.KindSpecificState.HumanNpcState;

			// Note:
			// - head, neck, shoulder, crotch, feet: based on current dipole length
			// - arms, legs : based on ideal (incl. adjustment)
			//

			vec2f const headPosition = mParticles.GetPosition(npc.DipoleState->SecondaryParticleState.ParticleIndex);
			vec2f const feetPosition = mParticles.GetPosition(npc.PrimaryParticleState.ParticleIndex);
			vec2f const actualBodyVector = feetPosition - headPosition; // From head to feet
			float const actualBodyLength = actualBodyVector.length();
			vec2f const actualBodyVDir = actualBodyVector.normalise(actualBodyLength);
			vec2f const actualBodyHDir = actualBodyVDir.to_perpendicular(); // Points R (of the screen)
			vec2f const neckPosition = headPosition + actualBodyVector * GameParameters::HumanNpcGeometry::HeadLengthFraction;
			vec2f const shoulderPosition = neckPosition + actualBodyVector * GameParameters::HumanNpcGeometry::ArmDepthFraction / 2.0f;
			vec2f const crotchPosition = headPosition + actualBodyVector * (GameParameters::HumanNpcGeometry::HeadLengthFraction + GameParameters::HumanNpcGeometry::TorsoLengthFraction);

			float const cosLeftArmAngle = std::cos(humanNpcState.LeftArmAngle);
			float const sinLeftArmAngle = std::sin(humanNpcState.LeftArmAngle);
			float const cosRightArmAngle = std::cos(humanNpcState.RightArmAngle);
			float const sinRightArmAngle = std::sin(humanNpcState.RightArmAngle);
			float const cosLeftLegAngle = std::cos(humanNpcState.LeftLegAngle);
			float const sinLeftLegAngle = std::sin(humanNpcState.LeftLegAngle);
			float const cosRightLegAngle = std::cos(humanNpcState.RightLegAngle);
			float const sinRightLegAngle = std::sin(humanNpcState.RightLegAngle);

			// Arm and Leg lengths are relative to ideal
			float const adjustedIdealHumanHeight = GameParameters::HumanNpcGeometry::BodyLength * mCurrentHumanNpcBodyLengthAdjustment;
			float const leftArmLength = adjustedIdealHumanHeight * GameParameters::HumanNpcGeometry::ArmLengthFraction * humanNpcState.LeftArmLengthMultiplier;
			float const rightArmLength = adjustedIdealHumanHeight * GameParameters::HumanNpcGeometry::ArmLengthFraction * humanNpcState.RightArmLengthMultiplier;
			float const leftLegLength = adjustedIdealHumanHeight * GameParameters::HumanNpcGeometry::LegLengthFraction * humanNpcState.LeftLegLengthMultiplier;
			float const rightLegLength = adjustedIdealHumanHeight * GameParameters::HumanNpcGeometry::LegLengthFraction * humanNpcState.RightLegLengthMultiplier;

			if (humanNpcState.CurrentFaceOrientation != 0.0f)
			{
				//
				// Front-back
				//

				float const halfHeadW = (adjustedIdealHumanHeight * GameParameters::HumanNpcGeometry::HeadWidthFraction) / 2.0f;
				float const halfTorsoW = (adjustedIdealHumanHeight * GameParameters::HumanNpcGeometry::TorsoWidthFraction) / 2.0f;
				float const halfArmW = (adjustedIdealHumanHeight * GameParameters::HumanNpcGeometry::ArmWidthFraction) / 2.0f;
				float const halfLegW = (adjustedIdealHumanHeight * GameParameters::HumanNpcGeometry::LegWidthFraction) / 2.0f;

				// Head

				shipRenderContext.UploadNpcQuad(
					Quadf(
						headPosition - actualBodyHDir * halfHeadW,
						headPosition + actualBodyHDir * halfHeadW,
						neckPosition - actualBodyHDir * halfHeadW,
						neckPosition + actualBodyHDir * halfHeadW),
					humanNpcState.CurrentFaceOrientation,
					humanNpcState.CurrentFaceDirectionX);

				// Arms and legs

				vec2f const leftArmJointPosition = shoulderPosition - actualBodyHDir * (halfTorsoW - halfArmW);
				vec2f const rightArmJointPosition = shoulderPosition + actualBodyHDir * (halfTorsoW - halfArmW);

				vec2f const leftLegJointPosition = crotchPosition - actualBodyHDir * (halfTorsoW - halfLegW);
				vec2f const rightLegJointPosition = crotchPosition + actualBodyHDir * (halfTorsoW - halfLegW);

				if (humanNpcState.CurrentFaceOrientation > 0.0f)
				{
					// Front

					// Left arm (on left side of the screen)
					vec2f const leftArmVector = actualBodyVDir.rotate(cosLeftArmAngle, sinLeftArmAngle) * leftArmLength;
					vec2f const leftArmTraverseDir = leftArmVector.normalise().to_perpendicular();
					shipRenderContext.UploadNpcQuad(
						Quadf(
							leftArmJointPosition - leftArmTraverseDir * halfArmW,
							leftArmJointPosition + leftArmTraverseDir * halfArmW,
							leftArmJointPosition + leftArmVector - leftArmTraverseDir * halfArmW,
							leftArmJointPosition + leftArmVector + leftArmTraverseDir * halfArmW),
						humanNpcState.CurrentFaceOrientation,
						humanNpcState.CurrentFaceDirectionX);

					// Right arm (on right side of the screen)
					vec2f const rightArmVector = actualBodyVDir.rotate(cosRightArmAngle, sinRightArmAngle) * rightArmLength;
					vec2f const rightArmTraverseDir = rightArmVector.normalise().to_perpendicular();
					shipRenderContext.UploadNpcQuad(
						Quadf(
							rightArmJointPosition - rightArmTraverseDir * halfArmW,
							rightArmJointPosition + rightArmTraverseDir * halfArmW,
							rightArmJointPosition + rightArmVector - rightArmTraverseDir * halfArmW,
							rightArmJointPosition + rightArmVector + rightArmTraverseDir * halfArmW),
						humanNpcState.CurrentFaceOrientation,
						humanNpcState.CurrentFaceDirectionX);

					// Left leg (on left side of the screen)
					vec2f const leftLegVector = actualBodyVDir.rotate(cosLeftLegAngle, sinLeftLegAngle) * leftLegLength;
					vec2f const leftLegTraverseDir = leftLegVector.normalise().to_perpendicular();
					shipRenderContext.UploadNpcQuad(
						Quadf(
							leftLegJointPosition - leftLegTraverseDir * halfLegW,
							leftLegJointPosition + leftLegTraverseDir * halfLegW,
							leftLegJointPosition + leftLegVector - leftLegTraverseDir * halfLegW,
							leftLegJointPosition + leftLegVector + leftLegTraverseDir * halfLegW),
						humanNpcState.CurrentFaceOrientation,
						humanNpcState.CurrentFaceDirectionX);

					// Right leg (on right side of the screen)
					vec2f const rightLegVector = actualBodyVDir.rotate(cosRightLegAngle, sinRightLegAngle) * rightLegLength;
					vec2f const rightLegTraverseDir = rightLegVector.normalise().to_perpendicular();
					shipRenderContext.UploadNpcQuad(
						Quadf(
							rightLegJointPosition - rightLegTraverseDir * halfLegW,
							rightLegJointPosition + rightLegTraverseDir * halfLegW,
							rightLegJointPosition + rightLegVector - rightLegTraverseDir * halfLegW,
							rightLegJointPosition + rightLegVector + rightLegTraverseDir * halfLegW),
						humanNpcState.CurrentFaceOrientation,
						humanNpcState.CurrentFaceDirectionX);
				}
				else
				{
					// Back

					// Left arm (on right side of screen)
					vec2f const leftArmVector = actualBodyVDir.rotate(cosLeftArmAngle, -sinLeftArmAngle) * leftArmLength;
					vec2f const leftArmTraverseDir = leftArmVector.normalise().to_perpendicular();
					shipRenderContext.UploadNpcQuad(
						Quadf(
							rightArmJointPosition - leftArmTraverseDir * halfArmW,
							rightArmJointPosition + leftArmTraverseDir * halfArmW,
							rightArmJointPosition + leftArmVector - leftArmTraverseDir * halfArmW,
							rightArmJointPosition + leftArmVector + leftArmTraverseDir * halfArmW),
						humanNpcState.CurrentFaceOrientation,
						humanNpcState.CurrentFaceDirectionX);

					// Right arm (on left side of the screen)
					vec2f const rightArmVector = actualBodyVDir.rotate(cosRightArmAngle, -sinRightArmAngle) * rightArmLength;
					vec2f const rightArmTraverseDir = rightArmVector.normalise().to_perpendicular();
					shipRenderContext.UploadNpcQuad(
						Quadf(
							leftArmJointPosition - rightArmTraverseDir * halfArmW,
							leftArmJointPosition + rightArmTraverseDir * halfArmW,
							leftArmJointPosition + rightArmVector - rightArmTraverseDir * halfArmW,
							leftArmJointPosition + rightArmVector + rightArmTraverseDir * halfArmW),
						humanNpcState.CurrentFaceOrientation,
						humanNpcState.CurrentFaceDirectionX);

					// Left leg (on right side of the screen)
					vec2f const leftLegVector = actualBodyVDir.rotate(cosLeftLegAngle, -sinLeftLegAngle) * leftLegLength;
					vec2f const leftLegTraverseDir = leftLegVector.normalise().to_perpendicular();
					shipRenderContext.UploadNpcQuad(
						Quadf(
							rightLegJointPosition - leftLegTraverseDir * halfLegW,
							rightLegJointPosition + leftLegTraverseDir * halfLegW,
							rightLegJointPosition + leftLegVector - leftLegTraverseDir * halfLegW,
							rightLegJointPosition + leftLegVector + leftLegTraverseDir * halfLegW),
						humanNpcState.CurrentFaceOrientation,
						humanNpcState.CurrentFaceDirectionX);

					// Right leg (on left side of the screen)
					vec2f const rightLegVector = actualBodyVDir.rotate(cosRightLegAngle, -sinRightLegAngle) * rightLegLength;
					vec2f const rightLegTraverseDir = rightLegVector.normalise().to_perpendicular();
					shipRenderContext.UploadNpcQuad(
						Quadf(
							leftLegJointPosition - rightLegTraverseDir * halfLegW,
							leftLegJointPosition + rightLegTraverseDir * halfLegW,
							leftLegJointPosition + rightLegVector - rightLegTraverseDir * halfLegW,
							leftLegJointPosition + rightLegVector + rightLegTraverseDir * halfLegW),
						humanNpcState.CurrentFaceOrientation,
						humanNpcState.CurrentFaceDirectionX);
				}

				// Torso

				shipRenderContext.UploadNpcQuad(
					Quadf(
						neckPosition - actualBodyHDir * halfTorsoW,
						neckPosition + actualBodyHDir * halfTorsoW,
						crotchPosition - actualBodyHDir * halfTorsoW,
						crotchPosition + actualBodyHDir * halfTorsoW),
					humanNpcState.CurrentFaceOrientation,
					humanNpcState.CurrentFaceDirectionX);
			}
			else
			{
				//
				// Left-Right
				//

				float const halfHeadD = (adjustedIdealHumanHeight * GameParameters::HumanNpcGeometry::HeadDepthFraction) / 2.0f;
				float const halfTorsoD = (adjustedIdealHumanHeight * GameParameters::HumanNpcGeometry::TorsoDepthFraction) / 2.0f;
				float const halfArmD = (adjustedIdealHumanHeight * GameParameters::HumanNpcGeometry::ArmDepthFraction) / 2.0f;
				float const halfLegD = (adjustedIdealHumanHeight * GameParameters::HumanNpcGeometry::LegDepthFraction) / 2.0f;

				// Note: angles are with vertical, regardless of L/R

				vec2f const leftArmVector = actualBodyVDir.rotate(cosLeftArmAngle, sinLeftArmAngle) * leftArmLength;
				vec2f const leftArmTraverseDir = leftArmVector.normalise().to_perpendicular();
				Quadf leftArmQuad(
					shoulderPosition - leftArmTraverseDir * halfArmD,
					shoulderPosition + leftArmTraverseDir * halfArmD,
					shoulderPosition + leftArmVector - leftArmTraverseDir * halfArmD,
					shoulderPosition + leftArmVector + leftArmTraverseDir * halfArmD);

				vec2f const rightArmVector = actualBodyVDir.rotate(cosRightArmAngle, sinRightArmAngle) * rightArmLength;
				vec2f const rightArmTraverseDir = rightArmVector.normalise().to_perpendicular();
				Quadf rightArmQuad(
					shoulderPosition - rightArmTraverseDir * halfArmD,
					shoulderPosition + rightArmTraverseDir * halfArmD,
					shoulderPosition + rightArmVector - rightArmTraverseDir * halfArmD,
					shoulderPosition + rightArmVector + rightArmTraverseDir * halfArmD);

				vec2f const leftLegVector = actualBodyVDir.rotate(cosLeftLegAngle, sinLeftLegAngle) * leftLegLength;
				vec2f const leftLegTraverseDir = leftLegVector.normalise().to_perpendicular();
				Quadf leftLegQuad(
					crotchPosition - leftLegTraverseDir * halfLegD,
					crotchPosition + leftLegTraverseDir * halfLegD,
					crotchPosition + leftLegVector - leftLegTraverseDir * halfLegD,
					crotchPosition + leftLegVector + leftLegTraverseDir * halfLegD);

				vec2f const rightLegVector = actualBodyVDir.rotate(cosRightLegAngle, sinRightLegAngle) * rightLegLength;
				vec2f const rightLegTraverseDir = rightLegVector.normalise().to_perpendicular();
				Quadf rightLegQuad(
					crotchPosition - rightLegTraverseDir * halfLegD,
					crotchPosition + rightLegTraverseDir * halfLegD,
					crotchPosition + rightLegVector - rightLegTraverseDir * halfLegD,
					crotchPosition + rightLegVector + rightLegTraverseDir * halfLegD);

				// Arm and legs far

				if (humanNpcState.CurrentFaceDirectionX > 0.0f)
				{
					// Left arm
					shipRenderContext.UploadNpcQuad(
						leftArmQuad,
						humanNpcState.CurrentFaceOrientation,
						humanNpcState.CurrentFaceDirectionX);

					// Left leg
					shipRenderContext.UploadNpcQuad(
						leftLegQuad,
						humanNpcState.CurrentFaceOrientation,
						humanNpcState.CurrentFaceDirectionX);
				}
				else
				{
					// Right arm
					shipRenderContext.UploadNpcQuad(
						rightArmQuad,
						humanNpcState.CurrentFaceOrientation,
						humanNpcState.CurrentFaceDirectionX);

					// Right leg
					shipRenderContext.UploadNpcQuad(
						rightLegQuad,
						humanNpcState.CurrentFaceOrientation,
						humanNpcState.CurrentFaceDirectionX);
				}

				// Head

				shipRenderContext.UploadNpcQuad(
					Quadf(
						headPosition - actualBodyHDir * halfHeadD,
						headPosition + actualBodyHDir * halfHeadD,
						neckPosition - actualBodyHDir * halfHeadD,
						neckPosition + actualBodyHDir * halfHeadD),
					humanNpcState.CurrentFaceOrientation,
					humanNpcState.CurrentFaceDirectionX);

				// Torso

				shipRenderContext.UploadNpcQuad(
					Quadf(
						neckPosition - actualBodyHDir * halfTorsoD,
						neckPosition + actualBodyHDir * halfTorsoD,
						crotchPosition - actualBodyHDir * halfTorsoD,
						crotchPosition + actualBodyHDir * halfTorsoD),
					humanNpcState.CurrentFaceOrientation,
					humanNpcState.CurrentFaceDirectionX);

				// Arms and legs near

				if (humanNpcState.CurrentFaceDirectionX > 0.0f)
				{
					// Right arm
					shipRenderContext.UploadNpcQuad(
						rightArmQuad,
						humanNpcState.CurrentFaceOrientation,
						humanNpcState.CurrentFaceDirectionX);

					// Right leg
					shipRenderContext.UploadNpcQuad(
						rightLegQuad,
						humanNpcState.CurrentFaceOrientation,
						humanNpcState.CurrentFaceDirectionX);
				}
				else
				{
					// Left arm
					shipRenderContext.UploadNpcQuad(
						leftArmQuad,
						humanNpcState.CurrentFaceOrientation,
						humanNpcState.CurrentFaceDirectionX);

					// Left leg
					shipRenderContext.UploadNpcQuad(
						leftLegQuad,
						humanNpcState.CurrentFaceOrientation,
						humanNpcState.CurrentFaceDirectionX);
				}
			}

			break;
		}

		case NpcKindType::Furniture:
		{
			shipRenderContext.UploadNpcParticle(
				mParticles.GetPosition(npc.PrimaryParticleState.ParticleIndex),
				mParticles.GetRenderColor(npc.PrimaryParticleState.ParticleIndex),
				1.0f);

			if (npc.DipoleState.has_value())
			{
				shipRenderContext.UploadNpcParticle(
					mParticles.GetPosition(npc.DipoleState->SecondaryParticleState.ParticleIndex),
					mParticles.GetRenderColor(npc.DipoleState->SecondaryParticleState.ParticleIndex),
					1.0f);
			}

			break;
		}
	}
}

void Npcs::Publish() const
{
#ifdef IN_BARYLAB
	std::optional<ConstrainedRegimeParticleProbe> constrainedRegimeParticleProbe;
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

				if (state.PrimaryParticleState.ParticleIndex == *mCurrentlySelectedParticle)
				{
					if (state.PrimaryParticleState.ConstrainedState.has_value())
					{
						constrainedRegimeParticleProbe.emplace(
							state.PrimaryParticleState.ConstrainedState->CurrentTriangle,
							state.PrimaryParticleState.ConstrainedState->CurrentTriangleBarycentricCoords);

						if (mCurrentOriginTriangle.has_value())
						{
							subjectParticleBarycentricCoordinatesWrtOriginTriangleChanged = shipMesh.ShipTriangles.ToBarycentricCoordinates(
								mParticles.GetPosition(state.PrimaryParticleState.ParticleIndex),
								*mCurrentOriginTriangle,
								shipMesh.ShipPoints);
						}
					}

					physicsParticleProbe.emplace(mParticles.GetVelocity(state.PrimaryParticleState.ParticleIndex));
				}
				else if (state.DipoleState.has_value()
					&& state.DipoleState->SecondaryParticleState.ParticleIndex == *mCurrentlySelectedParticle)
				{
					if (state.DipoleState->SecondaryParticleState.ConstrainedState.has_value())
					{
						constrainedRegimeParticleProbe.emplace(
							state.DipoleState->SecondaryParticleState.ConstrainedState->CurrentTriangle,
							state.DipoleState->SecondaryParticleState.ConstrainedState->CurrentTriangleBarycentricCoords);

						if (mCurrentOriginTriangle.has_value())
						{
							subjectParticleBarycentricCoordinatesWrtOriginTriangleChanged = shipMesh.ShipTriangles.ToBarycentricCoordinates(
								mParticles.GetPosition(state.DipoleState->SecondaryParticleState.ParticleIndex),
								*mCurrentOriginTriangle,
								shipMesh.ShipPoints);
						}
					}

					physicsParticleProbe.emplace(mParticles.GetVelocity(state.DipoleState->SecondaryParticleState.ParticleIndex));
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
				case StateType::KindSpecificStateType::HumanNpcStateType::BehaviorType::Constrained_KnockedOut:
				{
					mGameEventHandler->OnHumanNpcBehaviorChanged("Constrained_KnockedOut");

					break;
				}

				case StateType::KindSpecificStateType::HumanNpcStateType::BehaviorType::Free_Aerial:
				{
					mGameEventHandler->OnHumanNpcBehaviorChanged("Free_Aerial");

					break;
				}

				default:
				{
					mGameEventHandler->OnHumanNpcBehaviorChanged("?");

					break;
				}
			}
		}
	}
#endif
}

}