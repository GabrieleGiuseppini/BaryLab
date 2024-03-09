/***************************************************************************************
* Original Author:      Gabriele Giuseppini
* Created:              2023-10-06
* Copyright:            Gabriele Giuseppini  (https://github.com/GabrieleGiuseppini)
***************************************************************************************/
#include "Physics.h"

#include "Colors.h"

#include <cassert>

namespace Physics {

void Npcs::Add(
	NpcType npcType,
	vec2f primaryPosition,
	std::optional<vec2f> secondaryPosition,
	float currentSimulationTime,
	Ship const & ship,
	GameParameters const & gameParameters)
{
	switch (npcType)
	{
		case NpcType::Furniture:
		{
			assert(mParticles.GetRemainingParticlesCount() >= 1);

			auto const & material = mMaterialDatabase.GetNpcMaterial(NpcMaterial::KindType::Furniture);

			auto const particleIndex = mParticles.Add(
				material.Mass,
				material.StaticFriction,
				material.KineticFriction,
				material.Elasticity,
				material.BuoyancyVolumeFill,
				primaryPosition,
				material.RenderColor);

			StateType::NpcParticleStateType primaryParticleState = StateType::NpcParticleStateType(
				particleIndex,
				CalculateParticleConstrainedState(
					primaryPosition,
					ship));

			mStateBuffer.emplace_back(
				npcType,
				primaryParticleState.ConstrainedState.has_value() ? StateType::RegimeType::Constrained : StateType::RegimeType::Free,
				std::move(primaryParticleState),
				std::nullopt,
				std::nullopt);

			return;
		}

		case NpcType::Human:
		{
			assert(mParticles.GetRemainingParticlesCount() >= 2);

			float const bodyLength = GameParameters::HumanNpcGeometry::BodyLength * gameParameters.HumanNpcBodyLengthAdjustment;

			// Feet (primary)

			auto const & feetMaterial = mMaterialDatabase.GetNpcMaterial(NpcMaterial::KindType::HumanFeet);

			auto const primaryParticleIndex = mParticles.Add(
				feetMaterial.Mass,
				feetMaterial.StaticFriction,
				feetMaterial.KineticFriction,
				feetMaterial.Elasticity,
				feetMaterial.BuoyancyVolumeFill,
				primaryPosition,
				feetMaterial.RenderColor);

			StateType::NpcParticleStateType primaryParticleState = StateType::NpcParticleStateType(
				primaryParticleIndex,
				CalculateParticleConstrainedState(
					primaryPosition,
					ship));

			// Head (secondary)

			auto const & headMaterial = mMaterialDatabase.GetNpcMaterial(NpcMaterial::KindType::HumanHead);

			if (!secondaryPosition)
			{
				secondaryPosition = primaryPosition + vec2f(0.0f, 1.0f) * bodyLength;
			}

			auto const secondaryParticleIndex = mParticles.Add(
				headMaterial.Mass,
				headMaterial.StaticFriction,
				headMaterial.KineticFriction,
				headMaterial.Elasticity,
				headMaterial.BuoyancyVolumeFill,
				*secondaryPosition,
				headMaterial.RenderColor);

			StateType::NpcParticleStateType secondaryParticleState = StateType::NpcParticleStateType(
				secondaryParticleIndex,
				CalculateParticleConstrainedState(
					*secondaryPosition,
					ship));

			float const massFactor =
				(feetMaterial.Mass * headMaterial.Mass)
				/ (feetMaterial.Mass + headMaterial.Mass);

			StateType::DipoleStateType dipoleState = StateType::DipoleStateType(
				std::move(secondaryParticleState),
				StateType::DipolePropertiesType(
					bodyLength,
					massFactor,
					1.0f));

			StateType::HumanNpcStateType humanState = InitializeHuman(
				primaryParticleState,
				secondaryParticleState,
				currentSimulationTime,
				ship);

			auto const regime = primaryParticleState.ConstrainedState.has_value()
				? StateType::RegimeType::Constrained
				: (dipoleState.SecondaryParticleState.ConstrainedState.has_value() ? StateType::RegimeType::Constrained : StateType::RegimeType::Free);

			mStateBuffer.emplace_back(
				npcType,
				regime,
				std::move(primaryParticleState),
				std::move(dipoleState),
				std::move(humanState));

			return;
		}
	}
}

void Npcs::SetPanicLevelForAllHumans(float panicLevel)
{
	for (auto & npc : mStateBuffer)
	{
		if (npc.Type == NpcType::Human)
		{
			assert(npc.HumanNpcState.has_value());
			npc.HumanNpcState->PanicLevel = panicLevel;
		}
	}
}

void Npcs::FlipHumanWalk(int npcIndex)
{
	if (npcIndex < mStateBuffer.size()
		&& mStateBuffer[npcIndex].HumanNpcState.has_value()
		&& mStateBuffer[npcIndex].HumanNpcState->CurrentBehavior == StateType::HumanNpcStateType::BehaviorType::Constrained_Walking)
	{
		FlipHumanWalk(*mStateBuffer[npcIndex].HumanNpcState, StrongTypedTrue<_DoImmediate>);
	}
}

void Npcs::FlipHumanFrontBack(int npcIndex)
{
	if (npcIndex < mStateBuffer.size()
		&& mStateBuffer[npcIndex].HumanNpcState.has_value())
	{
		mStateBuffer[npcIndex].HumanNpcState->CurrentFaceOrientation *= -1.0f;
	}
}

void Npcs::MoveParticleBy(
	ElementIndex particleIndex,
	vec2f const & offset,
	float currentSimulationTime,
	Ship const & ship)
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
	// Initialize NPC state
	//

	for (auto const n : *this)
	{
		auto & state = mStateBuffer[n];

		if (state.PrimaryParticleState.ParticleIndex == particleIndex
			|| (state.DipoleState.has_value() && state.DipoleState->SecondaryParticleState.ParticleIndex == particleIndex))
		{
			state = MaterializeNpcState(n, currentSimulationTime, ship);
			break;
		}
	}

	//
	// Select particle
	//

	SelectParticle(particleIndex, ship);

	//
	// Reset trajectories
	//

	mCurrentParticleTrajectory.reset();
	mCurrentParticleTrajectoryNotification.reset();
}

void Npcs::RotateParticlesWithShip(
	vec2f const & centerPos,
	float cosAngle,
	float sinAngle,
	Ship const & ship)
{
	//
	// Rotate particles
	//

	for (auto const n : *this)
	{
		auto & state = mStateBuffer[n];

		RotateParticleWithShip(
			state.PrimaryParticleState,
			centerPos,
			cosAngle,
			sinAngle,
			ship);

		if (state.DipoleState.has_value())
		{
			RotateParticleWithShip(
				state.DipoleState->SecondaryParticleState,
				centerPos,
				cosAngle,
				sinAngle,
				ship);
		}
	}
}

void Npcs::OnPointMoved(
	float currentSimulationTime,
	Ship const & ship)
{
	//
	// Recalculate state of all NPCs
	//

	for (auto const n : *this)
	{
		mStateBuffer[n] = MaterializeNpcState(n, currentSimulationTime, ship);
	}
}

void Npcs::Update(
	float currentSimulationTime,
	Ship const & ship,
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

	UpdateNpcs(currentSimulationTime, ship, gameParameters);

	//
	// Publish
	//

	Publish(ship);
}

void Npcs::Render(RenderContext & renderContext)
{
	//
	// Particles & limbs
	//

	float const adjustedStandardHumanHeight = GameParameters::HumanNpcGeometry::BodyLength * mCurrentHumanNpcBodyLengthAdjustment;

	renderContext.UploadParticlesStart();
	renderContext.UploadNpcHumanLimbsStart();

	for (auto const i : *this)
	{
		auto const & state = mStateBuffer[i];

		switch (mNpcRenderMode)
		{
			case NpcRenderMode::Limbs:
			{
				if (state.HumanNpcState.has_value())
				{
					assert(state.DipoleState.has_value());

					// Note:
					// - head, neck, shoulder, crotch, feet: based on current dipole length
					// - arms, legs : based on ideal * adjustment
					//

					vec2f const headPosition = mParticles.GetPosition(state.DipoleState->SecondaryParticleState.ParticleIndex);
					vec2f const feetPosition = mParticles.GetPosition(state.PrimaryParticleState.ParticleIndex);
					vec2f const actualBodyVector = feetPosition - headPosition; // From head to feet
					float const actualBodyLength = actualBodyVector.length();
					vec2f const actualBodyVDir = actualBodyVector.normalise(actualBodyLength);
					vec2f const actualBodyHDir = actualBodyVDir.to_perpendicular(); // Points R (of the screen)
					vec2f const neckPosition = headPosition + actualBodyVector * GameParameters::HumanNpcGeometry::HeadLengthFraction;
					vec2f const shoulderPosition = neckPosition + actualBodyVector * GameParameters::HumanNpcGeometry::ArmDepthFraction / 2.0f;
					vec2f const crotchPosition = headPosition + actualBodyVector * (GameParameters::HumanNpcGeometry::HeadLengthFraction + GameParameters::HumanNpcGeometry::TorsoLengthFraction);

					float const cosLeftArmAngle = std::cos(state.HumanNpcState->LeftArmAngle);
					float const sinLeftArmAngle = std::sin(state.HumanNpcState->LeftArmAngle);
					float const cosRightArmAngle = std::cos(state.HumanNpcState->RightArmAngle);
					float const sinRightArmAngle = std::sin(state.HumanNpcState->RightArmAngle);
					float const cosLeftLegAngle = std::cos(state.HumanNpcState->LeftLegAngle);
					float const sinLeftLegAngle = std::sin(state.HumanNpcState->LeftLegAngle);
					float const cosRightLegAngle = std::cos(state.HumanNpcState->RightLegAngle);
					float const sinRightLegAngle = std::sin(state.HumanNpcState->RightLegAngle);

					float const leftArmLength = adjustedStandardHumanHeight * GameParameters::HumanNpcGeometry::ArmLengthFraction * state.HumanNpcState->LeftArmLengthMultiplier;
					float const rightArmLength = adjustedStandardHumanHeight * GameParameters::HumanNpcGeometry::ArmLengthFraction * state.HumanNpcState->RightArmLengthMultiplier;
					float const leftLegLength = adjustedStandardHumanHeight * GameParameters::HumanNpcGeometry::LegLengthFraction * state.HumanNpcState->LeftLegLengthMultiplier;
					float const rightLegLength = adjustedStandardHumanHeight * GameParameters::HumanNpcGeometry::LegLengthFraction * state.HumanNpcState->RightLegLengthMultiplier;

					if (state.HumanNpcState->CurrentFaceOrientation != 0.0f)
					{
						//
						// Front-back
						//

						float const halfHeadW = (adjustedStandardHumanHeight * GameParameters::HumanNpcGeometry::HeadWidthFraction) / 2.0f;
						float const halfTorsoW = (adjustedStandardHumanHeight * GameParameters::HumanNpcGeometry::TorsoWidthFraction) / 2.0f;
						float const halfArmW = (adjustedStandardHumanHeight * GameParameters::HumanNpcGeometry::ArmWidthFraction) / 2.0f;
						float const halfLegW = (adjustedStandardHumanHeight * GameParameters::HumanNpcGeometry::LegWidthFraction) / 2.0f;

						// Head

						renderContext.UploadNpcHumanLimb(
							Quadf(
								headPosition - actualBodyHDir * halfHeadW,
								headPosition + actualBodyHDir * halfHeadW,
								neckPosition - actualBodyHDir * halfHeadW,
								neckPosition + actualBodyHDir * halfHeadW),
							state.HumanNpcState->CurrentFaceOrientation,
							state.HumanNpcState->CurrentFaceDirectionX);

						// Arms and legs

						vec2f const leftArmJointPosition = shoulderPosition - actualBodyHDir * (halfTorsoW - halfArmW);
						vec2f const rightArmJointPosition = shoulderPosition + actualBodyHDir * (halfTorsoW - halfArmW);

						vec2f const leftLegJointPosition = crotchPosition - actualBodyHDir * (halfTorsoW - halfLegW);
						vec2f const rightLegJointPosition = crotchPosition + actualBodyHDir * (halfTorsoW - halfLegW);

						if (state.HumanNpcState->CurrentFaceOrientation > 0.0f)
						{
							// Front

							// Left arm (on left side of the screen)
							vec2f const leftArmVector = actualBodyVDir.rotate(cosLeftArmAngle, sinLeftArmAngle) * leftArmLength;
							vec2f const leftArmTraverseDir = leftArmVector.normalise().to_perpendicular();
							renderContext.UploadNpcHumanLimb(
								Quadf(
									leftArmJointPosition - leftArmTraverseDir * halfArmW,
									leftArmJointPosition + leftArmTraverseDir * halfArmW,
									leftArmJointPosition + leftArmVector - leftArmTraverseDir * halfArmW,
									leftArmJointPosition + leftArmVector + leftArmTraverseDir * halfArmW),
								state.HumanNpcState->CurrentFaceOrientation,
								state.HumanNpcState->CurrentFaceDirectionX);

							// Right arm (on right side of the screen)
							vec2f const rightArmVector = actualBodyVDir.rotate(cosRightArmAngle, sinRightArmAngle) * rightArmLength;
							vec2f const rightArmTraverseDir = rightArmVector.normalise().to_perpendicular();
							renderContext.UploadNpcHumanLimb(
								Quadf(
									rightArmJointPosition - rightArmTraverseDir * halfArmW,
									rightArmJointPosition + rightArmTraverseDir * halfArmW,
									rightArmJointPosition + rightArmVector - rightArmTraverseDir * halfArmW,
									rightArmJointPosition + rightArmVector + rightArmTraverseDir * halfArmW),
								state.HumanNpcState->CurrentFaceOrientation,
								state.HumanNpcState->CurrentFaceDirectionX);

							// Left leg (on left side of the screen)
							vec2f const leftLegVector = actualBodyVDir.rotate(cosLeftLegAngle, sinLeftLegAngle) * leftLegLength;
							vec2f const leftLegTraverseDir = leftLegVector.normalise().to_perpendicular();
							renderContext.UploadNpcHumanLimb(
								Quadf(
									leftLegJointPosition - leftLegTraverseDir * halfLegW,
									leftLegJointPosition + leftLegTraverseDir * halfLegW,
									leftLegJointPosition + leftLegVector - leftLegTraverseDir * halfLegW,
									leftLegJointPosition + leftLegVector + leftLegTraverseDir * halfLegW),
								state.HumanNpcState->CurrentFaceOrientation,
								state.HumanNpcState->CurrentFaceDirectionX);

							// Right leg (on right side of the screen)
							vec2f const rightLegVector = actualBodyVDir.rotate(cosRightLegAngle, sinRightLegAngle) * rightLegLength;
							vec2f const rightLegTraverseDir = rightLegVector.normalise().to_perpendicular();
							renderContext.UploadNpcHumanLimb(
								Quadf(
									rightLegJointPosition - rightLegTraverseDir * halfLegW,
									rightLegJointPosition + rightLegTraverseDir * halfLegW,
									rightLegJointPosition + rightLegVector - rightLegTraverseDir * halfLegW,
									rightLegJointPosition + rightLegVector + rightLegTraverseDir * halfLegW),
								state.HumanNpcState->CurrentFaceOrientation,
								state.HumanNpcState->CurrentFaceDirectionX);
						}
						else
						{
							// Back

							// Left arm (on right side of screen)
							vec2f const leftArmVector = actualBodyVDir.rotate(cosLeftArmAngle, -sinLeftArmAngle) * leftArmLength;
							vec2f const leftArmTraverseDir = leftArmVector.normalise().to_perpendicular();
							renderContext.UploadNpcHumanLimb(
								Quadf(
									rightArmJointPosition - leftArmTraverseDir * halfArmW,
									rightArmJointPosition + leftArmTraverseDir * halfArmW,
									rightArmJointPosition + leftArmVector - leftArmTraverseDir * halfArmW,
									rightArmJointPosition + leftArmVector + leftArmTraverseDir * halfArmW),
								state.HumanNpcState->CurrentFaceOrientation,
								state.HumanNpcState->CurrentFaceDirectionX);

							// Right arm (on left side of the screen)
							vec2f const rightArmVector = actualBodyVDir.rotate(cosRightArmAngle, -sinRightArmAngle) * rightArmLength;
							vec2f const rightArmTraverseDir = rightArmVector.normalise().to_perpendicular();
							renderContext.UploadNpcHumanLimb(
								Quadf(
									leftArmJointPosition - rightArmTraverseDir * halfArmW,
									leftArmJointPosition + rightArmTraverseDir * halfArmW,
									leftArmJointPosition + rightArmVector - rightArmTraverseDir * halfArmW,
									leftArmJointPosition + rightArmVector + rightArmTraverseDir * halfArmW),
								state.HumanNpcState->CurrentFaceOrientation,
								state.HumanNpcState->CurrentFaceDirectionX);

							// Left leg (on right side of the screen)
							vec2f const leftLegVector = actualBodyVDir.rotate(cosLeftLegAngle, -sinLeftLegAngle) * leftLegLength;
							vec2f const leftLegTraverseDir = leftLegVector.normalise().to_perpendicular();
							renderContext.UploadNpcHumanLimb(
								Quadf(
									rightLegJointPosition - leftLegTraverseDir * halfLegW,
									rightLegJointPosition + leftLegTraverseDir * halfLegW,
									rightLegJointPosition + leftLegVector - leftLegTraverseDir * halfLegW,
									rightLegJointPosition + leftLegVector + leftLegTraverseDir * halfLegW),
								state.HumanNpcState->CurrentFaceOrientation,
								state.HumanNpcState->CurrentFaceDirectionX);

							// Right leg (on left side of the screen)
							vec2f const rightLegVector = actualBodyVDir.rotate(cosRightLegAngle, -sinRightLegAngle) * rightLegLength;
							vec2f const rightLegTraverseDir = rightLegVector.normalise().to_perpendicular();
							renderContext.UploadNpcHumanLimb(
								Quadf(
									leftLegJointPosition - rightLegTraverseDir * halfLegW,
									leftLegJointPosition + rightLegTraverseDir * halfLegW,
									leftLegJointPosition + rightLegVector - rightLegTraverseDir * halfLegW,
									leftLegJointPosition + rightLegVector + rightLegTraverseDir * halfLegW),
								state.HumanNpcState->CurrentFaceOrientation,
								state.HumanNpcState->CurrentFaceDirectionX);
						}

						// Torso

						renderContext.UploadNpcHumanLimb(
							Quadf(
								neckPosition - actualBodyHDir * halfTorsoW,
								neckPosition + actualBodyHDir * halfTorsoW,
								crotchPosition - actualBodyHDir * halfTorsoW,
								crotchPosition + actualBodyHDir * halfTorsoW),
							state.HumanNpcState->CurrentFaceOrientation,
							state.HumanNpcState->CurrentFaceDirectionX);
					}
					else
					{
						//
						// Left-Right
						//

						float const halfHeadD = (adjustedStandardHumanHeight * GameParameters::HumanNpcGeometry::HeadDepthFraction) / 2.0f;
						float const halfTorsoD = (adjustedStandardHumanHeight * GameParameters::HumanNpcGeometry::TorsoDepthFraction) / 2.0f;
						float const halfArmD = (adjustedStandardHumanHeight * GameParameters::HumanNpcGeometry::ArmDepthFraction) / 2.0f;
						float const halfLegD = (adjustedStandardHumanHeight * GameParameters::HumanNpcGeometry::LegDepthFraction) / 2.0f;

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

						if (state.HumanNpcState->CurrentFaceDirectionX > 0.0f)
						{
							// Left arm
							renderContext.UploadNpcHumanLimb(
								leftArmQuad,
								state.HumanNpcState->CurrentFaceOrientation,
								state.HumanNpcState->CurrentFaceDirectionX);

							// Left leg
							renderContext.UploadNpcHumanLimb(
								leftLegQuad,
								state.HumanNpcState->CurrentFaceOrientation,
								state.HumanNpcState->CurrentFaceDirectionX);
						}
						else
						{
							// Right arm
							renderContext.UploadNpcHumanLimb(
								rightArmQuad,
								state.HumanNpcState->CurrentFaceOrientation,
								state.HumanNpcState->CurrentFaceDirectionX);

							// Right leg
							renderContext.UploadNpcHumanLimb(
								rightLegQuad,
								state.HumanNpcState->CurrentFaceOrientation,
								state.HumanNpcState->CurrentFaceDirectionX);
						}

						// Head

						renderContext.UploadNpcHumanLimb(
							Quadf(
								headPosition - actualBodyHDir * halfHeadD,
								headPosition + actualBodyHDir * halfHeadD,
								neckPosition - actualBodyHDir * halfHeadD,
								neckPosition + actualBodyHDir * halfHeadD),
							state.HumanNpcState->CurrentFaceOrientation,
							state.HumanNpcState->CurrentFaceDirectionX);

						// Torso

						renderContext.UploadNpcHumanLimb(
							Quadf(
								neckPosition - actualBodyHDir * halfTorsoD,
								neckPosition + actualBodyHDir * halfTorsoD,
								crotchPosition - actualBodyHDir * halfTorsoD,
								crotchPosition + actualBodyHDir * halfTorsoD),
							state.HumanNpcState->CurrentFaceOrientation,
							state.HumanNpcState->CurrentFaceDirectionX);

						// Arms and legs near

						if (state.HumanNpcState->CurrentFaceDirectionX > 0.0f)
						{
							// Right arm
							renderContext.UploadNpcHumanLimb(
								rightArmQuad,
								state.HumanNpcState->CurrentFaceOrientation,
								state.HumanNpcState->CurrentFaceDirectionX);

							// Right leg
							renderContext.UploadNpcHumanLimb(
								rightLegQuad,
								state.HumanNpcState->CurrentFaceOrientation,
								state.HumanNpcState->CurrentFaceDirectionX);
						}
						else
						{
							// Left arm
							renderContext.UploadNpcHumanLimb(
								leftArmQuad,
								state.HumanNpcState->CurrentFaceOrientation,
								state.HumanNpcState->CurrentFaceDirectionX);

							// Left leg
							renderContext.UploadNpcHumanLimb(
								leftLegQuad,
								state.HumanNpcState->CurrentFaceOrientation,
								state.HumanNpcState->CurrentFaceDirectionX);
						}
					}
				}
				else
				{
					RenderParticle(state.PrimaryParticleState, renderContext);

					if (state.DipoleState.has_value())
					{
						RenderParticle(state.DipoleState->SecondaryParticleState, renderContext);
					}
				}

				break;
			}

			case NpcRenderMode::Physical:
			{
				RenderParticle(state.PrimaryParticleState, renderContext);

				if (state.DipoleState.has_value())
				{
					RenderParticle(state.DipoleState->SecondaryParticleState, renderContext);
				}
			}
		}
	}

	renderContext.UploadNpcHumanLimbsEnd();
	renderContext.UploadParticlesEnd();

	if (mNpcRenderMode == NpcRenderMode::Physical)
	{
		//
		// Springs
		//

		renderContext.UploadSpringsStart();

		for (auto const i : *this)
		{
			auto const & state = mStateBuffer[i];

			if (state.DipoleState.has_value())
			{
				renderContext.UploadSpring(
					mParticles.GetPosition(state.PrimaryParticleState.ParticleIndex),
					mParticles.GetPosition(state.DipoleState->SecondaryParticleState.ParticleIndex),
					rgbaColor(0x4a, 0x4a, 0x4a, 0xff));
			}
		}

		renderContext.UploadSpringsEnd();
	}

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
}

bool Npcs::IsTriangleConstrainingCurrentlySelectedParticle(ElementIndex triangleIndex) const
{
	if (mCurrentlySelectedParticle.has_value())
	{
		for (auto const i : *this)
		{
			auto const & state = mStateBuffer[i];

			if (state.PrimaryParticleState.ParticleIndex == *mCurrentlySelectedParticle
				&& state.PrimaryParticleState.ConstrainedState.has_value()
				&& triangleIndex == state.PrimaryParticleState.ConstrainedState->CurrentTriangle)
			{
				return true;
			}

			if (state.DipoleState.has_value()
				&& state.DipoleState->SecondaryParticleState.ParticleIndex == *mCurrentlySelectedParticle
				&& state.DipoleState->SecondaryParticleState.ConstrainedState.has_value()
				&& triangleIndex == state.DipoleState->SecondaryParticleState.ConstrainedState->CurrentTriangle)
			{
				return true;
			}
		}
	}

	return false;
}

bool Npcs::IsSpringHostingCurrentlySelectedParticle(
	ElementIndex springIndex,
	Ship const & ship) const
{
	if (mCurrentlySelectedParticle.has_value())
	{
		for (auto const i : *this)
		{
			auto const & state = mStateBuffer[i];

			if (state.PrimaryParticleState.ParticleIndex == *mCurrentlySelectedParticle ||
				(state.DipoleState.has_value() && state.DipoleState->SecondaryParticleState.ParticleIndex == *mCurrentlySelectedParticle))
			{
				if (state.PrimaryParticleState.ConstrainedState.has_value()
					&& state.PrimaryParticleState.ConstrainedState->CurrentVirtualEdgeOrdinal >= 0
					&& springIndex == ship.GetTriangles().GetSubSprings(state.PrimaryParticleState.ConstrainedState->CurrentTriangle).SpringIndices[state.PrimaryParticleState.ConstrainedState->CurrentVirtualEdgeOrdinal])
				{
					return true;
				}

				if (state.DipoleState.has_value()
					&& state.DipoleState->SecondaryParticleState.ConstrainedState.has_value()
					&& state.DipoleState->SecondaryParticleState.ConstrainedState->CurrentVirtualEdgeOrdinal >= 0
					&& springIndex == ship.GetTriangles().GetSubSprings(state.DipoleState->SecondaryParticleState.ConstrainedState->CurrentTriangle).SpringIndices[state.DipoleState->SecondaryParticleState.ConstrainedState->CurrentVirtualEdgeOrdinal])
				{
					return true;
				}
			}
		}
	}

	return false;
}

////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////

void Npcs::RotateParticleWithShip(
	StateType::NpcParticleStateType const & npcParticleState,
	vec2f const & centerPos,
	float cosAngle,
	float sinAngle,
	Ship const & ship)
{
	vec2f newPosition;

	if (npcParticleState.ConstrainedState.has_value())
	{
		// Simply set position from current bary coords

		newPosition = ship.GetTriangles().FromBarycentricCoordinates(
			npcParticleState.ConstrainedState->CurrentTriangleBarycentricCoords,
			npcParticleState.ConstrainedState->CurrentTriangle,
			ship.GetPoints());
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

void Npcs::RenderParticle(
	StateType::NpcParticleStateType const & particleState,
	RenderContext & renderContext)
{
	renderContext.UploadParticle(
		mParticles.GetPosition(particleState.ParticleIndex),
		mParticles.GetRenderColor(particleState.ParticleIndex),
		1.0f);
}

void Npcs::Publish(Ship const & ship)
{
	std::optional<ConstrainedRegimeParticleProbe> constrainedRegimeParticleProbe;
	std::optional<bcoords3f> subjectParticleBarycentricCoordinatesWrtOriginTriangleChanged;
	std::optional<PhysicsParticleProbe> physicsParticleProbe;

	if (mCurrentlySelectedParticle.has_value())
	{
		for (auto const n : *this)
		{
			auto const & state = mStateBuffer[n];

			if (state.PrimaryParticleState.ParticleIndex == *mCurrentlySelectedParticle)
			{
				if (state.PrimaryParticleState.ConstrainedState.has_value())
				{
					constrainedRegimeParticleProbe.emplace(
						state.PrimaryParticleState.ConstrainedState->CurrentTriangle,
						state.PrimaryParticleState.ConstrainedState->CurrentTriangleBarycentricCoords);

					if (mCurrentOriginTriangle.has_value())
					{
						subjectParticleBarycentricCoordinatesWrtOriginTriangleChanged = ship.GetTriangles().ToBarycentricCoordinates(
							mParticles.GetPosition(state.PrimaryParticleState.ParticleIndex),
							*mCurrentOriginTriangle,
							ship.GetPoints());
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
						subjectParticleBarycentricCoordinatesWrtOriginTriangleChanged = ship.GetTriangles().ToBarycentricCoordinates(
							mParticles.GetPosition(state.DipoleState->SecondaryParticleState.ParticleIndex),
							*mCurrentOriginTriangle,
							ship.GetPoints());
					}
				}

				physicsParticleProbe.emplace(mParticles.GetVelocity(state.DipoleState->SecondaryParticleState.ParticleIndex));
			}
		}
	}

	mGameEventDispatcher.OnSubjectParticleConstrainedRegimeUpdated(constrainedRegimeParticleProbe);
	mGameEventDispatcher.OnSubjectParticleBarycentricCoordinatesWrtOriginTriangleChanged(subjectParticleBarycentricCoordinatesWrtOriginTriangleChanged);
	mGameEventDispatcher.OnSubjectParticlePhysicsUpdated(physicsParticleProbe);
}

}