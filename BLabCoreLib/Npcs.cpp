/***************************************************************************************
* Original Author:      Gabriele Giuseppini
* Created:              2023-10-06
* Copyright:            Gabriele Giuseppini  (https://github.com/GabrieleGiuseppini)
***************************************************************************************/
#include "Npcs.h"

#include "Colors.h"

#include <cassert>

void Npcs::Add(
	NpcType npcType,
	vec2f const & primaryPosition,
	Mesh const & mesh)
{
	assert(mParticles.GetParticleCount() < LabParameters::MaxNpcs);

	//
	// Add primary particle
	//

	ElementIndex const primaryParticleIndex = mParticles.GetParticleCount();
	mParticles.Add(primaryPosition, rgbaColor(0x60, 0x60, 0x60, 0xff));
	auto primaryParticleState = MaterializeParticleState(
		primaryPosition,
		primaryParticleIndex,
		mesh);

	//
	// Add secondary particle
	//

	std::optional<StateType::NpcParticleStateType> secondaryParticleState;

	if (npcType != NpcType::Furniture)
	{
		// TODO
		assert(false);
	}
	
	//
	// Calculate regime
	//

	auto const regime = primaryParticleState.ConstrainedState.has_value()
		? StateType::RegimeType::Constrained
		: (secondaryParticleState.has_value() && secondaryParticleState->ConstrainedState.has_value()) ? StateType::RegimeType::Constrained : StateType::RegimeType::Free;

	mStateBuffer.emplace_back(
		regime,
		std::move(primaryParticleState),
		std::move(secondaryParticleState));
}

void Npcs::MoveParticleBy(
	ElementIndex particleIndex,
	vec2f const & offset,
	Mesh const & mesh)
{
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
			|| (state.SecondaryParticleState.has_value() && state.SecondaryParticleState->ParticleIndex == particleIndex))
		{
			state = MaterializeNpcState(n, mesh);
			break;
		}
	}

	//
	// Reset simulation
	//

	ResetSimulationStepState();

	//
	// Select particle
	//

	SelectParticle(particleIndex);

	//
	// Reset trajectories
	//

	mCurrentParticleTrajectory.reset();
	mCurrentParticleTrajectoryNotification.reset();
}

void Npcs::RotateParticlesWithMesh(
	vec2f const & centerPos,
	float cosAngle,
	float sinAngle,
	Mesh const & mesh)
{
	//
	// Rotate particles
	//

	for (auto const n : *this)
	{
		auto & state = mStateBuffer[n];

		RotateParticleWithMesh(
			state.PrimaryParticleState,
			centerPos,
			cosAngle,
			sinAngle,
			mesh);

		if (state.SecondaryParticleState.has_value())
		{
			RotateParticleWithMesh(
				*state.SecondaryParticleState,
				centerPos,
				cosAngle,
				sinAngle,
				mesh);
		}
	}

	//
	// Reset simulation
	//

	ResetSimulationStepState();
}

void Npcs::OnVertexMoved(Mesh const & mesh)
{
	//
	// Recalculate state of all NPCs
	//

	for (auto const n : *this)
	{
		mStateBuffer[n] = MaterializeNpcState(n, mesh);
	}

	//
	// Reset simulation
	//

	ResetSimulationStepState();
}

void Npcs::Update(
	Mesh const & mesh,
	LabParameters const & labParameters)
{
	//
	// Update NPCs' state
	//

	UpdateNpcs(mesh, labParameters);

	//
	// Publish
	//

	std::optional<ConstrainedRegimeParticleProbe> constrainedRegimeParticleProbe;
	std::optional<vec3f> subjectParticleBarycentricCoordinatesWrtOriginTriangleChanged;
	std::optional<PhysicsParticleProbe> physicsParticleProbe;

	if (mCurrentlySelectedParticle.has_value())
	{
		for (auto const n : *this)
		{
			auto const & state = mStateBuffer[n];

			if (state.PrimaryParticleState.ParticleIndex == *mCurrentlySelectedParticle
				&& state.PrimaryParticleState.ConstrainedState.has_value())
			{
				constrainedRegimeParticleProbe.emplace(
					state.PrimaryParticleState.ConstrainedState->CurrentTriangle,
					state.PrimaryParticleState.ConstrainedState->CurrentTriangleBarycentricCoords);
				
				if (mCurrentOriginTriangle.has_value())
				{
					subjectParticleBarycentricCoordinatesWrtOriginTriangleChanged = mesh.GetTriangles().ToBarycentricCoordinates(
						mParticles.GetPosition(state.PrimaryParticleState.ParticleIndex),
						*mCurrentOriginTriangle,
						mesh.GetVertices());
				}

				physicsParticleProbe.emplace(mParticles.GetVelocity(state.PrimaryParticleState.ParticleIndex));
			}

			if (state.SecondaryParticleState.has_value()
				&& state.SecondaryParticleState->ParticleIndex == *mCurrentlySelectedParticle
				&& state.SecondaryParticleState->ConstrainedState.has_value())
			{
				constrainedRegimeParticleProbe.emplace(
					state.SecondaryParticleState->ConstrainedState->CurrentTriangle,
					state.SecondaryParticleState->ConstrainedState->CurrentTriangleBarycentricCoords);

				if (mCurrentOriginTriangle.has_value())
				{
					subjectParticleBarycentricCoordinatesWrtOriginTriangleChanged = mesh.GetTriangles().ToBarycentricCoordinates(
						mParticles.GetPosition(state.SecondaryParticleState->ParticleIndex),
						*mCurrentOriginTriangle,
						mesh.GetVertices());
				}

				physicsParticleProbe.emplace(mParticles.GetVelocity(state.SecondaryParticleState->ParticleIndex));
			}
		}
	}

	mEventDispatcher.OnSubjectParticleConstrainedRegimeUpdated(constrainedRegimeParticleProbe);
	mEventDispatcher.OnSubjectParticleBarycentricCoordinatesWrtOriginTriangleChanged(subjectParticleBarycentricCoordinatesWrtOriginTriangleChanged);
	mEventDispatcher.OnSubjectParticlePhysicsUpdated(physicsParticleProbe);
}

void Npcs::Render(RenderContext & renderContext)
{
	//
	// Particles
	//

	renderContext.UploadParticlesStart();

	for (auto const i : *this)
	{
		auto const & state = mStateBuffer[i];

		RenderParticle(state.PrimaryParticleState, renderContext);

		if (state.SecondaryParticleState.has_value())
		{
			RenderParticle(*state.SecondaryParticleState, renderContext);
		}
	}

	renderContext.UploadParticlesEnd();

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
		vec2f sourcePosition;

		// If the trajectory's particle is being ray-traced, use its current position in its quest to the trajectory;
		// otherwise, use its real position
		if (IsParticleBeingRayTraced(mCurrentParticleTrajectory->ParticleIndex))
		{
			sourcePosition = mSimulationStepState.TrajectoryState->CurrentPosition;
		}
		else
		{
			sourcePosition = mParticles.GetPosition(mCurrentParticleTrajectory->ParticleIndex);
		}

		renderContext.UploadParticleTrajectory(
			sourcePosition,
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

			if (state.SecondaryParticleState.has_value()
				&& state.SecondaryParticleState->ParticleIndex == *mCurrentlySelectedParticle
				&& state.SecondaryParticleState->ConstrainedState.has_value()
				&& triangleIndex == state.SecondaryParticleState->ConstrainedState->CurrentTriangle)
			{
				return true;
			}
		}
	}

	return false;
}

////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////

void Npcs::RotateParticleWithMesh(
	StateType::NpcParticleStateType const & npcParticleState,
	vec2f const & centerPos,
	float cosAngle,
	float sinAngle,
	Mesh const & mesh)
{
	vec2f newPosition;

	if (npcParticleState.ConstrainedState.has_value())
	{
		// Simply set position from current bary coords

		newPosition = mesh.GetTriangles().FromBarycentricCoordinates(
			npcParticleState.ConstrainedState->CurrentTriangleBarycentricCoords,
			npcParticleState.ConstrainedState->CurrentTriangle,
			mesh.GetVertices());
	}
	else
	{
		// Rotate particle

		vec2f const centeredPos = mParticles.GetPosition(npcParticleState.ParticleIndex) - centerPos;
		vec2f const rotatedPos = vec2f(
			centeredPos.x * cosAngle - centeredPos.y * sinAngle,
			centeredPos.x * sinAngle + centeredPos.y * cosAngle);

		newPosition = rotatedPos + centeredPos;
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

	// Render shadow at end of trajectory, it this particle is being ray-traced
	if (IsParticleBeingRayTraced(particleState.ParticleIndex))
	{
		renderContext.UploadParticle(
			mSimulationStepState.TrajectoryState->CurrentPosition,
			mParticles.GetRenderColor(particleState.ParticleIndex),
			0.5f);
	}
}
