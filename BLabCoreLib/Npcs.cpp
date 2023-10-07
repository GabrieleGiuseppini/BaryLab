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
	assert(mParticles.GetElementCount() < LabParameters::MaxNpcs);

	//
	// Add primary particle
	//

	ElementIndex const primaryParticleIndex = mParticles.GetElementCount();
	mParticles.Add(primaryPosition, rgbaColor(0x60, 0x60, 0x60, 0xff));
	auto primaryParticleState = CalculateInitialParticleState(
		primaryPosition,
		primaryParticleIndex,
		mesh);

	auto regime = primaryParticleState.ConstrainedState ? StateType::RegimeType::Constrained : StateType::RegimeType::Free;

	//
	// Add secondary particle
	//

	std::optional<StateType::NpcParticleStateType> secondaryParticleState;

	if (npcType != NpcType::Furniture)
	{
		// TODO
		assert(false);
	}
	
	mStateBuffer.emplace_back(
		regime,
		std::move(primaryParticleState),
		std::move(secondaryParticleState));
}

void Npcs::Update(Mesh const & mesh)
{
	// TODO
	(void)mesh;

	//
	// Publish
	//

	std::optional<ConstrainedRegimeParticleProbe> constrainedRegimeParticleProbe;
	std::optional<vec3f> subjectParticleBarycentricCoordinatesWrtOriginTriangleChanged;
	std::optional<PhysicsParticleProbe> physicsParticleProbe;

	if (mCurrentlySelectedParticle.has_value())
	{
		for (auto const i : *this)
		{
			auto const & state = mStateBuffer[i];

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
		// TODOHERE: from mSimulationState, is current particle is particle of mCurrentParticleTrajectory
		if (mModel->GetParticles().GetState(mCurrentParticleTrajectory->ParticleIndex).TrajectoryState.has_value())
		{
			sourcePosition = mModel->GetParticles().GetState(mCurrentParticleTrajectory->ParticleIndex).TrajectoryState->CurrentPosition;
		}
		else
		{
			sourcePosition = mModel->GetParticles().GetPosition(mCurrentParticleTrajectory->ParticleIndex);
		}

		renderContext.UploadParticleTrajectory(
			sourcePosition,
			mCurrentParticleTrajectory->TargetPosition,
			rgbaColor(0x99, 0x99, 0x99, 0xff));
	}

	renderContext.UploadParticleTrajectoriesEnd();

}

////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////

void Npcs::RenderParticle(
	StateType::NpcParticleStateType const & particleState,
	RenderContext & renderContext)
{
	renderContext.UploadParticle(		
		mParticles.GetPosition(particleState.ParticleIndex),
		mParticles.GetRenderColor(particleState.ParticleIndex),
		1.0f);

	// Render end of trajectory
	if (particleState.TrajectoryState.has_value())
	{
		renderContext.UploadParticle(
			particleState.TrajectoryState->CurrentPosition,
			mParticles.GetRenderColor(particleState.ParticleIndex),
			0.5f);
	}
}

Npcs::StateType::NpcParticleStateType Npcs::CalculateInitialParticleState(
	vec2f const & position,
	ElementIndex particleIndex,
	Mesh const & mesh)
{
	std::optional<StateType::NpcParticleStateType::ConstrainedStateType> constrainedState;

	ElementIndex const triangleIndex = mesh.GetTriangles().FindContaining(position, mesh.GetVertices());
	if (triangleIndex != NoneElementIndex)
	{
		vec3f const barycentricCoords = mesh.GetTriangles().ToBarycentricCoordinatesFromWithinTriangle(
			position,
			triangleIndex,
			mesh.GetVertices());

		{
			assert(barycentricCoords[0] >= 0.0f && barycentricCoords[0] <= 1.0f);
			assert(barycentricCoords[1] >= 0.0f && barycentricCoords[1] <= 1.0f);
			assert(barycentricCoords[2] >= 0.0f && barycentricCoords[2] <= 1.0f);
		}

		constrainedState.emplace(
			triangleIndex,
			barycentricCoords);
	}

	return StateType::NpcParticleStateType(
		particleIndex,
		std::move(constrainedState));
}

bool Npcs::IsTriangleConstrainingCurrentlySelectedParticle(
	ElementIndex triangleIndex,
	Mesh const & mesh)
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
