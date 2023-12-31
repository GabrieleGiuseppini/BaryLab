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
	StructuralMaterialDatabase const & materialDatabase,
	Mesh const & mesh)
{
	assert(mParticles.GetParticleCount() < LabParameters::MaxNpcs);

	// Primary particle state

	ElementIndex const primaryParticleIndex = mParticles.GetParticleCount();

	StateType::NpcParticleStateType primaryParticleState = StateType::NpcParticleStateType(
		primaryParticleIndex,
		CalculateParticleConstrainedState(
			primaryPosition,
			mesh));

	// Add particles and eventual other states

	std::optional<StateType::DipoleStateType> dipoleState;
	std::optional<StateType::HumanNpcStateType> humanNpcState;

	switch (npcType)
	{
		case NpcType::Furniture:
		{
			auto const & material = materialDatabase.GetStructuralMaterial(StructuralMaterialDatabase::UniqueMaterialKeyType::Furniture);

			mParticles.Add(
				material.Mass,
				material.StaticFriction,
				material.KineticFriction,
				primaryPosition,
				material.RenderColor);

			break;
		}

		case NpcType::Human:
		{
			// Feet (primary)

			auto const & feetMaterial = materialDatabase.GetStructuralMaterial(StructuralMaterialDatabase::UniqueMaterialKeyType::HumanFeet);

			mParticles.Add(
				feetMaterial.Mass,
				feetMaterial.StaticFriction,
				feetMaterial.KineticFriction,
				primaryPosition,
				feetMaterial.RenderColor);

			// Head (secondary)

			ElementIndex const headParticleIndex = mParticles.GetParticleCount();

			auto const & headMaterial = materialDatabase.GetStructuralMaterial(StructuralMaterialDatabase::UniqueMaterialKeyType::HumanHead);

			vec2f const secondaryPosition = primaryPosition + vec2f(0.0f, 1.0f) * LabParameters::HumanNpcLength;

			mParticles.Add(
				headMaterial.Mass,
				headMaterial.StaticFriction,
				headMaterial.KineticFriction,
				secondaryPosition,
				headMaterial.RenderColor);

			StateType::NpcParticleStateType secondaryParticleState = StateType::NpcParticleStateType(
				headParticleIndex,
				CalculateParticleConstrainedState(
					secondaryPosition,
					mesh));

			humanNpcState = InitializeHuman(
				primaryParticleState,
				secondaryParticleState);

			float const massFactor =
				(feetMaterial.Mass * headMaterial.Mass)
				/ (feetMaterial.Mass + headMaterial.Mass);

			dipoleState.emplace(
				std::move(secondaryParticleState),
				StateType::DipolePropertiesType(
					LabParameters::HumanNpcLength,
					massFactor,
					1.0f));

			break;
		}
	}
	
	//
	// Calculate regime
	//

	auto const regime = primaryParticleState.ConstrainedState.has_value()
		? StateType::RegimeType::Constrained
		: (dipoleState.has_value() && dipoleState->SecondaryParticleState.ConstrainedState.has_value()) ? StateType::RegimeType::Constrained : StateType::RegimeType::Free;

	mStateBuffer.emplace_back(
		npcType,
		regime,
		std::move(primaryParticleState),
		std::move(dipoleState),
		std::move(humanNpcState));
}

void Npcs::MoveParticleBy(
	ElementIndex particleIndex,
	vec2f const & offset,
	Mesh const & mesh)
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
			state = MaterializeNpcState(n, mesh);
			break;
		}
	}

	//
	// Select particle
	//

	SelectParticle(particleIndex, mesh);

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

		if (state.DipoleState.has_value())
		{
			RotateParticleWithMesh(
				state.DipoleState->SecondaryParticleState,
				centerPos,
				cosAngle,
				sinAngle,
				mesh);
		}
	}
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

	Publish(mesh);
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

		if (state.DipoleState.has_value())
		{
			RenderParticle(state.DipoleState->SecondaryParticleState, renderContext);
		}
	}

	renderContext.UploadParticlesEnd();

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

void Npcs::Publish(Mesh const & mesh)
{
	std::optional<ConstrainedRegimeParticleProbe> constrainedRegimeParticleProbe;
	std::optional<vec3f> subjectParticleBarycentricCoordinatesWrtOriginTriangleChanged;
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
						subjectParticleBarycentricCoordinatesWrtOriginTriangleChanged = mesh.GetTriangles().ToBarycentricCoordinates(
							mParticles.GetPosition(state.PrimaryParticleState.ParticleIndex),
							*mCurrentOriginTriangle,
							mesh.GetVertices());
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
						subjectParticleBarycentricCoordinatesWrtOriginTriangleChanged = mesh.GetTriangles().ToBarycentricCoordinates(
							mParticles.GetPosition(state.DipoleState->SecondaryParticleState.ParticleIndex),
							*mCurrentOriginTriangle,
							mesh.GetVertices());
					}
				}

				physicsParticleProbe.emplace(mParticles.GetVelocity(state.DipoleState->SecondaryParticleState.ParticleIndex));
			}
		}
	}

	mEventDispatcher.OnSubjectParticleConstrainedRegimeUpdated(constrainedRegimeParticleProbe);
	mEventDispatcher.OnSubjectParticleBarycentricCoordinatesWrtOriginTriangleChanged(subjectParticleBarycentricCoordinatesWrtOriginTriangleChanged);
	mEventDispatcher.OnSubjectParticlePhysicsUpdated(physicsParticleProbe);
}