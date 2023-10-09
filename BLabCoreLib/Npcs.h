/***************************************************************************************
* Original Author:      Gabriele Giuseppini
* Created:              2023-10-06
* Copyright:            Gabriele Giuseppini  (https://github.com/GabrieleGiuseppini)
***************************************************************************************/
#pragma once

#include "BLabTypes.h"
#include "ElementIndexRangeIterator.h"
#include "EventDispatcher.h"
#include "LabParameters.h"
#include "Mesh.h"
#include "NpcParticles.h"
#include "RenderContext.h"
#include "Vectors.h"

#include <optional>
#include <vector>

class Npcs final
{
public:

	enum class NpcType
	{
		Furniture,
		Human
	};

	struct StateType final
	{
		enum class RegimeType
		{
			Constrained,
			Free
		};

		RegimeType Regime;

		struct NpcParticleStateType final
		{
			ElementIndex ParticleIndex;

			struct ConstrainedStateType final
			{
				ElementIndex CurrentTriangle;

				vec3f CurrentTriangleBarycentricCoords;

				// Not in FS - only used because of step-by-step ray tracing
				vec2f MeshRelativeVelocity; // Velocity of particle (as in velocity buffer), but relative to mesh at the moment velocity was calculated

				ConstrainedStateType(
					ElementIndex currentTriangle,
					vec3f const & currentTriangleBarycentricCoords)
					: CurrentTriangle(currentTriangle)
					, CurrentTriangleBarycentricCoords(currentTriangleBarycentricCoords)
					, MeshRelativeVelocity(vec2f::zero())
				{}
			};

			std::optional<ConstrainedStateType> ConstrainedState;

			NpcParticleStateType(
				ElementIndex particleIndex,
				std::optional<ConstrainedStateType> && constrainedState)
				: ParticleIndex(particleIndex)
				, ConstrainedState(std::move(constrainedState))
			{}
		};

		NpcParticleStateType PrimaryParticleState; // e.g. feet
		std::optional<NpcParticleStateType> SecondaryParticleState; // e.g. head

		StateType(
			RegimeType regime,
			NpcParticleStateType && primaryParticleState,
			std::optional<NpcParticleStateType> && secondaryParticleState)
			: Regime(regime)
			, PrimaryParticleState(std::move(primaryParticleState))
			, SecondaryParticleState(std::move(secondaryParticleState))
		{}
	};

	Npcs(
		EventDispatcher & eventDispatcher,
		bool isGravityEnabled)
		: mEventDispatcher(eventDispatcher)
		// Container
		, mStateBuffer()
		, mParticles(LabParameters::MaxNpcs * LabParameters::MaxParticlesPerNpc)
		// Parameters
		, mGravityGate(isGravityEnabled ? 1.0f : 0.0f)
	{}

	void Add(
		NpcType npcType,
		vec2f const & primaryPosition,
		Mesh const & mesh);

	void MoveParticleBy(
		ElementIndex particleIndex,
		vec2f const & offset,
		Mesh const & mesh);

	void RotateParticlesWithMesh(
		vec2f const & centerPos,
		float cosAngle,
		float sinAngle,
		Mesh const & mesh);

	void OnVertexMoved(Mesh const & mesh);

	void Update(
		Mesh const & mesh,
		LabParameters const & labParameters);

	void Render(RenderContext & renderContext);

	inline element_index_range_iterator begin() const noexcept
	{
		return element_index_range_iterator(0u);
	}

	inline element_index_range_iterator end() const noexcept
	{
		return element_index_range_iterator(static_cast<ElementCount>(mStateBuffer.size()));
	}

	StateType const & GetState(ElementIndex npcIndex) const
	{
		return mStateBuffer[npcIndex];
	}

	StateType & GetState(ElementIndex npcIndex)
	{
		return mStateBuffer[npcIndex];
	}

	NpcParticles const & GetParticles() const
	{
		return mParticles;
	}

	NpcParticles & GetParticles()
	{
		return mParticles;
	}

	bool IsGravityEnabled() const
	{
		return mGravityGate != 0.0f;
	}

	void SetGravityEnabled(bool isEnabled)
	{
		mGravityGate = isEnabled ? 1.0f : 0.0f;
	}

	//
	// Probing
	//

	void SelectParticle(ElementIndex particleIndex)
	{
		mCurrentlySelectedParticle = particleIndex;
	}	

	std::optional<ElementIndex> GetCurrentOriginTriangle() const
	{
		return mCurrentOriginTriangle;
	}

	void SelectOriginTriangle(ElementIndex triangleIndex)
	{
		mCurrentOriginTriangle = triangleIndex;
	}

	void ResetOriginTriangle()
	{
		mCurrentOriginTriangle.reset();
	}

	void NotifyParticleTrajectory(
		ElementIndex particleIndex,
		vec2f const & targetPosition)
	{
		mCurrentParticleTrajectoryNotification.emplace(
			particleIndex,
			targetPosition);

		mCurrentParticleTrajectory.reset();
	}

	void SetParticleTrajectory(
		ElementIndex particleIndex, 
		vec2f const & targetPosition)
	{
		mCurrentParticleTrajectory.emplace(
			particleIndex,
			targetPosition);

		mCurrentParticleTrajectoryNotification.reset();
	}

private:

	void RotateParticleWithMesh(
		StateType::NpcParticleStateType const & npcParticleState,
		vec2f const & centerPos,
		float cosAngle,
		float sinAngle,
		Mesh const & mesh);

	void RenderParticle(
		StateType::NpcParticleStateType const & particleState,
		RenderContext & renderContext);

	bool IsTriangleConstrainingCurrentlySelectedParticle(ElementIndex triangleIndex) const;

private:

	//
	// Simulation
	//

	// This struct maintains the state of a simulation step, so that we can run the whole step a sub-step at a time
	struct SimulationStepStateType final
	{
		// The NPC we are simulating
		ElementIndex CurrentNpcIndex;
		
		// Whether the current NPC's particle being simulated is the primary or the secondary one
		bool CurrentIsPrimaryParticle; 

		struct TrajectoryStateType
		{
			vec2f SourcePosition;
			vec2f TargetPosition;

			struct ConstrainedStateType
			{
				vec3f TargetPositionCurrentTriangleBarycentricCoords;
				vec2f MeshDisplacement;

				ConstrainedStateType(
					vec3f const & targetPositionCurrentTriangleBarycentricCoords,
					vec2f const & meshDisplacement)
					: TargetPositionCurrentTriangleBarycentricCoords(targetPositionCurrentTriangleBarycentricCoords)
					, MeshDisplacement(meshDisplacement)
				{}
			};

			std::optional<ConstrainedStateType> ConstrainedState; // Always set when in constrained state; updated when current triangle changes

			vec2f CurrentPosition;

			TrajectoryStateType(
				vec2f const & sourcePosition,
				vec2f const & targetPosition,
				std::optional<ConstrainedStateType> constrainedState)
				: SourcePosition(sourcePosition)
				, TargetPosition(targetPosition)
				, ConstrainedState(std::move(constrainedState))
				, CurrentPosition(sourcePosition)
			{}
		};

		// When set, we have a trajectory target for the particle; when not set, we have to calculate a target
		std::optional<TrajectoryStateType> TrajectoryState; 

		SimulationStepStateType()
			: CurrentNpcIndex(0)
			, CurrentIsPrimaryParticle(true)
			, TrajectoryState()
		{}
	};

	SimulationStepStateType mSimulationStepState;

	// TODOHERE

	struct CalculatedTrajectoryTarget final
	{
		vec2f Position;

		std::optional<SimulationStepStateType::TrajectoryStateType::ConstrainedStateType> ConstrainedStateInfo; // Returned when in constrained state

		CalculatedTrajectoryTarget(
			vec2f const & position,
			std::optional<SimulationStepStateType::TrajectoryStateType::ConstrainedStateType> constrainedStateInfo)
			: Position(position)
			, ConstrainedStateInfo(std::move(constrainedStateInfo))
		{}
	};

	CalculatedTrajectoryTarget CalculateTrajectoryTarget(
		ElementIndex particleIndex,
		LabParameters const & labParameters) const;

	struct FinalParticleState final
	{
		vec2f Position;
		std::optional<vec2f> Velocity;

		FinalParticleState(
			vec2f const & position,
			std::optional<vec2f> velocity)
			: Position(position)
			, Velocity(std::move(velocity))
		{}
	};

	std::optional<FinalParticleState> UpdateParticleTrajectoryTrace(
		Npcs::StateType::NpcParticleStateType & particleState,
		LabParameters const & labParameters);

	StateType MaterializeNpcState(
		ElementIndex npcIndex,
		Mesh const & mesh) const;

	StateType::NpcParticleStateType MaterializeParticleState(
		vec2f const & position,
		ElementIndex particleIndex,
		Mesh const & mesh) const;

	void ResetSimulationStepState();

	bool IsParticleBeingRayTraced(ElementIndex particleIndex) const
	{
		return mSimulationStepState.TrajectoryState.has_value()
			&& ((mSimulationStepState.CurrentIsPrimaryParticle && mStateBuffer[mSimulationStepState.CurrentNpcIndex].PrimaryParticleState.ParticleIndex == particleIndex)
				|| (!mSimulationStepState.CurrentIsPrimaryParticle && mStateBuffer[mSimulationStepState.CurrentNpcIndex].SecondaryParticleState->ParticleIndex == particleIndex));
	}

private:

	EventDispatcher & mEventDispatcher;

	//
	// Container
	//

	std::vector<StateType> mStateBuffer;

	NpcParticles mParticles;

	//
	// Probing
	//

	std::optional<ElementIndex> mCurrentlySelectedParticle;
	std::optional<ElementIndex> mCurrentOriginTriangle;

	std::optional<ParticleTrajectory> mCurrentParticleTrajectory;
	std::optional<ParticleTrajectory> mCurrentParticleTrajectoryNotification;


	//
	// Simulation parameters
	//

	float mGravityGate;
};
