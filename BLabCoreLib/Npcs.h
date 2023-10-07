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

			// Not in FS - only used because of step-by-step ray tracing
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

			// Not in FS - only used because of step-by-step ray tracing
			std::optional<TrajectoryStateType> TrajectoryState; // When set, we have a trajectory target; when not set, we have to calculate a target

			NpcParticleStateType(
				ElementIndex particleIndex,
				std::optional<ConstrainedStateType> && constrainedState)
				: ParticleIndex(particleIndex)
				, ConstrainedState(std::move(constrainedState))
				, TrajectoryState()
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
		, mGravityGate(0.0f)
	{}

	void Add(
		NpcType npcType,
		vec2f const & primaryPosition,
		Mesh const & mesh);

	void Update(Mesh const & mesh);

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

	bool IsTriangleConstrainingCurrentlySelectedParticle(
		ElementIndex triangleIndex,
		Mesh const & mesh);

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

private:

	void RenderParticle(
		StateType::NpcParticleStateType const & particleState,
		RenderContext & renderContext);

	//
	// Simulation
	//

	// This struct maintains the state of a simulation step, so that we can run the whole step a sub-step at a time
	struct SimulationStepState final
	{
		// The NPC we are simulating
		ElementIndex CurrentNpcIndex;

		// Whether the current NPC's particle being simulated is the primary or the secondary one
		bool CurrentIsPrimaryParticle;
	};

	StateType::NpcParticleStateType CalculateInitialParticleState(		
		vec2f const & position,
		ElementIndex particleIndex,
		Mesh const & mesh);

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
