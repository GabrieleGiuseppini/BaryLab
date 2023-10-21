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

		struct NpcParticleStateType final
		{
			ElementIndex ParticleIndex;

			struct ConstrainedStateType final
			{
				ElementIndex CurrentTriangle;

				vec3f CurrentTriangleBarycentricCoords;

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

		struct DipolePropertiesType final
		{
			float DipoleLength;
			float MassFactor; // Net of MassAdjustment
			float BaseStiffnessCoefficient;

			DipolePropertiesType(
				float dipoleLength,
				float massFactor,
				float baseStiffnessCoefficient)
				: DipoleLength(dipoleLength)
				, MassFactor(massFactor)
				, BaseStiffnessCoefficient(baseStiffnessCoefficient)
			{}
		};

		struct DipoleStateType
		{
			NpcParticleStateType SecondaryParticleState; // e.g. head
			DipolePropertiesType DipoleProperties;

			DipoleStateType(
				NpcParticleStateType && secondaryParticleState,
				DipolePropertiesType const & dipoleProperties)
				: SecondaryParticleState(std::move(secondaryParticleState))
				, DipoleProperties(dipoleProperties)
			{}
		};

		struct HumanNpcStateType final
		{
			enum class BehaviorType
			{
				Constrained_KnockedOut,
				Constrained_Rising,
				Constrained_Equilibrium,
				Constrained_Walking,

				Free_KnockedOut
			};

			BehaviorType CurrentBehavior;

			float CurrentStateValue;
			float TargetStateValue;

			HumanNpcStateType(
				BehaviorType initialBehavior,
				float currentStateValue,
				float targetStateValue)
				: CurrentBehavior(initialBehavior)
				, CurrentStateValue(currentStateValue)
				, TargetStateValue(targetStateValue)
			{}
		};

		NpcType Type;

		RegimeType Regime;

		NpcParticleStateType PrimaryParticleState; // e.g. feet
		std::optional<DipoleStateType> DipoleState;

		std::optional<HumanNpcStateType> HumanNpcState;

		StateType(
			NpcType type,
			RegimeType regime,
			NpcParticleStateType && primaryParticleState,
			std::optional<DipoleStateType> && dipoleState,
			std::optional<HumanNpcStateType> && humanNpcState)
			: Type(type)
			, Regime(regime)
			, PrimaryParticleState(std::move(primaryParticleState))
			, DipoleState(std::move(dipoleState))
			, HumanNpcState(std::move(humanNpcState))
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
		, mIsStepByStepMode(false)
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

	bool IsAtBeginningOfSimulationStep() const
	{
		return mSimulationStepState.IsInitial();
	}

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

	bool GetIsStepByStepMode() const
	{
		return mIsStepByStepMode;
	}

	void SetIsStepByStepMode(bool isStepByStepMode)
	{
		mIsStepByStepMode = isStepByStepMode;
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

	bool IsTriangleConstrainingCurrentlySelectedParticle(ElementIndex triangleIndex) const;

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

			vec2f CurrentPosition; // This is the only position we modify during the sub-steps; actual position is calculated when sub-steps are complete

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

		bool IsInitial() const
		{
			return CurrentNpcIndex == 0
				&& CurrentIsPrimaryParticle == true
				&& !TrajectoryState.has_value();
		}
	};

	SimulationStepStateType mSimulationStepState;

	void UpdateNpcs(
		Mesh const & mesh,
		LabParameters const & labParameters);

	struct TrajectoryTargetDipoleArg final
	{
		StateType::NpcParticleStateType & OtherParticle;
		StateType::DipolePropertiesType & DipoleProperties;

		TrajectoryTargetDipoleArg(
			StateType::NpcParticleStateType & otherParticle,
			StateType::DipolePropertiesType & dipoleProperties)
			: OtherParticle(otherParticle)
			, DipoleProperties(dipoleProperties)
		{}
	};

	struct CalculatedTrajectoryTargetRetVal final
	{
		vec2f Position;
		std::optional<SimulationStepStateType::TrajectoryStateType::ConstrainedStateType> ConstrainedStateInfo; // Returned when in constrained state

		CalculatedTrajectoryTargetRetVal(
			vec2f const & position,
			std::optional<SimulationStepStateType::TrajectoryStateType::ConstrainedStateType> constrainedStateInfo)
			: Position(position)
			, ConstrainedStateInfo(std::move(constrainedStateInfo))
		{}
	};

	CalculatedTrajectoryTargetRetVal CalculateTrajectoryTarget(
		StateType::NpcParticleStateType & particle,
		std::optional<TrajectoryTargetDipoleArg> const & dipoleArg,
		bool isPrimaryParticle,
		Mesh const & mesh,
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

	// When returns a final particle state, the simulation of this particle has completed
	std::optional<FinalParticleState> UpdateParticleTrajectoryTrace(
		Npcs::StateType::NpcParticleStateType & particleState,
		std::optional<TrajectoryTargetDipoleArg> const & dipoleArg,
		bool isPrimaryParticle,
		SimulationStepStateType::TrajectoryStateType & trajectoryState,
		Mesh const & mesh,
		LabParameters const & labParameters);

	StateType MaterializeNpcState(
		ElementIndex npcIndex,
		Mesh const & mesh) const;

	std::optional<StateType::NpcParticleStateType::ConstrainedStateType> CalculateParticleConstrainedState(
		vec2f const & position,
		Mesh const & mesh) const;

	void ResetSimulationStepState();

	bool IsParticleBeingRayTraced(ElementIndex particleIndex) const
	{
		return mSimulationStepState.TrajectoryState.has_value()
			&& ((mSimulationStepState.CurrentIsPrimaryParticle && mStateBuffer[mSimulationStepState.CurrentNpcIndex].PrimaryParticleState.ParticleIndex == particleIndex)
				|| (!mSimulationStepState.CurrentIsPrimaryParticle && mStateBuffer[mSimulationStepState.CurrentNpcIndex].DipoleState->SecondaryParticleState.ParticleIndex == particleIndex));
	}

	bool IsEdgeFloorToParticle(
		ElementIndex edgeElementIndex,
		ElementIndex triangleElementIndex,
		Mesh const & mesh) const
	{
		//
		// An edge is a floor for a given (constrained) particle if:
		// - It is a floor; AND
		// - The triangle is _not_ sealed, OR it _is_ sealed bu crossing the edge would make the particle free
		//

		if (mesh.GetEdges().GetSurfaceType(edgeElementIndex) != SurfaceType::Floor)
		{
			// Not even a floor
			return false;
		}

		bool const isSealedTriangle = 
			mesh.GetEdges().GetSurfaceType(mesh.GetTriangles().GetSubEdgeAIndex(triangleElementIndex)) == SurfaceType::Floor
			&& mesh.GetEdges().GetSurfaceType(mesh.GetTriangles().GetSubEdgeBIndex(triangleElementIndex)) == SurfaceType::Floor
			&& mesh.GetEdges().GetSurfaceType(mesh.GetTriangles().GetSubEdgeCIndex(triangleElementIndex)) == SurfaceType::Floor;

		if (!isSealedTriangle)
		{
			return true;
		}

		ElementIndex const oppositeTriangle = mesh.GetEdges().GetOppositeTriangle(edgeElementIndex, triangleElementIndex);
		if (oppositeTriangle == NoneElementIndex)
		{
			// Crossing this floor makes the particle free
			return true;
		}

		return false;
	}

	bool DoesFloorSeparateFromPrimaryParticle(
		ElementIndex primaryParticleIndex,
		ElementIndex secondaryParticleIndex,
		ElementIndex edgeElementIndex,
		Mesh const & mesh) const
	{
		vec2f const aPos = mesh.GetEdges().GetEndpointAPosition(edgeElementIndex, mesh.GetVertices());
		vec2f const bPos = mesh.GetEdges().GetEndpointBPosition(edgeElementIndex, mesh.GetVertices());
		vec2f const p1Pos = mParticles.GetPosition(primaryParticleIndex);
		vec2f const p2Pos = mParticles.GetPosition(secondaryParticleIndex);

		// ((y1−y2)(ax−x1)+(x2−x1)(ay−y1)) * ((y1−y2)(bx−x1)+(x2−x1)(by−y1)) < 0
		float const magic = ((aPos.y - bPos.y) * (p1Pos.x - aPos.x) + (bPos.x - aPos.x) * (p1Pos.y - aPos.y))
			* ((aPos.y - bPos.y) * (p2Pos.x - aPos.x) + (bPos.x - aPos.x) * (p2Pos.y - aPos.y));

		return magic < -0.0001f;
	}


	//
	// Human simulation
	//

	StateType::HumanNpcStateType InitializeHuman(
		StateType::NpcParticleStateType const & primaryParticleState,
		StateType::NpcParticleStateType const & secondaryParticleState,
		Mesh const & mesh) const;

	void UpdateHuman(
		StateType::HumanNpcStateType & humanState,
		StateType::NpcParticleStateType const & primaryParticleState,
		StateType::NpcParticleStateType const & secondaryParticleState,
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

	bool mIsStepByStepMode;

	float mGravityGate;
};
