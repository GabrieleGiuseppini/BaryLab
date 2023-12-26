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
#include "StructuralMaterialDatabase.h"
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
				Constrained_KnockedOut, // Does nothing
				Constrained_Rising, // Tries to stand up (appliying torque)
				Constrained_Equilibrium, // Stands up; continues to adjust alignment with torque
				Constrained_Walking, // Walks; continues to adjust alignment with torque

				Free_KnockedOut // Does nothing
			};

			BehaviorType CurrentBehavior;

			float CurrentFaceDirectionX; // [-1.0f, 1.0f]

			float CurrentStateValue;
			float TargetStateValue;

			float CurrentWalkingMagnitude; // [0.0f -> 1.0f]

			HumanNpcStateType(
				BehaviorType initialBehavior,
				float currentStateValue,
				float targetStateValue)
				: CurrentBehavior(initialBehavior)
				, CurrentFaceDirectionX(1.0f) // Futurework: randomize
				, CurrentStateValue(currentStateValue)
				, TargetStateValue(targetStateValue)
				, CurrentWalkingMagnitude(0.0f)
			{}

			void TransitionToState(
				BehaviorType behavior,
				float currentStateValue,
				float targetStateValue)
			{
				CurrentBehavior = behavior;
				CurrentStateValue = currentStateValue;
				TargetStateValue = targetStateValue;
			}
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
		, mGravityGate(isGravityEnabled ? 1.0f : 0.0f)
	{}

	void Add(
		NpcType npcType,
		vec2f const & primaryPosition,
		StructuralMaterialDatabase const & materialDatabase,
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

	void Publish(Mesh const & mesh);

private:

	//
	// Simulation 
	//

	void UpdateNpcs(
		Mesh const & mesh,
		LabParameters const & labParameters);

	struct DipoleArg final
	{
		StateType::NpcParticleStateType & OtherParticle;
		StateType::DipolePropertiesType & DipoleProperties;

		DipoleArg(
			StateType::NpcParticleStateType & otherParticle,
			StateType::DipolePropertiesType & dipoleProperties)
			: OtherParticle(otherParticle)
			, DipoleProperties(dipoleProperties)
		{}
	};

	void UpdateNpcParticle(
		StateType::NpcParticleStateType & particle,
		std::optional<DipoleArg> const & dipoleArg,
		bool isPrimaryParticle,
		StateType & npc,
		Mesh const & mesh,
		LabParameters const & labParameters);

	vec2f CalculateNpcParticlePhysicalForces(
		StateType::NpcParticleStateType & particle,
		float particleMass,
		LabParameters const & labParameters) const;

	void UpdateNpcParticle_Free(
		StateType::NpcParticleStateType & particle,
		vec2f const & startPosition,
		vec2f const & endPosition,
		NpcParticles & particles) const;

	std::optional<float> UpdateNpcParticle_ConstrainedNonInertial(
		StateType::NpcParticleStateType & particle,
		std::optional<DipoleArg> const & dipoleArg,
		bool const isPrimaryParticle,
		int edgeOrdinal,
		vec2f const & edgeDir,
		vec2f const & particleStartAbsolutePosition,				
		vec2f const & trajectoryStartAbsolutePosition,
		vec3f trajectoryEndBarycentricCoords,
		vec2f const & flattenedTrajectory,
		float edgePhysicalTraveledPlanned,
		float edgeWalkedPlanned,
		vec2f const meshVelocity,
		float dt,
		NpcParticles & particles,
		Mesh const & mesh,
		LabParameters const & labParameters) const;

	void UpdateNpcParticle_ConstrainedInertial(
		StateType::NpcParticleStateType & particle,
		std::optional<DipoleArg> const & dipoleArg,
		bool const isPrimaryParticle,
		vec3f const trajectoryStartBarycentricCoords,
		vec3f trajectoryEndBarycentricCoords,
		vec2f const & totalEdgeWalkedActual,
		vec2f const meshVelocity,
		float dt,
		NpcParticles & particles,
		Mesh const & mesh,
		LabParameters const & labParameters) const;

	inline void BounceConstrainedNpcParticle(
		ElementIndex particleIndex,
		StateType::NpcParticleStateType::ConstrainedStateType & particleConstrainedState,
		vec2f const & trajectory,
		vec2f const & bouncePosition,
		vec2f const & bounceEdgeNormal,
		vec2f const meshVelocity,
		float dt,
		NpcParticles & particles,
		LabParameters const & labParameters) const;

private:

	//
	// Helpers
	//

	StateType MaterializeNpcState(
		ElementIndex npcIndex,
		NpcParticles & particles,
		Mesh const & mesh) const;

	std::optional<StateType::NpcParticleStateType::ConstrainedStateType> CalculateParticleConstrainedState(
		vec2f const & position,
		Mesh const & mesh) const;

	bool IsEdgeFloorToParticle(
		ElementIndex edgeElementIndex,
		ElementIndex triangleElementIndex,
		Mesh const & mesh) const
	{
		//
		// An edge is a floor for a given (constrained) particle if:
		// - It is a floor; AND
		// - The triangle is _not_ sealed, OR it _is_ sealed but crossing the edge would make the particle free
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
		vec2f const & primaryParticlePosition,
		vec2f const & secondaryParticlePosition,
		ElementIndex edgeElementIndex,
		Mesh const & mesh) const
	{
		vec2f const aPos = mesh.GetEdges().GetEndpointAPosition(edgeElementIndex, mesh.GetVertices());
		vec2f const bPos = mesh.GetEdges().GetEndpointBPosition(edgeElementIndex, mesh.GetVertices());
		vec2f const & p1Pos = primaryParticlePosition;
		vec2f const & p2Pos = secondaryParticlePosition;

		// ((y1−y2)(ax−x1)+(x2−x1)(ay−y1)) * ((y1−y2)(bx−x1)+(x2−x1)(by−y1)) < 0
		float const magic = ((aPos.y - bPos.y) * (p1Pos.x - aPos.x) + (bPos.x - aPos.x) * (p1Pos.y - aPos.y))
			* ((aPos.y - bPos.y) * (p2Pos.x - aPos.x) + (bPos.x - aPos.x) * (p2Pos.y - aPos.y));

		LogMessage("TODOTEST: Magic=", magic, " (ppos=", p1Pos, " spos=", p2Pos, ")");

		return magic < -0.0001f;
	}

private:

	//
	// Human simulation
	//

	StateType::HumanNpcStateType InitializeHuman(
		StateType::NpcParticleStateType const & primaryParticleState,
		StateType::NpcParticleStateType const & secondaryParticleState,
		NpcParticles & particles) const;

	void UpdateHuman(
		StateType::HumanNpcStateType & humanState,
		StateType::NpcParticleStateType const & primaryParticleState,
		StateType::NpcParticleStateType const & secondaryParticleState,
		Mesh const & mesh,
		LabParameters const & labParameters);

	bool MaintainAndCheckHumanEquilibrium(
		ElementIndex primaryParticleIndex,
		ElementIndex secondaryParticleIndex,
		NpcParticles & particles,
		LabParameters const & labParameters);

	void RunWalkingHumanStateMachine(
		StateType::HumanNpcStateType & humanState,
		StateType::NpcParticleStateType const & primaryParticleState,
		StateType::NpcParticleStateType const & secondaryParticleState,
		Mesh const & mesh,
		LabParameters const & labParameters);

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
