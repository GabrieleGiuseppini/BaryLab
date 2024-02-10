/***************************************************************************************
* Original Author:      Gabriele Giuseppini
* Created:              2023-10-06
* Copyright:            Gabriele Giuseppini  (https://github.com/GabrieleGiuseppini)
***************************************************************************************/
#pragma once

#include "BLabTypes.h"
#include "BarycentricCoords.h"
#include "ElementIndexRangeIterator.h"
#include "EventDispatcher.h"
#include "LabParameters.h"
#include "Mesh.h"
#include "NpcParticles.h"
#include "RenderContext.h"
#include "StrongTypeDef.h"
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

				bcoords3f CurrentTriangleBarycentricCoords;

				ElementIndex CurrentVirtualEdgeElementIndex; // When set, we are "conceptually" along this edge - might not be really the case e.g. when we're at a vertex

				vec2f MeshRelativeVelocity; // Velocity of particle (as in velocity buffer), but relative to mesh at the moment velocity was calculated

				ConstrainedStateType(
					ElementIndex currentTriangle,
					bcoords3f const & currentTriangleBarycentricCoords)
					: CurrentTriangle(currentTriangle)
					, CurrentTriangleBarycentricCoords(currentTriangleBarycentricCoords)
					, CurrentVirtualEdgeElementIndex(NoneElementIndex)
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

			float CurrentEquilibriumSoftTerminationDecision; // [0.0f, 1.0f]

			float CurrentWalkMagnitude; // [0.0f, 1.0f]
			float TargetWalkMagnitude; // [0.0f, 1.0f]

			float CurrentWalkFlipDecision; // [0.0f, 1.0f]
			float TargetWalkFlipDecision; // [0.0f, 1.0f]

			// Updated in any state, but reset when starting to walk
			// This is *mesh-relative* when we're in constrained state
			float TotalEdgeTraveledSinceWalkStart; // [0.0f, +INF]

			// Animation

			float LegRightAngle;
			float LegLeftAngle;

			vec2f TopPoint;
			vec2f NeckPoint;
			vec2f CrotchPoint;
			vec2f LegRightPoint;			
			vec2f LegLeftPoint;

			HumanNpcStateType(
				BehaviorType initialBehavior,
				float currentStateValue,
				float targetStateValue)
				: CurrentBehavior(initialBehavior)
				, CurrentFaceDirectionX(1.0f) // Futurework: randomize
				, CurrentStateValue(currentStateValue)
				, TargetStateValue(targetStateValue)
				, CurrentEquilibriumSoftTerminationDecision(0.0f)
				, CurrentWalkMagnitude(0.0f)
				, TargetWalkMagnitude(0.0f)
				, CurrentWalkFlipDecision(0.0f)
				, TargetWalkFlipDecision(0.0f)
				, TotalEdgeTraveledSinceWalkStart(0.0f)
				// Animation
				, LegRightAngle(0.0f)
				, LegLeftAngle(0.0f)
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

public:

	Npcs(
		EventDispatcher & eventDispatcher,
		bool isGravityEnabled)
		: mEventDispatcher(eventDispatcher)
		// Container
		, mStateBuffer()
		, mParticles(LabParameters::MaxNpcs * LabParameters::MaxParticlesPerNpc)
		// Parameters
		, mGravityGate(isGravityEnabled ? 1.0f : 0.0f)
		, mNpcRenderMode(NpcRenderMode::Limbs)
	{}

	void Add(
		NpcType npcType,
		vec2f primaryPosition,
		std::optional<vec2f> secondaryPosition,
		StructuralMaterialDatabase const & materialDatabase,
		Mesh const & mesh,
		LabParameters const & labParameters);

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

	NpcRenderMode GetNpcRenderMode() const
	{
		return mNpcRenderMode;
	}

	void SetNpcRenderMode(NpcRenderMode value)
	{
		mNpcRenderMode = value;
	}

	//
	// Probing
	//

	void SelectParticle(
		ElementIndex particleIndex,
		Mesh const & mesh)
	{
		mCurrentlySelectedParticle = particleIndex;
		Publish(mesh);
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

	bool IsEdgeHostingCurrentlySelectedParticle(ElementIndex edgeIndex) const;

public:

	static bool IsEdgeFloorToParticle(
		ElementIndex edgeElementIndex,
		ElementIndex triangleElementIndex,
		Mesh const & mesh)
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

	static bool DoesFloorSeparateFromPrimaryParticle(
		vec2f const & primaryParticlePosition,
		vec2f const & secondaryParticlePosition,
		ElementIndex edgeElementIndex,
		Mesh const & mesh)
	{
		vec2f const aPos = mesh.GetEdges().GetEndpointAPosition(edgeElementIndex, mesh.GetVertices());
		vec2f const bPos = mesh.GetEdges().GetEndpointBPosition(edgeElementIndex, mesh.GetVertices());
		vec2f const & p1Pos = primaryParticlePosition;
		vec2f const & p2Pos = secondaryParticlePosition;

		// ((y1−y2)(ax−x1)+(x2−x1)(ay−y1)) * ((y1−y2)(bx−x1)+(x2−x1)(by−y1)) < 0
		float const magic = ((aPos.y - bPos.y) * (p1Pos.x - aPos.x) + (bPos.x - aPos.x) * (p1Pos.y - aPos.y))
			* ((aPos.y - bPos.y) * (p2Pos.x - aPos.x) + (bPos.x - aPos.x) * (p2Pos.y - aPos.y));

		return magic < -0.0001f;
	}

	static bool IsAtTarget(float currentValue, float targetValue)
	{
		return std::abs(targetValue - currentValue) < 0.01f;
	}

	static bool IsOnFloorEdge(
		Npcs::StateType::NpcParticleStateType::ConstrainedStateType const & constrainedState,
		Mesh const & mesh)
	{
		auto const & baryCoords = constrainedState.CurrentTriangleBarycentricCoords;
		auto const & triangleIndex = constrainedState.CurrentTriangle;

		if (baryCoords[0] == 0.0f
			&& IsEdgeFloorToParticle(
				mesh.GetTriangles().GetSubEdges(triangleIndex).EdgeIndices[1],
				triangleIndex,
				mesh))
		{
			return true;
		}

		if (baryCoords[1] == 0.0f
			&& IsEdgeFloorToParticle(
				mesh.GetTriangles().GetSubEdges(triangleIndex).EdgeIndices[2],
				triangleIndex,
				mesh))
		{
			return true;
		}

		if (baryCoords[2] == 0.0f
			&& IsEdgeFloorToParticle(
				mesh.GetTriangles().GetSubEdges(triangleIndex).EdgeIndices[0],
				triangleIndex,
				mesh))
		{
			return true;
		}

		return false;
	}

	// Head->Feet
	static vec2f CalculateHumanVector(ElementIndex primaryParticleIndex, ElementIndex secondaryParticleIndex, NpcParticles const & particles)
	{
		return particles.GetPosition(primaryParticleIndex) - particles.GetPosition(secondaryParticleIndex);
	}

	static float CalculateVerticalAlignment(vec2f const & humanVector)
	{
		return humanVector.normalise().dot(LabParameters::GravityDir);
	}

	static float CalculateVerticalAlignment(ElementIndex primaryParticleIndex, ElementIndex secondaryParticleIndex, NpcParticles const & particles)
	{
		return CalculateVerticalAlignment(CalculateHumanVector(primaryParticleIndex, secondaryParticleIndex, particles));
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

	void Publish(Mesh const & mesh);

private:

	//
	// Simulation 
	//

	void UpdateNpcs(
		Mesh const & mesh,
		LabParameters const & labParameters);

	void UpdateNpcParticle(
		StateType & npc,
		bool isPrimaryParticle,
		Mesh const & mesh,
		LabParameters const & labParameters);

	void CalculateNpcParticlePreliminaryForces(
		StateType const & npc,
		bool isPrimaryParticle,
		LabParameters const & labParameters);

	vec2f CalculateNpcParticleDefinitiveForces(
		StateType const & npc,
		bool isPrimaryParticle,
		float particleMass,
		LabParameters const & labParameters) const;

	void UpdateNpcParticle_Free(
		StateType::NpcParticleStateType & particle,
		vec2f const & startPosition,
		vec2f const & endPosition,
		NpcParticles & particles,
		LabParameters const & labParameters) const;

	// Returns total edge traveled (in step), and isStop
	std::tuple<float, bool> UpdateNpcParticle_ConstrainedNonInertial(
		StateType & npc,
		bool isPrimaryParticle,
		int edgeOrdinal,
		vec2f const & edgeDir,
		vec2f const & particleStartAbsolutePosition,				
		vec2f const & trajectoryStartAbsolutePosition,
		vec2f const & flattenedTrajectoryEndAbsolutePosition,
		bcoords3f flattenedTrajectoryEndBarycentricCoords,
		vec2f const & flattenedTrajectory,
		float edgeTraveledPlanned,
		vec2f const meshVelocity,
		float dt,
		NpcParticles & particles,
		Mesh const & mesh,
		LabParameters const & labParameters) const;

	float UpdateNpcParticle_ConstrainedInertial(
		StateType & npc,
		bool isPrimaryParticle,
		vec2f const & particleStartAbsolutePosition,
		bcoords3f const segmentTrajectoryStartBarycentricCoords,
		vec2f const & segmentTrajectoryEndAbsolutePosition,
		bcoords3f segmentTrajectoryEndBarycentricCoords,
		vec2f const meshVelocity,
		float segmentDt,
		NpcParticles & particles,
		Mesh const & mesh,
		LabParameters const & labParameters) const;

	struct NavigateVertexOutcome
	{
		enum class OutcomeType
		{
			CompletedNavigation,
			EncounteredFloor,
			ConvertedToFree
		};

		OutcomeType Type;

		int EncounteredFloorEdgeOrdinal; // In particle's current triangle

		static NavigateVertexOutcome MakeCompletedNavigationOutcome()
		{
			return NavigateVertexOutcome(OutcomeType::CompletedNavigation, -1);
		}

		static NavigateVertexOutcome MakeEncounteredFloorOutcome(int encounteredFloorEdgeOrdinal)
		{
			return NavigateVertexOutcome(OutcomeType::EncounteredFloor, encounteredFloorEdgeOrdinal);
		}

		static NavigateVertexOutcome MakeConvertedToFreeOutcome()
		{
			return NavigateVertexOutcome(OutcomeType::ConvertedToFree, -1);
		}

	private:

		NavigateVertexOutcome(
			OutcomeType type,
			int encounteredFloorEdgeOrdinal)
			: Type(type)
			, EncounteredFloorEdgeOrdinal(encounteredFloorEdgeOrdinal)
		{}
	};

	NavigateVertexOutcome NavigateVertex(
		StateType & npc,
		bool isPrimaryParticle,
		int vertexOrdinal,
		vec2f const & particleStartAbsolutePosition,
		vec2f const & trajectoryStartAbsolutePosition,
		vec2f const & trajectoryEndAbsolutePosition,
		bcoords3f trajectoryEndBarycentricCoords,
		bool isInitialStateUnknown,
		NpcParticles & particles,
		Mesh const & mesh,
		LabParameters const & labParameters) const;

	inline void BounceConstrainedNpcParticle(
		StateType & npc,
		bool isPrimaryParticle,
		vec2f const & trajectory,
		vec2f const & bouncePosition,
		vec2f const & bounceEdgeNormal,
		vec2f const meshVelocity,
		float dt,
		NpcParticles & particles,
		LabParameters const & labParameters) const;

	void OnImpact(
		vec2f const & impactVector,
		vec2f const & bounceEdgeNormal,
		StateType & npc,
		bool isPrimaryParticle) const;

	void UpdateNpcAnimation(
		StateType & npc,
		bool isPrimaryParticle,
		Mesh const & mesh,
		LabParameters const & labParameters);

private:

	//
	// Helpers
	//

	StateType MaterializeNpcState(
		ElementIndex npcIndex,
		Mesh const & mesh) const;

	std::optional<StateType::NpcParticleStateType::ConstrainedStateType> CalculateParticleConstrainedState(
		vec2f const & position,
		Mesh const & mesh) const;

private:

	//
	// Human simulation
	//

	StateType::HumanNpcStateType InitializeHuman(
		StateType::NpcParticleStateType const & primaryParticleState,
		StateType::NpcParticleStateType const & secondaryParticleState) const;

	void UpdateHuman(
		StateType::HumanNpcStateType & humanState,
		StateType::NpcParticleStateType const & primaryParticleState,
		StateType::NpcParticleStateType const & secondaryParticleState,
		Mesh const & mesh,
		LabParameters const & labParameters);

	bool CheckAndMaintainHumanEquilibrium(
		ElementIndex primaryParticleIndex,
		ElementIndex secondaryParticleIndex,
		bool isRisingState,
		bool isOnEdge,
		NpcParticles & particles,
		LabParameters const & labParameters);

	void RunWalkingHumanStateMachine(
		StateType::HumanNpcStateType & humanState,
		StateType::NpcParticleStateType const & primaryParticleState,
		Mesh const & mesh,
		LabParameters const & labParameters);

	using DoImmediate = StrongTypedBool<struct _DoImmediate>;

	void FlipHumanWalk(
		StateType::HumanNpcStateType & humanState, 
		DoImmediate doImmediate) const;

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
	NpcRenderMode mNpcRenderMode;
};
