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
#include "Log.h"
#include "Physics.h"
#include "RenderContext.h"
#include "StrongTypeDef.h"
#include "StructuralMaterialDatabase.h"
#include "Vectors.h"

#include <optional>
#include <vector>


template<typename... TArgs>
void LogNpcDebug(TArgs&&... args)
{
#ifdef IN_BARYLAB
	LogDebug(std::forward<TArgs>(args)...);
#else
	Logger::Instance.LogToNothing(std::forward<TArgs>(args)...);
#endif
}

namespace Physics {

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

				int CurrentVirtualEdgeOrdinal; // When set, we are "conceptually" along this edge - might not be really the case e.g. when we're at a vertex

				vec2f MeshRelativeVelocity; // Velocity of particle (as in velocity buffer), but relative to mesh (ship) at the moment velocity was calculated

				ConstrainedStateType(
					ElementIndex currentTriangle,
					bcoords3f const & currentTriangleBarycentricCoords)
					: CurrentTriangle(currentTriangle)
					, CurrentTriangleBarycentricCoords(currentTriangleBarycentricCoords)
					, CurrentVirtualEdgeOrdinal(-1)
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

			vec2f const & GetApplicableVelocity(NpcParticles const & particles) const
			{
				if (ConstrainedState.has_value())
				{
					return ConstrainedState->MeshRelativeVelocity;
				}
				else
				{
					return particles.GetVelocity(ParticleIndex);
				}
			}
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
				Constrained_Falling, // Clueless, does nothing; with feet on edge
				Constrained_Aerial, // Clueless, does nothing; with feet in air
				Constrained_KnockedOut, // Clueless, does nothing, not relevant where; waits until can rise

				Constrained_Rising, // Tries to stand up (appliying torque)
				Constrained_Equilibrium, // Stands up; continues to adjust alignment with torque
				Constrained_Walking, // Walks; continues to adjust alignment with torque

				Free_Aerial, // Does nothing
				Free_InWater, // Does nothing, but waits to swim
				Free_Swimming // Swims
			};

			BehaviorType CurrentBehavior;

			union BehaviorStateType
			{
				struct Constrained_FallingStateType
				{
					float ProgressToAerial;
					float ProgressToKnockedOut;

					void Reset()
					{
						ProgressToAerial = 0.0f;
						ProgressToKnockedOut = 0.0f;
					}
				} Constrained_Falling;

				struct Constrained_AerialStateType
				{
					float ProgressToFalling;
					float ProgressToKnockedOut;

					void Reset()
					{
						ProgressToFalling = 0.0f;
						ProgressToKnockedOut = 0.0f;
					}
				} Constrained_Aerial;

				struct Constrained_KnockedOutStateType
				{
					float ProgressToRising;
					float ProgressToAerial;

					void Reset()
					{
						ProgressToRising = 0.0f;
						ProgressToAerial = 0.0f;
					}
				} Constrained_KnockedOut;

				struct Constrained_RisingStateType
				{
					float CurrentSoftTerminationDecision; // [0.0f, 1.0f]

					void Reset()
					{
						CurrentSoftTerminationDecision = 0.0f;
					}
				} Constrained_Rising;

				struct Constrained_EquilibriumStateType
				{
					float ProgressToWalking;

					void Reset()
					{
						ProgressToWalking = 0.0f;
					}
				} Constrained_Equilibrium;

				struct Constrained_WalkingStateType
				{
					float CurrentWalkMagnitude; // [0.0f, 1.0f]

					float CurrentFlipDecision; // [0.0f, 1.0f]
					float TargetFlipDecision; // [0.0f, 1.0f]

					void Reset()
					{
						CurrentWalkMagnitude = 0.0f;
						CurrentFlipDecision = 0.0f;
						TargetFlipDecision = 0.0f;
					}
				} Constrained_Walking;

				struct Free_AerialStateType
				{
					void Reset()
					{
					}
				} Free_Aerial;

				struct Free_InWaterType
				{
					float ProgressToSwimming;

					void Reset()
					{
						ProgressToSwimming = 0.0f;
					}
				} Free_InWater;

				struct Free_SwimmingType
				{
					void Reset()
					{
					}
				} Free_Swimming;

			} CurrentBehaviorState;

			float CurrentStateTransitionSimulationTimestamp;
			float TotalDistanceTraveledOnEdgeSinceStateTransition; // [0.0f, +INF] - when we're constrained on an edge (e.g. walking)
			float TotalDistanceTraveledOffEdgeSinceStateTransition; // [0.0f, +INF] - when we're constrained off an edge or free

			float CurrentEquilibriumSoftTerminationDecision; // Cross-state

			float CurrentFaceOrientation; // [-1.0f, 0.0f, 1.0f]
			float CurrentFaceDirectionX; // [-1.0f, 0.0f, 1.0f]

			float PanicLevel; // [0.0f ... +INF)

			// Animation

			// Angles are CCW relative to vertical, regardless of where the NPC is looking towards (L/R)
			// (when we flip we pretend immediate mirroring of limbs from the point of view of the human,
			// so angles are independent from direction, and animation is smoother)

			// "Left" and "Right" are relative to screen when the NPC is looking at us
			// (so "right arm" is really its left arm)
			float RightLegAngle;
			float RightLegLengthMultiplier;
			float LeftLegAngle;
			float LeftLegLengthMultiplier;
			float RightArmAngle;
			float RightArmLengthMultiplier;
			float LeftArmAngle;
			float LeftArmLengthMultiplier;

			static float constexpr InitialArmAngle = Pi<float> / 2.0f * 0.3f;
			static float constexpr InitialLegAngle = 0.2f;

			HumanNpcStateType(
				BehaviorType initialBehavior,
				float currentSimulationTime)
				: CurrentEquilibriumSoftTerminationDecision(0.0f)
				, CurrentFaceOrientation(1.0f)
				, CurrentFaceDirectionX(0.0f)
				, PanicLevel(0.0f)
				// Animation
				, RightLegAngle(InitialLegAngle)
				, RightLegLengthMultiplier(1.0f)
				, LeftLegAngle(-InitialLegAngle)
				, LeftLegLengthMultiplier(1.0f)
				, RightArmAngle(InitialArmAngle)
				, RightArmLengthMultiplier(1.0f)
				, LeftArmAngle(-InitialArmAngle)
				, LeftArmLengthMultiplier(1.0f)
			{
				TransitionToState(initialBehavior, currentSimulationTime);
			}

			void TransitionToState(
				BehaviorType behavior,
				float currentSimulationTime)
			{
				CurrentBehavior = behavior;

				switch (behavior)
				{
					case BehaviorType::Constrained_Aerial:
					{
						CurrentBehaviorState.Constrained_Aerial.Reset();
						break;
					}

					case BehaviorType::Constrained_Equilibrium:
					{
						CurrentBehaviorState.Constrained_Equilibrium.Reset();
						break;
					}

					case BehaviorType::Constrained_Falling:
					{
						CurrentBehaviorState.Constrained_Falling.Reset();
						break;
					}

					case BehaviorType::Constrained_KnockedOut:
					{
						CurrentBehaviorState.Constrained_KnockedOut.Reset();
						break;
					}

					case BehaviorType::Constrained_Rising:
					{
						CurrentBehaviorState.Constrained_Rising.Reset();
						CurrentEquilibriumSoftTerminationDecision = 0.0f; // Start clean
						break;
					}

					case BehaviorType::Constrained_Walking:
					{
						CurrentBehaviorState.Constrained_Walking.Reset();
						break;
					}

					case BehaviorType::Free_Aerial:
					{
						CurrentBehaviorState.Free_Aerial.Reset();
						break;
					}

					case BehaviorType::Free_InWater:
					{
						CurrentBehaviorState.Free_InWater.Reset();
						break;
					}

					case BehaviorType::Free_Swimming:
					{
						CurrentBehaviorState.Free_Swimming.Reset();
						break;
					}
				}

				CurrentStateTransitionSimulationTimestamp = currentSimulationTime;
				TotalDistanceTraveledOnEdgeSinceStateTransition = 0.0f;
				TotalDistanceTraveledOffEdgeSinceStateTransition = 0.0f;
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
		Physics::World & parentWorld,
		EventDispatcher & eventDispatcher,
		LabParameters const & labParameters,
		bool isGravityEnabled)
		: mParentWorld(parentWorld)
		, mEventDispatcher(eventDispatcher)
		// Container
		, mStateBuffer()
		, mParticles(LabParameters::MaxNpcs * LabParameters::MaxParticlesPerNpc)
		// Parameters
		, mGravityGate(isGravityEnabled ? 1.0f : 0.0f)
		, mNpcRenderMode(NpcRenderMode::Limbs)
		, mCurrentHumanNpcBodyLengthAdjustment(labParameters.HumanNpcBodyLengthAdjustment)
	{}

	void Add(
		NpcType npcType,
		vec2f primaryPosition,
		std::optional<vec2f> secondaryPosition,
		float currentSimulationTime,
		StructuralMaterialDatabase const & materialDatabase,
		Ship const & ship,
		LabParameters const & labParameters);

	void SetPanicLevelForAllHumans(float panicLevel);

	void FlipHumanWalk(int npcIndex);

	void FlipHumanFrontBack(int npcIndex);

	void MoveParticleBy(
		ElementIndex particleIndex,
		vec2f const & offset,
		float currentSimulationTime,
		Ship const & ship);

	void RotateParticlesWithShip(
		vec2f const & centerPos,
		float cosAngle,
		float sinAngle,
		Ship const & ship);

	void OnVertexMoved(
		float currentSimulationTime,
		Ship const & ship);

	void Update(
		float currentSimulationTime,
		Ship const & ship,
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
		Ship const & ship)
	{
		mCurrentlySelectedParticle = particleIndex;
		Publish(ship);
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

	bool IsEdgeHostingCurrentlySelectedParticle(
		ElementIndex edgeIndex,
		Ship const & ship) const;

public:

	static bool IsEdgeFloorToParticle(
		ElementIndex triangleElementIndex,
		int edgeOrdinal,
		Ship const & ship)
	{
		//
		// An edge is a floor for a given (constrained) particle if:
		// - It is a floor; AND
		// - The triangle is _not_ sealed, OR it _is_ sealed but crossing the edge would make the particle free
		//

		if (ship.GetTriangles().GetSubEdgeSurfaceType(triangleElementIndex, edgeOrdinal) != SurfaceType::Floor)
		{
			// Not even a floor
			return false;
		}

		bool const isSealedTriangle =
			ship.GetTriangles().GetSubEdgeSurfaceType(triangleElementIndex, 0) == SurfaceType::Floor
			&& ship.GetTriangles().GetSubEdgeSurfaceType(triangleElementIndex, 1) == SurfaceType::Floor
			&& ship.GetTriangles().GetSubEdgeSurfaceType(triangleElementIndex, 2) == SurfaceType::Floor;

		if (!isSealedTriangle)
		{
			return true;
		}

		auto const & oppositeTriangleInfo = ship.GetTriangles().GetOppositeTriangle(triangleElementIndex, edgeOrdinal);
		if (oppositeTriangleInfo.TriangleElementIndex == NoneElementIndex || ship.GetTriangles().IsDeleted(oppositeTriangleInfo.TriangleElementIndex))
		{
			// Crossing this floor makes the particle free
			return true;
		}

		return false;
	}

	static bool DoesFloorSeparateFromPrimaryParticle(
		vec2f const & primaryParticlePosition,
		vec2f const & secondaryParticlePosition,
		ElementIndex triangleElementIndex,
		int edgeOrdinal,
		Ship const & ship)
	{
		ElementIndex const edgeElementIndex = ship.GetTriangles().GetSubEdges(triangleElementIndex).EdgeIndices[edgeOrdinal];
		vec2f const aPos = ship.GetEdges().GetEndpointAPosition(edgeElementIndex, ship.GetVertices());
		vec2f const bPos = ship.GetEdges().GetEndpointBPosition(edgeElementIndex, ship.GetVertices());
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
		Ship const & ship)
	{
		auto const & baryCoords = constrainedState.CurrentTriangleBarycentricCoords;
		auto const & triangleIndex = constrainedState.CurrentTriangle;

		if (baryCoords[0] == 0.0f
			&& IsEdgeFloorToParticle(
				triangleIndex,
				1,
				ship))
		{
			return true;
		}

		if (baryCoords[1] == 0.0f
			&& IsEdgeFloorToParticle(
				triangleIndex,
				2,
				ship))
		{
			return true;
		}

		if (baryCoords[2] == 0.0f
			&& IsEdgeFloorToParticle(
				triangleIndex,
				0,
				ship))
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

	void RotateParticleWithShip(
		StateType::NpcParticleStateType const & npcParticleState,
		vec2f const & centerPos,
		float cosAngle,
		float sinAngle,
		Ship const & ship);

	void RenderParticle(
		StateType::NpcParticleStateType const & particleState,
		RenderContext & renderContext);

	void Publish(Ship const & ship);

private:

	//
	// Simulation
	//

	void UpdateNpcs(
		float currentSimulationTime,
		Ship const & ship,
		LabParameters const & labParameters);

	void UpdateNpcParticle(
		StateType & npc,
		bool isPrimaryParticle,
		Ship const & ship,
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
		Ship const & ship,
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
		Ship const & ship,
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
		Ship const & ship,
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
		float currentSimulationTime,
		StateType & npc,
		bool isPrimaryParticle,
		Ship const & ship,
		LabParameters const & labParameters);

private:

	//
	// Helpers
	//

	StateType MaterializeNpcState(
		ElementIndex npcIndex,
		float currentSimulationTime,
		Ship const & ship) const;

	std::optional<StateType::NpcParticleStateType::ConstrainedStateType> CalculateParticleConstrainedState(
		vec2f const & position,
		Ship const & ship) const;

private:

	//
	// Human simulation
	//

	StateType::HumanNpcStateType InitializeHuman(
		StateType::NpcParticleStateType const & primaryParticleState,
		StateType::NpcParticleStateType const & secondaryParticleState,
		float currentSimulationTime,
		Ship const & ship) const;

	void UpdateHuman(
		float currentSimulationTime,
		StateType & npc,
		Ship const & ship,
		LabParameters const & labParameters);

	bool CheckAndMaintainHumanEquilibrium(
		ElementIndex primaryParticleIndex,
		ElementIndex secondaryParticleIndex,
		bool isRisingState,
		bool doMaintainEquilibrium,
		NpcParticles & particles,
		LabParameters const & labParameters);

	void RunWalkingHumanStateMachine(
		StateType::HumanNpcStateType & humanState,
		StateType::NpcParticleStateType const & primaryParticleState,
		Ship const & ship,
		LabParameters const & labParameters);

	void OnHumanImpact(
		vec2f const & impactVector,
		vec2f const & bounceEdgeNormal,
		StateType & npc,
		bool isPrimaryParticle) const;

	using DoImmediate = StrongTypedBool<struct _DoImmediate>;

	void FlipHumanWalk(
		StateType::HumanNpcStateType & humanState,
		DoImmediate doImmediate) const;

	void TransitionHumanToFree(
		float currentSimulationTime,
		StateType & npc);

	float CalculateActualHumanWalkingAbsoluteSpeed(
		StateType::HumanNpcStateType & humanState,
		LabParameters const & labParameters) const;

private:

	World & mParentWorld;
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

	// Cached from game parameters
	float mCurrentHumanNpcBodyLengthAdjustment;
};

}