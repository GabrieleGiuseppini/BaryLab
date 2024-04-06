/***************************************************************************************
* Original Author:      Gabriele Giuseppini
* Created:              2023-10-06
* Copyright:            Gabriele Giuseppini  (https://github.com/GabrieleGiuseppini)
***************************************************************************************/
#pragma once

#include "GameEventDispatcher.h"
#include "GameParameters.h"
#include "MaterialDatabase.h"
#include "PerfStats.h"
#include "Physics.h"
#include "RenderContext.h"

#include <GameCore/BarycentricCoords.h>
#include <GameCore/ElementIndexRangeIterator.h>
#include <GameCore/GameRandomEngine.h>
#include <GameCore/GameTypes.h>
#include <GameCore/Log.h>
#include <GameCore/StrongTypeDef.h>
#include <GameCore/Vectors.h>

#include <memory>
#include <optional>
#include <vector>

#ifdef IN_BARYLAB
#ifdef _DEBUG
#define IN_BARYLAB_DEBUG
#endif
// TODOTEST
//#define BARYLAB_PROBING
#endif

template<typename... TArgs>
void LogNpcDebug(TArgs&&... args)
{
#ifdef IN_BARYLAB_DEBUG
	LogDebug(std::forward<TArgs>(args)...);
#else
	Logger::Instance.LogToNothing(std::forward<TArgs>(args)...);
#endif
}

namespace Physics {

class Npcs final
{
private:

#pragma pack(push)

	struct LimbVector
	{
		float RightLeg;
		float LeftLeg;
		float RightArm;
		float LeftArm;

		float const * fptr() const
		{
			return &(RightLeg);
		}

		float * fptr()
		{
			return &(RightLeg);
		}
	};

#pragma pack(pop)

	//
	// The NPC state.
	//

	struct StateType final
	{
		enum class RegimeType
		{
			BeingPlaced,
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

				int CurrentVirtualEdgeOrdinal; // When set, we are "conceptually" along this edge of the current triangle - might not be really the case e.g. when we're at a vertex

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
			float MassFactor; // Purely from materials

			// Calculated
			float SpringStiffnessCoefficient;
			float SpringDampingCoefficient;

			DipolePropertiesType(
				float dipoleLength,
				float massFactor)
				: DipoleLength(dipoleLength)
				, MassFactor(massFactor)
				, SpringStiffnessCoefficient(0.0f) // Will be recalculated
				, SpringDampingCoefficient(0.0f) // Will be recalculated
			{}
		};

		struct DipoleStateType final
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

		union KindSpecificStateType
		{
			struct FurnitureNpcStateType final
			{
				FurnitureNpcKindType const Kind;

				FurnitureNpcStateType(FurnitureNpcKindType kind)
					: Kind(kind)
				{}
			} FurnitureNpcState;

			struct HumanNpcStateType final
			{
				HumanNpcKindType const Kind;
				float const Height; // This is "ideal"; dipole length is "real"
				float const WalkingSpeedBase;

				enum class BehaviorType
				{
					BeingPlaced, // Initial state, just monkeying

					Constrained_Falling, // Clueless, does nothing; with feet on edge
					Constrained_Aerial, // Clueless, does nothing; with feet in air
					Constrained_KnockedOut, // Clueless, does nothing, not relevant where; waits until can rise

					Constrained_Rising, // Tries to stand up (appliying torque)
					Constrained_Equilibrium, // Stands up; continues to adjust alignment with torque
					Constrained_Walking, // Walks; continues to adjust alignment with torque

					Free_Aerial, // Does nothing
					Free_InWater, // Does nothing, but waits to swim
					Free_Swimming_Style1, // Swims
					Free_Swimming_Style2, // Swims
					Free_Swimming_Style3  // Swims
				};

				BehaviorType CurrentBehavior;

				union BehaviorStateType
				{
					struct BeingPlacedStateType
					{
						void Reset()
						{
						}
					} BeingPlaced;

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
						float ProgressToRising;

						void Reset()
						{
							ProgressToFalling = 0.0f;
							ProgressToRising = 0.0f;
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

				// Panic levels
				float ShipOnFirePanicLevel; // [0.0f ... +1.0f]
				float GeneralizedPanicLevel; // [0.0f ... +1.0f]
				float ResultantPanicLevel; // [0.0f ... +INF)

				// Animation

				struct AnimationStateType
				{
					// Angles are CCW relative to vertical, regardless of where the NPC is looking towards (L/R)
					// (when we flip we pretend immediate mirroring of limbs from the point of view of the human,
					// so angles are independent from direction, and animation is smoother)

					// "Left" and "Right" are relative to screen when the NPC is looking at us
					// (so "right arm" is really its left arm)
					LimbVector LimbAngles;
					LimbVector LimbAnglesCos;
					LimbVector LimbAnglesSin;

					static float constexpr InitialArmAngle = Pi<float> / 2.0f * 0.3f;
					static float constexpr InitialLegAngle = 0.2f;

					LimbVector LimbLengthMultipliers;
					float LowerExtremityLengthMultiplier; // Multiplier for the part of the body from the crotch down to the feet

					AnimationStateType()
						: LimbAngles({ InitialLegAngle, -InitialLegAngle, InitialArmAngle, -InitialArmAngle })
						, LimbAnglesCos({std::cosf(LimbAngles.RightLeg), std::cosf(LimbAngles.LeftLeg), std::cosf(LimbAngles.RightArm), std::cosf(LimbAngles.LeftArm)})
						, LimbAnglesSin({ std::sinf(LimbAngles.RightLeg), std::sinf(LimbAngles.LeftLeg), std::sinf(LimbAngles.RightArm), std::sinf(LimbAngles.LeftArm) })
						, LimbLengthMultipliers({1.0f, 1.0f, 1.0f, 1.0f})
						, LowerExtremityLengthMultiplier(1.0f)
					{}
				} AnimationState;

				HumanNpcStateType(
					HumanNpcKindType kind,
					float height,
					float walkingSpeedBase,
					BehaviorType initialBehavior,
					float currentSimulationTime)
					: Kind(kind)
					, Height(height)
					, WalkingSpeedBase(walkingSpeedBase)
					, CurrentEquilibriumSoftTerminationDecision(0.0f)
					, CurrentFaceOrientation(1.0f)
					, CurrentFaceDirectionX(0.0f)
					, ShipOnFirePanicLevel(0.0f)
					, GeneralizedPanicLevel(0.0f)
					, ResultantPanicLevel(0.0f)
					// Animation
					, AnimationState()
				{
					TransitionToState(initialBehavior, currentSimulationTime);
				}

				void TransitionToState(
					BehaviorType behavior,
					float currentSimulationTime)
				{
					LogNpcDebug("  HumanBehaviorTransition: ", int(CurrentBehavior), " -> ", int(behavior));

					CurrentBehavior = behavior;

					switch (behavior)
					{
						case BehaviorType::BeingPlaced:
						{
							CurrentBehaviorState.BeingPlaced.Reset();
							break;
						}

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

						case BehaviorType::Free_Swimming_Style1:
						case BehaviorType::Free_Swimming_Style2:
						case BehaviorType::Free_Swimming_Style3:
						{
							CurrentBehaviorState.Free_Swimming.Reset();
							break;
						}
					}

					CurrentStateTransitionSimulationTimestamp = currentSimulationTime;
					TotalDistanceTraveledOnEdgeSinceStateTransition = 0.0f;
					TotalDistanceTraveledOffEdgeSinceStateTransition = 0.0f;
				}
			} HumanNpcState;

			explicit KindSpecificStateType(FurnitureNpcStateType && furnitureState)
				: FurnitureNpcState(std::move(furnitureState))
			{}

			explicit KindSpecificStateType(HumanNpcStateType && humanState)
				: HumanNpcState(std::move(humanState))
			{}
		};

		//
		// Members
		//

		// The ID of this NPC.
		NpcId const Id;

		// The type of this NPC.
		NpcKindType const Kind;

		// The current ship that this NPC belongs to.
		// NPCs always belong to a ship, and can change ships during the
		// course of their lives.
		ShipId CurrentShipId;

		// The current plane ID. When not set, it stands-in for the topmost
		// plane ID in the current ship.
		std::optional<PlaneId> CurrentPlaneId;

		// The current regime.
		RegimeType CurrentRegime;

		// The state of the first (mandatory) particle (e.g. the feet of a human).
		NpcParticleStateType PrimaryParticleState;

		// The state of the dipole, when this NPC is a dipole.
		std::optional<DipoleStateType> DipoleState;

		// The additional state specific to the type of this NPC.
		KindSpecificStateType KindSpecificState;

		// The current highlight state of this NPC.
		NpcHighlightType Highlight;

		// Randomness specific to this NPC.
		float RandomNormalizedUniformSeed;

		StateType(
			NpcId id,
			NpcKindType kind,
			ShipId initialShipId,
			std::optional<PlaneId> initialPlaneId,
			RegimeType initialRegime,
			NpcParticleStateType && primaryParticleState,
			std::optional<DipoleStateType> && dipoleState,
			KindSpecificStateType && kindSpecificState)
			: Id(id)
			, Kind(kind)
			, CurrentShipId(initialShipId)
			, CurrentPlaneId(std::move(initialPlaneId))
			, CurrentRegime(initialRegime)
			, PrimaryParticleState(std::move(primaryParticleState))
			, DipoleState(std::move(dipoleState))
			, KindSpecificState(std::move(kindSpecificState))
			, Highlight(NpcHighlightType::None)
			, RandomNormalizedUniformSeed(GameRandomEngine::GetInstance().GenerateNormalizedUniformReal())
		{}
	};

	//
	// The information heading the list of NPCs in a ship.
	//

	struct ShipNpcsType final
	{
		Ship const & ShipMesh;
		std::vector<NpcId> Npcs;

		ShipNpcsType(Ship const & shipMesh)
			: ShipMesh(shipMesh)
			, Npcs()
		{}
	};

public:

	Npcs(
		Physics::World & parentWorld,
		MaterialDatabase const & materialDatabase,
		std::shared_ptr<GameEventDispatcher> gameEventHandler,
		GameParameters const & gameParameters)
		: mParentWorld(parentWorld)
		, mMaterialDatabase(materialDatabase)
		, mGameEventHandler(std::move(gameEventHandler))
		// Container
		, mStateBuffer()
		, mShips()
		, mParticles(GameParameters::MaxNpcs * GameParameters::MaxParticlesPerNpc)
		// Stats
		, mHumanNpcCount(0)
		, mFurnitureNpcCount(0)
		, mFreeRegimeHumanNpcCount(0)
		, mConstrainedRegimeHumanNpcCount(0)
		// Simulation Parameters
		, mCurrentHumanNpcBodyLengthAdjustment(gameParameters.HumanNpcBodyLengthAdjustment)
		, mCurrentSpringReductionFraction(gameParameters.NpcSpringReductionFraction)
		, mCurrentSpringDampingCoefficient(gameParameters.NpcSpringReductionFraction)
	{}

	void Update(
		float currentSimulationTime,
		GameParameters const & gameParameters);

	void RenderUpload(
		RenderContext & renderContext,
		PerfStats & perfStats);

	///////////////////////////////

	NpcParticles const & GetParticles() const
	{
		return mParticles;
	}

	void OnShipAdded(Ship const & ship);

	void OnShipRemoved(ShipId shipId);

	std::optional<PickedObjectId<NpcId>> BeginPlaceNewFurnitureNpc(
		FurnitureNpcKindType furnitureKind,
		vec2f const & worldCoordinates,
		float currentSimulationTime);

	std::optional<PickedObjectId<NpcId>> BeginPlaceNewHumanNpc(
		HumanNpcKindType humanKind,
		vec2f const & worldCoordinates,
		float currentSimulationTime);

	std::optional<PickedObjectId<NpcId>> ProbeNpcAt(
		vec2f const & position,
		GameParameters const & gameParameters) const;

	void BeginMoveNpc(
		NpcId id,
		float currentSimulationTime);

	void MoveNpcTo(
		NpcId id,
		vec2f const & position,
		vec2f const & offset);

	void EndMoveNpc(
		NpcId id,
		float currentSimulationTime);

	void CompleteNewNpc(
		NpcId id,
		float currentSimulationTime);

	void RemoveNpc(NpcId id);

	void AbortNewNpc(NpcId id);

	void HighlightNpc(
		NpcId id,
		NpcHighlightType highlight);

	void SetGeneralizedPanicLevelForAllHumans(float panicLevel);

public:

#ifdef IN_BARYLAB

	//
	// Barylab-specific
	//

	bool AddHumanNpc(
		HumanNpcKindType humanKind,
		vec2f const & worldCoordinates,
		float currentSimulationTime);

	void FlipHumanWalk(int npcIndex);

	void FlipHumanFrontBack(int npcIndex);

	void MoveParticleBy(
		ElementIndex particleIndex,
		vec2f const & offset,
		float currentSimulationTime);

	void RotateParticlesWithShip(
		vec2f const & centerPos,
		float cosAngle,
		float sinAngle);

	void RotateParticleWithShip(
		StateType::NpcParticleStateType const & npcParticleState,
		vec2f const & centerPos,
		float cosAngle,
		float sinAngle,
		Ship const & shipMesh);

	void OnPointMoved(float currentSimulationTime);

	void OnMassAdjustmentChanged(float massAdjustment)
	{
		if (massAdjustment != mCurrentMassAdjustment)
		{
			mCurrentMassAdjustment = massAdjustment;
			RecalculateSpringForceParameters();
		}
	}

	void OnGravityAdjustmentChanged(float gravityAdjustment)
	{
		if (gravityAdjustment != mCurrentGravityAdjustment)
		{
			mCurrentGravityAdjustment = gravityAdjustment;
			RecalculateSpringForceParameters();
		}
	}

	//
	// Probing
	//

	void SelectNpc(NpcId npcId)
	{
		mCurrentlySelectedNpc = npcId;
		Publish();
	}

	std::optional<NpcId> GetCurrentlySelectedNpc() const
	{
		return mCurrentlySelectedNpc;
	}

	void DeselectNpc()
	{
		mCurrentlySelectedNpc.reset();
		Publish();
	}

	void SelectParticle(ElementIndex particleIndex)
	{
		mCurrentlySelectedParticle = particleIndex;
		Publish();
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

	bool IsSpringHostingCurrentlySelectedParticle(ElementIndex springIndex) const;

	void Publish() const;

#endif

private:

	NpcId GetNewNpcId();

	ShipId GetTopmostShipId() const;

	std::optional<ElementId> FindTopmostTriangleContaining(vec2f const & position) const;

	static ElementIndex FindTriangleContaining(
		vec2f const & position,
		Ship const & shipMesh);

	void TransferNpcToShip(
		StateType & npc,
		ShipId newShip);

	void RenderNpc(
		StateType const & npc,
		ShipRenderContext & shipRenderContext);

	void PublishHumanNpcStats();

private:

	//
	// Simulation
	//

	void ResetNpcStateToWorld(
		StateType & npc,
		float currentSimulationTime);

	void ResetNpcStateToWorld(
		StateType & npc,
		float currentSimulationTime,
		Ship const & shipMesh,
		std::optional<ElementIndex> primaryParticleTriangleIndex) const;

	void TransitionParticleToConstrainedState(
		StateType & npc,
		bool isPrimaryParticle,
		StateType::NpcParticleStateType::ConstrainedStateType constrainedState);

	void TransitionParticleToFreeState(
		StateType & npc,
		bool isPrimaryParticle);

	static std::optional<StateType::NpcParticleStateType::ConstrainedStateType> CalculateParticleConstrainedState(
		vec2f const & position,
		Ship const & shipMesh,
		std::optional<ElementIndex> triangleIndex);

	void OnMayBeNpcRegimeChanged(
		StateType::RegimeType oldRegime,
		StateType & npc);

	static StateType::RegimeType CalculateRegime(StateType const & npc);

	void UpdateNpcs(
		float currentSimulationTime,
		GameParameters const & gameParameters);

	void UpdateNpcParticlePhysics(
		StateType & npc,
		bool isPrimaryParticle,
		Ship const & shipMesh,
		float currentSimulationTime,
		GameParameters const & gameParameters);

	void CalculateNpcParticlePreliminaryForces(
		StateType const & npc,
		bool isPrimaryParticle,
		GameParameters const & gameParameters);

	vec2f CalculateNpcParticleDefinitiveForces(
		StateType const & npc,
		bool isPrimaryParticle,
		float particleMass,
		GameParameters const & gameParameters) const;

	void RecalculateSpringForceParameters();

	void RecalculateSpringForceParameters(StateType::DipolePropertiesType & dipoleProperties) const;

	void RecalculateHumanNpcBodyLengths();

	void UpdateNpcParticle_Free(
		StateType::NpcParticleStateType & particle,
		vec2f const & startPosition,
		vec2f const & endPosition,
		NpcParticles & particles,
		GameParameters const & gameParameters) const;

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
		Ship const & shipMesh,
		NpcParticles & particles,
		float currentSimulationTime,
		GameParameters const & gameParameters);

	float UpdateNpcParticle_ConstrainedInertial(
		StateType & npc,
		bool isPrimaryParticle,
		vec2f const & particleStartAbsolutePosition,
		bcoords3f const segmentTrajectoryStartBarycentricCoords,
		vec2f const & segmentTrajectoryEndAbsolutePosition,
		bcoords3f segmentTrajectoryEndBarycentricCoords,
		vec2f const meshVelocity,
		float segmentDt,
		Ship const & shipMesh,
		NpcParticles & particles,
		float currentSimulationTime,
		GameParameters const & gameParameters);

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
		Ship const & shipMesh,
		NpcParticles & particles,
		GameParameters const & gameParameters);

	inline void BounceConstrainedNpcParticle(
		StateType & npc,
		bool isPrimaryParticle,
		vec2f const & trajectory,
		vec2f const & bouncePosition,
		vec2f const & bounceEdgeNormal,
		vec2f const meshVelocity,
		float dt,
		NpcParticles & particles,
		float currentSimulationTime,
		GameParameters const & gameParameters) const;

	void OnImpact(
		StateType & npc,
		bool isPrimaryParticle,
		vec2f const & normalResponse,
		vec2f const & bounceEdgeNormal,
		float currentSimulationTime) const;

	void UpdateNpcAnimation(
		StateType & npc,
		bool isPrimaryParticle,
		float currentSimulationTime,
		Ship const & shipMesh,
		GameParameters const & gameParameters);

	static bool IsEdgeFloorToParticle(
		ElementIndex triangleElementIndex,
		int edgeOrdinal,
		StateType const & npc,
		bool isPrimaryParticle,
		NpcParticles const & npcParticles,
		Ship const & shipMesh)
	{
		// First off: if not a floor, it's not a floor

		if (shipMesh.GetTriangles().GetSubSpringNpcSurfaceType(triangleElementIndex, edgeOrdinal) == NpcSurfaceType::Open)
		{
			return false;
		}

		// Ok, it's a floor

		// If it's a primary, then every floor is a floor

		if (isPrimaryParticle)
		{
			return true;
		}

		// Ok, it's a floor and this is a secondary particle

		assert(npc.DipoleState.has_value());

		// If the primary is not on the other side of this edge, then every floor is a floor

		vec2f const & primaryPosition = npcParticles.GetPosition(npc.PrimaryParticleState.ParticleIndex);
		bcoords3f const primaryBaryCoords = shipMesh.GetTriangles().ToBarycentricCoordinates(primaryPosition, triangleElementIndex, shipMesh.GetPoints());

		// It's on the other side of the edge if its "edge's" b-coord is negative
		if (primaryBaryCoords[(edgeOrdinal + 2) % 3] >= 0.0f)
		{
			return true;
		}

		// Ok, it's a floor and it's separating this secondary particle from the primary

		// Now a bit of a hack: at this moment we're hurting because of the "hanging head" problem, i.e.
		// a human NPC ending with its head on an edge and it feet hanging underneath. To prevent this,
		// we consider this as a floor only if the human is not "quite vertical"

		vec2f const & secondaryPosition = npcParticles.GetPosition(npc.DipoleState->SecondaryParticleState.ParticleIndex);
		vec2f const humanDir = (primaryPosition - secondaryPosition).normalise(); // Pointing down to feet

		// It's vertical when y is -1.0 (cos of angle)
		return humanDir.y > -0.8f;
	}

	static bool IsOnFloorEdge(
		Npcs::StateType::NpcParticleStateType::ConstrainedStateType const & constrainedState,
		Ship const & shipMesh)
	{
		auto const & baryCoords = constrainedState.CurrentTriangleBarycentricCoords;
		auto const & triangleIndex = constrainedState.CurrentTriangle;

		if (baryCoords[0] == 0.0f
			&& shipMesh.GetTriangles().GetSubSpringNpcSurfaceType(triangleIndex, 1) != NpcSurfaceType::Open)
		{
			return true;
		}

		if (baryCoords[1] == 0.0f
			&& shipMesh.GetTriangles().GetSubSpringNpcSurfaceType(triangleIndex, 2) != NpcSurfaceType::Open)
		{
			return true;
		}

		if (baryCoords[2] == 0.0f
			&& shipMesh.GetTriangles().GetSubSpringNpcSurfaceType(triangleIndex, 0) != NpcSurfaceType::Open)
		{
			return true;
		}

		return false;
	}

	static vec2f CalculateDipoleVector(ElementIndex primaryParticleIndex, ElementIndex secondaryParticleIndex, NpcParticles const & particles)
	{
		return particles.GetPosition(primaryParticleIndex) - particles.GetPosition(secondaryParticleIndex);
	}

	static float CalculateVerticalAlignment(vec2f const & vector)
	{
		return vector.normalise().dot(GameParameters::GravityDir);
	}

	static float CalculateDipoleVerticalAlignment(ElementIndex primaryParticleIndex, ElementIndex secondaryParticleIndex, NpcParticles const & particles)
	{
		return CalculateVerticalAlignment(CalculateDipoleVector(primaryParticleIndex, secondaryParticleIndex, particles));
	}

	//
	// Human simulation
	//

	static StateType::KindSpecificStateType::HumanNpcStateType::BehaviorType CalculateHumanBehavior(StateType & npc);

	void UpdateHuman(
		StateType & npc,
		float currentSimulationTime,
		Ship const & shipMesh,
		GameParameters const & gameParameters);

	bool CheckAndMaintainHumanEquilibrium(
		ElementIndex primaryParticleIndex,
		ElementIndex secondaryParticleIndex,
		bool isRisingState,
		bool doMaintainEquilibrium,
		NpcParticles & particles,
		GameParameters const & gameParameters);

	void RunWalkingHumanStateMachine(
		StateType::KindSpecificStateType::HumanNpcStateType & humanState,
		StateType::NpcParticleStateType const & primaryParticleState,
		Ship const & shipMesh,
		GameParameters const & gameParameters);

	void OnHumanImpact(
		StateType & npc,
		bool isPrimaryParticle,
		vec2f const & normalResponse,
		vec2f const & bounceEdgeNormal,
		float currentSimulationTime) const;

	using DoImmediate = StrongTypedBool<struct _DoImmediate>;

	void FlipHumanWalk(
		StateType::KindSpecificStateType::HumanNpcStateType & humanState,
		DoImmediate doImmediate) const;

	void TransitionHumanBehaviorToFree(
		StateType & npc,
		float currentSimulationTime);

	float CalculateActualHumanWalkingAbsoluteSpeed(
		StateType::KindSpecificStateType::HumanNpcStateType & humanState,
		GameParameters const & gameParameters) const;

	float CalculateHumanWalkingSpeedAdjustment(
		StateType::KindSpecificStateType::HumanNpcStateType & humanState,
		GameParameters const & gameParameters) const;

private:

	World & mParentWorld;
	MaterialDatabase const & mMaterialDatabase;
	std::shared_ptr<GameEventDispatcher> mGameEventHandler;

	//
	// Container
	//
	// Use cases:
	//	1. Reaching all NPCs of a specific ship(e.g.because of ship - wide interactions, such as electrical tool, alarm, deleting ship, etc.)
	//	2. Allow an NPC to move ships(e.g.the one "being placed")
	//	3. Reaching an NPC by its ID
	//

	// The actual container of NPC states, indexed by NPC ID.
	// Indices are stable; elements are null'ed when removed.
	std::vector<std::optional<StateType>> mStateBuffer;

	// All the ships - together with their NPCs - indexed by Ship ID.
	// Indices are stable; elements are null'ed when removed.
	std::vector<std::optional<ShipNpcsType>> mShips;

	// All of the NPC particles.
	NpcParticles mParticles;

	//
	// Stats
	//

	ElementCount mHumanNpcCount;
	ElementCount mFurnitureNpcCount;
	ElementCount mFreeRegimeHumanNpcCount;
	ElementCount mConstrainedRegimeHumanNpcCount;

	//
	// Simulation parameters
	//

	// Cached from game parameters
	float mCurrentHumanNpcBodyLengthAdjustment;
	float mCurrentSpringReductionFraction;
	float mCurrentSpringDampingCoefficient;

#ifdef IN_BARYLAB

	// Cached from LabController
	float mCurrentMassAdjustment{ 1.0f };
	float mCurrentGravityAdjustment{ 1.0f };

#endif

#ifdef IN_BARYLAB

	//
	// Probing
	//

	std::optional<NpcId> mCurrentlySelectedNpc;
	std::optional<ElementIndex> mCurrentlySelectedParticle;
	std::optional<ElementIndex> mCurrentOriginTriangle;

	std::optional<ParticleTrajectory> mCurrentParticleTrajectory;
	std::optional<ParticleTrajectory> mCurrentParticleTrajectoryNotification;

#endif
};

}