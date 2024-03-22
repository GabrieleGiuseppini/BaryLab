/***************************************************************************************
 * Original Author:		Gabriele Giuseppini
 * Created:				2023-07-23
 * Copyright:			Gabriele Giuseppini  (https://github.com/GabrieleGiuseppini)
 ***************************************************************************************/
#include "Physics.h"

#include <GameCore/GameMath.h>

#include <algorithm>
#include <limits>

namespace Physics {

void Npcs::ResetNpcStateToWorld(
    StateType & npc,
    float currentSimulationTime)
{
    //
    // Find topmost triangle - among all ships - which contains this NPC
    //

    // Take the position of the primary particle as representative of the NPC
    auto const & position = mParticles.GetPosition(npc.PrimaryParticleState.ParticleIndex);

    auto const topmostTriangle = FindTopmostTriangleContaining(position);
    if (topmostTriangle.has_value())
    {
        // Primary is in a triangle!

        assert(mShips[topmostTriangle->GetShipId()].has_value());

        TransferNpcToShip(
            npc,
            topmostTriangle->GetShipId());

        ResetNpcStateToWorld(
            npc,
            currentSimulationTime,
            mShips[topmostTriangle->GetShipId()]->ShipMesh,
            topmostTriangle->GetLocalObjectId());
    }
    else
    {
        // No luck; means we're free, and pick topmost ship for that

        auto const topmostShipId = GetTopmostShipId();

        assert(mShips[topmostShipId].has_value());

        TransferNpcToShip(
            npc,
            topmostShipId);

        ResetNpcStateToWorld(
            npc,
            currentSimulationTime,
            mShips[GetTopmostShipId()]->ShipMesh,
            std::nullopt);
    }
}

void Npcs::ResetNpcStateToWorld(
    StateType & npc,
    float currentSimulationTime,
    Ship const & shipMesh,
    std::optional<ElementIndex> primaryParticleTriangleIndex) const
{
    // Plane ID

    if (primaryParticleTriangleIndex.has_value())
    {
        // Use the plane ID of this triangle
        ElementIndex const trianglePointIndex = shipMesh.GetTriangles().GetPointAIndex(*primaryParticleTriangleIndex);
        npc.CurrentPlaneId = shipMesh.GetPoints().GetPlaneId(trianglePointIndex);
    }
    else
    {
        // Primary is free, hence this NPC is on the topmost plane ID of its current ship
        npc.CurrentPlaneId = std::nullopt;
    }

    // Primary particle

    npc.PrimaryParticleState.ConstrainedState = CalculateParticleConstrainedState(
        mParticles.GetPosition(npc.PrimaryParticleState.ParticleIndex),
        shipMesh,
        primaryParticleTriangleIndex);

    // Secondary particle

    if (npc.DipoleState.has_value())
    {
        npc.DipoleState->SecondaryParticleState.ConstrainedState = CalculateParticleConstrainedState(
            mParticles.GetPosition(npc.DipoleState->SecondaryParticleState.ParticleIndex),
            shipMesh,
            std::nullopt);
    }

    // Regime

    npc.CurrentRegime = CalculateRegime(npc);

    // Kind specific

    switch (npc.Kind)
    {
        case NpcKindType::Furniture:
        {
            // Reset state
            npc.KindSpecificState.FurnitureNpcState = StateType::KindSpecificStateType::FurnitureNpcStateType();

            break;
        }

        case NpcKindType::Human:
        {
            // Change behavior
            npc.KindSpecificState.HumanNpcState.TransitionToState(
                CalculateHumanBehavior(npc),
                currentSimulationTime);

            break;
        }
    }

    Publish();
}

void Npcs::TransitionParticleToConstrainedState(
    StateType & npc,
    bool isPrimaryParticle,
    StateType::NpcParticleStateType::ConstrainedStateType constrainedState)
{
    // Transition

    if (isPrimaryParticle)
    {
        npc.PrimaryParticleState.ConstrainedState = constrainedState;
    }
    else
    {
        assert(npc.DipoleState.has_value());
        npc.DipoleState->SecondaryParticleState.ConstrainedState = constrainedState;
    }

    // Regime

    auto const oldRegime = npc.CurrentRegime;

    npc.CurrentRegime = CalculateRegime(npc);

    OnMayBeNpcRegimeChanged(oldRegime, npc);
}

void Npcs::TransitionParticleToFreeState(
    StateType & npc,
    bool isPrimaryParticle)
{
    // Transition

    if (isPrimaryParticle)
    {
        npc.PrimaryParticleState.ConstrainedState.reset();
    }
    else
    {
        assert(npc.DipoleState.has_value());
        npc.DipoleState->SecondaryParticleState.ConstrainedState.reset();
    }

    // Regime

    auto const oldRegime = npc.CurrentRegime;

    npc.CurrentRegime = CalculateRegime(npc);

    OnMayBeNpcRegimeChanged(oldRegime, npc);
}

std::optional<Npcs::StateType::NpcParticleStateType::ConstrainedStateType> Npcs::CalculateParticleConstrainedState(
    vec2f const & position,
    Ship const & shipMesh,
    std::optional<ElementIndex> triangleIndex)
{
    std::optional<StateType::NpcParticleStateType::ConstrainedStateType> constrainedState;

    if (!triangleIndex.has_value())
    {
        triangleIndex = FindTriangleContaining(position, shipMesh);
    }

    assert(triangleIndex.has_value());

    if (triangleIndex != NoneElementIndex && !shipMesh.GetTriangles().IsDeleted(*triangleIndex))
    {
        bcoords3f const barycentricCoords = shipMesh.GetTriangles().ToBarycentricCoordinatesFromWithinTriangle(
            position,
            *triangleIndex,
            shipMesh.GetPoints());

        assert(barycentricCoords.is_on_edge_or_internal());

        return Npcs::StateType::NpcParticleStateType::ConstrainedStateType(
            *triangleIndex,
            barycentricCoords);
    }

    return std::nullopt;
}

void Npcs::OnMayBeNpcRegimeChanged(
    StateType::RegimeType oldRegime,
    StateType & npc)
{
    if (oldRegime == npc.CurrentRegime)
    {
        // Nothing to do
        return;
    }

    if (npc.Kind == NpcKindType::Human)
    {
        //
        // Update stats
        //

        bool doPublishStats = false;

        if (oldRegime == StateType::RegimeType::Constrained)
        {
            assert(mConstrainedRegimeHumanNpcCount > 0);
            --mConstrainedRegimeHumanNpcCount;
            doPublishStats = true;
        }
        else if (oldRegime == StateType::RegimeType::Free)
        {
            assert(mFreeRegimeHumanNpcCount > 0);
            --mFreeRegimeHumanNpcCount;
            doPublishStats = true;
        }

        if (npc.CurrentRegime == StateType::RegimeType::Constrained)
        {
            ++mConstrainedRegimeHumanNpcCount;
            doPublishStats = true;
        }
        else if (npc.CurrentRegime == StateType::RegimeType::Free)
        {
            ++mFreeRegimeHumanNpcCount;
            doPublishStats = true;
        }

        if (doPublishStats)
        {
            PublishHumanNpcStats();
        }
    }
}

Npcs::StateType::RegimeType Npcs::CalculateRegime(StateType const & npc)
{
    // Constrained if at least one is constrained;
    // conversely: Free only if both are free
    return npc.PrimaryParticleState.ConstrainedState.has_value()
        ? StateType::RegimeType::Constrained
        : (npc.DipoleState.has_value() && npc.DipoleState->SecondaryParticleState.ConstrainedState.has_value()) ? StateType::RegimeType::Constrained : StateType::RegimeType::Free;
}

void Npcs::UpdateNpcs(
    float currentSimulationTime,
    GameParameters const & gameParameters)
{
    LogNpcDebug("----------------------------------");
    LogNpcDebug("----------------------------------");
    LogNpcDebug("----------------------------------");

    //
    // 1. Reset buffers
    //

    // Note: no need to reset PreliminaryForces as we'll recalculate all of them

    //
    // 2. Check if a free secondary particle should become constrained
    // 3. Calculate preliminary forces
    //

    for (auto & npcState : mStateBuffer)
    {
        if (npcState.has_value())
        {
            // Secondary free becoming constrained

            if (npcState->DipoleState.has_value()
                && !npcState->DipoleState->SecondaryParticleState.ConstrainedState.has_value() // Secondary is free
                && npcState->PrimaryParticleState.ConstrainedState.has_value()) // And primary is constrained
            {
                assert(mShips[npcState->CurrentShipId].has_value());

                auto newConstrainedState = CalculateParticleConstrainedState(
                    mParticles.GetPosition(npcState->DipoleState->SecondaryParticleState.ParticleIndex),
                    mShips[npcState->CurrentShipId]->ShipMesh,
                    std::nullopt);

                if (newConstrainedState.has_value())
                {
                    TransitionParticleToConstrainedState(*npcState, false, std::move(*newConstrainedState));
                }
            }

            // Preliminary Forces

            CalculateNpcParticlePreliminaryForces(
                *npcState,
                true,
                gameParameters);

            if (npcState->DipoleState.has_value())
            {
                CalculateNpcParticlePreliminaryForces(
                    *npcState,
                    false,
                    gameParameters);
            }
        }
    }

    //
    // 4. Update physical state
    //

    for (auto & npcState : mStateBuffer)
    {
        if (npcState.has_value())
        {
            LogNpcDebug("NPC ", npcState->Id);

            assert(mShips[npcState->CurrentShipId].has_value());
            auto const & shipMesh = mShips[npcState->CurrentShipId]->ShipMesh;

            UpdateNpcParticlePhysics(
                *npcState,
                true,
                shipMesh,
                gameParameters);

            if (npcState->DipoleState.has_value())
            {
                UpdateNpcParticlePhysics(
                    *npcState,
                    false,
                    shipMesh,
                    gameParameters);
            }
        }
    }

    //
    // 5. Update behavioral state machines
    //

    mParticles.ResetEquilibriumTorque();

    for (auto & npcState : mStateBuffer)
    {
        if (npcState.has_value())
        {
            if (npcState->Kind == NpcKindType::Human)
            {
                assert(mShips[npcState->CurrentShipId].has_value());
                auto const & shipMesh = mShips[npcState->CurrentShipId]->ShipMesh;

                UpdateHuman(
                    *npcState,
                    currentSimulationTime,
                    shipMesh,
                    gameParameters);
            }
        }
    }

    //
    // 6. Update animation
    //

    for (auto & npcState : mStateBuffer)
    {
        if (npcState.has_value())
        {
            LogNpcDebug("NPC ", npcState->Id);

            assert(mShips[npcState->CurrentShipId].has_value());
            auto const & shipMesh = mShips[npcState->CurrentShipId]->ShipMesh;

            UpdateNpcAnimation(
                *npcState,
                true,
                currentSimulationTime,
                shipMesh,
                gameParameters);

            if (npcState->DipoleState.has_value())
            {
                UpdateNpcAnimation(
                    *npcState,
                    false,
                    currentSimulationTime,
                    shipMesh,
                    gameParameters);
            }
        }
    }
}

void Npcs::UpdateNpcParticlePhysics(
    StateType & npc,
    bool isPrimaryParticle,
    Ship const & shipMesh,
    GameParameters const & gameParameters)
{
    //
    // Here be dragons!
    //

    auto & npcParticle = isPrimaryParticle ? npc.PrimaryParticleState : npc.DipoleState->SecondaryParticleState;

    LogNpcDebug("----------------------------------");
    LogNpcDebug("  ", isPrimaryParticle ? "Primary" : "Secondary");

    float const particleMass = mParticles.GetPhysicalProperties(npcParticle.ParticleIndex).Mass * gameParameters.MassAdjustment;
    float const dt = GameParameters::SimulationTimeStepDuration;

    vec2f const particleStartAbsolutePosition = mParticles.GetPosition(npcParticle.ParticleIndex);

    // Calculate physical displacement - once and for all, as whole loop
    // will attempt to move to trajectory that always ends here

    vec2f physicsDeltaPos;
#ifdef IN_BARYLAB
    if (mCurrentParticleTrajectory.has_value() && npcParticle.ParticleIndex == mCurrentParticleTrajectory->ParticleIndex)
    {
        // Consume externally-supplied trajectory

        physicsDeltaPos = mCurrentParticleTrajectory->TargetPosition - particleStartAbsolutePosition;
        mCurrentParticleTrajectory.reset();
    }
    else
#endif
    {
        // Integrate forces

        vec2f const physicalForces = CalculateNpcParticleDefinitiveForces(
            npc,
            isPrimaryParticle,
            particleMass,
            gameParameters);

        physicsDeltaPos = mParticles.GetVelocity(npcParticle.ParticleIndex) * dt + (physicalForces / particleMass) * dt * dt;
    }

    if (npc.CurrentRegime == StateType::RegimeType::BeingPlaced
        && ((npc.Kind == NpcKindType::Human && !isPrimaryParticle) || (npc.Kind != NpcKindType::Human)))
    {
        //
        // Particle is being placed
        //

        LogNpcDebug("    Being placed");
    }
    else if (!npcParticle.ConstrainedState.has_value())
    {
        //
        // Particle is free
        //

        LogNpcDebug("    Free: velocity=", mParticles.GetVelocity(npcParticle.ParticleIndex), " prelimF=", mParticles.GetPreliminaryForces(npcParticle.ParticleIndex), " physicsDeltaPos=", physicsDeltaPos);
        LogNpcDebug("    StartPosition=", particleStartAbsolutePosition, " StartVelocity=", mParticles.GetVelocity(npcParticle.ParticleIndex));

        UpdateNpcParticle_Free(
            npcParticle,
            particleStartAbsolutePosition,
            particleStartAbsolutePosition + physicsDeltaPos,
            mParticles,
            gameParameters);

        LogNpcDebug("    EndPosition=", mParticles.GetPosition(npcParticle.ParticleIndex), " EndVelocity=", mParticles.GetVelocity(npcParticle.ParticleIndex));

        // Update total distance traveled
        if (npc.Kind == NpcKindType::Human
            && isPrimaryParticle) // Human is represented by primary particle
        {
            npc.KindSpecificState.HumanNpcState.TotalDistanceTraveledOffEdgeSinceStateTransition += physicsDeltaPos.length();
        }

        // We're done
    }
    else
    {
        //
        // Constrained
        //

        assert(npcParticle.ConstrainedState.has_value());

        LogNpcDebug("    Constrained: velocity=", mParticles.GetVelocity(npcParticle.ParticleIndex), " prelimF=", mParticles.GetPreliminaryForces(npcParticle.ParticleIndex), " physicsDeltaPos=", physicsDeltaPos);

        // Loop tracing trajectory from TrajectoryStart (== current bary coords) to TrajectoryEnd (== start absolute pos + deltaPos);
        // each step moves the next TrajectoryStart a bit ahead.
        // Each iteration of the loop either exits (completes), or moves current bary coords (and calcs remaining dt) when it wants
        // to "continue" an impact while on edge-moving-against-it, i.e. when it wants to recalculate a new flattened traj.
        //    - In this case, the iteration doesn't change current absolute position nor velocity; it only updates current bary coords
        //      to what will become the next TrajectoryStart
        //    - In this case, at next iteration:
        //          - TrajectoryStart (== current bary coords) is new
        //          - TrajectoryEnd (== start absolute pos + physicsDeltaPos) is same as before
        //
        // Each iteration of the loop performs either an "inertial ray-tracing step" - i.e. with the particle free to move
        // inside of a triangle - or a "non-inertial step" - i.e. with the particle pushed against an edge, and thus
        // moving in a non-inertial frame as the mesh acceleration spawns the appearance of apparent forces.
        // - After an inertial iteration we won't enter a non-inertial iteration
        // - A non-inertial iteration might be followed by an inertial one

        // Initialize trajectory start (wrt current triangle) as the absolute pos of the particle as if it just
        // moved with the mesh, staying in its position wrt its triangle; in other words, it's the new theoretical
        // position after just mesh displacement
        vec2f trajectoryStartAbsolutePosition = shipMesh.GetTriangles().FromBarycentricCoordinates(
            npcParticle.ConstrainedState->CurrentTriangleBarycentricCoords,
            npcParticle.ConstrainedState->CurrentTriangle,
            shipMesh.GetPoints());

        // Calculate mesh velocity for the whole loop as the pure displacement of the triangle containing this particle
        vec2f const meshVelocity = (particleStartAbsolutePosition - trajectoryStartAbsolutePosition) / GameParameters::SimulationTimeStepDuration;

        // Machinery to detect 2- or 3-iteration paths that don't move particle (positional well,
        // aka gravity well)

        struct PastBarycentricPosition
        {
            ElementIndex TriangleElementIndex;
            bcoords3f BarycentricCoords;

            PastBarycentricPosition(
                ElementIndex triangleElementIndex,
                bcoords3f barycentricCoords)
                : TriangleElementIndex(triangleElementIndex)
                , BarycentricCoords(barycentricCoords)
            {}

            bool operator==(PastBarycentricPosition const & other) const
            {
                return this->TriangleElementIndex == other.TriangleElementIndex
                    && this->BarycentricCoords == other.BarycentricCoords;
            }
        };

        std::array<std::optional<PastBarycentricPosition>, 2> pastBarycentricPositions = { // A queue, insert at zero - initialized with now
            PastBarycentricPosition(npcParticle.ConstrainedState->CurrentTriangle, npcParticle.ConstrainedState->CurrentTriangleBarycentricCoords),
            std::nullopt
        };

        // Total displacement walked along the edge - as a sum of the vectors from the individual steps
        //
        // We will consider the particle's starting position to incrementally move by this,
        // in order to keep trajectory directory invariant with walk
        vec2f totalEdgeWalkedActual = vec2f::zero();

        // And here's for something seemingly obscure.
        //
        // At our first trajectory flattening for a non-inertial step we take the calculated
        // absolute flattened trajectory length (including walk) as the maximum (absolute)
        // distance that we're willing to travel in the (remaining) time quantum.
        // After all this is really the projection of the real and apparent forces acting on
        // the particle onto the (first) edge, and the one we should keep as the particle's
        // movement vector is now tied to the edge.
        //
        // At times an iteration might want to travel more than what we had decided is the
        // max we're willing to, for example because we've traveled a lot almost orthogonally
        // to the theoretical trajectory, and while there's little left along the flattened
        // trajectory, the trajectory end point might still be quite far from where we are.
        // If we stop being constrained by the edge that causes the travel to be almost
        // orthogonal to the trajectory, we might become subject to an abnormal quantity of
        // displacement - yielding also an abnormal velocity calculation.
        //
        // We thus clamp the magnitude of flattened trajectory vectors so to never exceed
        // this max distance we're willing to travel - hence the nickname "budget".

        std::optional<float> edgeDistanceToTravelMax;
        float edgeDistanceTraveledTotal = 0.0f; // To keep track of total distance

        for (float remainingDt = dt; ; )
        {
            assert(remainingDt > 0.0f);

            LogNpcDebug("    ------------------------");
            LogNpcDebug("    Triangle=", npcParticle.ConstrainedState->CurrentTriangle, " B-Coords=", npcParticle.ConstrainedState->CurrentTriangleBarycentricCoords, " RemainingDt=", remainingDt);

            //
            // We ray-trace the particle along a trajectory that starts at the position at which the particle
            // would be if it moved with the mesh (i.e. staying at its position wrt its triangle) and ends
            // at the position at which the particle would be if it were only subject to physical forces
            // and to walking
            //

            //
            // Calculate trajectory
            //
            // Trajectory is the apparent displacement, i.e. the displacement of the particle
            // *from the point of view of the mesh* (or of the particle itself), thus made up of
            // both the physical displacement and the mesh displacement.
            //
            // Given that trajectory discounts physics move, trajectory is the displacement
            // caused by the apparent forces. In fact, we've verified that when the particle has
            // the same velocity as the mesh, trajectory is zero.
            //
            // Note: when there's no physics displacement, this amounts to pure mesh move
            // (PartPos - FromBary(PartPosBary)). On the other hand, when the mesh is at rest, this
            // amounts to pure physical displacement.
            //

            // Absolute position of particle if it only moved due to physical forces and walking
            vec2f const trajectoryEndAbsolutePosition = particleStartAbsolutePosition + physicsDeltaPos + totalEdgeWalkedActual;

            // Trajectory
            // Note: on first iteration this is the same as physicsDeltaPos + meshDisplacement
            // TrajectoryStartAbsolutePosition changes at each iteration to be the absolute translation of the current particle bary coords
            vec2f const trajectory = trajectoryEndAbsolutePosition - trajectoryStartAbsolutePosition;

            LogNpcDebug("    TrajectoryStartAbsolutePosition=", trajectoryStartAbsolutePosition, " PhysicsDeltaPos=", physicsDeltaPos, " TotalEdgeWalkedActual=", totalEdgeWalkedActual,
                " => TrajectoryEndAbsolutePosition=", trajectoryEndAbsolutePosition, " Trajectory=", trajectory);

            //
            // Check if we are at a corner; if we are, travel through corners - according to trajectory - until
            // any of the following:
            //  - Trajectory points towards a floor
            //  - Trajectory points *inside* triangle
            //  - Particle becomes free
            //

            std::optional<int> vertexOrdinal = npcParticle.ConstrainedState->CurrentTriangleBarycentricCoords.try_get_vertex();
            if (vertexOrdinal.has_value())
            {
                bcoords3f const trajectoryEndBarycentricCoords = shipMesh.GetTriangles().ToBarycentricCoordinates(
                    trajectoryEndAbsolutePosition,
                    npcParticle.ConstrainedState->CurrentTriangle,
                    shipMesh.GetPoints());

                auto const outcome = NavigateVertex(
                    npc,
                    isPrimaryParticle,
                    *vertexOrdinal,
                    particleStartAbsolutePosition,
                    trajectoryStartAbsolutePosition,
                    trajectoryEndAbsolutePosition,
                    trajectoryEndBarycentricCoords,
                    true, // Check immediately whether we're directed towards the interior of the triangle
                    shipMesh,
                    mParticles,
                    gameParameters);

                switch (outcome.Type)
                {
                    case NavigateVertexOutcome::OutcomeType::CompletedNavigation:
                    {
                        // Continue with ray-tracing
                        break;
                    }

                    case NavigateVertexOutcome::OutcomeType::ConvertedToFree:
                    {
                        // Done
                        remainingDt = 0.0f;
                        break;
                    }

                    case NavigateVertexOutcome::OutcomeType::EncounteredFloor:
                    {
                        // Continue with ray-tracing, handling impact using "flattening" algorithm
                        break;
                    }
                }
            }

            if (remainingDt == 0.0f)
            {
                break;
            }

            //
            // Check which of two cases applies at this moment:
            // 1. On edge and moving against (or along) it (i.e. in a non-inertial frame)
            // 2. Not on edge, or on edge but not moving against (nor along) it
            //

            LogNpcDebug("    Checking edges");

            // Check all edges, stop at first one that is floor and against which we're moving
            bool hasFoundNonInertial = false;
            for (int vi = 0; vi < 3; ++vi)
            {
                if (npcParticle.ConstrainedState->CurrentTriangleBarycentricCoords[vi] == 0.0f)
                {
                    // We are on this edge

                    int const edgeOrdinal = (vi + 1) % 3;

                    assert(isPrimaryParticle || npc.DipoleState.has_value());

                    LogNpcDebug("      edge ", edgeOrdinal, ": isFloor=", IsEdgeFloorToParticle(npcParticle.ConstrainedState->CurrentTriangle, edgeOrdinal, shipMesh));

                    // Check if this is really a floor to this particle
                    if (IsEdgeFloorToParticle(npcParticle.ConstrainedState->CurrentTriangle, edgeOrdinal, shipMesh)
                        && (isPrimaryParticle || !DoesFloorSeparateFromPrimaryParticle(
                            mParticles.GetPosition(npc.DipoleState->SecondaryParticleState.ParticleIndex),
                            trajectoryStartAbsolutePosition, // Current (virtual, not yet real) position of this (secondary) particle
                            npcParticle.ConstrainedState->CurrentTriangle,
                            edgeOrdinal,
                            shipMesh)))
                    {
                        // On floor edge - so potentially in a non-inertial frame

                        // Check now whether we're moving *against* the floor

                        vec2f const edgeVector = shipMesh.GetTriangles().GetSubSpringVector(
                            npcParticle.ConstrainedState->CurrentTriangle,
                            edgeOrdinal,
                            shipMesh.GetPoints());

                        vec2f const edgeDir = edgeVector.normalise();
                        vec2f const edgeNormal = edgeDir.to_perpendicular(); // Points outside of triangle (i.e. towards floor)

                        if (trajectory.dot(edgeNormal) >= 0.0f) // If 0, no normal force - hence no friction; however we want to take this codepath anyways for consistency
                        {
                            //
                            // Case 1: Non-inertial: on edge and moving against (or along) it, pushed by it
                            //

                            LogNpcDebug("    ConstrainedNonInertial: triangle=", npcParticle.ConstrainedState->CurrentTriangle, " edgeOrdinal=", edgeOrdinal, " bCoords=", npcParticle.ConstrainedState->CurrentTriangleBarycentricCoords, " trajectory=", trajectory);
                            LogNpcDebug("    StartPosition=", mParticles.GetPosition(npcParticle.ParticleIndex), " StartVelocity=", mParticles.GetVelocity(npcParticle.ParticleIndex), " MeshVelocity=", meshVelocity, " StartMRVelocity=", npcParticle.ConstrainedState->MeshRelativeVelocity);

                            npcParticle.ConstrainedState->CurrentVirtualEdgeOrdinal = edgeOrdinal;

                            //
                            // We're moving against the floor, hence we are in a non-inertial frame...
                            // ...take friction into account and flaten trajectory
                            //

                            //
                            // Calculate magnitude of flattened trajectory - i.e. component of trajectory
                            // along (i.e. tangent) edge, positive when in the direction of the edge
                            //

                            float trajectoryT = trajectory.dot(edgeDir);

                            //
                            // Update tangential component of trajectory with friction
                            //

                            {
                                // Normal trajectory: apparent (integrated) force against the floor;
                                // positive when against the floor

                                float const trajectoryN = trajectory.dot(edgeDir.to_perpendicular());

                                //
                                // Choose between kinetic and static friction
                                //
                                // Note that we want to check actual velocity here, not
                                // *intention* to move (which is trajectory)
                                //

                                float frictionCoefficient;
                                if (std::abs(npcParticle.ConstrainedState->MeshRelativeVelocity.dot(edgeDir)) > 0.01f) // Magic number
                                {
                                    // Kinetic friction
                                    frictionCoefficient = mParticles.GetPhysicalProperties(npcParticle.ParticleIndex).KineticFriction * gameParameters.KineticFrictionAdjustment;
                                }
                                else
                                {
                                    // Static friction
                                    frictionCoefficient = mParticles.GetPhysicalProperties(npcParticle.ParticleIndex).StaticFriction * gameParameters.StaticFrictionAdjustment;
                                }

                                // Calculate friction (integrated) force magnitude (along edgeDir),
                                // which is the same as apparent force, up to max friction threshold
                                float tFriction = std::min(std::abs(trajectoryT), frictionCoefficient * std::max(trajectoryN, 0.0f));
                                if (trajectoryT >= 0.0f)
                                {
                                    tFriction *= -1.0f;
                                }

                                LogNpcDebug("        friction: trajectoryN=", trajectoryN, " relVel=", npcParticle.ConstrainedState->MeshRelativeVelocity,
                                    " trajectoryT=", trajectoryT, " tFriction=", tFriction);

                                // Update trajectory with friction
                                trajectoryT += tFriction;
                            }

                            // The (signed) edge length that we plan to travel independently from walking
                            float const edgePhysicalTraveledPlanned = trajectoryT;

                            //
                            // Calculate displacement due to walking, if any
                            //

                            // The (signed) edge length that we plan to travel exclusively via walking
                            float edgeWalkedPlanned = 0.0f;

                            if (npc.Kind == NpcKindType::Human
                                && npc.KindSpecificState.HumanNpcState.CurrentBehavior == StateType::KindSpecificStateType::HumanNpcStateType::BehaviorType::Constrained_Walking
                                && isPrimaryParticle)
                            {
                                //
                                // Walking displacement projected along edge - what's needed to reach total desired walking
                                // displacement, together with physical
                                // - Never reduces physical
                                // - Never adds to physical so much as to cause resultant to be faster than walking speed
                                //

                                vec2f const idealWalkDir = vec2f(npc.KindSpecificState.HumanNpcState.CurrentFaceDirectionX, 0.0f);
                                assert(idealWalkDir.length() == 1.0f);

                                float const idealWalkMagnitude = CalculateActualHumanWalkingAbsoluteSpeed(npc.KindSpecificState.HumanNpcState, gameParameters) * remainingDt;

                                vec2f walkDir; // Actual absolute direction of walk - along the edge
                                if (idealWalkDir.dot(edgeDir) >= 0.0f)
                                {
                                    // Same direction as edge (ahead is towards larger)

                                    walkDir = edgeDir;
                                    edgeWalkedPlanned = Clamp(idealWalkMagnitude - edgePhysicalTraveledPlanned, 0.0f, idealWalkMagnitude);
                                }
                                else
                                {
                                    // Opposite direction as edge (ahead is towards smaller)

                                    walkDir = -edgeDir;
                                    edgeWalkedPlanned = Clamp(-idealWalkMagnitude - edgePhysicalTraveledPlanned, -idealWalkMagnitude, 0.0f);
                                }

                                // Apply gravity resistance: too steep slopes (wrt vertical) are gently clamped to zero,
                                // to prevent walking on floors that are too steep
                                float constexpr ResistanceSlopeStart = 0.50f; // Start slightly before expected 45-degree ramp
                                float constexpr ResistanceSlopeEnd = 0.85f; // Max slope we're willing to climb
                                float const x = walkDir.dot(-GameParameters::GravityDir); // TODO: perf: this is just y
                                if (x >= ResistanceSlopeStart)
                                {
                                    float const x2 = (x - ResistanceSlopeStart) / (ResistanceSlopeEnd - ResistanceSlopeStart);
                                    float const gravityResistance = std::max(1.0f - (x2 * x2), 0.0f);

                                    edgeWalkedPlanned *= gravityResistance;
                                }

                                if (npc.KindSpecificState.HumanNpcState.CurrentBehaviorState.Constrained_Walking.CurrentWalkMagnitude != 0.0f)
                                {
                                    LogNpcDebug("        idealWalkMagnitude=", idealWalkMagnitude, " => edgeWalkedPlanned=", edgeWalkedPlanned, " (@", npc.KindSpecificState.HumanNpcState.CurrentBehaviorState.Constrained_Walking.CurrentWalkMagnitude, ")");
                                }
                            }

                            //
                            // Calculate total (signed) displacement we plan on undergoing
                            //

                            float const edgeTraveledPlanned = edgePhysicalTraveledPlanned + edgeWalkedPlanned; // Resultant

                            if (!edgeDistanceToTravelMax)
                            {
                                edgeDistanceToTravelMax = std::abs(edgeTraveledPlanned);

                                LogNpcDebug("        initialized distance budget: edgeDistanceToTravelMax=", *edgeDistanceToTravelMax);
                            }

                            // Make sure we don't travel more than what we're willing to

                            float adjustedEdgeTraveledPlanned;

                            float const remainingDistanceBudget = *edgeDistanceToTravelMax - edgeDistanceTraveledTotal;
                            assert(remainingDistanceBudget >= 0.0f);
                            if (std::abs(edgeTraveledPlanned) > remainingDistanceBudget)
                            {
                                if (edgeTraveledPlanned >= 0.0f)
                                {
                                    adjustedEdgeTraveledPlanned = std::min(edgeTraveledPlanned, remainingDistanceBudget);
                                }
                                else
                                {
                                    adjustedEdgeTraveledPlanned = std::max(edgeTraveledPlanned, -remainingDistanceBudget);
                                }

                                LogNpcDebug("        travel exceeds budget (edgeTraveledPlanned=", edgeTraveledPlanned, " budget=", remainingDistanceBudget,
                                    " => adjustedEdgeTraveledPlanned=", adjustedEdgeTraveledPlanned);
                            }
                            else
                            {
                                adjustedEdgeTraveledPlanned = edgeTraveledPlanned;
                            }

                            //
                            // Recover flattened trajectory as a vector
                            //

                            vec2f const flattenedTrajectory = edgeDir * adjustedEdgeTraveledPlanned;

                            //
                            // Calculate trajectory target
                            //

                            vec2f flattenedTrajectoryEndAbsolutePosition = trajectoryStartAbsolutePosition + flattenedTrajectory;

                            bcoords3f flattenedTrajectoryEndBarycentricCoords = shipMesh.GetTriangles().ToBarycentricCoordinates(
                                flattenedTrajectoryEndAbsolutePosition,
                                npcParticle.ConstrainedState->CurrentTriangle,
                                shipMesh.GetPoints());

                            //
                            // Due to numerical slack, ensure target barycentric coords are still along edge
                            //

                            int const edgeVertexOrdinal = (edgeOrdinal + 2) % 3;
                            flattenedTrajectoryEndBarycentricCoords[edgeVertexOrdinal] = 0.0f;
                            flattenedTrajectoryEndBarycentricCoords[(edgeVertexOrdinal + 1) % 3] = 1.0f - flattenedTrajectoryEndBarycentricCoords[(edgeVertexOrdinal + 2) % 3];

                            LogNpcDebug("        flattenedTrajectory=", flattenedTrajectory, " flattenedTrajectoryEndAbsolutePosition=", flattenedTrajectoryEndAbsolutePosition, " flattenedTrajectoryEndBarycentricCoords=", flattenedTrajectoryEndBarycentricCoords);

                            //
                            // Ray-trace using non-inertial physics;
                            // will return when completed or when current edge is over
                            //
                            // If needs to continue, returns the (signed) actual edge traveled, which is implicitly the
                            // (signed) actual edge physically traveled plus the (signed) actual edge walked during the
                            // consumed dt portion of the remaning dt
                            //
                            // Fact: at each iteration, the actual movement of the particle will be the result of phys traj and imposed walk displacement
                            // Fact: phys traj displacement (planned) is itself dependant from remaining_dt, because of advancement of particle's current bary coords
                            // Fat: walk displacement (planned) is also dependant from remaining_dt, because we use walk velocity
                            // Fact: so, the actual movement includes the consumed_dt's portion (fraction) of both phys traj and imposed walk
                            //

                            LogNpcDebug("        edgePhysicalTraveledPlanned=", edgePhysicalTraveledPlanned, " edgeWalkedPlanned=", edgeWalkedPlanned);

                            auto const [edgeTraveledActual, doStop] = UpdateNpcParticle_ConstrainedNonInertial(
                                npc,
                                isPrimaryParticle,
                                edgeOrdinal,
                                edgeDir,
                                particleStartAbsolutePosition,
                                trajectoryStartAbsolutePosition,
                                flattenedTrajectoryEndAbsolutePosition,
                                flattenedTrajectoryEndBarycentricCoords,
                                flattenedTrajectory,
                                adjustedEdgeTraveledPlanned,
                                meshVelocity,
                                remainingDt,
                                shipMesh,
                                mParticles,
                                gameParameters);

                            LogNpcDebug("    Actual edge traveled in non-inertial step: ", edgeTraveledActual);

                            if (doStop)
                            {
                                // Completed
                                remainingDt = 0.0f;
                            }
                            else
                            {
                                if (edgeTraveledActual == 0.0f)
                                {
                                    // No movement

                                    // Check if we're in a well
                                    if (PastBarycentricPosition tmp(npcParticle.ConstrainedState->CurrentTriangle, npcParticle.ConstrainedState->CurrentTriangleBarycentricCoords);
                                        pastBarycentricPositions[0] == tmp || pastBarycentricPositions[1] == tmp)
                                    {
                                        //
                                        // Well - stop here
                                        //

                                        LogNpcDebug("    Detected well - stopping here");

                                        // Update particle's physics, considering that we are in a well and thus still (wrt mesh)

                                        vec2f const particleEndAbsolutePosition = shipMesh.GetTriangles().FromBarycentricCoordinates(
                                            npcParticle.ConstrainedState->CurrentTriangleBarycentricCoords,
                                            npcParticle.ConstrainedState->CurrentTriangle,
                                            shipMesh.GetPoints());

                                        mParticles.SetPosition(npcParticle.ParticleIndex, particleEndAbsolutePosition);

                                        // No (relative) velocity
                                        mParticles.SetVelocity(npcParticle.ParticleIndex, -meshVelocity);
                                        npcParticle.ConstrainedState->MeshRelativeVelocity = vec2f::zero();

                                        // Consume the whole time quantum
                                        remainingDt = 0.0f;
                                    }
                                    else
                                    {
                                        // Resultant of physical and walked is zero => no movement at all
                                        // Note: either physical and/or walked could individually be non-zero
                                        //
                                        // For now, assume no dt consumed and continue
                                    }
                                }
                                else
                                {
                                    // We have moved

                                    // Calculate consumed dt
                                    assert(edgeTraveledActual * adjustedEdgeTraveledPlanned >= 0.0f); // Should have same sign
                                    float const dtFractionConsumed = adjustedEdgeTraveledPlanned != 0.0f
                                        ? std::min(edgeTraveledActual / adjustedEdgeTraveledPlanned, 1.0f) // Signs should agree anyway
                                        : 1.0f; // If we were planning no travel, any movement is a whole consumption
                                    LogNpcDebug("        dtFractionConsumed=", dtFractionConsumed);
                                    remainingDt *= (1.0f - dtFractionConsumed);

                                    // Reset well detection machinery
                                    std::for_each(
                                        pastBarycentricPositions.begin(),
                                        pastBarycentricPositions.end(),
                                        [](auto & bp) { bp.reset(); });
                                }
                            }

                            // Update total (absolute) distance traveled along (an) edge
                            edgeDistanceTraveledTotal += std::abs(edgeTraveledActual);

                            // If we haven't completed, there is still some distance remaining in the budget
                            //
                            // Note: if this doesn't hold, at the next iteration we'll move by zero and we'll reset
                            // velocity to zero, even though we have moved in this step, thus yielding an erroneous
                            // zero velocity
                            assert(remainingDt == 0.0f || edgeDistanceTraveledTotal < *edgeDistanceToTravelMax);

                            // Update total human distance traveled
                            if (npc.Kind == NpcKindType::Human
                                && isPrimaryParticle) // Human is represented by primary particle
                            {
                                // Note: does not include eventual disance traveled after becoming free; fine because we will transition and wipe out total traveled
                                npc.KindSpecificState.HumanNpcState.TotalDistanceTraveledOnEdgeSinceStateTransition += std::abs(edgeTraveledActual);
                            }

                            // Update total vector walked along edge
                            // Note: we use unadjusted edge traveled planned, as edge walked planned is also unadjusted,
                            // and we are only interested in the ratio anyway
                            float const edgeWalkedActual = edgeTraveledPlanned != 0.0f
                                ? edgeTraveledActual * (edgeWalkedPlanned / edgeTraveledPlanned)
                                : 0.0f; // Unlikely, but read above for rationale behind 0.0
                            totalEdgeWalkedActual += edgeDir * edgeWalkedActual;
                            LogNpcDebug("        edgeWalkedActual=", edgeWalkedActual, " totalEdgeWalkedActual=", totalEdgeWalkedActual);

                            // Update well detection machinery
                            if (npcParticle.ConstrainedState.has_value()) // We might have left constrained state (not to return to it anymore)
                            {
                                pastBarycentricPositions[1] = pastBarycentricPositions[0];
                                pastBarycentricPositions[0].emplace(npcParticle.ConstrainedState->CurrentTriangle, npcParticle.ConstrainedState->CurrentTriangleBarycentricCoords);
                            }

                            if (npcParticle.ConstrainedState.has_value())
                            {
                                LogNpcDebug("    EndPosition=", mParticles.GetPosition(npcParticle.ParticleIndex), " EndVelocity=", mParticles.GetVelocity(npcParticle.ParticleIndex), " EndMRVelocity=", npcParticle.ConstrainedState->MeshRelativeVelocity);
                            }
                            else
                            {
                                LogNpcDebug("    EndPosition=", mParticles.GetPosition(npcParticle.ParticleIndex), " EndVelocity=", mParticles.GetVelocity(npcParticle.ParticleIndex));
                            }

                            hasFoundNonInertial = true;

                            break;
                        }
                        else
                        {
                            LogNpcDebug("      traj.edgeN=", trajectory.dot(edgeNormal));
                        }
                    }
                }
            } // for (vertex)

            if (!hasFoundNonInertial)
            {
                //
                // Case 2: Inertial: not on edge or on edge but not moving against (nor along) it
                //

                LogNpcDebug("    ConstrainedInertial: triangle=", npcParticle.ConstrainedState->CurrentTriangle, " bCoords=", npcParticle.ConstrainedState->CurrentTriangleBarycentricCoords, " physicsDeltaPos=", physicsDeltaPos);
                LogNpcDebug("    StartPosition=", mParticles.GetPosition(npcParticle.ParticleIndex), " StartVelocity=", mParticles.GetVelocity(npcParticle.ParticleIndex), " MeshVelocity=", meshVelocity, " StartMRVelocity=", npcParticle.ConstrainedState->MeshRelativeVelocity);

                npcParticle.ConstrainedState->CurrentVirtualEdgeOrdinal = -1;

                //
                // Calculate target barycentric coords
                //

                bcoords3f const trajectoryEndBarycentricCoords = shipMesh.GetTriangles().ToBarycentricCoordinates(
                    trajectoryEndAbsolutePosition,
                    npcParticle.ConstrainedState->CurrentTriangle,
                    shipMesh.GetPoints());

                //
                // Move towards target bary coords
                //

                float totalTraveled = UpdateNpcParticle_ConstrainedInertial(
                    npc,
                    isPrimaryParticle,
                    particleStartAbsolutePosition,
                    npcParticle.ConstrainedState->CurrentTriangleBarycentricCoords, // segmentTrajectoryStartBarycentricCoords
                    trajectoryEndAbsolutePosition,
                    trajectoryEndBarycentricCoords,
                    meshVelocity,
                    remainingDt,
                    shipMesh,
                    mParticles,
                    gameParameters);

                // Update total traveled
                if (npc.Kind == NpcKindType::Human
                    && isPrimaryParticle) // Human is represented by primary particle
                {
                    // Note: does not include eventual disance traveled after becoming free; fine because we will transition and wipe out total traveled
                    npc.KindSpecificState.HumanNpcState.TotalDistanceTraveledOffEdgeSinceStateTransition += std::abs(totalTraveled);
                }

                if (npcParticle.ConstrainedState.has_value())
                {
                    LogNpcDebug("    EndPosition=", mParticles.GetPosition(npcParticle.ParticleIndex), " EndVelocity=", mParticles.GetVelocity(npcParticle.ParticleIndex), " EndMRVelocity=", npcParticle.ConstrainedState->MeshRelativeVelocity);
                }
                else
                {
                    LogNpcDebug("    EndPosition=", mParticles.GetPosition(npcParticle.ParticleIndex), " EndVelocity=", mParticles.GetVelocity(npcParticle.ParticleIndex));
                }

                // Consume whole time quantum and stop
                remainingDt = 0.0f;
            }

            if (remainingDt <= 0.0f)
            {
                assert(remainingDt > -0.0001f); // If negative, it's only because of numerical slack

                // Consumed whole time quantum, loop completed
                break;
            }

            //
            // Update trajectory start position for next iteration
            //

            // Current (virtual, not yet real) position of this particle
            trajectoryStartAbsolutePosition = shipMesh.GetTriangles().FromBarycentricCoordinates(
                npcParticle.ConstrainedState->CurrentTriangleBarycentricCoords,
                npcParticle.ConstrainedState->CurrentTriangle,
                shipMesh.GetPoints());
        }
    }

#ifdef IN_BARYLAB
    if (mCurrentlySelectedParticle == npcParticle.ParticleIndex)
    {
        // Publish final velocities

        vec2f const particleVelocity = (mParticles.GetPosition(npcParticle.ParticleIndex) - particleStartAbsolutePosition) / GameParameters::SimulationTimeStepDuration;

        mGameEventHandler->OnCustomProbe("VelX", particleVelocity.x);
        mGameEventHandler->OnCustomProbe("VelY", particleVelocity.y);
    }
#endif
}

void Npcs::CalculateNpcParticlePreliminaryForces(
    StateType const & npc,
    bool isPrimaryParticle,
    GameParameters const & gameParameters)
{
    auto & npcParticle = isPrimaryParticle ? npc.PrimaryParticleState : npc.DipoleState->SecondaryParticleState;

    //
    // Calculate world forces
    //

    float const particleMass = mParticles.GetPhysicalProperties(npcParticle.ParticleIndex).Mass * gameParameters.MassAdjustment;

    // 1. World forces - gravity

    vec2f preliminaryForces = GameParameters::Gravity * gameParameters.GravityAdjustment * mGravityGate * particleMass;

    if (!npcParticle.ConstrainedState.has_value())
    {
        // Check whether we are underwater

        float constexpr BuoyancyInterfaceWidth = 0.4f;

        vec2f testParticlePosition = mParticles.GetPosition(npcParticle.ParticleIndex);
        if (npc.Kind == NpcKindType::Human && !isPrimaryParticle)
        {
            // Head - use an offset
            assert(npc.DipoleState.has_value());
            testParticlePosition += (mParticles.GetPosition(npc.PrimaryParticleState.ParticleIndex) - testParticlePosition) * BuoyancyInterfaceWidth * 2.0f / 3.0f;
        }

        float const particleDepth = mParentWorld.GetOceanSurface().GetDepth(testParticlePosition);
        float const uwCoefficient = Clamp(particleDepth, 0.0f, BuoyancyInterfaceWidth) / BuoyancyInterfaceWidth;
        if (uwCoefficient > 0.0f)
        {
            // Underwater

            // 2. World forces - buoyancy

            preliminaryForces.y +=
                GameParameters::GravityMagnitude * 1000.0f
                * mParticles.GetPhysicalProperties(npcParticle.ParticleIndex).BuoyancyVolumeFill
                * gameParameters.BuoyancyAdjustment
                * uwCoefficient;

            // 3. World forces - water drag

            preliminaryForces +=
                -mParticles.GetVelocity(npcParticle.ParticleIndex)
                * GameParameters::WaterFrictionDragCoefficient
                * gameParameters.WaterFrictionDragCoefficientAdjustment;
        }
    }

    // 3. Spring forces

    if (npc.DipoleState.has_value())
    {
        ElementIndex const otherParticleIndex = isPrimaryParticle ? npc.DipoleState->SecondaryParticleState.ParticleIndex : npc.PrimaryParticleState.ParticleIndex;

        float const dt = GameParameters::SimulationTimeStepDuration;

        vec2f const springDisplacement = mParticles.GetPosition(otherParticleIndex) - mParticles.GetPosition(npcParticle.ParticleIndex); // Towards other
        float const springDisplacementLength = springDisplacement.length();
        vec2f const springDir = springDisplacement.normalise_approx(springDisplacementLength);

        //
        // 3a. Hooke's law
        //

        float const springStiffnessCoefficient =
            gameParameters.SpringReductionFraction
            * npc.DipoleState->DipoleProperties.MassFactor * gameParameters.MassAdjustment
            / (dt * dt);

        // Calculate spring force on this particle
        float const fSpring =
            (springDisplacementLength - npc.DipoleState->DipoleProperties.DipoleLength)
            * springStiffnessCoefficient;

        //
        // 3b. Damper forces
        //
        // Damp the velocities of each endpoint pair, as if the points were also connected by a damper
        // along the same direction as the spring
        //

        float const springDampingCoefficient =
            gameParameters.SpringDampingCoefficient
            * npc.DipoleState->DipoleProperties.MassFactor * gameParameters.MassAdjustment
            / dt;

        // Calculate damp force on this particle
        vec2f const relVelocity = mParticles.GetVelocity(otherParticleIndex) - mParticles.GetVelocity(npcParticle.ParticleIndex);
        float const fDamp =
            relVelocity.dot(springDir)
            * springDampingCoefficient;

        //
        // Apply forces
        //

        preliminaryForces += springDir * (fSpring + fDamp);
    }

    // 4. External forces

    preliminaryForces += mParticles.GetExternalForces(npcParticle.ParticleIndex);


    mParticles.SetPreliminaryForces(npcParticle.ParticleIndex, preliminaryForces);
}

vec2f Npcs::CalculateNpcParticleDefinitiveForces(
    StateType const & npc,
    bool isPrimaryParticle,
    float particleMass,
    GameParameters const & gameParameters) const
{
    auto & npcParticle = isPrimaryParticle ? npc.PrimaryParticleState : npc.DipoleState->SecondaryParticleState;

    vec2f definitiveForces = mParticles.GetPreliminaryForces(npcParticle.ParticleIndex);

    //
    // Human Equlibrium Torque
    //

    if (npc.Kind == NpcKindType::Human
        && (npc.KindSpecificState.HumanNpcState.CurrentBehavior == StateType::KindSpecificStateType::HumanNpcStateType::BehaviorType::Constrained_Walking
            || npc.KindSpecificState.HumanNpcState.CurrentBehavior == StateType::KindSpecificStateType::HumanNpcStateType::BehaviorType::Constrained_Equilibrium
            || npc.KindSpecificState.HumanNpcState.CurrentBehavior == StateType::KindSpecificStateType::HumanNpcStateType::BehaviorType::Constrained_Rising)
        && !isPrimaryParticle)
    {
        assert(npc.DipoleState.has_value());

        ElementIndex const primaryParticleIndex = npc.PrimaryParticleState.ParticleIndex;
        ElementIndex const secondaryParticleIndex = npc.DipoleState->SecondaryParticleState.ParticleIndex;

        // Given that we apply torque onto the secondary particle *after* the primary has been simulated
        // (so that we take into account the primary's new position), and thus we see the primary where it
        // is at the end of the step - possibly far away if mesh velocity is high, we want to _predict_
        // where the secondary will be by its own velocity

        vec2f const secondaryPredictedPosition =
            mParticles.GetPosition(secondaryParticleIndex)
            + mParticles.GetVelocity(secondaryParticleIndex) * GameParameters::SimulationTimeStepDuration;
            ;
        vec2f const humanVector = secondaryPredictedPosition - mParticles.GetPosition(primaryParticleIndex);

        // Calculate CW angle between head and vertical (pointing up);
        // positive when human is CW wrt vertical
        //
        // |   H
        // |  /
        // |-/
        // |/
        //
        float const staticDisplacementAngleCW = (-GameParameters::GravityDir).angleCw(humanVector);

        // Calculate CW angle that head would rotate by (relative to feet) due to relative velocity alone;
        // positive when new position is CW wrt old
        //
        // |   H
        // |  /
        // | /\
        // |/__L___H'
        //
        vec2f const relativeVelocityDisplacement =
            (mParticles.GetVelocity(secondaryParticleIndex) - mParticles.GetVelocity(primaryParticleIndex))
            * GameParameters::SimulationTimeStepDuration;
        float const relativeVelocityAngleCW = humanVector.angleCw(humanVector + relativeVelocityDisplacement);

        //
        // Calculate the torque on the secondary (head)
        // required to maintain alignment with vertical
        //

        // Calculate angle that we want to enforce with this torque
        float const stiffness = gameParameters.HumanNpcEquilibriumTorqueStiffnessCoefficient + std::min(npc.KindSpecificState.HumanNpcState.ResultantPanicLevel, 1.0f) * 0.0005f;
        float const totalTorqueAngleCW =
            staticDisplacementAngleCW * stiffness
            + relativeVelocityAngleCW * gameParameters.HumanNpcEquilibriumTorqueDampingCoefficient;

        // Calculate (linear) force that generates this rotation
        vec2f const torqueDisplacement = humanVector.rotate(totalTorqueAngleCW) - humanVector;
        vec2f const equilibriumTorqueForce =
            torqueDisplacement
            * particleMass / (GameParameters::SimulationTimeStepDuration * GameParameters::SimulationTimeStepDuration)
            * mParticles.GetEquilibriumTorque(secondaryParticleIndex);

        definitiveForces += equilibriumTorqueForce;
    }

    return definitiveForces;
}

void Npcs::UpdateNpcParticle_Free(
    StateType::NpcParticleStateType & particle,
    vec2f const & startPosition,
    vec2f const & endPosition,
    NpcParticles & particles,
    GameParameters const & gameParameters) const
{
    assert(!particle.ConstrainedState.has_value());

    // Update position
    particles.SetPosition(
        particle.ParticleIndex,
        endPosition);

    // Update velocity
    // Use whole time quantum for velocity, as communicated start/end positions are those planned for whole dt
    particles.SetVelocity(
        particle.ParticleIndex,
        (endPosition - startPosition) / GameParameters::SimulationTimeStepDuration * (1.0f - gameParameters.GlobalDamping));
}

std::tuple<float, bool> Npcs::UpdateNpcParticle_ConstrainedNonInertial(
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
    GameParameters const & gameParameters)
{
    auto & npcParticle = isPrimaryParticle ? npc.PrimaryParticleState : npc.DipoleState->SecondaryParticleState;
    assert(npcParticle.ConstrainedState.has_value());
    auto & npcParticleConstrainedState = *npcParticle.ConstrainedState;

    //
    // Ray-trace along the flattened trajectory, ending only at one of the following three conditions:
    // 1. Reached destination: terminate
    // 2. Becoming free: do free movement and terminate
    // 3. Impact:
    //      - With bounce: impart bounce velocity and terminate
    //      - Without bounce: advance simulation by how much walked, re-flatten trajectory, and continue
    //

    //
    // If target is on/in triangle, we move to target
    //

    if (flattenedTrajectoryEndBarycentricCoords.is_on_edge_or_internal())
    {
        LogNpcDebug("      Target is on/in triangle, moving to target");

        //
        // Update particle and exit - consuming whole time quantum
        //

        npcParticleConstrainedState.CurrentTriangleBarycentricCoords = flattenedTrajectoryEndBarycentricCoords;

        vec2f const particleEndAbsolutePosition = shipMesh.GetTriangles().FromBarycentricCoordinates(
            flattenedTrajectoryEndBarycentricCoords,
            npcParticleConstrainedState.CurrentTriangle,
            shipMesh.GetPoints());

        particles.SetPosition(npcParticle.ParticleIndex, particleEndAbsolutePosition);

        //
        // Velocity: given that we've completed *along the edge*, then we can calculate
        // our (relative) velocity based on the distance traveled along this (remaining) time quantum
        //
        // We take into account only the edge traveled at this moment, divided by the length of this time quantum:
        // V = signed_edge_traveled_actual * edgeDir / this_dt
        //
        // Now: consider that at this moment we've reached the planned end of the iteration's sub-trajectory;
        // we can then assume that signed_edge_traveled_actual == signed_edge_traveled_planned (verified via assert)
        //

        assert(std::abs((particleEndAbsolutePosition - trajectoryStartAbsolutePosition).dot(edgeDir) - edgeTraveledPlanned) < 0.001f);

        vec2f const relativeVelocity =
            edgeDir
            * edgeTraveledPlanned
            / dt;

        particles.SetVelocity(npcParticle.ParticleIndex, relativeVelocity - meshVelocity);
        npcParticleConstrainedState.MeshRelativeVelocity = relativeVelocity;

        LogNpcDebug("        edgeTraveleded (==planned)=", edgeTraveledPlanned, " absoluteVelocity=", particles.GetVelocity(npcParticle.ParticleIndex));

        // Complete
        return std::make_tuple(
            edgeTraveledPlanned, // We moved by how much we planned
            true); // Stop here
    }

    //
    // Target is outside triangle
    //

    LogNpcDebug("      Target is outside triangle");

    //
    // Find closest intersection in the direction of the trajectory, which is
    // a vertex of this triangle
    //
    // Guaranteed to exist and within trajectory, because target is outside
    // of triangle and we're on an edge
    //

    // Here we take advantage of the fact that we know we're on an edge

    assert(npcParticleConstrainedState.CurrentTriangleBarycentricCoords[(edgeOrdinal + 2) % 3] == 0.0f);

    int intersectionVertexOrdinal;
    if (flattenedTrajectoryEndBarycentricCoords[edgeOrdinal] < 0.0f)
    {
        // It's the vertex at the end of our edge
        assert(flattenedTrajectoryEndBarycentricCoords[(edgeOrdinal + 1) % 3] > 0.0f); // Because target is outside of triangle
        intersectionVertexOrdinal = (edgeOrdinal + 1) % 3;
    }
    else
    {
        // It's the vertex before our edge
        assert(flattenedTrajectoryEndBarycentricCoords[(edgeOrdinal + 1) % 3] < 0.0f);
        intersectionVertexOrdinal = edgeOrdinal;
    }

    //
    // Move to intersection (which is a vertex by now)
    //
    // Note that except for the very first iteration, any other iteration will travel
    // zero distance at this moment
    //

    bcoords3f intersectionBarycentricCoords = bcoords3f::zero();
    intersectionBarycentricCoords[intersectionVertexOrdinal] = 1.0f;

    LogNpcDebug("      Moving to intersection vertex ", intersectionVertexOrdinal, ": ", intersectionBarycentricCoords);

    npcParticleConstrainedState.CurrentTriangleBarycentricCoords = intersectionBarycentricCoords;

    // Update (signed) edge traveled, assuming we're always moving along initial edge
    // (other edges we're only passing through cuspids)
    vec2f const intersectionAbsolutePosition = shipMesh.GetPoints().GetPosition(shipMesh.GetTriangles().GetPointIndices(npcParticleConstrainedState.CurrentTriangle)[intersectionVertexOrdinal]);
    assert(intersectionAbsolutePosition == shipMesh.GetTriangles().FromBarycentricCoordinates(
        intersectionBarycentricCoords,
        npcParticleConstrainedState.CurrentTriangle,
        shipMesh.GetPoints()));

    float const edgeTraveled = (intersectionAbsolutePosition - trajectoryStartAbsolutePosition).dot(edgeDir);

    LogNpcDebug("        edgeTraveled=", edgeTraveled);

    //
    // Navigate this vertex now, until any of these:
    // - We are directed inside a triangle
    // - We hit a floor
    // - We become free
    //

    auto const navigationOutcome = NavigateVertex(
        npc,
        isPrimaryParticle,
        intersectionVertexOrdinal,
        particleStartAbsolutePosition,
        trajectoryStartAbsolutePosition,
        flattenedTrajectoryEndAbsolutePosition,
        flattenedTrajectoryEndBarycentricCoords,
        false, // No need to check whether we are directed into _this_ triangle
        shipMesh,
        particles,
        gameParameters);

    switch (navigationOutcome.Type)
    {
        case NavigateVertexOutcome::OutcomeType::CompletedNavigation:
        {
            return std::make_tuple(edgeTraveled, false);
        }

        case NavigateVertexOutcome::OutcomeType::ConvertedToFree:
        {
            return std::make_tuple(edgeTraveled, true);
        }

        case NavigateVertexOutcome::OutcomeType::EncounteredFloor:
        {
            //
            // We might have hit a tiny bump (e.g. because of triangles slightly bent); in this case we don't want to bounce
            //

            float const flattenedTrajectoryLength = flattenedTrajectory.length();

            vec2f const floorEdgeDir =
                shipMesh.GetTriangles().GetSubSpringVector(npcParticleConstrainedState.CurrentTriangle, navigationOutcome.EncounteredFloorEdgeOrdinal, shipMesh.GetPoints())
                .normalise();
            vec2f const floorEdgeNormal = floorEdgeDir.to_perpendicular();

            // Check angle between desired (original) trajectory and edge
            float const trajProjOntoEdgeNormal = flattenedTrajectory.normalise().dot(floorEdgeNormal);
            if (trajProjOntoEdgeNormal <= 0.71f) // PI/4+
            {
                //
                // Impact continuation (no bounce)
                //
                // Stop here and then check trajectory in new situation
                //

                LogNpcDebug("      Impact continuation (trajProjOntoEdgeNormal=", trajProjOntoEdgeNormal, ")");

                return std::make_tuple(edgeTraveled, false);
            }
            else
            {
                //
                // Bounce - calculate bounce response, using the *apparent* (trajectory)
                // velocity - since this one includes the mesh velocity
                //

                LogNpcDebug("      Bounce (trajProjOntoEdgeNormal=", trajProjOntoEdgeNormal, ")");

                BounceConstrainedNpcParticle(
                    npc,
                    isPrimaryParticle,
                    flattenedTrajectory,
                    intersectionAbsolutePosition,
                    floorEdgeNormal,
                    meshVelocity,
                    dt,
                    particles,
                    gameParameters);

                // Terminate
                return std::make_tuple(edgeTraveled, true);
            }
        }
    }

    assert(false);
    return std::make_tuple(edgeTraveled, false);
}

float Npcs::UpdateNpcParticle_ConstrainedInertial(
    StateType & npc,
    bool isPrimaryParticle,
    vec2f const & particleStartAbsolutePosition, // Since beginning of whole time quantum, not just this step
    bcoords3f const segmentTrajectoryStartBarycentricCoords, // In current triangle
    vec2f const & segmentTrajectoryEndAbsolutePosition,
    bcoords3f segmentTrajectoryEndBarycentricCoords, // In current triangle; mutable
    vec2f const meshVelocity,
    float segmentDt,
    Ship const & shipMesh,
    NpcParticles & particles,
    GameParameters const & gameParameters)
{
    auto & npcParticle = isPrimaryParticle ? npc.PrimaryParticleState : npc.DipoleState->SecondaryParticleState;
    assert(npcParticle.ConstrainedState.has_value());
    auto & npcParticleConstrainedState = *npcParticle.ConstrainedState;

    vec2f const segmentTrajectoryStartAbsolutePosition = shipMesh.GetTriangles().FromBarycentricCoordinates(
        segmentTrajectoryStartBarycentricCoords,
        npcParticle.ConstrainedState->CurrentTriangle,
        shipMesh.GetPoints());

    //
    // Ray-trace along the specified trajectory, ending only at one of the following three conditions:
    // 1. Reached destination: terminate
    // 2. Becoming free: do free movement and terminate
    // 3. Impact with bounce: impart bounce velocity and terminate
    //

    for (int iIter = 0; ; ++iIter)
    {
        assert(npcParticleConstrainedState.CurrentTriangleBarycentricCoords.is_on_edge_or_internal());

        LogNpcDebug("    SegmentTrace ", iIter);
        LogNpcDebug("      triangle=", npcParticleConstrainedState.CurrentTriangle, " bCoords=", npcParticleConstrainedState.CurrentTriangleBarycentricCoords,
            " segmentTrajStartBCoords=", segmentTrajectoryStartBarycentricCoords, " segmentTrajEndBCoords=", segmentTrajectoryEndBarycentricCoords);

        //
        // If target is on/in triangle, we move to target
        //

        if (segmentTrajectoryEndBarycentricCoords.is_on_edge_or_internal())
        {
            LogNpcDebug("      Target is on/in triangle, moving to target");

            //
            // Update particle and exit - consuming whole time quantum
            //

            // Move particle to end of trajectory
            npcParticleConstrainedState.CurrentTriangleBarycentricCoords = segmentTrajectoryEndBarycentricCoords;
            particles.SetPosition(npcParticle.ParticleIndex, segmentTrajectoryEndAbsolutePosition);

            // Use whole time quantum for velocity, as particleStartAbsolutePosition is fixed at t0
            vec2f const totalAbsoluteTraveledVector = segmentTrajectoryEndAbsolutePosition - particleStartAbsolutePosition;
            vec2f const absoluteVelocity = totalAbsoluteTraveledVector / GameParameters::SimulationTimeStepDuration;
            particles.SetVelocity(npcParticle.ParticleIndex, absoluteVelocity);
            npcParticleConstrainedState.MeshRelativeVelocity = absoluteVelocity + meshVelocity;

            LogNpcDebug("        totalAbsoluteTraveledVector=", totalAbsoluteTraveledVector, " absoluteVelocity=", particles.GetVelocity(npcParticle.ParticleIndex));

            // Return (mesh-relative) distance traveled with this move
            return (segmentTrajectoryEndAbsolutePosition - segmentTrajectoryStartAbsolutePosition).length();
        }

        //
        // We're inside or on triangle, and target is outside triangle;
        // if we're on edge, trajectory is along this edge
        //

        LogNpcDebug("      Target is outside triangle");

        //
        // Find closest intersection in the direction of the trajectory
        //
        // Guaranteed to exist and within trajectory, because target is outside
        // of triangle and we're inside or on an edge
        //

#ifdef _DEBUG

        struct EdgeIntersectionDiag
        {
            float Den;
            float T;
            bcoords3f IntersectionPoint;

            EdgeIntersectionDiag(
                float den,
                float t)
                : Den(den)
                , T(t)
                , IntersectionPoint()
            {}
        };

        std::array<std::optional<EdgeIntersectionDiag>, 3> diags;

#endif // DEBUG

        int intersectionVertexOrdinal = -1;
        float minIntersectionT = std::numeric_limits<float>::max();

        for (int vi = 0; vi < 3; ++vi)
        {
            int const edgeOrdinal = (vi + 1) % 3;

            // Only consider edges that we ancounter ahead along the trajectory
            if (segmentTrajectoryEndBarycentricCoords[vi] < 0.0f)
            {
                float const den = npcParticleConstrainedState.CurrentTriangleBarycentricCoords[vi] - segmentTrajectoryEndBarycentricCoords[vi];
                float const t = npcParticleConstrainedState.CurrentTriangleBarycentricCoords[vi] / den;

#ifdef _DEBUG
                diags[vi].emplace(den, t);
                diags[vi]->IntersectionPoint =
                    npcParticleConstrainedState.CurrentTriangleBarycentricCoords
                    + (segmentTrajectoryEndBarycentricCoords - npcParticleConstrainedState.CurrentTriangleBarycentricCoords) * t;
#endif

                assert(t > -Epsilon<float>); // Some numeric slack, trajectory is here guaranteed to be pointing into this edge

                LogNpcDebug("        t[v", vi, " e", edgeOrdinal, "] = ", t);

                if (t < minIntersectionT)
                {
                    intersectionVertexOrdinal = vi;
                    minIntersectionT = t;
                }
            }
        }

        assert(intersectionVertexOrdinal >= 0); // Guaranteed to exist
        assert(minIntersectionT > -Epsilon<float> && minIntersectionT <= 1.0f); // Guaranteed to exist, and within trajectory

        int const intersectionEdgeOrdinal = (intersectionVertexOrdinal + 1) % 3;

        // Calculate intersection barycentric coordinates

        bcoords3f intersectionBarycentricCoords;
        intersectionBarycentricCoords[intersectionVertexOrdinal] = 0.0f;
        float const lNext = Clamp( // Barycentric coord of next vertex at intersection; enforcing it's within triangle
            npcParticleConstrainedState.CurrentTriangleBarycentricCoords[(intersectionVertexOrdinal + 1) % 3] * (1.0f - minIntersectionT)
            + segmentTrajectoryEndBarycentricCoords[(intersectionVertexOrdinal + 1) % 3] * minIntersectionT,
            0.0f,
            1.0f);
        intersectionBarycentricCoords[(intersectionVertexOrdinal + 1) % 3] = lNext;
        intersectionBarycentricCoords[(intersectionVertexOrdinal + 2) % 3] = 1.0f - lNext;

        assert(intersectionBarycentricCoords.is_on_edge_or_internal());

        //
        // Move to intersection, by moving barycentric coords
        //

        LogNpcDebug("      Moving bary coords to intersection with edge ", intersectionEdgeOrdinal, " ", intersectionBarycentricCoords);

        npcParticleConstrainedState.CurrentTriangleBarycentricCoords = intersectionBarycentricCoords;

        //
        // Check if impacted with floor
        //

        vec2f const intersectionAbsolutePosition = shipMesh.GetTriangles().FromBarycentricCoordinates(
            intersectionBarycentricCoords,
            npcParticleConstrainedState.CurrentTriangle,
            shipMesh.GetPoints());

        assert(isPrimaryParticle || npc.DipoleState.has_value());

        if (IsEdgeFloorToParticle(npcParticleConstrainedState.CurrentTriangle, intersectionEdgeOrdinal, shipMesh)
            && (isPrimaryParticle || !DoesFloorSeparateFromPrimaryParticle(
                mParticles.GetPosition(npc.DipoleState->SecondaryParticleState.ParticleIndex),
                intersectionAbsolutePosition,
                npcParticleConstrainedState.CurrentTriangle,
                intersectionEdgeOrdinal,
                shipMesh)))
        {
            //
            // Impact and bounce
            //

            LogNpcDebug("      Impact and bounce");

            //
            // Calculate bounce response, using the *apparent* (trajectory)
            // velocity - since this one includes the mesh velocity
            //

            vec2f const trajectory = segmentTrajectoryEndAbsolutePosition - segmentTrajectoryStartAbsolutePosition;

            vec2f const intersectionEdgeDir =
                shipMesh.GetTriangles().GetSubSpringVector(npcParticleConstrainedState.CurrentTriangle, intersectionEdgeOrdinal, shipMesh.GetPoints())
                .normalise();
            vec2f const intersectionEdgeNormal = intersectionEdgeDir.to_perpendicular();

            BounceConstrainedNpcParticle(
                npc,
                isPrimaryParticle,
                trajectory,
                intersectionAbsolutePosition,
                intersectionEdgeNormal,
                meshVelocity,
                segmentDt,
                particles,
                gameParameters);

            // Return (mesh-relative) distance traveled with this move
            return (segmentTrajectoryEndAbsolutePosition - segmentTrajectoryStartAbsolutePosition).length();
        }
        else
        {
            //
            // Not floor, climb over edge
            //

            LogNpcDebug("      Climbing over non-floor edge");

            // Find opposite triangle
            auto const & oppositeTriangleInfo = shipMesh.GetTriangles().GetOppositeTriangle(npcParticleConstrainedState.CurrentTriangle, intersectionEdgeOrdinal);
            if (oppositeTriangleInfo.TriangleElementIndex == NoneElementIndex || shipMesh.GetTriangles().IsDeleted(oppositeTriangleInfo.TriangleElementIndex))
            {
                //
                // Become free
                //

                LogNpcDebug("      No opposite triangle found, becoming free");

                //
                // Move to endpoint and exit, consuming whole quantum
                //

                TransitionParticleToFreeState(npc, isPrimaryParticle);

                UpdateNpcParticle_Free(
                    npcParticle,
                    particleStartAbsolutePosition,
                    segmentTrajectoryEndAbsolutePosition,
                    particles,
                    gameParameters);

                vec2f const totalTraveledVector = intersectionAbsolutePosition - particleStartAbsolutePosition; // We consider constrained only
                return totalTraveledVector.length();
            }
            else
            {
                //
                // Move to edge of opposite triangle
                //

                LogNpcDebug("      Moving to edge ", oppositeTriangleInfo.EdgeOrdinal, " of opposite triangle ", oppositeTriangleInfo.TriangleElementIndex);

                // Move to triangle
                npcParticleConstrainedState.CurrentTriangle = oppositeTriangleInfo.TriangleElementIndex;

                // Calculate new current barycentric coords (wrt new triangle)
                bcoords3f newBarycentricCoords; // In new triangle
                newBarycentricCoords[(oppositeTriangleInfo.EdgeOrdinal + 2) % 3] = 0.0f;
                newBarycentricCoords[oppositeTriangleInfo.EdgeOrdinal] = npcParticleConstrainedState.CurrentTriangleBarycentricCoords[(intersectionEdgeOrdinal + 1) % 3];
                newBarycentricCoords[(oppositeTriangleInfo.EdgeOrdinal + 1) % 3] = npcParticleConstrainedState.CurrentTriangleBarycentricCoords[intersectionEdgeOrdinal];

                LogNpcDebug("      B-Coords: ", npcParticleConstrainedState.CurrentTriangleBarycentricCoords, " -> ", newBarycentricCoords);

                assert(newBarycentricCoords.is_on_edge_or_internal());

                npcParticleConstrainedState.CurrentTriangleBarycentricCoords = newBarycentricCoords;

                // Translate target coords to this triangle, for next iteration
                auto const oldSegmentTrajectoryEndBarycentricCoords = segmentTrajectoryEndBarycentricCoords; // For logging
                // Note: here we introduce a lot of error - the target bary coords are not anymore
                // guaranteed to lie exactly on the (continuation of the) edge
                segmentTrajectoryEndBarycentricCoords = shipMesh.GetTriangles().ToBarycentricCoordinates(
                    segmentTrajectoryEndAbsolutePosition,
                    oppositeTriangleInfo.TriangleElementIndex,
                    shipMesh.GetPoints());

                LogNpcDebug("      TrajEndB-Coords: ", oldSegmentTrajectoryEndBarycentricCoords, " -> ", segmentTrajectoryEndBarycentricCoords);

                // Continue
            }
        }
    }
}

Npcs::NavigateVertexOutcome Npcs::NavigateVertex(
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
    GameParameters const & gameParameters)
{
    //
    // Note: we don't move here
    //

    auto & npcParticle = isPrimaryParticle ? npc.PrimaryParticleState : npc.DipoleState->SecondaryParticleState;
    assert(npcParticle.ConstrainedState.has_value());

    for (int iIter = 0; ; ++iIter)
    {
        LogNpcDebug("    NavigateVertex: iter=", iIter);

        assert(iIter < 5); // Detect and break on infinite loops

        // The two vertices around the vertex we are on - seen in clockwise order
        int const nextVertexOrdinal = (vertexOrdinal + 1) % 3;
        int const prevVertexOrdinal = (vertexOrdinal + 2) % 3;

        // Pre-conditions: we are at this vertex
        assert(npcParticle.ConstrainedState->CurrentTriangleBarycentricCoords[vertexOrdinal] == 1.0f);
        assert(npcParticle.ConstrainedState->CurrentTriangleBarycentricCoords[nextVertexOrdinal] == 0.0f);
        assert(npcParticle.ConstrainedState->CurrentTriangleBarycentricCoords[prevVertexOrdinal] == 0.0f);

        LogNpcDebug("      Triangle=", npcParticle.ConstrainedState->CurrentTriangle, " Vertex=", vertexOrdinal, " TrajectoryEndBarycentricCoords=", trajectoryEndBarycentricCoords);

        if (isInitialStateUnknown)
        {
            // Check whether we are directed towards the *interior* of this triangle, or elsewhere
            if (trajectoryEndBarycentricCoords[prevVertexOrdinal] >= 0.0f
                && trajectoryEndBarycentricCoords[nextVertexOrdinal] >= 0.0f)
            {
                //
                // We go inside this triangle - stop where we are, we'll then check trajectory in new situation
                //

                LogNpcDebug("      Trajectory extends inside triangle - CompletedNavigation");

                return NavigateVertexOutcome::MakeCompletedNavigationOutcome();
            }
        }

        //
        // Find next edge that we intersect at this vertex
        //

        int crossedEdgeOrdinal;
        if (trajectoryEndBarycentricCoords[prevVertexOrdinal] <= trajectoryEndBarycentricCoords[nextVertexOrdinal])
        {
            // Clockwise - next edge
            crossedEdgeOrdinal = vertexOrdinal;
        }
        else
        {
            // Anti-clockwise - prev edge
            crossedEdgeOrdinal = (vertexOrdinal + 2) % 3;
        }

        LogNpcDebug("      Trajectory crosses triangle: crossedEdgeOrdinal=", crossedEdgeOrdinal);

        //
        // Check whether this new edge is floor
        //

        if (IsEdgeFloorToParticle(npcParticle.ConstrainedState->CurrentTriangle, crossedEdgeOrdinal, shipMesh)
            && (isPrimaryParticle || !DoesFloorSeparateFromPrimaryParticle(
                mParticles.GetPosition(npc.DipoleState->SecondaryParticleState.ParticleIndex),
                trajectoryStartAbsolutePosition, // Current (virtual, not yet real) position of this (secondary) particle
                npcParticle.ConstrainedState->CurrentTriangle,
                crossedEdgeOrdinal,
                shipMesh)))
        {
            //
            // Encountered floor
            //

            LogNpcDebug("      Crossed edge is floor - EncounteredFloor");

            return NavigateVertexOutcome::MakeEncounteredFloorOutcome(crossedEdgeOrdinal);
        }

        //
        // Not floor, climb over edge
        //

        // Find opposite triangle
        auto const & oppositeTriangleInfo = shipMesh.GetTriangles().GetOppositeTriangle(npcParticle.ConstrainedState->CurrentTriangle, crossedEdgeOrdinal);
        if (oppositeTriangleInfo.TriangleElementIndex == NoneElementIndex || shipMesh.GetTriangles().IsDeleted(oppositeTriangleInfo.TriangleElementIndex))
        {
            //
            // Become free
            //

            LogNpcDebug("      No opposite triangle found, becoming free - ConvertedToFree");

            TransitionParticleToFreeState(npc, isPrimaryParticle);

            UpdateNpcParticle_Free(
                npcParticle,
                particleStartAbsolutePosition,
                trajectoryEndAbsolutePosition,
                particles,
                gameParameters);

            return NavigateVertexOutcome::MakeConvertedToFreeOutcome();
        }

        //
        // Move to triangle
        //

        LogNpcDebug("      Moving to edge ", oppositeTriangleInfo.EdgeOrdinal, " of opposite triangle ", oppositeTriangleInfo.TriangleElementIndex);

        npcParticle.ConstrainedState->CurrentTriangle = oppositeTriangleInfo.TriangleElementIndex;

        // Calculate new current barycentric coords (wrt opposite triangle - note that we haven't moved)
        bcoords3f newBarycentricCoords; // In new triangle
        newBarycentricCoords[(oppositeTriangleInfo.EdgeOrdinal + 2) % 3] = 0.0f;
        newBarycentricCoords[oppositeTriangleInfo.EdgeOrdinal] = npcParticle.ConstrainedState->CurrentTriangleBarycentricCoords[(crossedEdgeOrdinal + 1) % 3];
        newBarycentricCoords[(oppositeTriangleInfo.EdgeOrdinal + 1) % 3] = npcParticle.ConstrainedState->CurrentTriangleBarycentricCoords[crossedEdgeOrdinal];

        LogNpcDebug("      B-Coords: ", npcParticle.ConstrainedState->CurrentTriangleBarycentricCoords, " -> ", newBarycentricCoords);

        assert(newBarycentricCoords.is_on_edge_or_internal());

        npcParticle.ConstrainedState->CurrentTriangleBarycentricCoords = newBarycentricCoords;

        // New vertex: we know that coord of vertex opposite of crossed edge (i.e. vertex with ordinal crossed_edge+2) is 0.0
        if (newBarycentricCoords[oppositeTriangleInfo.EdgeOrdinal] == 0.0f)
        {
            // Between edge and edge+1
            vertexOrdinal = (oppositeTriangleInfo.EdgeOrdinal + 1) % 3;
        }
        else
        {
            // Between edge and edge-1
            assert(newBarycentricCoords[(oppositeTriangleInfo.EdgeOrdinal + 1) % 3] == 0.0f);
            vertexOrdinal = oppositeTriangleInfo.EdgeOrdinal;
        }

        //
        // Translate target bary coords
        //

        trajectoryEndBarycentricCoords = shipMesh.GetTriangles().ToBarycentricCoordinates(
            trajectoryEndAbsolutePosition,
            npcParticle.ConstrainedState->CurrentTriangle,
            shipMesh.GetPoints(),
            1.0e-07f); // Be strict with roundings - there is a tiny region where the trajectory would oscillate between two triangles (e.g. {0.998306274, 0.00169372594, -3.49245965e-10} -> {-1.61526414e-09, 0.00169372361, 0.998306274})

        LogNpcDebug("      TrajEndB-Coords: ", trajectoryEndBarycentricCoords);

        // At next iteration initial state will be unknown
        isInitialStateUnknown = true;
    }
}

void Npcs::BounceConstrainedNpcParticle(
    StateType & npc,
    bool isPrimaryParticle,
    vec2f const & trajectory,
    vec2f const & bouncePosition,
    vec2f const & bounceEdgeNormal,
    vec2f const meshVelocity,
    float dt,
    NpcParticles & particles,
    GameParameters const & gameParameters) const
{
    auto & npcParticle = isPrimaryParticle ? npc.PrimaryParticleState : npc.DipoleState->SecondaryParticleState;
    assert(npcParticle.ConstrainedState.has_value());

    // Decompose apparent (physical, not walked) particle velocity into normal and tangential

    vec2f const apparentParticleVelocity = trajectory / dt;
    float const apparentParticleVelocityAlongNormal = apparentParticleVelocity.dot(bounceEdgeNormal);
    vec2f const normalVelocity = bounceEdgeNormal * apparentParticleVelocityAlongNormal;
    vec2f const tangentialVelocity = apparentParticleVelocity - normalVelocity;

    // Calculate normal reponse: Vn' = -e*Vn (e = elasticity, [0.0 - 1.0])
    vec2f const normalResponse =
        -normalVelocity
        * particles.GetPhysicalProperties(npcParticle.ParticleIndex).Elasticity
        * gameParameters.ElasticityAdjustment;

    // Calculate tangential response: Vt' = a*Vt (a = (1.0-friction), [0.0 - 1.0])
    vec2f const tangentialResponse =
        tangentialVelocity
        * std::max(0.0f, 1.0f - particles.GetPhysicalProperties(npcParticle.ParticleIndex).KineticFriction * gameParameters.KineticFrictionAdjustment);

    // Given that we've been working in *apparent* space (we've calc'd the collision response to *trajectory* which is apparent displacement),
    // we need to transform velocity to absolute particle velocity
    vec2f const resultantAbsoluteVelocity = (normalResponse + tangentialResponse) - meshVelocity;

    LogNpcDebug("        trajectory=", trajectory, " apparentParticleVelocity=", apparentParticleVelocity, " nr=", normalResponse, " tr=", tangentialResponse, " rr=", resultantAbsoluteVelocity);

    //
    // Set position and velocity
    //

    particles.SetPosition(npcParticle.ParticleIndex, bouncePosition);

    particles.SetVelocity(npcParticle.ParticleIndex, resultantAbsoluteVelocity);
    npcParticle.ConstrainedState->MeshRelativeVelocity = resultantAbsoluteVelocity + meshVelocity;

    //
    // Publish impact
    //

    OnImpact(
        npc,
        isPrimaryParticle,
        normalVelocity,
        bounceEdgeNormal);
}

void Npcs::OnImpact(
    StateType & npc,
    bool isPrimaryParticle,
    vec2f const & impactVector,
    vec2f const & bounceEdgeNormal) const // Pointing outside of triangle
{
    LogNpcDebug("    OnImpact(", impactVector.length(), ", ", bounceEdgeNormal, ")");

    // Human state machine
    if (npc.Kind == NpcKindType::Human)
    {
        OnHumanImpact(
            npc,
            isPrimaryParticle,
            impactVector,
            bounceEdgeNormal);
    }
}

void Npcs::UpdateNpcAnimation(
    StateType & npc,
    bool isPrimaryParticle,
    float currentSimulationTime,
    Ship const & shipMesh,
    GameParameters const & gameParameters)
{
    if (npc.Kind == NpcKindType::Human && isPrimaryParticle) // Take the primary as the only representative of a human
    {
        assert(npc.DipoleState.has_value());
        auto & humanNpcState = npc.KindSpecificState.HumanNpcState;
        using HumanNpcStateType = StateType::KindSpecificStateType::HumanNpcStateType;

        ElementIndex const primaryParticleIndex = npc.PrimaryParticleState.ParticleIndex;
        auto const & primaryContrainedState = npc.PrimaryParticleState.ConstrainedState;
        ElementIndex const secondaryParticleIndex = npc.DipoleState->SecondaryParticleState.ParticleIndex;
        auto const & secondaryConstrainedState = npc.DipoleState->SecondaryParticleState.ConstrainedState;

        float targetRightArmAngle = humanNpcState.RightArmAngle;
        float targetRightArmLengthMultiplier = 1.0f;
        float targetLeftArmAngle = humanNpcState.LeftArmAngle;
        float targetLeftArmLengthMultiplier = 1.0f;

        float targetRightLegAngle = humanNpcState.RightLegAngle;
        float targetRightLegLengthMultiplier = 1.0f;
        float targetLeftLegAngle = humanNpcState.LeftLegAngle;
        float targetLeftLegLengthMultiplier = 1.0f;

        float convergenceRate = 0.0f;

        switch (humanNpcState.CurrentBehavior)
        {
            case HumanNpcStateType::BehaviorType::BeingPlaced:
            {
                float const arg =
                    (currentSimulationTime - humanNpcState.CurrentStateTransitionSimulationTimestamp) * 1.0f
                    + humanNpcState.TotalDistanceTraveledOffEdgeSinceStateTransition * 0.2f;

                float const yArms = std::sin(arg * Pi<float> * 2.0f);
                targetRightArmAngle = Pi<float> / 2.0f + Pi<float> / 2.0f * 0.8f * yArms;
                targetLeftArmAngle = -targetRightArmAngle;

                float const yLegs = std::sin(arg * Pi<float> * 2.0f + npc.RandomNormalizedUniformSeed * Pi<float> * 2.0f);
                targetRightLegAngle = (1.0f + yLegs) / 2.0f * Pi<float> / 2.0f * 0.45f;
                targetLeftLegAngle = -targetRightLegAngle;

                convergenceRate = 0.3f;

                break;
            }

            case HumanNpcStateType::BehaviorType::Constrained_Rising:
            {
                // Arms slightly open
                targetRightArmAngle = HumanNpcStateType::InitialArmAngle;
                targetLeftArmAngle = -HumanNpcStateType::InitialArmAngle;

                // Legs slightly open
                targetRightLegAngle = HumanNpcStateType::InitialLegAngle * 2.0f;
                targetLeftLegAngle = -HumanNpcStateType::InitialLegAngle * 2.0f;

                // Leg and arm that are against floor "help"
                if (primaryContrainedState.has_value() && primaryContrainedState->CurrentVirtualEdgeOrdinal >= 0)
                {
                    vec2f const edgeVector = shipMesh.GetTriangles().GetSubSpringVector(primaryContrainedState->CurrentTriangle, primaryContrainedState->CurrentVirtualEdgeOrdinal, shipMesh.GetPoints());
                    vec2f const head = mParticles.GetPosition(secondaryParticleIndex);
                    vec2f const feet = mParticles.GetPosition(primaryParticleIndex);

                    float humanEdgeAngle = edgeVector.angleCw(head - feet);
                    if (humanEdgeAngle < 0.0f)
                    {
                        humanEdgeAngle += Pi<float>;
                    }

                    // Angle at which left arm is perpendicular to body and touching edge
                    float constexpr MaxAngle = 0.40489178628508342331207292900944f;
                    //static_assert(MaxAngle == std::atan(GameParameters::HumanNpcGeometry::ArmLengthFraction / (1.0f - GameParameters::HumanNpcGeometry::HeadLengthFraction)));

                    if (humanEdgeAngle <= Pi<float> / 2.0f)
                    {
                        // *  0 --> PI/2
                        //  \
                        //   \
                        // ----

                        // Left arm grows up to PI/2 until max angle, then goes down to rest

                        if (humanEdgeAngle <= MaxAngle)
                        {
                            // Quadratically
                            targetLeftArmAngle = -(1.0f - (MaxAngle - humanEdgeAngle) * (MaxAngle - humanEdgeAngle) / (MaxAngle * MaxAngle)) * Pi<float> / 2.0f;
                        }
                        else
                        {
                            targetLeftArmAngle = -Pi<float> / 2.0f - (MaxAngle - humanEdgeAngle) / (MaxAngle - Pi<float> / 2.0f) * (HumanNpcStateType::InitialArmAngle - Pi<float> / 2.0f);
                        }

                        // Right arm helps
                        targetRightArmAngle = targetLeftArmAngle + 0.15f;

                        // Right leg slightly more open
                        targetRightLegAngle *= 0.2f;
                    }
                    else
                    {
                        // <--  *   PI/2 <-- PI
                        //     /
                        //    /
                        //  ------

                        // Right arm grows up to PI/2 until PI - max angle, then goes down to rest

                        if (humanEdgeAngle >= Pi<float> - MaxAngle)
                        {
                            // Quadratically
                            targetRightArmAngle = Pi<float> / 2.0f - (Pi<float> - MaxAngle - humanEdgeAngle) * (Pi<float> -MaxAngle - humanEdgeAngle) / (MaxAngle * MaxAngle) * (Pi<float> / 2.0f);
                        }
                        else
                        {
                            targetRightArmAngle = HumanNpcStateType::InitialArmAngle + (humanEdgeAngle - Pi<float> / 2.0f) / (Pi<float> -MaxAngle - Pi<float> / 2.0f) * (Pi<float> / 2.0f - HumanNpcStateType::InitialArmAngle);
                        }

                        // Left arm helps
                        targetLeftArmAngle = targetRightArmAngle - 0.15f;

                        // Left leg slightly more open
                        targetLeftLegAngle *= 0.2f;
                    }
                }

                convergenceRate = 0.1f;

                break;
            }

            case HumanNpcStateType::BehaviorType::Constrained_Equilibrium:
            {
                // Just small arms angle

                float constexpr ArmsAngle = HumanNpcStateType::InitialArmAngle;

                targetRightArmAngle = ArmsAngle;
                targetLeftArmAngle = -ArmsAngle;

                targetRightLegAngle = 0.0f;
                targetLeftLegAngle = 0.0f;

                convergenceRate = 0.1f;

                break;
            }

            case HumanNpcStateType::BehaviorType::Constrained_Walking:
            {
                //
                // Calculate leg angle based on distance traveled
                //

                float const actualWalkingSpeed = CalculateActualHumanWalkingAbsoluteSpeed(humanNpcState, gameParameters);
                float const actualHalfStepLengthFraction = (GameParameters::HumanNpcGeometry::StepLengthFraction * std::sqrt(actualWalkingSpeed) / 2.0f);
                float const maxLegAngle = std::atan(actualHalfStepLengthFraction / GameParameters::HumanNpcGeometry::LegLengthFraction);

                float const adjustedStandardHumanHeight = GameParameters::HumanNpcGeometry::BodyLength * gameParameters.HumanNpcBodyLengthAdjustment;
                float const stepLength = GameParameters::HumanNpcGeometry::StepLengthFraction * adjustedStandardHumanHeight;
                float const distance =
                    humanNpcState.TotalDistanceTraveledOnEdgeSinceStateTransition
                    + 0.3f * humanNpcState.TotalDistanceTraveledOffEdgeSinceStateTransition;
                float const distanceInTwoSteps = std::fmod(distance + 3.0f * stepLength / 2.0f, stepLength * 2.0f);
                LogNpcDebug("distanceInTwoSteps=", distanceInTwoSteps);

                float const legAngle = std::abs(stepLength - distanceInTwoSteps) / stepLength * 2.0f * maxLegAngle - maxLegAngle;

                targetRightLegAngle = legAngle;
                targetLeftLegAngle = -legAngle;

                // Arms depend on panic
                if (humanNpcState.ResultantPanicLevel < 0.0001f)
                {
                    // Arms aperture depends on speed

                    // At base speed (1m/s): 1.4
                    // Swing more
                    float const apertureMultiplier = 1.4f + (actualWalkingSpeed - 1.0f) * 0.4f;
                    targetRightArmAngle = targetLeftLegAngle * apertureMultiplier;
                }
                else
                {
                    // Arms raised up

                    float const elapsed = currentSimulationTime - humanNpcState.CurrentStateTransitionSimulationTimestamp;
                    float const halfPeriod = 1.0f - 0.6f * std::min(humanNpcState.ResultantPanicLevel, 4.0f) / 4.0f;
                    float const inPeriod = std::fmod(elapsed, halfPeriod * 2.0f);

                    float constexpr MaxAngle = Pi<float> / 2.0f;
                    float const angle = std::abs(halfPeriod - inPeriod) / halfPeriod * 2.0f * MaxAngle - MaxAngle;

                    // PanicMultiplier: p=0.0 => 1.0 p=2.0 => 0.4
                    float const panicMultiplier = 0.4f + 0.6f * (1.0f - std::min(humanNpcState.ResultantPanicLevel, 2.0f) / 2.0f);
                    targetRightArmAngle = Pi<float> -angle * panicMultiplier;
                }
                targetLeftArmAngle = -targetRightArmAngle;

                convergenceRate = 0.25f;

                if (primaryContrainedState.has_value() && primaryContrainedState->CurrentVirtualEdgeOrdinal >= 0)
                {
                    //
                    // We are walking on an edge - make sure feet don't look weird on sloped edges
                    //

                    ElementIndex const edgeElementIndex = shipMesh.GetTriangles().GetSubSprings(primaryContrainedState->CurrentTriangle).SpringIndices[primaryContrainedState->CurrentVirtualEdgeOrdinal];
                    // Note: we do not care if not in CW order
                    vec2f const edg1 = shipMesh.GetSprings().GetEndpointAPosition(edgeElementIndex, shipMesh.GetPoints());
                    vec2f const edg2 = shipMesh.GetSprings().GetEndpointBPosition(edgeElementIndex, shipMesh.GetPoints());
                    vec2f const edgVector = edg2 - edg1;
                    vec2f const edgDir = edgVector.normalise();

                    //
                    // 1. Limit leg angles if on slope
                    //

                    vec2f const headPosition = mParticles.GetPosition(secondaryParticleIndex);
                    vec2f const feetPosition = mParticles.GetPosition(primaryParticleIndex);
                    vec2f const actualBodyVector = feetPosition - headPosition; // From head to feet
                    vec2f const actualBodyDir = actualBodyVector.normalise();

                    float const bodyToVirtualEdgeAlignment = std::abs(edgDir.dot(actualBodyDir.to_perpendicular()));
                    float const angleLimitFactor = bodyToVirtualEdgeAlignment * bodyToVirtualEdgeAlignment * bodyToVirtualEdgeAlignment;
                    targetRightLegAngle *= angleLimitFactor;
                    targetLeftLegAngle *= angleLimitFactor;

                    //
                    // 2. Constrain feet onto edge - i.e. adjust leg lengths
                    //

                    //
                    // leg1 + tl * (leg2 - leg1) = edg1 + te * (edg2 - edg1)
                    // =>
                    // tl = (edg1.y - leg1.y) * (edg2.x - edg1.x) + (leg1.x - edg1.x) * (edg2.y - edg1.y)
                    //      -----------------------------------------------------------------------------
                    //                                  edg X leg
                    //

                    float const adjustedStandardLegLength = GameParameters::HumanNpcGeometry::LegLengthFraction * adjustedStandardHumanHeight;
                    vec2f const crotchPosition = headPosition + actualBodyVector * (GameParameters::HumanNpcGeometry::HeadLengthFraction + GameParameters::HumanNpcGeometry::TorsoLengthFraction);

                    // PERF: the below may be optimized: leg*2 is not needed and leg*1 is crotch => no dependency on legs
                    {
                        vec2f const legrVector = actualBodyDir.rotate(targetRightLegAngle) * adjustedStandardLegLength;
                        float const edgCrossRightLeg = edgVector.cross(legrVector);
                        if (std::abs(edgCrossRightLeg) > 0.0000001f)
                        {
                            vec2f const legr1 = crotchPosition;
                            vec2f const legr2 = crotchPosition + legrVector;

                            targetRightLegLengthMultiplier = ((edg1.y - legr1.y) * (edg2.x - edg1.x) + (legr1.x - edg1.x) * (edg2.y - edg1.y)) / edgCrossRightLeg;
                        }
                    }

                    {
                        vec2f const leglVector = actualBodyDir.rotate(targetLeftLegAngle) * adjustedStandardLegLength;
                        float const edgCrossLeftLeg = edgVector.cross(leglVector);
                        if (std::abs(edgCrossLeftLeg) > 0.0000001f)
                        {
                            vec2f const legl1 = crotchPosition;
                            vec2f const legl2 = crotchPosition + leglVector;

                            targetLeftLegLengthMultiplier = ((edg1.y - legl1.y) * (edg2.x - edg1.x) + (legl1.x - edg1.x) * (edg2.y - edg1.y)) / edgCrossLeftLeg;
                        }
                    }
                }

                break;
            }

            case HumanNpcStateType::BehaviorType::Constrained_Falling:
            {
                // Both arms in direction of face, depending on head velocity in that direction

                vec2f const headPosition = mParticles.GetPosition(secondaryParticleIndex);
                vec2f const feetPosition = mParticles.GetPosition(primaryParticleIndex);
                vec2f const actualBodyVector = feetPosition - headPosition; // From head to feet
                vec2f const actualBodyDir = actualBodyVector.normalise();

                // The extent to which we move arms depends on the avg velocity or head+feet

                vec2f const & headVelocity = npc.DipoleState->SecondaryParticleState.GetApplicableVelocity(mParticles);
                vec2f const & feetVelocity = npc.PrimaryParticleState.GetApplicableVelocity(mParticles);

                float const avgVelocityAlongBodyPerp = (headVelocity + feetVelocity).dot(actualBodyDir.to_perpendicular());
                float const targetDepth = LinearStep(0.0f, 3.0f, std::abs(avgVelocityAlongBodyPerp));

                if (humanNpcState.CurrentFaceDirectionX >= 0.0f)
                {
                    // We want to send arms to the right...
                    // ...but not against face direction
                    if (humanNpcState.CurrentFaceDirectionX >= 0.0f)
                    {
                        targetRightArmAngle = Pi<float> / 2.0f * targetDepth + 0.09f;
                        targetLeftArmAngle = targetRightArmAngle - 0.18f;
                    }
                    else
                    {
                        targetRightArmAngle = 0.0f;
                        targetLeftArmAngle = 0.0f;
                    }
                }
                else
                {
                    // We want to send arms to the left...
                   // ...but not against face direction
                    if (humanNpcState.CurrentFaceDirectionX <= 0.0f)
                    {
                        targetLeftArmAngle = -Pi<float> / 2.0f * targetDepth - 0.09f;
                        targetRightArmAngle = targetLeftArmAngle + 0.18f;
                    }
                    else
                    {
                        targetRightArmAngle = 0.0f;
                        targetLeftArmAngle = 0.0f;
                    }
                }

                convergenceRate = 0.08f * (0.2f + targetDepth);

                // Close legs
                targetRightLegAngle = 0.05f;
                targetLeftLegAngle = -0.05f;

                break;
            }

            case HumanNpcStateType::BehaviorType::Constrained_KnockedOut:
            {
                // Check if both head and feet are on a floor
                bool const isHeadOnEdge = primaryContrainedState.has_value() && IsOnFloorEdge(*primaryContrainedState, shipMesh);
                bool const areFootOnEdge = secondaryConstrainedState.has_value() && IsOnFloorEdge(*secondaryConstrainedState, shipMesh);
                if (isHeadOnEdge && areFootOnEdge)
                {
                    // Arms: +/- PI or 0, depending on where they are now

                    if (humanNpcState.RightArmAngle >= -Pi<float> / 2.0f
                        && humanNpcState.RightArmAngle <= Pi<float> / 2.0f)
                    {
                        targetRightArmAngle = 0.0f;
                    }
                    else
                    {
                        targetRightArmAngle = Pi<float>;
                    }

                    if (humanNpcState.LeftArmAngle >= -Pi<float> / 2.0f
                        && humanNpcState.LeftArmAngle <= Pi<float> / 2.0f)
                    {
                        targetLeftArmAngle = 0.0f;
                    }
                    else
                    {
                        targetLeftArmAngle = Pi<float>;
                    }

                    // Legs: 0

                    targetRightLegAngle = 0.0f;
                    targetLeftLegAngle = 0.0f;

                    convergenceRate = 0.05f; // Quite slow
                }

                break;
            }

            case HumanNpcStateType::BehaviorType::Constrained_Aerial:
            case HumanNpcStateType::BehaviorType::Free_Aerial:
            case HumanNpcStateType::BehaviorType::Free_InWater:
            {
                //
                // Rag doll
                //

                vec2f const headPosition = mParticles.GetPosition(secondaryParticleIndex);
                vec2f const feetPosition = mParticles.GetPosition(primaryParticleIndex);
                vec2f const actualBodyDir = (feetPosition - headPosition).normalise();

                // Arms: always up, unless horizontal

                float const horizontality = std::abs(actualBodyDir.dot(GameParameters::GravityDir));
                float const armAngle = Pi<float> - (Pi<float> / 2.0f) / std::exp(horizontality * 2.2f);
                targetRightArmAngle = armAngle;
                targetLeftArmAngle = -targetRightArmAngle;

                if (humanNpcState.CurrentBehavior != HumanNpcStateType::BehaviorType::Constrained_KnockedOut)
                {
                    // Legs: when arms far from rest, tight; when arms close, at fixed angle

                    float constexpr LegRestAngle = HumanNpcStateType::InitialLegAngle;
                    float const legAngle = LegRestAngle;

                    // Legs inclined in direction opposite of relvel, by an amount proportional to relvel itself
                    vec2f const relativeVelocity = mParticles.GetVelocity(primaryParticleIndex) - mParticles.GetVelocity(secondaryParticleIndex);
                    float const relVelPerpToBody = relativeVelocity.dot(actualBodyDir.to_perpendicular());
                    float const legAngleOffset = -SmoothStep(0.0f, 3.0f, std::abs(relVelPerpToBody)) * (LegRestAngle + 0.3f) * (relVelPerpToBody < 0.0f ? -1.0f : 1.0f);
                    targetRightLegAngle = legAngle - legAngleOffset;
                    targetLeftLegAngle = -legAngle - legAngleOffset;
                }
                else
                {
                    // Leave legs as-is
                }

                convergenceRate = 0.1f;

                break;
            }

            case HumanNpcStateType::BehaviorType::Free_Swimming:
            {
                //
                // 1 period:
                //
                //  _----|         1.0
                // /     \
                // |      \_____|  0.0
                //              |
                //

                float constexpr Period1 = 3.00f;
                float constexpr Period2 = 1.00f;

                float const elapsed = currentSimulationTime - humanNpcState.CurrentStateTransitionSimulationTimestamp;
                float const panicAccelerator = 1.0f + std::min(humanNpcState.ResultantPanicLevel, 2.0f) / 2.0f * 4.0f;

                float const arg =
                    Period1 / 2.0f // Start some-halfway-through to avoid sudden extreme angles
                    + (elapsed) * 2.6f * panicAccelerator
                    + humanNpcState.TotalDistanceTraveledOffEdgeSinceStateTransition * 0.7f;

                float const inPeriod = fmod(arg, (Period1 + Period2));
                float const y = (inPeriod < Period1)
                    ? std::sqrt(inPeriod / Period1)
                    : ((inPeriod - Period1) - Period2) * ((inPeriod - Period1) - Period2) / std::sqrt(Period2);

                // 0: 0, 2: 1, >+ INF: 1
                float const depthDamper = Clamp(mParentWorld.GetOceanSurface().GetDepth(mParticles.GetPosition(secondaryParticleIndex)) / 1.5f, 0.0f, 1.0f);

                // Arms: flapping around PI/2, with amplitude depending on depth
                float constexpr ArmAngleAmplitude = 2.9f; // Half of this on each side of center angle
                float const armCenterAngle = Pi<float> / 2.0f;
                float const armAngle =
                    armCenterAngle
                    + (y * 2.0f - 1.0f) * ArmAngleAmplitude / 2.0f * (depthDamper * 0.75f + 0.25f);
                targetRightArmAngle = armAngle;
                targetLeftArmAngle = -targetRightArmAngle;

                // Legs:flapping around a (small) angle, which becomes even smaller
                // width depth amplitude depending on depth
                float constexpr LegAngleAmplitude = 0.25f * 2.0f; // Half of this on each side of center angle
                float const legCenterAngle = 0.25f * (depthDamper * 0.5f + 0.5f);
                float const legAngle =
                    legCenterAngle
                    + (y * 2.0f - 1.0f) * LegAngleAmplitude / 2.0f * (depthDamper * 0.35f + 0.65f);
                targetRightLegAngle = legAngle;
                targetLeftLegAngle = -targetRightLegAngle;

                // Converge rate depends on how long we've been in this state
                float const MaxConvergenceWait = 3.5f;
                convergenceRate = 0.01f + Clamp(elapsed, 0.0f, MaxConvergenceWait) / MaxConvergenceWait * (0.25f - 0.01f);

                break;
            }
        }

        humanNpcState.RightArmAngle += (targetRightArmAngle - humanNpcState.RightArmAngle) * convergenceRate;
        humanNpcState.RightArmLengthMultiplier += (targetRightArmLengthMultiplier - humanNpcState.RightArmLengthMultiplier) * convergenceRate;
        humanNpcState.LeftArmAngle += (targetLeftArmAngle - humanNpcState.LeftArmAngle) * convergenceRate;
        humanNpcState.LeftArmLengthMultiplier += (targetLeftArmLengthMultiplier - humanNpcState.LeftArmLengthMultiplier) * convergenceRate;

        humanNpcState.RightLegAngle += (targetRightLegAngle - humanNpcState.RightLegAngle) * convergenceRate;
        humanNpcState.RightLegLengthMultiplier += (targetRightLegLengthMultiplier - humanNpcState.RightLegLengthMultiplier) * convergenceRate;
        humanNpcState.LeftLegAngle += (targetLeftLegAngle - humanNpcState.LeftLegAngle) * convergenceRate;
        humanNpcState.LeftLegLengthMultiplier += (targetLeftLegLengthMultiplier - humanNpcState.LeftLegLengthMultiplier) * convergenceRate;
    }
}

}