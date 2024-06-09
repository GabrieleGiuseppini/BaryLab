/***************************************************************************************
 * Original Author:		Gabriele Giuseppini
 * Created:				2023-07-23
 * Copyright:			Gabriele Giuseppini  (https://github.com/GabrieleGiuseppini)
 ***************************************************************************************/
#include "Physics.h"

#include <GameCore/GameMath.h>

#include <array>
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
    auto const & position = mParticles.GetPosition(npc.ParticleMesh.Particles[0].ParticleIndex);

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

    // Particles

    for (size_t p = 0; p < npc.ParticleMesh.Particles.size(); ++p)
    {
        if (p == 0)
        {
            // Primary
            npc.ParticleMesh.Particles[p].ConstrainedState = CalculateParticleConstrainedState(
                mParticles.GetPosition(npc.ParticleMesh.Particles[p].ParticleIndex),
                shipMesh,
                primaryParticleTriangleIndex);
        }
        else
        {
            // Secondaries

            if (!npc.ParticleMesh.Particles[0].ConstrainedState.has_value())
            {
                // When primary is free, also secondary is free (and thus whole NPC)
                if (npc.ParticleMesh.Particles[p].ConstrainedState.has_value())
                {
                    npc.ParticleMesh.Particles[p].ConstrainedState.reset();
                }
            }
            else
            {
                npc.ParticleMesh.Particles[p].ConstrainedState = CalculateParticleConstrainedState(
                    mParticles.GetPosition(npc.ParticleMesh.Particles[p].ParticleIndex),
                    shipMesh,
                    std::nullopt);
            }
        }
    }

    // Regime

    npc.CurrentRegime = CalculateRegime(npc);

    // Kind specific

    switch (npc.Kind)
    {
        case NpcKindType::Furniture:
        {
            // Nop

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

#ifdef BARYLAB_PROBING
    Publish();
#endif
}

void Npcs::TransitionParticleToConstrainedState(
    StateType & npc,
    int npcParticleOrdinal,
    StateType::NpcParticleStateType::ConstrainedStateType constrainedState)
{
    // Transition
    npc.ParticleMesh.Particles[npcParticleOrdinal].ConstrainedState = constrainedState;

    // Regime
    auto const oldRegime = npc.CurrentRegime;
    npc.CurrentRegime = CalculateRegime(npc);
    OnMayBeNpcRegimeChanged(oldRegime, npc);
}

void Npcs::TransitionParticleToFreeState(
    StateType & npc,
    int npcParticleOrdinal)
{
    // Transition
    npc.ParticleMesh.Particles[npcParticleOrdinal].ConstrainedState.reset();
    if (npcParticleOrdinal == 0)
    {
        // When primary is free, also secondaries are free (and thus whole NPC)
        for (size_t p = 1; p < npc.ParticleMesh.Particles.size(); ++p)
        {
            if (npc.ParticleMesh.Particles[p].ConstrainedState.has_value())
            {
                npc.ParticleMesh.Particles[p].ConstrainedState.reset();
            }
        }
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
    // Constrained iff primary is constrained
    assert(npc.ParticleMesh.Particles.size() > 0);
    return npc.ParticleMesh.Particles[0].ConstrainedState.has_value()
        ? StateType::RegimeType::Constrained
        : StateType::RegimeType::Free;
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

            // TODO: move the inner primary check before this loop - for tidyness
            for (auto p = 1; p < npcState->ParticleMesh.Particles.size(); ++p)
            {
                if (!npcState->ParticleMesh.Particles[p].ConstrainedState.has_value() // Secondary is free
                    && npcState->ParticleMesh.Particles[0].ConstrainedState.has_value()) // And primary is constrained
                {
                    assert(mShips[npcState->CurrentShipId].has_value());

                    auto newConstrainedState = CalculateParticleConstrainedState(
                        mParticles.GetPosition(npcState->ParticleMesh.Particles[p].ParticleIndex),
                        mShips[npcState->CurrentShipId]->ShipMesh,
                        std::nullopt);

                    if (newConstrainedState.has_value())
                    {
                        // Make this secondary constrained
                        TransitionParticleToConstrainedState(*npcState, static_cast<int>(p), std::move(*newConstrainedState));
                    }
                }
            }

            // Preliminary Forces

            for (auto p = 0; p < npcState->ParticleMesh.Particles.size(); ++p)
            {
                CalculateNpcParticlePreliminaryForces(
                    *npcState,
                    static_cast<int>(p),
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

            for (auto p = 0; p < npcState->ParticleMesh.Particles.size(); ++p)
            {
                UpdateNpcParticlePhysics(
                    *npcState,
                    static_cast<int>(p),
                    shipMesh,
                    currentSimulationTime,
                    gameParameters);
            }
        }
    }

    //
    // 5. Update behavioral state machines
    //

    LogNpcDebug("----------------------------------");

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
            assert(mShips[npcState->CurrentShipId].has_value());
            auto const & shipMesh = mShips[npcState->CurrentShipId]->ShipMesh;

            for (auto p = 0; p < npcState->ParticleMesh.Particles.size(); ++p)
            {
                UpdateNpcAnimation(
                    *npcState,
                    static_cast<int>(p),
                    currentSimulationTime,
                    shipMesh);
            }
        }
    }
}

void Npcs::UpdateNpcParticlePhysics(
    StateType & npc,
    int npcParticleOrdinal,
    Ship const & shipMesh,
    float currentSimulationTime,
    GameParameters const & gameParameters)
{
    //
    // Here be dragons!
    //

    auto & npcParticle = npc.ParticleMesh.Particles[npcParticleOrdinal];

    LogNpcDebug("----------------------------------");
    LogNpcDebug("  Particle ", npcParticleOrdinal);

    float const particleMass =
        mParticles.GetPhysicalProperties(npcParticle.ParticleIndex).Mass
#ifdef IN_BARYLAB
        * mCurrentMassAdjustment
#endif
        ;
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
            npcParticleOrdinal,
            particleMass,
            gameParameters);

        physicsDeltaPos = mParticles.GetVelocity(npcParticle.ParticleIndex) * dt + (physicalForces / particleMass) * dt * dt;
    }

    if (npc.CurrentRegime == StateType::RegimeType::BeingPlaced
        && ((npc.Kind == NpcKindType::Human && npcParticleOrdinal != 0) || (npc.Kind != NpcKindType::Human)))
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
            && npcParticleOrdinal == 0) // Human is represented by primary particle
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

        // Loop tracing trajectory from TrajectoryStart (== current bary coords in new mesh state) to TrajectoryEnd (== start absolute pos + deltaPos);
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
        // around the inside of a triangle - or a "non-inertial step" - i.e. with the particle pushed against an edge, and
        // thus moving in a non-inertial frame as the mesh acceleration spawns the appearance of apparent forces.
        // - After an inertial iteration we won't enter a non-inertial iteration
        // - A non-inertial iteration might be followed by an inertial one

        // Initialize trajectory start (wrt current triangle) as the absolute pos of the particle as if it just
        // moved with the mesh, staying in its position wrt its triangle; in other words, it's the new theoretical
        // position after just mesh displacement
        vec2f trajectoryStartAbsolutePosition = shipMesh.GetTriangles().FromBarycentricCoordinates(
            npcParticle.ConstrainedState->CurrentBCoords.BCoords,
            npcParticle.ConstrainedState->CurrentBCoords.TriangleElementIndex,
            shipMesh.GetPoints());

        // Calculate mesh velocity for the whole loop as the pure displacement of the triangle containing this particle
        // (end in-triangle position - start absolute position)
        vec2f const meshVelocity = (trajectoryStartAbsolutePosition - particleStartAbsolutePosition) / GameParameters::SimulationTimeStepDuration;

        // Machinery to detect 2- or 3-iteration paths that don't move particle (positional well,
        // aka gravity well)

        std::array<std::optional<AbsoluteTriangleBCoords>, 2> pastBarycentricPositions = { // A queue, insert at zero - initialized with now
            npcParticle.ConstrainedState->CurrentBCoords,
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
        float edgeDistanceTraveledTotal = 0.0f; // To keep track of total distance, and of whether we have moved or not

        // The edge - wrt the current particle's triangle - that we are currently traveling on.
        // Meaningful only during Constrained-NonInertial phases.
        //
        // We begin by not knowing where we are on; when we don't know, we figure this out
        // either by means of vertex navigation (if we are on a vertex) - eventually involving
        // a floor choice in case we are a walking human NPC - or by checking on which floors
        // we are incident (if we are not on a vertex, or navigation of the vertex yields
        // a non-definite edge).
        // From them on, as long as we determine we're surely on an edge (i.e. non-inertial), we remember
        // the edge.

        std::optional<int> currentNonInertialFloorEdgeOrdinal;

        for (float remainingDt = dt; ; )
        {
            assert(remainingDt > 0.0f);

            LogNpcDebug("    ------------------------");
            LogNpcDebug("    New iter: CurrentBCoords=", npcParticle.ConstrainedState->CurrentBCoords, " RemainingDt=", remainingDt);

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
            // Determine whether we are insisting on an edge - i.e. whether we're non-inertial
            //
            // If this is the first time, we determine the edge via vertex navigation (if we
            // are at a vertex) or via floor normals (if we are on an edge, or at a vertex
            // but unsure of the floor as the trajectory might hint to the inside of a triangle).
            //
            // If it's not the first time, then we know we;re at one of these two:
            //  - We cannot determine the edge (because the trajectory might hint to the inside of a triangle)
            //  - We know we are on an edge
            //

            // The edge we determine now - of the triangle @ npcParticle.ConstrainedState->CurrentBCoords; this is a floor
            int nonInertialEdgeOrdinal = -1;

            if (!currentNonInertialFloorEdgeOrdinal.has_value())
            {
                // First time
                if (std::optional<int> const vertexOrdinal = npcParticle.ConstrainedState->CurrentBCoords.BCoords.try_get_vertex();
                    vertexOrdinal.has_value())
                {
                    //
                    // We are on a vertex (*two* edges)
                    //
                    // There are two conditions that bring us here:
                    //  1. After a bounce or an impact continuation
                    //  2. Initial (e.g. placed by chance on a vertex)
                    //

                    LogNpcDebug("    On a vertex: ", *vertexOrdinal, " of triangle ", npcParticle.ConstrainedState->CurrentBCoords.TriangleElementIndex);

                    // Check what comes next

                    bcoords3f const trajectoryEndBarycentricCoords = shipMesh.GetTriangles().ToBarycentricCoordinates(
                        trajectoryEndAbsolutePosition,
                        npcParticle.ConstrainedState->CurrentBCoords.TriangleElementIndex,
                        shipMesh.GetPoints());

                    auto const outcome = NavigateVertex(
                        npc,
                        npcParticleOrdinal,
                        npcParticle.ConstrainedState->CurrentVirtualFloor,
                        *vertexOrdinal,
                        trajectory,
                        trajectoryEndAbsolutePosition,
                        trajectoryEndBarycentricCoords,
                        shipMesh,
                        mParticles);

                    switch (outcome.Type)
                    {
                        case NavigateVertexOutcome::OutcomeType::BecomeFree:
                        {
                            // Transition to free immediately

                            TransitionParticleToFreeState(npc, npcParticleOrdinal);

                            UpdateNpcParticle_Free(
                                npcParticle,
                                particleStartAbsolutePosition,
                                trajectoryEndAbsolutePosition,
                                mParticles,
                                gameParameters);

                            remainingDt = 0.0f;

                            break;
                        }

                        case NavigateVertexOutcome::OutcomeType::ContinueAlongFloor:
                        {
                            //
                            // This either tells us that our trajectory points us towards this floor (if not walking),
                            // or that we're walking and decided for this floor: in either case we are non-inertial on this floor
                            //

                            assert(outcome.FloorEdgeOrdinal >= 0);
                            nonInertialEdgeOrdinal = outcome.FloorEdgeOrdinal;

                            // Move to NavigationOutcome
                            npcParticle.ConstrainedState->CurrentBCoords = outcome.TriangleBCoords;

                            break;
                        }

                        case NavigateVertexOutcome::OutcomeType::ContinueToInterior:
                        {
                            //
                            // This tells us that the trajectory appears to tell that we're
                            // moving into a triangle or along an edge of it, but this could
                            // be a numerical illusion.
                            // Need to use floor normals to see if we are
                            // incident to a triangle
                            //

                            assert(outcome.FloorEdgeOrdinal == -1);
                            nonInertialEdgeOrdinal = -1; // Determine via floor-normal checks

                            // Move to NavigationOutcome
                            npcParticle.ConstrainedState->CurrentBCoords = outcome.TriangleBCoords;

#ifdef _DEBUG

                            // Assert dot-with-normal is negative or very small for all edges we touch
                            // (after all the numerical illusion is just lack of precision)
                            for (int vi = 0; vi < 3; ++vi)
                            {
                                if (npcParticle.ConstrainedState->CurrentBCoords.BCoords[vi] == 0.0f)
                                {
                                    int const edgeOrdinal = (vi + 1) % 3;

                                    vec2f const edgeVector = shipMesh.GetTriangles().GetSubSpringVector(
                                        npcParticle.ConstrainedState->CurrentBCoords.TriangleElementIndex,
                                        edgeOrdinal,
                                        shipMesh.GetPoints());
                                    vec2f const edgeDir = edgeVector.normalise();
                                    vec2f const edgeNormal = edgeDir.to_perpendicular(); // Points outside of triangle (i.e. towards floor)
                                    assert(trajectory.dot(edgeNormal) < 0.01f);
                                }
                            }
#endif
                            break;
                        }

                        case NavigateVertexOutcome::OutcomeType::ImpactOnFloor:
                        {
                            // Let the subsequent non-inertial ray tracing take care of this

                            assert(outcome.FloorEdgeOrdinal >= 0);
                            nonInertialEdgeOrdinal = outcome.FloorEdgeOrdinal;

                            // Do not move to NavigationOutcome, we'll get there via ray tracing

                            break;
                        }
                    }

                    // Follow-up free conversion
                    if (remainingDt == 0.0f)
                    {
                        break;
                    }
                }
            }
            else
            {
                LogNpcDebug("    We remember from earlier that we are on a floor edge/inside triangle: ", *currentNonInertialFloorEdgeOrdinal, " of triangle ", npcParticle.ConstrainedState->CurrentBCoords.TriangleElementIndex);

                nonInertialEdgeOrdinal = *currentNonInertialFloorEdgeOrdinal;
            }

            if (nonInertialEdgeOrdinal == -1)
            {
                LogNpcDebug("    Determining floor edge with floor normals");

                // Find edge that trajectory is against the most
                float bestHeadsOnNess = std::numeric_limits<float>::lowest();
                for (int vi = 0; vi < 3; ++vi)
                {
                    if (npcParticle.ConstrainedState->CurrentBCoords.BCoords[vi] == 0.0f)
                    {
                        // We are on this edge

                        int const edgeOrdinal = (vi + 1) % 3;

                        LogNpcDebug("    On an edge: ", edgeOrdinal, " of triangle ", npcParticle.ConstrainedState->CurrentBCoords.TriangleElementIndex);

                        // Check if this is really a floor to this particle
                        if (IsEdgeFloorToParticle(npcParticle.ConstrainedState->CurrentBCoords.TriangleElementIndex, edgeOrdinal, npc, npcParticleOrdinal, mParticles, shipMesh))
                        {
                            // The edge is a floor

                            LogNpcDebug("      The edge is a floor");

                            // Check now whether we're moving *against* the floor

                            vec2f const edgeVector = shipMesh.GetTriangles().GetSubSpringVector(
                                npcParticle.ConstrainedState->CurrentBCoords.TriangleElementIndex,
                                edgeOrdinal,
                                shipMesh.GetPoints());
                            vec2f const edgeNormal = edgeVector.to_perpendicular(); // Points outside of triangle (i.e. towards floor)
                            float const headsOnNess = trajectory.dot(edgeNormal);
                            if (trajectory.dot(edgeNormal) >= 0.0f) // If 0, no normal force - hence no friction; however we want to take this codepath anyways for consistency
                            {
                                //
                                // We are insisting against this edge floor
                                //

                                LogNpcDebug("      We are moving against this floor (alignment=", headsOnNess, ")");

                                if (headsOnNess > bestHeadsOnNess)
                                {
                                    // Winner
                                    nonInertialEdgeOrdinal = edgeOrdinal;
                                    bestHeadsOnNess = headsOnNess;
                                }
                            }
                            else
                            {
                                LogNpcDebug("      We are not moving against this floor (alignment=", headsOnNess, ")");
                            }
                        }
                        else
                        {
                            LogNpcDebug("      The edge is not a floor");
                        }
                    }
                }
            }

            if (nonInertialEdgeOrdinal >= 0)
            {
                //
                // Case 1: Non-inertial: on edge and moving against (or along) it, pushed by it
                //

                LogNpcDebug("    ConstrainedNonInertial: triangle=", npcParticle.ConstrainedState->CurrentBCoords.TriangleElementIndex, " nonInertialEdgeOrdinal=", nonInertialEdgeOrdinal, " bCoords=", npcParticle.ConstrainedState->CurrentBCoords.BCoords, " trajectory=", trajectory);
                LogNpcDebug("    StartPosition=", mParticles.GetPosition(npcParticle.ParticleIndex), " StartVelocity=", mParticles.GetVelocity(npcParticle.ParticleIndex), " MeshVelocity=", meshVelocity, " StartMRVelocity=", npcParticle.ConstrainedState->MeshRelativeVelocity);

                // Set now our current edge-ness for the rest of this simulation step, which is
                // mostly for animation purposes.
                //
                // In fact, the subsequent UpdateNpcParticle_ConstrainedNonInertial call has 5 outcomes:
                //  - We end up inside the triangle - along this same edge;
                //  - We travel up to a vertex, and then:
                //      - End towards interior of triangle: in this case we might have reached a new triangle,
                //        but we'll go through this again and become inertial (resetting this edge-ness);
                //      - Ready to continue on a new edge: we could update this edge-ness to the new edge, but we'll go
                //        through this again, and thus we save the effort;
                //      - Become free: no more constrained state anyways;
                //      - Hit floor: in this case we _do_ want the current virtual floor to be the one _before_ the bounce.

                npcParticle.ConstrainedState->CurrentVirtualFloor.emplace(
                    npcParticle.ConstrainedState->CurrentBCoords.TriangleElementIndex,
                    nonInertialEdgeOrdinal);

                //
                // We're moving against the floor, hence we are in a non-inertial frame...
                // ...take friction into account and flaten trajectory
                //

                //
                // Calculate magnitude of flattened trajectory - i.e. component of trajectory
                // along (i.e. tangent) edge, positive when in the direction of the edge
                //

                vec2f const edgeDir = shipMesh.GetTriangles().GetSubSpringVector(
                    npcParticle.ConstrainedState->CurrentBCoords.TriangleElementIndex,
                    nonInertialEdgeOrdinal,
                    shipMesh.GetPoints()).normalise();

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
                    && npcParticleOrdinal == 0)
                {
                    //
                    // Walking displacement projected along edge - what's needed to reach total desired walking
                    // displacement, together with physical
                    // - Never reduces physical
                    // - Never adds to physical so much as to cause resultant to be faster than walking speed
                    //

                    vec2f const idealWalkDir = vec2f(npc.KindSpecificState.HumanNpcState.CurrentFaceDirectionX, 0.0f);
                    assert(idealWalkDir.length() == 1.0f);

                    float const idealWalkMagnitude = CalculateActualHumanWalkingAbsoluteSpeed(npc.KindSpecificState.HumanNpcState) * remainingDt;

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
                    //
                    // Note: walkDir.y is sin(slope angle between horiz and dir)
                    float constexpr NeighborhoodWidth = 1.0f - GameParameters::MaxHumanNpcWalkSinSlope;
                    float constexpr ResistanceSinSlopeStart = GameParameters::MaxHumanNpcWalkSinSlope - NeighborhoodWidth / 2.0f;
                    if (walkDir.y >= ResistanceSinSlopeStart) // walkDir.y is component along vertical, pointing up
                    {
                        float const y2 = (walkDir.y - ResistanceSinSlopeStart) / NeighborhoodWidth;
                        float const gravityResistance = std::max(1.0f - y2, 0.0f);

                        LogNpcDebug("        gravityResistance=", gravityResistance, " (walkDir.y=", walkDir.y, ")");

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
                    npcParticle.ConstrainedState->CurrentBCoords.TriangleElementIndex,
                    shipMesh.GetPoints());

                //
                // Due to numerical slack, ensure target barycentric coords are still along edge
                //

                int const edgeVertexOrdinal = (nonInertialEdgeOrdinal + 2) % 3;
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

                auto const [edgeTraveledActual, doStop, newFloorEdgeOrdinal] = UpdateNpcParticle_ConstrainedNonInertial(
                    npc,
                    npcParticleOrdinal,
                    nonInertialEdgeOrdinal,
                    edgeDir,
                    particleStartAbsolutePosition,
                    trajectoryStartAbsolutePosition,
                    flattenedTrajectoryEndAbsolutePosition,
                    flattenedTrajectoryEndBarycentricCoords,
                    flattenedTrajectory,
                    adjustedEdgeTraveledPlanned,
                    (edgeDistanceTraveledTotal > 0.0f), // hasMovedInStep
                    meshVelocity,
                    remainingDt,
                    shipMesh,
                    mParticles,
                    currentSimulationTime,
                    gameParameters);

                if (npcParticle.ConstrainedState.has_value())
                {
                    LogNpcDebug("    EndBCoords=", npcParticle.ConstrainedState->CurrentBCoords);
                }
                else
                {
                    LogNpcDebug("    Became free");
                }

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
                        if (pastBarycentricPositions[0] == npcParticle.ConstrainedState->CurrentBCoords
                            || pastBarycentricPositions[1] == npcParticle.ConstrainedState->CurrentBCoords)
                        {
                            //
                            // Well - stop here
                            //

                            LogNpcDebug("    Detected well - stopping here");

                            // Update particle's physics, considering that we are in a well and thus still (wrt mesh)

                            vec2f const particleEndAbsolutePosition = shipMesh.GetTriangles().FromBarycentricCoordinates(
                                npcParticle.ConstrainedState->CurrentBCoords.BCoords,
                                npcParticle.ConstrainedState->CurrentBCoords.TriangleElementIndex,
                                shipMesh.GetPoints());

                            mParticles.SetPosition(npcParticle.ParticleIndex, particleEndAbsolutePosition);

                            // No (relative) velocity (so just mesh velocity, and no global damping)
                            mParticles.SetVelocity(npcParticle.ParticleIndex, meshVelocity);
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
                        static_assert(pastBarycentricPositions.size() == 2);
                        pastBarycentricPositions[0].reset();
                        pastBarycentricPositions[1].reset();
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
                    && npcParticleOrdinal == 0) // Human is represented by primary particle
                {
                    // Note: does not include eventual distance traveled after becoming free; fine because we will transition and wipe out total traveled
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
                    pastBarycentricPositions[0] = npcParticle.ConstrainedState->CurrentBCoords;
                }

                // Update current floor edge we're walking on, or inside triangle
                currentNonInertialFloorEdgeOrdinal = newFloorEdgeOrdinal;

                if (npcParticle.ConstrainedState.has_value())
                {
                    LogNpcDebug("    EndPosition=", mParticles.GetPosition(npcParticle.ParticleIndex), " EndVelocity=", mParticles.GetVelocity(npcParticle.ParticleIndex), " EndMRVelocity=", npcParticle.ConstrainedState->MeshRelativeVelocity);
                }
                else
                {
                    LogNpcDebug("    EndPosition=", mParticles.GetPosition(npcParticle.ParticleIndex), " EndVelocity=", mParticles.GetVelocity(npcParticle.ParticleIndex));
                }
            }
            else
            {
                //
                // Case 2: Inertial: not on edge or on edge but not moving against (nor along) it
                //

                LogNpcDebug("    ConstrainedInertial: CurrentBCoords=", npcParticle.ConstrainedState->CurrentBCoords, " physicsDeltaPos=", physicsDeltaPos);
                LogNpcDebug("    StartPosition=", mParticles.GetPosition(npcParticle.ParticleIndex), " StartVelocity=", mParticles.GetVelocity(npcParticle.ParticleIndex), " MeshVelocity=", meshVelocity, " StartMRVelocity=", npcParticle.ConstrainedState->MeshRelativeVelocity);

                // We are not on an edge floor (nor we'll got on it now)
                npcParticle.ConstrainedState->CurrentVirtualFloor.reset();
                currentNonInertialFloorEdgeOrdinal.reset();

                //
                // Calculate target barycentric coords
                //

                bcoords3f const trajectoryEndBarycentricCoords = shipMesh.GetTriangles().ToBarycentricCoordinates(
                    trajectoryEndAbsolutePosition,
                    npcParticle.ConstrainedState->CurrentBCoords.TriangleElementIndex,
                    shipMesh.GetPoints());

                //
                // Move towards target bary coords
                //

                float totalTraveled = UpdateNpcParticle_ConstrainedInertial(
                    npc,
                    npcParticleOrdinal,
                    particleStartAbsolutePosition,
                    trajectoryStartAbsolutePosition, // segmentTrajectoryStartAbsolutePosition
                    trajectoryEndAbsolutePosition,
                    trajectoryEndBarycentricCoords,
                    (edgeDistanceTraveledTotal > 0.0f),
                    meshVelocity,
                    remainingDt,
                    shipMesh,
                    mParticles,
                    currentSimulationTime,
                    gameParameters);

                if (npcParticle.ConstrainedState.has_value())
                {
                    LogNpcDebug("    EndBCoords=", npcParticle.ConstrainedState->CurrentBCoords);
                }
                else
                {
                    LogNpcDebug("    Became free");
                }

                // Update total traveled
                if (npc.Kind == NpcKindType::Human
                    && npcParticleOrdinal == 0) // Human is represented by primary particle
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
                npcParticle.ConstrainedState->CurrentBCoords.BCoords,
                npcParticle.ConstrainedState->CurrentBCoords.TriangleElementIndex,
                shipMesh.GetPoints());
        }
    }

#ifdef BARYLAB_PROBING
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
    int npcParticleOrdinal,
    GameParameters const & gameParameters)
{
    auto & npcParticle = npc.ParticleMesh.Particles[npcParticleOrdinal];

    //
    // Calculate world forces
    //

    float const particleMass =
        mParticles.GetPhysicalProperties(npcParticle.ParticleIndex).Mass
#ifdef IN_BARYLAB
        * mCurrentMassAdjustment
#endif
        ;

    // 1. World forces - gravity

    vec2f preliminaryForces =
        GameParameters::Gravity
#ifdef IN_BARYLAB
        * mCurrentGravityAdjustment
#endif
        * particleMass;

    if (!npcParticle.ConstrainedState.has_value() && npc.CurrentRegime != StateType::RegimeType::BeingPlaced)
    {
        // Check whether we are underwater

        float constexpr BuoyancyInterfaceWidth = 0.4f;

        vec2f testParticlePosition = mParticles.GetPosition(npcParticle.ParticleIndex);
        if (npc.Kind == NpcKindType::Human && npcParticleOrdinal > 0)
        {
            // Head - use an offset
            testParticlePosition += (mParticles.GetPosition(npc.ParticleMesh.Particles[0].ParticleIndex) - testParticlePosition) * BuoyancyInterfaceWidth * 2.0f / 3.0f;
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

    for (auto const & connectedSpring : npcParticle.ConnectedSprings)
    {
        float const dt = GameParameters::SimulationTimeStepDuration;

        vec2f const springDisplacement = mParticles.GetPosition(connectedSpring.OtherEndpointIndex) - mParticles.GetPosition(npcParticle.ParticleIndex); // Towards other
        float const springDisplacementLength = springDisplacement.length();
        vec2f const springDir = springDisplacement.normalise_approx(springDisplacementLength);

        //
        // 3a. Hooke's law
        //

        auto const & spring = npc.ParticleMesh.Springs[connectedSpring.SpringOrdinal];

        // Calculate spring force on this particle
        float const fSpring =
            (springDisplacementLength - spring.DipoleLength)
            * spring.SpringStiffnessCoefficient;

        //
        // 3b. Damper forces
        //
        // Damp the velocities of each endpoint pair, as if the points were also connected by a damper
        // along the same direction as the spring
        //

        // Calculate damp force on this particle
        vec2f const relVelocity = mParticles.GetVelocity(connectedSpring.OtherEndpointIndex) - mParticles.GetVelocity(npcParticle.ParticleIndex);
        float const fDamp =
            relVelocity.dot(springDir)
            * spring.SpringDampingCoefficient;

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
    int npcParticleOrdinal,
    float particleMass,
    GameParameters const & gameParameters) const
{
    auto & npcParticle = npc.ParticleMesh.Particles[npcParticleOrdinal];

    vec2f definitiveForces = mParticles.GetPreliminaryForces(npcParticle.ParticleIndex);

    //
    // Human Equlibrium Torque
    //

    if (npc.Kind == NpcKindType::Human
        && npc.KindSpecificState.HumanNpcState.EquilibriumTorque != 0.0f
        && npcParticleOrdinal > 0)
    {
        ElementIndex const primaryParticleIndex = npc.ParticleMesh.Particles[0].ParticleIndex;
        ElementIndex const secondaryParticleIndex = npcParticle.ParticleIndex;

        vec2f const feetPosition = mParticles.GetPosition(primaryParticleIndex);

        // Given that we apply torque onto the secondary particle *after* the primary has been simulated
        // (so that we take into account the primary's new position), and thus we see the primary where it
        // is at the end of the step - possibly far away if mesh velocity is high, we want to _predict_
        // where the secondary will be by its own velocity

        vec2f const headPredictedPosition =
            mParticles.GetPosition(secondaryParticleIndex)
            + mParticles.GetVelocity(secondaryParticleIndex) * GameParameters::SimulationTimeStepDuration;

        vec2f const humanDir = (headPredictedPosition - feetPosition).normalise_approx();

        // Calculate radial force direction
        vec2f radialDir = humanDir.to_perpendicular(); // CCW
        if (humanDir.x < 0.0f)
        {
            // Head is to the left of the vertical
            radialDir *= -1.0f;
        }

        //
        // First component: Hookean force proportional to length of arc to destination, directed towards the ideal head position
        //  - But we approximate the arc with the chord, i.e.the distance between source and destination
        //

        vec2f const idealHeadPosition = feetPosition + vec2f(0.0f, npc.KindSpecificState.HumanNpcState.Height * mCurrentHumanNpcBodyLengthAdjustment);

        float const stiffnessCoefficient =
            gameParameters.HumanNpcEquilibriumTorqueStiffnessCoefficient
            + std::min(npc.KindSpecificState.HumanNpcState.ResultantPanicLevel, 1.0f) * 0.0005f;

        float const force1Magnitude =
            (idealHeadPosition - headPredictedPosition).length()
            * stiffnessCoefficient;

        //
        // Second component: damp force proportional to component of relative velocity that is orthogonal to human vector, opposite that velocity
        //

        vec2f const relativeVelocity = mParticles.GetVelocity(secondaryParticleIndex) - mParticles.GetVelocity(primaryParticleIndex);
        float const orthoRelativeVelocity = relativeVelocity.dot(radialDir);

        float const dampCoefficient = gameParameters.HumanNpcEquilibriumTorqueDampingCoefficient;

        float const force2Magnitude =
            -orthoRelativeVelocity
            * dampCoefficient;

        //
        // Combine
        //

        vec2f const equilibriumTorqueForce =
            radialDir
            * (force1Magnitude + force2Magnitude)
            / mCurrentHumanNpcBodyLengthAdjustment // Note: we divide by human length adjustment to maintain torque independent from lever length
            * particleMass / (GameParameters::SimulationTimeStepDuration * GameParameters::SimulationTimeStepDuration);

        definitiveForces += equilibriumTorqueForce;
    }

    return definitiveForces;
}

void Npcs::RecalculateSpringForceParameters()
{
    // Visit all springs

    for (auto & npc : mStateBuffer)
    {
        if (npc.has_value())
        {
            for (auto & spring : npc->ParticleMesh.Springs)
            {
                RecalculateSpringForceParameters(spring);
            }
        }
    }
}

void Npcs::RecalculateSpringForceParameters(StateType::NpcSpringStateType & spring) const
{
    spring.SpringStiffnessCoefficient =
        mCurrentSpringReductionFraction
        * spring.MassFactor
#ifdef IN_BARYLAB
        * mCurrentMassAdjustment
#endif
        / (GameParameters::SimulationTimeStepDuration * GameParameters::SimulationTimeStepDuration);

    spring.SpringDampingCoefficient =
        mCurrentSpringDampingCoefficient
        * spring.MassFactor
#ifdef IN_BARYLAB
        * mCurrentMassAdjustment
#endif
        / GameParameters::SimulationTimeStepDuration;
}

void Npcs::RecalculateHumanNpcDipoleLengths()
{
    for (auto & state : mStateBuffer)
    {
        if (state.has_value())
        {
            if (state->Kind == NpcKindType::Human)
            {
                assert(state->ParticleMesh.Springs.size() == 1);
                state->ParticleMesh.Springs[0].DipoleLength = CalculateHumanNpcDipoleLength(state->KindSpecificState.HumanNpcState.Height);
            }
        }
    }
}

float Npcs::CalculateHumanNpcDipoleLength(float baseHeight) const
{
    return baseHeight * mCurrentHumanNpcBodyLengthAdjustment;
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

Npcs::ConstrainedNonInertialOutcome Npcs::UpdateNpcParticle_ConstrainedNonInertial(
    StateType & npc,
    int npcParticleOrdinal,
    int edgeOrdinal,
    vec2f const & edgeDir,
    vec2f const & particleStartAbsolutePosition,
    vec2f const & trajectoryStartAbsolutePosition,
    vec2f const & flattenedTrajectoryEndAbsolutePosition,
    bcoords3f flattenedTrajectoryEndBarycentricCoords,
    vec2f const & flattenedTrajectory,
    float edgeTraveledPlanned,
    bool hasMovedInStep,
    vec2f const meshVelocity,
    float dt,
    Ship const & shipMesh,
    NpcParticles & particles,
    float currentSimulationTime,
    GameParameters const & gameParameters)
{
    auto & npcParticle = npc.ParticleMesh.Particles[npcParticleOrdinal];
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

        npcParticleConstrainedState.CurrentBCoords.BCoords = flattenedTrajectoryEndBarycentricCoords;
        particles.SetPosition(npcParticle.ParticleIndex, flattenedTrajectoryEndAbsolutePosition);

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

        assert(std::abs((flattenedTrajectoryEndAbsolutePosition - trajectoryStartAbsolutePosition).dot(edgeDir) - edgeTraveledPlanned) < 0.001f);

        vec2f const relativeVelocity =
            edgeDir
            * edgeTraveledPlanned
            / dt;

        vec2f const absoluteVelocity =
            // Do not damp velocity if we're trying to maintain equilibrium
            relativeVelocity * ((npc.Kind != NpcKindType::Human || npc.KindSpecificState.HumanNpcState.EquilibriumTorque == 0.0f) ? (1.0f - gameParameters.GlobalDamping) : 1.0f)
            + meshVelocity;

        particles.SetVelocity(npcParticle.ParticleIndex, absoluteVelocity);
        npcParticleConstrainedState.MeshRelativeVelocity = relativeVelocity;

        LogNpcDebug("        edgeTraveleded (==planned)=", edgeTraveledPlanned, " absoluteVelocity=", particles.GetVelocity(npcParticle.ParticleIndex));

        // Complete
        return {
            edgeTraveledPlanned,    // We moved by how much we planned
            true,                   // Stop here
            -1 };                   // No edge (but we're stopping anyways)
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

    assert(npcParticleConstrainedState.CurrentBCoords.BCoords[(edgeOrdinal + 2) % 3] == 0.0f);

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
    // Move to intersection vertex
    //
    // Note that except for the very first iteration, any other iteration will travel
    // zero distance at this moment
    //

    bcoords3f intersectionBarycentricCoords = bcoords3f::zero();
    intersectionBarycentricCoords[intersectionVertexOrdinal] = 1.0f;

    LogNpcDebug("      Moving to intersection vertex ", intersectionVertexOrdinal, ": ", intersectionBarycentricCoords);

    // Calculate (signed) edge traveled
    float edgeTraveled;
    if (intersectionBarycentricCoords == npcParticleConstrainedState.CurrentBCoords.BCoords)
    {
        // We haven't moved - ensure pure zero
        edgeTraveled = 0.0f;
    }
    else
    {
        vec2f const intersectionAbsolutePosition = shipMesh.GetPoints().GetPosition(shipMesh.GetTriangles().GetPointIndices(npcParticleConstrainedState.CurrentBCoords.TriangleElementIndex)[intersectionVertexOrdinal]);
        assert(intersectionAbsolutePosition == shipMesh.GetTriangles().FromBarycentricCoordinates(
            intersectionBarycentricCoords,
            npcParticleConstrainedState.CurrentBCoords.TriangleElementIndex,
            shipMesh.GetPoints()));

        edgeTraveled = (intersectionAbsolutePosition - trajectoryStartAbsolutePosition).dot(edgeDir);
    }

    LogNpcDebug("        edgeTraveled=", edgeTraveled);

    // Move
    npcParticleConstrainedState.CurrentBCoords.BCoords = intersectionBarycentricCoords;

    //
    // Navigate this vertex now
    //

    auto const navigationOutcome = NavigateVertex(
        npc,
        npcParticleOrdinal,
        TriangleAndEdge(npcParticleConstrainedState.CurrentBCoords.TriangleElementIndex, edgeOrdinal),
        intersectionVertexOrdinal,
        flattenedTrajectory,
        flattenedTrajectoryEndAbsolutePosition,
        flattenedTrajectoryEndBarycentricCoords,
        shipMesh,
        particles);

    switch (navigationOutcome.Type)
    {
        case NavigateVertexOutcome::OutcomeType::BecomeFree:
        {
            // Become free

            TransitionParticleToFreeState(npc, npcParticleOrdinal);

            UpdateNpcParticle_Free(
                npcParticle,
                particleStartAbsolutePosition,
                flattenedTrajectoryEndAbsolutePosition,
                particles,
                gameParameters);

            // Terminate
            return {
                edgeTraveled,
                true,
                std::nullopt };
        }

        case NavigateVertexOutcome::OutcomeType::ContinueAlongFloor:
        {
            // Impact continuation, continue

            // Move to NavigationOutcome
            npcParticleConstrainedState.CurrentBCoords = navigationOutcome.TriangleBCoords;

            // Continue
            return {
                edgeTraveled,
                false,
                navigationOutcome.FloorEdgeOrdinal }; // This is the edge we want to be on, we may bypass the next NavigateVertex
        }

        case NavigateVertexOutcome::OutcomeType::ContinueToInterior:
        {
            // Continue

            // Move to NavigationOutcome
            npcParticleConstrainedState.CurrentBCoords = navigationOutcome.TriangleBCoords;

            // Continue
            return {
                edgeTraveled,
                false,
                -1 }; // Our next edge is...no edge
        }

        case NavigateVertexOutcome::OutcomeType::ImpactOnFloor:
        {
            // Bounce

            vec2f const floorEdgeDir =
                shipMesh.GetTriangles().GetSubSpringVector(
                    navigationOutcome.TriangleBCoords.TriangleElementIndex,
                    navigationOutcome.FloorEdgeOrdinal,
                    shipMesh.GetPoints())
                .normalise();
            vec2f const floorEdgeNormal = floorEdgeDir.to_perpendicular();

            vec2f const bounceAbsolutePosition = shipMesh.GetTriangles().FromBarycentricCoordinates(
                navigationOutcome.TriangleBCoords.BCoords,
                navigationOutcome.TriangleBCoords.TriangleElementIndex,
                shipMesh.GetPoints());

            BounceConstrainedNpcParticle(
                npc,
                npcParticleOrdinal,
                flattenedTrajectory,
                hasMovedInStep || (edgeTraveled != 0.0f),
                bounceAbsolutePosition,
                floorEdgeNormal,
                meshVelocity,
                dt,
                particles,
                currentSimulationTime,
                gameParameters);

            // Move to bounce position - so that we can eventually make a fresh new choice
            // now that we have acquired a velocity in a different direction
            npcParticleConstrainedState.CurrentBCoords = navigationOutcome.TriangleBCoords;

            // Terminate
            return {
                edgeTraveled,
                true,
                std::nullopt }; // Certainly the edge we bounce on is not the one we want to be on - and moreover, we are done with the loop now
        }
    }

    assert(false);
    return { edgeTraveled, false, -1 };
}

float Npcs::UpdateNpcParticle_ConstrainedInertial(
    StateType & npc,
    int npcParticleOrdinal,
    vec2f const & particleStartAbsolutePosition, // Since beginning of whole time quantum, not just this step
    vec2f const & segmentTrajectoryStartAbsolutePosition,
    vec2f const & segmentTrajectoryEndAbsolutePosition,
    bcoords3f segmentTrajectoryEndBarycentricCoords, // In current triangle; mutable
    bool hasMovedInStep,
    vec2f const meshVelocity,
    float segmentDt,
    Ship const & shipMesh,
    NpcParticles & particles,
    float currentSimulationTime,
    GameParameters const & gameParameters)
{
    auto & npcParticle = npc.ParticleMesh.Particles[npcParticleOrdinal];
    assert(npcParticle.ConstrainedState.has_value());
    auto & npcParticleConstrainedState = *npcParticle.ConstrainedState;

    //
    // Ray-trace along the specified trajectory, ending only at one of the following three conditions:
    // 1. Reached destination: terminate
    // 2. Becoming free: do free movement and terminate
    // 3. Impact with bounce: impart bounce velocity and terminate
    //

    for (int iIter = 0; ; ++iIter)
    {
        assert(iIter < 8); // Detect and debug-break on infinite loops

        assert(npcParticleConstrainedState.CurrentBCoords.BCoords.is_on_edge_or_internal());

        LogNpcDebug("    SegmentTrace ", iIter);
        LogNpcDebug("      triangle=", npcParticleConstrainedState.CurrentBCoords.TriangleElementIndex, " bCoords=", npcParticleConstrainedState.CurrentBCoords.BCoords,
            " segmentTrajEndBCoords=", segmentTrajectoryEndBarycentricCoords);

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
            npcParticleConstrainedState.CurrentBCoords.BCoords = segmentTrajectoryEndBarycentricCoords;
            particles.SetPosition(npcParticle.ParticleIndex, segmentTrajectoryEndAbsolutePosition);

            // Use whole time quantum for velocity, as particleStartAbsolutePosition is fixed at t0
            vec2f const totalAbsoluteTraveledVector = segmentTrajectoryEndAbsolutePosition - particleStartAbsolutePosition;
            vec2f const relativeVelocity = totalAbsoluteTraveledVector / GameParameters::SimulationTimeStepDuration - meshVelocity;

            // Do not damp velocity if we're trying to maintain equilibrium
            vec2f const absoluteVelocity =
                relativeVelocity * ((npc.Kind != NpcKindType::Human || npc.KindSpecificState.HumanNpcState.EquilibriumTorque == 0.0f) ? (1.0f - gameParameters.GlobalDamping) : 1.0f)
                + meshVelocity;

            particles.SetVelocity(npcParticle.ParticleIndex, absoluteVelocity);
            npcParticleConstrainedState.MeshRelativeVelocity = relativeVelocity;

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
                float const den = npcParticleConstrainedState.CurrentBCoords.BCoords[vi] - segmentTrajectoryEndBarycentricCoords[vi];
                float const t = npcParticleConstrainedState.CurrentBCoords.BCoords[vi] / den;

#ifdef _DEBUG
                diags[vi].emplace(den, t);
                diags[vi]->IntersectionPoint =
                    npcParticleConstrainedState.CurrentBCoords.BCoords
                    + (segmentTrajectoryEndBarycentricCoords - npcParticleConstrainedState.CurrentBCoords.BCoords) * t;
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
            npcParticleConstrainedState.CurrentBCoords.BCoords[(intersectionVertexOrdinal + 1) % 3] * (1.0f - minIntersectionT)
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

        npcParticleConstrainedState.CurrentBCoords.BCoords = intersectionBarycentricCoords;

        //
        // Check if impacted with floor
        //

        vec2f const intersectionAbsolutePosition = shipMesh.GetTriangles().FromBarycentricCoordinates(
            intersectionBarycentricCoords,
            npcParticleConstrainedState.CurrentBCoords.TriangleElementIndex,
            shipMesh.GetPoints());

        if (IsEdgeFloorToParticle(npcParticleConstrainedState.CurrentBCoords.TriangleElementIndex, intersectionEdgeOrdinal, npc, npcParticleOrdinal, mParticles, shipMesh))
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
                shipMesh.GetTriangles().GetSubSpringVector(
                    npcParticleConstrainedState.CurrentBCoords.TriangleElementIndex,
                    intersectionEdgeOrdinal,
                    shipMesh.GetPoints())
                .normalise();
            vec2f const intersectionEdgeNormal = intersectionEdgeDir.to_perpendicular();

            BounceConstrainedNpcParticle(
                npc,
                npcParticleOrdinal,
                trajectory,
                hasMovedInStep || ((intersectionAbsolutePosition - segmentTrajectoryStartAbsolutePosition).length() != 0.0f),
                intersectionAbsolutePosition,
                intersectionEdgeNormal,
                meshVelocity,
                segmentDt,
                particles,
                currentSimulationTime,
                gameParameters);

            // Remember that - at least for this frame - we are non-inertial on this floor
            //
            // We do this to interrupt long streaks off-floor while we bounce multiple times on the floor
            // while walking; without interruptions the multiple bounces, with no intervening non-inertial
            // step, would cause the NPC to stop walking. Afte rall we're not lying as we're really
            // non-inertial on this floor
            npcParticleConstrainedState.CurrentVirtualFloor.emplace(npcParticleConstrainedState.CurrentBCoords.TriangleElementIndex, intersectionEdgeOrdinal);

            // Return (mesh-relative) distance traveled with this move
            return (segmentTrajectoryEndAbsolutePosition - segmentTrajectoryStartAbsolutePosition).length();
        }

        //
        // Not floor, climb over edge
        //

        LogNpcDebug("      Climbing over non-floor edge");

        // Find opposite triangle
        auto const & oppositeTriangleInfo = shipMesh.GetTriangles().GetOppositeTriangle(npcParticleConstrainedState.CurrentBCoords.TriangleElementIndex, intersectionEdgeOrdinal);
        if (oppositeTriangleInfo.TriangleElementIndex == NoneElementIndex || shipMesh.GetTriangles().IsDeleted(oppositeTriangleInfo.TriangleElementIndex))
        {
            //
            // Become free
            //

            LogNpcDebug("      No opposite triangle found, becoming free");

            //
            // Move to endpoint and exit, consuming whole quantum
            //

            TransitionParticleToFreeState(npc, npcParticleOrdinal);

            UpdateNpcParticle_Free(
                npcParticle,
                particleStartAbsolutePosition,
                segmentTrajectoryEndAbsolutePosition,
                particles,
                gameParameters);

            vec2f const totalTraveledVector = intersectionAbsolutePosition - particleStartAbsolutePosition; // We consider constrained portion only
            return totalTraveledVector.length();
        }
        else
        {
            //
            // Move to edge of opposite triangle
            //

            LogNpcDebug("      Moving to edge ", oppositeTriangleInfo.EdgeOrdinal, " of opposite triangle ", oppositeTriangleInfo.TriangleElementIndex);

            // Calculate new current barycentric coords (wrt new triangle)
            bcoords3f newBarycentricCoords; // In new triangle
            newBarycentricCoords[(oppositeTriangleInfo.EdgeOrdinal + 2) % 3] = 0.0f;
            newBarycentricCoords[oppositeTriangleInfo.EdgeOrdinal] = npcParticleConstrainedState.CurrentBCoords.BCoords[(intersectionEdgeOrdinal + 1) % 3];
            newBarycentricCoords[(oppositeTriangleInfo.EdgeOrdinal + 1) % 3] = npcParticleConstrainedState.CurrentBCoords.BCoords[intersectionEdgeOrdinal];

            LogNpcDebug("      B-Coords: ", npcParticleConstrainedState.CurrentBCoords.BCoords, " -> ", newBarycentricCoords);

            assert(newBarycentricCoords.is_on_edge_or_internal());

            // Move to triangle and coords
            npcParticleConstrainedState.CurrentBCoords = AbsoluteTriangleBCoords(oppositeTriangleInfo.TriangleElementIndex, newBarycentricCoords);

            // Translate target coords to this triangle, for next iteration
            auto const oldSegmentTrajectoryEndBarycentricCoords = segmentTrajectoryEndBarycentricCoords; // For logging
            // Note: here we introduce a lot of error - the target bary coords are not anymore
            // guaranteed to lie exactly on the (continuation of the) edge
            segmentTrajectoryEndBarycentricCoords = shipMesh.GetTriangles().ToBarycentricCoordinatesInsideEdge(
                segmentTrajectoryEndAbsolutePosition,
                oppositeTriangleInfo.TriangleElementIndex,
                shipMesh.GetPoints(),
                oppositeTriangleInfo.EdgeOrdinal);

            LogNpcDebug("      TrajEndB-Coords: ", oldSegmentTrajectoryEndBarycentricCoords, " -> ", segmentTrajectoryEndBarycentricCoords);

            // Continue
        }
    }
}

inline Npcs::NavigateVertexOutcome Npcs::NavigateVertex(
    StateType const & npc,
    int npcParticleOrdinal,
    std::optional<TriangleAndEdge> const & walkedEdge,
    int vertexOrdinal, // Mutable
    vec2f const & trajectory,
    vec2f const & trajectoryEndAbsolutePosition,
    bcoords3f trajectoryEndBarycentricCoords, // Mutable
    Ship const & shipMesh,
    NpcParticles const & particles)
{
    // Rules of the code:
    // - We communicate the particle's final state (triangle & bcoords) upon leaving
    //      - And most importantly, we don't move
    // - We take into account *actual* (resultant physical) movement (i.e. trajectory), rather than *intended* (walkdir) movement

    //
    // See whether this particle is the primary (feet) of a walking human
    //

    auto const & npcParticle = npc.ParticleMesh.Particles[npcParticleOrdinal];
    assert(npcParticle.ConstrainedState.has_value());

    if (npc.Kind == NpcKindType::Human
        && npc.KindSpecificState.HumanNpcState.CurrentBehavior == StateType::KindSpecificStateType::HumanNpcStateType::BehaviorType::Constrained_Walking
        && walkedEdge.has_value()
        && npcParticleOrdinal == 0)
    {
        //
        // Navigate around this vertex according to CW/CCW and choose floors to walk on; we return any of these:
        // - We are directed inside a triangle (iff no floors have been found)
        // - We have chosen a floor
        // - We have detected an impact against a floor (iff no floors have been found)
        // - We become free
        //

        //
        // When in walking state and arriving at a vertex at which there are more than two *viable* floors there (incl.incoming, so >= 2 + 1), choose which one to take
        //    - It's like "not seeing" certain floors
        // Note: *viable* == with right slope for walking on it
        //    - We only consider those floors that are in a sector centered around walk(face) dir, up amplitude equal to = / -MaxSlopeForWalking, and down amplitude less than vertical
        //    - This allows us to take a down "stair" that is almost vertical
        //    - We only choose in our direction because the simulation is still very much physical, i.e.informed by trajectory
        //
        // We divide floor candidates in two groups: easy slope and hard slope (i.e. almost falling); we only choose among hard slopes if there
        // are no easy slopes. Rationale: NPC going about S and hitting wall-floor conjunction; we don't want it to take wall going down
        //      |
        //   ---|
        //    */|
        //    / |
        //

        struct AbsoluteTriangleBCoordsAndEdge
        {
            AbsoluteTriangleBCoords TriangleBCoords;
            int EdgeOrdinal;
        };

        std::array<AbsoluteTriangleBCoordsAndEdge, GameParameters::MaxSpringsPerPoint> floorCandidatesEasySlope;
        size_t floorCandidatesEasySlopeCount = 0;
        std::array<AbsoluteTriangleBCoordsAndEdge, GameParameters::MaxSpringsPerPoint> floorCandidatesHardSlope;
        size_t floorCandidatesHardSlopeCount = 0;
        std::optional<AbsoluteTriangleBCoordsAndEdge> firstBounceableFloor;
        std::optional<AbsoluteTriangleBCoords> firstTriangleInterior;

        AbsoluteTriangleBCoords currentAbsoluteBCoords = npcParticle.ConstrainedState->CurrentBCoords;

        vec2f const vertexAbsolutePosition = shipMesh.GetPoints().GetPosition(shipMesh.GetTriangles().GetPointIndices(currentAbsoluteBCoords.TriangleElementIndex)[vertexOrdinal]);

        // The two vertices around the vertex we are on - seen in clockwise order
        int nextVertexOrdinal = (vertexOrdinal + 1) % 3;
        int prevVertexOrdinal = (vertexOrdinal + 2) % 3;

        // Determine orientation of our visit (CW vs CCW)
        //
        // Assumption is that trajectoryEndBarycentricCoords is *outside* triangle,
        // and we are at the last vertex possible in this triangle before leaving it
        //
        //   0    * b[0] > b[2]
        //   |\  /
        //   | \
        //   | *\
        //  2----1
        //      /
        //     *  b[0] < b[2]

        RotationDirectionType const orientation = (trajectoryEndBarycentricCoords[prevVertexOrdinal] <= trajectoryEndBarycentricCoords[nextVertexOrdinal])
            ? RotationDirectionType::CounterClockwise
            : RotationDirectionType::Clockwise;

        // Remember floor geometry of starting edge
        NpcFloorGeometryType const initialFloorGeometry = shipMesh.GetTriangles().GetSubSpringNpcFloorGeometry(walkedEdge->TriangleElementIndex, walkedEdge->EdgeOrdinal);

        for (int iIter = 0; ; ++iIter)
        {
            LogNpcDebug("    NavigateVertex_Walking: iter=", iIter, " orientation=", orientation == RotationDirectionType::Clockwise ? "CW" : "CCW",
                " nCandidatesEasy=", floorCandidatesEasySlopeCount, " nCandidatesHard=", floorCandidatesHardSlopeCount, " hasBounceableFloor=", firstBounceableFloor.has_value() ? "T" : "F",
                " hasFirstTriangleInterior=", firstTriangleInterior.has_value() ? "T" : "F");

            assert(iIter < GameParameters::MaxSpringsPerPoint); // Detect and debug-break on infinite loops

            // Pre-conditions: we are at this vertex
            assert(currentAbsoluteBCoords.BCoords[vertexOrdinal] == 1.0f);
            assert(currentAbsoluteBCoords.BCoords[nextVertexOrdinal] == 0.0f);
            assert(currentAbsoluteBCoords.BCoords[prevVertexOrdinal] == 0.0f);

            LogNpcDebug("      Triangle=", currentAbsoluteBCoords.TriangleElementIndex, " Vertex=", vertexOrdinal, " BCoords=", currentAbsoluteBCoords.BCoords,
                " TrajectoryEndBarycentricCoords=", trajectoryEndBarycentricCoords);

            //
            // Check whether we are directed towards the *interior* of this triangle, if we don't
            // know yet which triangle we'd be going inside
            //

            if (!firstTriangleInterior.has_value()
                && trajectoryEndBarycentricCoords[prevVertexOrdinal] >= 0.0f
                && trajectoryEndBarycentricCoords[nextVertexOrdinal] >= 0.0f)
            {
                LogNpcDebug("      Trajectory extends inside triangle, remembering it");

                // Remember these absolute BCoords as we'll go there if we don't have any candidates nor we bounce on a floor
                firstTriangleInterior = currentAbsoluteBCoords;

                // Continue, so that we may find candidates at a lower slope
            }

            //
            // Find next edge that we cross
            //

            int const crossedEdgeOrdinal = (orientation == RotationDirectionType::CounterClockwise)
                ? vertexOrdinal // Next edge
                : (vertexOrdinal + 2) % 3; // Previous edge

            LogNpcDebug("      Next crossed edge: ", crossedEdgeOrdinal);

            //
            // Check whether we've gone too far around
            //

            auto const & nextVertexPos = shipMesh.GetPoints().GetPosition(shipMesh.GetTriangles().GetPointIndices(currentAbsoluteBCoords.TriangleElementIndex)[nextVertexOrdinal]);
            auto const & prevVertexPos = shipMesh.GetPoints().GetPosition(shipMesh.GetTriangles().GetPointIndices(currentAbsoluteBCoords.TriangleElementIndex)[prevVertexOrdinal]);
            if (prevVertexPos.x <= vertexAbsolutePosition.x && nextVertexPos.x > vertexAbsolutePosition.x)
            {
                LogNpcDebug("      Gone too far - may stop search here");
                break;
            }

            //
            // Check whether this new edge is floor
            //

            if (IsEdgeFloorToParticle(currentAbsoluteBCoords.TriangleElementIndex, crossedEdgeOrdinal, npc, npcParticleOrdinal, particles, shipMesh))
            {
                //
                // Encountered floor
                //

                auto const crossedEdgeFloorGeometry = shipMesh.GetTriangles().GetSubSpringNpcFloorGeometry(currentAbsoluteBCoords.TriangleElementIndex, crossedEdgeOrdinal);

                LogNpcDebug("        Crossed edge is floor (geometry:", int(crossedEdgeFloorGeometry), ")");

                //
                // Check whether it's a viable floor, i.e. whether its direction is:
                //    - x: in direction of movement (including close to 0, i.e. almost vertical)
                //    - y: lower than MaxHumanNpcWalkSinSlope
                //
                // Notes:
                //  - Here we check viability wrt *actual* (resultant physical) movement, rather than *intended* (walkdir) movement
                //  - Viability is based on actual gravity direction - not on *apparent* gravity
                //

                vec2f const crossedEdgeDir = shipMesh.GetTriangles().GetSubSpringVector(
                    currentAbsoluteBCoords.TriangleElementIndex,
                    crossedEdgeOrdinal,
                    shipMesh.GetPoints()).normalise();

                bool const isViable =
                    crossedEdgeDir.x < 0.0f
                    && (
                        (orientation == RotationDirectionType::CounterClockwise)
                        ? crossedEdgeDir.y <= GameParameters::MaxHumanNpcWalkSinSlope
                        : -crossedEdgeDir.y <= GameParameters::MaxHumanNpcWalkSinSlope
                        );

                LogNpcDebug("          Edge is ", isViable ? "viable" : "non-viable", " (edgeDir=", crossedEdgeDir, ")");

                if (isViable)
                {
                    //
                    // Viable floor - add to candidates
                    //

                    if (std::abs(crossedEdgeDir.y) < 0.98f) // Magic slope
                    {
                        floorCandidatesEasySlope[floorCandidatesEasySlopeCount++] = { currentAbsoluteBCoords, crossedEdgeOrdinal };

                        LogNpcDebug("          Added to easy candidates: new count=", floorCandidatesEasySlopeCount);
                    }
                    else
                    {
                        floorCandidatesHardSlope[floorCandidatesHardSlopeCount++] = { currentAbsoluteBCoords, crossedEdgeOrdinal };

                        LogNpcDebug("          Added to hard candidates: new count=", floorCandidatesHardSlopeCount);
                    }
                }
                else
                {
                    //
                    // Not viable - see if it's bounceable
                    //
                    // Bounceable: iff we hit it according to our current direction
                    //   - Note: all edges are in our direction until we've found a triangle interior
                    //

                    if (!firstTriangleInterior.has_value())
                    {
                        LogNpcDebug("          Edge is bounceable");

                        //
                        // Bounceable: store it so we may bounce on it in case there are no candidates
                        //
                        // Note: we give same-depth priority over other-depth, because:
                        //  - If no other walls - nor candidates - exist, we're ok with bouncing on S at --> _\
                        //  - But when given a choice with vertical - e.g. --> _\| - we prefer bouncing on vertical, for better physics
                        //      - Thus honoring semi-invisible nature of S
                        //
                        // Also: we want to enforce a concept of entering depth 1 areas from depth 2 areas, hence,
                        // if we're on S, we want to remember the _last_ of H/V we encounter (which would be the
                        // second at most, as there cannot be any third bounceable H/V if we're on S), as that is
                        // consistent with being "inside a depth 1 area"
                        //

                        bool doRememberBounceableFloor = false;

                        if (!firstBounceableFloor.has_value())
                        {
                            doRememberBounceableFloor = true;
                        }
                        else
                        {
                            auto const currentBounceableFloorGeometry = shipMesh.GetTriangles().GetSubSpringNpcFloorGeometry(firstBounceableFloor->TriangleBCoords.TriangleElementIndex, firstBounceableFloor->EdgeOrdinal);
                            if ((NpcFloorGeometryDepth(currentBounceableFloorGeometry) != NpcFloorGeometryDepth(initialFloorGeometry) && NpcFloorGeometryDepth(crossedEdgeFloorGeometry) == NpcFloorGeometryDepth(initialFloorGeometry))
                                ||
                                (
                                    NpcFloorGeometryDepth(initialFloorGeometry) == NpcFloorGeometryDepthType::Depth2
                                    && NpcFloorGeometryDepth(currentBounceableFloorGeometry) == NpcFloorGeometryDepthType::Depth1
                                    && NpcFloorGeometryDepth(crossedEdgeFloorGeometry) == NpcFloorGeometryDepthType::Depth1
                                    ))
                            {
                                doRememberBounceableFloor = true;
                            }
                        }

                        if (doRememberBounceableFloor)
                        {
                            LogNpcDebug("            Remembering bounceable floor");
                            firstBounceableFloor = { currentAbsoluteBCoords, crossedEdgeOrdinal };
                        }
                    }
                }

                //
                // Check now if it's an impenetrable wall
                //
                // Impenetrable: iff HonV or VonH
                //  - We allow SonS to be penetrable, so that while we go down ladder and encounter ladder going up, we may go beyond it if there's a floor behind it
                //

                if ((crossedEdgeFloorGeometry == NpcFloorGeometryType::Depth1H && initialFloorGeometry == NpcFloorGeometryType::Depth1V)
                    || (crossedEdgeFloorGeometry == NpcFloorGeometryType::Depth1V && initialFloorGeometry == NpcFloorGeometryType::Depth1H))
                {
                    LogNpcDebug("          Impenetrable, stopping here");

                    // If this was viable, we'll choose it; if it was not viable, we'll bounce on it if it was in direction
                    //  - Note: if it's not in direction (i.e. if we have an interior), and if we have no candidates and no bounceables - then we will take the interior

                    break;
                }
            }

            //
            // Climb over edge
            //

            // Find opposite triangle
            auto const & oppositeTriangleInfo = shipMesh.GetTriangles().GetOppositeTriangle(currentAbsoluteBCoords.TriangleElementIndex, crossedEdgeOrdinal);
            if (oppositeTriangleInfo.TriangleElementIndex == NoneElementIndex || shipMesh.GetTriangles().IsDeleted(oppositeTriangleInfo.TriangleElementIndex))
            {
                //
                // Found a free region
                //

                LogNpcDebug("      No opposite triangle found, found free region");

                // If we've found this free region before anything else, become free; otherwise,
                // stop search now and proceed with what we've found
                if (floorCandidatesEasySlopeCount + floorCandidatesHardSlopeCount != 0
                    || firstBounceableFloor.has_value()
                    || firstTriangleInterior.has_value())
                {
                    LogNpcDebug("        Free region after other options, stopping here");

                    break;
                }
                else
                {
                    // Detected free

                    LogNpcDebug("        Free region before anything else, detected free");

                    return NavigateVertexOutcome::MakeBecomeFreeOutcome();
                }
            }

            LogNpcDebug("      Opposite triangle found: ", oppositeTriangleInfo.TriangleElementIndex);

            // See whether we've gone around
            if (oppositeTriangleInfo.TriangleElementIndex == npcParticle.ConstrainedState->CurrentBCoords.TriangleElementIndex)
            {
                // Time to stop

                LogNpcDebug("        Opposite triangle is self (done full round), stopping search");

                // We must have found a triangle into which we go
                assert(firstTriangleInterior.has_value());

                break;
            }

            //
            // Move to triangle
            //

            LogNpcDebug("      Continuing search from edge ", oppositeTriangleInfo.EdgeOrdinal, " of opposite triangle ", oppositeTriangleInfo.TriangleElementIndex);

            // Calculate new current barycentric coords (wrt opposite triangle - note that we haven't moved)
            bcoords3f newBarycentricCoords; // In new triangle
            newBarycentricCoords[(oppositeTriangleInfo.EdgeOrdinal + 2) % 3] = 0.0f;
            newBarycentricCoords[oppositeTriangleInfo.EdgeOrdinal] = currentAbsoluteBCoords.BCoords[(crossedEdgeOrdinal + 1) % 3];
            newBarycentricCoords[(oppositeTriangleInfo.EdgeOrdinal + 1) % 3] = currentAbsoluteBCoords.BCoords[crossedEdgeOrdinal];

            LogNpcDebug("        B-Coords: ", currentAbsoluteBCoords.BCoords, " -> ", newBarycentricCoords);

            assert(newBarycentricCoords.is_on_edge_or_internal());

            // Move to triangle and b-coords
            currentAbsoluteBCoords = AbsoluteTriangleBCoords(oppositeTriangleInfo.TriangleElementIndex, newBarycentricCoords);

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

            // The two vertices around the vertex we are on - seen in clockwise order
            nextVertexOrdinal = (vertexOrdinal + 1) % 3;
            prevVertexOrdinal = (vertexOrdinal + 2) % 3;

            //
            // Translate target bary coords - if we still need them
            //

            if (!firstTriangleInterior.has_value())
            {
                trajectoryEndBarycentricCoords = shipMesh.GetTriangles().ToBarycentricCoordinatesInsideEdge(
                    trajectoryEndAbsolutePosition,
                    oppositeTriangleInfo.TriangleElementIndex,
                    shipMesh.GetPoints(),
                    oppositeTriangleInfo.EdgeOrdinal);

                LogNpcDebug("        New TrajEndB-Coords: ", trajectoryEndBarycentricCoords);
            }
        }

        //
        // Process results
        //

        if (floorCandidatesEasySlopeCount > 0 || floorCandidatesHardSlopeCount > 0)
        {
            // Only choose among hard ones if there are no easy ones

            AbsoluteTriangleBCoordsAndEdge chosenFloor;
            if (floorCandidatesEasySlopeCount > 0)
            {
                size_t chosenCandidateIndex = (floorCandidatesEasySlopeCount == 1)
                    ? 0
                    : GameRandomEngine::GetInstance().Choose(floorCandidatesEasySlopeCount);
                chosenFloor = floorCandidatesEasySlope[chosenCandidateIndex];

                LogNpcDebug("    Chosen easy candidate ", chosenCandidateIndex, " (", chosenFloor.TriangleBCoords.TriangleElementIndex, ":", chosenFloor.TriangleBCoords.BCoords,
                    ", ", chosenFloor.EdgeOrdinal, ") out of ", floorCandidatesEasySlopeCount);
            }
            else
            {
                size_t chosenCandidateIndex = (floorCandidatesHardSlopeCount == 1)
                    ? 0
                    : GameRandomEngine::GetInstance().Choose(floorCandidatesHardSlopeCount);
                chosenFloor = floorCandidatesHardSlope[chosenCandidateIndex];

                LogNpcDebug("    Chosen hard candidate ", chosenCandidateIndex, " (", chosenFloor.TriangleBCoords.TriangleElementIndex, ":", chosenFloor.TriangleBCoords.BCoords,
                    ", ", chosenFloor.EdgeOrdinal, ") out of ", floorCandidatesEasySlopeCount);
            }

            return NavigateVertexOutcome::MakeContinueAlongFloorOutcome(
                chosenFloor.TriangleBCoords,
                chosenFloor.EdgeOrdinal);
        }
        else if (firstBounceableFloor.has_value())
        {
            // Impact on this floor

            LogNpcDebug("    Impact on floor ", firstBounceableFloor->TriangleBCoords.TriangleElementIndex, ":", firstBounceableFloor->EdgeOrdinal);

            return NavigateVertexOutcome::MakeImpactOnFloorOutcome(
                firstBounceableFloor->TriangleBCoords,
                firstBounceableFloor->EdgeOrdinal);
        }
        else
        {
            // Inside triangle

            assert(firstTriangleInterior.has_value());

            LogNpcDebug("    Going inside triangle ", *firstTriangleInterior);

            return NavigateVertexOutcome::MakeContinueToInteriorOutcome(*firstTriangleInterior);
        }
    }
    else
    {
        //
        // Navigate this vertex along the trajectory direction, until any of these:
        // - We are directed inside a triangle
        // - We have encountered a floor
        // - We have detected an impact against a floor
        // - We become free
        //

        AbsoluteTriangleBCoords currentAbsoluteBCoords = npcParticle.ConstrainedState->CurrentBCoords;

        for (int iIter = 0; ; ++iIter)
        {
            LogNpcDebug("    NavigateVertex_NonWalking: iter=", iIter);

            assert(iIter < 8); // Detect and debug-break on infinite loops

            // The two vertices around the vertex we are on - seen in clockwise order
            int const nextVertexOrdinal = (vertexOrdinal + 1) % 3;
            int const prevVertexOrdinal = (vertexOrdinal + 2) % 3;

            // Pre-conditions: we are at this vertex
            assert(currentAbsoluteBCoords.BCoords[vertexOrdinal] == 1.0f);
            assert(currentAbsoluteBCoords.BCoords[nextVertexOrdinal] == 0.0f);
            assert(currentAbsoluteBCoords.BCoords[prevVertexOrdinal] == 0.0f);

            LogNpcDebug("      Triangle=", currentAbsoluteBCoords.TriangleElementIndex, " Vertex=", vertexOrdinal, " TrajectoryEndBarycentricCoords=", trajectoryEndBarycentricCoords);

            //
            // Check whether we are directed towards the *interior* of this triangle - including its edges
            //

            if (trajectoryEndBarycentricCoords[prevVertexOrdinal] >= 0.0f
                && trajectoryEndBarycentricCoords[nextVertexOrdinal] >= 0.0f)
            {
                //
                // We go inside (or on) this triangle - stop where we are, we'll then check trajectory in new situation
                //

                LogNpcDebug("      Going inside triangle ", currentAbsoluteBCoords);

                return NavigateVertexOutcome::MakeContinueToInteriorOutcome(currentAbsoluteBCoords);
            }

            //
            // Find next edge that we intersect at this vertex
            //

            RotationDirectionType const orientation = (trajectoryEndBarycentricCoords[prevVertexOrdinal] <= trajectoryEndBarycentricCoords[nextVertexOrdinal])
                ? RotationDirectionType::CounterClockwise
                : RotationDirectionType::Clockwise;

            int const crossedEdgeOrdinal = (orientation == RotationDirectionType::CounterClockwise)
                ? vertexOrdinal // Next edge
                : (vertexOrdinal + 2) % 3; // Prev edge

            LogNpcDebug("      Trajectory crosses triangle: crossedEdgeOrdinal=", crossedEdgeOrdinal);

            //
            // Check whether this new edge is floor
            //

            if (IsEdgeFloorToParticle(currentAbsoluteBCoords.TriangleElementIndex, crossedEdgeOrdinal, npc, npcParticleOrdinal, particles, shipMesh))
            {
                //
                // Encountered floor
                //

                LogNpcDebug("      Crossed edge is floor");

                // Determine viability of floor: check angle between desired (original) trajectory and edge

                vec2f const floorEdgeDir = shipMesh.GetTriangles().GetSubSpringVector(
                    currentAbsoluteBCoords.TriangleElementIndex,
                    crossedEdgeOrdinal,
                    shipMesh.GetPoints()).normalise();
                vec2f const floorEdgeNormal = floorEdgeDir.to_perpendicular();

                float const trajProjOntoEdgeNormal = trajectory.normalise().dot(floorEdgeNormal);
                bool const isViable = (trajProjOntoEdgeNormal <= 0.71f); // PI/4+

                LogNpcDebug("        Edge is ", isViable ? "viable" : "non-viable", " (floorEdgeNormal=", floorEdgeNormal, ")");

                if (isViable)
                {
                    return NavigateVertexOutcome::MakeContinueAlongFloorOutcome(currentAbsoluteBCoords, crossedEdgeOrdinal);
                }
                else
                {
                    return NavigateVertexOutcome::MakeImpactOnFloorOutcome(currentAbsoluteBCoords, crossedEdgeOrdinal);
                }
            }

            //
            // Not floor, climb over edge
            //

            // Find opposite triangle
            auto const & oppositeTriangleInfo = shipMesh.GetTriangles().GetOppositeTriangle(currentAbsoluteBCoords.TriangleElementIndex, crossedEdgeOrdinal);
            if (oppositeTriangleInfo.TriangleElementIndex == NoneElementIndex || shipMesh.GetTriangles().IsDeleted(oppositeTriangleInfo.TriangleElementIndex))
            {
                //
                // Become free
                //

                LogNpcDebug("      No opposite triangle found, becoming free");

                return NavigateVertexOutcome::MakeBecomeFreeOutcome();
            }

            //
            // Move to triangle
            //

            LogNpcDebug("      Moving to edge ", oppositeTriangleInfo.EdgeOrdinal, " of opposite triangle ", oppositeTriangleInfo.TriangleElementIndex);

            // Calculate new current barycentric coords (wrt opposite triangle - note that we haven't moved)
            bcoords3f newBarycentricCoords; // In new triangle
            newBarycentricCoords[(oppositeTriangleInfo.EdgeOrdinal + 2) % 3] = 0.0f;
            newBarycentricCoords[oppositeTriangleInfo.EdgeOrdinal] = currentAbsoluteBCoords.BCoords[(crossedEdgeOrdinal + 1) % 3];
            newBarycentricCoords[(oppositeTriangleInfo.EdgeOrdinal + 1) % 3] = currentAbsoluteBCoords.BCoords[crossedEdgeOrdinal];

            LogNpcDebug("      B-Coords: ", currentAbsoluteBCoords.BCoords, " -> ", newBarycentricCoords);

            assert(newBarycentricCoords.is_on_edge_or_internal());

            // Move to triangle and b-coords
            currentAbsoluteBCoords = AbsoluteTriangleBCoords(oppositeTriangleInfo.TriangleElementIndex, newBarycentricCoords);

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

            trajectoryEndBarycentricCoords = shipMesh.GetTriangles().ToBarycentricCoordinatesInsideEdge(
                trajectoryEndAbsolutePosition,
                oppositeTriangleInfo.TriangleElementIndex,
                shipMesh.GetPoints(),
                oppositeTriangleInfo.EdgeOrdinal);

            LogNpcDebug("      TrajEndB-Coords: ", trajectoryEndBarycentricCoords);
        }
    }
}

void Npcs::BounceConstrainedNpcParticle(
    StateType & npc,
    int npcParticleOrdinal,
    vec2f const & trajectory,
    bool hasMovedInStep,
    vec2f const & bouncePosition,
    vec2f const & bounceEdgeNormal,
    vec2f const meshVelocity,
    float dt,
    NpcParticles & particles,
    float currentSimulationTime,
    GameParameters const & gameParameters) const
{
    auto & npcParticle = npc.ParticleMesh.Particles[npcParticleOrdinal];
    assert(npcParticle.ConstrainedState.has_value());

    // Calculate apparent velocity:
    //  - If we have made any movement in this step, then this is an actual collision
    //    after having traveled part of the (segment) trajectory; we calculate then the
    //    collision velocity by means of the (segment) trajectory we wanted to travel during
    //    the (segment) dt
    //      - Note: here we take into account new forces that sprung up just before this iteration
    //  - Otherwise, the previous iteration stopped just short of the bounce, either heads-on
    //    (in which case we have MRVelocity) or because we are resting here (in which case we
    //    have no MRVelocity); we thus infer whether this is a real collision or not by means
    //    of the particle's mesh-relative velocity
    //      - Note: at this moment the particle's mesh-relative velocity does _not_ include
    //        yet an eventual mesh displacement that took place between the previous step and this
    //        step, nor any new forces that sprung up just before this iteration
    //          - We assume new forces that sprung up just before this iteration have no effect on
    //            bounce, because particle is already on floor
    //
    // If we end up with zero apparent velocity, it means we're just resting on this edge,
    // at this position

    vec2f const apparentParticleVelocity = hasMovedInStep
        ? trajectory / dt
        : npcParticle.ConstrainedState->MeshRelativeVelocity;
    float const apparentParticleVelocityAlongNormal = apparentParticleVelocity.dot(bounceEdgeNormal);

    LogNpcDebug("      BounceConstrainedNpcParticle: apparentParticleVelocity=", apparentParticleVelocity, " (hasMovedInStep=", hasMovedInStep,
        " meshRelativeVelocity=", npcParticle.ConstrainedState->MeshRelativeVelocity, ")");

    if (apparentParticleVelocityAlongNormal != 0.0f)
    {
        // Decompose apparent particle velocity into normal and tangential
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

        // Calculate whole response (which, given that we've been working in *apparent* space (we've calc'd the collision response to *trajectory* which is apparent displacement)),
        // is a relative velocity (relative to mesh)
        vec2f const resultantRelativeVelocity = (normalResponse + tangentialResponse);

        // Do not damp velocity if we're trying to maintain equilibrium
        vec2f const resultantAbsoluteVelocity =
            resultantRelativeVelocity * ((npc.Kind != NpcKindType::Human || npc.KindSpecificState.HumanNpcState.EquilibriumTorque == 0.0f) ? (1.0f - gameParameters.GlobalDamping) : 1.0f)
            + meshVelocity;

        LogNpcDebug("        Impact: trajectory=", trajectory, " apparentParticleVelocity=", apparentParticleVelocity, " nr=", normalResponse, " tr=", tangentialResponse, " rr=", resultantRelativeVelocity);

        //
        // Set position and velocity
        //

        particles.SetPosition(npcParticle.ParticleIndex, bouncePosition);

        particles.SetVelocity(npcParticle.ParticleIndex, resultantAbsoluteVelocity);
        npcParticle.ConstrainedState->MeshRelativeVelocity = resultantRelativeVelocity;

        //
        // Publish impact
        //

        OnImpact(
            npc,
            npcParticleOrdinal,
            normalVelocity,
            bounceEdgeNormal,
            currentSimulationTime);
    }
    else
    {
        LogNpcDebug("        Not an impact");

        //
        // Set position and velocity
        //

        particles.SetPosition(npcParticle.ParticleIndex, bouncePosition);

        particles.SetVelocity(npcParticle.ParticleIndex, meshVelocity);
        npcParticle.ConstrainedState->MeshRelativeVelocity = vec2f::zero();
    }
}

void Npcs::OnImpact(
    StateType & npc,
    int npcParticleOrdinal,
    vec2f const & normalResponse,
    vec2f const & bounceEdgeNormal, // Pointing outside of triangle
    float currentSimulationTime) const
{
    LogNpcDebug("    OnImpact(mag=", normalResponse.length(), ", bounceEdgeNormal=", bounceEdgeNormal, ")");

    // Human state machine
    if (npc.Kind == NpcKindType::Human)
    {
        OnHumanImpact(
            npc,
            npcParticleOrdinal,
            normalResponse,
            bounceEdgeNormal,
            currentSimulationTime);
    }
}

void Npcs::UpdateNpcAnimation(
    StateType & npc,
    int npcParticleOrdinal,
    float currentSimulationTime,
    Ship const & shipMesh)
{
    if (npc.Kind == NpcKindType::Human && npcParticleOrdinal == 0) // Take the primary as the only representative of a human
    {
        auto & humanNpcState = npc.KindSpecificState.HumanNpcState;
        auto & animationState = humanNpcState.AnimationState;
        using HumanNpcStateType = StateType::KindSpecificStateType::HumanNpcStateType;

        assert(npc.ParticleMesh.Particles.size() == 2);
        ElementIndex const primaryParticleIndex = npc.ParticleMesh.Particles[0].ParticleIndex;
        auto const & primaryContrainedState = npc.ParticleMesh.Particles[0].ConstrainedState;
        ElementIndex const secondaryParticleIndex = npc.ParticleMesh.Particles[1].ParticleIndex;
        auto const & secondaryConstrainedState = npc.ParticleMesh.Particles[1].ConstrainedState;

        //
        // Angles and thigh
        //

        // Target: begin with current
        FS_ALIGN16_BEG LimbVector targetAngles(animationState.LimbAngles) FS_ALIGN16_END;

        float convergenceRate = 0.0f;

        // Stuff we calc in some cases and which we need again later for lengths
        float humanEdgeAngle = 0.0f;
        float adjustedStandardHumanHeight = 0.0f;
        vec2f edg1, edg2, edgVector, edgDir;
        vec2f feetPosition, actualBodyVector, actualBodyDir;
        float periodicValue = 0.0f;

        float targetUpperLegLengthFraction = 1.0f;

        // Angle of human wrt edge until which arm is angled to the max
        // (extent of early stage during rising)
        float constexpr MaxHumanEdgeAngleForArms = 0.40489178628508342331207292900944f;
        //static_assert(MaxHumanEdgeAngleForArms == std::atan(GameParameters::HumanNpcGeometry::ArmLengthFraction / (1.0f - GameParameters::HumanNpcGeometry::HeadLengthFraction)));

        switch (humanNpcState.CurrentBehavior)
        {
            case HumanNpcStateType::BehaviorType::BeingPlaced:
            {
                float const arg =
                    (
                        (currentSimulationTime - humanNpcState.CurrentStateTransitionSimulationTimestamp) * 1.0f
                        + humanNpcState.TotalDistanceTraveledOffEdgeSinceStateTransition * 0.2f
                    ) * (1.0f + humanNpcState.ResultantPanicLevel * 0.2f)
                    * (Pi<float> * 2.0f + npc.RandomNormalizedUniformSeed * 4.0f);

                float const yArms = std::sin(arg);
                targetAngles.RightArm = Pi<float> / 2.0f + Pi<float> / 2.0f * 0.7f * yArms;
                targetAngles.LeftArm = -targetAngles.RightArm;

                float const yLegs = std::sin(arg + npc.RandomNormalizedUniformSeed * Pi<float> * 2.0f);
                targetAngles.RightLeg = (1.0f + yLegs) / 2.0f * Pi<float> / 2.0f * 0.3f;
                targetAngles.LeftLeg = -targetAngles.RightLeg;

                convergenceRate = 0.3f;

                break;
            }

            case HumanNpcStateType::BehaviorType::Constrained_PreRising:
            {
                // Move arms against floor (PI/2 wrt body)

                if (primaryContrainedState.has_value() && primaryContrainedState->CurrentVirtualFloor.has_value())
                {
                    vec2f const edgeVector = shipMesh.GetTriangles().GetSubSpringVector(
                        primaryContrainedState->CurrentVirtualFloor->TriangleElementIndex,
                        primaryContrainedState->CurrentVirtualFloor->EdgeOrdinal,
                        shipMesh.GetPoints());
                    vec2f const head = mParticles.GetPosition(secondaryParticleIndex);
                    vec2f const feet = mParticles.GetPosition(primaryParticleIndex);

                    float const humanFloorAlignment = (head - feet).dot(edgeVector);

                    float constexpr MaxArmAngle = Pi<float> / 2.0f;
                    float constexpr OtherArmDeltaAngle = 0.3f;

                    if (humanFloorAlignment >= 0.0f)
                    {
                        targetAngles.LeftArm = -MaxArmAngle;
                        targetAngles.RightArm = -MaxArmAngle + OtherArmDeltaAngle;
                    }
                    else
                    {
                        targetAngles.RightArm = MaxArmAngle;
                        targetAngles.LeftArm = MaxArmAngle - OtherArmDeltaAngle;
                    }
                }

                // Legs stay

                convergenceRate = 0.09f;

                break;
            }

            case HumanNpcStateType::BehaviorType::Constrained_Rising:
            {
                //
                // Leg and arm that are against floor "help"
                //

                if (primaryContrainedState.has_value() && primaryContrainedState->CurrentVirtualFloor.has_value())
                {
                    // Remember the virtual edge that we're rising against, so we can survive
                    // small bursts of being off the edge
                    humanNpcState.CurrentBehaviorState.Constrained_Rising.VirtualEdgeRisingAgainst = *primaryContrainedState->CurrentVirtualFloor;
                }

                if (humanNpcState.CurrentBehaviorState.Constrained_Rising.VirtualEdgeRisingAgainst.TriangleElementIndex != NoneElementIndex)
                {

                    vec2f const edgeVector = shipMesh.GetTriangles().GetSubSpringVector(
                        humanNpcState.CurrentBehaviorState.Constrained_Rising.VirtualEdgeRisingAgainst.TriangleElementIndex,
                        humanNpcState.CurrentBehaviorState.Constrained_Rising.VirtualEdgeRisingAgainst.EdgeOrdinal,
                        shipMesh.GetPoints());
                    vec2f const head = mParticles.GetPosition(secondaryParticleIndex);
                    vec2f const feet = mParticles.GetPosition(primaryParticleIndex);

                    // First off, we calculate the max possible human-edge vector, considering that
                    // human converges towards full vertical (gravity-only :-( )
                    float maxHumanEdgeAngle = edgeVector.angleCw(vec2f(0.0f, 1.0f)); // Also angle between edge and vertical

                    // Calculate angle between human and edge (angle that we need to rotate human CW to get onto edge)
                    humanEdgeAngle = edgeVector.angleCw(head - feet); // [0.0 ... PI]
                    if (humanEdgeAngle < 0.0f)
                    {
                        // Two possible inaccuracies here:
                        // o -8.11901e-06: this is basically 0.0
                        // o -3.14159: this is basically +PI

                        if (humanEdgeAngle >= -Pi<float> / 2.0f) // Just sentinel for side of inaccuracy
                        {
                            humanEdgeAngle = 0.0f;
                        }
                        else
                        {
                            humanEdgeAngle = Pi<float>;
                        }
                    }

                    bool isOnLeftSide; // Of screen - i.e. head to the left side of the edge (exploiting CWness of edge)
                    if (humanEdgeAngle <= maxHumanEdgeAngle)
                    {
                        isOnLeftSide = true;
                    }
                    else
                    {
                        isOnLeftSide = false;

                        // Normalize to simplify math below
                        humanEdgeAngle = Pi<float> - humanEdgeAngle;
                        maxHumanEdgeAngle = Pi<float> - maxHumanEdgeAngle;
                    }

                    // Max angle of arm wrt body - kept until MaxAngle
                    float constexpr MaxArmAngle = Pi<float> / 2.0f;

                    // Rest angle of arm wrt body - reached when fully erect
                    float constexpr RestArmAngle = HumanNpcStateType::AnimationStateType::InitialArmAngle * 0.3f;

                    // DeltaAngle of other arm
                    float constexpr OtherArmDeltaAngle = 0.3f;

                    // AngleMultiplier of other leg when closing knees
                    float constexpr OtherLegAlphaAngle = 0.7f;

                    // Shortening of angle path for legs becoming straight when closing knees
                    float const AnglePathShorteningForLegsInLateStage = 0.9f;

                    //

                    //  *  0 --> maxHumanEdgeAngle (which is PI/2 when edge is flat)
                    //   \
                    //   |\
                    // -----

                    // Arm: at MaxArmAngle until MaxHumanEdgeAngleForArms, then goes down to rest

                    float targetArm;
                    float targetLeg = 0.0f; // Start with legs closed - we'll change if we're in the early stage of rising and we're L/R

                    if (humanEdgeAngle <= MaxHumanEdgeAngleForArms)
                    {
                        // Early stage

                        // Arms: leave them where they are (MaxArmAngle)
                        targetArm = MaxArmAngle;

                        // Legs: we want a knee (iff we're facing L/R)

                        if (humanNpcState.CurrentFaceOrientation == 0.0f)
                        {
                            targetLeg = MaxArmAngle;
                            targetUpperLegLengthFraction = humanEdgeAngle / MaxHumanEdgeAngleForArms * 0.5f; // 0.0 @ 0.0 -> 0.5 @ MaxHumanEdgeAngleForArms
                        }
                    }
                    else
                    {
                        // Late stage: -> towards maxHumanEdgeAngle

                        // Arms: MaxArmAngle -> RestArmAngle

                        targetArm = MaxArmAngle + (MaxHumanEdgeAngleForArms - humanEdgeAngle) / (MaxHumanEdgeAngleForArms - maxHumanEdgeAngle) * (RestArmAngle - MaxArmAngle); // MaxArmAngle @ MaxHumanEdgeAngleForArms -> RestArmAngle @ maxHumanEdgeAngle

                        // Legs: towards zero

                        if (humanNpcState.CurrentFaceOrientation == 0.0f)
                        {
                            targetLeg = std::max(
                                MaxArmAngle - (MaxHumanEdgeAngleForArms - humanEdgeAngle) / (MaxHumanEdgeAngleForArms - maxHumanEdgeAngle * AnglePathShorteningForLegsInLateStage) * MaxArmAngle, // MaxArmAngle @ MaxHumanEdgeAngleForArms -> 0 @ maxHumanEdgeAngle-e
                                0.0f);
                            targetUpperLegLengthFraction = 0.5f;
                        }
                    }

                    // Knees cannot bend backwards!
                    if ((humanNpcState.CurrentFaceDirectionX > 0.0f && isOnLeftSide)
                        || (humanNpcState.CurrentFaceDirectionX < 0.0f && !isOnLeftSide))
                    {
                        targetLeg *= -1.0f;
                    }

                    if (isOnLeftSide)
                    {
                        targetAngles.LeftArm = -targetArm;
                        targetAngles.RightArm = targetAngles.LeftArm + OtherArmDeltaAngle;

                        targetAngles.LeftLeg = -targetLeg;
                        targetAngles.RightLeg = targetAngles.LeftLeg * OtherLegAlphaAngle;
                    }
                    else
                    {
                        targetAngles.RightArm = targetArm;
                        targetAngles.LeftArm = targetAngles.RightArm - OtherArmDeltaAngle;

                        targetAngles.RightLeg = targetLeg;
                        targetAngles.LeftLeg = targetAngles.RightLeg * OtherLegAlphaAngle;
                    }
                }
                else
                {
                    // Let's survive small bursts and keep current angles; after all we'll lose
                    // this state very quickly if the burst is too long

                    targetUpperLegLengthFraction = animationState.UpperLegLengthFraction;
                }

                convergenceRate = 0.45f;

                break;
            }

            case HumanNpcStateType::BehaviorType::Constrained_Equilibrium:
            {
                // Just small arms angle

                float constexpr ArmsAngle = HumanNpcStateType::AnimationStateType::InitialArmAngle;

                targetAngles.RightArm = ArmsAngle;
                targetAngles.LeftArm = -ArmsAngle;

                targetAngles.RightLeg = 0.0f;
                targetAngles.LeftLeg = 0.0f;

                convergenceRate = 0.1f;

                break;
            }

            case HumanNpcStateType::BehaviorType::Constrained_Walking:
            {
                //
                // Calculate leg angle based on distance traveled
                //

                // Add some dependency on walking speed
                float const actualWalkingSpeed = CalculateActualHumanWalkingAbsoluteSpeed(humanNpcState);
                float const MaxLegAngle =
                    0.41f // std::atan((GameParameters::HumanNpcGeometry::StepLengthFraction / 2.0f) / GameParameters::HumanNpcGeometry::LegLengthFraction)
                    * std::sqrt(actualWalkingSpeed * 0.9f);

                adjustedStandardHumanHeight = humanNpcState.Height * mCurrentHumanNpcBodyLengthAdjustment;
                float const stepLength = GameParameters::HumanNpcGeometry::StepLengthFraction * adjustedStandardHumanHeight;
                float const distance =
                    humanNpcState.TotalDistanceTraveledOnEdgeSinceStateTransition
                    + 0.3f * humanNpcState.TotalDistanceTraveledOffEdgeSinceStateTransition;
                float const distanceInTwoSteps = std::fmod(distance + 3.0f * stepLength / 2.0f, stepLength * 2.0f);

                float const legAngle = std::abs(stepLength - distanceInTwoSteps) / stepLength * 2.0f * MaxLegAngle - MaxLegAngle;

                targetAngles.RightLeg = legAngle;
                targetAngles.LeftLeg = -legAngle;

                // Arms depend on panic
                if (humanNpcState.ResultantPanicLevel < 0.0001f)
                {
                    // Arms aperture depends on speed

                    // At base speed (1m/s): 1.4
                    // Swing more
                    float const apertureMultiplier = 1.4f + (actualWalkingSpeed - 1.0f) * 0.4f;
                    targetAngles.RightArm = targetAngles.LeftLeg * apertureMultiplier;
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
                    targetAngles.RightArm = Pi<float> -angle * panicMultiplier;
                }
                targetAngles.LeftArm = -targetAngles.RightArm;

                convergenceRate = 0.25f;

                if (primaryContrainedState.has_value() && primaryContrainedState->CurrentVirtualFloor.has_value())
                {
                    //
                    // We are walking on an edge - make sure feet don't look weird on sloped edges
                    //

                    ElementIndex const edgeElementIndex =
                        shipMesh.GetTriangles().GetSubSprings(primaryContrainedState->CurrentVirtualFloor->TriangleElementIndex)
                        .SpringIndices[primaryContrainedState->CurrentVirtualFloor->EdgeOrdinal];
                    // Note: we do not care if not in CW order
                    edg1 = shipMesh.GetSprings().GetEndpointAPosition(edgeElementIndex, shipMesh.GetPoints());
                    edg2 = shipMesh.GetSprings().GetEndpointBPosition(edgeElementIndex, shipMesh.GetPoints());
                    edgVector = edg2 - edg1;
                    edgDir = edgVector.normalise_approx();

                    //
                    // 1. Limit leg angles if on slope
                    //

                    vec2f const headPosition = mParticles.GetPosition(secondaryParticleIndex);
                    feetPosition = mParticles.GetPosition(primaryParticleIndex);
                    actualBodyVector = feetPosition - headPosition; // From head to feet
                    actualBodyDir = actualBodyVector.normalise_approx();

                    float const bodyToVirtualEdgeAlignment = std::abs(edgDir.dot(actualBodyDir.to_perpendicular()));
                    float const angleLimitFactor = bodyToVirtualEdgeAlignment * bodyToVirtualEdgeAlignment * bodyToVirtualEdgeAlignment;
                    targetAngles.RightLeg *= angleLimitFactor;
                    targetAngles.LeftLeg *= angleLimitFactor;
                }

                break;
            }

            case HumanNpcStateType::BehaviorType::Constrained_Falling:
            {
                // Both arms in direction of face, depending on head velocity in that direction

                vec2f const headPosition = mParticles.GetPosition(secondaryParticleIndex);
                feetPosition = mParticles.GetPosition(primaryParticleIndex);
                actualBodyVector = feetPosition - headPosition; // From head to feet
                actualBodyDir = actualBodyVector.normalise_approx();

                // The extent to which we move arms depends on the avg velocity or head+feet

                vec2f const & headVelocity = npc.ParticleMesh.Particles[1].GetApplicableVelocity(mParticles);
                vec2f const & feetVelocity = npc.ParticleMesh.Particles[0].GetApplicableVelocity(mParticles);

                float const avgVelocityAlongBodyPerp = (headVelocity + feetVelocity).dot(actualBodyDir.to_perpendicular());
                float const targetDepth = LinearStep(0.0f, 3.0f, std::abs(avgVelocityAlongBodyPerp));

                if (humanNpcState.CurrentFaceDirectionX >= 0.0f)
                {
                    // We want to send arms to the right...
                    // ...but not against face direction
                    if (humanNpcState.CurrentFaceDirectionX >= 0.0f)
                    {
                        targetAngles.RightArm = Pi<float> / 2.0f * targetDepth + 0.09f;
                        targetAngles.LeftArm = targetAngles.RightArm - 0.18f;
                    }
                    else
                    {
                        targetAngles.RightArm = 0.0f;
                        targetAngles.LeftArm = 0.0f;
                    }
                }
                else
                {
                    // We want to send arms to the left...
                   // ...but not against face direction
                    if (humanNpcState.CurrentFaceDirectionX <= 0.0f)
                    {
                        targetAngles.LeftArm = -Pi<float> / 2.0f * targetDepth - 0.09f;
                        targetAngles.RightArm = targetAngles.LeftArm + 0.18f;
                    }
                    else
                    {
                        targetAngles.RightArm = 0.0f;
                        targetAngles.LeftArm = 0.0f;
                    }
                }

                convergenceRate = 0.08f * (0.2f + targetDepth);

                // Close legs
                targetAngles.RightLeg = 0.05f;
                targetAngles.LeftLeg = -0.05f;

                break;
            }

            case HumanNpcStateType::BehaviorType::Constrained_KnockedOut:
            {
                // Check if both head and feet are on a floor
                bool const isHeadOnEdge = primaryContrainedState.has_value() && primaryContrainedState->CurrentVirtualFloor.has_value();
                bool const areFootOnEdge = secondaryConstrainedState.has_value() && secondaryConstrainedState->CurrentVirtualFloor.has_value();
                if (isHeadOnEdge && areFootOnEdge)
                {
                    // Arms: +/- PI or 0, depending on where they are now

                    if (animationState.LimbAngles.RightArm >= -Pi<float> / 2.0f
                        && animationState.LimbAngles.RightArm <= Pi<float> / 2.0f)
                    {
                        targetAngles.RightArm = 0.0f;
                    }
                    else
                    {
                        targetAngles.RightArm = Pi<float>;
                    }

                    if (animationState.LimbAngles.LeftArm >= -Pi<float> / 2.0f
                        && animationState.LimbAngles.LeftArm <= Pi<float> / 2.0f)
                    {
                        targetAngles.LeftArm = 0.0f;
                    }
                    else
                    {
                        targetAngles.LeftArm = Pi<float>;
                    }
                }

                // Legs: 0

                targetAngles.RightLeg = 0.0f;
                targetAngles.LeftLeg = 0.0f;

                convergenceRate = 0.05f; // Quite slow

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
                feetPosition = mParticles.GetPosition(primaryParticleIndex);
                actualBodyVector = feetPosition - headPosition;
                actualBodyDir = actualBodyVector.normalise_approx();

                // Arms: always up, unless horizontal or foot on the floor

                float const horizontality = std::abs(actualBodyDir.dot(GameParameters::GravityDir));

                float const armAngle = (primaryContrainedState.has_value() && primaryContrainedState->CurrentVirtualFloor.has_value())
                    ? Pi<float> / 2.0f
                    : Pi<float> - (Pi<float> / 2.0f) / std::exp(horizontality * 2.2f);
                targetAngles.RightArm = armAngle;
                targetAngles.LeftArm = -targetAngles.RightArm;

                if (humanNpcState.CurrentBehavior != HumanNpcStateType::BehaviorType::Constrained_KnockedOut)
                {
                    // Legs: when arms far from rest, tight; when arms close, at fixed angle

                    float constexpr LegRestAngle = HumanNpcStateType::AnimationStateType::InitialLegAngle;
                    float const legAngle = LegRestAngle;

                    // Legs inclined in direction opposite of relvel, by an amount proportional to relvel itself
                    vec2f const relativeVelocity = mParticles.GetVelocity(primaryParticleIndex) - mParticles.GetVelocity(secondaryParticleIndex);
                    float const relVelPerpToBody = relativeVelocity.dot(actualBodyDir.to_perpendicular());
                    float const legAngleOffset = -SmoothStep(0.0f, 3.0f, std::abs(relVelPerpToBody)) * (LegRestAngle + 0.3f) * (relVelPerpToBody < 0.0f ? -1.0f : 1.0f);
                    targetAngles.RightLeg = legAngle - legAngleOffset;
                    targetAngles.LeftLeg = -legAngle - legAngleOffset;
                }
                else
                {
                    // Leave legs as-is
                }

                convergenceRate = 0.1f;

                break;
            }

            case HumanNpcStateType::BehaviorType::Free_Swimming_Style1:
            {
                //
                // Arms and legs up<->down
                //

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

                float elapsed = currentSimulationTime - humanNpcState.CurrentStateTransitionSimulationTimestamp;
                // Prolong first period
                float constexpr ActualLeadInTime = 6.0f;
                if (elapsed < ActualLeadInTime)
                {
                    elapsed = elapsed / ActualLeadInTime * Period1;
                }
                else
                {
                    elapsed = elapsed - Period1;
                }

                float const panicAccelerator = 1.0f + std::min(humanNpcState.ResultantPanicLevel, 2.0f) / 2.0f * 4.0f;

                float const arg =
                    Period1 / 2.0f // Start some-halfway-through to avoid sudden extreme angles
                    + elapsed * 2.6f * panicAccelerator
                    + humanNpcState.TotalDistanceTraveledOffEdgeSinceStateTransition * 0.7f;

                float const inPeriod = fmod(arg, (Period1 + Period2));
                // y: [0.0 ... 1.0]
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
                targetAngles.RightArm = armAngle;
                targetAngles.LeftArm = -targetAngles.RightArm;

                // Legs:flapping around a (small) angle, which becomes even smaller
                // width depth amplitude depending on depth
                float constexpr LegAngleAmplitude = 0.25f * 2.0f; // Half of this on each side of center angle
                float const legCenterAngle = 0.25f * (depthDamper * 0.5f + 0.5f);
                float const legAngle =
                    legCenterAngle
                    + (y * 2.0f - 1.0f) * LegAngleAmplitude / 2.0f * (depthDamper * 0.35f + 0.65f);
                targetAngles.RightLeg = legAngle;
                targetAngles.LeftLeg = -targetAngles.RightLeg;

                // Convergence rate depends on how long we've been in this state
                float const MaxConvergenceWait = 3.5f;
                convergenceRate = 0.01f + Clamp(elapsed, 0.0f, MaxConvergenceWait) / MaxConvergenceWait * (0.25f - 0.01f);

                break;
            }

            case HumanNpcStateType::BehaviorType::Free_Swimming_Style2:
            {
                //
                // Trappelen
                //

                float constexpr Period = 2.00f;

                float const elapsed = currentSimulationTime - humanNpcState.CurrentStateTransitionSimulationTimestamp;
                float const panicAccelerator = 1.0f + std::min(humanNpcState.ResultantPanicLevel, 2.0f) / 2.0f * 1.0f;

                float const arg =
                    elapsed * 2.6f * panicAccelerator
                    + humanNpcState.TotalDistanceTraveledOffEdgeSinceStateTransition * 0.7f;

                float const inPeriod = fmod(arg, Period);
                // periodicValue: [0.0 ... 1.0]
                periodicValue = (inPeriod < Period / 2.0f)
                    ? inPeriod / (Period / 2.0f)
                    : 1.0f - (inPeriod - (Period / 2.0f)) / (Period / 2.0f);

                // Arms: around a small angle
                targetAngles.RightArm = StateType::KindSpecificStateType::HumanNpcStateType::AnimationStateType::InitialArmAngle + (periodicValue - 0.5f) * Pi<float>/8.0f;
                targetAngles.LeftArm = -targetAngles.RightArm;

                // Legs: perfectly vertical
                targetAngles.RightLeg = 0.0f;
                targetAngles.LeftLeg = 0.0f;

                // Convergence rate depends on how long we've been in this state
                float const MaxConvergenceWait = 3.5f;
                convergenceRate = 0.01f + Clamp(elapsed, 0.0f, MaxConvergenceWait) / MaxConvergenceWait * (0.25f - 0.01f);

                break;
            }

            case HumanNpcStateType::BehaviorType::Free_Swimming_Style3:
            {
                //
                // Trappelen
                //

                float constexpr Period = 2.00f;

                float const elapsed = currentSimulationTime - humanNpcState.CurrentStateTransitionSimulationTimestamp;
                float const panicAccelerator = 1.0f + std::min(humanNpcState.ResultantPanicLevel, 2.0f) / 2.0f * 2.0f;

                float const arg =
                    elapsed * 2.6f * panicAccelerator
                    + humanNpcState.TotalDistanceTraveledOffEdgeSinceStateTransition * 0.7f;

                float const inPeriod = fmod(arg, Period);
                // periodicValue: [0.0 ... 1.0]
                periodicValue = (inPeriod < Period / 2.0f)
                    ? inPeriod / (Period / 2.0f)
                    : 1.0f - (inPeriod - (Period / 2.0f)) / (Period / 2.0f);

                // Arms: one arm around around a large angle; the other fixed around a small angle
                float const angle1 = (Pi<float> -StateType::KindSpecificStateType::HumanNpcStateType::AnimationStateType::InitialArmAngle) + (periodicValue - 0.5f) * Pi<float> / 8.0f;
                float const angle2 = -StateType::KindSpecificStateType::HumanNpcStateType::AnimationStateType::InitialArmAngle;
                if (npc.RandomNormalizedUniformSeed >= 0.0f)
                {
                    targetAngles.RightArm = angle1;
                    targetAngles.LeftArm = angle2;
                }
                else
                {
                    targetAngles.RightArm = -angle2;
                    targetAngles.LeftArm = -angle1;
                }

                // Legs: perfectly vertical
                targetAngles.RightLeg = 0.0f;
                targetAngles.LeftLeg = 0.0f;

                // Convergence rate depends on how long we've been in this state
                float const MaxConvergenceWait = 3.5f;
                convergenceRate = 0.01f + Clamp(elapsed, 0.0f, MaxConvergenceWait) / MaxConvergenceWait * (0.25f - 0.01f);

                break;
            }
        }

        // Converge
        animationState.LimbAngles.ConvergeTo(targetAngles, convergenceRate);
        animationState.UpperLegLengthFraction = targetUpperLegLengthFraction;

        // Calculate sins and coss
        SinCos4(animationState.LimbAngles.fptr(), animationState.LimbAnglesSin.fptr(), animationState.LimbAnglesCos.fptr());

        //
        // Length Multipliers
        //

        FS_ALIGN16_BEG LimbVector targetLengthMultipliers({ 1.0f, 1.0f, 1.0f, 1.0f }) FS_ALIGN16_END;
        float limbLengthConvergenceRate = convergenceRate;

        float targetLowerExtremityLengthMultiplier = 1.0f;

        float constexpr MinPrerisingArmLengthMultiplier = 0.35f;

        switch (humanNpcState.CurrentBehavior)
        {
            case HumanNpcStateType::BehaviorType::Constrained_PreRising:
            {
                // Retract arms
                targetLengthMultipliers.RightArm = MinPrerisingArmLengthMultiplier;
                targetLengthMultipliers.LeftArm = MinPrerisingArmLengthMultiplier;

                break;
            }

            case HumanNpcStateType::BehaviorType::Constrained_Rising:
            {
                if (humanNpcState.CurrentBehaviorState.Constrained_Rising.VirtualEdgeRisingAgainst.TriangleElementIndex != NoneElementIndex) // Locals guaranteed to be calc'd
                {
                    // Recoil arms

                    // For such a small angle, tan(x) ~= x
                    float const targetArmLengthMultiplier =
                        MinPrerisingArmLengthMultiplier
                        + Clamp(humanEdgeAngle / MaxHumanEdgeAngleForArms, 0.0f, 1.0f) * (1.0f - MinPrerisingArmLengthMultiplier);

                    targetLengthMultipliers.RightArm = targetArmLengthMultiplier;
                    targetLengthMultipliers.LeftArm = targetArmLengthMultiplier;
                }
                else
                {
                    // Survive small bursts of losing the edge
                    targetLengthMultipliers.RightArm = animationState.LimbLengthMultipliers.RightArm;
                    targetLengthMultipliers.LeftArm = animationState.LimbLengthMultipliers.LeftArm;
                }

                break;
            }

            case HumanNpcStateType::BehaviorType::Constrained_Walking:
            {
                // Take into account that crotch is lower
                targetLowerExtremityLengthMultiplier = animationState.LimbAnglesCos.RightLeg;

                if (primaryContrainedState.has_value() && primaryContrainedState->CurrentVirtualFloor.has_value())
                {
                    //
                    // We are walking on an edge - make sure feet don't look weird on sloped edges
                    //

                    //
                    // 2. Constrain feet onto edge - i.e. adjust leg lengths
                    //

                    //
                    // Using parametric eq's (tl=scalar from leg1 to leg2, te=scalar from edg1 to edg2):
                    //
                    // leg1 + tl * (leg2 - leg1) = edg1 + te * (edg2 - edg1)
                    // =>
                    // tl = (edg1.y - leg1.y) * (edg2.x - edg1.x) + (leg1.x - edg1.x) * (edg2.y - edg1.y)
                    //      -----------------------------------------------------------------------------
                    //                                  edg X leg
                    //

                    float constexpr MaxLengthMultiplier = 1.4f;

                    float const adjustedStandardLegLength = GameParameters::HumanNpcGeometry::LegLengthFraction * adjustedStandardHumanHeight;
                    vec2f const crotchPosition = feetPosition - actualBodyVector * (GameParameters::HumanNpcGeometry::LegLengthFraction * targetLowerExtremityLengthMultiplier);

                    // leg*1 is crotchPosition
                    float const numerator = (edg1.y - crotchPosition.y) * (edg2.x - edg1.x) + (crotchPosition.x - edg1.x) * (edg2.y - edg1.y);

                    {
                        vec2f const legrVector = actualBodyDir.rotate(animationState.LimbAnglesCos.RightLeg, animationState.LimbAnglesSin.RightLeg) * adjustedStandardLegLength;
                        float const edgCrossRightLeg = edgVector.cross(legrVector);
                        if (std::abs(edgCrossRightLeg) > 0.0000001f)
                        {
                            //targetRightLegLengthMultiplier = numerator / edgCrossRightLeg;
                            float const candidate = numerator / edgCrossRightLeg;
                            if (candidate > 0.01f)
                            {
                                targetLengthMultipliers.RightLeg = std::min(candidate, MaxLengthMultiplier);
                            }
                        }
                    }

                    {
                        vec2f const leglVector = actualBodyDir.rotate(animationState.LimbAnglesCos.LeftLeg, animationState.LimbAnglesSin.LeftLeg) * adjustedStandardLegLength;
                        float const edgCrossLeftLeg = edgVector.cross(leglVector);
                        if (std::abs(edgCrossLeftLeg) > 0.0000001f)
                        {
                            //targetLeftLegLengthMultiplier = numerator / edgCrossLeftLeg;
                            float const candidate = numerator / edgCrossLeftLeg;
                            if (candidate > 0.01f)
                            {
                                targetLengthMultipliers.LeftLeg = std::min(candidate, MaxLengthMultiplier);
                            }
                        }
                    }

                    limbLengthConvergenceRate = 0.09f;
                }

                break;
            }

            case HumanNpcStateType::BehaviorType::Free_Swimming_Style2:
            {
                //
                // Trappelen lengths
                //

                float constexpr TrappelenExtent = 0.3f;
                targetLengthMultipliers.RightLeg = 1.0f - (1.0f - periodicValue) * TrappelenExtent;
                targetLengthMultipliers.LeftLeg = 1.0f - periodicValue * TrappelenExtent;

                break;
            }

            case HumanNpcStateType::BehaviorType::Free_Swimming_Style3:
            {
                //
                // Trappelen lengths
                //

                float constexpr TrappelenExtent = 0.3f;
                targetLengthMultipliers.RightLeg = 1.0f - (1.0f - periodicValue) * TrappelenExtent;
                targetLengthMultipliers.LeftLeg = 1.0f - periodicValue * TrappelenExtent;

                break;
            }

            case HumanNpcStateType::BehaviorType::BeingPlaced:
            case HumanNpcStateType::BehaviorType::Constrained_Equilibrium:
            case HumanNpcStateType::BehaviorType::Constrained_Falling:
            case HumanNpcStateType::BehaviorType::Constrained_KnockedOut:
            case HumanNpcStateType::BehaviorType::Constrained_Aerial:
            case HumanNpcStateType::BehaviorType::Free_Aerial:
            case HumanNpcStateType::BehaviorType::Free_InWater:
            case HumanNpcStateType::BehaviorType::Free_Swimming_Style1:
            {
                // Nop
                break;
            }
        }

        // Converge
        animationState.LimbLengthMultipliers.ConvergeTo(targetLengthMultipliers, limbLengthConvergenceRate);
        animationState.LowerExtremityLengthMultiplier += (targetLowerExtremityLengthMultiplier - animationState.LowerExtremityLengthMultiplier) * convergenceRate;
    }
}

}