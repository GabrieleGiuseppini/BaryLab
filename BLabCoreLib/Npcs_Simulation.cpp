/***************************************************************************************
 * Original Author:		Gabriele Giuseppini
 * Created:				2023-07-23
 * Copyright:			Gabriele Giuseppini  (https://github.com/GabrieleGiuseppini)
 ***************************************************************************************/
#include "Npcs.h"

#include "BLabMath.h"
#include "Log.h"

#include <limits>

void Npcs::UpdateNpcs(
    float currentSimulationTime,
    Mesh const & mesh,
    LabParameters const & labParameters)
{
    LogMessage("----------------------------------");
    LogMessage("----------------------------------");
    LogMessage("----------------------------------");

    //
    // 1. Reset buffers
    //

    // Note: no need to reset PreliminaryForces as we'll recalculate all of them

    //
    // 2. Check if a free secondary particle should become constrained
    // 3. Calculate preliminary forces 
    //

    for (auto const n : *this)
    {
        auto & npcState = mStateBuffer[n];

        // Secondary free becoming constrained

        if (npcState.DipoleState.has_value()
            && !npcState.DipoleState->SecondaryParticleState.ConstrainedState.has_value() // Secondary is free
            && npcState.PrimaryParticleState.ConstrainedState.has_value()) // And primary is constrained
        {
            npcState.DipoleState->SecondaryParticleState.ConstrainedState = CalculateParticleConstrainedState(
                mParticles.GetPosition(npcState.DipoleState->SecondaryParticleState.ParticleIndex),
                mesh);
        }

        // Preliminary Forces

        CalculateNpcParticlePreliminaryForces(
            npcState,
            true,
            labParameters);

        if (npcState.DipoleState.has_value())
        {
            CalculateNpcParticlePreliminaryForces(
                npcState,
                false,
                labParameters);
        }
    }

    //
    // 4. Update state
    //

    for (auto const n : *this)
    {
        LogMessage("NPC ", n);

        auto & npcState = mStateBuffer[n];

        UpdateNpcParticle(
            npcState,
            true,
            mesh,
            labParameters);

        if (npcState.DipoleState.has_value())
        {
            UpdateNpcParticle(
                npcState,
                false,
                mesh,
                labParameters);
        }
    }

    //
    // 5. Update behavioral state machines
    //

    mParticles.ResetEquilibriumTorque();

    for (auto const n : *this)
    {
        auto & npcState = mStateBuffer[n];

        if (npcState.Type == NpcType::Human)
        {
            assert(npcState.DipoleState.has_value());
            assert(npcState.HumanNpcState.has_value());

            UpdateHuman(
                currentSimulationTime,
                *npcState.HumanNpcState,
                npcState.PrimaryParticleState,
                npcState.DipoleState->SecondaryParticleState,
                mesh,
                labParameters);
        }
    }

    //
    // 6. Update animation
    //

    for (auto const n : *this)
    {
        auto & npcState = mStateBuffer[n];

        UpdateNpcAnimation(
            npcState,
            true,
            mesh,
            labParameters);

        if (npcState.DipoleState.has_value())
        {
            UpdateNpcAnimation(
                npcState,
                false,
                mesh,
                labParameters);
        }
    }

    // TODOTEST
    if (mStateBuffer[0].PrimaryParticleState.ConstrainedState.has_value())
        mEventDispatcher.OnCustomProbe("CurrentVirtEdge", float(mStateBuffer[0].PrimaryParticleState.ConstrainedState->CurrentVirtualEdgeElementIndex));
    else
        mEventDispatcher.OnCustomProbe("CurrentVirtEdge", -1.0f);
}

void Npcs::UpdateNpcParticle(
    StateType & npc,
    bool isPrimaryParticle,
    Mesh const & mesh,
    LabParameters const & labParameters)
{
    auto & npcParticle = isPrimaryParticle ? npc.PrimaryParticleState : npc.DipoleState->SecondaryParticleState;

    LogMessage("----------------------------------");
    LogMessage("  ", isPrimaryParticle ? "Primary" : "Secondary");

    float const particleMass = mParticles.GetPhysicalProperties(npcParticle.ParticleIndex).Mass * labParameters.MassAdjustment;
    float const dt = LabParameters::SimulationTimeStepDuration;

    vec2f const particleStartAbsolutePosition = mParticles.GetPosition(npcParticle.ParticleIndex);

    // Calculate physical displacement - once and for all, as whole loop
    // will attempt to move to trajectory that always ends here

    vec2f physicsDeltaPos;
    if (mCurrentParticleTrajectory.has_value() && npcParticle.ParticleIndex == mCurrentParticleTrajectory->ParticleIndex)
    {
        // Consume externally-supplied trajectory

        physicsDeltaPos = mCurrentParticleTrajectory->TargetPosition - particleStartAbsolutePosition;
        mCurrentParticleTrajectory.reset();
    }
    else
    {
        // Integrate forces

        vec2f const physicalForces = CalculateNpcParticleDefinitiveForces(
            npc,
            isPrimaryParticle,
            particleMass,
            labParameters);

        physicsDeltaPos = mParticles.GetVelocity(npcParticle.ParticleIndex) * dt + (physicalForces / particleMass) * dt * dt;
    }

    if (!npcParticle.ConstrainedState.has_value())
    {
        // 
        // Particle is free
        //

        LogMessage("    Free: velocity=", mParticles.GetVelocity(npcParticle.ParticleIndex), " prelimF=", mParticles.GetPreliminaryForces(npcParticle.ParticleIndex), " physicsDeltaPos=", physicsDeltaPos);
        LogMessage("    StartPosition=", particleStartAbsolutePosition, " StartVelocity=", mParticles.GetVelocity(npcParticle.ParticleIndex));

        UpdateNpcParticle_Free(
            npcParticle,
            particleStartAbsolutePosition,
            particleStartAbsolutePosition + physicsDeltaPos,
            mParticles,
            labParameters);

        LogMessage("    EndPosition=", mParticles.GetPosition(npcParticle.ParticleIndex), " EndVelocity=", mParticles.GetVelocity(npcParticle.ParticleIndex));

        // We're done
    }
    else
    {
        //
        // Constrained
        //

        assert(npcParticle.ConstrainedState.has_value());

        LogMessage("    Constrained: velocity=", mParticles.GetVelocity(npcParticle.ParticleIndex), " prelimF=", mParticles.GetPreliminaryForces(npcParticle.ParticleIndex), " physicsDeltaPos=", physicsDeltaPos);

        // Loop tracing trajectory from TrajectoryStart (== current bary coords) to TrajectoryEnd (== start absolute pos + deltaPos); 
        // each step moves the next TrajectoryStart a bit ahead.
        // Each iteration of the loop either exits (completes), or moves current bary coords (and calcs remaining dt) when it wants 
        // to "continue" an impact while on edge-moving-against-it, i.e. when it wants to recalculate a new flattened traj.
        //    - In this case, the iteration doesn't change current absolute position, nor velocity; it only updates current bary coords 
        //      to what will become the next TrajectoryStart
        //    - In this case, at next iteration:
        //          - TrajectoryStart (== current bary coords) is new
        //          - TrajectoryEnd (== start absolute pos + physicsDeltaPos) is same as before

        // Initialize absolute position of particle (wrt current triangle) as if it moved with the mesh, staying in its position wrt its triangle;
        // it's the new theoretical position after mesh displacement
        vec2f trajectoryStartAbsolutePosition = mesh.GetTriangles().FromBarycentricCoordinates(
            npcParticle.ConstrainedState->CurrentTriangleBarycentricCoords,
            npcParticle.ConstrainedState->CurrentTriangle,
            mesh.GetVertices());

        // Calculate mesh displacement for the whole loop via pure displacement of triangle containing particle
        vec2f const meshVelocity = (particleStartAbsolutePosition - trajectoryStartAbsolutePosition) / LabParameters::SimulationTimeStepDuration;

        // Machinery to detect 3-iteration paths that doesn't move particle (positional well)

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
        };

        std::optional<PastBarycentricPosition> pastPastBarycentricPosition;
        std::optional<PastBarycentricPosition> pastBarycentricPosition;

        // Total displacement walked along the edge - as a sum of the vectors from the individual steps
        vec2f totalEdgeWalkedActual = vec2f::zero();

        for (float remainingDt = dt; ; )
        {
            assert(remainingDt > 0.0f);

            LogMessage("    ------------------------");
            LogMessage("    Triangle=", npcParticle.ConstrainedState->CurrentTriangle, " B-Coords=", npcParticle.ConstrainedState->CurrentTriangleBarycentricCoords, " RemainingDt=", remainingDt);

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

            // Absolute position of particle if it only moved due to physical forces and walking;
            vec2f const trajectoryEndAbsolutePosition = particleStartAbsolutePosition + physicsDeltaPos + totalEdgeWalkedActual;

            // Trajectory
            // Note: on first iteration this is the same as physicsDeltaPos + meshDisplacement
            vec2f const trajectory = trajectoryEndAbsolutePosition - trajectoryStartAbsolutePosition;

            LogMessage("    TrajectoryStartAbsolutePosition=", trajectoryStartAbsolutePosition, " PhysicsDeltaPos=", physicsDeltaPos, " TotalEdgeWalkedActual=", totalEdgeWalkedActual,
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
                bcoords3f const trajectoryEndBarycentricCoords = mesh.GetTriangles().ToBarycentricCoordinates(
                    trajectoryEndAbsolutePosition,
                    npcParticle.ConstrainedState->CurrentTriangle,
                    mesh.GetVertices());

                auto const outcome = NavigateVertex(
                    npc,
                    isPrimaryParticle,
                    *vertexOrdinal,
                    particleStartAbsolutePosition,
                    trajectoryStartAbsolutePosition,
                    trajectoryEndAbsolutePosition,
                    trajectoryEndBarycentricCoords,
                    true, // Check immediately whether we're directed towards the interior of the triangle
                    mParticles,
                    mesh,
                    labParameters);

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
                        // Handle impact with "flattening" algorithm
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

            LogMessage("    Checking edges");

            // Check all edges, stop at first one that is floor and against which we're moving
            bool hasFoundNonInertial = false;
            for (int vi = 0; vi < 3; ++vi)
            {
                if (npcParticle.ConstrainedState->CurrentTriangleBarycentricCoords[vi] == 0.0f)
                {
                    // We are on this edge

                    int const edgeOrdinal = (vi + 1) % 3;
                    ElementIndex const currentEdgeElementIndex = mesh.GetTriangles().GetSubEdges(npcParticle.ConstrainedState->CurrentTriangle).EdgeIndices[edgeOrdinal];

                    assert(isPrimaryParticle || npc.DipoleState.has_value());

                    LogMessage("      edge ", edgeOrdinal, ": isFloor=", IsEdgeFloorToParticle(currentEdgeElementIndex, npcParticle.ConstrainedState->CurrentTriangle, mesh));

                    // Check if this is really a floor to this particle
                    if (IsEdgeFloorToParticle(currentEdgeElementIndex, npcParticle.ConstrainedState->CurrentTriangle, mesh)
                        && (isPrimaryParticle || !DoesFloorSeparateFromPrimaryParticle(
                            mParticles.GetPosition(npc.DipoleState->SecondaryParticleState.ParticleIndex),
                            trajectoryStartAbsolutePosition, // Current (virtual, not yet real) position of this (secondary) particle
                            currentEdgeElementIndex,
                            mesh)))
                    {
                        // On floor edge - so potentially in a non-inertial frame

                        // Check now whether we're moving *against* the floor

                        vec2f const edgeVector = mesh.GetTriangles().GetSubEdgeVector(
                            npcParticle.ConstrainedState->CurrentTriangle,
                            edgeOrdinal,
                            mesh.GetVertices());

                        vec2f const edgeDir = edgeVector.normalise();
                        vec2f const edgeNormal = edgeDir.to_perpendicular(); // Points outside of triangle (i.e. towards floor)

                        if (trajectory.dot(edgeNormal) >= 0.0f) // If 0, no normal force - hence no friction; however we want to take this codepath anyways for consistency
                        {
                            //
                            // Case 1: Non-inertial: on edge and moving against (or along) it, pushed by it
                            //

                            LogMessage("    ConstrainedNonInertial: triangle=", npcParticle.ConstrainedState->CurrentTriangle, " edgeOrdinal=", edgeOrdinal, " bCoords=", npcParticle.ConstrainedState->CurrentTriangleBarycentricCoords, " trajectory=", trajectory);
                            LogMessage("    StartPosition=", mParticles.GetPosition(npcParticle.ParticleIndex), " StartVelocity=", mParticles.GetVelocity(npcParticle.ParticleIndex), " MeshVelocity=", meshVelocity, " StartMRVelocity=", npcParticle.ConstrainedState->MeshRelativeVelocity);

                            npcParticle.ConstrainedState->CurrentVirtualEdgeElementIndex = currentEdgeElementIndex;

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
                                    frictionCoefficient = mParticles.GetPhysicalProperties(npcParticle.ParticleIndex).KineticFriction * labParameters.KineticFrictionAdjustment;
                                }
                                else
                                {
                                    // Static friction
                                    frictionCoefficient = mParticles.GetPhysicalProperties(npcParticle.ParticleIndex).StaticFriction * labParameters.StaticFrictionAdjustment;
                                }

                                // Calculate friction (integrated) force magnitude (along edgeDir),
                                // which is the same as apparent force, up to max friction threshold
                                float tFriction = std::min(std::abs(trajectoryT), frictionCoefficient * std::max(trajectoryN, 0.0f));
                                if (trajectoryT >= 0.0f)
                                {
                                    tFriction *= -1.0f;
                                }

                                LogMessage("        friction: trajectoryN=", trajectoryN, " relVel=", npcParticle.ConstrainedState->MeshRelativeVelocity,
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

                            if (npc.HumanNpcState.has_value() 
                                && npc.HumanNpcState->CurrentBehavior == StateType::HumanNpcStateType::BehaviorType::Constrained_Walking
                                && isPrimaryParticle)
                            {
                                //
                                // Walking displacement projected along edge - what's needed to reach total desired walking
                                // displacement, together with physical
                                // - Never reduces physical
                                // - Never adds to physical so much as to cause resultant to be faster than walking speed
                                //

                                vec2f const idealWalkDir = vec2f(npc.HumanNpcState->CurrentFaceDirectionX, 0.0f);
                                assert(idealWalkDir.length() == 1.0f);

                                float const idealWalkMagnitude =
                                    labParameters.HumanNpcWalkingSpeed
                                    * remainingDt
                                    * npc.HumanNpcState->CurrentBehaviorState.Constrained_Walking.CurrentWalkMagnitude;

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
                                float const gravityResistance = LinearStep(
                                    0.60f, // Start slightly before expected 45-degree ramp
                                    0.79f,
                                    walkDir.dot(-LabParameters::GravityDir)); // TODO: perf: this is just y

                                edgeWalkedPlanned *= (1.0f - gravityResistance);

                                if (npc.HumanNpcState->CurrentBehaviorState.Constrained_Walking.CurrentWalkMagnitude != 0.0f)
                                {
                                    LogMessage("        idealWalkMagnitude=", idealWalkMagnitude, " gravityResistance=", gravityResistance, " => edgeWalkedPlanned=", edgeWalkedPlanned, " (@", npc.HumanNpcState->CurrentBehaviorState.Constrained_Walking.CurrentWalkMagnitude, ")");
                                }
                            }

                            //
                            // Recover flattened trajectory as a vector
                            //

                            float const edgeTraveledPlanned = edgePhysicalTraveledPlanned + edgeWalkedPlanned; // Resultant

                            vec2f const flattenedTrajectory = edgeDir * edgeTraveledPlanned;

                            //
                            // Calculate trajectory target
                            //

                            vec2f flattenedTrajectoryEndAbsolutePosition = trajectoryStartAbsolutePosition + flattenedTrajectory;

                            bcoords3f flattenedTrajectoryEndBarycentricCoords = mesh.GetTriangles().ToBarycentricCoordinates(
                                flattenedTrajectoryEndAbsolutePosition,
                                npcParticle.ConstrainedState->CurrentTriangle,
                                mesh.GetVertices());

                            //
                            // Due to numerical slack, ensure target barycentric coords are still along edge
                            //

                            int const edgeVertexOrdinal = (edgeOrdinal + 2) % 3;
                            flattenedTrajectoryEndBarycentricCoords[edgeVertexOrdinal] = 0.0f;
                            flattenedTrajectoryEndBarycentricCoords[(edgeVertexOrdinal + 1) % 3] = 1.0f - flattenedTrajectoryEndBarycentricCoords[(edgeVertexOrdinal + 2) % 3];

                            LogMessage("        flattenedTrajectory=", flattenedTrajectory, " flattenedTrajectoryEndAbsolutePosition=", flattenedTrajectoryEndAbsolutePosition, " flattenedTrajectoryEndBarycentricCoords=", flattenedTrajectoryEndBarycentricCoords);

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

                            LogMessage("        edgePhysicalTraveledPlanned=", edgePhysicalTraveledPlanned, " edgeWalkedPlanned=", edgeWalkedPlanned);

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
                                edgeTraveledPlanned,
                                meshVelocity,
                                remainingDt,
                                mParticles,
                                mesh,
                                labParameters);

                            LogMessage("    Actual edge traveled in non-inertial step: ", edgeTraveledActual);

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
                                    if (pastPastBarycentricPosition.has_value()
                                        && pastPastBarycentricPosition->TriangleElementIndex == npcParticle.ConstrainedState->CurrentTriangle
                                        && pastPastBarycentricPosition->BarycentricCoords == npcParticle.ConstrainedState->CurrentTriangleBarycentricCoords)
                                    {
                                        //
                                        // Well - stop here
                                        //

                                        LogMessage("    Detected well - stopping here");

                                        // Update particle's physics, considering that we are in a well and thus still (wrt mesh)

                                        vec2f const particleEndAbsolutePosition = mesh.GetTriangles().FromBarycentricCoordinates(
                                            npcParticle.ConstrainedState->CurrentTriangleBarycentricCoords,
                                            npcParticle.ConstrainedState->CurrentTriangle,
                                            mesh.GetVertices());

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
                                    assert(edgeTraveledActual * edgeTraveledPlanned >= 0.0f); // Should have same sign
                                    float const dtFractionConsumed = edgeTraveledPlanned != 0.0f
                                        ? std::min(edgeTraveledActual / edgeTraveledPlanned, 1.0f) // Signs should agree anyway
                                        : 1.0f; // If we were planning no travel, any movement is a whole consumption
                                    LogMessage("        dtFractionConsumed=", dtFractionConsumed);
                                    remainingDt *= (1.0f - dtFractionConsumed);

                                    // Reset well detection machinery
                                    pastPastBarycentricPosition.reset();
                                    pastBarycentricPosition.reset();
                                }

                            }

                            // Update total (edge) traveled
                            if (npc.HumanNpcState.has_value() && isPrimaryParticle)
                            {
                                npc.HumanNpcState->TotalDistanceTraveledSinceStateTransition += std::abs(edgeTraveledActual);
                            }

                            // Update total vector walked along edge
                            float const edgeWalkedActual = edgeTraveledPlanned != 0.0f
                                ? edgeTraveledActual * (edgeWalkedPlanned / edgeTraveledPlanned)
                                : 0.0f; // Unlikely, but read above for rationale behind 0.0
                            totalEdgeWalkedActual += edgeDir * edgeWalkedActual;
                            LogMessage("        edgeWalkedActual=", edgeWalkedActual, " totalEdgeWalkedActual=", totalEdgeWalkedActual);

                            if (npcParticle.ConstrainedState.has_value())
                            {
                                LogMessage("    EndPosition=", mParticles.GetPosition(npcParticle.ParticleIndex), " EndVelocity=", mParticles.GetVelocity(npcParticle.ParticleIndex), " EndMRVelocity=", npcParticle.ConstrainedState->MeshRelativeVelocity);

                                // Update well detection machinery
                                pastPastBarycentricPosition = pastBarycentricPosition;
                                pastBarycentricPosition.emplace(npcParticle.ConstrainedState->CurrentTriangle, npcParticle.ConstrainedState->CurrentTriangleBarycentricCoords);
                            }
                            else
                            {
                                LogMessage("    EndPosition=", mParticles.GetPosition(npcParticle.ParticleIndex), " EndVelocity=", mParticles.GetVelocity(npcParticle.ParticleIndex));
                            }

                            hasFoundNonInertial = true;

                            break;
                        }
                        else
                        {
                            LogMessage("      traj.edgeN=", trajectory.dot(edgeNormal));
                        }
                    }
                }
            } // for (vertex)

            if (!hasFoundNonInertial)
            {
                //
                // Case 2: Inertial: not on edge or on edge but not moving against (nor along) it
                //

                LogMessage("    ConstrainedInertial: triangle=", npcParticle.ConstrainedState->CurrentTriangle, " bCoords=", npcParticle.ConstrainedState->CurrentTriangleBarycentricCoords, " physicsDeltaPos=", physicsDeltaPos);
                LogMessage("    StartPosition=", mParticles.GetPosition(npcParticle.ParticleIndex), " StartVelocity=", mParticles.GetVelocity(npcParticle.ParticleIndex), " MeshVelocity=", meshVelocity, " StartMRVelocity=", npcParticle.ConstrainedState->MeshRelativeVelocity);

                npcParticle.ConstrainedState->CurrentVirtualEdgeElementIndex = NoneElementIndex;

                //
                // Calculate target barycentric coords
                //

                bcoords3f const trajectoryEndBarycentricCoords = mesh.GetTriangles().ToBarycentricCoordinates(
                    trajectoryEndAbsolutePosition,
                    npcParticle.ConstrainedState->CurrentTriangle,
                    mesh.GetVertices());

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
                    mParticles,
                    mesh,
                    labParameters);

                // Update total traveled
                if (npc.HumanNpcState.has_value() && isPrimaryParticle)
                {
                    npc.HumanNpcState->TotalDistanceTraveledSinceStateTransition += std::abs(totalTraveled);
                }

                if (npcParticle.ConstrainedState.has_value())
                {
                    LogMessage("    EndPosition=", mParticles.GetPosition(npcParticle.ParticleIndex), " EndVelocity=", mParticles.GetVelocity(npcParticle.ParticleIndex), " EndMRVelocity=", npcParticle.ConstrainedState->MeshRelativeVelocity);
                }
                else
                {
                    LogMessage("    EndPosition=", mParticles.GetPosition(npcParticle.ParticleIndex), " EndVelocity=", mParticles.GetVelocity(npcParticle.ParticleIndex));
                }

                // Consume whole time quantum
                remainingDt = 0.0f;
            }

            if (remainingDt <= 0.0f)
            {
                assert(remainingDt > -0.0001f); // If negative, only because of numerical slack

                // Consumed whole time quantum, loop completed
                break;
            }

            //
            // Update trajectory start position for next iteration
            //

            // Current (virtual, not yet real) position of this particle
            trajectoryStartAbsolutePosition = mesh.GetTriangles().FromBarycentricCoordinates(
                npcParticle.ConstrainedState->CurrentTriangleBarycentricCoords,
                npcParticle.ConstrainedState->CurrentTriangle,
                mesh.GetVertices());
        }
    }

    if (mCurrentlySelectedParticle == npcParticle.ParticleIndex)
    {
        // Publish final velocities

        vec2f const particleVelocity = (mParticles.GetPosition(npcParticle.ParticleIndex) - particleStartAbsolutePosition) / LabParameters::SimulationTimeStepDuration;

        mEventDispatcher.OnCustomProbe("VelX", particleVelocity.x);
        mEventDispatcher.OnCustomProbe("VelY", particleVelocity.y);
    }
}

void Npcs::CalculateNpcParticlePreliminaryForces(
    StateType const & npc,
    bool isPrimaryParticle,
    LabParameters const & labParameters)
{
    auto & npcParticle = isPrimaryParticle ? npc.PrimaryParticleState : npc.DipoleState->SecondaryParticleState;

    //
    // Calculate world forces
    //

    float const particleMass = mParticles.GetPhysicalProperties(npcParticle.ParticleIndex).Mass * labParameters.MassAdjustment;

    // 1. World forces - gravity

    vec2f preliminaryForces = LabParameters::Gravity * labParameters.GravityAdjustment * mGravityGate * particleMass;

    if (!npcParticle.ConstrainedState.has_value())
    {
        // Check whether we are underwater

        float constexpr BuoyancyInterfaceWidth = 0.4f;

        vec2f testParticlePosition = mParticles.GetPosition(npcParticle.ParticleIndex);
        if (npc.Type == NpcType::Human && !isPrimaryParticle)
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
                LabParameters::GravityMagnitude * 1000.0f
                * mParticles.GetPhysicalProperties(npcParticle.ParticleIndex).BuoyancyVolumeFill
                * labParameters.BuoyancyAdjustment
                * uwCoefficient;

            // 3. World forces - water drag

            preliminaryForces +=
                -mParticles.GetVelocity(npcParticle.ParticleIndex)
                * LabParameters::WaterFrictionDragCoefficient
                * labParameters.WaterFrictionDragCoefficientAdjustment;
        }
    }

    // 3. Spring forces

    if (npc.DipoleState.has_value())
    {
        ElementIndex const otherParticleIndex = isPrimaryParticle ? npc.DipoleState->SecondaryParticleState.ParticleIndex : npc.PrimaryParticleState.ParticleIndex;

        float const dt = LabParameters::SimulationTimeStepDuration;

        vec2f const springDisplacement = mParticles.GetPosition(otherParticleIndex) - mParticles.GetPosition(npcParticle.ParticleIndex); // Towards other
        float const springDisplacementLength = springDisplacement.length();
        vec2f const springDir = springDisplacement.normalise_approx(springDisplacementLength);

        //
        // 3a. Hooke's law
        //

        float const springStiffnessCoefficient =
            labParameters.SpringReductionFraction
            * npc.DipoleState->DipoleProperties.MassFactor * labParameters.MassAdjustment
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
            labParameters.SpringDampingCoefficient
            * npc.DipoleState->DipoleProperties.MassFactor * labParameters.MassAdjustment
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
    LabParameters const & labParameters) const
{
    auto & npcParticle = isPrimaryParticle ? npc.PrimaryParticleState : npc.DipoleState->SecondaryParticleState;
    
    vec2f definitiveForces = mParticles.GetPreliminaryForces(npcParticle.ParticleIndex);

    //
    // Human Equlibrium Torque
    //

    if (npc.HumanNpcState.has_value() 
        && (npc.HumanNpcState->CurrentBehavior == StateType::HumanNpcStateType::BehaviorType::Constrained_Walking 
            || npc.HumanNpcState->CurrentBehavior == StateType::HumanNpcStateType::BehaviorType::Constrained_Equilibrium
            || npc.HumanNpcState->CurrentBehavior == StateType::HumanNpcStateType::BehaviorType::Constrained_Rising)
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
            + mParticles.GetVelocity(secondaryParticleIndex) * LabParameters::SimulationTimeStepDuration;
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
        float const staticDisplacementAngleCW = (-LabParameters::GravityDir).angleCw(humanVector);

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
            * LabParameters::SimulationTimeStepDuration;
        float const relativeVelocityAngleCW = humanVector.angleCw(humanVector + relativeVelocityDisplacement);

        //
        // Calculate the torque on the secondary (head)
        // required to maintain alignment with vertical
        //

        // Calculate angle that we want to enforce with this torque
        float const totalTorqueAngleCW =
            staticDisplacementAngleCW * labParameters.HumanNpcEquilibriumTorqueStiffnessCoefficient
            + relativeVelocityAngleCW * labParameters.HumanNpcEquilibriumTorqueDampingCoefficient;

        // Calculate (linear) force that generates this rotation
        vec2f const torqueDisplacement = humanVector.rotate(totalTorqueAngleCW) - humanVector;
        vec2f const equilibriumTorqueForce =
            torqueDisplacement
            * particleMass / (LabParameters::SimulationTimeStepDuration * LabParameters::SimulationTimeStepDuration)
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
    LabParameters const & labParameters) const
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
        (endPosition - startPosition) / LabParameters::SimulationTimeStepDuration * (1.0f - labParameters.GlobalDamping));
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
    NpcParticles & particles,
    Mesh const & mesh,
    LabParameters const & labParameters) const
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
        LogMessage("      Target is on/in triangle, moving to target");

        //
        // Update particle and exit - consuming whole time quantum
        //            

        npcParticleConstrainedState.CurrentTriangleBarycentricCoords = flattenedTrajectoryEndBarycentricCoords;

        vec2f const particleEndAbsolutePosition = mesh.GetTriangles().FromBarycentricCoordinates(
            flattenedTrajectoryEndBarycentricCoords,
            npcParticleConstrainedState.CurrentTriangle,
            mesh.GetVertices());

        particles.SetPosition(npcParticle.ParticleIndex, particleEndAbsolutePosition);

        //
        // Velocity: given that we've completed *along the edge*, then we can calculate
        // our (relative) velocity based on the distance traveled along this time quantum
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

        LogMessage("        edgeTraveledPlanned=", edgeTraveledPlanned, " absoluteVelocity=", particles.GetVelocity(npcParticle.ParticleIndex));

        // Complete
        return std::make_tuple(
            edgeTraveledPlanned, // We moved by how much we planned
            true); // Stop here
    }

    //
    // Target is outside triangle
    //

    LogMessage("      Target is outside triangle");

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

    LogMessage("      Moving to intersection vertex ", intersectionVertexOrdinal, ": ", intersectionBarycentricCoords);

    npcParticleConstrainedState.CurrentTriangleBarycentricCoords = intersectionBarycentricCoords;

    // Update (signed) edge traveled, assuming we're always moving along initial edge
    // (other edges we're only passing through cuspids)
    vec2f const intersectionAbsolutePosition = mesh.GetVertices().GetPosition(mesh.GetTriangles().GetVertexIndices(npcParticleConstrainedState.CurrentTriangle)[intersectionVertexOrdinal]);
    assert(intersectionAbsolutePosition == mesh.GetTriangles().FromBarycentricCoordinates(
        intersectionBarycentricCoords,
        npcParticleConstrainedState.CurrentTriangle,
        mesh.GetVertices()));

    float const edgeTraveled = (intersectionAbsolutePosition - trajectoryStartAbsolutePosition).dot(edgeDir);

    LogMessage("        edgeTraveled=", edgeTraveled);

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
        particles,
        mesh,
        labParameters);

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
                mesh.GetTriangles().GetSubEdgeVector(npcParticleConstrainedState.CurrentTriangle, navigationOutcome.EncounteredFloorEdgeOrdinal, mesh.GetVertices())
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

                LogMessage("      Impact continuation (trajProjOntoEdgeNormal=", trajProjOntoEdgeNormal, ")");

                return std::make_tuple(edgeTraveled, false);
            }
            else
            {
                //
                // Bounce - calculate bounce response, using the *apparent* (trajectory) 
                // velocity - since this one includes the mesh velocity
                //

                LogMessage("      Bounce (trajProjOntoEdgeNormal=", trajProjOntoEdgeNormal, ")");

                BounceConstrainedNpcParticle(
                    npc,
                    isPrimaryParticle,
                    flattenedTrajectory,
                    intersectionAbsolutePosition,
                    floorEdgeNormal,
                    meshVelocity,
                    dt,
                    particles,
                    labParameters);

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
    NpcParticles & particles,
    Mesh const & mesh,
    LabParameters const & labParameters) const
{
    auto & npcParticle = isPrimaryParticle ? npc.PrimaryParticleState : npc.DipoleState->SecondaryParticleState;
    assert(npcParticle.ConstrainedState.has_value());
    auto & npcParticleConstrainedState = *npcParticle.ConstrainedState;

    vec2f const segmentTrajectoryStartAbsolutePosition = mesh.GetTriangles().FromBarycentricCoordinates(
        segmentTrajectoryStartBarycentricCoords,
        npcParticle.ConstrainedState->CurrentTriangle,
        mesh.GetVertices());

    //
    // Ray-trace along the specified trajectory, ending only at one of the following three conditions:
    // 1. Reached destination: terminate
    // 2. Becoming free: do free movement and terminate
    // 3. Impact with bounce: impart bounce velocity and terminate
    //

    for (int iIter = 0; ; ++iIter)
    {
        assert(npcParticleConstrainedState.CurrentTriangleBarycentricCoords.is_on_edge_or_internal());

        LogMessage("    SegmentTrace ", iIter);
        LogMessage("      triangle=", npcParticleConstrainedState.CurrentTriangle, " bCoords=", npcParticleConstrainedState.CurrentTriangleBarycentricCoords,
            " segmentTrajStartBCoords=", segmentTrajectoryStartBarycentricCoords, " segmentTrajEndBCoords=", segmentTrajectoryEndBarycentricCoords);

        //
        // If target is on/in triangle, we move to target
        //

        if (segmentTrajectoryEndBarycentricCoords.is_on_edge_or_internal())
        {
            LogMessage("      Target is on/in triangle, moving to target");

            //
            // Update particle and exit - consuming whole time quantum
            //            

            // Move particle to end of trajectory
            npcParticleConstrainedState.CurrentTriangleBarycentricCoords = segmentTrajectoryEndBarycentricCoords;
            particles.SetPosition(npcParticle.ParticleIndex, segmentTrajectoryEndAbsolutePosition);

            // Use whole time quantum for velocity, as particleStartAbsolutePosition is fixed at t0
            vec2f const totalAbsoluteTraveledVector = segmentTrajectoryEndAbsolutePosition - particleStartAbsolutePosition;
            vec2f const absoluteVelocity = totalAbsoluteTraveledVector / LabParameters::SimulationTimeStepDuration;
            particles.SetVelocity(npcParticle.ParticleIndex, absoluteVelocity);
            npcParticleConstrainedState.MeshRelativeVelocity = absoluteVelocity + meshVelocity;

            LogMessage("        totalAbsoluteTraveledVector=", totalAbsoluteTraveledVector, " absoluteVelocity=", particles.GetVelocity(npcParticle.ParticleIndex));

            // Return (mesh-relative) distance traveled with this move
            return (segmentTrajectoryEndAbsolutePosition - segmentTrajectoryStartAbsolutePosition).length();
        }

        //
        // We're inside or on triangle, and target is outside triangle;
        // if we're on edge, trajectory is along this edge
        //

        LogMessage("      Target is outside triangle");

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

                LogMessage("        t[v", vi, " e", edgeOrdinal, "] = ", t);

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
        ElementIndex const intersectionEdgeElementIndex = mesh.GetTriangles().GetSubEdges(npcParticleConstrainedState.CurrentTriangle).EdgeIndices[intersectionEdgeOrdinal];

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

        LogMessage("      Moving bary coords to intersection with edge ", intersectionEdgeOrdinal, " ", intersectionBarycentricCoords);

        npcParticleConstrainedState.CurrentTriangleBarycentricCoords = intersectionBarycentricCoords;

        //
        // Check if impacted with floor
        //

        vec2f const intersectionAbsolutePosition = mesh.GetTriangles().FromBarycentricCoordinates(
            intersectionBarycentricCoords,
            npcParticleConstrainedState.CurrentTriangle,
            mesh.GetVertices());

        assert(isPrimaryParticle || npc.DipoleState.has_value());

        if (IsEdgeFloorToParticle(intersectionEdgeElementIndex, npcParticleConstrainedState.CurrentTriangle, mesh)
            && (isPrimaryParticle || !DoesFloorSeparateFromPrimaryParticle(
                mParticles.GetPosition(npc.DipoleState->SecondaryParticleState.ParticleIndex),
                intersectionAbsolutePosition,
                intersectionEdgeElementIndex,
                mesh)))
        {
            //
            // Impact and bounce
            //

            LogMessage("      Impact and bounce");

            //
            // Calculate bounce response, using the *apparent* (trajectory) 
            // velocity - since this one includes the mesh velocity
            //

            vec2f const trajectory = segmentTrajectoryEndAbsolutePosition - segmentTrajectoryStartAbsolutePosition;

            vec2f const intersectionEdgeDir =
                mesh.GetTriangles().GetSubEdgeVector(npcParticleConstrainedState.CurrentTriangle, intersectionEdgeOrdinal, mesh.GetVertices())
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
                labParameters);

            // Return (mesh-relative) distance traveled with this move
            return (segmentTrajectoryEndAbsolutePosition - segmentTrajectoryStartAbsolutePosition).length();
        }
        else
        {
            //
            // Not floor, climb over edge
            //

            LogMessage("      Climbing over non-floor edge");

            // Find opposite triangle
            ElementIndex const oppositeTriangle = mesh.GetEdges().GetOppositeTriangle(intersectionEdgeElementIndex, npcParticleConstrainedState.CurrentTriangle);
            if (oppositeTriangle == NoneElementIndex || mesh.GetTriangles().IsDeleted(oppositeTriangle))
            {
                //
                // Become free
                //

                LogMessage("      No opposite triangle found, becoming free");

                //
                // Move to endpoint and exit, consuming whole quantum
                //

                npcParticle.ConstrainedState.reset();

                UpdateNpcParticle_Free(
                    npcParticle,
                    particleStartAbsolutePosition,
                    segmentTrajectoryEndAbsolutePosition,
                    particles,
                    labParameters);

                vec2f const totalTraveledVector = intersectionAbsolutePosition - particleStartAbsolutePosition; // We consider constrained only
                return totalTraveledVector.length();
            }
            else
            {
                //
                // Move to edge of opposite triangle 
                //

                int const oppositeTriangleEdgeOrdinal = mesh.GetTriangles().GetSubEdgeOrdinal(oppositeTriangle, intersectionEdgeElementIndex);

                LogMessage("      Moving to edge ", oppositeTriangleEdgeOrdinal, " of opposite triangle ", oppositeTriangle);

                // Move to triangle
                npcParticleConstrainedState.CurrentTriangle = oppositeTriangle;

                // Calculate new current barycentric coords (wrt new triangle)
                bcoords3f newBarycentricCoords; // In new triangle
                newBarycentricCoords[(oppositeTriangleEdgeOrdinal + 2) % 3] = 0.0f;
                newBarycentricCoords[oppositeTriangleEdgeOrdinal] = npcParticleConstrainedState.CurrentTriangleBarycentricCoords[(intersectionEdgeOrdinal + 1) % 3];
                newBarycentricCoords[(oppositeTriangleEdgeOrdinal + 1) % 3] = npcParticleConstrainedState.CurrentTriangleBarycentricCoords[intersectionEdgeOrdinal];

                LogMessage("      B-Coords: ", npcParticleConstrainedState.CurrentTriangleBarycentricCoords, " -> ", newBarycentricCoords);

                assert(newBarycentricCoords.is_on_edge_or_internal());

                npcParticleConstrainedState.CurrentTriangleBarycentricCoords = newBarycentricCoords;

                // Translate target coords to this triangle, for next iteration
                auto const oldSegmentTrajectoryEndBarycentricCoords = segmentTrajectoryEndBarycentricCoords; // For logging
                // Note: here we introduce a lot of error - the target bary coords are not anymore
                // guaranteed to lie exactly on the (continuation of the) edge
                segmentTrajectoryEndBarycentricCoords = mesh.GetTriangles().ToBarycentricCoordinates(
                    segmentTrajectoryEndAbsolutePosition,
                    oppositeTriangle,
                    mesh.GetVertices());

                LogMessage("      TrajEndB-Coords: ", oldSegmentTrajectoryEndBarycentricCoords, " -> ", segmentTrajectoryEndBarycentricCoords);

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
    NpcParticles & particles,
    Mesh const & mesh,
    LabParameters const & labParameters) const
{
    //
    // Note: we don't move here
    //

    auto & npcParticle = isPrimaryParticle ? npc.PrimaryParticleState : npc.DipoleState->SecondaryParticleState;
    assert(npcParticle.ConstrainedState.has_value());

    for (int iIter = 0; ; ++iIter)
    {
        LogDebug("    NavigateVertex: iter=", iIter);

        // The two vertices around the vertex we are on - seen in clockwise order
        int const nextVertexOrdinal = (vertexOrdinal + 1) % 3;
        int const prevVertexOrdinal = (vertexOrdinal + 2) % 3;

        // Pre-conditions: we are at this vertex
        assert(npcParticle.ConstrainedState->CurrentTriangleBarycentricCoords[vertexOrdinal] == 1.0f);
        assert(npcParticle.ConstrainedState->CurrentTriangleBarycentricCoords[nextVertexOrdinal] == 0.0f);
        assert(npcParticle.ConstrainedState->CurrentTriangleBarycentricCoords[prevVertexOrdinal] == 0.0f);

        LogDebug("      Triangle=", npcParticle.ConstrainedState->CurrentTriangle, " Vertex=", vertexOrdinal, " TrajectoryEndBarycentricCoords=", trajectoryEndBarycentricCoords);

        if (isInitialStateUnknown)
        {
            // Check whether we are directed towards the *interior* of this triangle, or elsewhere
            if (trajectoryEndBarycentricCoords[prevVertexOrdinal] >= 0.0f
                && trajectoryEndBarycentricCoords[nextVertexOrdinal] >= 0.0f)
            {
                //
                // We go inside this triangle - stop where we are, we'll then check trajectory in new situation
                //

                LogDebug("      Trajectory extends inside triangle - CompletedNavigation");

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

        LogDebug("      Trajectory crosses triangle: crossedEdgeOrdinal=", crossedEdgeOrdinal);

        //
        // Check whether this new edge is floor
        //

        ElementIndex const crossedEdgeElementIndex = mesh.GetTriangles().GetSubEdges(npcParticle.ConstrainedState->CurrentTriangle).EdgeIndices[crossedEdgeOrdinal];

        if (IsEdgeFloorToParticle(crossedEdgeElementIndex, npcParticle.ConstrainedState->CurrentTriangle, mesh)
            && (isPrimaryParticle || !DoesFloorSeparateFromPrimaryParticle(
                mParticles.GetPosition(npc.DipoleState->SecondaryParticleState.ParticleIndex),
                trajectoryStartAbsolutePosition, // Current (virtual, not yet real) position of this (secondary) particle
                crossedEdgeElementIndex,
                mesh)))
        {
            //
            // Encountered floor
            //

            LogDebug("      Crossed edge is floor - EncounteredFloor");

            return NavigateVertexOutcome::MakeEncounteredFloorOutcome(crossedEdgeOrdinal);
        }

        //
        // Not floor, climb over edge
        //

        // Find opposite triangle
        ElementIndex const oppositeTriangle = mesh.GetEdges().GetOppositeTriangle(crossedEdgeElementIndex, npcParticle.ConstrainedState->CurrentTriangle);
        if (oppositeTriangle == NoneElementIndex || mesh.GetTriangles().IsDeleted(oppositeTriangle))
        {
            //
            // Become free
            //

            LogDebug("      No opposite triangle found, becoming free - ConvertedToFree");

            npcParticle.ConstrainedState.reset();

            UpdateNpcParticle_Free(
                npcParticle,
                particleStartAbsolutePosition,
                trajectoryEndAbsolutePosition,
                particles,
                labParameters);

            return NavigateVertexOutcome::MakeConvertedToFreeOutcome();
        }

        //
        // Move to triangle
        //

        int const oppositeTriangleCrossedEdgeOrdinal = mesh.GetTriangles().GetSubEdgeOrdinal(oppositeTriangle, crossedEdgeElementIndex);

        LogDebug("      Moving to edge ", oppositeTriangleCrossedEdgeOrdinal, " of opposite triangle ", oppositeTriangle);

        npcParticle.ConstrainedState->CurrentTriangle = oppositeTriangle;

        // Calculate new current barycentric coords (wrt opposite triangle - note that we haven't moved)
        bcoords3f newBarycentricCoords; // In new triangle
        newBarycentricCoords[(oppositeTriangleCrossedEdgeOrdinal + 2) % 3] = 0.0f;
        newBarycentricCoords[oppositeTriangleCrossedEdgeOrdinal] = npcParticle.ConstrainedState->CurrentTriangleBarycentricCoords[(crossedEdgeOrdinal + 1) % 3];
        newBarycentricCoords[(oppositeTriangleCrossedEdgeOrdinal + 1) % 3] = npcParticle.ConstrainedState->CurrentTriangleBarycentricCoords[crossedEdgeOrdinal];

        LogDebug("      B-Coords: ", npcParticle.ConstrainedState->CurrentTriangleBarycentricCoords, " -> ", newBarycentricCoords);

        assert(newBarycentricCoords.is_on_edge_or_internal());

        npcParticle.ConstrainedState->CurrentTriangleBarycentricCoords = newBarycentricCoords;

        // New vertex: we know that coord of vertex opposite of crossed edge (i.e. vertex with ordinal crossed_edge+2) is 0.0
        if (newBarycentricCoords[oppositeTriangleCrossedEdgeOrdinal] == 0.0f)
        {
            // Between edge and edge+1
            vertexOrdinal = (oppositeTriangleCrossedEdgeOrdinal + 1) % 3;
        }
        else
        {
            // Between edge and edge-1
            assert(newBarycentricCoords[(oppositeTriangleCrossedEdgeOrdinal + 1) % 3] == 0.0f);
            vertexOrdinal = oppositeTriangleCrossedEdgeOrdinal;
        }

        //
        // Translate target bary coords
        //

        trajectoryEndBarycentricCoords = mesh.GetTriangles().ToBarycentricCoordinates(
            trajectoryEndAbsolutePosition,
            npcParticle.ConstrainedState->CurrentTriangle,
            mesh.GetVertices());

        LogMessage("      TrajEndB-Coords: ", trajectoryEndBarycentricCoords);

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
    LabParameters const & labParameters) const
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
        * labParameters.ElasticityAdjustment;

    // Calculate tangential response: Vt' = a*Vt (a = (1.0-friction), [0.0 - 1.0])
    vec2f const tangentialResponse =
        tangentialVelocity
        * std::max(0.0f, 1.0f - particles.GetPhysicalProperties(npcParticle.ParticleIndex).KineticFriction * labParameters.KineticFrictionAdjustment);

    // Given that we've been working in *apparent* space (we've calc'd the collision response to *trajectory* which is apparent displacement),
    // we need to transform velocity to absolute particle velocity
    vec2f const resultantAbsoluteVelocity = (normalResponse + tangentialResponse) - meshVelocity;

    LogMessage("        trajectory=", trajectory, " apparentParticleVelocity=", apparentParticleVelocity, " nr=", normalResponse, " tr=", tangentialResponse, " rr=", resultantAbsoluteVelocity);

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
        normalVelocity,
        bounceEdgeNormal,
        npc,
        isPrimaryParticle);
}

void Npcs::OnImpact(
    vec2f const & impactVector,
    vec2f const & bounceEdgeNormal, // Pointing outside of triangle
    StateType & npc,
    bool /*isPrimaryParticle*/) const
{
    LogMessage("    OnImpact(", impactVector, ", ", bounceEdgeNormal, ")");

    // Human state machine
    if (npc.HumanNpcState.has_value())
    {
        switch (npc.HumanNpcState->CurrentBehavior)
        {
            case StateType::HumanNpcStateType::BehaviorType::Constrained_Walking:
            {
                // Check alignment of impact with walking direction; if hit => flip
                if (bounceEdgeNormal.dot(vec2f(npc.HumanNpcState->CurrentFaceDirectionX, 0.0f)) > 0.0f
                    && npc.HumanNpcState->CurrentBehaviorState.Constrained_Walking.CurrentWalkMagnitude != 0.0f)
                {
                    // Flip now
                    FlipHumanWalk(*npc.HumanNpcState, StrongTypedTrue<_DoImmediate>);
                }

                break;
            }

            default:
            {
                break;
            }
        }
    }
}

void Npcs::UpdateNpcAnimation(
    StateType & npc,
    bool isPrimaryParticle,
    Mesh const & mesh,
    LabParameters const & labParameters)
{
    if (npc.Type == NpcType::Human && isPrimaryParticle)
    {
        assert(npc.DipoleState.has_value());
        assert(npc.HumanNpcState.has_value());

        ElementIndex const primaryParticleIndex = npc.PrimaryParticleState.ParticleIndex;
        ElementIndex const secondaryParticleIndex = npc.DipoleState->SecondaryParticleState.ParticleIndex;

        float const humanHeight = LabParameters::HumanNpcGeometry::BodyLength * labParameters.HumanNpcBodyLengthAdjustment;
        vec2f const headPosition = mParticles.GetPosition(secondaryParticleIndex);
        vec2f const feetPosition = mParticles.GetPosition(primaryParticleIndex);
        float const legLength = LabParameters::HumanNpcGeometry::LegLengthFraction * humanHeight;
        vec2f const crotchPosition = headPosition + (feetPosition - headPosition) * (LabParameters::HumanNpcGeometry::HeadLengthFraction + LabParameters::HumanNpcGeometry::TorsoLengthFraction);

        float targetLegRightAngle = 0.0f;
        float targetLegLeftAngle = 0.0f;
        float angleConvergenceRate = 0.0f;

        std::optional<vec2f> rightFootPosition;
        std::optional<vec2f> leftFootPosition;
        
        switch (npc.HumanNpcState->CurrentBehavior)
        {
            case StateType::HumanNpcStateType::BehaviorType::Constrained_KnockedOut:
            case StateType::HumanNpcStateType::BehaviorType::Constrained_Rising:
            case StateType::HumanNpcStateType::BehaviorType::Constrained_Equilibrium:
            {
                targetLegRightAngle = 0.0f;
                targetLegLeftAngle = 0.0f;
                angleConvergenceRate = 0.3f;

                break;
            }

            case StateType::HumanNpcStateType::BehaviorType::Constrained_Walking:
            {
                float const MaxLegAngle = std::atan((LabParameters::HumanNpcGeometry::StepLengthFraction / 2.0f) / LabParameters::HumanNpcGeometry::LegLengthFraction);

                float const stepLength = LabParameters::HumanNpcGeometry::StepLengthFraction * humanHeight;
                float const distanceInTwoSteps = std::fmod(npc.HumanNpcState->TotalDistanceTraveledSinceStateTransition + 3.0f * stepLength / 2.0f, stepLength * 2.0f);
                LogMessage("distanceInTwoSteps=", distanceInTwoSteps);

                targetLegRightAngle = std::abs(stepLength - distanceInTwoSteps) / stepLength * 2.0f * MaxLegAngle - MaxLegAngle;
                targetLegLeftAngle = -targetLegRightAngle;
                angleConvergenceRate = 1.0f;
                
                if (npc.PrimaryParticleState.ConstrainedState->CurrentVirtualEdgeElementIndex != NoneElementIndex)
                {
                    // Constrain feet onto current virtual edge

                    vec2f const e1 = mesh.GetEdges().GetEndpointAPosition(npc.PrimaryParticleState.ConstrainedState->CurrentVirtualEdgeElementIndex, mesh.GetVertices());
                    vec2f const e2 = mesh.GetEdges().GetEndpointBPosition(npc.PrimaryParticleState.ConstrainedState->CurrentVirtualEdgeElementIndex, mesh.GetVertices());

                    vec2f const r1 = crotchPosition;
                    vec2f const r2 = feetPosition + vec2f(1.0f, 0.0f) * std::tan(targetLegRightAngle) * legLength;

                    vec2f const l1 = crotchPosition;
                    vec2f const l2 = feetPosition + vec2f(1.0f, 0.0f) * std::tan(targetLegLeftAngle) * legLength;

                    // ((x1 - x0) * v1 - u1 * (y1 - y0)) / (u0 * v1 - u1 * v0)
                    // x0, y0 == e1
                    // x1, y1 == r1
                    // u0, v0 == (e2 - e1) // TODO: norm?
                    // u1, v1 == (r2 - r1) // TODO: norm?

                    vec2f const uv0 = (e2 - e1).normalise();
                    vec2f const uvr1 = (r2 - r1).normalise();
                    vec2f const uvl1 = (l2 - l1).normalise();

                    float const tr = ((r1.x - e1.x) * uvr1.y - (r1.y - e1.y) * uvr1.x) / (uv0.x * uvr1.y - uvr1.x * uv0.y);
                    float const tl = ((r1.x - e1.x) * uvl1.y - (r1.y - e1.y) * uvl1.x) / (uv0.x * uvl1.y - uvl1.x * uv0.y);

                    rightFootPosition = e1 + uv0 * tr;
                    leftFootPosition = e1 + uv0 * tl;
                }


                break;
            }

            case StateType::HumanNpcStateType::BehaviorType::Free_Aerial:
            {
                targetLegRightAngle = 0.0f;
                targetLegLeftAngle = 0.0f;
                angleConvergenceRate = 0.3f;

                break;
            }
        }

        npc.HumanNpcState->LegRightAngle += (targetLegRightAngle - npc.HumanNpcState->LegRightAngle) * angleConvergenceRate;
        npc.HumanNpcState->LegLeftAngle += (targetLegLeftAngle - npc.HumanNpcState->LegLeftAngle) * angleConvergenceRate;

        npc.HumanNpcState->TopPoint = headPosition;
        npc.HumanNpcState->NeckPoint = headPosition + (feetPosition - headPosition) * LabParameters::HumanNpcGeometry::HeadLengthFraction;
        npc.HumanNpcState->CrotchPoint = headPosition + (feetPosition - headPosition) * (LabParameters::HumanNpcGeometry::HeadLengthFraction + LabParameters::HumanNpcGeometry::TorsoLengthFraction);

        if (rightFootPosition.has_value())
            npc.HumanNpcState->LegRightPoint = *rightFootPosition;
        else
            npc.HumanNpcState->LegRightPoint = feetPosition + vec2f(1.0f, 0.0f) * std::tan(npc.HumanNpcState->LegRightAngle) * legLength;

        if (leftFootPosition.has_value())
            npc.HumanNpcState->LegLeftPoint = *leftFootPosition;
        else
            npc.HumanNpcState->LegLeftPoint = feetPosition + vec2f(1.0f, 0.0f) * std::tan(npc.HumanNpcState->LegLeftAngle) * legLength;
    }
}

Npcs::StateType Npcs::MaterializeNpcState(
    ElementIndex npcIndex,
    float currentSimulationTime,
    Mesh const & mesh) const
{
    auto const & state = mStateBuffer[npcIndex];

    // Primary particle

    StateType::NpcParticleStateType primaryParticleState = StateType::NpcParticleStateType(
        state.PrimaryParticleState.ParticleIndex,
        CalculateParticleConstrainedState(
            mParticles.GetPosition(state.PrimaryParticleState.ParticleIndex),
            mesh));

    // Secondary particle

    std::optional<StateType::DipoleStateType> dipoleState;

    if (state.DipoleState.has_value())
    {
        StateType::NpcParticleStateType secondaryParticleState = StateType::NpcParticleStateType(
            state.DipoleState->SecondaryParticleState.ParticleIndex,
            CalculateParticleConstrainedState(
                mParticles.GetPosition(state.DipoleState->SecondaryParticleState.ParticleIndex),
                mesh));

        dipoleState.emplace(
            std::move(secondaryParticleState),
            state.DipoleState->DipoleProperties);
    }

    // Human NPC state

    std::optional<StateType::HumanNpcStateType> humanNpcState;

    if (state.Type == NpcType::Human)
    {
        assert(dipoleState.has_value());

        humanNpcState = InitializeHuman(
            primaryParticleState,
            dipoleState->SecondaryParticleState,
            currentSimulationTime);
    }

    // Regime

    auto const regime = primaryParticleState.ConstrainedState.has_value()
        ? StateType::RegimeType::Constrained
        : (dipoleState.has_value() && dipoleState->SecondaryParticleState.ConstrainedState.has_value()) ? StateType::RegimeType::Constrained : StateType::RegimeType::Free;

    return StateType(
        state.Type,
        regime,
        std::move(primaryParticleState),
        std::move(dipoleState),
        std::move(humanNpcState));
}

std::optional<Npcs::StateType::NpcParticleStateType::ConstrainedStateType> Npcs::CalculateParticleConstrainedState(
    vec2f const & position,
    Mesh const & mesh) const
{
    std::optional<StateType::NpcParticleStateType::ConstrainedStateType> constrainedState;

    ElementIndex const triangleIndex = mesh.GetTriangles().FindContaining(position, mesh.GetVertices());
    if (triangleIndex != NoneElementIndex && !mesh.GetTriangles().IsDeleted(triangleIndex))
    {
        bcoords3f const barycentricCoords = mesh.GetTriangles().ToBarycentricCoordinatesFromWithinTriangle(
            position,
            triangleIndex,
            mesh.GetVertices());

        assert(barycentricCoords.is_on_edge_or_internal());

        return Npcs::StateType::NpcParticleStateType::ConstrainedStateType(
            triangleIndex,
            barycentricCoords);
    }

    return std::nullopt;
}
