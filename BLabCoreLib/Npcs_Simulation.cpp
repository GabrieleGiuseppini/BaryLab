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
    Mesh const & mesh,
    LabParameters const & labParameters)
{
    LogMessage("----------------------------------");
    LogMessage("----------------------------------");
    LogMessage("----------------------------------");

    //
    // 1. Reset buffers
    //

    mParticles.ResetEquilibriumTorque();
    // Note: no need to reset PreliminaryForces as we'll recalculate all of them

    //
    // 2. Check if a free secondary particle should become constrained
    // 3. Update behavioral state machines
    // 4. Calculate preliminary forces 
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

        // Behavior

        if (npcState.Type == NpcType::Human)
        {
            assert(npcState.DipoleState.has_value());
            assert(npcState.HumanNpcState.has_value());

            UpdateHuman(
                *npcState.HumanNpcState,
                npcState.PrimaryParticleState,
                npcState.DipoleState->SecondaryParticleState,
                mesh,
                labParameters);
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
    // 5. Update state
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

        ElementIndex currentTriangleElementIndex = npcParticle.ConstrainedState->CurrentTriangle;

        // Initialize absolute position of particle (wrt current triangle) as if it moved with the mesh, staying in its position wrt its triangle;
        // it's the new theoretical position after mesh displacement
        vec2f trajectoryStartAbsolutePosition = mesh.GetTriangles().FromBarycentricCoordinates(
            npcParticle.ConstrainedState->CurrentTriangleBarycentricCoords,
            currentTriangleElementIndex,
            mesh.GetVertices());

        // Calculate mesh displacement for the whole loop via pure displacement of triangle containing particle
        vec2f const meshVelocity = (particleStartAbsolutePosition - trajectoryStartAbsolutePosition) / LabParameters::SimulationTimeStepDuration;

        // Machinery to detect 3-iteration paths that doesn't move particle (positional well)

        struct PastBarycentricPosition
        {
            ElementIndex TriangleElementIndex;
            vec3f BarycentricCoords;

            PastBarycentricPosition(
                ElementIndex triangleElementIndex,
                vec3f barycentricCoords)
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
            LogMessage("    Triangle=", currentTriangleElementIndex, " B-Coords=", npcParticle.ConstrainedState->CurrentTriangleBarycentricCoords, " RemainingDt=", remainingDt);

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

            int nonZeroVertex = -1;
            int zeroCoordsCount = 0;

            if (npcParticle.ConstrainedState->CurrentTriangleBarycentricCoords[0] == 0.0f)
            {
                ++zeroCoordsCount;
            }
            else
            {
                nonZeroVertex = 0;
            }

            if (npcParticle.ConstrainedState->CurrentTriangleBarycentricCoords[1] == 0.0f)
            {
                ++zeroCoordsCount;
            }
            else
            {
                nonZeroVertex = 1;
            }

            if (npcParticle.ConstrainedState->CurrentTriangleBarycentricCoords[2] == 0.0f)
            {
                ++zeroCoordsCount;
            }
            else
            {
                nonZeroVertex = 2;
            }

            while (zeroCoordsCount == 2)
            {
                //
                // We are at a corner
                //

                assert(nonZeroVertex >= 0);
                int const nextNonZeroVertexOrdinal = (nonZeroVertex + 1) % 3;
                int const prevNonZeroVertexOrdinal = (nonZeroVertex + 2) % 3;

                vec3f const trajectoryEndBarycentricCoords = mesh.GetTriangles().ToBarycentricCoordinates(
                    trajectoryEndAbsolutePosition,
                    currentTriangleElementIndex,
                    mesh.GetVertices());

                LogMessage("    At a corner: vertex=", nonZeroVertex, " TrajectoryEndBarycentricCoords=", trajectoryEndBarycentricCoords);

                if (trajectoryEndBarycentricCoords[prevNonZeroVertexOrdinal] >= 0.0f
                    && trajectoryEndBarycentricCoords[nextNonZeroVertexOrdinal] >= 0.0f)
                {
                    //
                    // We go inside this triangle - stop where we are, we'll then check trajectory in new situation
                    //

                    LogMessage("        Trajectory extends inside new triangle - exiting and continuing");

                    break;
                }

                //
                // Find next edge that we intersect at this vertex
                //

                int crossedEdgeOrdinal;
                if (trajectoryEndBarycentricCoords[prevNonZeroVertexOrdinal] <= trajectoryEndBarycentricCoords[nextNonZeroVertexOrdinal])
                {
                    // Clockwise - next edge
                    crossedEdgeOrdinal = nonZeroVertex;
                }
                else
                {
                    // Anti-clockwise - prev edge
                    crossedEdgeOrdinal = (nonZeroVertex + 2) % 3;
                }

                LogMessage("        Trajectory crosses triangle: crossedEdgeOrdinal=", crossedEdgeOrdinal);

                //
                // Check whether edge is floor
                //

                ElementIndex const crossedEdgeElementIndex = mesh.GetTriangles().GetSubEdges(currentTriangleElementIndex).EdgeIndices[crossedEdgeOrdinal];

                if (IsEdgeFloorToParticle(crossedEdgeElementIndex, currentTriangleElementIndex, mesh)
                    && (isPrimaryParticle || !DoesFloorSeparateFromPrimaryParticle(
                        mParticles.GetPosition(npc.DipoleState->SecondaryParticleState.ParticleIndex),
                        trajectoryStartAbsolutePosition, // Current (virtual, not yet real) position of this (secondary) particle
                        crossedEdgeElementIndex,
                        mesh)))
                {
                    LogMessage("        Crossed edge is floor - exiting and continuing");
                    break;
                }

                //
                // Not floor, climb over edge
                //

                // Find opposite triangle
                ElementIndex const oppositeTriangle = mesh.GetEdges().GetOppositeTriangle(crossedEdgeElementIndex, currentTriangleElementIndex);
                if (oppositeTriangle == NoneElementIndex)
                {
                    //
                    // Become free
                    //

                    LogMessage("        No opposite triangle found, becoming free");

                    //
                    // Move to endpoint and exit, consuming whole quantum
                    //

                    npcParticle.ConstrainedState.reset();

                    UpdateNpcParticle_Free(
                        npcParticle,
                        particleStartAbsolutePosition,
                        trajectoryEndAbsolutePosition,
                        mParticles,
                        labParameters);

                    remainingDt = 0.0f;
                    break;
                }

                //
                // Move to triangle
                //

                int const oppositeTriangleCrossedEdgeOrdinal = mesh.GetTriangles().GetSubEdgeOrdinal(oppositeTriangle, crossedEdgeElementIndex);

                LogMessage("        Moving to edge ", oppositeTriangleCrossedEdgeOrdinal, " of opposite triangle ", oppositeTriangle);

                npcParticle.ConstrainedState->CurrentTriangle = oppositeTriangle;
                currentTriangleElementIndex = npcParticle.ConstrainedState->CurrentTriangle; // TODO: get rid of this

                // New barycentric coords (we haven't moved but just changed triangle)
                // Calculate new current barycentric coords (wrt opposite triangle)
                vec3f newBarycentricCoords; // In new triangle
                newBarycentricCoords[(oppositeTriangleCrossedEdgeOrdinal + 2) % 3] = 0.0f;
                newBarycentricCoords[oppositeTriangleCrossedEdgeOrdinal] = npcParticle.ConstrainedState->CurrentTriangleBarycentricCoords[(crossedEdgeOrdinal + 1) % 3];
                newBarycentricCoords[(oppositeTriangleCrossedEdgeOrdinal + 1) % 3] = npcParticle.ConstrainedState->CurrentTriangleBarycentricCoords[crossedEdgeOrdinal];

                LogMessage("          B-Coords: ", npcParticle.ConstrainedState->CurrentTriangleBarycentricCoords, " -> ", newBarycentricCoords);

                {
                    assert(newBarycentricCoords[0] >= 0.0f && newBarycentricCoords[0] <= 1.0f);
                    assert(newBarycentricCoords[1] >= 0.0f && newBarycentricCoords[1] <= 1.0f);
                    assert(newBarycentricCoords[2] >= 0.0f && newBarycentricCoords[2] <= 1.0f);
                }

                npcParticle.ConstrainedState->CurrentTriangleBarycentricCoords = newBarycentricCoords;

                // New vertex: we know that coord of vertex opposite of crossed edge (i.e. vertex with ordinal crossed_edge+2) is 0.0
                if (newBarycentricCoords[oppositeTriangleCrossedEdgeOrdinal] == 0.0f)
                {
                    // Between edge and edge+1
                    nonZeroVertex = (oppositeTriangleCrossedEdgeOrdinal + 1) % 3;
                }
                else
                {
                    // Between edge and edge-1
                    assert(newBarycentricCoords[(oppositeTriangleCrossedEdgeOrdinal + 1) % 3] == 0.0f);
                    nonZeroVertex = oppositeTriangleCrossedEdgeOrdinal;
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
                    ElementIndex const currentEdgeElementIndex = mesh.GetTriangles().GetSubEdges(currentTriangleElementIndex).EdgeIndices[edgeOrdinal];

                    assert(isPrimaryParticle || npc.DipoleState.has_value());

                    LogMessage("      edge ", edgeOrdinal, ": isFloor=", IsEdgeFloorToParticle(currentEdgeElementIndex, currentTriangleElementIndex, mesh));

                    // Check if this is really a floor to this particle
                    if (IsEdgeFloorToParticle(currentEdgeElementIndex, currentTriangleElementIndex, mesh)
                        && (isPrimaryParticle || !DoesFloorSeparateFromPrimaryParticle(
                            mParticles.GetPosition(npc.DipoleState->SecondaryParticleState.ParticleIndex),
                            trajectoryStartAbsolutePosition, // Current (virtual, not yet real) position of this (secondary) particle
                            currentEdgeElementIndex,
                            mesh)))
                    {
                        // On floor edge - so potentially in a non-inertial frame

                        // Check now whether we're moving *against* the floor

                        vec2f const edgeVector = mesh.GetTriangles().GetSubEdgeVector(
                            currentTriangleElementIndex,
                            edgeOrdinal,
                            mesh.GetVertices());

                        vec2f const edgeDir = edgeVector.normalise();
                        vec2f const edgeNormal = edgeDir.to_perpendicular(); // Points outside of triangle (i.e. towards floor)

                        if (trajectory.dot(edgeNormal) >= 0.0f) // If 0, no normal force - hence no friction; however we want to take this codepath anyways for consistency
                        {
                            //
                            // Case 1: Non-inertial: on edge and moving against (or along) it, pushed by it
                            //

                            LogMessage("    ConstrainedNonInertial: triangle=", currentTriangleElementIndex, " edgeOrdinal=", edgeOrdinal, " bCoords=", npcParticle.ConstrainedState->CurrentTriangleBarycentricCoords, " trajectory=", trajectory);
                            LogMessage("    StartPosition=", mParticles.GetPosition(npcParticle.ParticleIndex), " StartVelocity=", mParticles.GetVelocity(npcParticle.ParticleIndex), " MeshVelocity=", meshVelocity, " StartMRVelocity=", npcParticle.ConstrainedState->MeshRelativeVelocity);

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

                            if (npc.HumanNpcState.has_value() && isPrimaryParticle)
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
                                    * npc.HumanNpcState->CurrentWalkMagnitude;

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

                                if (npc.HumanNpcState->CurrentWalkMagnitude != 0.0f)
                                {
                                    LogMessage("        idealWalkMagnitude=", idealWalkMagnitude, " gravityResistance=", gravityResistance, " => edgeWalkedPlanned=", edgeWalkedPlanned, " (@", npc.HumanNpcState->CurrentWalkMagnitude, ")");
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

                            vec3f flattenedTrajectoryEndBarycentricCoords = mesh.GetTriangles().ToBarycentricCoordinates(
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

                            std::optional<float> const edgeTraveledActual = UpdateNpcParticle_ConstrainedNonInertial(
                                npc,
                                isPrimaryParticle,
                                edgeOrdinal,
                                edgeDir,
                                particleStartAbsolutePosition,
                                trajectoryStartAbsolutePosition,
                                flattenedTrajectoryEndBarycentricCoords,
                                flattenedTrajectory,
                                edgeTraveledPlanned,
                                meshVelocity,
                                remainingDt,
                                mParticles,
                                mesh,
                                labParameters);

                            LogMessage("    Actual edge traveled in non-inertial step: ", edgeTraveledActual.has_value() ? std::to_string(*edgeTraveledActual) : "N/A");

                            float edgeWalkedActual;

                            if (!edgeTraveledActual.has_value())
                            {
                                // Completed
                                remainingDt = 0.0f;

                                // Walked is the planned one
                                edgeWalkedActual = edgeWalkedPlanned;
                            }
                            else
                            {
                                if (*edgeTraveledActual == 0.0f)
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

                                    // Calculated walked actual
                                    //
                                    // If actual resultant is zero, it's either because planned resultant was zero, 
                                    // or because we've cut short a trajectory.
                                    // If planned resultant was zero then we've reached target - this step should have completed;
                                    // however, for safety, we may assume this happens (i.e. that it may return zero when planned was zero).
                                    // In this case do nothing, as we'll complete later on.
                                    // If on the other hand we've cut short, since planned was non-zero and actual is zero, then walked now is also zero
                                    edgeWalkedActual = 0.0f;
                                }
                                else
                                {
                                    // We have moved

                                    assert(edgeTraveledActual.has_value());

                                    // Calculate consumed dt
                                    assert(*edgeTraveledActual * edgeTraveledPlanned >= 0.0f); // Should have same sign
                                    float const dtFractionConsumed = edgeTraveledPlanned != 0.0f
                                        ? std::min(*edgeTraveledActual / edgeTraveledPlanned, 1.0f) // Signs should agree anyway
                                        : 1.0f; // If we were planning no travel, any movement is a whole consumption
                                    LogMessage("        dtFractionConsumed=", dtFractionConsumed);
                                    remainingDt *= (1.0f - dtFractionConsumed);

                                    // Calculate actual walked via planned walked's proportion with total
                                    edgeWalkedActual = edgeTraveledPlanned != 0.0f
                                        ? *edgeTraveledActual * (edgeWalkedPlanned / edgeTraveledPlanned)
                                        : 0.0f; // Unlikely, but read above for rationale behind 0.0

                                    // Reset well detection machinery
                                    pastPastBarycentricPosition.reset();
                                    pastBarycentricPosition.reset();
                                }
                            }

                            // Update total vector walked along edge
                            totalEdgeWalkedActual += edgeDir * edgeWalkedActual;
                            LogMessage("        totalEdgeWalkedActual=", totalEdgeWalkedActual);

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

                LogMessage("    ConstrainedInertial: triangle=", currentTriangleElementIndex, " bCoords=", npcParticle.ConstrainedState->CurrentTriangleBarycentricCoords, " physicsDeltaPos=", physicsDeltaPos);
                LogMessage("    StartPosition=", mParticles.GetPosition(npcParticle.ParticleIndex), " StartVelocity=", mParticles.GetVelocity(npcParticle.ParticleIndex), " MeshVelocity=", meshVelocity, " StartMRVelocity=", npcParticle.ConstrainedState->MeshRelativeVelocity);

                //
                // Calculate target barycentric coords
                //

                vec3f const trajectoryEndBarycentricCoords = mesh.GetTriangles().ToBarycentricCoordinates(
                    trajectoryEndAbsolutePosition,
                    npcParticle.ConstrainedState->CurrentTriangle,
                    mesh.GetVertices());

                //
                // Move towards target bary coords
                //

                UpdateNpcParticle_ConstrainedInertial(
                    npc,
                    isPrimaryParticle,
                    particleStartAbsolutePosition,
                    npcParticle.ConstrainedState->CurrentTriangleBarycentricCoords, // trajectoryStartBarycentricCoords
                    trajectoryEndBarycentricCoords,
                    meshVelocity,
                    remainingDt,
                    mParticles,
                    mesh,
                    labParameters);

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

            currentTriangleElementIndex = npcParticle.ConstrainedState->CurrentTriangle;

            // Current (virtual, not yet real) position of this particle
            trajectoryStartAbsolutePosition = mesh.GetTriangles().FromBarycentricCoordinates(
                npcParticle.ConstrainedState->CurrentTriangleBarycentricCoords,
                currentTriangleElementIndex,
                mesh.GetVertices());
        }

        //
        // Store total edge walked, if this is the primary of a human
        //

        if (npc.HumanNpcState.has_value() && isPrimaryParticle)
        {
            // TODO: add to new HumanNpcState member
            // TODO: decied if we need total edge *walked* or total edge *traveled* for animation
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

        float const uwCoefficient = Clamp(labParameters.SeaLevel - testParticlePosition.y, 0.0f, BuoyancyInterfaceWidth) / BuoyancyInterfaceWidth;
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

std::optional<float> Npcs::UpdateNpcParticle_ConstrainedNonInertial(
    StateType & npc,
    bool isPrimaryParticle,
    int edgeOrdinal,
    vec2f const & edgeDir,
    vec2f const & particleStartAbsolutePosition,
    vec2f const & trajectoryStartAbsolutePosition,
    vec3f trajectoryEndBarycentricCoords,
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

    if (trajectoryEndBarycentricCoords.x >= 0.0f && trajectoryEndBarycentricCoords.x <= 1.0f
        && trajectoryEndBarycentricCoords.y >= 0.0f && trajectoryEndBarycentricCoords.y <= 1.0f
        && trajectoryEndBarycentricCoords.z >= 0.0f && trajectoryEndBarycentricCoords.z <= 1.0f)
    {
        LogMessage("      Target is on/in triangle, moving to target");

        //
        // Update particle and exit - consuming whole time quantum
        //            

        npcParticleConstrainedState.CurrentTriangleBarycentricCoords = trajectoryEndBarycentricCoords;

        vec2f const particleEndAbsolutePosition = mesh.GetTriangles().FromBarycentricCoordinates(
            trajectoryEndBarycentricCoords,
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
        return std::nullopt;
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

    // Figure out if we're moving clockwise or counter-clockwise (fron the point of view of inside
    // the triangle); this will hold for all the triangles traveled while going through a cuspid
    bool isClockwise;
    if (trajectoryEndBarycentricCoords[edgeOrdinal] < 0.0f)
    {
        // It's the vertex at the end of our edge
        assert(trajectoryEndBarycentricCoords[(edgeOrdinal + 1) % 3] > 0.0f); // Because target is outside of triangle
        isClockwise = true;
    }
    else
    {
        // It's the vertex before our edge
        assert(trajectoryEndBarycentricCoords[(edgeOrdinal + 1) % 3] < 0.0f);
        isClockwise = false;
    }

    int intersectionEdgeOrdinal = (isClockwise)
        ? (edgeOrdinal + 1) % 3 // Next edge
        : (edgeOrdinal + 2) % 3; // Previous edge

    // Vertex we end at
    int intersectionVertexOrdinal = (isClockwise)
        ? intersectionEdgeOrdinal // Initial vertex
        : (intersectionEdgeOrdinal + 1) % 3; // Final vertex

    // In a loop:
    //  - Move to edge intersection
    //  - Check floor & bumps
    //  - Climb triangle
    //  - Find edge intersection

    float totalEdgeTraveled = 0.0f; // Signed

    for (vec2f lastTrajectoryStartAbsolutePosition = trajectoryStartAbsolutePosition; ;)
    {
        //
        // Move to intersection (which is a vertex by now)
        //
        // Note that except for the very first iteration, any other iteration will travel
        // zero distance at this moment
        //

        vec3f intersectionBarycentricCoords = vec3f::zero();
        intersectionBarycentricCoords[intersectionVertexOrdinal] = 1.0f;

        LogMessage("      Moving to intersection with edge ordinal ", intersectionEdgeOrdinal, ": ", intersectionBarycentricCoords);

        npcParticleConstrainedState.CurrentTriangleBarycentricCoords = intersectionBarycentricCoords;

        // Update (signed) edge traveled, assuming we're always moving along initial edge
        // (other edges we're only passing through cuspids)
        vec2f const intersectionAbsolutePosition = mesh.GetVertices().GetPosition(mesh.GetTriangles().GetVertexIndices(npcParticleConstrainedState.CurrentTriangle)[intersectionVertexOrdinal]);
        assert(intersectionAbsolutePosition == mesh.GetTriangles().FromBarycentricCoordinates(
            intersectionBarycentricCoords,
            npcParticleConstrainedState.CurrentTriangle,
            mesh.GetVertices()));

        float const edgeTraveled = (intersectionAbsolutePosition - lastTrajectoryStartAbsolutePosition).dot(edgeDir);
        totalEdgeTraveled += edgeTraveled;

        LogMessage("        totalEdgeTraveled=", totalEdgeTraveled);

        lastTrajectoryStartAbsolutePosition = intersectionAbsolutePosition;

        //
        // Check if impacted with floor
        //

        ElementIndex const intersectionEdgeElementIndex = mesh.GetTriangles().GetSubEdges(npcParticleConstrainedState.CurrentTriangle).EdgeIndices[intersectionEdgeOrdinal];

        assert(isPrimaryParticle || npc.DipoleState.has_value());

        if (IsEdgeFloorToParticle(intersectionEdgeElementIndex, npcParticleConstrainedState.CurrentTriangle, mesh)
            && (isPrimaryParticle || !DoesFloorSeparateFromPrimaryParticle(
                mParticles.GetPosition(npc.DipoleState->SecondaryParticleState.ParticleIndex),
                intersectionAbsolutePosition,
                intersectionEdgeElementIndex,
                mesh)))
        {
            //
            // Impact
            //

            LogMessage("      Impact");

            //
            // We might have hit a tiny bump (e.g. because of triangles slightly bent); in this case we don't want to bounce
            //

            float const flattenedTrajectoryLength = flattenedTrajectory.length();

            vec2f const intersectionEdgeDir =
                mesh.GetTriangles().GetSubEdgeVector(npcParticleConstrainedState.CurrentTriangle, intersectionEdgeOrdinal, mesh.GetVertices())
                .normalise();
            vec2f const intersectionEdgeNormal = intersectionEdgeDir.to_perpendicular();

            // Check angle between desired (original) trajectory and edge
            float const trajProjOntoEdgeNormal = flattenedTrajectory.normalise_approx(flattenedTrajectoryLength).dot(intersectionEdgeNormal);
            if (trajProjOntoEdgeNormal <= 0.71f) // PI/4+
            {
                //
                // Impact continuation (no bounce)
                //
                // Stop here and then check trajectory in new situation
                //

                LogMessage("      Impact continuation (trajProjOntoEdgeNormal=", trajProjOntoEdgeNormal, ")");

                return totalEdgeTraveled;
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
                    intersectionEdgeNormal,
                    meshVelocity,
                    dt,
                    particles,
                    labParameters);

                // Terminate
                return std::nullopt;
            }
        }

        //
        // Not floor, climb over edge
        //

        LogMessage("      Climbing over non-floor edge");

        // Find opposite triangle
        ElementIndex const oppositeTriangle = mesh.GetEdges().GetOppositeTriangle(intersectionEdgeElementIndex, npcParticleConstrainedState.CurrentTriangle);
        if (oppositeTriangle == NoneElementIndex)
        {
            //
            // Become free
            //

            LogMessage("      No opposite triangle found, becoming free");

            //
            // Move to endpoint and exit, consuming whole quantum
            //

            vec2f const endPosition = mesh.GetTriangles().FromBarycentricCoordinates(
                trajectoryEndBarycentricCoords,
                npcParticleConstrainedState.CurrentTriangle,
                mesh.GetVertices());

            npcParticle.ConstrainedState.reset();

            UpdateNpcParticle_Free(
                npcParticle,
                particleStartAbsolutePosition,
                endPosition,
                particles,
                labParameters);

            return std::nullopt;
        }

        //
        // Move to edge of opposite triangle 
        //

        int const oppositeTriangleIntersectionEdgeOrdinal = mesh.GetTriangles().GetSubEdgeOrdinal(oppositeTriangle, intersectionEdgeElementIndex);

        LogMessage("      Moving to edge ", oppositeTriangleIntersectionEdgeOrdinal, " of opposite triangle ", oppositeTriangle);

        // Save current absolute trajectory end
        vec2f const oldTrajectoryEndAbsolutePosition = mesh.GetTriangles().FromBarycentricCoordinates(
            trajectoryEndBarycentricCoords,
            npcParticleConstrainedState.CurrentTriangle,
            mesh.GetVertices());

        npcParticleConstrainedState.CurrentTriangle = oppositeTriangle;

        // Calculate new current barycentric coords (wrt opposite triangle)
        vec3f newBarycentricCoords; // In new triangle
        newBarycentricCoords[(oppositeTriangleIntersectionEdgeOrdinal + 2) % 3] = 0.0f;
        newBarycentricCoords[oppositeTriangleIntersectionEdgeOrdinal] = npcParticleConstrainedState.CurrentTriangleBarycentricCoords[(intersectionEdgeOrdinal + 1) % 3];
        newBarycentricCoords[(oppositeTriangleIntersectionEdgeOrdinal + 1) % 3] = npcParticleConstrainedState.CurrentTriangleBarycentricCoords[intersectionEdgeOrdinal];

        LogMessage("        B-Coords: ", npcParticleConstrainedState.CurrentTriangleBarycentricCoords, " -> ", newBarycentricCoords);

        {
            assert(newBarycentricCoords[0] >= 0.0f && newBarycentricCoords[0] <= 1.0f);
            assert(newBarycentricCoords[1] >= 0.0f && newBarycentricCoords[1] <= 1.0f);
            assert(newBarycentricCoords[2] >= 0.0f && newBarycentricCoords[2] <= 1.0f);
        }

        npcParticleConstrainedState.CurrentTriangleBarycentricCoords = newBarycentricCoords;

        //
        // Translate trajectory end coords to this triangle
        //

        auto const oldTrajectoryEndBarycentricCoords = trajectoryEndBarycentricCoords;

        // Note: here we introduce a lot of error - the target bary coords are not anymore
        // guaranteed to lie exactly on the (continuation of the) edge
        trajectoryEndBarycentricCoords = mesh.GetTriangles().ToBarycentricCoordinates(
            oldTrajectoryEndAbsolutePosition,
            oppositeTriangle,
            mesh.GetVertices());

        LogMessage("        TrajEndB-Coords: ", oldTrajectoryEndBarycentricCoords, " -> ", trajectoryEndBarycentricCoords);

        //
        // Now move on to this triangle
        //

        //
        // See if trajectory requires another edge intersection & crossing
        //

        // Vertex we end at
        intersectionVertexOrdinal = (isClockwise)
            ? (oppositeTriangleIntersectionEdgeOrdinal + 1) % 3 // Final vertex of intersection edge
            : oppositeTriangleIntersectionEdgeOrdinal; // First vertex of intersection edge

        LogMessage("        Vertex=", intersectionVertexOrdinal);

        if (trajectoryEndBarycentricCoords[(intersectionVertexOrdinal + 1) % 3] >= 0.0f
            && trajectoryEndBarycentricCoords[(intersectionVertexOrdinal + 2) % 3] >= 0.0f)
        {
            //
            // It's inside this triangle - stop where we are, we'll then check trajectory in new situation
            //

            LogMessage("      Trajectory extends inside new triangle - exiting and continuing");

            return totalEdgeTraveled;
        }

        //
        // Find next edge that we intersect at this cuspid
        //

        LogMessage("      Trajectory crosses new triangle - finding next edge intersected at cuspid");

        intersectionEdgeOrdinal = (isClockwise)
            ? intersectionVertexOrdinal // Edge following vertex
            : (intersectionVertexOrdinal + 2) % 3; // Edge preceding vertex
    }
}

void Npcs::UpdateNpcParticle_ConstrainedInertial(
    StateType & npc,
    bool isPrimaryParticle,
    vec2f const & particleStartAbsolutePosition, // Since beginning of whole time quantum, not just this step
    vec3f const trajectoryStartBarycentricCoords,
    vec3f trajectoryEndBarycentricCoords, // Mutable
    vec2f const meshVelocity,
    float dt,
    NpcParticles & particles,
    Mesh const & mesh,
    LabParameters const & labParameters) const
{
    auto & npcParticle = isPrimaryParticle ? npc.PrimaryParticleState : npc.DipoleState->SecondaryParticleState;
    assert(npcParticle.ConstrainedState.has_value());
    auto & npcParticleConstrainedState = *npcParticle.ConstrainedState;

    vec2f const trajectoryStartAbsolutePosition = mesh.GetTriangles().FromBarycentricCoordinates(
        trajectoryStartBarycentricCoords,
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
        {
            assert(npcParticleConstrainedState.CurrentTriangleBarycentricCoords[0] >= 0.0f && npcParticleConstrainedState.CurrentTriangleBarycentricCoords[0] <= 1.0f);
            assert(npcParticleConstrainedState.CurrentTriangleBarycentricCoords[1] >= 0.0f && npcParticleConstrainedState.CurrentTriangleBarycentricCoords[1] <= 1.0f);
            assert(npcParticleConstrainedState.CurrentTriangleBarycentricCoords[2] >= 0.0f && npcParticleConstrainedState.CurrentTriangleBarycentricCoords[2] <= 1.0f);
        }

        LogMessage("    SegmentTrace ", iIter);
        LogMessage("      triangle=", npcParticleConstrainedState.CurrentTriangle, " bCoords=", npcParticleConstrainedState.CurrentTriangleBarycentricCoords,
            " trajStartBCoords=", trajectoryStartBarycentricCoords, " trajEndBCoords=", trajectoryEndBarycentricCoords);

        //
        // If target is on/in triangle, we move to target
        //

        if (trajectoryEndBarycentricCoords.x >= 0.0f && trajectoryEndBarycentricCoords.x <= 1.0f
            && trajectoryEndBarycentricCoords.y >= 0.0f && trajectoryEndBarycentricCoords.y <= 1.0f
            && trajectoryEndBarycentricCoords.z >= 0.0f && trajectoryEndBarycentricCoords.z <= 1.0f)
        {
            LogMessage("      Target is on/in triangle, moving to target");

            //
            // Update particle and exit - consuming whole time quantum
            //            

            npcParticleConstrainedState.CurrentTriangleBarycentricCoords = trajectoryEndBarycentricCoords;

            vec2f const particleEndAbsolutePosition = mesh.GetTriangles().FromBarycentricCoordinates(
                trajectoryEndBarycentricCoords,
                npcParticleConstrainedState.CurrentTriangle,
                mesh.GetVertices());

            particles.SetPosition(npcParticle.ParticleIndex, particleEndAbsolutePosition);

            // Use whole time quantum for velocity, as particleStartAbsolutePosition is fixed at t0
            vec2f const absoluteVelocity = (particleEndAbsolutePosition - particleStartAbsolutePosition) / LabParameters::SimulationTimeStepDuration;
            particles.SetVelocity(npcParticle.ParticleIndex, absoluteVelocity);
            npcParticleConstrainedState.MeshRelativeVelocity = absoluteVelocity + meshVelocity;

            LogMessage("        traveledActual=", (particleEndAbsolutePosition - particleStartAbsolutePosition), " absoluteVelocity=", particles.GetVelocity(npcParticle.ParticleIndex));

            // We have consumed the whole time quantum
            return;
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
            vec3f IntersectionPoint;

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
            if (trajectoryEndBarycentricCoords[vi] < 0.0f)
            {
                float const den = npcParticleConstrainedState.CurrentTriangleBarycentricCoords[vi] - trajectoryEndBarycentricCoords[vi];
                float const t = npcParticleConstrainedState.CurrentTriangleBarycentricCoords[vi] / den;

#ifdef _DEBUG
                diags[vi].emplace(den, t);
                diags[vi]->IntersectionPoint =
                    npcParticleConstrainedState.CurrentTriangleBarycentricCoords
                    + (trajectoryEndBarycentricCoords - npcParticleConstrainedState.CurrentTriangleBarycentricCoords) * t;
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

        vec3f intersectionBarycentricCoords;
        intersectionBarycentricCoords[intersectionVertexOrdinal] = 0.0f;
        float const lNext = Clamp( // Barycentric coord of next vertex at intersection; enforcing it's within triangle
            npcParticleConstrainedState.CurrentTriangleBarycentricCoords[(intersectionVertexOrdinal + 1) % 3] * (1.0f - minIntersectionT)
            + trajectoryEndBarycentricCoords[(intersectionVertexOrdinal + 1) % 3] * minIntersectionT,
            0.0f,
            1.0f);
        intersectionBarycentricCoords[(intersectionVertexOrdinal + 1) % 3] = lNext;
        intersectionBarycentricCoords[(intersectionVertexOrdinal + 2) % 3] = 1.0f - lNext;

        {
            assert(intersectionBarycentricCoords[0] >= 0.0f && intersectionBarycentricCoords[0] <= 1.0f);
            assert(intersectionBarycentricCoords[1] >= 0.0f && intersectionBarycentricCoords[1] <= 1.0f);
            assert(intersectionBarycentricCoords[2] >= 0.0f && intersectionBarycentricCoords[2] <= 1.0f);
        }

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

            vec2f const trajectoryEndAbsolutePosition = mesh.GetTriangles().FromBarycentricCoordinates(
                trajectoryEndBarycentricCoords,
                npcParticleConstrainedState.CurrentTriangle,
                mesh.GetVertices());

            vec2f const trajectory = trajectoryEndAbsolutePosition - trajectoryStartAbsolutePosition;

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
                dt,
                particles,
                labParameters);

            // Consume the entire time quantum
            return;
        }
        else
        {
            //
            // Not floor, climb over edge
            //

            LogMessage("      Climbing over non-floor edge");

            // Find opposite triangle
            ElementIndex const oppositeTriangle = mesh.GetEdges().GetOppositeTriangle(intersectionEdgeElementIndex, npcParticleConstrainedState.CurrentTriangle);
            if (oppositeTriangle == NoneElementIndex)
            {
                //
                // Become free
                //

                LogMessage("      No opposite triangle found, becoming free");

                //
                // Move to endpoint and exit, consuming whole quantum
                //

                vec2f const trajectoryEndAbsolutePosition = mesh.GetTriangles().FromBarycentricCoordinates(
                    trajectoryEndBarycentricCoords,
                    npcParticleConstrainedState.CurrentTriangle,
                    mesh.GetVertices());

                npcParticle.ConstrainedState.reset();

                UpdateNpcParticle_Free(
                    npcParticle,
                    particleStartAbsolutePosition,
                    trajectoryEndAbsolutePosition,
                    particles,
                    labParameters);

                return;
            }
            else
            {
                //
                // Move to edge of opposite triangle 
                //

                int const oppositeTriangleEdgeOrdinal = mesh.GetTriangles().GetSubEdgeOrdinal(oppositeTriangle, intersectionEdgeElementIndex);

                LogMessage("      Moving to edge ", oppositeTriangleEdgeOrdinal, " of opposite triangle ", oppositeTriangle);

                // Save current absolute trajectory end
                vec2f const oldTrajectoryEndAbsolutePosition = mesh.GetTriangles().FromBarycentricCoordinates(
                    trajectoryEndBarycentricCoords,
                    npcParticleConstrainedState.CurrentTriangle,
                    mesh.GetVertices());

                npcParticleConstrainedState.CurrentTriangle = oppositeTriangle;

                // Calculate new current barycentric coords (wrt opposite triangle)
                vec3f newBarycentricCoords; // In new triangle
                newBarycentricCoords[(oppositeTriangleEdgeOrdinal + 2) % 3] = 0.0f;
                newBarycentricCoords[oppositeTriangleEdgeOrdinal] = npcParticleConstrainedState.CurrentTriangleBarycentricCoords[(intersectionEdgeOrdinal + 1) % 3];
                newBarycentricCoords[(oppositeTriangleEdgeOrdinal + 1) % 3] = npcParticleConstrainedState.CurrentTriangleBarycentricCoords[intersectionEdgeOrdinal];

                LogMessage("      B-Coords: ", npcParticleConstrainedState.CurrentTriangleBarycentricCoords, " -> ", newBarycentricCoords);

                {
                    assert(newBarycentricCoords[0] >= 0.0f && newBarycentricCoords[0] <= 1.0f);
                    assert(newBarycentricCoords[1] >= 0.0f && newBarycentricCoords[1] <= 1.0f);
                    assert(newBarycentricCoords[2] >= 0.0f && newBarycentricCoords[2] <= 1.0f);
                }

                npcParticleConstrainedState.CurrentTriangleBarycentricCoords = newBarycentricCoords;

                // Translate target coords to this triangle, for next iteration

                auto const oldTrajectoryEndBarycentricCoords = trajectoryEndBarycentricCoords;

                // Note: here we introduce a lot of error - the target bary coords are not anymore
                // guaranteed to lie exactly on the (continuation of the) edge
                trajectoryEndBarycentricCoords = mesh.GetTriangles().ToBarycentricCoordinates(
                    oldTrajectoryEndAbsolutePosition,
                    oppositeTriangle,
                    mesh.GetVertices());

                LogMessage("      TrajEndB-Coords: ", oldTrajectoryEndBarycentricCoords, " -> ", trajectoryEndBarycentricCoords);

                // Continue
            }
        }
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
                    && npc.HumanNpcState->CurrentWalkMagnitude != 0.0f)
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
    Mesh const & /*mesh*/,
    LabParameters const & /*labParameters*/)
{
    if (npc.Type == NpcType::Human && isPrimaryParticle) // TODO: should we invoke it directly on UpdateHumanNpcAnimation?
    {
        assert(npc.DipoleState.has_value());
        assert(npc.HumanNpcState.has_value());

        switch (npc.HumanNpcState->CurrentBehavior)
        {
            case StateType::HumanNpcStateType::BehaviorType::Constrained_Walking:
            {
                ElementIndex const primaryParticleIndex = npc.PrimaryParticleState.ParticleIndex;
                ElementIndex const secondaryParticleIndex = npc.DipoleState->SecondaryParticleState.ParticleIndex;

                npc.HumanNpcState->TopPoint = mParticles.GetPosition(secondaryParticleIndex);
                npc.HumanNpcState->CrotchPoint = mParticles.GetPosition(secondaryParticleIndex) + (mParticles.GetPosition(primaryParticleIndex) - mParticles.GetPosition(secondaryParticleIndex)) * 0.5f;
                npc.HumanNpcState->RightLegPoint = npc.HumanNpcState->CrotchPoint + vec2f(0.5f, -0.5f) * LabParameters::HumanNpcLength;
                npc.HumanNpcState->LeftLegPoint = npc.HumanNpcState->CrotchPoint + vec2f(-0.5f, -0.5f) * LabParameters::HumanNpcLength;

                // TODOHERE

                break;
            }

            default:
            {
                break;
            }
        }
    }
}

Npcs::StateType Npcs::MaterializeNpcState(
    ElementIndex npcIndex,
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
            dipoleState->SecondaryParticleState);
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
    if (triangleIndex != NoneElementIndex)
    {
        vec3f const barycentricCoords = mesh.GetTriangles().ToBarycentricCoordinatesFromWithinTriangle(
            position,
            triangleIndex,
            mesh.GetVertices());

        {
            assert(barycentricCoords[0] >= 0.0f && barycentricCoords[0] <= 1.0f);
            assert(barycentricCoords[1] >= 0.0f && barycentricCoords[1] <= 1.0f);
            assert(barycentricCoords[2] >= 0.0f && barycentricCoords[2] <= 1.0f);
        }

        return Npcs::StateType::NpcParticleStateType::ConstrainedStateType(
            triangleIndex,
            barycentricCoords);
    }

    return std::nullopt;
}
