/***************************************************************************************
 * Original Author:		Gabriele Giuseppini
 * Created:				2023-07-23
 * Copyright:			Gabriele Giuseppini  (https://github.com/GabrieleGiuseppini)
 ***************************************************************************************/
#include "Npcs.h"

#include "BLabMath.h"
#include "Log.h"

#include <limits>

void Npcs::UpdateNpcs2(
    Mesh const & mesh,
    LabParameters const & labParameters)
{
    LogMessage("----------------------------------");
    LogMessage("----------------------------------");
    LogMessage("----------------------------------");

    //
    // 1. Check if a free secondary particle should become constrained
    // 2. Update behavioral state machines
    // 3. Calculate spring forces 
    //

    for (auto const n : *this)
    {
        auto & state = mStateBuffer[n];

        // Secondary free becoming constrained
        
        if (state.DipoleState.has_value()
            && !state.DipoleState->SecondaryParticleState.ConstrainedState.has_value() // Secondary is free
            && state.PrimaryParticleState.ConstrainedState.has_value()) // And primary is constrained
        {
            state.DipoleState->SecondaryParticleState.ConstrainedState = CalculateParticleConstrainedState(
                mParticles.GetPosition(state.DipoleState->SecondaryParticleState.ParticleIndex),
                mesh);
        }

        // Behavior

        if (state.Type == NpcType::Human)
        {
            assert(state.DipoleState.has_value());
            assert(state.HumanNpcState.has_value());

            UpdateHuman(
                *state.HumanNpcState,
                state.PrimaryParticleState,
                state.DipoleState->SecondaryParticleState,
                mesh,
                labParameters);
        }

        // Spring forces

        if (state.DipoleState.has_value())
        {
            ElementIndex const primaryParticleIndex = state.PrimaryParticleState.ParticleIndex;
            ElementIndex const secondaryParticleIndex = state.DipoleState->SecondaryParticleState.ParticleIndex;

            float const dt = LabParameters::SimulationTimeStepDuration;

            vec2f const springDisplacement = mParticles.GetPosition(secondaryParticleIndex) - mParticles.GetPosition(primaryParticleIndex); // Towards secondary particle
            float const springDisplacementLength = springDisplacement.length();
            vec2f const springDir = springDisplacement.normalise_approx(springDisplacementLength);

            //
            // 1. Hooke's law
            //

            float const springStiffnessCoefficient =
                labParameters.SpringReductionFraction
                * state.DipoleState->DipoleProperties.MassFactor * labParameters.MassAdjustment
                / (dt * dt);

            // Calculate spring force on this particle
            float const fSpring =
                (springDisplacementLength - state.DipoleState->DipoleProperties.DipoleLength)
                * springStiffnessCoefficient;

            //
            // 2. Damper forces
            //
            // Damp the velocities of each endpoint pair, as if the points were also connected by a damper
            // along the same direction as the spring
            //

            float const springDampingCoefficient =
                labParameters.SpringDampingCoefficient
                * state.DipoleState->DipoleProperties.MassFactor * labParameters.MassAdjustment
                / dt;

            // Calculate damp force on this particle
            vec2f const relVelocity = mParticles.GetVelocity(secondaryParticleIndex) - mParticles.GetVelocity(primaryParticleIndex);
            float const fDamp =
                relVelocity.dot(springDir)
                * springDampingCoefficient;

            //
            // 3. Apply forces
            //

            vec2f const primaryParticlesSpringForces = springDir * (fSpring + fDamp);

            mParticles.SetSpringForces(
                primaryParticleIndex,
                primaryParticlesSpringForces);

            mParticles.SetSpringForces(
                secondaryParticleIndex,
                -primaryParticlesSpringForces);
        }
    }

    //
    // 3. Update state
    //

    for (auto const n : *this)
    {
        LogMessage("NPC ", n);
        LogMessage("----------------------------------");

        auto & npcState = mStateBuffer[n];

        UpdateNpcParticle2(
            npcState.PrimaryParticleState,
            npcState.DipoleState.has_value()
                ? DipoleArg(npcState.DipoleState->SecondaryParticleState, npcState.DipoleState->DipoleProperties)
                : std::optional<DipoleArg>(),
            true,
            npcState,
            mesh,
            labParameters);

        if (npcState.DipoleState.has_value())
        {
            UpdateNpcParticle2(
                npcState.DipoleState->SecondaryParticleState,
                DipoleArg(npcState.PrimaryParticleState, npcState.DipoleState->DipoleProperties),
                false,
                npcState,
                mesh,
                labParameters);
        }
    }
}

void Npcs::UpdateNpcParticle2(
    StateType::NpcParticleStateType & particle,
    std::optional<DipoleArg> const & dipoleArg,
    bool isPrimaryParticle,
    StateType & npc,
    Mesh const & mesh,
    LabParameters const & labParameters)
{
    // TODOHERE: might be needed later for walking
    (void)npc;

    LogMessage("  ", isPrimaryParticle ? "Primary" : "Secondary");

    float const particleMass = LabParameters::ParticleMass * labParameters.MassAdjustment;
    float const dt = LabParameters::SimulationTimeStepDuration;

    vec2f const particleStartAbsolutePosition = mParticles.GetPosition(particle.ParticleIndex);

    // Calculate physical displacement - once and for all, as whole loop
    // will attempt to move to trajectory that always ends here

    vec2f const physicalForces =
        mParticles.GetWorldForce(particle.ParticleIndex)
        + LabParameters::Gravity * labParameters.GravityAdjustment * mGravityGate * particleMass
        + mParticles.GetSpringForces(particle.ParticleIndex);

    vec2f const physicsDeltaPos =
        mParticles.GetVelocity(particle.ParticleIndex) * dt
        + (physicalForces / particleMass) * dt * dt;

    // Absolute position of particle if it only moved due to physical forces; constant
    // during the whole loop
    vec2f const trajectoryEndAbsolutePosition = particleStartAbsolutePosition + physicsDeltaPos;

    if (!particle.ConstrainedState.has_value())
    {
        // 
        // Free
        //

        LogMessage("    Free: physicsDeltaPos=", physicsDeltaPos);
        LogMessage("    StartPosition=", particleStartAbsolutePosition, " StartVelocity=", mParticles.GetVelocity(particle.ParticleIndex));

        UpdateNpcParticle_Free2(
            particle,
            particleStartAbsolutePosition,
            trajectoryEndAbsolutePosition,
            mParticles);

        LogMessage("    EndPosition=", mParticles.GetPosition(particle.ParticleIndex), " EndVelocity=", mParticles.GetVelocity(particle.ParticleIndex));

        // We're done
        return;
    }

    //
    // Constrained
    //

    assert(particle.ConstrainedState.has_value());

    // Loop tracing trajectory from TrajectoryStart (== current bary coords) to TrajectoryEnd (== start absolute pos + deltaPos); 
    // each step moves the next TrajectoryStart a bit ahead.
    // Each iteration of the loop either exits (completes), or moves currenty bary coords (and calcs remaining dt) when it wants 
    // to "continue" an impact while on edge-moving-against-it, i.e. when it wants to recalculate a new flattened traj.
    //    - In this case, the iteration doesn't change current absolute position, nor velocity; it only updates current bary coords 
    //      to what will become the next TrajectoryStart
    //    - In this case, at next iteration:
    //          - TrajectoryStart (== current bary coords) is new
    //          - TrajectoryEnd (== start absolute pos + physicsDeltaPos) is same as before

    ElementIndex currentTriangleElementIndex = particle.ConstrainedState->CurrentTriangle;

    // Initialize absolute position of particle (wrt current triangle) as if it moved with the mesh, staying in its position wrt its triangle;
    // it's the new theoretical position after mesh displacement
    vec2f trajectoryStartAbsolutePosition = mesh.GetTriangles().FromBarycentricCoordinates(
        particle.ConstrainedState->CurrentTriangleBarycentricCoords,
        currentTriangleElementIndex,
        mesh.GetVertices());

    // Calculate mesh displacement for the whole loop
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

    for (float remainingDt = dt; ; )
    {
        assert(remainingDt > 0.0f);

        LogMessage("    ------------------------");
        LogMessage("    remainingDt=", remainingDt);

        //
        // We ray-trace the particle along a trajectory that starts at the position at which the particle
        // would be if it moved with the mesh (i.e. staying at its position wrt its triangle) and ends
        // at the position at which the particle would be if it were only subject to physical forces
        //

        //
        // Calculate trajectory
        //
        // Trajectory is the apparent displacement, i.e. the displacement of the particle
        // *from the point of view of the mesh* (or of the particle itself).
        //
        // Given that trajectory discounts physics move, trajectory is the displacement
        // caused by the apparent forces. In fact, we've verified that when the particle has
        // the same velocity as the mesh, trajectory is zero.
        //
        // Note: when there's no physics displacement, this amounts to pure mesh move
        // (PartPos - FromBary(PartPosBary)). On the other hand, when the mesh is at rest, this
        // amounts to pure physical displacement.
        //

        // Note: on first iteration this is the same as physicsDeltaPos + meshDisplacement
        vec2f const trajectory = trajectoryEndAbsolutePosition - trajectoryStartAbsolutePosition;

        //
        // Check which of two cases applies at this moment:
        // 1. On edge and moving against (or along) it (i.e. in a non-inertial frame)
        // 2. Not on edge or on edge but not moving against (nor along) it
        //
        
        // Check all edges, stop at first one that is floor and against which we're moving
        bool hasFoundNonInertial = false;
        for (int vi = 0; vi < 3; ++vi)
        {
            if (particle.ConstrainedState->CurrentTriangleBarycentricCoords[vi] == 0.0f)
            {
                // We are on this edge

                int const edgeOrdinal = (vi + 1) % 3;
                ElementIndex const currentEdgeElementIndex = mesh.GetTriangles().GetSubEdges(currentTriangleElementIndex).EdgeIndices[edgeOrdinal];

                assert(isPrimaryParticle || dipoleArg.has_value());

                LogMessage("      edge ", edgeOrdinal, ": isFloor=", IsEdgeFloorToParticle(currentEdgeElementIndex, currentTriangleElementIndex, mesh));

                // Check if this is really a floor to this particle
                if (IsEdgeFloorToParticle(currentEdgeElementIndex, currentTriangleElementIndex, mesh)
                    && (isPrimaryParticle || !DoesFloorSeparateFromPrimaryParticle(
                        mParticles.GetPosition(dipoleArg->OtherParticle.ParticleIndex),
                        trajectoryStartAbsolutePosition, // == New theoretical position after mesh displacement
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

                        LogMessage("    ConstrainedNonInertial: triangle=", currentTriangleElementIndex, " edgeOrdinal=", edgeOrdinal, " bCoords=", particle.ConstrainedState->CurrentTriangleBarycentricCoords, " trajectory=", trajectory);
                        LogMessage("    StartPosition=", mParticles.GetPosition(particle.ParticleIndex), " StartVelocity=", mParticles.GetVelocity(particle.ParticleIndex), " MeshVelocity=", meshVelocity, " StartMRVelocity=", particle.ConstrainedState->MeshRelativeVelocity);

                        // Ray-trace using non-inertial physics;
                        // will return when completed or when current edge is over
                        float const newRemainingDt = UpdateNpcParticle_ConstrainedNonInertial2(
                            particle,
                            dipoleArg,
                            isPrimaryParticle,
                            edgeOrdinal,
                            edgeDir,
                            particleStartAbsolutePosition,
                            trajectoryStartAbsolutePosition,
                            trajectory,
                            meshVelocity,
                            remainingDt,
                            mParticles,
                            mesh,
                            labParameters);

                        if (newRemainingDt == remainingDt)
                        {
                            // No movement

                            // Check if we're in a well
                            if (pastPastBarycentricPosition.has_value()
                                && pastPastBarycentricPosition->TriangleElementIndex == particle.ConstrainedState->CurrentTriangle
                                && pastPastBarycentricPosition->BarycentricCoords == particle.ConstrainedState->CurrentTriangleBarycentricCoords)
                            {
                                //
                                // Well - stop here
                                //

                                LogMessage("    Detected well - stopping here");

                                // Update particle's physics, considering that we are in a well and thus still (wrt mesh)

                                vec2f const particleEndAbsolutePosition = mesh.GetTriangles().FromBarycentricCoordinates(
                                    particle.ConstrainedState->CurrentTriangleBarycentricCoords,
                                    particle.ConstrainedState->CurrentTriangle,
                                    mesh.GetVertices());

                                mParticles.SetPosition(particle.ParticleIndex, particleEndAbsolutePosition);

                                // No (relative) velocity
                                mParticles.SetVelocity(particle.ParticleIndex, -meshVelocity);
                                particle.ConstrainedState->MeshRelativeVelocity = vec2f::zero();

                                // Consume the whole time quantum
                                remainingDt = 0.0f;
                            }
                        }
                        else
                        {
                            // We have moved
                            remainingDt = newRemainingDt;

                            // Reset well detection machinery
                            pastPastBarycentricPosition.reset();
                            pastBarycentricPosition.reset();
                        }
                            
                        // Update well detection machinery
                        pastPastBarycentricPosition = pastBarycentricPosition;
                        pastBarycentricPosition.emplace(particle.ConstrainedState->CurrentTriangle, particle.ConstrainedState->CurrentTriangleBarycentricCoords);

                        if (particle.ConstrainedState.has_value())
                        {
                            LogMessage("    EndPosition=", mParticles.GetPosition(particle.ParticleIndex), " EndVelocity=", mParticles.GetVelocity(particle.ParticleIndex), " EndMRVelocity=", particle.ConstrainedState->MeshRelativeVelocity);
                        }
                        else
                        {
                            LogMessage("    EndPosition=", mParticles.GetPosition(particle.ParticleIndex), " EndVelocity=", mParticles.GetVelocity(particle.ParticleIndex));
                        }

                        hasFoundNonInertial = true;

                        break;
                    }
                    else
                    {
                        LogMessage("      Traj.EdgeN=", trajectory.dot(edgeNormal));
                    }
                }
            }
        } // for (vertex)

        if (!hasFoundNonInertial)
        {
            //
            // Case 2: Inertial: not on edge or on edge but not moving against (nor along) it
            //

            LogMessage("    ConstrainedInertial: triangle=", currentTriangleElementIndex, " bCoords=", particle.ConstrainedState->CurrentTriangleBarycentricCoords, " physicsDeltaPos=", physicsDeltaPos);
            LogMessage("    StartPosition=", mParticles.GetPosition(particle.ParticleIndex), " StartVelocity=", mParticles.GetVelocity(particle.ParticleIndex), " MeshVelocity=", meshVelocity, " StartMRVelocity=", particle.ConstrainedState->MeshRelativeVelocity);

            //
            // Calculate target barycentric coords
            //

            vec3f const trajectoryEndBarycentricCoords = mesh.GetTriangles().ToBarycentricCoordinates(
                trajectoryEndAbsolutePosition,
                particle.ConstrainedState->CurrentTriangle,
                mesh.GetVertices());

            //
            // Move towards target bary coords
            //

            UpdateNpcParticle_ConstrainedInertial2(
                particle,
                dipoleArg,
                isPrimaryParticle,
                particle.ConstrainedState->CurrentTriangleBarycentricCoords, // trajectoryStartBarycentricCoords
                trajectoryEndBarycentricCoords,
                meshVelocity,
                remainingDt,
                mParticles,
                mesh,
                labParameters);

            if (particle.ConstrainedState.has_value())
            {
                LogMessage("    EndPosition=", mParticles.GetPosition(particle.ParticleIndex), " EndVelocity=", mParticles.GetVelocity(particle.ParticleIndex), " EndMRVelocity=", particle.ConstrainedState->MeshRelativeVelocity);
            }
            else
            {
                LogMessage("    EndPosition=", mParticles.GetPosition(particle.ParticleIndex), " EndVelocity=", mParticles.GetVelocity(particle.ParticleIndex));
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

        currentTriangleElementIndex = particle.ConstrainedState->CurrentTriangle;

        trajectoryStartAbsolutePosition = mesh.GetTriangles().FromBarycentricCoordinates(
            particle.ConstrainedState->CurrentTriangleBarycentricCoords,
            currentTriangleElementIndex,
            mesh.GetVertices());
    }
    
    // Publish final velocities

    vec2f const particleVelocity = (mParticles.GetPosition(particle.ParticleIndex) - particleStartAbsolutePosition) / LabParameters::SimulationTimeStepDuration;

    mEventDispatcher.OnCustomProbe("VelX", particleVelocity.x);
    mEventDispatcher.OnCustomProbe("VelY", particleVelocity.y);
}

void Npcs::UpdateNpcParticle_Free2(
    StateType::NpcParticleStateType & particle,
    vec2f const & startPosition,
    vec2f const & endPosition,
    NpcParticles & particles) const
{
    assert(!particle.ConstrainedState.has_value());

    // Update position
    particles.SetPosition(
        particle.ParticleIndex,
        endPosition);

    // Update velocity
    // Use whole time quantum for velocity, as particleStartAbsolutePosition is fixed at t0
    particles.SetVelocity(
        particle.ParticleIndex,
        (endPosition - startPosition) / LabParameters::SimulationTimeStepDuration);
}

float Npcs::UpdateNpcParticle_ConstrainedNonInertial2(
    StateType::NpcParticleStateType & particle,
    std::optional<DipoleArg> const & dipoleArg,
    bool const isPrimaryParticle,
    int edgeOrdinal,
    vec2f const & edgeDir,
    vec2f const & particleStartAbsolutePosition,
    vec2f const & trajectoryStartAbsolutePosition,
    vec2f const & theoreticalTrajectory,
    vec2f const meshVelocity,
    float dt,
    NpcParticles & particles,
    Mesh const & mesh,
    LabParameters const & labParameters) const
{
    assert(particle.ConstrainedState.has_value());

    //
    // We're moving against the floor, hence we are in a non-inertial frame...
    // ...take friction into account and flaten trajectory
    //

    //
    // Calculate magnitude of flattened trajectory - i.e. component of trajectory
    // along (tangent) edge, positive when in the direction of the edge
    //

    float trajectoryT = theoreticalTrajectory.dot(edgeDir);

    //
    // Update tangential component of trajectory with friction
    //

    {
        // Normal trajectory: apparent (integrated) force against the floor;
        // positive when against the floor

        float const trajectoryN = theoreticalTrajectory.dot(edgeDir.to_perpendicular());

        //
        // Choose between kinetic and static friction
        //
        // Note that we want to check actual velocity here, not
        // *intention* to move (which is trajectory)
        //

        float frictionCoefficient;
        if (std::abs(particle.ConstrainedState->MeshRelativeVelocity.dot(edgeDir)) > 0.01f) // Magic number
        {
            // Kinetic friction
            frictionCoefficient = labParameters.KineticFriction;
        }
        else
        {
            // Static friction
            frictionCoefficient = labParameters.StaticFriction;
        }

        // Calculate friction (integrated) force magnitude (along edgeDir),
        // which is the same as apparent force, up to max friction threshold
        float tFriction = std::min(std::abs(trajectoryT), frictionCoefficient * std::max(trajectoryN, 0.0f));
        if (trajectoryT >= 0.0f)
        {
            tFriction *= -1.0f;
        }

        LogMessage("        friction: trajectoryN=", trajectoryN, " relVel=", particle.ConstrainedState->MeshRelativeVelocity, 
            " trajectoryT=", trajectoryT, " tFriction=", tFriction);

        // Update trajectory with friction
        trajectoryT += tFriction;
    }

    //
    // Recovered flattened trajectory as a vector
    //

    vec2f const flattenedTrajectory = edgeDir * trajectoryT;

    //
    // Calculate trajectory target
    //

    vec2f trajectoryEndAbsolutePosition = trajectoryStartAbsolutePosition + flattenedTrajectory;

    ////// TODOHERE: walking
    ////// Adjust target pos with superimposed displacement - which is NOT taken into account for velocity
    ////std::optional<vec2f> secondaryVoluntarySuperimposedDisplacement;
    ////if (npc.HumanNpcState.has_value())
    ////{
    ////    vec2f const walkVector = edgeDir * edgeDir.dot(vec2f(1.0f, 0.0f)) * mParticles.GetVoluntarySuperimposedDisplacement(particle.ParticleIndex);
    ////    if (walkVector != vec2f::zero())
    ////    {
    ////        LogMessage("WalkVector@1: ", walkVector, " ", npc.HumanNpcState->CurrentWalkingMagnitude);
    ////    }

    ////    trajectoryEndAbsolutePosition += walkVector;

    ////    if (isPrimaryParticle)
    ////    {
    ////        // Impart same displacement to secondary particle
    ////        secondaryVoluntarySuperimposedDisplacement = walkVector;
    ////    }
    ////}

    vec3f trajectoryEndBarycentricCoords = mesh.GetTriangles().ToBarycentricCoordinates(
        trajectoryEndAbsolutePosition,
        particle.ConstrainedState->CurrentTriangle,
        mesh.GetVertices());

    //
    // Due to numerical slack, ensure target barycentric coords are still along edge
    //

    int const edgeVertexOrdinal = (edgeOrdinal + 2) % 3;
    trajectoryEndBarycentricCoords[edgeVertexOrdinal] = 0.0f;
    trajectoryEndBarycentricCoords[(edgeVertexOrdinal + 1) % 3] = 1.0f - trajectoryEndBarycentricCoords[(edgeVertexOrdinal + 2) % 3];

    LogMessage("        trajectoryEndAbsolutePosition=", trajectoryEndAbsolutePosition, " trajectoryEndBarycentricCoords=", trajectoryEndBarycentricCoords);

    //
    // Ray-trace along the flattened trajectory, ending only at one of the following three conditions:
    // 1. Reached destination: terminate
    // 2. Becoming free: do free movement and terminate
    // 3. Impact:
    //      - With bounce: impart bounce velocity and terminate
    //      - Without bounce: advance simulation by how much walked, re-flatten trajectory, and continue
    //

    float distanceTraveledAlongEdge = 0.0f;

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

        particle.ConstrainedState->CurrentTriangleBarycentricCoords = trajectoryEndBarycentricCoords;

        vec2f const particleEndAbsolutePosition = mesh.GetTriangles().FromBarycentricCoordinates(
            trajectoryEndBarycentricCoords,
            particle.ConstrainedState->CurrentTriangle,
            mesh.GetVertices());

        particles.SetPosition(particle.ParticleIndex, particleEndAbsolutePosition);

        // Use whole time quantum for velocity, as particleStartAbsolutePosition is fixed at t0        
        vec2f const absoluteVelocity = (particleEndAbsolutePosition - particleStartAbsolutePosition) / LabParameters::SimulationTimeStepDuration;
        particles.SetVelocity(particle.ParticleIndex, absoluteVelocity);
        particle.ConstrainedState->MeshRelativeVelocity = absoluteVelocity + meshVelocity;

        // We have consumed the whole time quantum
        return 0.0f;
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

    assert(particle.ConstrainedState->CurrentTriangleBarycentricCoords[(edgeOrdinal + 2) % 3] == 0.0f);

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

    for (vec2f lastTrajectoryStartAbsolutePosition = trajectoryStartAbsolutePosition; ;)
    {
        //
        // Move to intersection (which is a vertex by now)
        //

        vec3f intersectionBarycentricCoords = vec3f::zero();
        intersectionBarycentricCoords[intersectionVertexOrdinal] = 1.0f;

        LogMessage("      Moving to intersection with edge ordinal ", intersectionEdgeOrdinal, ": ", intersectionBarycentricCoords);

        particle.ConstrainedState->CurrentTriangleBarycentricCoords = intersectionBarycentricCoords;

        // Update distance traveled

        vec2f const intersectionAbsolutePosition = mesh.GetVertices().GetPosition(mesh.GetTriangles().GetVertexIndices(particle.ConstrainedState->CurrentTriangle)[intersectionVertexOrdinal]);
        assert(intersectionAbsolutePosition == mesh.GetTriangles().FromBarycentricCoordinates(
            intersectionBarycentricCoords,
            particle.ConstrainedState->CurrentTriangle,
            mesh.GetVertices()));

        distanceTraveledAlongEdge += (intersectionAbsolutePosition - lastTrajectoryStartAbsolutePosition).length();

        LogMessage("        distanceTraveledAlongEdge=", distanceTraveledAlongEdge);

        lastTrajectoryStartAbsolutePosition = intersectionAbsolutePosition;

        //
        // Check if impacted with floor
        //

        ElementIndex const intersectionEdgeElementIndex = mesh.GetTriangles().GetSubEdges(particle.ConstrainedState->CurrentTriangle).EdgeIndices[intersectionEdgeOrdinal];

        assert(isPrimaryParticle || dipoleArg.has_value());

        if (IsEdgeFloorToParticle(intersectionEdgeElementIndex, particle.ConstrainedState->CurrentTriangle, mesh)
            && (isPrimaryParticle || !DoesFloorSeparateFromPrimaryParticle(
                mParticles.GetPosition(dipoleArg->OtherParticle.ParticleIndex),
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

            float const trajectoryLength = flattenedTrajectory.length();

            vec2f const intersectionEdgeDir =
                mesh.GetTriangles().GetSubEdgeVector(particle.ConstrainedState->CurrentTriangle, intersectionEdgeOrdinal, mesh.GetVertices())
                .normalise();
            vec2f const intersectionEdgeNormal = intersectionEdgeDir.to_perpendicular();

            // Check angle between desired (original) trajectory and edge
            float const trajProjOntoEdgeNormal = flattenedTrajectory.normalise_approx(trajectoryLength).dot(intersectionEdgeNormal);
            if (trajProjOntoEdgeNormal <= 0.71f) // PI/4+
            {
                //
                // Impact continuation (no bounce)
                //
                // Stop here and then check trajectory in new situation
                //

                LogMessage("      Impact continuation (trajProjOntoEdgeNormal=", trajProjOntoEdgeNormal, ")");

                // Calculate consumed time quantum
                float consumedDt;
                if (trajectoryLength != 0.0f)
                {
                    consumedDt = dt * distanceTraveledAlongEdge / trajectoryLength;
                }
                else
                {
                    // If no distance was to be traveled in this time quantum, then we've consumed it all
                    consumedDt = dt;
                }

                LogMessage("        consumedDt=", consumedDt);

                return dt - consumedDt;
            }
            else
            {
                //
                // Bounce - calculate bounce response, using the *apparent* (trajectory) 
                // velocity - since this one includes the mesh velocity
                //

                LogMessage("      Bounce (trajProjOntoEdgeNormal=", trajProjOntoEdgeNormal, ")");

                BounceConstrainedNpcParticle(
                    particle.ParticleIndex,
                    *particle.ConstrainedState,
                    flattenedTrajectory,
                    intersectionAbsolutePosition,
                    intersectionEdgeNormal,
                    meshVelocity,
                    dt,
                    particles,
                    labParameters);

                // Consume the entire time quantum
                return 0.0f;
            }
        }

        //
        // Not floor, climb over edge
        //

        LogMessage("      Climbing over non-floor edge");

        // Find opposite triangle
        ElementIndex const oppositeTriangle = mesh.GetEdges().GetOppositeTriangle(intersectionEdgeElementIndex, particle.ConstrainedState->CurrentTriangle);
        if (oppositeTriangle == NoneElementIndex)
        {
            //
            // Become free
            //

            LogMessage("      No opposite triangle found, becoming free");

            //
            // Move to endpoint and exit, consuming whole quantum
            //

            particle.ConstrainedState.reset();

            UpdateNpcParticle_Free2(
                particle,
                particleStartAbsolutePosition,
                mesh.GetTriangles().FromBarycentricCoordinates(
                    trajectoryEndBarycentricCoords,
                    particle.ConstrainedState->CurrentTriangle,
                    mesh.GetVertices()),
                particles);

            return 0.0f;
        }

        //
        // Move to edge of opposite triangle 
        //

        int const oppositeTriangleIntersectionEdgeOrdinal = mesh.GetTriangles().GetSubEdgeOrdinal(oppositeTriangle, intersectionEdgeElementIndex);

        LogMessage("      Moving to edge ", oppositeTriangleIntersectionEdgeOrdinal, " of opposite triangle ", oppositeTriangle);

        // Save current absolute trajectory end
        vec2f const oldTrajectoryEndAbsolutePosition = mesh.GetTriangles().FromBarycentricCoordinates(
            trajectoryEndBarycentricCoords,
            particle.ConstrainedState->CurrentTriangle,
            mesh.GetVertices());

        particle.ConstrainedState->CurrentTriangle = oppositeTriangle;

        // Calculate new current barycentric coords (wrt opposite triangle)
        vec3f newBarycentricCoords; // In new triangle
        newBarycentricCoords[(oppositeTriangleIntersectionEdgeOrdinal + 2) % 3] = 0.0f;
        newBarycentricCoords[oppositeTriangleIntersectionEdgeOrdinal] = particle.ConstrainedState->CurrentTriangleBarycentricCoords[(intersectionEdgeOrdinal + 1) % 3];
        newBarycentricCoords[(oppositeTriangleIntersectionEdgeOrdinal + 1) % 3] = particle.ConstrainedState->CurrentTriangleBarycentricCoords[intersectionEdgeOrdinal];

        LogMessage("        B-Coords: ", particle.ConstrainedState->CurrentTriangleBarycentricCoords, " -> ", newBarycentricCoords);

        {
            assert(newBarycentricCoords[0] >= 0.0f && newBarycentricCoords[0] <= 1.0f);
            assert(newBarycentricCoords[1] >= 0.0f && newBarycentricCoords[1] <= 1.0f);
            assert(newBarycentricCoords[2] >= 0.0f && newBarycentricCoords[2] <= 1.0f);
        }

        particle.ConstrainedState->CurrentTriangleBarycentricCoords = newBarycentricCoords;

        //
        // Translate trajectory end coords to this triangle, for next iteration
        //

        auto const oldTrajectoryEndBarycentricCoords = trajectoryEndBarycentricCoords;

        // TODOHERE: here we introduce a lot of error
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

            // Calculate consumed time quantum
            float consumedDt;
            float const trajectoryLength = flattenedTrajectory.length();
            if (trajectoryLength != 0.0f)
            {
                consumedDt = dt * distanceTraveledAlongEdge / trajectoryLength;
            }
            else
            {
                // If no distance was to be traveled in this time quantum, then we've consumed it all
                consumedDt = dt;
            }

            LogMessage("        consumedDt=", consumedDt);

            return dt - consumedDt;
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

void Npcs::UpdateNpcParticle_ConstrainedInertial2(
    StateType::NpcParticleStateType & particle,
    std::optional<DipoleArg> const & dipoleArg,
    bool const isPrimaryParticle,
    vec3f const trajectoryStartBarycentricCoords,
    vec3f trajectoryEndBarycentricCoords, // Mutable    
    vec2f const meshVelocity,
    float dt,
    NpcParticles & particles,
    Mesh const & mesh,
    LabParameters const & labParameters) const
{
    vec2f const particleStartAbsolutePosition = particles.GetPosition(particle.ParticleIndex);

    vec2f const trajectoryStartAbsolutePosition = mesh.GetTriangles().FromBarycentricCoordinates(
        trajectoryStartBarycentricCoords,
        particle.ConstrainedState->CurrentTriangle,
        mesh.GetVertices());

    //
    // Ray-trace along the specified trajectory, ending only at one of the following three conditions:
    // 1. Reached destination: terminate
    // 2. Becoming free: do free movement and terminate
    // 3. Impact with bounce: impart bounce velocity and terminate
    //

    for (int iIter = 0; ; ++iIter)
    {
        assert(particle.ConstrainedState.has_value());

        {
            assert(particle.ConstrainedState->CurrentTriangleBarycentricCoords[0] >= 0.0f && particle.ConstrainedState->CurrentTriangleBarycentricCoords[0] <= 1.0f);
            assert(particle.ConstrainedState->CurrentTriangleBarycentricCoords[1] >= 0.0f && particle.ConstrainedState->CurrentTriangleBarycentricCoords[1] <= 1.0f);
            assert(particle.ConstrainedState->CurrentTriangleBarycentricCoords[2] >= 0.0f && particle.ConstrainedState->CurrentTriangleBarycentricCoords[2] <= 1.0f);
        }

        LogMessage("    SegmentTrace ", iIter);
        LogMessage("      triangle=", particle.ConstrainedState->CurrentTriangle, " bCoords=", particle.ConstrainedState->CurrentTriangleBarycentricCoords,
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

            particle.ConstrainedState->CurrentTriangleBarycentricCoords = trajectoryEndBarycentricCoords;

            vec2f const particleEndAbsolutePosition = mesh.GetTriangles().FromBarycentricCoordinates(
                trajectoryEndBarycentricCoords,
                particle.ConstrainedState->CurrentTriangle,
                mesh.GetVertices());

            particles.SetPosition(particle.ParticleIndex, particleEndAbsolutePosition);

            // Use whole time quantum for velocity, as particleStartAbsolutePosition is fixed at t0
            vec2f const absoluteVelocity = (particleEndAbsolutePosition - particleStartAbsolutePosition) / LabParameters::SimulationTimeStepDuration;
            particles.SetVelocity(particle.ParticleIndex, absoluteVelocity);
            particle.ConstrainedState->MeshRelativeVelocity = absoluteVelocity + meshVelocity;

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
                float const den = particle.ConstrainedState->CurrentTriangleBarycentricCoords[vi] - trajectoryEndBarycentricCoords[vi];
                float const t = particle.ConstrainedState->CurrentTriangleBarycentricCoords[vi] / den;

#ifdef _DEBUG
                diags[vi].emplace(den, t);
                diags[vi]->IntersectionPoint =
                    particle.ConstrainedState->CurrentTriangleBarycentricCoords
                    + (trajectoryEndBarycentricCoords - particle.ConstrainedState->CurrentTriangleBarycentricCoords) * t;
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
        ElementIndex const intersectionEdgeElementIndex = mesh.GetTriangles().GetSubEdges(particle.ConstrainedState->CurrentTriangle).EdgeIndices[intersectionEdgeOrdinal];

        // Calculate intersection barycentric coordinates

        vec3f intersectionBarycentricCoords;
        intersectionBarycentricCoords[intersectionVertexOrdinal] = 0.0f;
        float const lNext = Clamp( // Barycentric coord of next vertex at intersection; enforcing it's within triangle
            particle.ConstrainedState->CurrentTriangleBarycentricCoords[(intersectionVertexOrdinal + 1) % 3] * (1.0f - minIntersectionT)
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

        particle.ConstrainedState->CurrentTriangleBarycentricCoords = intersectionBarycentricCoords;

        //
        // Check if impacted with floor
        //

        vec2f const intersectionAbsolutePosition = mesh.GetTriangles().FromBarycentricCoordinates(
            intersectionBarycentricCoords,
            particle.ConstrainedState->CurrentTriangle,
            mesh.GetVertices());

        assert(isPrimaryParticle || dipoleArg.has_value());

        if (IsEdgeFloorToParticle(intersectionEdgeElementIndex, particle.ConstrainedState->CurrentTriangle, mesh)
            && (isPrimaryParticle || !DoesFloorSeparateFromPrimaryParticle(
                mParticles.GetPosition(dipoleArg->OtherParticle.ParticleIndex),
                intersectionAbsolutePosition,
                intersectionEdgeElementIndex,
                mesh)))
        {
            //
            // Impact and bounce
            //

            LogMessage("      Impact and bounce");

            // Move to intersection, by moving barycentric coords

            LogMessage("      Moving bary coords to intersection with edge ", intersectionEdgeOrdinal, " ", intersectionBarycentricCoords);

            particle.ConstrainedState->CurrentTriangleBarycentricCoords = intersectionBarycentricCoords;

            //
            // Calculate bounce response, using the *apparent* (trajectory) 
            // velocity - since this one includes the mesh velocity
            //

            vec2f const trajectoryEndAbsolutePosition = mesh.GetTriangles().FromBarycentricCoordinates(
                trajectoryEndBarycentricCoords,
                particle.ConstrainedState->CurrentTriangle,
                mesh.GetVertices());

            vec2f const trajectory = trajectoryEndAbsolutePosition - trajectoryStartAbsolutePosition;

            vec2f const intersectionEdgeDir =
                mesh.GetTriangles().GetSubEdgeVector(particle.ConstrainedState->CurrentTriangle, intersectionEdgeOrdinal, mesh.GetVertices())
                .normalise();
            vec2f const intersectionEdgeNormal = intersectionEdgeDir.to_perpendicular();

            BounceConstrainedNpcParticle(
                particle.ParticleIndex,
                *particle.ConstrainedState,
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
            ElementIndex const oppositeTriangle = mesh.GetEdges().GetOppositeTriangle(intersectionEdgeElementIndex, particle.ConstrainedState->CurrentTriangle);
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
                    particle.ConstrainedState->CurrentTriangle,
                    mesh.GetVertices());

                particle.ConstrainedState.reset();

                UpdateNpcParticle_Free2(
                    particle,
                    particleStartAbsolutePosition,
                    trajectoryEndAbsolutePosition,
                    particles);

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
                    particle.ConstrainedState->CurrentTriangle,
                    mesh.GetVertices());

                particle.ConstrainedState->CurrentTriangle = oppositeTriangle;

                // Calculate new current barycentric coords (wrt opposite triangle)
                vec3f newBarycentricCoords; // In new triangle
                newBarycentricCoords[(oppositeTriangleEdgeOrdinal + 2) % 3] = 0.0f;
                newBarycentricCoords[oppositeTriangleEdgeOrdinal] = particle.ConstrainedState->CurrentTriangleBarycentricCoords[(intersectionEdgeOrdinal + 1) % 3];
                newBarycentricCoords[(oppositeTriangleEdgeOrdinal + 1) % 3] = particle.ConstrainedState->CurrentTriangleBarycentricCoords[intersectionEdgeOrdinal];

                LogMessage("      B-Coords: ", particle.ConstrainedState->CurrentTriangleBarycentricCoords, " -> ", newBarycentricCoords);

                {
                    assert(newBarycentricCoords[0] >= 0.0f && newBarycentricCoords[0] <= 1.0f);
                    assert(newBarycentricCoords[1] >= 0.0f && newBarycentricCoords[1] <= 1.0f);
                    assert(newBarycentricCoords[2] >= 0.0f && newBarycentricCoords[2] <= 1.0f);
                }

                particle.ConstrainedState->CurrentTriangleBarycentricCoords = newBarycentricCoords;

                // Translate target coords to this triangle, for next iteration

                auto const oldTrajectoryEndBarycentricCoords = trajectoryEndBarycentricCoords;

                // TODOHERE: here we introduce a lot of error
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
    ElementIndex particleIndex,
    StateType::NpcParticleStateType::ConstrainedStateType & particleConstrainedState,
    vec2f const & trajectory,
    vec2f const & bouncePosition,
    vec2f const & bounceEdgeNormal,
    vec2f const meshVelocity,
    float dt,
    NpcParticles & particles,
    LabParameters const & labParameters) const
{
    // Decompose apparent particle velocity into normal and tangential

    vec2f const apparentParticleVelocity = trajectory / dt;

    LogMessage("        apparentParticleVelocity=", apparentParticleVelocity);

    float const apparentParticleVelocityAlongNormal = apparentParticleVelocity.dot(bounceEdgeNormal);
    vec2f const normalVelocity = bounceEdgeNormal * apparentParticleVelocityAlongNormal;
    vec2f const tangentialVelocity = apparentParticleVelocity - normalVelocity;

    // Calculate normal reponse: Vn' = -e*Vn (e = elasticity, [0.0 - 1.0])
    vec2f const normalResponse =
        -normalVelocity
        * labParameters.Elasticity;

    // Calculate tangential response: Vt' = a*Vt (a = (1.0-friction), [0.0 - 1.0])
    vec2f const tangentialResponse =
        tangentialVelocity
        * (1.0f - labParameters.KineticFriction);

    // Given that we've been working in *apparent* space (we've calc'd the collision response to *trajectory* which is apparent displacement),
    // we need to transform velocity to absolute particle velocity
    vec2f const resultantAbsoluteVelocity = (normalResponse + tangentialResponse) - meshVelocity;

    LogMessage("        nr=", normalResponse, " tr=", tangentialResponse, " rr=", resultantAbsoluteVelocity);

    //
    // Exit, consuming whole time quantum
    //
    // Note: we consume whole time quantum as a subsequent bounce during the same simulation step wouldn't see mesh moving
    //

    particles.SetPosition(particleIndex, bouncePosition);

    particles.SetVelocity(particleIndex, resultantAbsoluteVelocity);
    particleConstrainedState.MeshRelativeVelocity = resultantAbsoluteVelocity + meshVelocity;
}
