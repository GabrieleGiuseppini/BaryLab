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
            dt,
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

    for (float remainingDt = dt; ; )
    {
        assert(remainingDt > 0.0f);

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

        // Note: this is the same as physicsDeltaPos + meshDisplacement
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

                        remainingDt = UpdateNpcParticle_ConstrainedNonInertial2(
                            particle,
                            dipoleArg,
                            isPrimaryParticle,
                            edgeOrdinal,
                            edgeDir,
                            trajectoryStartAbsolutePosition,
                            trajectory,
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

                        hasFoundNonInertial = true;

                        break;
                    }
                }
            }
        }

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

            // TODOHERE: will it ever be remainingDt != 0, here that we are in inertial case?
            remainingDt = UpdateNpcParticle_ConstrainedTraceSegment2(
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
    float dt,
    NpcParticles & particles) const
{
    assert(!particle.ConstrainedState.has_value());

    // Update position
    particles.SetPosition(
        particle.ParticleIndex,
        endPosition);

    // Update velocity
    particles.SetVelocity(
        particle.ParticleIndex,
        (endPosition - startPosition) / dt);
}

float Npcs::UpdateNpcParticle_ConstrainedNonInertial2(
    StateType::NpcParticleStateType & particle,
    std::optional<DipoleArg> const & dipoleArg,
    bool isPrimaryParticle,
    int edgeOrdinal,
    vec2f const & edgeDir,
    vec2f const & trajectoryStartAbsolutePosition,
    vec2f const & trajectory,
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

        LogMessage("      friction: trajectoryN=", trajectoryN, " relVel=", particle.ConstrainedState->MeshRelativeVelocity, 
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

    vec2f targetAbsolutePosition = trajectoryStartAbsolutePosition + flattenedTrajectory;

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

    ////    targetAbsolutePosition += walkVector;

    ////    if (isPrimaryParticle)
    ////    {
    ////        // Impart same displacement to secondary particle
    ////        secondaryVoluntarySuperimposedDisplacement = walkVector;
    ////    }
    ////}

    vec3f targetBarycentricCoords = mesh.GetTriangles().ToBarycentricCoordinates(
        targetAbsolutePosition,
        particle.ConstrainedState->CurrentTriangle,
        mesh.GetVertices());

    //
    // Due to numerical slack, ensure target barycentric coords are still along edge
    //

    int const vertexOrdinal = (edgeOrdinal + 2) % 3;
    targetBarycentricCoords[vertexOrdinal] = 0.0f;
    targetBarycentricCoords[(vertexOrdinal + 1) % 3] = 1.0f - targetBarycentricCoords[(vertexOrdinal + 2) % 3];

    LogMessage("      targetAbsolutePosition=", targetAbsolutePosition, " targetBarycentricCoords=", targetBarycentricCoords);

    return UpdateNpcParticle_ConstrainedTraceSegment2(
        particle,
        dipoleArg,
        isPrimaryParticle,
        particle.ConstrainedState->CurrentTriangleBarycentricCoords, // trajectoryStartBarycentricCoords
        targetBarycentricCoords,
        meshVelocity,
        dt,
        particles,
        mesh,
        labParameters);
}

float Npcs::UpdateNpcParticle_ConstrainedTraceSegment2(
    StateType::NpcParticleStateType & particle,
    std::optional<DipoleArg> const & dipoleArg,
    bool isPrimaryParticle,
    vec3f const trajectoryStartBarycentricCoords,
    vec3f trajectoryEndBarycentricCoords,
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
    // 1. Reached destination: simulation step has completed
    // 2. Becoming free: simulation step is half-completed, a portion remains for free movement
    // 3. Impact with bounce: simulation step is half-completed, a portion remains for another constrained step
    //

    for (int i = 0; ; ++i)
    {
        assert(particle.ConstrainedState.has_value());

        {
            assert(particle.ConstrainedState->CurrentTriangleBarycentricCoords[0] >= 0.0f && particle.ConstrainedState->CurrentTriangleBarycentricCoords[0] <= 1.0f);
            assert(particle.ConstrainedState->CurrentTriangleBarycentricCoords[1] >= 0.0f && particle.ConstrainedState->CurrentTriangleBarycentricCoords[1] <= 1.0f);
            assert(particle.ConstrainedState->CurrentTriangleBarycentricCoords[2] >= 0.0f && particle.ConstrainedState->CurrentTriangleBarycentricCoords[2] <= 1.0f);
        }

        LogMessage("    SegmentTrace ", i, ":");
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

            vec2f const absoluteVelocity = (particleEndAbsolutePosition - particleStartAbsolutePosition) / dt;
            particles.SetVelocity(particle.ParticleIndex, absoluteVelocity);
            particle.ConstrainedState->MeshRelativeVelocity = absoluteVelocity + meshVelocity;

            // We have consumed the whole time quantum
            return 0.0f;
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

        //
        // Move to intersection
        //

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

        // Move barycentric coords to intersection

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
            // Impact - calculate bounce response, using the *apparent* (trajectory) 
            // velocity - since this one includes the mesh velocity
            //

            LogMessage("      Impact");

            // Decompose apparent particle velocity into normal and tangential

            vec2f const trajectoryEndAbsolutePosition = mesh.GetTriangles().FromBarycentricCoordinates(
                trajectoryEndBarycentricCoords,
                particle.ConstrainedState->CurrentTriangle,
                mesh.GetVertices());

            vec2f const trajectory = trajectoryEndAbsolutePosition - trajectoryStartAbsolutePosition;
            vec2f const apparentParticleVelocity = trajectory / dt;

            LogMessage("        apparentParticleVelocity=", apparentParticleVelocity);

            vec2f const edgeDir =
                mesh.GetTriangles().GetSubEdgeVector(particle.ConstrainedState->CurrentTriangle, intersectionEdgeOrdinal, mesh.GetVertices())
                .normalise();
            vec2f const edgeNormal = edgeDir.to_perpendicular();
            float const apparentParticleVelocityAlongNormal = apparentParticleVelocity.dot(edgeNormal);
            vec2f const normalVelocity = edgeNormal * apparentParticleVelocityAlongNormal;
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

            particles.SetPosition(particle.ParticleIndex, intersectionAbsolutePosition);

            particles.SetVelocity(particle.ParticleIndex, resultantAbsoluteVelocity);
            particle.ConstrainedState->MeshRelativeVelocity = resultantAbsoluteVelocity + meshVelocity;

            // Consume the entire time quantum
            return 0.0f;
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
                    dt,
                    particles);

                return 0.0f;
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

                // Translate target coords to this triangle

                trajectoryEndBarycentricCoords = mesh.GetTriangles().ToBarycentricCoordinates(
                    oldTrajectoryEndAbsolutePosition,
                    particle.ConstrainedState->CurrentTriangle,
                    mesh.GetVertices());

                // Continue
            }
        }
    }
}