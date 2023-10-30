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
    float const particleMass = LabParameters::ParticleMass * labParameters.MassAdjustment;

    for (float remainingDt = LabParameters::SimulationTimeStepDuration; remainingDt > 0.0f; /*reduced in loop*/)
    {
        LogMessage("  ", isPrimaryParticle ? "Primary" : "Secondary", ": remainingDt=", remainingDt);

        //
        // Calculate desired displacement
        //

        vec2f const physicalForces =
            mParticles.GetWorldForce(particle.ParticleIndex)
            + LabParameters::Gravity * labParameters.GravityAdjustment * mGravityGate * labParameters.ParticleMass * labParameters.MassAdjustment
            + mParticles.GetSpringForces(particle.ParticleIndex);

        vec2f const physicsDeltaPos =
            mParticles.GetVelocity(particle.ParticleIndex) * remainingDt
            + (physicalForces / particleMass) * remainingDt * remainingDt;

        //
        // Check which of three cases applies at this moment:
        // 1. Free
        // 2. Constrained, on edge and moving against (or along) it (i.e. in a non-inertial frame)
        // 3. Constrained, not on edge or on edge but not moving against (nor along) it
        //

        if (!particle.ConstrainedState.has_value())
        {
            // 
            // Case 1: Free
            //

            UpdateNpcParticle_Free2(
                particle,
                physicsDeltaPos,
                mParticles);

            // We've consumed the whole time quantum
            remainingDt = 0.0f;
        }
        else
        {
            //
            // Constrained
            //

            ElementIndex const currentTriangleElementIndex = particle.ConstrainedState->CurrentTriangle;

            // New position of particle in moved triangle (i.e. sourcePos + non-inertial displacement)
            // - does not include physics delta pos
            vec2f const newTheoreticalPositionAfterMeshDisplacement = mesh.GetTriangles().FromBarycentricCoordinates(
                particle.ConstrainedState->CurrentTriangleBarycentricCoords,
                currentTriangleElementIndex,
                mesh.GetVertices());

            // Calculate mesh displacement
            vec2f const meshDisplacement = mParticles.GetPosition(particle.ParticleIndex) - newTheoreticalPositionAfterMeshDisplacement;

            //
            // Check whether we're on a floor and wanting to move against it
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
                            newTheoreticalPositionAfterMeshDisplacement,
                            currentEdgeElementIndex,
                            mesh)))
                    {
                        // On floor edge - so potentially in a non-inertial frame

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
                        vec2f const trajectory = mParticles.GetPosition(particle.ParticleIndex) + physicsDeltaPos - newTheoreticalPositionAfterMeshDisplacement;

                        // Check whether we're moving *against* the floor

                        vec2f const edgeVector = mesh.GetTriangles().GetSubEdgeVector(
                            currentTriangleElementIndex,
                            edgeOrdinal,
                            mesh.GetVertices());

                        vec2f const edgeDir = edgeVector.normalise();
                        vec2f const edgeNormal = edgeDir.to_perpendicular(); // Points outside of triangle (i.e. towards floor)

                        if (trajectory.dot(edgeNormal) >= 0.0f) // If 0, no normal force - hence no friction; however we want to take this codepath anyways for consistency
                        {
                            //
                            // Case 2: Constrained, on edge and moving against (or along) it (i.e. in a non-inertial frame)
                            //

                            remainingDt = UpdateNpcParticle_ConstrainedNonInertial2(
                                particle,
                                dipoleArg,
                                isPrimaryParticle,
                                npc,
                                physicsDeltaPos,
                                mParticles,
                                mesh,
                                labParameters);

                            hasFoundNonInertial = true;

                            break;
                        }
                    }
                }
            }

            if (!hasFoundNonInertial)
            {
                //
                // Case 3: Constrained, not on edge or on edge but not moving against (nor along) it
                //

                remainingDt = UpdateNpcParticle_ConstrainedInertial2(
                    particle,
                    dipoleArg,
                    isPrimaryParticle,
                    npc,
                    physicsDeltaPos,
                    mParticles,
                    mesh,
                    labParameters);
            }
        }
    }
}

void Npcs::UpdateNpcParticle_Free2(
    StateType::NpcParticleStateType & particle,
    vec2f const & absoluteDisplacement,
    NpcParticles & particles) const
{
    LogMessage("    Free: absoluteDisplacement=", absoluteDisplacement);

    // Update position
    particles.SetPosition(
        particle.ParticleIndex,
        particles.GetPosition(particle.ParticleIndex) + absoluteDisplacement);

    // Update velocity
    particles.SetVelocity(
        particle.ParticleIndex,
        absoluteDisplacement / LabParameters::SimulationTimeStepDuration);
}

float Npcs::UpdateNpcParticle_ConstrainedNonInertial2(
    StateType::NpcParticleStateType & particle,
    std::optional<DipoleArg> const & dipoleArg,
    bool isPrimaryParticle,
    StateType const & npc,
    vec2f const & physicsDeltaPos,
    NpcParticles & particles,
    Mesh const & mesh,
    LabParameters const & labParameters) const
{
    LogMessage("    ConstrainedNonInertial: physicsDeltaPos=", physicsDeltaPos);

    // TODOHERE
    (void)particle;
    (void)dipoleArg;
    (void)isPrimaryParticle;
    (void)npc;
    (void)particles;
    (void)mesh;
    (void)labParameters;

    /*
                            //
                            // We're moving against the floor, hence we are in a non-inertial frame...
                            // ...take friction into account and flaten trajectory
                            //

                            LogMessage("    Particle is on floor edge ", edgeOrdinal, ", moving against it");

                            //
                            // Calculate magnitude of flattened trajectory - i.e. component of trajectory
                            // along (tangent) edge, appi.e. arent (integrated) force tangential to the floor;
                            // positive when in the direction of the edge
                            //

                            float trajectoryT = trajectory.dot(edgeDir);

                            //
                            // Update tangential component of trajectory with friction
                            //

                            {
                                // Normal trajectory: apparent (integrated) force against the floor;
                                // positive when against the floor

                                float const trajectoryN = trajectory.dot(edgeNormal);

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

                                LogMessage("      friction: trajectoryN=", trajectoryN, " relVel=", particle.ConstrainedState->MeshRelativeVelocity, " trajectoryT=", trajectoryT,
                                    " tFriction=", tFriction);

                                // Update trajectory with friction
                                trajectoryT += tFriction;
                            }

                            //
                            // Recovered flattened trajectory as a vector
                            //

                            vec2f const flattenedTrajectory = edgeDir * trajectoryT;

                            vec2f targetPosition = newTheoreticalPositionAfterMeshDisplacement + flattenedTrajectory;
                            vec2f const targetAbsoluteVelocity = (targetPosition - particlePosition) / LabParameters::SimulationTimeStepDuration;

                            // Adjust target pos with superimposed displacement - which is NOT taken into account for velocity
                            std::optional<vec2f> secondaryVoluntarySuperimposedDisplacement;
                            if (npc.HumanNpcState.has_value())
                            {
                                vec2f const walkVector = edgeDir * edgeDir.dot(vec2f(1.0f, 0.0f)) * mParticles.GetVoluntarySuperimposedDisplacement(particle.ParticleIndex);
                                if (walkVector != vec2f::zero())
                                {
                                    LogMessage("WalkVector@1: ", walkVector, " ", npc.HumanNpcState->CurrentWalkingMagnitude);
                                }

                                targetPosition += walkVector;

                                if (isPrimaryParticle)
                                {
                                    // Impart same displacement to secondary particle
                                    secondaryVoluntarySuperimposedDisplacement = walkVector;
                                }
                            }

                            //
                            // Due to numerical slack, ensure target barycentric coords are along edge
                            //

                            vec3f targetBarycentricCoords = triangles.ToBarycentricCoordinates(
                                targetPosition,
                                currentTriangleElementIndex,
                                vertices);

                            LogMessage("      targetPosition=", targetPosition, " targetBarycentricCoords before forcing=", targetBarycentricCoords);

                            // Force to be on edge
                            int const vertexOrdinal = (edgeOrdinal + 2) % 3;
                            targetBarycentricCoords[vertexOrdinal] = 0.0f;
                            targetBarycentricCoords[(vertexOrdinal + 1) % 3] = 1.0f - targetBarycentricCoords[(vertexOrdinal + 2) % 3];

                            vec2f const adjustedTargetPosition = triangles.FromBarycentricCoordinates(
                                targetBarycentricCoords,
                                currentTriangleElementIndex,
                                vertices);

                            LogMessage("      flattened traj=", flattenedTrajectory, " flattened target coords=", targetBarycentricCoords, " actual target pos=", adjustedTargetPosition);

                            return CalculatedTrajectoryTargetRetVal(
                                adjustedTargetPosition,
                                targetAbsoluteVelocity,
                                SimulationStepStateType::TrajectoryStateType::ConstrainedStateType(
                                    targetBarycentricCoords,
                                    meshDisplacement),
                                std::move(secondaryVoluntarySuperimposedDisplacement));
    */

    return 0.0f;
}

float Npcs::UpdateNpcParticle_ConstrainedInertial2(
    StateType::NpcParticleStateType & particle,
    std::optional<DipoleArg> const & dipoleArg,
    bool isPrimaryParticle,
    StateType const & npc,
    vec2f const & physicsDeltaPos,
    NpcParticles & particles,
    Mesh const & mesh,
    LabParameters const & labParameters) const
{
    LogMessage("    ConstrainedInertial: physicsDeltaPos=", physicsDeltaPos);

    // TODOHERE
    (void)particle;
    (void)dipoleArg;
    (void)isPrimaryParticle;
    (void)npc;
    (void)particles;
    (void)mesh;
    (void)labParameters;
    return 0.0f;
}