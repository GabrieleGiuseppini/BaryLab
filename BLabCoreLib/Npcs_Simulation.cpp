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
    if (IsAtBeginningOfSimulationStep())
    {
        LogMessage("----------------------------------");
        LogMessage("----------------------------------");
        LogMessage("----------------------------------");

        //
        // 1. Update behavioral state machines
        // 2. Calculate spring forces for the entire step
        //

        for (auto const n : *this)
        {
            auto & state = mStateBuffer[n];

            // Behavior

            if (state.Type == NpcType::Human)
            {
                assert(state.DipoleState.has_value());
                assert(state.HumanNpcState.has_value());

                UpdateHuman(
                    *state.HumanNpcState,
                    state.PrimaryParticleState,
                    state.DipoleState->SecondaryParticleState,
                    mesh);
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

                LogMessage("  springDisplacement[", n, "]: ", springDisplacement, " (Part1:", mParticles.GetPosition(primaryParticleIndex), 
                    " Part2:", mParticles.GetPosition(secondaryParticleIndex), ")");

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

        LogMessage("----------------------------------");
    }

    while (true)
    {
        auto & npcState = mStateBuffer[mSimulationStepState.CurrentNpcIndex];

        const auto [npcParticleState, dipoleArg, isPrimaryParticle] = [&]()
            {
                if (mSimulationStepState.CurrentIsPrimaryParticle)
                {
                    return std::tuple<StateType::NpcParticleStateType &, std::optional<TrajectoryTargetDipoleArg>, bool>(
                        npcState.PrimaryParticleState,
                        npcState.DipoleState.has_value()
                            ? TrajectoryTargetDipoleArg(npcState.DipoleState->SecondaryParticleState, npcState.DipoleState->DipoleProperties)
                            : std::optional<TrajectoryTargetDipoleArg>(),
                        true);
                }
                else
                {
                    assert(npcState.DipoleState.has_value());

                    return std::tuple<StateType::NpcParticleStateType &, std::optional<TrajectoryTargetDipoleArg>, bool>(
                        npcState.DipoleState->SecondaryParticleState,
                        TrajectoryTargetDipoleArg(
                            npcState.PrimaryParticleState,
                            npcState.DipoleState->DipoleProperties),
                        false);
                }
            }();

        ElementIndex const particleIndex = npcParticleState.ParticleIndex;

        LogMessage("UpdateNpcs - NpcIndex:", mSimulationStepState.CurrentNpcIndex, " IsPrimaryParticle:", mSimulationStepState.CurrentIsPrimaryParticle, 
            " ParticleIndex:", particleIndex);

        if (!mSimulationStepState.TrajectoryState)
        {
            //
            // Calculate a trajectory for this particle
            //

            LogMessage("  Trajectory calculation");

            // Calculate target position
            vec2f targetPosition;
            std::optional<SimulationStepStateType::TrajectoryStateType::ConstrainedStateType> trajectoryConstrainedState;
            if (mCurrentParticleTrajectory.has_value()
                && mCurrentParticleTrajectory->ParticleIndex == particleIndex)
            {
                // We have a user-imposed trajectory for this particle...
                // ...use the provided trajectory

                targetPosition = mCurrentParticleTrajectory->TargetPosition;

                // Calculate constrained state for the trajectory
                if (npcParticleState.ConstrainedState.has_value())
                {
                    // Target bary coords

                    vec3f const targetPositionCurrentTriangleBarycentricCoords = mesh.GetTriangles().ToBarycentricCoordinates(
                        targetPosition,
                        npcParticleState.ConstrainedState->CurrentTriangle,
                        mesh.GetVertices());

                    // Mesh displacement

                    vec2f const newTheoreticalPositionAfterMeshDisplacement = mesh.GetTriangles().FromBarycentricCoordinates(
                        npcParticleState.ConstrainedState->CurrentTriangleBarycentricCoords,
                        npcParticleState.ConstrainedState->CurrentTriangle,
                        mesh.GetVertices());

                    vec2f const meshDisplacement =
                        mParticles.GetPosition(particleIndex)
                        - newTheoreticalPositionAfterMeshDisplacement;

                    //

                    trajectoryConstrainedState.emplace(
                        targetPositionCurrentTriangleBarycentricCoords,
                        meshDisplacement);
                }
            }
            else
            {
                // Calculate a trajectory for this particle, based on physics

                auto const trajectoryTarget = CalculateTrajectoryTarget(
                    npcParticleState,
                    dipoleArg,
                    isPrimaryParticle,
                    mesh,
                    labParameters);
                
                //

                targetPosition = trajectoryTarget.Position;
                trajectoryConstrainedState = trajectoryTarget.ConstrainedStateInfo;
            }

            // Calculate source position
            vec2f sourcePosition;
            if (npcParticleState.ConstrainedState.has_value())
            {
                // Constrained state: trajectory is from current bary coords in cur triangle 
                // (i.e. in mesh's reference frame, after mesh movement) up to target position calculated by physics

                sourcePosition = mesh.GetTriangles().FromBarycentricCoordinates(
                    npcParticleState.ConstrainedState->CurrentTriangleBarycentricCoords,
                    npcParticleState.ConstrainedState->CurrentTriangle,
                    mesh.GetVertices());
            }
            else
            {
                // Free state: from current (absolute) position up to target position calculated by physics

                sourcePosition = mParticles.GetPosition(particleIndex);
            }

            mSimulationStepState.TrajectoryState.emplace(
                sourcePosition,
                targetPosition,
                trajectoryConstrainedState);

            //
            // Update trajectory notification
            //

            mCurrentParticleTrajectory.emplace(particleIndex, targetPosition);
            mCurrentParticleTrajectoryNotification.reset();

            if (mIsStepByStepMode)
            {
                // We are done
                break;
            }
        }

        assert(mSimulationStepState.TrajectoryState.has_value());

        //
        // Do a trajectory tracing step for the current particle
        //

        LogMessage("  Trajectory step");

        auto const finalParticleState = UpdateParticleTrajectoryTrace(npcParticleState, *mSimulationStepState.TrajectoryState, mesh, labParameters);
        if (finalParticleState)
        {
            LogMessage("  Particle ", particleIndex, " COMPLETED");

            //
            // Finalize particle's velocities
            //

            vec2f absoluteVelocity;
            if (finalParticleState->Velocity.has_value())
            {
                LogMessage("  FinalV: delivered: ", *(finalParticleState->Velocity));

                absoluteVelocity = *(finalParticleState->Velocity);
            }
            else
            {
                LogMessage("  FinalV: calculated: deltaV: ", (finalParticleState->Position - mParticles.GetPosition(particleIndex)) / LabParameters::SimulationTimeStepDuration);

                absoluteVelocity = (finalParticleState->Position - mParticles.GetPosition(particleIndex)) / LabParameters::SimulationTimeStepDuration;
            }

            mParticles.SetVelocity(
                particleIndex,
                absoluteVelocity);

            assert(mSimulationStepState.TrajectoryState.has_value());

            if (mSimulationStepState.TrajectoryState->ConstrainedState.has_value())
            {
                // We're constrained...

                // ...update mesh-relative velocity as well, for next iteration (friction)

                vec2f const meshRelativeVelocity =
                    absoluteVelocity
                    + mSimulationStepState.TrajectoryState->ConstrainedState->MeshDisplacement / LabParameters::SimulationTimeStepDuration;

                LogMessage("  FinalMRV: meshRelativeVelocity=", meshRelativeVelocity);

                assert(npcParticleState.ConstrainedState.has_value());
                npcParticleState.ConstrainedState->MeshRelativeVelocity = meshRelativeVelocity;

                // Make sure particle's regime is consistent
                mStateBuffer[mSimulationStepState.CurrentNpcIndex].Regime = StateType::RegimeType::Constrained;
            }
            else
            {
                // We're free

                assert(!npcParticleState.ConstrainedState.has_value());

                // Make sure particle's regime is consistent
                mStateBuffer[mSimulationStepState.CurrentNpcIndex].Regime = (dipoleArg.has_value() && dipoleArg->SecondaryParticle.ConstrainedState.has_value())
                    ? StateType::RegimeType::Constrained
                    : StateType::RegimeType::Free;
            }

            //
            // Finalize particle's position
            //

            LogMessage("  FinalP: calculated: ", finalParticleState->Position);

            mParticles.SetPosition(
                particleIndex, 
                finalParticleState->Position);

            //
            // Advance state
            //

            // Trajectories

            mSimulationStepState.TrajectoryState.reset();

            if (mCurrentParticleTrajectory.has_value() && mCurrentParticleTrajectory->ParticleIndex == npcParticleState.ParticleIndex)
            {
                mCurrentParticleTrajectory.reset();
            }

            // Particle

            if (mSimulationStepState.CurrentIsPrimaryParticle
                && dipoleArg.has_value())
            {
                // Move on to next particle
                mSimulationStepState.CurrentIsPrimaryParticle = false;
            }
            else
            {
                // Move on to next NPC
                mSimulationStepState.CurrentIsPrimaryParticle = true;
                ++mSimulationStepState.CurrentNpcIndex;
                if (mSimulationStepState.CurrentNpcIndex == static_cast<ElementIndex>(mStateBuffer.size()))
                {
                    // Simulation step complete
                    mSimulationStepState.CurrentNpcIndex = 0;

                    assert(mSimulationStepState.IsInitial());

                    // We're done
                    break;
                }
            }
        }
        else
        {
            // This particle requires more simulation
        }

        if (mIsStepByStepMode)
        {
            // We are done
            break;
        }
    }
}

Npcs::CalculatedTrajectoryTargetRetVal Npcs::CalculateTrajectoryTarget(
    StateType::NpcParticleStateType & particle,
    std::optional<TrajectoryTargetDipoleArg> const & dipole,
    bool isPrimaryParticle,
    Mesh const & mesh,
    LabParameters const & labParameters) const
{
    // TODO: we'll use them for forced spring displacement (first) and for ghost (second)
    (void)dipole;
    (void)isPrimaryParticle;

    float const dt = LabParameters::SimulationTimeStepDuration;
    float const particleMass = LabParameters::ParticleMass * labParameters.MassAdjustment;

    Vertices const & vertices = mesh.GetVertices();
    Edges const & edges = mesh.GetEdges();
    Triangles const & triangles = mesh.GetTriangles();

    vec2f const particlePosition = mParticles.GetPosition(particle.ParticleIndex);

    //
    // Integrate physical forces, including eventual bounce velocity from a previous impact
    //

    vec2f const physicalForces =
        mParticles.GetWorldForce(particle.ParticleIndex)
        + LabParameters::Gravity * labParameters.GravityAdjustment * mGravityGate * labParameters.ParticleMass * labParameters.MassAdjustment
        + mParticles.GetSpringForces(particle.ParticleIndex);

    vec2f const physicsDeltaPos =
        mParticles.GetVelocity(particle.ParticleIndex) * dt
        + physicalForces / particleMass * dt * dt;

    LogMessage("      springForces=", mParticles.GetSpringForces(particle.ParticleIndex), " physicsDeltaPos=", physicsDeltaPos);

    //
    // Check whether we're on a floor
    //

    // The apparent (non-inertial) displacement
    vec2f meshDisplacement;

    if (particle.ConstrainedState.has_value())
    {
        ElementIndex const currentTriangleElementIndex = particle.ConstrainedState->CurrentTriangle;

        // New position of particle in moved triangle (i.e. sourcePos + non-inertial displacement)
        // - does not include physics delta pos
        vec2f const newTheoreticalPositionAfterMeshDisplacement = triangles.FromBarycentricCoordinates(
            particle.ConstrainedState->CurrentTriangleBarycentricCoords,
            currentTriangleElementIndex,
            vertices);

        // Calculate mesh displacement
        meshDisplacement = particlePosition - newTheoreticalPositionAfterMeshDisplacement;

        // Calculate ghost condition - we don't consider a surface floor (and thus allow 
        // particle to ghost through it) if it's in a triangle made fully of floors
        bool const isGhost =
            edges.GetSurfaceType(triangles.GetSubEdgeAIndex(currentTriangleElementIndex)) == SurfaceType::Floor
            && edges.GetSurfaceType(triangles.GetSubEdgeBIndex(currentTriangleElementIndex)) == SurfaceType::Floor
            && edges.GetSurfaceType(triangles.GetSubEdgeCIndex(currentTriangleElementIndex)) == SurfaceType::Floor;

        // TODO: other ghost condition
        // TODO: encapsulate in helper method as it's used twice

        if (!isGhost)
        {
            // Check all edges, stop at first one that is floor
            for (int vi = 0; vi < 3; ++vi)
            {
                if (particle.ConstrainedState->CurrentTriangleBarycentricCoords[vi] == 0.0f)
                {
                    int const edgeOrdinal = (vi + 1) % 3;
                    ElementIndex const currentEdgeElementIndex = triangles.GetSubEdges(currentTriangleElementIndex).EdgeIndices[edgeOrdinal];

                    if (edges.GetSurfaceType(currentEdgeElementIndex) == SurfaceType::Floor)
                    {
                        //
                        // On floor edge - so potentially in a non-inertial frame
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
                        vec2f const trajectory = particlePosition + physicsDeltaPos - newTheoreticalPositionAfterMeshDisplacement;

                        LogMessage("      startPos=", particlePosition, " physicsDeltaPos=", physicsDeltaPos,
                            " meshDispl=", meshDisplacement, " traj=", trajectory);

                        //
                        // Check whether we're moving *against* the floor
                        //

                        vec2f const edgeVector = mesh.GetTriangles().GetSubEdgeVector(
                            currentTriangleElementIndex,
                            edgeOrdinal,
                            vertices);

                        vec2f const edgeDir = edgeVector.normalise();
                        vec2f const edgeNormal = edgeDir.to_perpendicular(); // Points outside of triangle (i.e. towards floor)

                        if (trajectory.dot(edgeNormal) >= 0.0f) // If 0, no normal force - hence no friction; however we want to take this codepath anyways for consistency
                        {
                            //
                            // We're moving against the floor, hence we are in a non-inertial frame
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

                            //
                            // Due to numerical slack, ensure target barycentric coords are along edge
                            //

                            vec2f const targetPos = newTheoreticalPositionAfterMeshDisplacement + flattenedTrajectory;
                            vec3f targetBarycentricCoords = triangles.ToBarycentricCoordinates(                                
                                targetPos,
                                currentTriangleElementIndex,
                                vertices);

                            LogMessage("      targetPos=", targetPos, " targetBarycentricCoords before forcing=", targetBarycentricCoords);

                            // Force to be on edge
                            int const vertexOrdinal = (edgeOrdinal + 2) % 3;
                            targetBarycentricCoords[vertexOrdinal] = 0.0f;
                            targetBarycentricCoords[(vertexOrdinal + 1) % 3] = 1.0f - targetBarycentricCoords[(vertexOrdinal + 2) % 3];

                            vec2f const targetCoords = triangles.FromBarycentricCoordinates(
                                targetBarycentricCoords,
                                currentTriangleElementIndex,
                                vertices);

                            LogMessage("      flattened traj=", flattenedTrajectory, " flattened target coords=", targetBarycentricCoords," actual target pos=", targetCoords);

                            return CalculatedTrajectoryTargetRetVal(
                                targetCoords,
                                SimulationStepStateType::TrajectoryStateType::ConstrainedStateType(
                                    targetBarycentricCoords,
                                    meshDisplacement));
                        }
                        else
                        {
                            LogMessage("    Particle is on floor edge, but not moving against it");

                            vec2f const targetPosition = mParticles.GetPosition(particle.ParticleIndex) + physicsDeltaPos;

                            return CalculatedTrajectoryTargetRetVal(
                                targetPosition,
                                SimulationStepStateType::TrajectoryStateType::ConstrainedStateType(
                                    triangles.ToBarycentricCoordinates(
                                        targetPosition,
                                        currentTriangleElementIndex,
                                        vertices),
                                    meshDisplacement));
                        }
                    }
                }
            }
        }

        // No floor, continue
    }
    else
    {
        // Inertial frame - no mesh displacement that we care about
        meshDisplacement = vec2f::zero();
    }

    //
    // Not on a floor edge - fully inertial frame
    //
    // Just use pure physical forces
    //

    LogMessage("    Particle not in constrained state or not on floor edge; using pure physical forces");

    vec2f const targetPosition = particlePosition + physicsDeltaPos;
    
    std::optional<SimulationStepStateType::TrajectoryStateType::ConstrainedStateType> constrainedState;
    if (particle.ConstrainedState.has_value())
    {
        constrainedState.emplace(
            triangles.ToBarycentricCoordinates(
                targetPosition,
                particle.ConstrainedState->CurrentTriangle,
                vertices),
            meshDisplacement);
    }

    return CalculatedTrajectoryTargetRetVal(
        targetPosition,
        constrainedState);
}

std::optional<Npcs::FinalParticleState> Npcs::UpdateParticleTrajectoryTrace(
    Npcs::StateType::NpcParticleStateType & particleState,
    SimulationStepStateType::TrajectoryStateType & trajectoryState,
    Mesh const & mesh,
    LabParameters const & labParameters)
{
    Vertices const & vertices = mesh.GetVertices();
    Edges const & edges = mesh.GetEdges();
    Triangles const & triangles = mesh.GetTriangles();

    //
    // If we are at target, we're done
    //

    if (trajectoryState.CurrentPosition == trajectoryState.TargetPosition)
    {
        // Reached destination

        LogMessage("    Reached destination");

        return FinalParticleState(
            trajectoryState.TargetPosition,
            std::nullopt);
    }

    //
    // If we're free, we move to target
    //

    if (!particleState.ConstrainedState)
    {
        LogMessage("    Particle is free, moving to target");

        // ...move to target position directly
        trajectoryState.CurrentPosition = trajectoryState.TargetPosition;

        return std::nullopt; // We'll end at next iteration
    }

    //
    // We are in constrained state
    //

    LogMessage("    Particle is in constrained state");

    assert(particleState.ConstrainedState.has_value());
    assert(trajectoryState.ConstrainedState.has_value());

    {
        assert(particleState.ConstrainedState->CurrentTriangleBarycentricCoords[0] >= 0.0f && particleState.ConstrainedState->CurrentTriangleBarycentricCoords[0] <= 1.0f);
        assert(particleState.ConstrainedState->CurrentTriangleBarycentricCoords[1] >= 0.0f && particleState.ConstrainedState->CurrentTriangleBarycentricCoords[1] <= 1.0f);
        assert(particleState.ConstrainedState->CurrentTriangleBarycentricCoords[2] >= 0.0f && particleState.ConstrainedState->CurrentTriangleBarycentricCoords[2] <= 1.0f);
    }

    ElementIndex const currentTriangle = particleState.ConstrainedState->CurrentTriangle;

    vec3f const & targetBarycentricCoords = trajectoryState.ConstrainedState->TargetPositionCurrentTriangleBarycentricCoords;

    //
    // If target is on/in triangle, we move to target
    //

    if (targetBarycentricCoords.x >= 0.0f && targetBarycentricCoords.x <= 1.0f
        && targetBarycentricCoords.y >= 0.0f && targetBarycentricCoords.y <= 1.0f
        && targetBarycentricCoords.z >= 0.0f && targetBarycentricCoords.z <= 1.0f)
    {
        LogMessage("    Target is on/in triangle (", targetBarycentricCoords, "), moving to target");

        // ...move to target position directly
        trajectoryState.CurrentPosition = trajectoryState.TargetPosition;
        particleState.ConstrainedState->CurrentTriangleBarycentricCoords = targetBarycentricCoords;

        return std::nullopt; //  // We'll end at next iteration
    }

    //
    // We're inside or on triangle, and target is outside triangle;
    // if we're on edge, trajectory is along this edge
    //

    LogMessage("    Target is outside triangle ", targetBarycentricCoords, "");

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
        if (targetBarycentricCoords[vi] < 0.0f)
        {
            float const den = particleState.ConstrainedState->CurrentTriangleBarycentricCoords[vi] - targetBarycentricCoords[vi];
            // TODOTEST
            ////float const t = IsAlmostZero(den)
            ////    ? std::numeric_limits<float>::max() // Parallel, meets at infinity
            ////    : particleState.ConstrainedState->CurrentTriangleBarycentricCoords[vi] / den;
            float const t = particleState.ConstrainedState->CurrentTriangleBarycentricCoords[vi] / den;

#ifdef _DEBUG
            diags[vi].emplace(den, t);
            diags[vi]->IntersectionPoint =
                particleState.ConstrainedState->CurrentTriangleBarycentricCoords
                + (targetBarycentricCoords - particleState.ConstrainedState->CurrentTriangleBarycentricCoords) * t;
#endif

            assert(t > -Epsilon<float>); // Some numeric slack, trajectory is here guaranteed to be pointing into this edge

            LogMessage("      t[v", vi, " e", edgeOrdinal, "] = ", t);            

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
    ElementIndex const intersectionEdgeElementIndex = triangles.GetSubEdges(currentTriangle).EdgeIndices[intersectionEdgeOrdinal];

    LogMessage("    Moving to intersection with edge ", intersectionEdgeOrdinal);

    // Calculate intersection barycentric coordinates
    vec3f intersectionBarycentricCoords;
    intersectionBarycentricCoords[intersectionVertexOrdinal] = 0.0f;
    float const lNext = Clamp( // Barycentric coord of next vertex at intersection; enforcing it's within triangle
        particleState.ConstrainedState->CurrentTriangleBarycentricCoords[(intersectionVertexOrdinal + 1) % 3] * (1.0f - minIntersectionT)
        + targetBarycentricCoords[(intersectionVertexOrdinal + 1) % 3] * minIntersectionT,
        0.0f,
        1.0f); 
    intersectionBarycentricCoords[(intersectionVertexOrdinal + 1) % 3] = lNext;
    intersectionBarycentricCoords[(intersectionVertexOrdinal + 2) % 3] = 1.0f - lNext;

    {
        assert(intersectionBarycentricCoords[0] >= 0.0f && intersectionBarycentricCoords[0] <= 1.0f);
        assert(intersectionBarycentricCoords[1] >= 0.0f && intersectionBarycentricCoords[1] <= 1.0f);
        assert(intersectionBarycentricCoords[2] >= 0.0f && intersectionBarycentricCoords[2] <= 1.0f);
    }

    // Move to intersection

    vec2f const intersectionPosition = triangles.FromBarycentricCoordinates(
        intersectionBarycentricCoords,
        currentTriangle,
        vertices);

    trajectoryState.CurrentPosition = intersectionPosition;
    particleState.ConstrainedState->CurrentTriangleBarycentricCoords = intersectionBarycentricCoords;

    //
    // Check if impacted with floor
    //

    // Calculate ghost condition - we don't consider a surface floor (and thus allow 
    // particle to ghost through it) if it's in a triangle made fully of floors
    bool const isGhost =
        edges.GetSurfaceType(triangles.GetSubEdgeAIndex(currentTriangle)) == SurfaceType::Floor
        && edges.GetSurfaceType(triangles.GetSubEdgeBIndex(currentTriangle)) == SurfaceType::Floor
        && edges.GetSurfaceType(triangles.GetSubEdgeCIndex(currentTriangle)) == SurfaceType::Floor;

    if (edges.GetSurfaceType(intersectionEdgeElementIndex) == SurfaceType::Floor
        && !isGhost)
    {
        //
        // Impact
        //

        LogMessage("    Impact");

        //
        // Update velocity with bounce response, using the *apparent* (trajectory) 
        // velocity - since this one includes the mesh velocity
        //

        // Decompose apparent particle velocity into normal and tangential
        vec2f const trajectory = trajectoryState.TargetPosition - trajectoryState.SourcePosition;
        vec2f const apparentParticleVelocity = trajectory / LabParameters::SimulationTimeStepDuration;
        vec2f const edgeDir =
            triangles.GetSubEdgeVector(currentTriangle, intersectionEdgeOrdinal, vertices)
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

        LogMessage("      traj=", trajectory, " apv=", apparentParticleVelocity, " nr=", normalResponse, " tr=", tangentialResponse);

        // Given that we calc the collision response to *trajectory* (i.e. apparent trajectory),
        // we need to transform it to absolute particle velocity
        vec2f const meshVelocity = trajectoryState.ConstrainedState->MeshDisplacement / LabParameters::SimulationTimeStepDuration;
        vec2f const resultantResponseVelocity = normalResponse + tangentialResponse - meshVelocity;

        LogMessage("      meshVelocity=", meshVelocity, " res=", resultantResponseVelocity);

        // 
        // Conclude here
        //

        return FinalParticleState(
            intersectionPosition,
            resultantResponseVelocity);
    }
    else
    {
        //
        // Not floor, climb over edge
        //

        LogMessage("    Climbing over non-floor edge");

        // Find opposite triangle
        ElementIndex const oppositeTriangle = edges.GetOppositeTriangle(intersectionEdgeElementIndex, currentTriangle);
        if (oppositeTriangle == NoneElementIndex)
        {
            //
            // Become free
            //

            LogMessage("    No opposite triangle found, becoming free");

            particleState.ConstrainedState.reset();
            trajectoryState.ConstrainedState.reset();

            return std::nullopt; // Will move and then stop
        }
        else
        {
            //
            // Move to edge of opposite triangle 
            //

            int const oppositeTriangleEdgeOrdinal = triangles.GetSubEdgeOrdinal(oppositeTriangle, intersectionEdgeElementIndex);

            LogMessage("    Moving to edge ", oppositeTriangleEdgeOrdinal, " of opposite triangle ", oppositeTriangle);

            particleState.ConstrainedState->CurrentTriangle = oppositeTriangle;

            // Calculate new target barycentric coords (wrt opposite triangle)
            trajectoryState.ConstrainedState->TargetPositionCurrentTriangleBarycentricCoords = triangles.ToBarycentricCoordinates(
                trajectoryState.TargetPosition,
                oppositeTriangle,
                vertices);

            // Calculate new current barycentric coords (wrt opposite triangle)
            vec3f newBarycentricCoords; // In new triangle
            newBarycentricCoords[(oppositeTriangleEdgeOrdinal + 2) % 3] = 0.0f;
            newBarycentricCoords[oppositeTriangleEdgeOrdinal] = particleState.ConstrainedState->CurrentTriangleBarycentricCoords[(intersectionEdgeOrdinal + 1) % 3];
            newBarycentricCoords[(oppositeTriangleEdgeOrdinal + 1) % 3] = particleState.ConstrainedState->CurrentTriangleBarycentricCoords[intersectionEdgeOrdinal];

            LogMessage("    B-Coords: ", particleState.ConstrainedState->CurrentTriangleBarycentricCoords, " -> ", newBarycentricCoords);

            {
                assert(newBarycentricCoords[0] >= 0.0f && newBarycentricCoords[0] <= 1.0f);
                assert(newBarycentricCoords[1] >= 0.0f && newBarycentricCoords[1] <= 1.0f);
                assert(newBarycentricCoords[2] >= 0.0f && newBarycentricCoords[2] <= 1.0f);
            }

            particleState.ConstrainedState->CurrentTriangleBarycentricCoords = newBarycentricCoords;

            return std::nullopt; // Continue
        }
    }
}

Npcs::StateType Npcs::MaterializeNpcState(
    ElementIndex npcIndex,
    Mesh const & mesh) const
{
    auto const & state = mStateBuffer[npcIndex];

    // Primary particle

    auto primaryParticleState = MaterializeParticleState(
        mParticles.GetPosition(state.PrimaryParticleState.ParticleIndex),
        state.PrimaryParticleState.ParticleIndex,
        mesh);

    // Secondary particle

    std::optional<StateType::DipoleStateType> dipoleState;

    if (state.DipoleState.has_value())
    {
        StateType::NpcParticleStateType secondaryParticleState = MaterializeParticleState(
            mParticles.GetPosition(state.DipoleState->SecondaryParticleState.ParticleIndex),
            state.DipoleState->SecondaryParticleState.ParticleIndex,
            mesh);

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
            mesh);
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

Npcs::StateType::NpcParticleStateType Npcs::MaterializeParticleState(
    vec2f const & position,
    ElementIndex particleIndex,
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

        constrainedState.emplace(
            triangleIndex,
            barycentricCoords);
    }

    return StateType::NpcParticleStateType(
        particleIndex,
        std::move(constrainedState));
}

void Npcs::ResetSimulationStepState()
{
    mSimulationStepState = SimulationStepStateType();
}
