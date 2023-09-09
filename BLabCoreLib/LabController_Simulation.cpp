/***************************************************************************************
 * Original Author:		Gabriele Giuseppini
 * Created:				2023-07-23
 * Copyright:			Gabriele Giuseppini  (https://github.com/GabrieleGiuseppini)
 ***************************************************************************************/
#include "LabController.h"

#include "BLabMath.h"
#include "Log.h"

#include <limits>

void LabController::InitializeParticleRegime(ElementIndex particleIndex)
{
    std::optional<Particles::StateType::ConstrainedStateType> constrainedState;

    ElementIndex const triangleIndex = FindTriangleContaining(mModel->GetParticles().GetPosition(particleIndex));
    if (triangleIndex != NoneElementIndex)
    {
        vec3f const barycentricCoords = mModel->GetMesh().GetTriangles().ToBarycentricCoordinatesFromWithinTriangle(
            mModel->GetParticles().GetPosition(particleIndex),
            triangleIndex,
            mModel->GetMesh().GetVertices());

        {
            assert(barycentricCoords[0] >= 0.0f && barycentricCoords[0] <= 1.0f);
            assert(barycentricCoords[1] >= 0.0f && barycentricCoords[1] <= 1.0f);
            assert(barycentricCoords[2] >= 0.0f && barycentricCoords[2] <= 1.0f);
        }

        constrainedState.emplace(
            triangleIndex,
            barycentricCoords);
    }

    mModel->GetParticles().GetState(particleIndex).ConstrainedState = constrainedState;
}

void LabController::UpdateSimulation(LabParameters const & labParameters)
{
    //
    // Update particles' state
    //

    auto & particles = mModel->GetParticles();

    LogMessage("----------------");
    LogMessage("----------------");

    for (auto const & p : particles)
    {
        LogMessage("----------------");

        auto & particleState = particles.GetState(p);

        if (!particleState.TrajectoryState.has_value())
        {
            //
            // We need a trajectory
            //

            LogMessage("Particle ", p, ": Trajectory state initialization");

            // Update mesh transformations
            UpdateMeshTransformations();

            // Calculate target position
            vec2f targetPosition;
            std::optional<Particles::StateType::TrajectoryStateType::ConstrainedStateType> trajectoryConstrainedState;
            if (mCurrentParticleTrajectory.has_value()
                && mCurrentParticleTrajectory->ParticleIndex == p)
            {
                // We have a user-imposed trajectory for this particle...
                // ...use the provided trajectory

                targetPosition = mCurrentParticleTrajectory->TargetPosition;

                if (particleState.ConstrainedState.has_value())
                {
                    vec3f const targetPositionCurrentTriangleBarycentricCoords = mModel->GetMesh().GetTriangles().ToBarycentricCoordinates(
                        targetPosition,
                        particleState.ConstrainedState->CurrentTriangle,
                        mModel->GetMesh().GetVertices());

                    vec2f const newTheoreticalPositionAfterMeshDisplacement = mModel->GetMesh().GetTriangles().FromBarycentricCoordinates(
                        particleState.ConstrainedState->CurrentTriangleBarycentricCoords,
                        particleState.ConstrainedState->CurrentTriangle,
                        mModel->GetMesh().GetVertices());

                    vec2f const meshDisplacement =
                        particles.GetPosition(p)
                        - newTheoreticalPositionAfterMeshDisplacement;

                    trajectoryConstrainedState.emplace(
                        targetPositionCurrentTriangleBarycentricCoords,
                        meshDisplacement);
                }
            }
            else
            {
                // Use physics to calculate trajectory

                auto const trajectoryTarget = CalculatePhysicsTarget(p, labParameters); 
                targetPosition = trajectoryTarget.Position;
                trajectoryConstrainedState = trajectoryTarget.ConstrainedStateInfo;
            }

            // Calculate source position
            vec2f sourcePosition;
            if (particleState.ConstrainedState.has_value())
            {
                // Constrained state: trajectory is from current bary coords in cur triangle 
                // (i.e. in mesh's reference frame) up to target position calculated by physics

                sourcePosition = mModel->GetMesh().GetTriangles().FromBarycentricCoordinates(
                    particleState.ConstrainedState->CurrentTriangleBarycentricCoords,
                    particleState.ConstrainedState->CurrentTriangle,
                    mModel->GetMesh().GetVertices());
            }
            else
            {
                // Free state: from current (absolute) position up to target position calculated by physics

                sourcePosition = particles.GetPosition(p);
            }

            particleState.TrajectoryState.emplace(
                sourcePosition,
                targetPosition,
                trajectoryConstrainedState);

            //
            // Update trajectory notification
            //

            mCurrentParticleTrajectory.emplace(p, targetPosition);
            mCurrentParticleTrajectoryNotification.reset();

            if (mRenderSimulationSteps)
            {
                continue;
            }
        }

        assert(particleState.TrajectoryState.has_value());

        LogMessage("Particle ", p, ": Trajectory state update");

        while (true)
        {
            auto const finalParticleState = UpdateParticleTrajectoryState(particles.GetState(p), labParameters);
            if (finalParticleState)
            {
                LogMessage("Particle ", p, " COMPLETED");

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
                    LogMessage("  FinalV: calculated: deltaV=", (finalParticleState->Position - particles.GetPosition(p)) / LabParameters::SimulationTimeStepDuration);

                    absoluteVelocity = (finalParticleState->Position - particles.GetPosition(p)) / LabParameters::SimulationTimeStepDuration;
                }

                particles.SetVelocity(
                    p,
                    absoluteVelocity);

                if (particleState.TrajectoryState->ConstrainedState.has_value())
                { 
                    vec2f const meshRelativeVelocity = 
                        absoluteVelocity
                        + particleState.TrajectoryState->ConstrainedState->MeshDisplacement / LabParameters::SimulationTimeStepDuration;

                    LogMessage("  FinalMRV: meshRelativeVelocity=", meshRelativeVelocity);

                    assert(particleState.ConstrainedState.has_value());
                    particleState.ConstrainedState->MeshRelativeVelocity = meshRelativeVelocity;
                }
                else
                {
                    // We're free
                    assert(!particleState.ConstrainedState.has_value());
                }                

                //
                // Finalize particle's position
                //

                particles.SetPosition(p, finalParticleState->Position);

                // Destroy trajectory state
                particleState.TrajectoryState.reset();
                mCurrentParticleTrajectory.reset();

                // We're done
                break;
            }

            if (mRenderSimulationSteps)
            {
                break;
            }
        }
    }
}

LabController::TrajectoryTarget LabController::CalculatePhysicsTarget(
    ElementIndex particleIndex,
    LabParameters const & labParameters) const
{
    float const dt = LabParameters::SimulationTimeStepDuration;
    float const particleMass = LabParameters::ParticleMass * labParameters.MassAdjustment;

    Particles const & particles = mModel->GetParticles();

    Vertices const & vertices = mModel->GetMesh().GetVertices();
    Edges const & edges = mModel->GetMesh().GetEdges();
    Triangles const & triangles = mModel->GetMesh().GetTriangles();

    //
    // Integrate physical forces, including eventual bounce velocity from a previous impact
    //

    vec2f const physicalForces = particles.GetWorldForce(particleIndex) * labParameters.GravityAdjustment;

    vec2f const physicsDeltaPos =
        particles.GetVelocity(particleIndex) * dt
        + physicalForces / particleMass * dt * dt;

    //
    // Check whether we're on a floor
    //

    // The apparent (non-inertial) displacement
    vec2f meshDisplacement;

    auto const & particleState = particles.GetState(particleIndex);
    if (particleState.ConstrainedState.has_value())
    {
        ElementIndex const currentTriangleElementIndex = particleState.ConstrainedState->CurrentTriangle;

        // New position of particle in moved triangle (i.e. sourcePos + non-inertial displacement)
        // - does not include physics delta pos
        vec2f const newTheoreticalPositionAfterMeshDisplacement = triangles.FromBarycentricCoordinates(
            particleState.ConstrainedState->CurrentTriangleBarycentricCoords,
            currentTriangleElementIndex,
            vertices);

        // Calculate mesh displacement
        meshDisplacement =
            particles.GetPosition(particleIndex)
            - newTheoreticalPositionAfterMeshDisplacement;

        // Calculate ghost condition - we don't consider a surface floor (and thus allow 
        // particle to ghost through it) if it's in a triangle made fully of floors
        bool const isGhost =
            edges.GetSurfaceType(triangles.GetSubEdgeAIndex(currentTriangleElementIndex)) == SurfaceType::Floor
            && edges.GetSurfaceType(triangles.GetSubEdgeBIndex(currentTriangleElementIndex)) == SurfaceType::Floor
            && edges.GetSurfaceType(triangles.GetSubEdgeCIndex(currentTriangleElementIndex)) == SurfaceType::Floor;

        if (!isGhost)
        {
            // Check all edges, stop at first one that is floor
            for (int vi = 0; vi < 3; ++vi)
            {
                if (particleState.ConstrainedState->CurrentTriangleBarycentricCoords[vi] == 0.0f)
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
                        vec2f trajectory = particles.GetPosition(particleIndex) 
                            + physicsDeltaPos
                            - newTheoreticalPositionAfterMeshDisplacement;

                        LogMessage("  startPos=", particles.GetPosition(particleIndex), " physicsDeltaPos=", physicsDeltaPos, " meshDispl=", meshDisplacement, " traj=", trajectory);

                        //
                        // Check whether we're moving *against* the floor
                        //

                        vec2f const edgeVector = mModel->GetMesh().GetTriangles().GetSubEdgeVector(
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

                            LogMessage("  Particle is on floor edge ", edgeOrdinal, ", moving against it");

                            //
                            // Update trajectory with friction
                            //

                            {
                                // Normal trajectory: apparent (integrated) force against the floor;
                                // positive when against the floor

                                float const tn = trajectory.dot(edgeNormal);

                                // Tangential trajectory: apparent (integrated) force tangential to the floor
                                // (flattened trajectory); positive when in the direction of
                                // the edge

                                float const tt = trajectory.dot(edgeDir);

                                //
                                // Choose between kinetic and static friction
                                //
                                // Note that we want to check actual velocity here, not
                                // *intention* to move (which is trajectory)
                                //
                                
                                float frictionCoefficient;
                                if (std::abs(particleState.ConstrainedState->MeshRelativeVelocity.dot(edgeDir)) > 0.01f) // Magic number
                                {
                                    // Kinetic friction
                                    frictionCoefficient = labParameters.KineticFriction;
                                }
                                else
                                {
                                    // Static friction
                                    frictionCoefficient = labParameters.StaticFriction;
                                }

                                // Calculate friction force magnitude (along edgeDir)
                                float tFriction = std::min(std::abs(tt), frictionCoefficient * std::max(tn, 0.0f));
                                if (tt >= 0.0f)
                                {
                                    tFriction *= -1.0f;
                                }

                                LogMessage("    friction: tn=", tn, " relVel=", particleState.ConstrainedState->MeshRelativeVelocity, " tt=", tt, 
                                    " tFriction=", tFriction);

                                // Update trajectory with friction
                                trajectory += edgeDir * tFriction;
                            }

                            //
                            // We're moving against the floor, so flatten the apparent move
                            // (i.e. take component of move along floor)
                            //

                            vec2f const flattenedTrajectory = edgeDir * trajectory.dot(edgeDir);

                            //
                            // Due to numerical slack, ensure target barycentric coords are along edge
                            //

                            vec2f const targetPos = newTheoreticalPositionAfterMeshDisplacement + flattenedTrajectory;
                            vec3f targetBarycentricCoords = triangles.ToBarycentricCoordinates(                                
                                targetPos,
                                currentTriangleElementIndex,
                                vertices);

                            LogMessage("    targetPos: ", targetPos, " targetBarycentricCoords before forcing : ", targetBarycentricCoords);

                            // Force to be on edge
                            int const vertexOrdinal = (edgeOrdinal + 2) % 3;
                            targetBarycentricCoords[vertexOrdinal] = 0.0f;
                            targetBarycentricCoords[(vertexOrdinal + 1) % 3] = 1.0f - targetBarycentricCoords[(vertexOrdinal + 2) % 3];

                            vec2f const targetCoords = triangles.FromBarycentricCoordinates(
                                targetBarycentricCoords,
                                currentTriangleElementIndex,
                                vertices);

                            LogMessage("    flattened traj: ", flattenedTrajectory, " flattened target coords : ", targetBarycentricCoords," actual target pos : ", targetCoords);

                            return TrajectoryTarget(
                                targetCoords,
                                Particles::StateType::TrajectoryStateType::ConstrainedStateType(
                                    targetBarycentricCoords,
                                    meshDisplacement));
                        }
                        else
                        {
                            LogMessage("  Particle is on floor edge, but not moving against it");

                            vec2f const targetPosition = particles.GetPosition(particleIndex) + physicsDeltaPos;

                            return TrajectoryTarget(
                                targetPosition,
                                Particles::StateType::TrajectoryStateType::ConstrainedStateType(
                                    triangles.ToBarycentricCoordinates(
                                        targetPosition,
                                        currentTriangleElementIndex,
                                        vertices),
                                    vec2f::zero()));
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

    LogMessage("  Particle not on floor edge; using pure physical forces");

    vec2f const targetPosition =
        particles.GetPosition(particleIndex)
        + physicsDeltaPos;
    
    std::optional<Particles::StateType::TrajectoryStateType::ConstrainedStateType> constrainedState;
    if (particleState.ConstrainedState.has_value())
    {
        constrainedState.emplace(
            triangles.ToBarycentricCoordinates(
                targetPosition,
                particleState.ConstrainedState->CurrentTriangle,
                vertices),
            meshDisplacement);
    }

    return TrajectoryTarget(
        targetPosition,
        constrainedState);
}

std::optional<LabController::FinalParticleState> LabController::UpdateParticleTrajectoryState(
    Particles::StateType & particleState,
    LabParameters const & labParameters)
{
    Particles const & particles = mModel->GetParticles();
    Vertices const & vertices = mModel->GetMesh().GetVertices();
    Edges const & edges = mModel->GetMesh().GetEdges();
    Triangles const & triangles = mModel->GetMesh().GetTriangles();

    assert(particleState.TrajectoryState.has_value());
    
    //
    // If we are at target, we're done
    //

    if (particleState.TrajectoryState->CurrentPosition == particleState.TrajectoryState->TargetPosition)
    {
        // Reached destination

        LogMessage("  Reached destination");

        return FinalParticleState(
            particleState.TrajectoryState->TargetPosition,
            std::nullopt);
    }

    //
    // If we're free, we move to target
    //

    if (!particleState.ConstrainedState)
    {
        LogMessage("  Particle is free, moving to target");

        // ...move to target position directly
        particleState.TrajectoryState->CurrentPosition = particleState.TrajectoryState->TargetPosition;

        return std::nullopt; // We'll end at next iteration
    }

    //
    // We are in constrained state
    //

    LogMessage("  Particle is in constrained state");

    assert(particleState.ConstrainedState.has_value());
    assert(particleState.TrajectoryState->ConstrainedState.has_value());

    {
        assert(particleState.ConstrainedState->CurrentTriangleBarycentricCoords[0] >= 0.0f && particleState.ConstrainedState->CurrentTriangleBarycentricCoords[0] <= 1.0f);
        assert(particleState.ConstrainedState->CurrentTriangleBarycentricCoords[1] >= 0.0f && particleState.ConstrainedState->CurrentTriangleBarycentricCoords[1] <= 1.0f);
        assert(particleState.ConstrainedState->CurrentTriangleBarycentricCoords[2] >= 0.0f && particleState.ConstrainedState->CurrentTriangleBarycentricCoords[2] <= 1.0f);
    }

    ElementIndex const currentTriangle = particleState.ConstrainedState->CurrentTriangle;

    vec3f const & targetBarycentricCoords = particleState.TrajectoryState->ConstrainedState->TargetPositionCurrentTriangleBarycentricCoords;

    //
    // If target is on/in triangle, we move to target
    //

    if (targetBarycentricCoords.x >= 0.0f && targetBarycentricCoords.x <= 1.0f
        && targetBarycentricCoords.y >= 0.0f && targetBarycentricCoords.y <= 1.0f
        && targetBarycentricCoords.z >= 0.0f && targetBarycentricCoords.z <= 1.0f)
    {
        LogMessage("  Target is on/in triangle (", targetBarycentricCoords, "), moving to target");

        // ...move to target position directly
        particleState.TrajectoryState->CurrentPosition = particleState.TrajectoryState->TargetPosition;
        particleState.ConstrainedState->CurrentTriangleBarycentricCoords = targetBarycentricCoords;

        return std::nullopt; //  // We'll end at next iteration
    }

    //
    // We're inside or on triangle, and target is outside triangle;
    // if we're on edge, trajectory is along this edge
    //

    LogMessage("  Target is outside triangle (", targetBarycentricCoords, ")");

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

            LogMessage("  t[v", vi, " e", edgeOrdinal, "] = ", t);            

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

    LogMessage("  Moving to intersection with edge ", intersectionEdgeOrdinal);

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

    particleState.TrajectoryState->CurrentPosition = intersectionPosition;
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

        LogMessage("  Impact");

        //
        // Update velocity with bounce response, using the *apparent* (trajectory) 
        // velocity - since this one includes the mesh velocity
        //

        // Decompose apparent particle velocity into normal and tangential
        vec2f const trajectory = particleState.TrajectoryState->TargetPosition - particleState.TrajectoryState->SourcePosition;
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

        LogMessage("    traj=", trajectory, " apv=", apparentParticleVelocity, " nr=", normalResponse, " tr=", tangentialResponse);

        // Given that we calc the collision response to *trajectory* (i.e. apparent trajectory),
        // we need to transform it to absolute particle velocity
        vec2f const meshVelocity = particleState.TrajectoryState->ConstrainedState->MeshDisplacement / LabParameters::SimulationTimeStepDuration;
        vec2f const resultantResponseVelocity = normalResponse + tangentialResponse - meshVelocity;

        LogMessage("    meshVelocity=", meshVelocity, " res=", resultantResponseVelocity);

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

        LogMessage("  Climbing over non-floor edge");

        // Find opposite triangle
        ElementIndex const oppositeTriangle = edges.GetOppositeTriangle(intersectionEdgeElementIndex, currentTriangle);
        if (oppositeTriangle == NoneElementIndex)
        {
            //
            // Become free
            //

            LogMessage("  No opposite triangle found, becoming free");

            particleState.ConstrainedState.reset();
            particleState.TrajectoryState->ConstrainedState.reset();

            return std::nullopt; // Will move and then stop
        }
        else
        {
            //
            // Move to edge of opposite triangle 
            //

            int const oppositeTriangleEdgeOrdinal = triangles.GetSubEdgeOrdinal(oppositeTriangle, intersectionEdgeElementIndex);

            LogMessage("  Moving to edge ", oppositeTriangleEdgeOrdinal, " of opposite triangle ", oppositeTriangle);

            particleState.ConstrainedState->CurrentTriangle = oppositeTriangle;

            // Calculate new target barycentric coords (wrt opposite triangle)
            particleState.TrajectoryState->ConstrainedState->TargetPositionCurrentTriangleBarycentricCoords = triangles.ToBarycentricCoordinates(
                particleState.TrajectoryState->TargetPosition,
                oppositeTriangle,
                vertices);

            // Calculate new current barycentric coords (wrt opposite triangle)
            vec3f newBarycentricCoords; // In new triangle
            newBarycentricCoords[(oppositeTriangleEdgeOrdinal + 2) % 3] = 0.0f;
            newBarycentricCoords[oppositeTriangleEdgeOrdinal] = particleState.ConstrainedState->CurrentTriangleBarycentricCoords[(intersectionEdgeOrdinal + 1) % 3];
            newBarycentricCoords[(oppositeTriangleEdgeOrdinal + 1) % 3] = particleState.ConstrainedState->CurrentTriangleBarycentricCoords[intersectionEdgeOrdinal];

            LogMessage("  B-Coords: ", particleState.ConstrainedState->CurrentTriangleBarycentricCoords, " -> ", newBarycentricCoords);

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