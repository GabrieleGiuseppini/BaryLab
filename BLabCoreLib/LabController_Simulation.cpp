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
            std::optional<vec3f> targetPositionCurrentTriangleBarycentricCoords;
            if (mCurrentParticleTrajectory.has_value()
                && mCurrentParticleTrajectory->ParticleIndex == p)
            {
                // We have a user-imposed trajectory for this particle...
                // ...use the provided trajectory

                targetPosition = mCurrentParticleTrajectory->TargetPosition;

                if (particleState.ConstrainedState.has_value())
                {
                    targetPositionCurrentTriangleBarycentricCoords = mModel->GetMesh().GetTriangles().ToBarycentricCoordinates(
                        targetPosition,
                        particleState.ConstrainedState->CurrentTriangle,
                        mModel->GetMesh().GetVertices());
                }
            }
            else
            {
                // Use physics to calculate trajectory

                auto const trajectoryTarget = CalculatePhysicsTarget(p, labParameters); 
                targetPosition = trajectoryTarget.Position;
                targetPositionCurrentTriangleBarycentricCoords = trajectoryTarget.CurrentTriangleBarycentricCoords;
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
                targetPositionCurrentTriangleBarycentricCoords);

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

                // Finalize particle
                if (finalParticleState->Velocity.has_value())
                {
                    particles.SetVelocity(
                        p, 
                        *(finalParticleState->Velocity));
                }
                else
                {
                    particles.SetVelocity(
                        p,
                        (finalParticleState->Position - particles.GetPosition(p)) / LabParameters::SimulationTimeStepDuration);
                }
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

    vec2f forces = particles.GetWorldForce(particleIndex) * labParameters.GravityAdjustment;

    //
    // Check whether we're on a floor
    //

    auto const & particleState = particles.GetState(particleIndex);
    if (particleState.ConstrainedState.has_value())
    {
        ElementIndex const currentTriangleElementIndex = particleState.ConstrainedState->CurrentTriangle;

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
                        // On floor edge
                        //

                        vec2f const edgeVector = mModel->GetMesh().GetTriangles().GetSubEdgeVector(
                            currentTriangleElementIndex,
                            edgeOrdinal,
                            vertices);

                        vec2f const edgeDir = edgeVector.normalise();

                        // New position of particle in moved triangle (i.e. sourcePos + non-inertial displacement)
                        vec2f const newTheoreticalPositionAfterMeshDisplacement = triangles.FromBarycentricCoordinates(
                            particleState.ConstrainedState->CurrentTriangleBarycentricCoords,
                            currentTriangleElementIndex,
                            vertices);

                        //
                        // Add friction
                        //
                        
                        {
                            //
                            // Calculate total apparent force against edge and along edge: sum of these: 
                            //  - world forces
                            //  - (apparent) force that generates mesh displacement
                            //      - TODO: should we consider relative velocity here?
                            //

                            vec2f const edgeNormal = edgeDir.to_perpendicular(); // Points outside of triangle (into floor)

                            // TODO: add forces first then project

                            float const worldForcesNormal = forces.dot(edgeNormal);
                            float const worldForcesTangent = forces.dot(edgeDir);

                            vec2f const meshDisplacement = particles.GetPosition(particleIndex) - newTheoreticalPositionAfterMeshDisplacement;
                            vec2f const meshApparentForce = meshDisplacement / (dt * dt) * particleMass;
                            float const meshApparentForceNormal = meshApparentForce.dot(edgeNormal);
                            float const meshApparentForceTangent = -meshApparentForce.dot(edgeDir);

                            float const fn = worldForcesNormal + meshApparentForceNormal;
                            float const ft = worldForcesTangent + meshApparentForceTangent;

                            LogMessage("  Friction: Fn=(", worldForcesNormal, " + ", meshApparentForceNormal, ")=", fn,
                                " Ft=(", worldForcesTangent, " + ", meshApparentForceTangent, ")=", ft);

                            //
                            // Calculate max static friction
                            //

                            float const fs = labParameters.StaticFriction * std::max(fn, 0.0f); // Friction exists only if incident

                            // TODOHERE

                            // TODO
                            //forces += frictionResponseForce;
                        }                        

                        //
                        // Integrate, including eventual bounce velocity from a previous impact
                        //

                        vec2f const deltaPos =
                            particles.GetVelocity(particleIndex) * dt
                            + forces / particleMass * dt * dt;

                        vec2f const targetPhysicsPosition = particles.GetPosition(particleIndex) + deltaPos;

                        //
                        // Calculate trajectory
                        //
                        // Trajectory is mesh move plus physics move, i.e., given that we are on a floor,
                        // it's the full, absolute move that would be traveled, if it were not for the floor
                        //
                        // Note: when there's no physics deltaPos, this amounts to pure mesh move in this frame
                        // (PartPos - FromBary(PartPosBary))
                        //

                        vec2f const meshDisplacement_test = 
                            particles.GetPosition(particleIndex)
                            - newTheoreticalPositionAfterMeshDisplacement;

                        vec2f const trajectory = 
                            targetPhysicsPosition
                            - newTheoreticalPositionAfterMeshDisplacement;

                        LogMessage("  startPos=", particles.GetPosition(particleIndex), " deltaPos=", deltaPos, " meshDispl=", meshDisplacement_test, " traj=", trajectory);

                        //
                        // Check whether we're moving *against* the floor
                        //

                        if (trajectory.dot(edgeVector.to_perpendicular()) > 0.0f) // Normal to edge is directed outside of triangle (i.e. towards floor)
                        {
                            //
                            // We're moving against the floor, so flatten the physical move
                            // (i.e. take component of move along floor)
                            //
                            // Note that staying on the edge is equivalent to having traveled
                            // the mesh move, so all that's left here is the *physical* move
                            //

                            vec2f const flattenedDeltaPos = edgeDir * deltaPos.dot(edgeDir);

                            // Due to numerical slack, ensure target barycentric coords are along edge

                            // TODOHERE: mistake? if deltaPos is zero, here we take target bary coords of current position which,
                            // if triangle is moving, would be fully outside of triangle. For example, if we start at (-1.5, -1.5) (b=(0.5, 0.0, 0.5))
                            // and mesh moves by (dm, 0.0), new target bary coords are (0.5+dm, -dm, 0.5), which if flattened by means
                            // of setting 0 coord to 0 makes new target bary coords as (0.5+dm, 0.0, 0.5-dm), yields a target pos of
                            // (-1.5, -1.484375), yielding v=(0.0, dm).
                            // If not a mistake, fix comment above.

                            vec3f targetBarycentricCoords = triangles.ToBarycentricCoordinates(
                                particles.GetPosition(particleIndex) + flattenedDeltaPos,
                                currentTriangleElementIndex,
                                vertices);

                            // Force to be on edge
                            int const vertexOrdinal = (edgeOrdinal + 2) % 3;
                            targetBarycentricCoords[vertexOrdinal] = 0.0f;
                            targetBarycentricCoords[(vertexOrdinal + 1) % 3] = 1.0f - targetBarycentricCoords[(vertexOrdinal + 2) % 3];

                            vec2f const targetCoords = triangles.FromBarycentricCoordinates(
                                targetBarycentricCoords,
                                currentTriangleElementIndex,
                                vertices);

                            LogMessage("  Particle is on floor edge ", edgeOrdinal, ", moving against it; flattened target coords: ", targetBarycentricCoords,
                                " actual target pos : ", targetCoords);

                            return TrajectoryTarget(
                                targetCoords,
                                targetBarycentricCoords);
                        }
                        else
                        {
                            LogMessage("  Particle is on floor edge, but not moving against it");

                            vec2f const targetPosition = particles.GetPosition(particleIndex) + deltaPos;

                            return TrajectoryTarget(
                                targetPosition,
                                triangles.ToBarycentricCoordinates(
                                    targetPosition,
                                    currentTriangleElementIndex,
                                    vertices));
                        }
                    }
                }
            }
        }

        // No floor, continue
    }

    //
    // Not on a floor edge, use pure force
    //

    LogMessage("  Particle not on floor edge; using pure force");

    //
    // Integrate simply, including eventual bounce velocity from a previous impact
    //

    vec2f const targetPosition =
        particles.GetPosition(particleIndex)
        + particles.GetVelocity(particleIndex) * dt
        + forces / particleMass * dt * dt;

    std::optional<vec3f> targetPositionCurrentTriangleBarycentricCoords;
    if (particleState.ConstrainedState.has_value())
    {
        targetPositionCurrentTriangleBarycentricCoords = triangles.ToBarycentricCoordinates(
            targetPosition,
            particleState.ConstrainedState->CurrentTriangle,
            vertices);
    }

    return TrajectoryTarget(
        targetPosition,
        targetPositionCurrentTriangleBarycentricCoords);
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
            // TODOTEST
            //(particleState.TrajectoryState->TargetPosition - particleState.TrajectoryState->SourcePosition) / LabParameters::SimulationTimeStepDuration);
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
    assert(particleState.TrajectoryState->TargetPositionCurrentTriangleBarycentricCoords.has_value());

    {
        assert(particleState.ConstrainedState->CurrentTriangleBarycentricCoords[0] >= 0.0f && particleState.ConstrainedState->CurrentTriangleBarycentricCoords[0] <= 1.0f);
        assert(particleState.ConstrainedState->CurrentTriangleBarycentricCoords[1] >= 0.0f && particleState.ConstrainedState->CurrentTriangleBarycentricCoords[1] <= 1.0f);
        assert(particleState.ConstrainedState->CurrentTriangleBarycentricCoords[2] >= 0.0f && particleState.ConstrainedState->CurrentTriangleBarycentricCoords[2] <= 1.0f);
    }

    ElementIndex const currentTriangle = particleState.ConstrainedState->CurrentTriangle;

    vec3f const & targetBarycentricCoords = *particleState.TrajectoryState->TargetPositionCurrentTriangleBarycentricCoords;

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
        // Update velocity with bounce response, using the *theoretical* (trajectory) 
        // velocity - since this one includes the mesh velocity
        //

        // Decompose theoretical particle velocity into normal and tangential
        vec2f const trajectory = particleState.TrajectoryState->TargetPosition - particleState.TrajectoryState->SourcePosition;
        vec2f const theoreticalParticleVelocity = trajectory / LabParameters::SimulationTimeStepDuration;
        vec2f const edgeDir =
            triangles.GetSubEdgeVector(currentTriangle, intersectionEdgeOrdinal, vertices)
            .normalise();
        vec2f const edgeNormal = edgeDir.to_perpendicular();
        float const theoreticalParticleVelocityAlongNormal = theoreticalParticleVelocity.dot(edgeNormal);
        vec2f const normalVelocity = edgeNormal * theoreticalParticleVelocityAlongNormal;
        vec2f const tangentialVelocity = theoreticalParticleVelocity - normalVelocity;

        // Calculate normal reponse: Vn' = -e*Vn (e = elasticity, [0.0 - 1.0])
        vec2f const normalResponse =
            -normalVelocity
            * labParameters.Elasticity;

        // Calculate tangential response: Vt' = a*Vt (a = (1.0-friction), [0.0 - 1.0])
        vec2f const tangentialResponse =
            tangentialVelocity
            * (1.0f - labParameters.KineticFriction);

        LogMessage("    traje=", trajectory, " pv = ", theoreticalParticleVelocity, " nr = ", normalResponse, " tr = ", tangentialResponse);

        vec2f const resultantResponseVelocity = normalResponse + tangentialResponse;

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
            particleState.TrajectoryState->TargetPositionCurrentTriangleBarycentricCoords = triangles.ToBarycentricCoordinates(
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