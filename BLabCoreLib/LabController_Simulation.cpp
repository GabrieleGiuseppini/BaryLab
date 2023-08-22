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
        vec3f const barycentricCoords = mModel->GetMesh().GetTriangles().ToBarycentricCoordinates(
            mModel->GetParticles().GetPosition(particleIndex),
            triangleIndex,
            mModel->GetMesh().GetVertices());

        constrainedState.emplace(
            triangleIndex,
            barycentricCoords);
    }

    mModel->GetParticles().GetState(particleIndex).ConstrainedState = constrainedState;
}

void LabController::UpdateSimulation(LabParameters const & labParameters)
{
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

            // Calculate target position
            vec2f targetPosition;
            if (mCurrentParticleTrajectory.has_value()
                && mCurrentParticleTrajectory->ParticleIndex == p)
            {
                // We have a user-imposed trajectory for this particle...
                // ...use the provided trajectory

                targetPosition = mCurrentParticleTrajectory->TargetPosition;
            }
            else
            {
                // Use physics to calculate trajectory

                targetPosition = 
                    particles.GetPosition(p) // Physics moves from current pos, before mesh move
                    + CalculatePhysicsDeltaPos(p, labParameters); 
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
                targetPosition);

            //
            // Update trajectory notification
            //

            mCurrentParticleTrajectory.emplace(p, targetPosition);
            mCurrentParticleTrajectoryNotification.reset();
        }
        else
        {
            assert(particleState.TrajectoryState.has_value());

            LogMessage("Particle ", p, ": Trajectory state update");

            auto const finalParticleState = UpdateParticleTrajectoryState(particles.GetState(p), labParameters);
            if (finalParticleState)
            {
                LogMessage("Particle ", p, " COMPLETED");

                // Finalize particle
                particles.SetPosition(p, finalParticleState->Position);
                particles.SetVelocity(p, finalParticleState->Velocity);

                // Destroy trajectory state
                particleState.TrajectoryState.reset();
                mCurrentParticleTrajectory.reset();
            }
        }
    }
}

vec2f LabController::CalculatePhysicsDeltaPos(
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

        int currentEdgeOrdinal = -1;
        if (particleState.ConstrainedState->CurrentTriangleBarycentricCoords.x == 0.0f)
        {
            currentEdgeOrdinal = 1;
        }
        else if (particleState.ConstrainedState->CurrentTriangleBarycentricCoords.y == 0.0f)
        {
            currentEdgeOrdinal = 2;
        }
        else if (particleState.ConstrainedState->CurrentTriangleBarycentricCoords.z == 0.0f)
        {
            currentEdgeOrdinal = 0;
        }

        if (currentEdgeOrdinal >= 0)
        {
            ElementIndex const currentEdgeElementIndex = triangles.GetSubEdges(currentTriangleElementIndex).EdgeIndices[currentEdgeOrdinal];

            // Calculate ghost condition - we don't consider a surface floor (and thus allow 
            // particle to ghost through it) if it's in a triangle made fully of floors
            bool const isGhost =
                edges.GetSurfaceType(triangles.GetSubEdgeAIndex(currentTriangleElementIndex)) == SurfaceType::Floor
                && edges.GetSurfaceType(triangles.GetSubEdgeBIndex(currentTriangleElementIndex)) == SurfaceType::Floor
                && edges.GetSurfaceType(triangles.GetSubEdgeCIndex(currentTriangleElementIndex)) == SurfaceType::Floor;

            if (edges.GetSurfaceType(currentEdgeElementIndex) == SurfaceType::Floor
                && !isGhost)
            {
                //
                // On floor edge
                //                

                //
                // Add friction
                //

                // TODO

                //
                // Integrate, including eventual bounce velocity from a previous impact
                //

                vec2f const deltaPos =
                    particles.GetVelocity(particleIndex) * dt
                    + forces / particleMass * dt * dt;

                //
                // If we're moving against the floor, flatten trajectory
                //

                vec2f const edgeVector = mModel->GetMesh().GetTriangles().GetSubEdgeVector(
                    currentTriangleElementIndex,
                    currentEdgeOrdinal,
                    vertices);

                if (deltaPos.dot(edgeVector.to_perpendicular()) > 0.0f) // Normal to edge is directed outside of triangle (i.e. towards floor)
                {
                    LogMessage("  Particle is on floor edge, moving against it");

                    //
                    // Flatten trajectory - i.e. take component of deltapos along floor
                    //

                    // TODOHERE: ...unless delta pos goes towards another edge that we're also on (i.e. we're in a vertex), in which case
                    // delta pos becomes zero
                    // TODO: is it really needed though? What happens if we're going towards another edge that we're also on (i.e. we're in vertex)?

                    vec2f const edgeDir = edgeVector.normalise();

                    return edgeDir * deltaPos.dot(edgeDir);
                }
                else
                {
                    LogMessage("  Particle is on floor edge, but not against it");

                    return deltaPos;
                }
            }
        }
    }    

    //
    // Not on a floor edge, use pure force
    //

    LogMessage("  Particle not on floor edge; using pure force");

    //
    // Integrate, including eventual bounce velocity from a previous impact
    //

    return
        particles.GetVelocity(particleIndex) * dt
        + forces / particleMass * dt * dt;

}

std::optional<LabController::FinalPerticleState> LabController::UpdateParticleTrajectoryState(
    Particles::StateType & particleState,
    LabParameters const & labParameters)
{
    Particles const & particles = mModel->GetParticles();
    Vertices const & vertices = mModel->GetMesh().GetVertices();
    Edges const & edges = mModel->GetMesh().GetEdges();
    Triangles const & triangles = mModel->GetMesh().GetTriangles();

    assert(particleState.TrajectoryState.has_value());
    vec2f const targetPosition = particleState.TrajectoryState->TargetPosition;
    vec2f const trajectory = targetPosition - particleState.TrajectoryState->SourcePosition;
    

    //
    // If we are at target, we're done
    //

    if (particleState.TrajectoryState->CurrentPosition == targetPosition)
    {
        // Reached destination

        LogMessage("  Reached destination");

        return FinalPerticleState(
            targetPosition,
            trajectory / LabParameters::SimulationTimeStepDuration);
    }

    //
    // If we're free, we move to target
    //

    if (!particleState.ConstrainedState)
    {
        LogMessage("  Particle is free, moving to target");

        // ...move to target position directly
        particleState.TrajectoryState->CurrentPosition = targetPosition;

        return std::nullopt; // We'll end at next iteration
    }

    //
    // We are in constrained state
    //

    LogMessage("  Particle is in constrained state");

    assert(particleState.ConstrainedState.has_value());

    ElementIndex const currentTriangle = particleState.ConstrainedState->CurrentTriangle;

    vec3f targetBarycentricCoords = triangles.ToBarycentricCoordinates(
        targetPosition,
        currentTriangle,
        vertices);

    //
    // If target is on/in triangle, we move to target
    //

    if (targetBarycentricCoords.x >= 0.0f && targetBarycentricCoords.x <= 1.0f
        && targetBarycentricCoords.y >= 0.0f && targetBarycentricCoords.y <= 1.0f
        && targetBarycentricCoords.z >= 0.0f && targetBarycentricCoords.z <= 1.0f)
    {
        LogMessage("  Target is on/in triangle, moving to target");

        // ...move to target position directly
        particleState.TrajectoryState->CurrentPosition = targetPosition;
        particleState.ConstrainedState->CurrentTriangleBarycentricCoords = targetBarycentricCoords;

        return std::nullopt; //  // We'll end at next iteration
    }

    //
    // We're inside or on triangle, and target is outside triangle;
    // if we're on edge, trajectory is along this edge
    //

    //
    // Find closest intersection in the direction of the trajectory
    //
    // Guaranteed to exist and within trajectory, because target is outside
    // of triangle and we're inside or on an edge
    //

    float constexpr TEpsilon = 0.0001f;    

    int intersectionVertexOrdinal = -1;
    float minIntersectionT = std::numeric_limits<float>::max();

    for (int vi = 0; vi < 3; ++vi)
    {
        int const edgeOrdinal = (vi + 1) % 3;

        // Only consider edges ahead of trajectory
        vec2f const edgeNormal = triangles.GetSubEdgeVector(currentTriangle, edgeOrdinal, vertices).to_perpendicular();
        if (trajectory.dot(edgeNormal) > 0.0f) // Stricly positive, hence not parallel
        {
            float const den = particleState.ConstrainedState->CurrentTriangleBarycentricCoords[vi] - targetBarycentricCoords[vi];
            float const t = IsAlmostZero(den)
                ? std::numeric_limits<float>::max() // Parallel, meets at infinity
                : particleState.ConstrainedState->CurrentTriangleBarycentricCoords[vi] / den;

            assert(t > -TEpsilon); // Some numeric slack, trajectory is here guaranteed to be pointing into this edge

            LogMessage("  t[v", vi, " e", edgeOrdinal, "] = ", t);            

            if (t < minIntersectionT)
            {
                intersectionVertexOrdinal = vi;
                minIntersectionT = t;
            }
        }
    }

    assert(intersectionVertexOrdinal >= 0); // Guaranteed to exist
    assert(minIntersectionT > -TEpsilon && minIntersectionT <= 1.0f); // Guaranteed to exist, and within trajectory

    //
    // Move to intersection
    //

    int const intersectionEdgeOrdinal = (intersectionVertexOrdinal + 1) % 3;
    ElementIndex const intersectionEdgeElementIndex = triangles.GetSubEdges(currentTriangle).EdgeIndices[intersectionEdgeOrdinal];

    LogMessage("  Moving to intersection with edge ", intersectionEdgeOrdinal);

    // Calculate intersection barycentric coordinates
    vec3f intersectionBarycentricCoords;
    intersectionBarycentricCoords[intersectionVertexOrdinal] = 0.0f;
    float const lNext = // Barycentric coord of next vertex at intersection
        particleState.ConstrainedState->CurrentTriangleBarycentricCoords[(intersectionVertexOrdinal + 1) % 3] * (1.0f - minIntersectionT)
        + targetBarycentricCoords[(intersectionVertexOrdinal + 1) % 3] * minIntersectionT;
    intersectionBarycentricCoords[(intersectionVertexOrdinal + 1) % 3] = lNext;
    intersectionBarycentricCoords[(intersectionVertexOrdinal + 2) % 3] = 1.0f - lNext;

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
        // Update velocity with bounce response
        //

        // Decompose particle velocity into normal and tangential
        vec2f const particleVelocity = trajectory / LabParameters::SimulationTimeStepDuration;
        vec2f const edgeDir =
            triangles.GetSubEdgeVector(currentTriangle, intersectionEdgeOrdinal, vertices)
            .normalise();
        vec2f const edgeNormal = edgeDir.to_perpendicular();
        float const particleVelocityAlongNormal = particleVelocity.dot(edgeNormal);
        vec2f const normalVelocity = edgeNormal * particleVelocityAlongNormal;
        vec2f const tangentialVelocity = particleVelocity - normalVelocity;

        // Calculate normal reponse: Vn' = -e*Vn (e = elasticity, [0.0 - 1.0])
        vec2f const normalResponse =
            -normalVelocity
            * labParameters.Elasticity;

        // Calculate tangential response: Vt' = a*Vt (a = (1.0-friction), [0.0 - 1.0])
        vec2f const tangentialResponse =
            tangentialVelocity
            * (1.0f - labParameters.Friction);

        vec2f const resultantVelocity = normalResponse + tangentialResponse;

        // 
        // Conclude here
        //

        return FinalPerticleState(
            intersectionPosition,
            resultantVelocity);
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

            // Calculate new barycentric coords (wrt opposite triangle)
            vec3f newBarycentricCoords; // In new triangle
            newBarycentricCoords[(oppositeTriangleEdgeOrdinal + 2) % 3] = 0.0f;
            newBarycentricCoords[oppositeTriangleEdgeOrdinal] = particleState.ConstrainedState->CurrentTriangleBarycentricCoords[(intersectionEdgeOrdinal + 1) % 3];
            newBarycentricCoords[(oppositeTriangleEdgeOrdinal + 1) % 3] = particleState.ConstrainedState->CurrentTriangleBarycentricCoords[intersectionEdgeOrdinal];

            LogMessage("  B-Coords: ", particleState.ConstrainedState->CurrentTriangleBarycentricCoords, " -> ", newBarycentricCoords);

            particleState.ConstrainedState->CurrentTriangleBarycentricCoords = newBarycentricCoords;

            return std::nullopt; // Continue
        }
    }
}