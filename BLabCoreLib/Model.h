/***************************************************************************************
 * Original Author:		Gabriele Giuseppini
 * Created:				2023-07-19
 * Copyright:			Gabriele Giuseppini  (https://github.com/GabrieleGiuseppini)
 ***************************************************************************************/
#pragma once

#include "EventDispatcher.h"
#include "Mesh.h"
#include "Particles.h"
#include "Vectors.h"

#include <memory>
#include <optional>

class Model
{
public:

	Model(		
		std::unique_ptr<Mesh> mesh,
		std::unique_ptr<Particles> particles,
		EventDispatcher & eventDispatcher)
		: mMesh(std::move(mesh))
		, mParticles(std::move(particles))
		, mOriginTriangle()
		, mTrajectoryDestination()
		, mEventDispatcher(eventDispatcher)
	{}

	Mesh const & GetMesh() const
	{
		return *mMesh;
	}

	Mesh & GetMesh()
	{
		return *mMesh;
	}

	Particles const & GetParticles() const
	{
		return *mParticles;
	}

	Particles & GetParticles()
	{
		return *mParticles;
	}

	void SetOriginTriangle(ElementIndex triangleIndex)
	{
		mOriginTriangle = triangleIndex;
	}

	void SetTrajectoryDestination(vec2f const & destinationPosition)
	{
		mTrajectoryDestination = destinationPosition;

		mEventDispatcher.OnTrajectoryToggled(true);
	}

private:

	EventDispatcher & mEventDispatcher;

	std::unique_ptr<Mesh> mMesh;
	std::unique_ptr<Particles> mParticles;
	
	// Simulation state
	std::optional<ElementIndex> mOriginTriangle;
	std::optional<vec2f> mTrajectoryDestination;
};
