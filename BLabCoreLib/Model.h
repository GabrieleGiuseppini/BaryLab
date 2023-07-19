/***************************************************************************************
 * Original Author:		Gabriele Giuseppini
 * Created:				2023-07-19
 * Copyright:			Gabriele Giuseppini  (https://github.com/GabrieleGiuseppini)
 ***************************************************************************************/
#pragma once

#include "EventDispatcher.h"
#include "Mesh.h"
#include "Vectors.h"

#include <memory>
#include <optional>

class Model
{
public:

	Model(		
		std::unique_ptr<Mesh> mesh,
		EventDispatcher & eventDispatcher)
		: mMesh(std::move(mesh))
		, mSubjectParticle(vec2f::zero())
		, mOriginTriangle()
		, mTrajectoryDestination()
		, mEventDispatcher(eventDispatcher)
	{}

	Mesh const & GetMesh() const
	{
		return *mMesh;
	}

	vec2f const & GetSubjectParticlePosition() const
	{
		return mSubjectParticle.Position;
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

	struct Particle
	{
		vec2f Position;

		Particle(vec2f const & position)
			: Position(position)
		{}
	};

private:

	EventDispatcher & mEventDispatcher;

	std::unique_ptr<Mesh> mMesh;
	Particle mSubjectParticle;
	std::optional<ElementIndex> mOriginTriangle;
	std::optional<vec2f> mTrajectoryDestination;
};
