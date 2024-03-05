/***************************************************************************************
 * Original Author:		Gabriele Giuseppini
 * Created:				2023-07-19
 * Copyright:			Gabriele Giuseppini  (https://github.com/GabrieleGiuseppini)
 ***************************************************************************************/
#pragma once

#include "EventDispatcher.h"
#include "Npcs.h"
#include "Physics.h"
#include "Vectors.h"

#include <memory>
#include <optional>

class Model
{
public:

	Model(
		std::unique_ptr<Physics::Ship> ship,
		std::unique_ptr<Npcs> npcs,
		EventDispatcher & eventDispatcher)
		: mShip(std::move(ship))
		, mNpcs(std::move(npcs))
		, mOriginTriangle()
		, mTrajectoryDestination()
		, mEventDispatcher(eventDispatcher)
	{}

	Physics::Ship const & GetShip() const
	{
		return *mShip;
	}

	Ship & GetShip()
	{
		return *mShip;
	}

	Npcs const & GetNpcs() const
	{
		return *mNpcs;
	}

	Npcs & GetNpcs()
	{
		return *mNpcs;
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

	void ResetNpcs(std::unique_ptr<Npcs> npcs)
	{
		mNpcs.reset();
		mNpcs = std::move(npcs);
	}

private:

	EventDispatcher & mEventDispatcher;

	std::unique_ptr<Physics::Ship> mShip;
	std::unique_ptr<Npcs> mNpcs;

	// Simulation state
	std::optional<ElementIndex> mOriginTriangle;
	std::optional<vec2f> mTrajectoryDestination;
};
