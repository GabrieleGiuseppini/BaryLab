/***************************************************************************************
 * Original Author:     Gabriele Giuseppini
 * Created:             2020-05-15
 * Copyright:           Gabriele Giuseppini  (https://github.com/GabrieleGiuseppini)
 ***************************************************************************************/
#pragma once

#include "Physics.h"

#include "GameParameters.h"
#include "MaterialDatabase.h"

#include <memory>
#include <vector>

 // Placeholder for real World

namespace Physics {

class World
{
public:

	World(
		MaterialDatabase const & materialDatabase,
		std::shared_ptr<GameEventDispatcher> gameEventHandler,
		GameParameters const & gameParameters,
		bool isGravityEnabled,
		float oceanDepth)
		: mGameEventHandler(std::move(gameEventHandler))
		, mAllShips()
		, mOceanSurface(oceanDepth)
		, mNpcs()
	{
		mNpcs = std::make_unique<Physics::Npcs>(
			*this,
			materialDatabase,
			mGameEventHandler,
			gameParameters,
			isGravityEnabled);
	}

	void AddShip(std::unique_ptr<Ship> ship)
	{
		mAllShips.emplace_back(std::move(ship));

		// Tell NPCs
		mNpcs->OnShipAdded(*mAllShips.back());
	}

	Ship const & GetShip() const
	{
		assert(mAllShips.size() > 0);
		return *mAllShips[0];
	}

	Ship & GetShip()
	{
		assert(mAllShips.size() > 0);
		return *mAllShips[0];
	}

	OceanSurface const & GetOceanSurface() const
	{
		return mOceanSurface;
	}

	OceanSurface & GetOceanSurface()
	{
		return mOceanSurface;
	}

	Physics::Npcs const & GetNpcs() const
	{
		return *mNpcs;
	}

	Physics::Npcs & GetNpcs()
	{
		return *mNpcs;
	}

	//

	void SetOriginTriangle(ElementIndex triangleIndex)
	{
		mOriginTriangle = triangleIndex;
	}

	void SetTrajectoryDestination(vec2f const & destinationPosition)
	{
		mTrajectoryDestination = destinationPosition;

		mGameEventHandler->OnTrajectoryToggled(true);
	}

	void Update(
		float currentSimulationTime,
		GameParameters const & gameParameters)
	{
		mNpcs->Update(
			currentSimulationTime,
			gameParameters);
	}

private:

	std::shared_ptr<GameEventDispatcher> mGameEventHandler;

	// Repository
	std::vector<std::unique_ptr<Ship>> mAllShips;
	OceanSurface mOceanSurface;
	std::unique_ptr<Npcs> mNpcs;

	// Simulation state
	std::optional<ElementIndex> mOriginTriangle;
	std::optional<vec2f> mTrajectoryDestination;
};

}