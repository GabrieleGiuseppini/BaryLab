/***************************************************************************************
* Original Author:		Gabriele Giuseppini
* Created:				2020-05-15
* Copyright:			Gabriele Giuseppini  (https://github.com/GabrieleGiuseppini)
***************************************************************************************/
#pragma once

#include <stdexcept>
#include <string>

class GameException : public std::runtime_error
{
public:

	GameException(std::string const & errorMessage)
		: std::runtime_error(errorMessage.c_str())
	{}
};