/***************************************************************************************
* Original Author:      Gabriele Giuseppini
* Created:              2020-05-16
* Copyright:            Gabriele Giuseppini  (https://github.com/GabrieleGiuseppini)
***************************************************************************************/
#pragma once

#include "BLabTypes.h"
#include "Mesh.h"
#include "MeshBuilderTypes.h"
#include "Points.h"
#include "Springs.h"
#include "Triangles.h"

#include <cstdint>
#include <filesystem>
#include <memory>
#include <optional>
#include <vector>

/*
 * This class contains all the logic for building a mesh out of a Definition.
 */
class MeshBuilder
{
public:

    static std::unique_ptr<Mesh> BuildMesh(std::filesystem::path const & meshDefinitionFilepath);

private:
};
