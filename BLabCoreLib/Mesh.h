/***************************************************************************************
 * Original Author:		Gabriele Giuseppini
 * Created:				2023-07-19
 * Copyright:			Gabriele Giuseppini  (https://github.com/GabrieleGiuseppini)
 ***************************************************************************************/
#pragma once

#include "Edges.h"
#include "Triangles.h"
#include "Vertices.h"

class Mesh
{
public:

    Mesh(
        Vertices && vertices,
        Edges && edges,
        Triangles && triangles)
        : mVertices(std::move(vertices))
        , mEdges(std::move(edges))
        , mTriangles(std::move(triangles))
    {}


    Vertices const & GetVertices() const
    {
        return mVertices;
    }

    Vertices & GetVertices()
    {
        return mVertices;
    }

    Edges const & GetEdges() const
    {
        return mEdges;
    }

    Edges & GetEdges()
    {
        return mEdges;
    }

    Triangles const & GetTriangles() const
    {
        return mTriangles;
    }

    Triangles & GetTriangles()
    {
        return mTriangles;
    }

private:

    Vertices mVertices;
    Edges mEdges;
    Triangles mTriangles;
};
