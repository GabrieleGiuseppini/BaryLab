/***************************************************************************************
* Original Author:		Gabriele Giuseppini
* Created:				2020-05-16
* Copyright:			Gabriele Giuseppini  (https://github.com/GabrieleGiuseppini)
***************************************************************************************/
#pragma once

#include "BLabOpenGL.h"
#include "BLabOpenGLMappedBuffer.h"
#include "LabParameters.h"
#include "ShaderManager.h"
#include "Vectors.h"
#include "ViewModel.h"

#include <optional>
#include <vector>

class RenderContext
{
public:

    RenderContext(
        int canvasWidth,
        int canvasHeight);

    ////////////////////////////////////////////////////////////////
    // View properties
    ////////////////////////////////////////////////////////////////

    float const & GetZoom() const
    {
        return mViewModel.GetZoom();
    }

    void SetZoom(float zoom)
    {
        mViewModel.SetZoom(zoom);
        mIsViewModelDirty = true;
    }

    vec2f const & GetCameraWorldPosition() const
    {
        return mViewModel.GetCameraWorldPosition();
    }

    void SetCameraWorldPosition(vec2f const & pos)
    {
        mViewModel.SetCameraWorldPosition(pos);
        mIsViewModelDirty = true;
    }

    int GetCanvasWidth() const
    {
        return mViewModel.GetCanvasWidth();
    }

    int GetCanvasHeight() const
    {
        return mViewModel.GetCanvasHeight();
    }

    void SetCanvasSize(int width, int height)
    {
        mViewModel.SetCanvasSize(width, height);
        mIsViewModelDirty = true;
        mIsCanvasSizeDirty = true;
    }

    float GetVisibleWorldWidth() const
    {
        return mViewModel.GetVisibleWorldWidth();
    }

    float GetVisibleWorldHeight() const
    {
        return mViewModel.GetVisibleWorldHeight();
    }

    float GetVisibleWorldLeft() const
    {
        return mViewModel.GetVisibleWorldTopLeft().x;
    }

    float GetVisibleWorldRight() const
    {
        return mViewModel.GetVisibleWorldBottomRight().x;
    }

    float GetVisibleWorldTop() const
    {
        return mViewModel.GetVisibleWorldTopLeft().y;
    }

    float GetVisibleWorldBottom() const
    {
        return mViewModel.GetVisibleWorldBottomRight().y;
    }

    float CalculateZoomForWorldWidth(float worldWidth) const
    {
        return mViewModel.CalculateZoomForWorldWidth(worldWidth);
    }

    float CalculateZoomForWorldHeight(float worldHeight) const
    {
        return mViewModel.CalculateZoomForWorldHeight(worldHeight);
    }

    vec2f ScreenToWorld(vec2f const & screenCoordinates) const
    {
        return mViewModel.ScreenToWorld(screenCoordinates);
    }

    vec2f ScreenOffsetToWorldOffset(vec2f const & screenOffset) const
    {
        return mViewModel.ScreenOffsetToWorldOffset(screenOffset);
    }

    vec2f WorldToScreen(vec2f const & worldCoordinates) const
    {
        return mViewModel.WorldToScreen(worldCoordinates);
    }

    bool IsGridEnabled() const
    {
        return mIsGridEnabled;
    }

    void SetGridEnabled(bool isGridEnabled)
    {
        mIsGridEnabled = isGridEnabled;
        mIsGridDirty = true;
    }

    ////////////////////////////////////////////////////////////////
    // Rendering
    ////////////////////////////////////////////////////////////////

    void RenderStart();

    void UploadVertices(
        size_t vertexCount,
        vec2f const * vertexPositions);

    void UploadEdgesStart();

    void UploadEdge(
        vec2f const & edgeEndpointAPosition,
        vec2f const & edgeEndpointBPosition,
        rgbaColor const & edgeColor);

    void UploadEdgesEnd();

    void UploadParticlesStart();

    void UploadParticle(
        vec2f const & particlePosition,
        rgbaColor const & particleColor);

    void UploadParticlesEnd();

    void UploadParticleTrajectoriesStart();

    void UploadParticleTrajectory(
        vec2f const & startPosition,
        vec2f const & endPosition,
        rgbaColor const & color);

    void UploadParticleTrajectoriesEnd();

    void UploadSelectedTrianglesStart();

    void UploadSelectedTriangle(
        vec2f const & endpointAPosition,
        vec2f const & endpointBPosition,
        vec2f const & endpointCPosition);

    void UploadSelectedTrianglesEnd();

    void RenderEnd();

private:

    std::unique_ptr<ShaderManager> mShaderManager;

    ViewModel mViewModel;

private:

    ////////////////////////////////////////////////////////////////
    // Settings
    ////////////////////////////////////////////////////////////////

    void ProcessSettingChanges();

    bool mIsCanvasSizeDirty;
    void OnCanvasSizeUpdated();

    bool mIsViewModelDirty;
    void OnViewModelUpdated();

    bool mIsGridDirty;
    void OnGridUpdated();

private:

    ////////////////////////////////////////////////////////////////
    // Vertices
    ////////////////////////////////////////////////////////////////

#pragma pack(push)

    struct VertexVertex
    {
        vec2f Position;
        vec2f VertexSpacePosition;

        VertexVertex(
            vec2f const & position,
            vec2f const & vertexSpacePosition)
            : Position(position)
            , VertexSpacePosition(vertexSpacePosition)
        {}
    };

#pragma pack(pop)

    size_t mVertexVertexCount;

    BLabOpenGLVAO mVertexVAO;

    BLabOpenGLMappedBuffer<VertexVertex, GL_ARRAY_BUFFER> mVertexVertexBuffer;
    BLabOpenGLVBO mVertexVertexVBO;

    ////////////////////////////////////////////////////////////////
    // Edges
    ////////////////////////////////////////////////////////////////

#pragma pack(push)

    struct EdgeVertex
    {
        vec2f Position;
        vec2f VertexSpacePosition;
        vec4f Color;

        EdgeVertex(
            vec2f const & position,
            vec2f const & vertexSpacePosition,
            vec4f const & color)
            : Position(position)
            , VertexSpacePosition(vertexSpacePosition)
            , Color(color)
        {}
    };

#pragma pack(pop)

    BLabOpenGLVAO mEdgeVAO;

    std::vector<EdgeVertex> mEdgeVertexBuffer;
    BLabOpenGLVBO mEdgeVertexVBO;

    ////////////////////////////////////////////////////////////////
    // Particles
    ////////////////////////////////////////////////////////////////

#pragma pack(push)

    struct ParticleVertex
    {
        vec2f Position;
        vec2f VertexSpacePosition;
        vec4f Color;

        ParticleVertex(
            vec2f const & position,
            vec2f const & vertexSpacePosition,
            vec4f const & color)
            : Position(position)
            , VertexSpacePosition(vertexSpacePosition)
            , Color(color)
        {}
    };

#pragma pack(pop)

    BLabOpenGLVAO mParticleVAO;

    std::vector<ParticleVertex> mParticleVertexBuffer;
    BLabOpenGLVBO mParticleVertexVBO;

    ////////////////////////////////////////////////////////////////
    // Particle trajectory
    ////////////////////////////////////////////////////////////////

#pragma pack(push)

    struct ParticleTrajectoryVertex
    {
        vec2f Position;
        vec2f VertexSpacePosition;
        vec4f Color;

        ParticleTrajectoryVertex(
            vec2f const & position,
            vec2f const & vertexSpacePosition,
            vec4f const & color)
            : Position(position)
            , VertexSpacePosition(vertexSpacePosition)
            , Color(color)
        {}
    };

#pragma pack(pop)

    BLabOpenGLVAO mParticleTrajectoryVAO;

    std::vector<ParticleTrajectoryVertex> mParticleTrajectoryVertexBuffer;
    BLabOpenGLVBO mParticleTrajectoryVertexVBO;

    ////////////////////////////////////////////////////////////////
    // Selected triangles
    ////////////////////////////////////////////////////////////////

#pragma pack(push)

    struct SelectedTriangleVertex
    {
        vec2f Position;

        SelectedTriangleVertex(
            vec2f const & position)
            : Position(position)
        {}
    };

#pragma pack(pop)

    BLabOpenGLVAO mSelectedTriangleVAO;

    std::vector<SelectedTriangleVertex> mSelectedTriangleVertexBuffer;
    BLabOpenGLVBO mSelectedTriangleVertexVBO;

    ////////////////////////////////////////////////////////////////
    // Grid
    ////////////////////////////////////////////////////////////////
    
#pragma pack(push)

    struct GridVertex
    {
        vec2f positionObject; // Object space

        GridVertex() = default;

        GridVertex(
            vec2f const & _positionObject)
            : positionObject(_positionObject)
        {}
    };

#pragma pack(pop)

    BLabOpenGLVAO mGridVAO;
    BLabOpenGLVBO mGridVBO;
    bool mIsGridEnabled;
};