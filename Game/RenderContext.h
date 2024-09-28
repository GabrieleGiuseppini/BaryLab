/***************************************************************************************
* Original Author:		Gabriele Giuseppini
* Created:				2020-05-16
* Copyright:			Gabriele Giuseppini  (https://github.com/GabrieleGiuseppini)
***************************************************************************************/
#pragma once

#include "RenderTypes.h"
#include "ViewModel.h"

#include <GameOpenGL/GameOpenGL.h>
#include <GameOpenGL/GameOpenGLMappedBuffer.h>
#include <GameOpenGL/ShaderManager.h>
#include <GameOpenGL/TriangleQuadElementArrayVBO.h>

#include <GameCore/BoundedVector.h>
#include <GameCore/GameTypes.h>
#include <GameCore/Vectors.h>

#include <optional>
#include <vector>

namespace Render {

class RenderContext;
using ShipRenderContext = RenderContext;

class RenderContext
{
public:

    RenderContext(
        int canvasWidth,
        int canvasHeight);

    ShipRenderContext & GetShipRenderContext(ShipId)
    {
        // Placeholder
        return *this;
    }

    NpcRenderModeType GetNpcRenderMode() const
    {
        return mNpcRenderMode;
    }

    void SetNpcRenderMode(NpcRenderModeType value)
    {
        mNpcRenderMode = value;
    }

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

    float ScreenOffsetToWorldOffset(float screenOffset) const
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

    void UploadSeaLevel(float value);

    void UploadVertices(
        size_t vertexCount,
        vec2f const * vertexPositions);

    void UploadEdgesStart();

    void UploadEdge(
        vec2f const & edgeEndpointAPosition,
        vec2f const & edgeEndpointBPosition,
        rgbaColor const & edgeColor,
        float thicknessAdjustment);

    void UploadEdgesEnd();

    void UploadNpcTextureQuadsStart(size_t maxQuadCount);

    Quad [[nodiscard]] & UploadNpcTextureQuadPosition()
    {
        return mNpcTextureQuadQuadBuffer.emplace_back_ghost();
    }

    void UploadNpcTextureQuadAttributes(
        float planeId,
        TextureCoordinatesQuad const & textureCoords,
        vec3f const & overlayColor)
    {
        mNpcTextureQuadTextureCoordsBuffer.emplace_back(
            Quad{
                vec2f(textureCoords.LeftX, textureCoords.TopY),
                vec2f(textureCoords.LeftX, textureCoords.BottomY),
                vec2f(textureCoords.RightX, textureCoords.TopY),
                vec2f(textureCoords.RightX, textureCoords.BottomY)
            });

        NpcTextureQuadAttributesVertex const a{
            planeId,
            overlayColor
        };

        auto * buf = &(mNpcTextureQuadAttributesVertexBuffer.emplace_back_ghost(4));        
        buf[0] = a;
        buf[1] = a;
        buf[2] = a;
        buf[3] = a;
    }

    void UploadNpcTextureQuadsEnd();

    void UploadNpcParticlesStart();

    void UploadNpcParticle(
        PlaneId planeId,
        vec2f const & particlePosition,
        rgbaColor const & particleColor,
        float alpha,
        rgbaColor const & overlayColor);

    void UploadNpcParticlesEnd();

    void UploadNpcSpringsStart();

    void UploadNpcSpring(
        PlaneId planeId,
        vec2f const & endpointAPosition,
        vec2f const & endpointBPosition,
        rgbaColor const & springColor);

    void UploadNpcSpringsEnd();

    void UploadNpcFlame(
        PlaneId /*planeId*/,
        vec2f const & /*baseCenterPosition*/,
        vec2f const & /*flameVector*/,
        float /*flameWindRotationAngle*/,
        float /*scale*/,
        float /*flamePersonalitySeed*/)
    {
        // Nop, placeholder
    }

    void UploadRectSelection(
        vec2f const & /*centerPosition*/,
        vec2f const & /*verticalDir*/,
        float /*width*/,
        float /*height*/,
        rgbColor const & /*color*/,
        float /*elapsedSimulationTime*/)
    {
        // Nop, placeholder
    }

    void UploadParticleTrajectoriesStart();

    void UploadParticleTrajectory(
        vec2f const & startPosition,
        vec2f const & endPosition,
        rgbaColor const & color);

    void UploadParticleTrajectoriesEnd();

    void UploadTrianglesStart();

    void UploadTriangle(
        vec2f const & endpointAPosition,
        vec2f const & endpointBPosition,
        vec2f const & endpointCPosition,
        rgbaColor const & particleColor);

    void UploadTrianglesEnd();

    void UploadShipVelocity(
        vec2f const & shipVelocity,
        float highlight);

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
    // Background
    ////////////////////////////////////////////////////////////////

#pragma pack(push)

    struct BackgroundVertex
    {
        vec2f PositionNdc;

        BackgroundVertex(
            vec2f const & positionNdc)
            : PositionNdc(positionNdc)
        {}
    };

#pragma pack(pop)

    GameOpenGLVAO mBackgroundVAO;

    std::vector<BackgroundVertex> mBackgroundVertexBuffer;
    GameOpenGLVBO mBackgroundVertexVBO;

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

    GameOpenGLVAO mVertexVAO;

    BLabOpenGLMappedBuffer<VertexVertex, GL_ARRAY_BUFFER> mVertexVertexBuffer;
    GameOpenGLVBO mVertexVertexVBO;

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

    GameOpenGLVAO mEdgeVAO;

    std::vector<EdgeVertex> mEdgeVertexBuffer;
    GameOpenGLVBO mEdgeVertexVBO;

    ////////////////////////////////////////////////////////////////
    // NPC Quads
    ////////////////////////////////////////////////////////////////

    struct NpcTextureQuadAttributesVertex
    {
        float PlaneId;
        vec3f OverlayColor;

        NpcTextureQuadAttributesVertex(
            float planeId,
            vec3f overlayColor)
            : PlaneId(planeId)
            , OverlayColor(overlayColor)
        {}
    };

    GameOpenGLVAO mNpcQuadVAO;

    BoundedVector<Quad> mNpcTextureQuadQuadBuffer; // 4 vertices
    GameOpenGLVBO mNpcTextureQuadQuadVBO;

    BoundedVector<Quad> mNpcTextureQuadTextureCoordsBuffer; // 4 vertices
    GameOpenGLVBO mNpcTextureQuadTextureCoordsVBO;

    BoundedVector<NpcTextureQuadAttributesVertex> mNpcTextureQuadAttributesVertexBuffer;
    GameOpenGLVBO mNpcTextureQuadAttributesVertexVBO;

    ////////////////////////////////////////////////////////////////
    // NPC Particles
    ////////////////////////////////////////////////////////////////

#pragma pack(push)

    struct NpcParticleVertex
    {
        vec2f Position;
        vec2f VertexSpacePosition;
        vec4f Color;
        vec4f OverlayColor;

        NpcParticleVertex(
            vec2f const & position,
            vec2f const & vertexSpacePosition,
            vec4f const & color,
            vec4f overlayColor)
            : Position(position)
            , VertexSpacePosition(vertexSpacePosition)
            , Color(color)
            , OverlayColor(overlayColor)
        {}
    };

#pragma pack(pop)

    GameOpenGLVAO mNpcParticleVAO;

    std::vector<NpcParticleVertex> mNpcParticleVertexBuffer;
    GameOpenGLVBO mNpcParticleVertexVBO;

    ////////////////////////////////////////////////////////////////
    // NPC Springs
    ////////////////////////////////////////////////////////////////

#pragma pack(push)

    struct NpcSpringVertex
    {
        vec2f Position;
        vec2f VertexSpacePosition;
        vec4f Color;

        NpcSpringVertex(
            vec2f const & position,
            vec2f const & vertexSpacePosition,
            vec4f const & color)
            : Position(position)
            , VertexSpacePosition(vertexSpacePosition)
            , Color(color)
        {}
    };

#pragma pack(pop)

    GameOpenGLVAO mNpcSpringVAO;

    std::vector<NpcSpringVertex> mNpcSpringVertexBuffer;
    GameOpenGLVBO mNpcSpringVertexVBO;

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

    GameOpenGLVAO mParticleTrajectoryVAO;

    std::vector<ParticleTrajectoryVertex> mParticleTrajectoryVertexBuffer;
    GameOpenGLVBO mParticleTrajectoryVertexVBO;

    ////////////////////////////////////////////////////////////////
    // Triangles
    ////////////////////////////////////////////////////////////////

#pragma pack(push)

    struct TriangleVertex
    {
        vec2f Position;
        vec4f Color;

        TriangleVertex(
            vec2f const & position,
            vec4f const & color)
            : Position(position)
            , Color(color)
        {}
    };

#pragma pack(pop)

    GameOpenGLVAO mTriangleVAO;

    std::vector<TriangleVertex> mTriangleVertexBuffer;
    GameOpenGLVBO mTriangleVertexVBO;

    ////////////////////////////////////////////////////////////////
    // Ship velocity
    ////////////////////////////////////////////////////////////////

#pragma pack(push)

    struct ShipVelocityVertex
    {
        vec2f PositionNdc;
        vec2f VertexSpacePosition;
        float Highlight;

        ShipVelocityVertex(
            vec2f const & positionNdc,
            vec2f const & vertexSpacePosition,
            float highlight)
            : PositionNdc(positionNdc)
            , VertexSpacePosition(vertexSpacePosition)
            , Highlight(highlight)
        {}
    };

#pragma pack(pop)

    GameOpenGLVAO mShipVelocityVAO;

    std::vector<ShipVelocityVertex> mShipVelocityVertexBuffer;
    GameOpenGLVBO mShipVelocityVertexVBO;

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

    GameOpenGLVAO mGridVAO;
    GameOpenGLVBO mGridVBO;
    bool mIsGridEnabled;

    ////////////////////////////////////////////////////////////////
    // Quad indices
    ////////////////////////////////////////////////////////////////

    std::unique_ptr<TriangleQuadElementArrayVBO> mElementIndices;

    ////////////////////////////////////////////////////////////////
    // Render parameters
    ////////////////////////////////////////////////////////////////

    NpcRenderModeType mNpcRenderMode;
};

}