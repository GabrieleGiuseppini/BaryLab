/***************************************************************************************
* Original Author:		Gabriele Giuseppini
* Created:				2020-05-16
* Copyright:			Gabriele Giuseppini  (https://github.com/GabrieleGiuseppini)
***************************************************************************************/
#include "RenderContext.h"

#include "GameParameters.h"
#include "ResourceLocator.h"

#include <cmath>

namespace Render {

RenderContext::RenderContext(
    int canvasWidth,
    int canvasHeight)
    : mShaderManager()
    , mViewModel(1.0f, vec2f::zero(), canvasWidth, canvasHeight)
    // Settings
    , mIsCanvasSizeDirty(true)
    , mIsViewModelDirty(true)
    , mIsGridDirty(true)
    ////
    , mVertexVertexCount(0)
    , mIsGridEnabled(false)
    ////
    , mElementIndices()
    , mNpcRenderMode(NpcRenderModeType::Limbs)
{
    GLuint tmpGLuint;

    //
    // Initialize OpenGL
    //

    try
    {
        GameOpenGL::InitOpenGL();
    }
    catch (std::exception const & e)
    {
        throw std::runtime_error("Error during OpenGL initialization: " + std::string(e.what()));
    }

    ////////////////////////////////////////////////////////////////
    // Initialize shaders, VAO's, and VBOs
    ////////////////////////////////////////////////////////////////

    mShaderManager = ShaderManager::CreateInstance(ResourceLocator::GetShadersRootFolderPath());

    mElementIndices = TriangleQuadElementArrayVBO::Create();

    //
    // Background
    //

    glGenVertexArrays(1, &tmpGLuint);
    mBackgroundVAO = tmpGLuint;
    glBindVertexArray(*mBackgroundVAO);

    glGenBuffers(1, &tmpGLuint);
    mBackgroundVertexVBO = tmpGLuint;
    glBindBuffer(GL_ARRAY_BUFFER, *mBackgroundVertexVBO);

    glEnableVertexAttribArray(static_cast<GLuint>(ShaderManager::VertexAttributeType::BackgroundAttributeGroup1));
    glVertexAttribPointer(static_cast<GLuint>(ShaderManager::VertexAttributeType::BackgroundAttributeGroup1), 2, GL_FLOAT, GL_FALSE, sizeof(BackgroundVertex), (void *)0);
    static_assert(sizeof(BackgroundVertex) == 2 * sizeof(float));

    {
        mBackgroundVertexBuffer.emplace_back(
            vec2f(-1.0f, 1.0f));

        mBackgroundVertexBuffer.emplace_back(
            vec2f(-1.0f, -1.0f));

        mBackgroundVertexBuffer.emplace_back(
            vec2f(1.0f, 1.0f));

        mBackgroundVertexBuffer.emplace_back(
            vec2f(1.0f, -1.0f));

        glBindBuffer(GL_ARRAY_BUFFER, *mBackgroundVertexVBO);
        glBufferData(GL_ARRAY_BUFFER, sizeof(BackgroundVertex) * mBackgroundVertexBuffer.size(), mBackgroundVertexBuffer.data(), GL_STATIC_DRAW);
        glBindBuffer(GL_ARRAY_BUFFER, 0);
    }

    glBindVertexArray(0);

    //
    // Vertices
    //

    mVertexVertexCount = 0;

    glGenVertexArrays(1, &tmpGLuint);
    mVertexVAO = tmpGLuint;
    glBindVertexArray(*mVertexVAO);

    glGenBuffers(1, &tmpGLuint);
    mVertexVertexVBO = tmpGLuint;
    glBindBuffer(GL_ARRAY_BUFFER, *mVertexVertexVBO);

    glEnableVertexAttribArray(static_cast<GLuint>(ShaderManager::VertexAttributeType::VertexAttributeGroup1));
    glVertexAttribPointer(static_cast<GLuint>(ShaderManager::VertexAttributeType::VertexAttributeGroup1), 4, GL_FLOAT, GL_FALSE, sizeof(VertexVertex), (void *)0);
    static_assert(sizeof(VertexVertex) == 4 * sizeof(float));

    glBindVertexArray(0);

    //
    // Edges
    //

    glGenVertexArrays(1, &tmpGLuint);
    mEdgeVAO = tmpGLuint;
    glBindVertexArray(*mEdgeVAO);

    glGenBuffers(1, &tmpGLuint);
    mEdgeVertexVBO = tmpGLuint;
    glBindBuffer(GL_ARRAY_BUFFER, *mEdgeVertexVBO);

    glEnableVertexAttribArray(static_cast<GLuint>(ShaderManager::VertexAttributeType::EdgeAttributeGroup1));
    glVertexAttribPointer(static_cast<GLuint>(ShaderManager::VertexAttributeType::EdgeAttributeGroup1), 4, GL_FLOAT, GL_FALSE, sizeof(EdgeVertex), (void *)0);
    glEnableVertexAttribArray(static_cast<GLuint>(ShaderManager::VertexAttributeType::EdgeAttributeGroup2));
    glVertexAttribPointer(static_cast<GLuint>(ShaderManager::VertexAttributeType::EdgeAttributeGroup2), 4, GL_FLOAT, GL_FALSE, sizeof(EdgeVertex), (void *)(4 * sizeof(float)));
    static_assert(sizeof(EdgeVertex) == 8 * sizeof(float));

    glBindVertexArray(0);

    //
    // NPC Quads
    //

    glGenVertexArrays(1, &tmpGLuint);
    mNpcQuadVAO = tmpGLuint;
    glBindVertexArray(*mNpcQuadVAO);

    glGenBuffers(1, &tmpGLuint);
    mNpcQuadVertexVBO = tmpGLuint;
    glBindBuffer(GL_ARRAY_BUFFER, *mNpcQuadVertexVBO);

    glEnableVertexAttribArray(static_cast<GLuint>(ShaderManager::VertexAttributeType::NpcQuadAttributeGroup1));
    glVertexAttribPointer(static_cast<GLuint>(ShaderManager::VertexAttributeType::NpcQuadAttributeGroup1), 4, GL_FLOAT, GL_FALSE, sizeof(NpcQuadVertex), (void *)0);
    glEnableVertexAttribArray(static_cast<GLuint>(ShaderManager::VertexAttributeType::NpcQuadAttributeGroup2));
    glVertexAttribPointer(static_cast<GLuint>(ShaderManager::VertexAttributeType::NpcQuadAttributeGroup2), 2, GL_FLOAT, GL_FALSE, sizeof(NpcQuadVertex), (void *)(4 * sizeof(float)));
    glEnableVertexAttribArray(static_cast<GLuint>(ShaderManager::VertexAttributeType::NpcQuadAttributeGroup3));
    glVertexAttribPointer(static_cast<GLuint>(ShaderManager::VertexAttributeType::NpcQuadAttributeGroup3), 4, GL_FLOAT, GL_FALSE, sizeof(NpcQuadVertex), (void *)((4 + 2) * sizeof(float)));
    static_assert(sizeof(NpcQuadVertex) == (4 + 2 + 4) * sizeof(float));

    glBindVertexArray(0);

    //
    // NPC Particles
    //

    glGenVertexArrays(1, &tmpGLuint);
    mNpcParticleVAO = tmpGLuint;
    glBindVertexArray(*mNpcParticleVAO);

    glGenBuffers(1, &tmpGLuint);
    mNpcParticleVertexVBO = tmpGLuint;
    glBindBuffer(GL_ARRAY_BUFFER, *mNpcParticleVertexVBO);

    glEnableVertexAttribArray(static_cast<GLuint>(ShaderManager::VertexAttributeType::NpcParticleAttributeGroup1));
    glVertexAttribPointer(static_cast<GLuint>(ShaderManager::VertexAttributeType::NpcParticleAttributeGroup1), 4, GL_FLOAT, GL_FALSE, sizeof(NpcParticleVertex), (void *)0);
    glEnableVertexAttribArray(static_cast<GLuint>(ShaderManager::VertexAttributeType::NpcParticleAttributeGroup2));
    glVertexAttribPointer(static_cast<GLuint>(ShaderManager::VertexAttributeType::NpcParticleAttributeGroup2), 2, GL_FLOAT, GL_FALSE, sizeof(NpcParticleVertex), (void *)(4 * sizeof(float)));
    glEnableVertexAttribArray(static_cast<GLuint>(ShaderManager::VertexAttributeType::NpcParticleAttributeGroup3));
    glVertexAttribPointer(static_cast<GLuint>(ShaderManager::VertexAttributeType::NpcParticleAttributeGroup3), 4, GL_FLOAT, GL_FALSE, sizeof(NpcParticleVertex), (void *)((4 + 2) * sizeof(float)));
    static_assert(sizeof(NpcParticleVertex) == (4 + 4 + 4) * sizeof(float));

    glBindVertexArray(0);

    //
    // NPC Springs
    //

    glGenVertexArrays(1, &tmpGLuint);
    mNpcSpringVAO = tmpGLuint;
    glBindVertexArray(*mNpcSpringVAO);

    glGenBuffers(1, &tmpGLuint);
    mNpcSpringVertexVBO = tmpGLuint;
    glBindBuffer(GL_ARRAY_BUFFER, *mNpcSpringVertexVBO);

    glEnableVertexAttribArray(static_cast<GLuint>(ShaderManager::VertexAttributeType::NpcSpringAttributeGroup1));
    glVertexAttribPointer(static_cast<GLuint>(ShaderManager::VertexAttributeType::NpcSpringAttributeGroup1), 4, GL_FLOAT, GL_FALSE, sizeof(NpcSpringVertex), (void *)0);
    glEnableVertexAttribArray(static_cast<GLuint>(ShaderManager::VertexAttributeType::NpcSpringAttributeGroup2));
    glVertexAttribPointer(static_cast<GLuint>(ShaderManager::VertexAttributeType::NpcSpringAttributeGroup2), 4, GL_FLOAT, GL_FALSE, sizeof(NpcSpringVertex), (void *)(4 * sizeof(float)));
    static_assert(sizeof(NpcSpringVertex) == 8 * sizeof(float));

    glBindVertexArray(0);

    //
    // Particle trajectories
    //

    glGenVertexArrays(1, &tmpGLuint);
    mParticleTrajectoryVAO = tmpGLuint;
    glBindVertexArray(*mParticleTrajectoryVAO);

    glGenBuffers(1, &tmpGLuint);
    mParticleTrajectoryVertexVBO = tmpGLuint;
    glBindBuffer(GL_ARRAY_BUFFER, *mParticleTrajectoryVertexVBO);

    glEnableVertexAttribArray(static_cast<GLuint>(ShaderManager::VertexAttributeType::ParticleTrajectoryAttributeGroup1));
    glVertexAttribPointer(static_cast<GLuint>(ShaderManager::VertexAttributeType::ParticleTrajectoryAttributeGroup1), 4, GL_FLOAT, GL_FALSE, sizeof(ParticleTrajectoryVertex), (void *)0);
    glEnableVertexAttribArray(static_cast<GLuint>(ShaderManager::VertexAttributeType::ParticleTrajectoryAttributeGroup2));
    glVertexAttribPointer(static_cast<GLuint>(ShaderManager::VertexAttributeType::ParticleTrajectoryAttributeGroup2), 4, GL_FLOAT, GL_FALSE, sizeof(ParticleTrajectoryVertex), (void *)(4 * sizeof(float)));
    static_assert(sizeof(ParticleTrajectoryVertex) == (4 + 4) * sizeof(float));

    glBindVertexArray(0);

    //
    // Triangles
    //

    glGenVertexArrays(1, &tmpGLuint);
    mTriangleVAO = tmpGLuint;
    glBindVertexArray(*mTriangleVAO);

    glGenBuffers(1, &tmpGLuint);
    mTriangleVertexVBO = tmpGLuint;
    glBindBuffer(GL_ARRAY_BUFFER, *mTriangleVertexVBO);

    glEnableVertexAttribArray(static_cast<GLuint>(ShaderManager::VertexAttributeType::TriangleAttributeGroup1));
    glVertexAttribPointer(static_cast<GLuint>(ShaderManager::VertexAttributeType::TriangleAttributeGroup1), 2, GL_FLOAT, GL_FALSE, sizeof(TriangleVertex), (void *)0);
    glEnableVertexAttribArray(static_cast<GLuint>(ShaderManager::VertexAttributeType::TriangleAttributeGroup2));
    glVertexAttribPointer(static_cast<GLuint>(ShaderManager::VertexAttributeType::TriangleAttributeGroup2), 4, GL_FLOAT, GL_FALSE, sizeof(TriangleVertex), (void *)(2 * sizeof(float)));
    static_assert(sizeof(TriangleVertex) ==(2 + 4) * sizeof(float));

    glBindVertexArray(0);

    //
    // Ship velocity
    //

    glGenVertexArrays(1, &tmpGLuint);
    mShipVelocityVAO = tmpGLuint;
    glBindVertexArray(*mShipVelocityVAO);

    glGenBuffers(1, &tmpGLuint);
    mShipVelocityVertexVBO = tmpGLuint;
    glBindBuffer(GL_ARRAY_BUFFER, *mShipVelocityVertexVBO);

    glEnableVertexAttribArray(static_cast<GLuint>(ShaderManager::VertexAttributeType::ShipVelocityAttributeGroup1));
    glVertexAttribPointer(static_cast<GLuint>(ShaderManager::VertexAttributeType::ShipVelocityAttributeGroup1), 4, GL_FLOAT, GL_FALSE, sizeof(ShipVelocityVertex), (void *)0);
    glEnableVertexAttribArray(static_cast<GLuint>(ShaderManager::VertexAttributeType::ShipVelocityAttributeGroup2));
    glVertexAttribPointer(static_cast<GLuint>(ShaderManager::VertexAttributeType::ShipVelocityAttributeGroup2), 1, GL_FLOAT, GL_FALSE, sizeof(ShipVelocityVertex), (void *)(4 * sizeof(float)));
    static_assert(sizeof(ShipVelocityVertex) == (4 + 1) * sizeof(float));

    glBindVertexArray(0);

    //
    // Grid
    //

    glGenVertexArrays(1, &tmpGLuint);
    mGridVAO = tmpGLuint;
    glBindVertexArray(*mGridVAO);

    glGenBuffers(1, &tmpGLuint);
    mGridVBO = tmpGLuint;
    glBindBuffer(GL_ARRAY_BUFFER, *mGridVBO);

    glEnableVertexAttribArray(static_cast<GLuint>(ShaderManager::VertexAttributeType::GridAttributeGroup1));
    glVertexAttribPointer(static_cast<GLuint>(ShaderManager::VertexAttributeType::GridAttributeGroup1), 2, GL_FLOAT, GL_FALSE, sizeof(GridVertex), (void *)0);

    glBindVertexArray(0);

    ////////////////////////////////////////////////////////////////
    // Initialize global settings
    ////////////////////////////////////////////////////////////////

    // Enable blend for alpha transparency
    glEnable(GL_BLEND);
    glBlendFunc(GL_SRC_ALPHA, GL_ONE_MINUS_SRC_ALPHA);

    // Disable depth test
    glDisable(GL_DEPTH_TEST);

    // Set polygon mode
    glPolygonMode(GL_FRONT_AND_BACK, GL_FILL);

    ////////////////////////////////////////////////////////////////
    // Set parameters in all shaders
    ////////////////////////////////////////////////////////////////

    ProcessSettingChanges();
}

void RenderContext::RenderStart()
{
    // Process setting changes
    ProcessSettingChanges();
}

void RenderContext::UploadSeaLevel(float value)
{
    mShaderManager->ActivateProgram<ShaderManager::ProgramType::Background>();
    mShaderManager->SetProgramParameter<ShaderManager::ProgramType::Background, ShaderManager::ProgramParameterType::SeaLevel>(
        mViewModel.WorldToNdc(vec2f(0.0f, value)).y);

}

void RenderContext::UploadVertices(
    size_t vertexCount,
    vec2f const * vertexPositions)
{
    //
    // Map buffer
    //

    glBindBuffer(GL_ARRAY_BUFFER, *mVertexVertexVBO);

    // Check whether we need to re-allocate the buffers
    if (vertexCount * 6 != mVertexVertexCount)
    {
        mVertexVertexCount = vertexCount * 6;

        glBufferData(GL_ARRAY_BUFFER, mVertexVertexCount * sizeof(VertexVertex), nullptr, GL_STREAM_DRAW);
        CheckOpenGLError();
    }

    mVertexVertexBuffer.map(mVertexVertexCount);

    //
    // Upload buffer
    //

    for (size_t v = 0; v < vertexCount; ++v)
    {
        vec2f const & vertexPosition = vertexPositions[v];

        float const xLeft = vertexPosition.x - GameParameters::VertexRadius;
        float const xRight = vertexPosition.x + GameParameters::VertexRadius;
        float const yTop = vertexPosition.y + GameParameters::VertexRadius;
        float const yBottom = vertexPosition.y - GameParameters::VertexRadius;

        // Left, bottom
        mVertexVertexBuffer.emplace_back(
            vec2f(xLeft, yBottom),
            vec2f(-1.0f, -1.0f));

        // Left, top
        mVertexVertexBuffer.emplace_back(
            vec2f(xLeft, yTop),
            vec2f(-1.0f, 1.0f));

        // Right, bottom
        mVertexVertexBuffer.emplace_back(
            vec2f(xRight, yBottom),
            vec2f(1.0f, -1.0f));

        // Left, top
        mVertexVertexBuffer.emplace_back(
            vec2f(xLeft, yTop),
            vec2f(-1.0f, 1.0f));

        // Right, bottom
        mVertexVertexBuffer.emplace_back(
            vec2f(xRight, yBottom),
            vec2f(1.0f, -1.0f));

        // Right, top
        mVertexVertexBuffer.emplace_back(
            vec2f(xRight, yTop),
            vec2f(1.0f, 1.0f));
    }

    //
    // Unmap buffer
    //

    assert(mVertexVertexBuffer.size() == mVertexVertexCount);

    mVertexVertexBuffer.unmap();

    glBindBuffer(GL_ARRAY_BUFFER, 0);
}

void RenderContext::UploadEdgesStart()
{
    //
    // Prepare buffer
    //

    mEdgeVertexBuffer.clear();
}

void RenderContext::UploadEdge(
    vec2f const & edgeEndpointAPosition,
    vec2f const & edgeEndpointBPosition,
    rgbaColor const & edgeColor,
    float thicknessAdjustment)
{
    vec2f const edgeVector = edgeEndpointBPosition - edgeEndpointAPosition;
    vec2f const edgeNormal = edgeVector.to_perpendicular().normalise() * thicknessAdjustment * GameParameters::EdgeThickness / 2.0f;

    vec2f const bottomLeft = edgeEndpointAPosition - edgeNormal;
    vec2f const bottomRight = edgeEndpointAPosition + edgeNormal;
    vec2f const topLeft = edgeEndpointBPosition - edgeNormal;
    vec2f const topRight = edgeEndpointBPosition + edgeNormal;

    vec4f const color = edgeColor.toVec4f();

    // Left, bottom
    mEdgeVertexBuffer.emplace_back(
        bottomLeft,
        vec2f(-1.0f, -1.0f),
        color);

    // Left, top
    mEdgeVertexBuffer.emplace_back(
        topLeft,
        vec2f(-1.0f, 1.0f),
        color);

    // Right, bottom
    mEdgeVertexBuffer.emplace_back(
        bottomRight,
        vec2f(1.0f, -1.0f),
        color);

    // Left, top
    mEdgeVertexBuffer.emplace_back(
        topLeft,
        vec2f(-1.0f, 1.0f),
        color);

    // Right, bottom
    mEdgeVertexBuffer.emplace_back(
        bottomRight,
        vec2f(1.0f, -1.0f),
        color);

    // Right, top
    mEdgeVertexBuffer.emplace_back(
        topRight,
        vec2f(1.0f, 1.0f),
        color);
}

void RenderContext::UploadEdgesEnd()
{
    //
    // Upload buffer, if needed
    //

    if (!mEdgeVertexBuffer.empty())
    {
        glBindBuffer(GL_ARRAY_BUFFER, *mEdgeVertexVBO);
        glBufferData(GL_ARRAY_BUFFER, sizeof(EdgeVertex) * mEdgeVertexBuffer.size(), mEdgeVertexBuffer.data(), GL_STREAM_DRAW);
        glBindBuffer(GL_ARRAY_BUFFER, 0);
    }
}

void RenderContext::UploadNpcTextureQuadsStart(size_t quadCount)
{
    //
    // Prepare buffer and indices
    //

    mNpcQuadVertexBuffer.reset(quadCount * 4);

    mElementIndices->EnsureSize(quadCount);
}

void RenderContext::UploadNpcTextureQuadsEnd()
{
    //
    // Upload buffer, if needed
    //

    if (!mNpcQuadVertexBuffer.empty())
    {
        glBindBuffer(GL_ARRAY_BUFFER, *mNpcQuadVertexVBO);
        glBufferData(GL_ARRAY_BUFFER, sizeof(NpcQuadVertex) * mNpcQuadVertexBuffer.size(), mNpcQuadVertexBuffer.data(), GL_STREAM_DRAW);
        glBindBuffer(GL_ARRAY_BUFFER, 0);
    }
}

void RenderContext::UploadNpcParticlesStart()
{
    //
    // Prepare buffer
    //

    mNpcParticleVertexBuffer.clear();
}

void RenderContext::UploadNpcParticle(
    PlaneId /*planeId*/,
    vec2f const & particlePosition,
    rgbaColor const & particleColor,
    float alpha,
    NpcHighlightType highlight)
{
    float const xLeft = particlePosition.x - GameParameters::ParticleRadius;
    float const xRight = particlePosition.x + GameParameters::ParticleRadius;
    float const yTop = particlePosition.y + GameParameters::ParticleRadius;
    float const yBottom = particlePosition.y - GameParameters::ParticleRadius;

    vec4f color = particleColor.toVec4f();
    color.w *= alpha;

    vec4f const overlayColor = NpcHighlightToOverlayColor(highlight);

    // Left, bottom
    mNpcParticleVertexBuffer.emplace_back(
        vec2f(xLeft, yBottom),
        vec2f(-1.0f, -1.0f),
        color,
        overlayColor);

    // Left, top
    mNpcParticleVertexBuffer.emplace_back(
        vec2f(xLeft, yTop),
        vec2f(-1.0f, 1.0f),
        color,
        overlayColor);

    // Right, bottom
    mNpcParticleVertexBuffer.emplace_back(
        vec2f(xRight, yBottom),
        vec2f(1.0f, -1.0f),
        color,
        overlayColor);

    // Left, top
    mNpcParticleVertexBuffer.emplace_back(
        vec2f(xLeft, yTop),
        vec2f(-1.0f, 1.0f),
        color,
        overlayColor);

    // Right, bottom
    mNpcParticleVertexBuffer.emplace_back(
        vec2f(xRight, yBottom),
        vec2f(1.0f, -1.0f),
        color,
        overlayColor);

    // Right, top
    mNpcParticleVertexBuffer.emplace_back(
        vec2f(xRight, yTop),
        vec2f(1.0f, 1.0f),
        color,
        overlayColor);
}

void RenderContext::UploadNpcParticlesEnd()
{
    //
    // Upload buffer, if needed
    //

    if (!mNpcParticleVertexBuffer.empty())
    {
        glBindBuffer(GL_ARRAY_BUFFER, *mNpcParticleVertexVBO);
        glBufferData(GL_ARRAY_BUFFER, sizeof(NpcParticleVertex) * mNpcParticleVertexBuffer.size(), mNpcParticleVertexBuffer.data(), GL_STREAM_DRAW);
        glBindBuffer(GL_ARRAY_BUFFER, 0);
    }
}


void RenderContext::UploadNpcSpringsStart()
{
    //
    // Prepare buffer
    //

    mNpcSpringVertexBuffer.clear();
}

void RenderContext::UploadNpcSpring(
    PlaneId /*planeId*/,
    vec2f const & endpointAPosition,
    vec2f const & endpointBPosition,
    rgbaColor const & springColor)
{
    vec2f const springVector = endpointBPosition - endpointAPosition;
    vec2f const springNormal = springVector.to_perpendicular().normalise() * GameParameters::SpringThickness / 2.0f;

    vec2f const bottomLeft = endpointAPosition - springNormal;
    vec2f const bottomRight = endpointAPosition + springNormal;
    vec2f const topLeft = endpointBPosition - springNormal;
    vec2f const topRight = endpointBPosition + springNormal;

    vec4f const color = springColor.toVec4f();

    // Left, bottom
    mNpcSpringVertexBuffer.emplace_back(
        bottomLeft,
        vec2f(-1.0f, -1.0f),
        color);

    // Left, top
    mNpcSpringVertexBuffer.emplace_back(
        topLeft,
        vec2f(-1.0f, 1.0f),
        color);

    // Right, bottom
    mNpcSpringVertexBuffer.emplace_back(
        bottomRight,
        vec2f(1.0f, -1.0f),
        color);

    // Left, top
    mNpcSpringVertexBuffer.emplace_back(
        topLeft,
        vec2f(-1.0f, 1.0f),
        color);

    // Right, bottom
    mNpcSpringVertexBuffer.emplace_back(
        bottomRight,
        vec2f(1.0f, -1.0f),
        color);

    // Right, top
    mNpcSpringVertexBuffer.emplace_back(
        topRight,
        vec2f(1.0f, 1.0f),
        color);
}

void RenderContext::UploadNpcSpringsEnd()
{
    //
    // Upload buffer, if needed
    //

    if (!mNpcSpringVertexBuffer.empty())
    {
        glBindBuffer(GL_ARRAY_BUFFER, *mNpcSpringVertexVBO);
        glBufferData(GL_ARRAY_BUFFER, sizeof(NpcSpringVertex) * mNpcSpringVertexBuffer.size(), mNpcSpringVertexBuffer.data(), GL_STREAM_DRAW);
        glBindBuffer(GL_ARRAY_BUFFER, 0);
    }
}

void RenderContext::UploadParticleTrajectoriesStart()
{
    //
    // Prepare buffer
    //

    mParticleTrajectoryVertexBuffer.clear();
}

void RenderContext::UploadParticleTrajectory(
    vec2f const & startPosition,
    vec2f const & endPosition,
    rgbaColor const & color)
{
    vec2f const trajectoryVector = endPosition - startPosition;
    vec2f const trajectoryNormal = trajectoryVector.to_perpendicular().normalise() * GameParameters::ParticleTrajectoryThickness / 2.0f;

    vec2f const bottomLeft = startPosition - trajectoryNormal;
    vec2f const bottomRight = startPosition + trajectoryNormal;
    vec2f const topLeft = endPosition - trajectoryNormal;
    vec2f const topRight = endPosition + trajectoryNormal;

    vec4f const colorf = color.toVec4f();

    // Left, bottom
    mParticleTrajectoryVertexBuffer.emplace_back(
        bottomLeft,
        vec2f(-1.0f, -1.0f),
        colorf);

    // Left, top
    mParticleTrajectoryVertexBuffer.emplace_back(
        topLeft,
        vec2f(-1.0f, 1.0f),
        colorf);

    // Right, bottom
    mParticleTrajectoryVertexBuffer.emplace_back(
        bottomRight,
        vec2f(1.0f, -1.0f),
        colorf);

    // Left, top
    mParticleTrajectoryVertexBuffer.emplace_back(
        topLeft,
        vec2f(-1.0f, 1.0f),
        colorf);

    // Right, bottom
    mParticleTrajectoryVertexBuffer.emplace_back(
        bottomRight,
        vec2f(1.0f, -1.0f),
        colorf);

    // Right, top
    mParticleTrajectoryVertexBuffer.emplace_back(
        topRight,
        vec2f(1.0f, 1.0f),
        colorf);
}

void RenderContext::UploadParticleTrajectoriesEnd()
{
    //
    // Upload buffer, if needed
    //

    if (!mParticleTrajectoryVertexBuffer.empty())
    {
        glBindBuffer(GL_ARRAY_BUFFER, *mParticleTrajectoryVertexVBO);
        glBufferData(GL_ARRAY_BUFFER, sizeof(ParticleTrajectoryVertex) * mParticleTrajectoryVertexBuffer.size(), mParticleTrajectoryVertexBuffer.data(), GL_STREAM_DRAW);
        glBindBuffer(GL_ARRAY_BUFFER, 0);
    }
}

void RenderContext::UploadTrianglesStart()
{
    //
    // Prepare buffer
    //

    mTriangleVertexBuffer.clear();
}

void RenderContext::UploadTriangle(
    vec2f const & endpointAPosition,
    vec2f const & endpointBPosition,
    vec2f const & endpointCPosition,
    rgbaColor const & color)
{
    auto const colorf = color.toVec4f();

    mTriangleVertexBuffer.emplace_back(endpointAPosition, colorf);
    mTriangleVertexBuffer.emplace_back(endpointBPosition, colorf);
    mTriangleVertexBuffer.emplace_back(endpointCPosition, colorf);
}

void RenderContext::UploadTrianglesEnd()
{
    //
    // Upload buffer, if needed
    //

    if (!mTriangleVertexBuffer.empty())
    {
        glBindBuffer(GL_ARRAY_BUFFER, *mTriangleVertexVBO);
        glBufferData(GL_ARRAY_BUFFER, sizeof(TriangleVertex) * mTriangleVertexBuffer.size(), mTriangleVertexBuffer.data(), GL_STREAM_DRAW);
        glBindBuffer(GL_ARRAY_BUFFER, 0);
    }
}

void RenderContext::UploadShipVelocity(
    vec2f const & shipVelocity,
    float highlight)
{
    //
    // Create vertices
    //

    mShipVelocityVertexBuffer.clear();

    if (shipVelocity != vec2f::zero())
    {
        // Dimensions - NDC
        float const quadrantWidth = 0.4f;
        float const quadrantHeight = quadrantWidth * mViewModel.GetAspectRatio();
        vec2f const halfQuadrantSize = vec2f(quadrantWidth / 2.0f, quadrantHeight / 2.0f);
        vec2f const center = vec2f(1.0f - quadrantWidth / 2.0f, -1.0f + quadrantHeight / 2.0f);

        // Vector - in normalized space
        float const vectorMagnitude = shipVelocity.length();
        vec2f const vectorDir = shipVelocity.normalise(vectorMagnitude);
        vec2f const vectorNormal = vectorDir.to_perpendicular();
        float constexpr ArrowHeadWidth = 0.3f; // To be kept inline with shader
        float constexpr MinVectorLength = ArrowHeadWidth + 0.04f; // Room for border & anti-aliasing, plus a portion of body
        float constexpr MaxVectorLength = 1.0f;
        float constexpr VectorWidth = ArrowHeadWidth;
        // 0.2 + (1.0-0.2)*(1.0 - e^(-x/5.0))  // Almost 1.0 at 30.0
        float const scaledVectorLength = MinVectorLength + (MaxVectorLength - MinVectorLength) * (1.0f - std::exp(-vectorMagnitude * 0.2f));

        //
        // A     J         L   MinLen    B
        // -------------------------------
        // |     |         |             |
        // |center         |             |
        // |     |         |             |
        // -------------------------------
        // C     K         M             D
        //

        // In NDC space

        vec2f const a = center + vectorNormal * VectorWidth / 2.0f * halfQuadrantSize;
        vec2f const c = center - vectorNormal * VectorWidth / 2.0f * halfQuadrantSize;

        // Left portion

        float constexpr LeftPortionWidth = MinVectorLength * 0.2f;
        vec2f const j = a + vectorDir * LeftPortionWidth * halfQuadrantSize;
        vec2f const k = c + vectorDir * LeftPortionWidth * halfQuadrantSize;

        mShipVelocityVertexBuffer.emplace_back(
            a,
            vec2f(-0.5f, 0.5f),
            highlight);

        mShipVelocityVertexBuffer.emplace_back(
            c,
            vec2f(-0.5f, -0.5f),
            highlight);

        mShipVelocityVertexBuffer.emplace_back(
            j,
            vec2f(-0.25f, 0.5f),
            highlight);

        mShipVelocityVertexBuffer.emplace_back(
            k,
            vec2f(-0.25f, -0.5f),
            highlight);

        // Middle portion

        float const middlePortionWidth = scaledVectorLength - LeftPortionWidth - MinVectorLength;
        vec2f const l = a + vectorDir * (LeftPortionWidth + middlePortionWidth) * halfQuadrantSize;
        vec2f const m = c + vectorDir * (LeftPortionWidth + middlePortionWidth) * halfQuadrantSize;

        mShipVelocityVertexBuffer.emplace_back(
            l,
            vec2f(0.0f, 0.5f),
            highlight);

        mShipVelocityVertexBuffer.emplace_back(
            m,
            vec2f(0.0f, -0.5f),
            highlight);

        // Right portion

        vec2f const b = a + vectorDir * scaledVectorLength * halfQuadrantSize;
        vec2f const d = c + vectorDir * scaledVectorLength * halfQuadrantSize;

        mShipVelocityVertexBuffer.emplace_back(
            b,
            vec2f(0.5f, 0.5f),
            highlight);

        mShipVelocityVertexBuffer.emplace_back(
            d,
            vec2f(0.5f, -0.5f),
            highlight);
    }

    //
    // Upload buffer, if needed
    //

    if (!mShipVelocityVertexBuffer.empty())
    {
        glBindBuffer(GL_ARRAY_BUFFER, *mShipVelocityVertexVBO);
        glBufferData(GL_ARRAY_BUFFER, sizeof(ShipVelocityVertex) * mShipVelocityVertexBuffer.size(), mShipVelocityVertexBuffer.data(), GL_STREAM_DRAW);
        glBindBuffer(GL_ARRAY_BUFFER, 0);
    }
}

void RenderContext::RenderEnd()
{
    ////////////////////////////////////////////////////////////////
    // Element Indices
    ////////////////////////////////////////////////////////////////

    if (mElementIndices->IsDirty())
    {
        mElementIndices->Upload();
    }

    ////////////////////////////////////////////////////////////////
    // Background
    ////////////////////////////////////////////////////////////////

    {
        glBindVertexArray(*mBackgroundVAO);

        mShaderManager->ActivateProgram<ShaderManager::ProgramType::Background>();

        assert(mBackgroundVertexBuffer.size() == 4);
        glDrawArrays(GL_TRIANGLE_STRIP, 0, static_cast<GLsizei>(mBackgroundVertexBuffer.size()));

        CheckOpenGLError();

        glBindVertexArray(0);
    }

    ////////////////////////////////////////////////////////////////
    // Triangles
    ////////////////////////////////////////////////////////////////

    if (!mTriangleVertexBuffer.empty())
    {
        glBindVertexArray(*mTriangleVAO);

        mShaderManager->ActivateProgram<ShaderManager::ProgramType::Triangles>();

        assert((mTriangleVertexBuffer.size() % 3) == 0);
        glDrawArrays(GL_TRIANGLES, 0, static_cast<GLsizei>(mTriangleVertexBuffer.size()));

        CheckOpenGLError();

        glBindVertexArray(0);
    }

    ////////////////////////////////////////////////////////////////
    // Ship velocity
    ////////////////////////////////////////////////////////////////

    if (!mShipVelocityVertexBuffer.empty())
    {
        glBindVertexArray(*mShipVelocityVAO);

        mShaderManager->ActivateProgram<ShaderManager::ProgramType::ShipVelocity>();

        glDrawArrays(GL_TRIANGLE_STRIP, 0, static_cast<GLsizei>(mShipVelocityVertexBuffer.size()));

        CheckOpenGLError();

        glBindVertexArray(0);
    }

    ////////////////////////////////////////////////////////////////
    // Edges
    ////////////////////////////////////////////////////////////////

    if (!mEdgeVertexBuffer.empty())
    {
        glBindVertexArray(*mEdgeVAO);

        mShaderManager->ActivateProgram<ShaderManager::ProgramType::Edges>();

        assert((mEdgeVertexBuffer.size() % 6) == 0);
        glDrawArrays(GL_TRIANGLES, 0, static_cast<GLsizei>(mEdgeVertexBuffer.size()));

        CheckOpenGLError();

        glBindVertexArray(0);
    }

    ////////////////////////////////////////////////////////////////
    // Vertices
    ////////////////////////////////////////////////////////////////

    glBindVertexArray(*mVertexVAO);

    mShaderManager->ActivateProgram<ShaderManager::ProgramType::Vertices>();

    assert((mVertexVertexCount % 6) == 0);
    glDrawArrays(GL_TRIANGLES, 0, static_cast<GLsizei>(mVertexVertexCount));

    CheckOpenGLError();

    glBindVertexArray(0);

    ////////////////////////////////////////////////////////////////
    // NPC Springs
    ////////////////////////////////////////////////////////////////

    if (!mNpcSpringVertexBuffer.empty())
    {
        glBindVertexArray(*mNpcSpringVAO);

        mShaderManager->ActivateProgram<ShaderManager::ProgramType::NpcSprings>();

        assert((mNpcSpringVertexBuffer.size() % 6) == 0);
        glDrawArrays(GL_TRIANGLES, 0, static_cast<GLsizei>(mNpcSpringVertexBuffer.size()));

        CheckOpenGLError();

        glBindVertexArray(0);
    }

    ////////////////////////////////////////////////////////////////
    // NPC Particles
    ////////////////////////////////////////////////////////////////

    if (!mNpcParticleVertexBuffer.empty())
    {
        glBindVertexArray(*mNpcParticleVAO);

        mShaderManager->ActivateProgram<ShaderManager::ProgramType::NpcParticles>();

        assert((mNpcParticleVertexBuffer.size() % 6) == 0);
        glDrawArrays(GL_TRIANGLES, 0, static_cast<GLsizei>(mNpcParticleVertexBuffer.size()));

        CheckOpenGLError();

        glBindVertexArray(0);
    }

    ////////////////////////////////////////////////////////////////
    // NPC Quads
    ////////////////////////////////////////////////////////////////

    if (!mNpcQuadVertexBuffer.empty())
    {
        glBindVertexArray(*mNpcQuadVAO);

        // Intel bug: cannot associate with VAO
        mElementIndices->Bind();

        mShaderManager->ActivateProgram<ShaderManager::ProgramType::NpcQuads>();

        assert((mNpcQuadVertexBuffer.size() % 4) == 0);

        glDrawElements(
            GL_TRIANGLES,
            static_cast<GLsizei>(mNpcQuadVertexBuffer.size() / 4 * 6),
            GL_UNSIGNED_INT,
            (GLvoid *)0);

        CheckOpenGLError();

        glBindVertexArray(0);
    }

    ////////////////////////////////////////////////////////////////
    // Particle trajectories
    ////////////////////////////////////////////////////////////////

    if (!mParticleTrajectoryVertexBuffer.empty())
    {
        glBindVertexArray(*mParticleTrajectoryVAO);

        mShaderManager->ActivateProgram<ShaderManager::ProgramType::NpcParticleTrajectories>();

        assert((mParticleTrajectoryVertexBuffer.size() % 6) == 0);
        glDrawArrays(GL_TRIANGLES, 0, static_cast<GLsizei>(mParticleTrajectoryVertexBuffer.size()));

        CheckOpenGLError();

        glBindVertexArray(0);
    }

    ////////////////////////////////////////////////////////////////
    // Grid
    ////////////////////////////////////////////////////////////////

    if (mIsGridEnabled)
    {
        glBindVertexArray(*mGridVAO);

        mShaderManager->ActivateProgram<ShaderManager::ProgramType::Grid>();

        glDrawArrays(GL_TRIANGLE_STRIP, 0, 4);

        CheckOpenGLError();

        glBindVertexArray(0);
    }

    ////////////////////////////////////////////////////////////////
    // Terminate
    ////////////////////////////////////////////////////////////////

    // Flush all pending commands (but not the GPU buffer)
    GameOpenGL::Flush();
}

//////////////////////////////////////////////////////////////////////////////////////////////

void RenderContext::ProcessSettingChanges()
{
    if (mIsCanvasSizeDirty)
    {
        OnCanvasSizeUpdated();
        mIsCanvasSizeDirty = false;
    }

    if (mIsViewModelDirty)
    {
        OnViewModelUpdated();
        mIsGridDirty = true;
        mIsViewModelDirty = false;
    }

    if (mIsGridDirty)
    {
        OnGridUpdated();
        mIsGridDirty = false;
    }
}

void RenderContext::OnCanvasSizeUpdated()
{
    glViewport(0, 0, mViewModel.GetCanvasWidth(), mViewModel.GetCanvasHeight());
}

void RenderContext::OnViewModelUpdated()
{
    //
    // Update ortho matrix
    //

    ViewModel::ProjectionMatrix const & orthoMatrix = mViewModel.GetOrthoMatrix();

    mShaderManager->ActivateProgram<ShaderManager::ProgramType::Edges>();
    mShaderManager->SetProgramParameter<ShaderManager::ProgramType::Edges, ShaderManager::ProgramParameterType::OrthoMatrix>(
        orthoMatrix);

    mShaderManager->ActivateProgram<ShaderManager::ProgramType::NpcParticles>();
    mShaderManager->SetProgramParameter<ShaderManager::ProgramType::NpcParticles, ShaderManager::ProgramParameterType::OrthoMatrix>(
        orthoMatrix);

    mShaderManager->ActivateProgram<ShaderManager::ProgramType::NpcQuads>();
    mShaderManager->SetProgramParameter<ShaderManager::ProgramType::NpcQuads, ShaderManager::ProgramParameterType::OrthoMatrix>(
        orthoMatrix);

    mShaderManager->ActivateProgram<ShaderManager::ProgramType::NpcSprings>();
    mShaderManager->SetProgramParameter<ShaderManager::ProgramType::NpcSprings, ShaderManager::ProgramParameterType::OrthoMatrix>(
        orthoMatrix);

    mShaderManager->ActivateProgram<ShaderManager::ProgramType::Vertices>();
    mShaderManager->SetProgramParameter<ShaderManager::ProgramType::Vertices, ShaderManager::ProgramParameterType::OrthoMatrix>(
        orthoMatrix);

    mShaderManager->ActivateProgram<ShaderManager::ProgramType::NpcParticleTrajectories>();
    mShaderManager->SetProgramParameter<ShaderManager::ProgramType::NpcParticleTrajectories, ShaderManager::ProgramParameterType::OrthoMatrix>(
        orthoMatrix);

    mShaderManager->ActivateProgram<ShaderManager::ProgramType::Triangles>();
    mShaderManager->SetProgramParameter<ShaderManager::ProgramType::Triangles, ShaderManager::ProgramParameterType::OrthoMatrix>(
        orthoMatrix);

    mShaderManager->ActivateProgram<ShaderManager::ProgramType::Grid>();
    mShaderManager->SetProgramParameter<ShaderManager::ProgramType::Grid, ShaderManager::ProgramParameterType::OrthoMatrix>(
        orthoMatrix);
}

void RenderContext::OnGridUpdated()
{
    //
    // Calculate vertex attributes
    //

    // Visible world coordinates
    vec2f const visibleWorldTopLeft = mViewModel.GetVisibleWorldTopLeft();
    vec2f const visibleWorldBottomRight = mViewModel.GetVisibleWorldBottomRight();

    // Vertices

    std::array<GridVertex, 4> vertexBuffer;

    // Bottom-left
    vertexBuffer[0] = GridVertex(
        vec2f(
            visibleWorldTopLeft.x,
            visibleWorldBottomRight.y));

    // Top left
    vertexBuffer[1] = GridVertex(
        visibleWorldTopLeft);

    // Bottom-right
    vertexBuffer[2] = GridVertex(
        visibleWorldBottomRight);

    // Top-right
    vertexBuffer[3] = GridVertex(
        vec2f(
            visibleWorldBottomRight.x,
            visibleWorldTopLeft.y));

    // Upload vertices
    glBindBuffer(GL_ARRAY_BUFFER, *mGridVBO);
    glBufferData(GL_ARRAY_BUFFER, vertexBuffer.size() * sizeof(GridVertex), vertexBuffer.data(), GL_STATIC_DRAW);
    CheckOpenGLError();
    glBindBuffer(GL_ARRAY_BUFFER, 0);

    //
    // Calculate aspect
    //

    float const pixelWorldWidth = mViewModel.ScreenOffsetToWorldOffset(vec2f(1.0f, -1.0f)).x; // x and y are the same here
    mShaderManager->ActivateProgram<ShaderManager::ProgramType::Grid>();
    mShaderManager->SetProgramParameter<ShaderManager::ProgramType::Grid, ShaderManager::ProgramParameterType::PixelWorldWidth>(
        pixelWorldWidth);

    int constexpr ExtraGridEnlargement = 2;
    float const worldStepSize = std::max(std::ldexp(1.0f, static_cast<int>(std::floor(std::log2f(pixelWorldWidth))) + 2 + ExtraGridEnlargement), 1.0f);

    mShaderManager->ActivateProgram<ShaderManager::ProgramType::Grid>();
    mShaderManager->SetProgramParameter<ShaderManager::ProgramType::Grid, ShaderManager::ProgramParameterType::WorldStep>(
        worldStepSize);
}

}