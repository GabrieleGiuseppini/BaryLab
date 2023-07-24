/***************************************************************************************
* Original Author:		Gabriele Giuseppini
* Created:				2020-05-16
* Copyright:			Gabriele Giuseppini  (https://github.com/GabrieleGiuseppini)
***************************************************************************************/
#include "RenderContext.h"

#include <cmath>

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
{
    GLuint tmpGLuint;

    //
    // Initialize OpenGL
    //

    try
    {
        BLabOpenGL::InitOpenGL();
    }
    catch (std::exception const & e)
    {
        throw std::runtime_error("Error during OpenGL initialization: " + std::string(e.what()));
    }

    ////////////////////////////////////////////////////////////////
    // Initialize shaders, VAO's, and VBOs
    ////////////////////////////////////////////////////////////////

    mShaderManager = ShaderManager::CreateInstance();

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
    // Particles
    //

    glGenVertexArrays(1, &tmpGLuint);
    mParticleVAO = tmpGLuint;
    glBindVertexArray(*mParticleVAO);

    glGenBuffers(1, &tmpGLuint);
    mParticleVertexVBO = tmpGLuint;
    glBindBuffer(GL_ARRAY_BUFFER, *mParticleVertexVBO);

    glEnableVertexAttribArray(static_cast<GLuint>(ShaderManager::VertexAttributeType::ParticleAttributeGroup1));
    glVertexAttribPointer(static_cast<GLuint>(ShaderManager::VertexAttributeType::ParticleAttributeGroup1), 4, GL_FLOAT, GL_FALSE, sizeof(ParticleVertex), (void *)0);
    glEnableVertexAttribArray(static_cast<GLuint>(ShaderManager::VertexAttributeType::ParticleAttributeGroup2));
    glVertexAttribPointer(static_cast<GLuint>(ShaderManager::VertexAttributeType::ParticleAttributeGroup2), 4, GL_FLOAT, GL_FALSE, sizeof(ParticleVertex), (void *)(4 * sizeof(float)));
    static_assert(sizeof(ParticleVertex) == (4 + 4) * sizeof(float));

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

    glEnableVertexAttribArray(static_cast<GLuint>(ShaderManager::VertexAttributeType::ParticleAttributeGroup1));
    glVertexAttribPointer(static_cast<GLuint>(ShaderManager::VertexAttributeType::ParticleAttributeGroup1), 4, GL_FLOAT, GL_FALSE, sizeof(ParticleTrajectoryVertex), (void *)0);
    static_assert(sizeof(ParticleTrajectoryVertex) == 4 * sizeof(float));

    glBindVertexArray(0);

    //
    // Selected triangles
    //

    glGenVertexArrays(1, &tmpGLuint);
    mSelectedTriangleVAO = tmpGLuint;
    glBindVertexArray(*mSelectedTriangleVAO);

    glGenBuffers(1, &tmpGLuint);
    mSelectedTriangleVertexVBO = tmpGLuint;
    glBindBuffer(GL_ARRAY_BUFFER, *mSelectedTriangleVertexVBO);

    glEnableVertexAttribArray(static_cast<GLuint>(ShaderManager::VertexAttributeType::SelectedTriangleAttributeGroup1));
    glVertexAttribPointer(static_cast<GLuint>(ShaderManager::VertexAttributeType::SelectedTriangleAttributeGroup1), 2, GL_FLOAT, GL_FALSE, sizeof(SelectedTriangleVertex), (void *)0);
    static_assert(sizeof(SelectedTriangleVertex) == 2 * sizeof(float));

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

    ////////////////////////////////////////////////////////////////
    // Set parameters in all shaders
    ////////////////////////////////////////////////////////////////

    ProcessSettingChanges();
}

void RenderContext::RenderStart()
{
    // Set polygon mode
    glPolygonMode(GL_FRONT_AND_BACK, GL_FILL);

    // Clear canvas - and depth buffer
    //vec3f constexpr ClearColor = rgbColor(0xca, 0xf4, 0xf4).toVec3f();
    vec3f constexpr ClearColor = rgbColor(0xff, 0xff, 0xff).toVec3f();
    glClearColor(ClearColor.x, ClearColor.y, ClearColor.z, 1.0f);
    glClear(GL_COLOR_BUFFER_BIT | GL_DEPTH_BUFFER_BIT);

    // Process setting changes
    ProcessSettingChanges();
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

        float const xLeft = vertexPosition.x - LabParameters::VertexRadius;
        float const xRight = vertexPosition.x + LabParameters::VertexRadius;
        float const yTop = vertexPosition.y + LabParameters::VertexRadius;
        float const yBottom = vertexPosition.y - LabParameters::VertexRadius;

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
    rgbaColor const & edgeColor)
{
    vec2f const edgeVector = edgeEndpointBPosition - edgeEndpointAPosition;
    vec2f const edgeNormal = edgeVector.to_perpendicular().normalise() * LabParameters::EdgeThickness / 2.0f;

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

void RenderContext::UploadParticlesStart()
{
    //
    // Prepare buffer
    //

    mParticleVertexBuffer.clear();
}

void RenderContext::UploadParticle(
    vec2f const & particlePosition,
    rgbaColor const & particleColor)
{
    float const xLeft = particlePosition.x - LabParameters::ParticleRadius;
    float const xRight = particlePosition.x + LabParameters::ParticleRadius;
    float const yTop = particlePosition.y + LabParameters::ParticleRadius;
    float const yBottom = particlePosition.y - LabParameters::ParticleRadius;

    vec4f const color = particleColor.toVec4f();

    // Left, bottom
    mParticleVertexBuffer.emplace_back(
        vec2f(xLeft, yBottom),
        vec2f(-1.0f, -1.0f),
        color);

    // Left, top
    mParticleVertexBuffer.emplace_back(
        vec2f(xLeft, yTop),
        vec2f(-1.0f, 1.0f),
        color);

    // Right, bottom
    mParticleVertexBuffer.emplace_back(
        vec2f(xRight, yBottom),
        vec2f(1.0f, -1.0f),
        color);

    // Left, top
    mParticleVertexBuffer.emplace_back(
        vec2f(xLeft, yTop),
        vec2f(-1.0f, 1.0f),
        color);

    // Right, bottom
    mParticleVertexBuffer.emplace_back(
        vec2f(xRight, yBottom),
        vec2f(1.0f, -1.0f),
        color);

    // Right, top
    mParticleVertexBuffer.emplace_back(
        vec2f(xRight, yTop),
        vec2f(1.0f, 1.0f),
        color);
}

void RenderContext::UploadParticlesEnd()
{
    //
    // Upload buffer, if needed
    //

    if (!mParticleVertexBuffer.empty())
    {
        glBindBuffer(GL_ARRAY_BUFFER, *mParticleVertexVBO);
        glBufferData(GL_ARRAY_BUFFER, sizeof(ParticleVertex) * mParticleVertexBuffer.size(), mParticleVertexBuffer.data(), GL_STREAM_DRAW);
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
    vec2f const & endPosition)
{
    vec2f const trajectoryVector = endPosition - startPosition;
    vec2f const edgeNormal = trajectoryVector.to_perpendicular().normalise() * LabParameters::ParticleTrajectoryThickness / 2.0f;

    vec2f const bottomLeft = startPosition - edgeNormal;
    vec2f const bottomRight = startPosition + edgeNormal;
    vec2f const topLeft = endPosition - edgeNormal;
    vec2f const topRight = endPosition + edgeNormal;

    // Left, bottom
    mParticleTrajectoryVertexBuffer.emplace_back(
        bottomLeft,
        vec2f(-1.0f, -1.0f));

    // Left, top
    mParticleTrajectoryVertexBuffer.emplace_back(
        topLeft,
        vec2f(-1.0f, 1.0f));

    // Right, bottom
    mParticleTrajectoryVertexBuffer.emplace_back(
        bottomRight,
        vec2f(1.0f, -1.0f));

    // Left, top
    mParticleTrajectoryVertexBuffer.emplace_back(
        topLeft,
        vec2f(-1.0f, 1.0f));

    // Right, bottom
    mParticleTrajectoryVertexBuffer.emplace_back(
        bottomRight,
        vec2f(1.0f, -1.0f));

    // Right, top
    mParticleTrajectoryVertexBuffer.emplace_back(
        topRight,
        vec2f(1.0f, 1.0f));
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

void RenderContext::UploadSelectedTrianglesStart()
{
    //
    // Prepare buffer
    //

    mSelectedTriangleVertexBuffer.clear();
}

void RenderContext::UploadSelectedTriangle(
    vec2f const & endpointAPosition,
    vec2f const & endpointBPosition,
    vec2f const & endpointCPosition)
{
    mSelectedTriangleVertexBuffer.emplace_back(endpointAPosition);
    mSelectedTriangleVertexBuffer.emplace_back(endpointBPosition);
    mSelectedTriangleVertexBuffer.emplace_back(endpointCPosition);
}

void RenderContext::UploadSelectedTrianglesEnd()
{
    //
    // Upload buffer, if needed
    //

    if (!mSelectedTriangleVertexBuffer.empty())
    {
        glBindBuffer(GL_ARRAY_BUFFER, *mSelectedTriangleVertexVBO);
        glBufferData(GL_ARRAY_BUFFER, sizeof(SelectedTriangleVertex) * mSelectedTriangleVertexBuffer.size(), mSelectedTriangleVertexBuffer.data(), GL_STREAM_DRAW);
        glBindBuffer(GL_ARRAY_BUFFER, 0);
    }
}

void RenderContext::RenderEnd()
{
    ////////////////////////////////////////////////////////////////
    // Render selected triangles
    ////////////////////////////////////////////////////////////////

    if (!mSelectedTriangleVertexBuffer.empty())
    {
        glBindVertexArray(*mSelectedTriangleVAO);

        mShaderManager->ActivateProgram<ShaderManager::ProgramType::SelectedTriangles>();

        assert((mSelectedTriangleVertexBuffer.size() % 3) == 0);
        glDrawArrays(GL_TRIANGLES, 0, static_cast<GLsizei>(mSelectedTriangleVertexBuffer.size()));

        CheckOpenGLError();

        glBindVertexArray(0);
    }

    ////////////////////////////////////////////////////////////////
    // Render edges
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
    // Render particle trajectories
    ////////////////////////////////////////////////////////////////

    if (!mParticleTrajectoryVertexBuffer.empty())
    {
        glBindVertexArray(*mParticleTrajectoryVAO);

        mShaderManager->ActivateProgram<ShaderManager::ProgramType::ParticleTrajectories>();

        assert((mParticleTrajectoryVertexBuffer.size() % 6) == 0);
        glDrawArrays(GL_TRIANGLES, 0, static_cast<GLsizei>(mParticleTrajectoryVertexBuffer.size()));

        CheckOpenGLError();

        glBindVertexArray(0);
    }

    ////////////////////////////////////////////////////////////////
    // Render vertices
    ////////////////////////////////////////////////////////////////

    glBindVertexArray(*mVertexVAO);

    mShaderManager->ActivateProgram<ShaderManager::ProgramType::Vertices>();

    assert((mVertexVertexCount % 6) == 0);
    glDrawArrays(GL_TRIANGLES, 0, static_cast<GLsizei>(mVertexVertexCount));

    CheckOpenGLError();

    glBindVertexArray(0);

    ////////////////////////////////////////////////////////////////
    // Render particles
    ////////////////////////////////////////////////////////////////

    if (!mParticleVertexBuffer.empty())
    {
        glBindVertexArray(*mParticleVAO);

        mShaderManager->ActivateProgram<ShaderManager::ProgramType::Particles>();

        assert((mParticleVertexBuffer.size() % 6) == 0);
        glDrawArrays(GL_TRIANGLES, 0, static_cast<GLsizei>(mParticleVertexBuffer.size()));

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
    BLabOpenGL::Flush();
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

    mShaderManager->ActivateProgram<ShaderManager::ProgramType::Particles>();
    mShaderManager->SetProgramParameter<ShaderManager::ProgramType::Particles, ShaderManager::ProgramParameterType::OrthoMatrix>(
        orthoMatrix);

    mShaderManager->ActivateProgram<ShaderManager::ProgramType::Vertices>();
    mShaderManager->SetProgramParameter<ShaderManager::ProgramType::Vertices, ShaderManager::ProgramParameterType::OrthoMatrix>(
        orthoMatrix);

    mShaderManager->ActivateProgram<ShaderManager::ProgramType::ParticleTrajectories>();
    mShaderManager->SetProgramParameter<ShaderManager::ProgramType::ParticleTrajectories, ShaderManager::ProgramParameterType::OrthoMatrix>(
        orthoMatrix);

    mShaderManager->ActivateProgram<ShaderManager::ProgramType::SelectedTriangles>();
    mShaderManager->SetProgramParameter<ShaderManager::ProgramType::SelectedTriangles, ShaderManager::ProgramParameterType::OrthoMatrix>(
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
