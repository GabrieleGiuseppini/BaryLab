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
    // Springs
    //

    glGenVertexArrays(1, &tmpGLuint);
    mSpringVAO = tmpGLuint;
    glBindVertexArray(*mSpringVAO);

    glGenBuffers(1, &tmpGLuint);
    mSpringVertexVBO = tmpGLuint;
    glBindBuffer(GL_ARRAY_BUFFER, *mSpringVertexVBO);

    glEnableVertexAttribArray(static_cast<GLuint>(ShaderManager::VertexAttributeType::SpringAttributeGroup1));
    glVertexAttribPointer(static_cast<GLuint>(ShaderManager::VertexAttributeType::SpringAttributeGroup1), 4, GL_FLOAT, GL_FALSE, sizeof(SpringVertex), (void *)0);
    glEnableVertexAttribArray(static_cast<GLuint>(ShaderManager::VertexAttributeType::SpringAttributeGroup2));
    glVertexAttribPointer(static_cast<GLuint>(ShaderManager::VertexAttributeType::SpringAttributeGroup2), 4, GL_FLOAT, GL_FALSE, sizeof(SpringVertex), (void *)(4 * sizeof(float)));
    static_assert(sizeof(SpringVertex) == 8 * sizeof(float));

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
    glEnableVertexAttribArray(static_cast<GLuint>(ShaderManager::VertexAttributeType::ParticleAttributeGroup2));
    glVertexAttribPointer(static_cast<GLuint>(ShaderManager::VertexAttributeType::ParticleAttributeGroup2), 4, GL_FLOAT, GL_FALSE, sizeof(ParticleTrajectoryVertex), (void *)(4 * sizeof(float)));
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
    // Mesh velocity
    //

    glGenVertexArrays(1, &tmpGLuint);
    mMeshVelocityVAO = tmpGLuint;
    glBindVertexArray(*mMeshVelocityVAO);

    glGenBuffers(1, &tmpGLuint);
    mMeshVelocityVertexVBO = tmpGLuint;
    glBindBuffer(GL_ARRAY_BUFFER, *mMeshVelocityVertexVBO);

    glEnableVertexAttribArray(static_cast<GLuint>(ShaderManager::VertexAttributeType::MeshVelocityAttributeGroup1));
    glVertexAttribPointer(static_cast<GLuint>(ShaderManager::VertexAttributeType::MeshVelocityAttributeGroup1), 4, GL_FLOAT, GL_FALSE, sizeof(MeshVelocityVertex), (void *)0);
    glEnableVertexAttribArray(static_cast<GLuint>(ShaderManager::VertexAttributeType::MeshVelocityAttributeGroup2));
    glVertexAttribPointer(static_cast<GLuint>(ShaderManager::VertexAttributeType::MeshVelocityAttributeGroup2), 1, GL_FLOAT, GL_FALSE, sizeof(MeshVelocityVertex), (void *)(4 * sizeof(float)));
    static_assert(sizeof(MeshVelocityVertex) == (4 + 1) * sizeof(float));

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

void RenderContext::UploadSpringsStart()
{
    //
    // Prepare buffer
    //

    mSpringVertexBuffer.clear();
}

void RenderContext::UploadSpring(
    vec2f const & endpointAPosition,
    vec2f const & endpointBPosition,
    rgbaColor const & springColor)
{
    vec2f const springVector = endpointBPosition - endpointAPosition;
    vec2f const springNormal = springVector.to_perpendicular().normalise() * LabParameters::SpringThickness / 2.0f;

    vec2f const bottomLeft = endpointAPosition - springNormal;
    vec2f const bottomRight = endpointAPosition + springNormal;
    vec2f const topLeft = endpointBPosition - springNormal;
    vec2f const topRight = endpointBPosition + springNormal;

    vec4f const color = springColor.toVec4f();

    // Left, bottom
    mSpringVertexBuffer.emplace_back(
        bottomLeft,
        vec2f(-1.0f, -1.0f),
        color);

    // Left, top
    mSpringVertexBuffer.emplace_back(
        topLeft,
        vec2f(-1.0f, 1.0f),
        color);

    // Right, bottom
    mSpringVertexBuffer.emplace_back(
        bottomRight,
        vec2f(1.0f, -1.0f),
        color);

    // Left, top
    mSpringVertexBuffer.emplace_back(
        topLeft,
        vec2f(-1.0f, 1.0f),
        color);

    // Right, bottom
    mSpringVertexBuffer.emplace_back(
        bottomRight,
        vec2f(1.0f, -1.0f),
        color);

    // Right, top
    mSpringVertexBuffer.emplace_back(
        topRight,
        vec2f(1.0f, 1.0f),
        color);
}

void RenderContext::UploadSpringsEnd()
{
    //
    // Upload buffer, if needed
    //

    if (!mSpringVertexBuffer.empty())
    {
        glBindBuffer(GL_ARRAY_BUFFER, *mSpringVertexVBO);
        glBufferData(GL_ARRAY_BUFFER, sizeof(SpringVertex) * mSpringVertexBuffer.size(), mSpringVertexBuffer.data(), GL_STREAM_DRAW);
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
    rgbaColor const & particleColor,
    float alpha)
{
    float const xLeft = particlePosition.x - LabParameters::ParticleRadius;
    float const xRight = particlePosition.x + LabParameters::ParticleRadius;
    float const yTop = particlePosition.y + LabParameters::ParticleRadius;
    float const yBottom = particlePosition.y - LabParameters::ParticleRadius;

    vec4f color = particleColor.toVec4f();
    color.w *= alpha;

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
    vec2f const & endPosition,
    rgbaColor const & color)
{
    vec2f const trajectoryVector = endPosition - startPosition;
    vec2f const trajectoryNormal = trajectoryVector.to_perpendicular().normalise() * LabParameters::ParticleTrajectoryThickness / 2.0f;

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

void RenderContext::UploadMeshVelocity(
    vec2f const & meshVelocity,
    float highlight)
{
    //
    // Create vertices
    //

    mMeshVelocityVertexBuffer.clear();

    if (meshVelocity != vec2f::zero())
    {
        // Dimensions - NDC
        float const quadrantWidth = 0.4f;
        float const quadrantHeight = quadrantWidth * mViewModel.GetAspectRatio();
        vec2f const halfQuadrantSize = vec2f(quadrantWidth / 2.0f, quadrantHeight / 2.0f);
        vec2f const center = vec2f(1.0f - quadrantWidth / 2.0f, -1.0f + quadrantHeight / 2.0f);

        // Vector - in normalized space
        float const vectorMagnitude = meshVelocity.length();
        vec2f const vectorDir = meshVelocity.normalise(vectorMagnitude);
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

        mMeshVelocityVertexBuffer.emplace_back(
            a,
            vec2f(-0.5f, 0.5f),
            highlight);

        mMeshVelocityVertexBuffer.emplace_back(
            c,
            vec2f(-0.5f, -0.5f),
            highlight);

        mMeshVelocityVertexBuffer.emplace_back(
            j,
            vec2f(-0.25f, 0.5f),
            highlight);

        mMeshVelocityVertexBuffer.emplace_back(
            k,
            vec2f(-0.25f, -0.5f),
            highlight);

        // Middle portion

        float const middlePortionWidth = scaledVectorLength - LeftPortionWidth - MinVectorLength;
        vec2f const l = a + vectorDir * (LeftPortionWidth + middlePortionWidth) * halfQuadrantSize;
        vec2f const m = c + vectorDir * (LeftPortionWidth + middlePortionWidth) * halfQuadrantSize;

        mMeshVelocityVertexBuffer.emplace_back(
            l,
            vec2f(0.0f, 0.5f),
            highlight);

        mMeshVelocityVertexBuffer.emplace_back(
            m,
            vec2f(0.0f, -0.5f),
            highlight);

        // Right portion

        vec2f const b = a + vectorDir * scaledVectorLength * halfQuadrantSize;
        vec2f const d = c + vectorDir * scaledVectorLength * halfQuadrantSize;

        mMeshVelocityVertexBuffer.emplace_back(
            b,
            vec2f(0.5f, 0.5f),
            highlight);

        mMeshVelocityVertexBuffer.emplace_back(
            d,
            vec2f(0.5f, -0.5f),
            highlight);
    }

    //
    // Upload buffer, if needed
    //

    if (!mMeshVelocityVertexBuffer.empty())
    {
        glBindBuffer(GL_ARRAY_BUFFER, *mMeshVelocityVertexVBO);
        glBufferData(GL_ARRAY_BUFFER, sizeof(MeshVelocityVertex) * mMeshVelocityVertexBuffer.size(), mMeshVelocityVertexBuffer.data(), GL_STREAM_DRAW);
        glBindBuffer(GL_ARRAY_BUFFER, 0);
    }
}

void RenderContext::RenderEnd()
{
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
    // Mesh velocity
    ////////////////////////////////////////////////////////////////

    if (!mMeshVelocityVertexBuffer.empty())
    {
        glBindVertexArray(*mMeshVelocityVAO);

        mShaderManager->ActivateProgram<ShaderManager::ProgramType::MeshVelocity>();

        glDrawArrays(GL_TRIANGLE_STRIP, 0, static_cast<GLsizei>(mMeshVelocityVertexBuffer.size()));

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
    // Springs
    ////////////////////////////////////////////////////////////////

    if (!mSpringVertexBuffer.empty())
    {
        glBindVertexArray(*mSpringVAO);

        mShaderManager->ActivateProgram<ShaderManager::ProgramType::Springs>();

        assert((mSpringVertexBuffer.size() % 6) == 0);
        glDrawArrays(GL_TRIANGLES, 0, static_cast<GLsizei>(mSpringVertexBuffer.size()));

        CheckOpenGLError();

        glBindVertexArray(0);
    }

    ////////////////////////////////////////////////////////////////
    // Particles
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
    // Particle trajectories
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

    mShaderManager->ActivateProgram<ShaderManager::ProgramType::Springs>();
    mShaderManager->SetProgramParameter<ShaderManager::ProgramType::Springs, ShaderManager::ProgramParameterType::OrthoMatrix>(
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
