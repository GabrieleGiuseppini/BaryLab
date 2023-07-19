/***************************************************************************************
* Original Author:		Gabriele Giuseppini
* Created:				2020-05-15
* Copyright:			Gabriele Giuseppini  (https://github.com/GabrieleGiuseppini)
***************************************************************************************/
#include "BLabOpenGL.h"

#include "SysSpecifics.h"

#include <algorithm>
#include <memory>
#include <numeric>

int BLabOpenGL::MaxVertexAttributes = 0;
int BLabOpenGL::MaxViewportWidth = 0;
int BLabOpenGL::MaxViewportHeight = 0;
int BLabOpenGL::MaxTextureSize = 0;
int BLabOpenGL::MaxRenderbufferSize = 0;

void BLabOpenGL::InitOpenGL()
{
    int status = gladLoadGL();
    if (!status)
    {
        throw BLabException("We are sorry, but this game requires OpenGL and it seems your graphics driver does not support it; the error is: failed to initialize GLAD");
    }

    //
    // Check OpenGL version
    //

    LogMessage("OpenGL version: ", GLVersion.major, ".", GLVersion.minor);

    if (GLVersion.major < MinOpenGLVersionMaj
        || (GLVersion.major == MinOpenGLVersionMaj && GLVersion.minor < MinOpenGLVersionMin))
    {
        throw BLabException(
            std::string("We are sorry, but this game requires at least OpenGL ")
            + std::to_string(MinOpenGLVersionMaj) + "." + std::to_string(MinOpenGLVersionMin)
            + ", while the version currently supported by your graphics driver is "
            + std::to_string(GLVersion.major) + "." + std::to_string(GLVersion.minor));
    }


    //
    // Init our extensions
    //

    InitOpenGLExt();


    //
    // Get some constants
    //

    GLint tmpConstant;

    glGetIntegerv(GL_MAX_VERTEX_ATTRIBS, &tmpConstant);
    MaxVertexAttributes = tmpConstant;
    LogMessage("GL_MAX_VERTEX_ATTRIBS=", MaxVertexAttributes);

    GLint maxViewportDims[2];
    glGetIntegerv(GL_MAX_VIEWPORT_DIMS, &(maxViewportDims[0]));
    MaxViewportWidth = maxViewportDims[0];
    MaxViewportHeight = maxViewportDims[1];
    LogMessage("GL_MAX_VIEWPORT_DIMS=", MaxViewportWidth, "x", MaxViewportHeight);

    glGetIntegerv(GL_MAX_TEXTURE_SIZE, &tmpConstant);
    MaxTextureSize = tmpConstant;
    LogMessage("GL_MAX_TEXTURE_SIZE=", MaxTextureSize);

    glGetIntegerv(GL_MAX_RENDERBUFFER_SIZE, &tmpConstant);
    MaxRenderbufferSize = tmpConstant;
    LogMessage("GL_MAX_RENDERBUFFER_SIZE=", MaxRenderbufferSize);
}

void BLabOpenGL::CompileShader(
    std::string const & shaderSource,
    GLenum shaderType,
    BLabOpenGLShaderProgram const & shaderProgram,
    std::string const & programName)
{
    char const * shaderSourceCString = shaderSource.c_str();
    std::string const shaderTypeName = (shaderType == GL_VERTEX_SHADER) ? "vertex" : "fragment";

    // Set source
    GLuint shader = glCreateShader(shaderType);
    glShaderSource(shader, 1, &shaderSourceCString, NULL);
    GLenum glError = glGetError();
    if (GL_NO_ERROR != glError)
    {
        throw BLabException("Error setting " + shaderTypeName + " shader source for program \"" + programName + "\"");
    }

    // Compile
    glCompileShader(shader);
    GLint success;
    glGetShaderiv(shader, GL_COMPILE_STATUS, &success);
    if (GL_FALSE == success)
    {
        char infoLog[1024];
        glGetShaderInfoLog(shader, sizeof(infoLog) - 1, NULL, infoLog);
        throw BLabException("Error compiling " + shaderTypeName + " shader: " + std::string(infoLog));
    }

    // Attach to program
    glAttachShader(*shaderProgram, shader);
    glError = glGetError();
    if (GL_NO_ERROR != glError)
    {
        throw BLabException("Error attaching compiled " + shaderTypeName + " shader to program \"" + programName + "\"");
    }

    // Delete shader
    glDeleteShader(shader);
}

void BLabOpenGL::LinkShaderProgram(
    BLabOpenGLShaderProgram const & shaderProgram,
    std::string const & programName)
{
    glLinkProgram(*shaderProgram);

    // Check
    int success;
    glGetProgramiv(*shaderProgram, GL_LINK_STATUS, &success);
    if (!success)
    {
        char infoLog[1024];
        glGetShaderInfoLog(*shaderProgram, sizeof(infoLog), NULL, infoLog);
        throw BLabException("Error linking " + programName + " shader program: " + std::string(infoLog));
    }
}

GLint BLabOpenGL::GetParameterLocation(
    BLabOpenGLShaderProgram const & shaderProgram,
    std::string const & parameterName)
{
    GLint parameterLocation = glGetUniformLocation(*shaderProgram, parameterName.c_str());

    GLenum glError = glGetError();
    if (parameterLocation == -1
        || GL_NO_ERROR != glError)
    {
        throw BLabException("Cannot retrieve location of parameter \"" + parameterName + "\"");
    }

    return parameterLocation;
}

void BLabOpenGL::BindAttributeLocation(
    BLabOpenGLShaderProgram const & shaderProgram,
    GLuint attributeIndex,
    std::string const & attributeName)
{
    glBindAttribLocation(
        *shaderProgram,
        attributeIndex,
        attributeName.c_str());

    GLenum glError = glGetError();
    if (GL_NO_ERROR != glError)
    {
        throw BLabException("Error binding attribute location for attribute \"" + attributeName + "\"");
    }
}

void BLabOpenGL::Flush()
{
    // We do it here to have this call in the stack, helping
    // performance profiling
    glFlush();
}