/***************************************************************************************
* Original Author:      Gabriele Giuseppini
* Created:              2020-05-16
* Copyright:            Gabriele Giuseppini  (https://github.com/GabrieleGiuseppini)
***************************************************************************************/
#pragma once

#include <Game/LabController.h>

#include <wx/image.h>
#include <wx/window.h>

#include <cassert>
#include <chrono>
#include <memory>
#include <optional>
#include <vector>

enum class ToolType
{
    MoveParticle = 0,
    MoveVertex = 1,
    RotateMeshByPosition = 2,
    RotateMeshByParticle = 3,
    SetParticleTrajectory = 4,
    SetOriginTriangle = 5,
    SelectParticle = 6
};

struct InputState
{
    bool IsLeftMouseDown;
    bool IsRightMouseDown;
    bool IsShiftKeyDown;
    vec2f MousePosition;
    vec2f PreviousMousePosition;

    InputState()
        : IsLeftMouseDown(false)
        , IsRightMouseDown(false)
        , IsShiftKeyDown(false)
        , MousePosition()
        , PreviousMousePosition()
    {
    }
};

/*
 * Base abstract class of all tools.
 */
class Tool
{
public:

    virtual ~Tool() {}

    ToolType GetToolType() const { return mToolType; }

    virtual void Initialize(InputState const & inputState) = 0;
    virtual void Deinitialize(InputState const & inputState) = 0;
    virtual void SetCurrentCursor() = 0;

    virtual void Update(InputState const & inputState) = 0;

protected:

    Tool(
        ToolType toolType,
        wxWindow * cursorWindow,
        std::shared_ptr<LabController> labController)
        : mCursorWindow(cursorWindow)
        , mLabController(std::move(labController))
        , mToolType(toolType)
    {}

    wxWindow * const mCursorWindow;
    std::shared_ptr<LabController> const mLabController;

private:

    ToolType const mToolType;
};


//////////////////////////////////////////////////////////////////////////////////////////
// Tools
//////////////////////////////////////////////////////////////////////////////////////////

class MoveParticleTool final : public Tool
{
public:

    MoveParticleTool(
        wxWindow * cursorWindow,
        std::shared_ptr<LabController> labController);

public:

    virtual void Initialize(InputState const & /*inputState*/) override
    {
        mCurrentEngagementState.reset();

        // Set cursor
        SetCurrentCursor();
    }

    virtual void Deinitialize(InputState const & /*inputState*/) override
    {
    }

    virtual void SetCurrentCursor() override
    {
        mCursorWindow->SetCursor(!!mCurrentEngagementState ? mDownCursor : mUpCursor);
    }

    virtual void Update(InputState const & inputState) override
    {
        bool const wasEngaged = !!mCurrentEngagementState;

        vec2f const newPosition = inputState.MousePosition;

        if (inputState.IsLeftMouseDown)
        {
            if (!mCurrentEngagementState)
            {
                //
                // Not engaged...
                // ...see if we're able to pick a particle and thus start engagement
                //

                auto const particleId = mLabController->TryPickNpcParticle(newPosition);
                if (particleId.has_value())
                {
                    //
                    // Engage!
                    //

                    mCurrentEngagementState.emplace(*particleId, newPosition);
                }
            }
            else
            {
                //
                // Engaged
                //

                vec2f const stride = newPosition - mCurrentEngagementState->LastPosition;

                mLabController->MoveNpcParticleBy(
                    mCurrentEngagementState->ParticleIndex,
                    stride,
                    vec2f::zero());

                // Update state
                mCurrentEngagementState->LastPosition = newPosition;
                mCurrentEngagementState->LastStride = stride;
            }
        }
        else
        {
            if (mCurrentEngagementState)
            {
                // Impart intertia
                mLabController->MoveNpcParticleBy(
                    mCurrentEngagementState->ParticleIndex,
                    vec2f::zero(),
                    mCurrentEngagementState->LastStride);

                // Disengage
                mCurrentEngagementState.reset();
            }
        }

        if (!!mCurrentEngagementState != wasEngaged)
        {
            // State change

            // Update cursor
            SetCurrentCursor();
        }
    }

private:

    // Our state

    struct EngagementState
    {
        ElementIndex ParticleIndex;
        vec2f LastPosition;
        vec2f LastStride;

        explicit EngagementState(
            ElementIndex particleIndex,
            vec2f const & currentPosition)
            : ParticleIndex(particleIndex)
            , LastPosition(currentPosition)
            , LastStride(vec2f::zero())
        {}
    };

    std::optional<EngagementState> mCurrentEngagementState; // When set, indicates it's engaged

    // The cursors
    wxCursor const mUpCursor;
    wxCursor const mDownCursor;
};

class MoveVertexTool final : public Tool
{
public:

    MoveVertexTool(
        wxWindow * cursorWindow,
        std::shared_ptr<LabController> labController);

public:

    virtual void Initialize(InputState const & /*inputState*/) override
    {
        mCurrentEngagementState.reset();

        // Set cursor
        SetCurrentCursor();
    }

    virtual void Deinitialize(InputState const & /*inputState*/) override
    {
    }

    virtual void SetCurrentCursor() override
    {
        mCursorWindow->SetCursor(!!mCurrentEngagementState ? mDownCursor : mUpCursor);
    }

    virtual void Update(InputState const & inputState) override
    {
        bool const wasEngaged = !!mCurrentEngagementState;

        vec2f const newPosition = inputState.MousePosition;

        if (inputState.IsLeftMouseDown)
        {
            if (!mCurrentEngagementState)
            {
                //
                // Not engaged...
                // ...see if we're able to pick a point and thus start engagement
                //

                auto const vertexId = mLabController->TryPickVertex(newPosition);
                if (vertexId.has_value())
                {
                    //
                    // Engage!
                    //

                    mCurrentEngagementState.emplace(*vertexId, newPosition);
                }
            }
            else
            {
                //
                // Engaged
                //

                vec2f const stride = newPosition - mCurrentEngagementState->LastPosition;

                mLabController->MoveVertexBy(
                    mCurrentEngagementState->VertexIndex,
                    stride);

                // Update state
                mCurrentEngagementState->LastPosition = newPosition;
            }
        }
        else
        {
            if (mCurrentEngagementState)
            {
                // Disengage
                mCurrentEngagementState.reset();
            }
        }

        if (!!mCurrentEngagementState != wasEngaged)
        {
            // State change

            // Update cursor
            SetCurrentCursor();
        }
    }

private:

    // Our state

    struct EngagementState
    {
        ElementIndex VertexIndex;
        vec2f LastPosition;

        explicit EngagementState(
            ElementIndex vertexIndex,
            vec2f const & currentPosition)
            : VertexIndex(vertexIndex)
            , LastPosition(currentPosition)
        {}
    };

    std::optional<EngagementState> mCurrentEngagementState; // When set, indicates it's engaged

    // The cursors
    wxCursor const mUpCursor;
    wxCursor const mDownCursor;
};

class RotateMeshByPositionTool final : public Tool
{
public:

    RotateMeshByPositionTool(
        wxWindow * cursorWindow,
        std::shared_ptr<LabController> labController);

public:

    virtual void Initialize(InputState const & /*inputState*/) override
    {
        // Reset state
        mCurrentEngagementState.reset();

        SetCurrentCursor();
    }

    virtual void Deinitialize(InputState const & /*inputState*/) override
    {
    }

    virtual void SetCurrentCursor() override
    {
        mCursorWindow->SetCursor(!!mCurrentEngagementState ? mDownCursor : mUpCursor);
    }

    virtual void Update(InputState const & inputState) override
    {
        if (mCurrentEngagementState.has_value())
        {
            if (inputState.IsLeftMouseDown)
            {
                // Rotate
                mLabController->RotateShipBy(
                    mCurrentEngagementState->RotationCenter,
                    inputState.MousePosition.y - mCurrentEngagementState->LastPosition.y);

                mCurrentEngagementState->LastPosition = inputState.MousePosition;
            }
            else
            {
                // Disengage
                mCurrentEngagementState.reset();
            }
        }
        else
        {
            if (inputState.IsLeftMouseDown)
            {
                // Engage
                mCurrentEngagementState.emplace(inputState.MousePosition);
            }
        }

        SetCurrentCursor();
    }

private:

    struct EngagementData
    {
        vec2f RotationCenter;
        vec2f LastPosition;

        EngagementData(vec2f const initialPosition)
            : RotationCenter(initialPosition)
            , LastPosition(initialPosition)
        {}
    };

    // When set, we're rotating
    std::optional<EngagementData> mCurrentEngagementState;

    // The cursors
    wxCursor const mUpCursor;
    wxCursor const mDownCursor;
};

class RotateMeshByParticleTool final : public Tool
{
public:

    RotateMeshByParticleTool(
        wxWindow * cursorWindow,
        std::shared_ptr<LabController> labController);

public:

    virtual void Initialize(InputState const & /*inputState*/) override
    {
        // Reset state
        mCurrentEngagementState.reset();

        SetCurrentCursor();
    }

    virtual void Deinitialize(InputState const & /*inputState*/) override
    {
    }

    virtual void SetCurrentCursor() override
    {
        mCursorWindow->SetCursor(!!mCurrentEngagementState ? mDownCursor : mUpCursor);
    }

    virtual void Update(InputState const & inputState) override
    {
        if (mCurrentEngagementState.has_value())
        {
            if (inputState.IsLeftMouseDown)
            {
                // Rotate
                mLabController->RotateShipBy(
                    0,
                    inputState.MousePosition.y - mCurrentEngagementState->LastPosition.y);

                mCurrentEngagementState->LastPosition = inputState.MousePosition;
            }
            else
            {
                // Disengage
                mCurrentEngagementState.reset();
            }
        }
        else
        {
            if (inputState.IsLeftMouseDown)
            {
                // Engage
                mCurrentEngagementState.emplace(inputState.MousePosition);
            }
        }

        SetCurrentCursor();
    }

private:

    struct EngagementData
    {
        vec2f LastPosition;

        EngagementData(vec2f const initialPosition)
            : LastPosition(initialPosition)
        {}
    };

    // When set, we're rotating
    std::optional<EngagementData> mCurrentEngagementState;

    // The cursors
    wxCursor const mUpCursor;
    wxCursor const mDownCursor;
};

class SetParticleTrajectoryTool final : public Tool
{
public:

    SetParticleTrajectoryTool(
        wxWindow * cursorWindow,
        std::shared_ptr<LabController> labController);

public:

    virtual void Initialize(InputState const & /*inputState*/) override
    {
        mCurrentSelectedParticle.reset();

        // Set cursor
        SetCurrentCursor();
    }

    virtual void Deinitialize(InputState const & /*inputState*/) override
    {
    }

    virtual void SetCurrentCursor() override
    {
        mCursorWindow->SetCursor(mCursor);
    }

    virtual void Update(InputState const & inputState) override
    {
        vec2f const mousePosition = inputState.MousePosition;

        if (inputState.IsLeftMouseDown)
        {
            if (!mCurrentSelectedParticle)
            {
                //
                // Not engaged...
                // ...see if we're able to pick a particle and thus start engagement
                //

                auto const particleId = mLabController->TryPickNpcParticle(mousePosition);
                if (particleId.has_value())
                {
                    //
                    // Engage!
                    //

                    mCurrentSelectedParticle.emplace(*particleId);
                }
            }
            else
            {
                //
                // Engaged
                //

                mLabController->NotifyNpcParticleTrajectory(
                    *mCurrentSelectedParticle,
                    mousePosition);
            }
        }
        else
        {
            if (mCurrentSelectedParticle)
            {
                // Set trajectory
                mLabController->SetNpcParticleTrajectory(
                    *mCurrentSelectedParticle,
                    mousePosition);

                // Disengage
                mCurrentSelectedParticle.reset();
            }
        }
    }

private:

    // Our state

    std::optional<ElementIndex> mCurrentSelectedParticle; // When set, we're setting a trajectory for a particle

    // The cursor
    wxCursor const mCursor;
};

class SetOriginTriangleTool final : public Tool
{
public:

    SetOriginTriangleTool(
        wxWindow * cursorWindow,
        std::shared_ptr<LabController> labController);

public:

    virtual void Initialize(InputState const & /*inputState*/) override
    {
        mIsLeftMouseDown = false;

        // Set cursor
        SetCurrentCursor();
    }

    virtual void Deinitialize(InputState const & /*inputState*/) override
    {
    }

    virtual void SetCurrentCursor() override
    {
        mCursorWindow->SetCursor(mCursor);
    }

    virtual void Update(InputState const & inputState) override
    {
        if (inputState.IsLeftMouseDown)
        {
            if (!mIsLeftMouseDown)
            {
                mLabController->TrySelectOriginTriangle(inputState.MousePosition);

                mIsLeftMouseDown = true;
            }
        }
        else
        {
            mIsLeftMouseDown = false;
        }
    }

private:

    // Our state

    bool mIsLeftMouseDown;

    // The cursor
    wxCursor const mCursor;
};

class SelectParticleTool final : public Tool
{
public:

    SelectParticleTool(
        wxWindow * cursorWindow,
        std::shared_ptr<LabController> labController);

public:

    virtual void Initialize(InputState const & /*inputState*/) override
    {
        mIsLeftMouseDown = false;

        // Set cursor
        SetCurrentCursor();
    }

    virtual void Deinitialize(InputState const & /*inputState*/) override
    {
    }

    virtual void SetCurrentCursor() override
    {
        mCursorWindow->SetCursor(mCursor);
    }

    virtual void Update(InputState const & inputState) override
    {
        if (inputState.IsLeftMouseDown)
        {
            if (!mIsLeftMouseDown)
            {
                mLabController->TrySelectParticle(inputState.MousePosition);

                mIsLeftMouseDown = true;
            }
        }
        else
        {
            mIsLeftMouseDown = false;
        }
    }

private:

    // Our state

    bool mIsLeftMouseDown;

    // The cursor
    wxCursor const mCursor;
};

