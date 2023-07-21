/***************************************************************************************
* Original Author:      Gabriele Giuseppini
* Created:              2020-05-16
* Copyright:            Gabriele Giuseppini  (https://github.com/GabrieleGiuseppini)
***************************************************************************************/
#pragma once

#include <BLabCoreLib/LabController.h>

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
    SetParticleTrajectory = 2,
    SetOriginTriangle = 3,
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

        if (inputState.IsLeftMouseDown)
        {
            if (!mCurrentEngagementState)
            {
                //
                // Not engaged...
                // ...see if we're able to pick a point and thus start engagement
                //

                vec2f const mousePosition = inputState.MousePosition;

                auto const vertexId = mLabController->TryPickVertex(mousePosition);
                if (vertexId.has_value())
                {
                    //
                    // Engage!
                    //

                    mCurrentEngagementState.emplace(*vertexId);
                }
            }
            else
            {
                //
                // Engaged
                //

                mLabController->MoveVertexTo(
                    mCurrentEngagementState->VertexIndex,
                    inputState.MousePosition);
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

        explicit EngagementState(ElementIndex vertexIndex)
            : VertexIndex(vertexIndex)
        {}
    };

    std::optional<EngagementState> mCurrentEngagementState; // When set, indicates it's engaged

    // The cursors
    wxCursor const mUpCursor;
    wxCursor const mDownCursor;
};

