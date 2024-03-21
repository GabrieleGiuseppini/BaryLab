/***************************************************************************************
* Original Author:      Gabriele Giuseppini
* Created:              2020-05-16
* Copyright:            Gabriele Giuseppini  (https://github.com/GabrieleGiuseppini)
***************************************************************************************/
#pragma once

#include "WxHelpers.h"

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
    SelectParticle = 6,

    AddHumanNpc = 7,
    MoveNpc = 8,
    RemoveNpc = 9
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

//////////////////////////////////////////////////////////////////////////////////////////////////////////////

class AddNpcToolBase : public Tool
{
    //
    // State machine:
    //  - Clean:
    //      - ButtonDown: BeginPlaceNewHumanNpc; -> Error or -> Engaged
    //  - Error:
    //      - ButtonUp: -> Clean
    //  - Engaged:
    //      - MouseMove: Move;
    //      - ButtonUp: Confirm; -> Clean
    //      - Reset: Abort;
    //

public:

    AddNpcToolBase(
        ToolType toolType,
        wxWindow * cursorWindow,
        std::shared_ptr<LabController> labController)
        : Tool(
            toolType,
            cursorWindow,
            std::move(labController))
        , mCurrentEngagementState(std::nullopt)
        , mClosedCursor(WxHelpers::MakeCursor("move_npc_cursor_down", 11, 29))
        , mOpenCursor(WxHelpers::MakeCursor("move_npc_cursor_up", 11, 29))
        , mErrorCursor(WxHelpers::MakeCursor("move_npc_cursor_error", 11, 29))
    {}

public:

    virtual void Initialize(InputState const & /*inputState*/) override
    {
        mCurrentEngagementState.reset();

        // Set cursor
        SetCurrentCursor();
    }

    virtual void Deinitialize(InputState const & /*inputState*/) override
    {
        if (mCurrentEngagementState.has_value()
            && mCurrentEngagementState->Npc.has_value())
        {
            // Abort
            mLabController->AbortNewNpc(mCurrentEngagementState->Npc->ObjectId);
        }
    }

    virtual void SetCurrentCursor() override
    {
        if (mCurrentEngagementState.has_value())
        {
            if (mCurrentEngagementState->Npc.has_value())
            {
                mCursorWindow->SetCursor(mClosedCursor);
            }
            else
            {
                mCursorWindow->SetCursor(mErrorCursor);
            }
        }
        else
        {
            mCursorWindow->SetCursor(mOpenCursor);
        }
    }

    virtual void Update(InputState const & inputState) override
    {
        if (!mCurrentEngagementState.has_value())
        {
            // Clean

            if (inputState.IsLeftMouseDown)
            {
                // -> Error or -> Engaged
                mCurrentEngagementState.emplace(InternalBeginPlaceNewNpc(inputState.MousePosition));
                SetCurrentCursor();
            }
        }
        else if (!mCurrentEngagementState->Npc.has_value())
        {
            // Error

            if (!inputState.IsLeftMouseDown)
            {
                // -> Clean
                mCurrentEngagementState.reset();
                SetCurrentCursor();
            }
        }
        else
        {
            // Engaged

            if (!inputState.IsLeftMouseDown)
            {
                // Confirm;
                mLabController->CompleteNewNpc(mCurrentEngagementState->Npc->ObjectId);

                // -> Clean
                mCurrentEngagementState.reset();
                SetCurrentCursor();
            }
            else
            {
                // Move;
                mLabController->MoveNpcTo(
                    mCurrentEngagementState->Npc->ObjectId,
                    inputState.MousePosition,
                    mCurrentEngagementState->Npc->WorldOffset);
            }
        }
    }

protected:

    virtual std::optional<PickedObjectId<NpcId>> InternalBeginPlaceNewNpc(vec2f const & screenCoordinates) = 0;

private:

    // Our state

    struct EngagementState
    {
        std::optional<PickedObjectId<NpcId>> Npc; // When not set, we're in error mode

        explicit EngagementState(std::optional<PickedObjectId<NpcId>> npc)
            : Npc(npc)
        {}
    };

    std::optional<EngagementState> mCurrentEngagementState; // When set, indicates it's engaged

    // The cursors
    wxCursor const mClosedCursor;
    wxCursor const mOpenCursor;
    wxCursor const mErrorCursor;
};

class AddHumanNpcTool final : public AddNpcToolBase
{
public:

    AddHumanNpcTool(
        wxWindow * cursorWindow,
        std::shared_ptr<LabController> labController);

    void SetHumanNpcKind(HumanNpcKindType humanNpcKind)
    {
        mHumanNpcKind = humanNpcKind;
    }

protected:

    std::optional<PickedObjectId<NpcId>> InternalBeginPlaceNewNpc(vec2f const & screenCoordinates) override
    {
        return mLabController->BeginPlaceNewHumanNpc(
            mHumanNpcKind,
            screenCoordinates);
    }

private:

    HumanNpcKindType mHumanNpcKind;
};

class MoveNpcTool final : public Tool
{
    //
    // State machine:
    //
    //  - Hovering (MouseUp): if we have Npc that's a candidate - and it's highlighted
    //  - Moving (MouseDown): if we have Npc we're moving it - and it's not highlighted
    //

public:

    MoveNpcTool(
        wxWindow * cursorWindow,
        std::shared_ptr<LabController> labController);

public:

    virtual void Initialize(InputState const & inputState) override
    {
        mNpc.reset();
        mIsMouseDown = inputState.IsLeftMouseDown;

        // Set cursor
        SetCurrentCursor();
    }

    virtual void Deinitialize(InputState const & /*inputState*/) override
    {
        if (mNpc.has_value() && !mIsMouseDown)
        {
            mLabController->HighlightNpc(mNpc->ObjectId, NpcHighlightType::None);
        }
    }

    virtual void SetCurrentCursor() override
    {
        if (mIsMouseDown)
        {
            mCursorWindow->SetCursor(mClosedCursor);
        }
        else
        {
            mCursorWindow->SetCursor(mOpenCursor);
        }
    }

    virtual void Update(InputState const & inputState) override
    {
        if (inputState.IsLeftMouseDown)
        {
            if (!mIsMouseDown)
            {
                // Clicked

                // If we have an NPC, it becomes the one we're moving
                if (mNpc.has_value())
                {
                    mLabController->BeginMoveNpc(mNpc->ObjectId);

                    // Now that it's moving, un-highlight it
                    mLabController->HighlightNpc(
                        mNpc->ObjectId,
                        NpcHighlightType::None);
                }

                mIsMouseDown = true;
                SetCurrentCursor();
            }
            else
            {
                // Moved while keeping down

                if (mNpc.has_value())
                {
                    mLabController->MoveNpcTo(
                        mNpc->ObjectId,
                        inputState.MousePosition,
                        mNpc->WorldOffset);
                }
            }
        }
        else
        {
            if (mIsMouseDown)
            {
                // Released

                // If we have an NPC, we've stopped moving it
                if (mNpc.has_value())
                {
                    mLabController->EndMoveNpc(mNpc->ObjectId);
                }

                mIsMouseDown = false;
                SetCurrentCursor();
            }

            // Probe at new position
            auto const probeOutcome = mLabController->ProbeNpcAt(inputState.MousePosition);
            if (probeOutcome)
            {
                if (mNpc.has_value() && mNpc->ObjectId != probeOutcome->ObjectId)
                {
                    mLabController->HighlightNpc(
                        mNpc->ObjectId,
                        NpcHighlightType::None);
                }

                if (!mNpc.has_value() || mNpc->ObjectId != probeOutcome->ObjectId)
                {
                    mLabController->HighlightNpc(
                        probeOutcome->ObjectId,
                        NpcHighlightType::Candidate);
                }
            }
            else
            {
                if (mNpc.has_value())
                {
                    mLabController->HighlightNpc(
                        mNpc->ObjectId,
                        NpcHighlightType::None);
                }
            }

            mNpc = probeOutcome; // Always update so to pick latest offset
        }
    }

private:

    // Our state
    std::optional<PickedObjectId<NpcId>> mNpc;
    bool mIsMouseDown;

    // The cursors
    wxCursor const mClosedCursor;
    wxCursor const mOpenCursor;
};

class RemoveNpcTool final : public Tool
{
    //
    // State machine:
    //
    //  - Hovering (MouseUp): if we have Npc that's a candidate - and it's highlighted
    //  - Moving (MouseDown): if we have Npc we're removing it - and it's not highlighted
    //

public:

    RemoveNpcTool(
        wxWindow * cursorWindow,
        std::shared_ptr<LabController> labController);

public:

    virtual void Initialize(InputState const & inputState) override
    {
        mNpc.reset();
        mIsMouseDown = inputState.IsLeftMouseDown;

        // Set cursor
        SetCurrentCursor();
    }

    virtual void Deinitialize(InputState const & /*inputState*/) override
    {
        if (mNpc.has_value() && !mIsMouseDown)
        {
            mLabController->HighlightNpc(*mNpc, NpcHighlightType::None);
        }
    }

    virtual void SetCurrentCursor() override
    {
        if (mIsMouseDown)
        {
            mCursorWindow->SetCursor(mClosedCursor);
        }
        else
        {
            mCursorWindow->SetCursor(mOpenCursor);
        }
    }

    virtual void Update(InputState const & inputState) override
    {
        if (inputState.IsLeftMouseDown)
        {
            if (!mIsMouseDown)
            {
                // Clicked

                // If we have an NPC, it becomes the one we're removing
                if (mNpc.has_value())
                {
                    mLabController->RemoveNpc(*mNpc);

                    mNpc = std::nullopt;
                }

                mIsMouseDown = true;
                SetCurrentCursor();
            }
            else
            {
                // Moved while keeping down

                // Nop
            }
        }
        else
        {
            if (mIsMouseDown)
            {
                // Released

                mIsMouseDown = false;
                SetCurrentCursor();
            }

            // Probe at new position
            auto const probeOutcome = mLabController->ProbeNpcAt(inputState.MousePosition);
            if (probeOutcome)
            {
                if (mNpc.has_value() && *mNpc != probeOutcome->ObjectId)
                {
                    mLabController->HighlightNpc(
                        *mNpc,
                        NpcHighlightType::None);
                }

                if (!mNpc.has_value() || *mNpc != probeOutcome->ObjectId)
                {
                    mLabController->HighlightNpc(
                        probeOutcome->ObjectId,
                        NpcHighlightType::Candidate);

                    mNpc = probeOutcome->ObjectId;
                }
            }
            else
            {
                if (mNpc.has_value())
                {
                    mLabController->HighlightNpc(
                        *mNpc,
                        NpcHighlightType::None);

                    mNpc = std::nullopt;
                }
            }
        }
    }

private:

    // Our state
    std::optional<NpcId> mNpc;
    bool mIsMouseDown;

    // The cursors
    wxCursor const mClosedCursor;
    wxCursor const mOpenCursor;
};

