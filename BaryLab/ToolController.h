/***************************************************************************************
* Original Author:      Gabriele Giuseppini
* Created:              2020-05-16
* Copyright:            Gabriele Giuseppini  (https://github.com/GabrieleGiuseppini)
***************************************************************************************/
#pragma once

#include "Tools.h"

#include <Game/LabController.h>

#include <wx/frame.h>

#include <cassert>
#include <memory>
#include <vector>

class ToolController
{
public:

    ToolController(
        ToolType initialToolType,
        wxWindow * cursorWindow,
        std::shared_ptr<LabController> labController);

    ToolType GetTool() const
    {
        if (nullptr != mCurrentTool)
            return mCurrentTool->GetToolType();
        else
            return ToolType::AddHumanNpc;
    }

    void SetTool(ToolType toolType)
    {
        assert(static_cast<size_t>(toolType) < mAllTools.size());

        // Notify old tool
        if (nullptr != mCurrentTool)
            mCurrentTool->Deinitialize(mInputState);

        // Switch tool
        mCurrentTool = mAllTools[static_cast<size_t>(toolType)].get();
        mCurrentTool->Initialize(mInputState);
    }

    void SetHumanNpcPlaceTool(NpcSubKindIdType subKind)
    {
        AddHumanNpcTool * tool = dynamic_cast<AddHumanNpcTool *>(mAllTools[static_cast<size_t>(ToolType::AddHumanNpc)].get());
        tool->SetHumanNpcKind(subKind);
        SetTool(ToolType::AddHumanNpc);
    }

    void UnsetTool()
    {
        // Notify old tool
        if (nullptr != mCurrentTool)
        {
            mCurrentTool->Deinitialize(mInputState);
        }

        mCurrentTool = nullptr;
    }

    void Update()
    {
        if (nullptr != mCurrentTool)
        {
            mCurrentTool->Update(mInputState);
        }
    }

    //
    // Getters
    //

    vec2f const & GetMouseScreenCoordinates() const
    {
        return mInputState.MousePosition;
    }

    bool IsShiftKeyDown() const
    {
        return mInputState.IsShiftKeyDown;
    }

    //
    // External event handlers
    //

    void OnMouseMove(
        int x,
        int y);

    void OnLeftMouseDown();

    void OnLeftMouseUp();

    void OnRightMouseDown();

    void OnRightMouseUp();

    void OnShiftKeyDown();

    void OnShiftKeyUp();

private:

    // Input state
    InputState mInputState;

    // Tool state
    Tool * mCurrentTool;
    std::vector<std::unique_ptr<Tool>> mAllTools; // Indexed by enum

private:

    wxWindow * const mCursorWindow;
    wxCursor mPanCursor;
    std::shared_ptr<LabController> const mLabController;
};
