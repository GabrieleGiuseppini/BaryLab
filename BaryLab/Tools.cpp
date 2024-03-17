/***************************************************************************************
* Original Author:      Gabriele Giuseppini
* Created:              2020-05-16
* Copyright:            Gabriele Giuseppini  (https://github.com/GabrieleGiuseppini)
***************************************************************************************/
#include "Tools.h"

MoveParticleTool::MoveParticleTool(
    wxWindow * cursorWindow,
    std::shared_ptr<LabController> labController)
    : Tool(
        ToolType::MoveParticle,
        cursorWindow,
        std::move(labController))
    , mCurrentEngagementState(std::nullopt)
    , mUpCursor(WxHelpers::MakeCursor("generic_cursor", 7, 8))
    , mDownCursor(WxHelpers::MakeCursor("generic_cursor", 7, 8))
{
}

MoveVertexTool::MoveVertexTool(
    wxWindow * cursorWindow,
    std::shared_ptr<LabController> labController)
    : Tool(
        ToolType::MoveVertex,
        cursorWindow,
        std::move(labController))
    , mCurrentEngagementState(std::nullopt)
    , mUpCursor(WxHelpers::MakeCursor("generic_cursor", 7, 8))
    , mDownCursor(WxHelpers::MakeCursor("generic_cursor", 7, 8))
{
}

RotateMeshByPositionTool::RotateMeshByPositionTool(
    wxWindow * cursorWindow,
    std::shared_ptr<LabController> labController)
    : Tool(
        ToolType::RotateMeshByPosition,
        cursorWindow,
        std::move(labController))
    , mCurrentEngagementState()
    , mUpCursor(WxHelpers::MakeCursor("rotate_mesh_cursor_up", 13, 5))
    , mDownCursor(WxHelpers::MakeCursor("rotate_mesh_cursor_down", 13, 5))
{
}

RotateMeshByParticleTool::RotateMeshByParticleTool(
    wxWindow * cursorWindow,
    std::shared_ptr<LabController> labController)
    : Tool(
        ToolType::RotateMeshByParticle,
        cursorWindow,
        std::move(labController))
    , mCurrentEngagementState()
    , mUpCursor(WxHelpers::MakeCursor("rotate_mesh_cursor_up", 13, 5))
    , mDownCursor(WxHelpers::MakeCursor("rotate_mesh_cursor_down", 13, 5))
{
}

SetParticleTrajectoryTool::SetParticleTrajectoryTool(
    wxWindow * cursorWindow,
    std::shared_ptr<LabController> labController)
    : Tool(
        ToolType::SetParticleTrajectory,
        cursorWindow,
        std::move(labController))
    , mCurrentSelectedParticle()
    , mCursor(WxHelpers::MakeCursor("generic_cursor", 7, 8))
{
}

SetOriginTriangleTool::SetOriginTriangleTool(
    wxWindow * cursorWindow,
    std::shared_ptr<LabController> labController)
    : Tool(
        ToolType::SetOriginTriangle,
        cursorWindow,
        std::move(labController))
    , mIsLeftMouseDown(false)
    , mCursor(WxHelpers::MakeCursor("generic_cursor", 7, 8))
{
}

SelectParticleTool::SelectParticleTool(
    wxWindow * cursorWindow,
    std::shared_ptr<LabController> labController)
    : Tool(
        ToolType::SelectParticle,
        cursorWindow,
        std::move(labController))
    , mIsLeftMouseDown(false)
    , mCursor(WxHelpers::MakeCursor("generic_cursor", 7, 8))
{
}

AddHumanNpcTool::AddHumanNpcTool(
    wxWindow * cursorWindow,
    std::shared_ptr<LabController> labController)
    : AddNpcToolBase(
        ToolType::AddHumanNpc,
        cursorWindow,
        std::move(labController))
{
}

MoveNpcTool::MoveNpcTool(
    wxWindow * cursorWindow,
    std::shared_ptr<LabController> labController)
    : Tool(
        ToolType::MoveNpc,
        cursorWindow,
        std::move(labController))
    , mClosedCursor(WxHelpers::MakeCursor("move_npc_cursor_down", 11, 29))
    , mOpenCursor(WxHelpers::MakeCursor("move_npc_cursor_up", 11, 29))
{
}
