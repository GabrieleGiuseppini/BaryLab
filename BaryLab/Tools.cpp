/***************************************************************************************
* Original Author:      Gabriele Giuseppini
* Created:              2020-05-16
* Copyright:            Gabriele Giuseppini  (https://github.com/GabrieleGiuseppini)
***************************************************************************************/
#include "Tools.h"

#include "WxHelpers.h"

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
