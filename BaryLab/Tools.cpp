/***************************************************************************************
* Original Author:      Gabriele Giuseppini
* Created:              2020-05-16
* Copyright:            Gabriele Giuseppini  (https://github.com/GabrieleGiuseppini)
***************************************************************************************/
#include "Tools.h"

#include "WxHelpers.h"

MoveVertexTool::MoveVertexTool(
    wxWindow * cursorWindow,
    std::shared_ptr<LabController> labController)
    : Tool(
        ToolType::MoveVertex,
        cursorWindow,
        std::move(labController))
    , mCurrentEngagementState(std::nullopt)
    , mUpCursor(WxHelpers::MakeCursor("move_vertex_up", 13, 5))
    , mDownCursor(WxHelpers::MakeCursor("move_vertex_down", 13, 5))
{
}
