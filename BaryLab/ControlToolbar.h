/***************************************************************************************
* Original Author:      Gabriele Giuseppini
* Created:              2020-05-21
* Copyright:            Gabriele Giuseppini  (https://github.com/GabrieleGiuseppini)
***************************************************************************************/
#pragma once

#include "Tools.h"

#include <wx/bmpbuttn.h>
#include <wx/panel.h>
#include <wx/tglbtn.h>

#include <memory>
#include <string>

class ControlToolbar final : public wxPanel
{
public:

    static wxEventType const wxEVT_TOOLBAR_ACTION;

    static long const ID_MOVE_PARTICLE;    
    static long const ID_SET_PARTICLE_GRAVITY;
    static long const ID_SET_PARTICLE_TRAJECTORY;
    static long const ID_MOVE_VERTEX;
    static long const ID_SET_ORIGIN_TRIANGLE;
    
    static long const ID_SIMULATION_CONTROL_PLAY;
    static long const ID_SIMULATION_CONTROL_PAUSE;
    static long const ID_SIMULATION_CONTROL_STEP;

    static long const ID_ACTION_RESET;
    static long const ID_ACTION_LOAD_MESH;
    static long const ID_ACTION_SETTINGS;

    static long const ID_VIEW_CONTROL_GRID;

public:

    ControlToolbar(wxWindow * parent);

    virtual ~ControlToolbar() = default;

    void ReconcialiteUI(bool isGravityEnabled);

    void ReconciliateUIWithTool(ToolType tool);

    bool ProcessKeyDown(
        int keyCode,
        int keyModifiers);

private:

    void OnSimulationControlButton(wxBitmapToggleButton * button);
    void OnSimulationControlStepButton();
    void OnActionResetButton();
    void OnActionLoadMeshButton();
    void OnActionSettingsButton();
    void OnViewControlButton(wxBitmapToggleButton * button);

private:

    wxBitmapToggleButton * mMoveParticleButton;
    wxBitmapToggleButton * mMoveVertexButton;
    wxBitmapToggleButton * mSetParticleTrajectoryButton;
    wxBitmapToggleButton * mSetOriginTriangleButton;
    wxBitmapToggleButton * mSetParticleGravityButton;

    wxBitmapToggleButton * mSimulationControlPlayButton;
    wxBitmapToggleButton * mSimulationControlPauseButton;
    wxBitmapButton * mSimulationControlStepButton;

    wxBitmapToggleButton * mViewControlGridButton;
};
