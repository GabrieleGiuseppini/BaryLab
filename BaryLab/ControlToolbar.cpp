/***************************************************************************************
* Original Author:      Gabriele Giuseppini
* Created:              2020-05-21
* Copyright:            Gabriele Giuseppini  (https://github.com/GabrieleGiuseppini)
***************************************************************************************/
#include "ControlToolbar.h"

#include "WxHelpers.h"

#include <BLabCoreLib/ResourceLocator.h>

#include <wx/bitmap.h>
#include <wx/gbsizer.h>
#include <wx/sizer.h>

#include <cassert>

wxEventType const ControlToolbar::wxEVT_TOOLBAR_ACTION = wxNewEventType();

long const ControlToolbar::ID_MOVE_PARTICLE = wxNewId();
long const ControlToolbar::ID_SET_PARTICLE_GRAVITY = wxNewId();
long const ControlToolbar::ID_SET_PARTICLE_TRAJECTORY = wxNewId();
long const ControlToolbar::ID_MOVE_VERTEX = wxNewId();
long const ControlToolbar::ID_SET_ORIGIN_TRIANGLE = wxNewId();

long const ControlToolbar::ID_SIMULATION_CONTROL_PLAY = wxNewId();
long const ControlToolbar::ID_SIMULATION_CONTROL_PAUSE = wxNewId();
long const ControlToolbar::ID_SIMULATION_CONTROL_STEP = wxNewId();

long const ControlToolbar::ID_ACTION_RESET = wxNewId();
long const ControlToolbar::ID_ACTION_LOAD_MESH = wxNewId();
long const ControlToolbar::ID_ACTION_SETTINGS = wxNewId();

long const ControlToolbar::ID_VIEW_CONTROL_GRID = wxNewId();

ControlToolbar::ControlToolbar(wxWindow* parent)
    : wxPanel(
        parent,
        wxID_ANY,
        wxDefaultPosition,
        wxDefaultSize,
        wxBORDER_SIMPLE)
{
    SetBackgroundColour(wxSystemSettings::GetColour(wxSYS_COLOUR_BTNFACE));

    wxBoxSizer * vSizer = new wxBoxSizer(wxVERTICAL);

    {
        wxGridSizer * gridSizer = new wxGridSizer(2, 2, 2);

        // Move particle
        {
            mMoveParticleButton = new wxBitmapToggleButton(
                this,
                ID_MOVE_PARTICLE,
                wxBitmap(
                    // TODO
                    (ResourceLocator::GetResourcesFolderPath() / "move_particle_icon.png").string(),
                    wxBITMAP_TYPE_PNG),
                wxDefaultPosition, wxDefaultSize, wxBU_EXACTFIT);

            mMoveParticleButton->Bind(wxEVT_TOGGLEBUTTON, 
                [this](wxCommandEvent & /*event*/) 
                { 
                    wxCommandEvent evt(wxEVT_TOOLBAR_ACTION, ID_MOVE_PARTICLE);
                    ProcessEvent(evt);
                });

            mMoveParticleButton->SetToolTip("Move a particle");

            gridSizer->Add(mMoveParticleButton);
        }

        // Move vertex
        {
            mMoveVertexButton = new wxBitmapToggleButton(
                this,
                ID_MOVE_VERTEX,
                wxBitmap(
                (ResourceLocator::GetResourcesFolderPath() / "move_vertex_icon.png").string(),
                wxBITMAP_TYPE_PNG),
                wxDefaultPosition, wxDefaultSize, wxBU_EXACTFIT);

            mMoveVertexButton->Bind(wxEVT_TOGGLEBUTTON,
                [this](wxCommandEvent & /*event*/)
                {
                    wxCommandEvent evt(wxEVT_TOOLBAR_ACTION, ID_MOVE_VERTEX);
                    ProcessEvent(evt);
                });

            mMoveVertexButton->SetToolTip("Move a mesh vertex");

            gridSizer->Add(mMoveVertexButton);
        }

        // Set particle trajectory
        {
            mSetParticleTrajectoryButton = new wxBitmapToggleButton(
                this,
                ID_SET_PARTICLE_TRAJECTORY,
                wxBitmap(
                    (ResourceLocator::GetResourcesFolderPath() / "set_particle_trajectory_icon.png").string(),
                    wxBITMAP_TYPE_PNG),
                wxDefaultPosition, wxDefaultSize, wxBU_EXACTFIT);

            mSetParticleTrajectoryButton->Bind(wxEVT_TOGGLEBUTTON,
                [this](wxCommandEvent & /*event*/)
                {
                    wxCommandEvent evt(wxEVT_TOOLBAR_ACTION, ID_SET_PARTICLE_TRAJECTORY);
                    ProcessEvent(evt);
                });

            mSetParticleTrajectoryButton->SetToolTip("Set a trajectory for a particle");

            gridSizer->Add(mSetParticleTrajectoryButton);
        }

        // Set origin triangle
        {
            mSetOriginTriangleButton = new wxBitmapToggleButton(
                this,
                ID_SET_ORIGIN_TRIANGLE,
                wxBitmap(
                    (ResourceLocator::GetResourcesFolderPath() / "set_origin_triangle_icon.png").string(),
                    wxBITMAP_TYPE_PNG),
                wxDefaultPosition, wxDefaultSize, wxBU_EXACTFIT);

            mSetOriginTriangleButton->Bind(wxEVT_TOGGLEBUTTON,
                [this](wxCommandEvent & /*event*/)
                {
                    wxCommandEvent evt(wxEVT_TOOLBAR_ACTION, ID_SET_ORIGIN_TRIANGLE);
                    ProcessEvent(evt);
                });

            mSetOriginTriangleButton->SetToolTip("Set a triangle as origin of the barycentric coordinates");

            gridSizer->Add(mSetOriginTriangleButton);
        }

        // Set particle gravity
        {
            mSetParticleGravityButton = new wxBitmapToggleButton(
                this,
                ID_SET_PARTICLE_GRAVITY,
                wxBitmap(
                    (ResourceLocator::GetResourcesFolderPath() / "gravity_icon.png").string(),
                    wxBITMAP_TYPE_PNG),
                wxDefaultPosition, wxDefaultSize, wxBU_EXACTFIT);

            mSetParticleGravityButton->Bind(wxEVT_TOGGLEBUTTON,
                [this](wxCommandEvent & /*event*/)
                {
                    wxCommandEvent evt(wxEVT_TOOLBAR_ACTION, ID_SET_PARTICLE_GRAVITY);
                    ProcessEvent(evt);
                });

            mSetParticleGravityButton->SetToolTip("Enable gravity for the particles");

            gridSizer->Add(mSetParticleGravityButton);
        }
        
        vSizer->Add(gridSizer, 0, wxALIGN_CENTER | wxALL, 5);
    }

    vSizer->AddSpacer(10);

    // Simulation control
    {
        wxGridSizer * gridSizer = new wxGridSizer(2, 2, 2);

        // Play
        {
            mSimulationControlPlayButton = new wxBitmapToggleButton(
                this,
                ID_SIMULATION_CONTROL_PLAY,
                wxBitmap(
                    (ResourceLocator::GetResourcesFolderPath() / "simcontrol_play.png").string(),
                    wxBITMAP_TYPE_PNG),
                wxDefaultPosition, wxDefaultSize, wxBU_EXACTFIT);

            mSimulationControlPlayButton->Bind(wxEVT_TOGGLEBUTTON, [this](wxCommandEvent & /*event*/) { OnSimulationControlButton(mSimulationControlPlayButton); });

            mSimulationControlPlayButton->SetToolTip("Start simulation auto-play");

            gridSizer->Add(mSimulationControlPlayButton);
        }

        // Pause
        {
            mSimulationControlPauseButton = new wxBitmapToggleButton(
                this,
                ID_SIMULATION_CONTROL_PAUSE,
                wxBitmap(
                    (ResourceLocator::GetResourcesFolderPath() / "simcontrol_pause.png").string(),
                    wxBITMAP_TYPE_PNG),
                wxDefaultPosition, wxDefaultSize, wxBU_EXACTFIT);

            mSimulationControlPauseButton->SetValue(true); // We start in pause

            mSimulationControlPauseButton->Bind(wxEVT_TOGGLEBUTTON, [this](wxCommandEvent & /*event*/) { OnSimulationControlButton(mSimulationControlPauseButton); });

            mSimulationControlPauseButton->SetToolTip("Pause simulation auto-play (SPACE)");

            gridSizer->Add(mSimulationControlPauseButton);
        }

        // Step
        {
            mSimulationControlStepButton = new wxBitmapButton(
                this,
                ID_SIMULATION_CONTROL_STEP,
                wxBitmap(
                    (ResourceLocator::GetResourcesFolderPath() / "simcontrol_step.png").string(),
                    wxBITMAP_TYPE_PNG),
                wxDefaultPosition, wxDefaultSize, wxBU_EXACTFIT);

            mSimulationControlStepButton->Enable(true); // We start in pause

            mSimulationControlStepButton->Bind(wxEVT_BUTTON, [this](wxCommandEvent & /*event*/) { OnSimulationControlStepButton(); });

            mSimulationControlStepButton->SetToolTip("Run a single simulation step (ENTER)");

            gridSizer->Add(mSimulationControlStepButton);
        }

        vSizer->Add(gridSizer, 0, wxALIGN_CENTER | wxALL, 5);
    }

    vSizer->AddSpacer(10);

    // Action
    {
        wxGridSizer * gridSizer = new wxGridSizer(2, 2, 2);

        // Reset
        {
            auto button = new wxButton(
                this,
                ID_ACTION_RESET,
                wxEmptyString, wxDefaultPosition, wxDefaultSize, wxBU_EXACTFIT);

            button->SetBitmap(
                wxBitmap(
                    (ResourceLocator::GetResourcesFolderPath() / "reset_icon.png").string(),
                    wxBITMAP_TYPE_PNG));

            button->Bind(
                wxEVT_BUTTON,
                [this](wxCommandEvent & /*event*/)
                {
                    wxCommandEvent evt(wxEVT_TOOLBAR_ACTION, ID_ACTION_RESET);
                    ProcessEvent(evt);
                });

            button->SetToolTip("Reset the simulation");

            gridSizer->Add(button);
        }

        // Load
        {
            auto button = new wxButton(
                this,
                ID_ACTION_LOAD_MESH,
                wxEmptyString, wxDefaultPosition, wxDefaultSize, wxBU_EXACTFIT);

            button->SetBitmap(
                wxBitmap(
                    (ResourceLocator::GetResourcesFolderPath() / "open_icon.png").string(),
                    wxBITMAP_TYPE_PNG));

            button->Bind(
                wxEVT_BUTTON,
                [this](wxCommandEvent & /*event*/)
                {
                    wxCommandEvent evt(wxEVT_TOOLBAR_ACTION, ID_ACTION_LOAD_MESH);
                    ProcessEvent(evt);
                });

            button->SetToolTip("Load a new object");

            gridSizer->Add(button);
        }

        // Settings
        {
            auto button = new wxButton(
                this,
                ID_ACTION_SETTINGS,
                wxEmptyString, wxDefaultPosition, wxDefaultSize, wxBU_EXACTFIT);

            button->SetBitmap(
                wxBitmap(
                    (ResourceLocator::GetResourcesFolderPath() / "settings_icon.png").string(),
                    wxBITMAP_TYPE_PNG));

            button->Bind(
                wxEVT_BUTTON,
                [this](wxCommandEvent & /*event*/)
                {
                    wxCommandEvent evt(wxEVT_TOOLBAR_ACTION, ID_ACTION_SETTINGS);
                    ProcessEvent(evt);
                });

            button->SetToolTip("Adjust the simulation's settings");

            gridSizer->Add(button);
        }

        vSizer->Add(gridSizer, 0, wxALIGN_CENTER | wxALL, 5);
    }

    vSizer->AddSpacer(10);

    // View Control
    {
        wxGridSizer * gridSizer = new wxGridSizer(2, 2, 2);

        // Grid
        {
            mViewControlGridButton = new wxBitmapToggleButton(
                this,
                ID_VIEW_CONTROL_GRID,
                wxBitmap(
                    (ResourceLocator::GetResourcesFolderPath() / "view_grid.png").string(),
                    wxBITMAP_TYPE_PNG),
                wxDefaultPosition, wxDefaultSize, wxBU_EXACTFIT);

            mViewControlGridButton->Bind(wxEVT_TOGGLEBUTTON, [this](wxCommandEvent & /*event*/) { OnViewControlButton(mViewControlGridButton); });

            mViewControlGridButton->SetToolTip("Enable or disable grid");

            gridSizer->Add(mViewControlGridButton);
        }

        vSizer->Add(gridSizer, 0, wxALIGN_CENTER | wxALL, 5);
    }

    vSizer->AddSpacer(10);

    this->SetSizer(vSizer);
}

void ControlToolbar::SetTool(ToolType tool)
{
    switch (tool)
    {
        case ToolType::MoveParticle:
        {
            mMoveParticleButton->SetValue(true);
            mMoveVertexButton->SetValue(false);
            mSetParticleTrajectoryButton->SetValue(false);
            mSetOriginTriangleButton->SetValue(false);

            break;
        }

        case ToolType::MoveVertex:
        {
            mMoveParticleButton->SetValue(false);
            mMoveVertexButton->SetValue(true);
            mSetParticleTrajectoryButton->SetValue(false);
            mSetOriginTriangleButton->SetValue(false);

            break;
        }

        case ToolType::SetOriginTriangle:
        {
            mMoveParticleButton->SetValue(false);
            mMoveVertexButton->SetValue(false);
            mSetParticleTrajectoryButton->SetValue(false);
            mSetOriginTriangleButton->SetValue(true);

            break;
        }

        case ToolType::SetParticleTrajectory:
        {
            mMoveParticleButton->SetValue(false);
            mMoveVertexButton->SetValue(false);
            mSetParticleTrajectoryButton->SetValue(true);
            mSetOriginTriangleButton->SetValue(false);

            break;
        }
    }
}

bool ControlToolbar::ProcessKeyDown(
    int keyCode,
    int keyModifiers)
{
    if (keyCode == WXK_SPACE)
    {
        // Pause
        if (!mSimulationControlPauseButton->GetValue())
        {
            mSimulationControlPauseButton->SetFocus();
            mSimulationControlPauseButton->SetValue(true);
            OnSimulationControlButton(mSimulationControlPauseButton);
            return true;
        }
    }
    else if (keyCode == WXK_RETURN)
    {
        // Step
        if (mSimulationControlPauseButton->GetValue())
        {
            OnSimulationControlStepButton();
            return true;
        }
    }
    else if (keyCode == 'R' && keyModifiers == wxMOD_CONTROL)
    {
        // Reset
        OnActionResetButton();
        return true;
    }
    else if (keyCode == 'O' && keyModifiers == wxMOD_CONTROL)
    {
        // Load
        OnActionLoadMeshButton();
        return true;
    }
    else if (keyCode == 'S' && keyModifiers == wxMOD_CONTROL)
    {
        // Settings
        OnActionSettingsButton();
        return true;
    }

    return false;
}

void ControlToolbar::OnSimulationControlButton(wxBitmapToggleButton * button)
{
    if (button->GetId() == ID_SIMULATION_CONTROL_PLAY)
    {
        if (button->GetValue())
        {
            // Set all others to off
            mSimulationControlPauseButton->SetValue(false);

            // Disable step
            mSimulationControlStepButton->Enable(false);

            // Fire event
            wxCommandEvent evt(wxEVT_TOOLBAR_ACTION, ID_SIMULATION_CONTROL_PLAY);
            ProcessEvent(evt);
        }
        else
        {
            // Set back to on
            mSimulationControlPlayButton->SetValue(true);
        }
    }
    else
    {
        assert(button->GetId() == ID_SIMULATION_CONTROL_PAUSE);

        if (button->GetValue())
        {
            // Set all others to off
            mSimulationControlPlayButton->SetValue(false);

            // Enable step
            mSimulationControlStepButton->Enable(true);

            // Fire event
            wxCommandEvent evt(wxEVT_TOOLBAR_ACTION, ID_SIMULATION_CONTROL_PAUSE);
            ProcessEvent(evt);
        }
        else
        {
            // Set back to on
            mSimulationControlPauseButton->SetValue(true);
        }
    }
}

void ControlToolbar::OnSimulationControlStepButton()
{
    // Fire event
    wxCommandEvent evt(wxEVT_TOOLBAR_ACTION, ID_SIMULATION_CONTROL_STEP);
    ProcessEvent(evt);
}

void ControlToolbar::OnActionResetButton()
{
    wxCommandEvent evt(wxEVT_TOOLBAR_ACTION, ID_ACTION_RESET);
    ProcessEvent(evt);
}

void ControlToolbar::OnActionLoadMeshButton()
{
    wxCommandEvent evt(wxEVT_TOOLBAR_ACTION, ID_ACTION_LOAD_MESH);
    ProcessEvent(evt);
}

void ControlToolbar::OnActionSettingsButton()
{
    wxCommandEvent evt(wxEVT_TOOLBAR_ACTION, ID_ACTION_SETTINGS);
    ProcessEvent(evt);
}

void ControlToolbar::OnViewControlButton(wxBitmapToggleButton * button)
{
    if (button->GetId() == ID_VIEW_CONTROL_GRID)
    {
        // Fire event
        wxCommandEvent evt(wxEVT_TOOLBAR_ACTION, ID_VIEW_CONTROL_GRID);
        evt.SetInt(button->GetValue() ? 1 : 0);
        ProcessEvent(evt);
    }
}