/***************************************************************************************
* Original Author:      Gabriele Giuseppini
* Created:              2020-05-21
* Copyright:            Gabriele Giuseppini  (https://github.com/GabrieleGiuseppini)
***************************************************************************************/
#include "ControlToolbar.h"

#include "WxHelpers.h"

#include "UIControls/LinearSliderCore.h"

#include <BLabCoreLib/ResourceLocator.h>

#include <wx/bitmap.h>
#include <wx/gbsizer.h>
#include <wx/sizer.h>

#include <cassert>
#include <memory>
#include <string>

wxEventType const ControlToolbar::wxEVT_TOOLBAR_ACTION = wxNewEventType();

long const ControlToolbar::ID_MOVE_PARTICLE = wxNewId();
long const ControlToolbar::ID_SET_PARTICLE_TRAJECTORY = wxNewId();
long const ControlToolbar::ID_SET_ORIGIN_TRIANGLE = wxNewId();
long const ControlToolbar::ID_MOVE_VERTEX = wxNewId();
long const ControlToolbar::ID_ROTATE_MESH_BY_POSITION = wxNewId();
long const ControlToolbar::ID_ROTATE_MESH_BY_PARTICLE = wxNewId();

long const ControlToolbar::ID_SET_PARTICLE_GRAVITY = wxNewId();

long const ControlToolbar::ID_SIMULATION_CONTROL_PLAY = wxNewId();
long const ControlToolbar::ID_SIMULATION_CONTROL_PAUSE = wxNewId();
long const ControlToolbar::ID_SIMULATION_CONTROL_STEP = wxNewId();

long const ControlToolbar::ID_ACTION_RESET = wxNewId();
long const ControlToolbar::ID_ACTION_LOAD_MESH = wxNewId();
long const ControlToolbar::ID_ACTION_SETTINGS = wxNewId();

long const ControlToolbar::ID_VIEW_CONTROL_GRID = wxNewId();
long const ControlToolbar::ID_RENDER_SIMULATION_STEPS = wxNewId();

wxDEFINE_EVENT(EVT_MESH_TRANSFORMATION_CHANGED, ControlToolbar::meshTransformationChangedEvent);

ControlToolbar::ControlToolbar(wxWindow * parent)
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

                    ReconciliateUIWithTool(ToolType::MoveParticle);
                });

            mMoveParticleButton->SetToolTip("Move a particle");

            gridSizer->Add(mMoveParticleButton);
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

                    ReconciliateUIWithTool(ToolType::SetParticleTrajectory);
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

                    ReconciliateUIWithTool(ToolType::SetOriginTriangle);
                });

            mSetOriginTriangleButton->SetToolTip("Set a triangle as origin of the barycentric coordinates");

            gridSizer->Add(mSetOriginTriangleButton);
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

                    ReconciliateUIWithTool(ToolType::MoveVertex);
                });

            mMoveVertexButton->SetToolTip("Move a mesh vertex");

            gridSizer->Add(mMoveVertexButton);
        }

        // Rotate mesh by position
        {
            mRotateMeshByPositionButton = new wxBitmapToggleButton(
                this,
                ID_ROTATE_MESH_BY_POSITION,
                wxBitmap(
                    (ResourceLocator::GetResourcesFolderPath() / "rotate_mesh_by_position_icon.png").string(),
                    wxBITMAP_TYPE_PNG),
                wxDefaultPosition, wxDefaultSize, wxBU_EXACTFIT);

            mRotateMeshByPositionButton->Bind(wxEVT_TOGGLEBUTTON,
                [this](wxCommandEvent & /*event*/)
                {
                    wxCommandEvent evt(wxEVT_TOOLBAR_ACTION, ID_ROTATE_MESH_BY_POSITION);
                    ProcessEvent(evt);

                    ReconciliateUIWithTool(ToolType::RotateMeshByPosition);
                });

            mMoveVertexButton->SetToolTip("Rotate the mesh around a point");

            gridSizer->Add(mRotateMeshByPositionButton);
        }

        // Rotate mesh by particle
        {
            mRotateMeshByParticleButton = new wxBitmapToggleButton(
                this,
                ID_ROTATE_MESH_BY_PARTICLE,
                wxBitmap(
                    (ResourceLocator::GetResourcesFolderPath() / "rotate_mesh_by_particle_icon.png").string(),
                    wxBITMAP_TYPE_PNG),
                wxDefaultPosition, wxDefaultSize, wxBU_EXACTFIT);

            mRotateMeshByParticleButton->Bind(wxEVT_TOGGLEBUTTON,
                [this](wxCommandEvent & /*event*/)
                {
                    wxCommandEvent evt(wxEVT_TOOLBAR_ACTION, ID_ROTATE_MESH_BY_PARTICLE);
                    ProcessEvent(evt);

                    ReconciliateUIWithTool(ToolType::RotateMeshByParticle);
                });

            mMoveVertexButton->SetToolTip("Rotate the mesh around a particle");

            gridSizer->Add(mRotateMeshByParticleButton);
        }


        vSizer->Add(gridSizer, 0, wxALIGN_CENTER | wxALL, 5);
    }

    vSizer->AddSpacer(10);

    {
        wxGridSizer * gridSizer = new wxGridSizer(2, 2, 2);

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
                    evt.SetInt(mSetParticleGravityButton->GetValue() ? 1 : 0);
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

        // Render simulation steps
        {
            mRenderSimulationStepsButton = new wxBitmapToggleButton(
                this,
                ID_RENDER_SIMULATION_STEPS,
                wxBitmap(
                    (ResourceLocator::GetResourcesFolderPath() / "render_simulation_steps.png").string(),
                    wxBITMAP_TYPE_PNG),
                wxDefaultPosition, wxDefaultSize, wxBU_EXACTFIT);

            mRenderSimulationStepsButton->Bind(wxEVT_TOGGLEBUTTON, [this](wxCommandEvent & /*event*/) { OnViewControlButton(mRenderSimulationStepsButton); });

            mRenderSimulationStepsButton->SetToolTip("Toggles rendering of simulation steps");

            gridSizer->Add(mRenderSimulationStepsButton);
        }

        vSizer->Add(gridSizer, 0, wxALIGN_CENTER | wxALL, 5);
    }

    vSizer->AddSpacer(10);

    // Mesh transformation
    {
        wxBoxSizer * hSizer = new wxBoxSizer(wxHORIZONTAL);

        static int constexpr SliderWidth = 50;
        static int constexpr SliderHeight = 140;

        // Horizontal velocity
        {
            mHorizontalMeshVelocitySlider = new SliderControl<float>(
                this,
                SliderWidth,
                SliderHeight,
                "H Vel",
                "The horizontal velocity of the mesh",
                [this](float value)
                {
                    meshTransformationChangedEvent evt(EVT_MESH_TRANSFORMATION_CHANGED, this->GetId(), 
                        vec2f(value, mVerticalMeshVelocitySlider->GetValue()));
                    ProcessWindowEvent(evt);
                },
                std::make_unique<LinearSliderCore>(
                   -30.0f,
                    30.0f));

            mHorizontalMeshVelocitySlider->SetValue(0.0f);

            hSizer->Add(
                mHorizontalMeshVelocitySlider,
                0, 
                wxALL, 
                1);
        }

        // Vertical velocity
        {
            mVerticalMeshVelocitySlider = new SliderControl<float>(
                this,
                SliderWidth,
                SliderHeight,
                "V Vel",
                "The vertical velocity of the mesh",
                [this](float value)
                {
                    meshTransformationChangedEvent evt(EVT_MESH_TRANSFORMATION_CHANGED, this->GetId(),
                        vec2f(mHorizontalMeshVelocitySlider->GetValue(), value));
                    ProcessWindowEvent(evt);
                },
                std::make_unique<LinearSliderCore>(
                    -30.0f,
                    30.0f));

            mVerticalMeshVelocitySlider->SetValue(0.0f);

            hSizer->Add(
                mVerticalMeshVelocitySlider,
                0,
                wxALL,
                1);
        }

        vSizer->Add(hSizer, 0, wxALIGN_CENTER | wxALL, 5);
    }

    vSizer->AddSpacer(10);

    this->SetSizer(vSizer);
}

void ControlToolbar::ReconcialiteUI(
    bool isGravityEnabled,
    bool isViewGridEnabled,
    bool isRenderSimulationStepsEnabled)
{
    mSetParticleGravityButton->SetValue(isGravityEnabled);
    mViewControlGridButton->SetValue(isViewGridEnabled);
    mRenderSimulationStepsButton->SetValue(isRenderSimulationStepsEnabled);
}

void ControlToolbar::ReconciliateUIWithTool(ToolType tool)
{
    switch (tool)
    {
        case ToolType::MoveParticle:
        {
            mMoveParticleButton->SetValue(true);
            mMoveVertexButton->SetValue(false);
            mRotateMeshByPositionButton->SetValue(false);
            mRotateMeshByParticleButton->SetValue(false);
            mSetParticleTrajectoryButton->SetValue(false);
            mSetOriginTriangleButton->SetValue(false);

            break;
        }

        case ToolType::MoveVertex:
        {
            mMoveParticleButton->SetValue(false);
            mMoveVertexButton->SetValue(true);
            mRotateMeshByPositionButton->SetValue(false);
            mRotateMeshByParticleButton->SetValue(false);
            mSetParticleTrajectoryButton->SetValue(false);
            mSetOriginTriangleButton->SetValue(false);

            break;
        }

        case ToolType::RotateMeshByPosition:
        {
            mMoveParticleButton->SetValue(false);
            mMoveVertexButton->SetValue(false);
            mRotateMeshByPositionButton->SetValue(true);
            mRotateMeshByParticleButton->SetValue(false);
            mSetParticleTrajectoryButton->SetValue(false);
            mSetOriginTriangleButton->SetValue(false);

            break;
        }

        case ToolType::RotateMeshByParticle:
        {
            mMoveParticleButton->SetValue(false);
            mMoveVertexButton->SetValue(false);
            mRotateMeshByPositionButton->SetValue(false);
            mRotateMeshByParticleButton->SetValue(true);
            mSetParticleTrajectoryButton->SetValue(false);
            mSetOriginTriangleButton->SetValue(false);

            break;
        }

        case ToolType::SetParticleTrajectory:
        {
            mMoveParticleButton->SetValue(false);
            mMoveVertexButton->SetValue(false);
            mRotateMeshByPositionButton->SetValue(false);
            mRotateMeshByParticleButton->SetValue(false);
            mSetParticleTrajectoryButton->SetValue(true);
            mSetOriginTriangleButton->SetValue(false);

            break;
        }

        case ToolType::SetOriginTriangle:
        {
            mMoveParticleButton->SetValue(false);
            mMoveVertexButton->SetValue(false);
            mRotateMeshByPositionButton->SetValue(false);
            mRotateMeshByParticleButton->SetValue(false);
            mSetParticleTrajectoryButton->SetValue(false);
            mSetOriginTriangleButton->SetValue(true);

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
    else if (button->GetId() == ID_RENDER_SIMULATION_STEPS)
    {
        // Fire event
        wxCommandEvent evt(wxEVT_TOOLBAR_ACTION, ID_RENDER_SIMULATION_STEPS);
        evt.SetInt(button->GetValue() ? 1 : 0);
        ProcessEvent(evt);
    }
}