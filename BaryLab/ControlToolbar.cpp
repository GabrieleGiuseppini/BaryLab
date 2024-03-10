/***************************************************************************************
* Original Author:      Gabriele Giuseppini
* Created:              2020-05-21
* Copyright:            Gabriele Giuseppini  (https://github.com/GabrieleGiuseppini)
***************************************************************************************/
#include "ControlToolbar.h"

#include "WxHelpers.h"

#include "UIControls/LinearSliderCore.h"

#include <Game/ResourceLocator.h>

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
long const ControlToolbar::ID_SELECT_PARTICLE = wxNewId();
long const ControlToolbar::ID_MOVE_VERTEX = wxNewId();
long const ControlToolbar::ID_ROTATE_MESH_BY_POSITION = wxNewId();
long const ControlToolbar::ID_ROTATE_MESH_BY_PARTICLE = wxNewId();

long const ControlToolbar::ID_ADD_HUMAN_NPC = wxNewId();
long const ControlToolbar::ID_MOVE_NPC = wxNewId();
long const ControlToolbar::ID_REMOVE_NPC = wxNewId();

long const ControlToolbar::ID_SET_PARTICLE_GRAVITY = wxNewId();

long const ControlToolbar::ID_SIMULATION_CONTROL_PLAY = wxNewId();
long const ControlToolbar::ID_SIMULATION_CONTROL_PAUSE = wxNewId();
long const ControlToolbar::ID_SIMULATION_CONTROL_STEP = wxNewId();

long const ControlToolbar::ID_ACTION_RESET = wxNewId();
long const ControlToolbar::ID_ACTION_LOAD_MESH = wxNewId();
long const ControlToolbar::ID_ACTION_SETTINGS = wxNewId();

long const ControlToolbar::ID_VIEW_CONTROL_GRID = wxNewId();

wxDEFINE_EVENT(EVT_MESH_TRANSFORMATION_CHANGED, ControlToolbar::meshTransformationChangedEvent);
wxDEFINE_EVENT(EVT_HUMAN_NPC_PANIC_LEVEL_CHANGED, ControlToolbar::humanNpcPanicLevelChangedEvent);

ControlToolbar::ControlToolbar(wxWindow * parent)
    : wxPanel(
        parent,
        wxID_ANY,
        wxDefaultPosition,
        wxDefaultSize,
        wxBORDER_SIMPLE)
{
    int constexpr SliderWidth = 50;
    int constexpr SliderHeight = 140;

    SetBackgroundColour(wxSystemSettings::GetColour(wxSYS_COLOUR_BTNFACE));

    wxBoxSizer * vSizer = new wxBoxSizer(wxVERTICAL);

    {
        wxGridSizer * gridSizer = new wxGridSizer(3, 2, 2);

        // Select particle
        {
            mSelectParticleButton = new wxBitmapToggleButton(
                this,
                ID_SELECT_PARTICLE,
                wxBitmap(
                    (ResourceLocator::GetResourcesFolderPath() / "select_particle_icon.png").string(),
                    wxBITMAP_TYPE_PNG),
                wxDefaultPosition, wxDefaultSize, wxBU_EXACTFIT);

            mSelectParticleButton->Bind(wxEVT_TOGGLEBUTTON,
                [this](wxCommandEvent & /*event*/)
                {
                    wxCommandEvent evt(wxEVT_TOOLBAR_ACTION, ID_SELECT_PARTICLE);
                    ProcessEvent(evt);

                    ReconciliateUIWithTool(ToolType::SelectParticle);
                });

            mSelectParticleButton->SetToolTip("Select a particle");

            gridSizer->Add(mSelectParticleButton);
        }

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

        /////////////////

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

        // Spacer
        {
            gridSizer->AddSpacer(0);
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

            mRotateMeshByPositionButton->SetToolTip("Rotate the mesh around a point");

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

            mRotateMeshByParticleButton->SetToolTip("Rotate the mesh around a particle");

            gridSizer->Add(mRotateMeshByParticleButton);
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
                    evt.SetInt(mSetParticleGravityButton->GetValue() ? 1 : 0);
                    ProcessEvent(evt);
                });

            mSetParticleGravityButton->SetToolTip("Enable gravity for the particles");

            gridSizer->Add(mSetParticleGravityButton);
        }

        vSizer->Add(gridSizer, 0, wxALIGN_CENTER | wxALL, 5);
    }

    vSizer->AddSpacer(10);

    // Tools
    {
        wxGridSizer * gridSizer = new wxGridSizer(3, 2, 2);

        // Add Human NPC
        {
            mAddHumanNpcButton = new wxBitmapToggleButton(
                this,
                ID_ADD_HUMAN_NPC,
                wxBitmap(
                    (ResourceLocator::GetResourcesFolderPath() / "add_human_npc_icon.png").string(),
                    wxBITMAP_TYPE_PNG),
                wxDefaultPosition, wxDefaultSize, wxBU_EXACTFIT);

            mAddHumanNpcButton->Bind(wxEVT_TOGGLEBUTTON,
                [this](wxCommandEvent & /*event*/)
                {
                    wxCommandEvent evt(wxEVT_TOOLBAR_ACTION, ID_ADD_HUMAN_NPC);
                    ProcessEvent(evt);

                    ReconciliateUIWithTool(ToolType::AddHumanNpc);
                });

            mAddHumanNpcButton->SetToolTip("Add a human NPC");

            gridSizer->Add(mAddHumanNpcButton);
        }

        // Move NPC
        {
            mMoveNpcButton = new wxBitmapToggleButton(
                this,
                ID_MOVE_NPC,
                wxBitmap(
                    (ResourceLocator::GetResourcesFolderPath() / "move_npc_icon.png").string(),
                    wxBITMAP_TYPE_PNG),
                wxDefaultPosition, wxDefaultSize, wxBU_EXACTFIT);

            mMoveNpcButton->Bind(wxEVT_TOGGLEBUTTON,
                [this](wxCommandEvent & /*event*/)
                {
                    wxCommandEvent evt(wxEVT_TOOLBAR_ACTION, ID_MOVE_NPC);
                    ProcessEvent(evt);

                    ReconciliateUIWithTool(ToolType::MoveNpc);
                });

            mMoveNpcButton->SetToolTip("Move an NPC");

            gridSizer->Add(mMoveNpcButton);
        }

        // Remove NPC
        {
            mRemoveNpcButton = new wxBitmapToggleButton(
                this,
                ID_REMOVE_NPC,
                wxBitmap(
                    (ResourceLocator::GetResourcesFolderPath() / "remove_npc_icon.png").string(),
                    wxBITMAP_TYPE_PNG),
                wxDefaultPosition, wxDefaultSize, wxBU_EXACTFIT);

            mRemoveNpcButton->Bind(wxEVT_TOGGLEBUTTON,
                [this](wxCommandEvent & /*event*/)
                {
                    wxCommandEvent evt(wxEVT_TOOLBAR_ACTION, ID_REMOVE_NPC);
                    ProcessEvent(evt);

                    ReconciliateUIWithTool(ToolType::RemoveNpc);
                });

            mRemoveNpcButton->SetToolTip("Remove an NPC");

            gridSizer->Add(mRemoveNpcButton);
        }

        vSizer->Add(gridSizer, 0, wxALIGN_CENTER | wxALL, 5);
    }

    vSizer->AddSpacer(10);

    // Simulation control
    {
        wxGridSizer * gridSizer = new wxGridSizer(3, 2, 2);

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
        wxGridSizer * gridSizer = new wxGridSizer(3, 2, 2);

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
        wxGridSizer * gridSizer = new wxGridSizer(1, 2, 2);

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

    // Mesh transformation
    {
        auto * gridSizer = new wxFlexGridSizer(2, 2, 2);

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

            gridSizer->Add(
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

            gridSizer->Add(
                mVerticalMeshVelocitySlider,
                0,
                wxALL,
                1);
        }

        // Horizontal velocity zero
        {
            auto * button = new wxButton(
                this,
                wxID_ANY,
                "Zero",
                wxDefaultPosition,
                wxSize(-1, -1));

            button->Bind(wxEVT_BUTTON,
                [this](wxCommandEvent & /*event*/)
                {
                    meshTransformationChangedEvent evt(EVT_MESH_TRANSFORMATION_CHANGED, this->GetId(),
                        vec2f(0.0f, mVerticalMeshVelocitySlider->GetValue()));
                    ProcessEvent(evt);

                    mHorizontalMeshVelocitySlider->SetValue(0.0f);
                });

            gridSizer->Add(
                button,
                0,
                wxALL,
                1);
        }

        // Vertical velocity zero
        {
            auto * button = new wxButton(
                this,
                wxID_ANY,
                "Zero",
                wxDefaultPosition,
                wxSize(-1, -1));

            button->Bind(wxEVT_BUTTON,
                [this](wxCommandEvent & /*event*/)
                {
                    meshTransformationChangedEvent evt(EVT_MESH_TRANSFORMATION_CHANGED, this->GetId(),
                        vec2f(mHorizontalMeshVelocitySlider->GetValue(), 0.0f));
                    ProcessEvent(evt);

                    mVerticalMeshVelocitySlider->SetValue(0.0f);
                });

            gridSizer->Add(
                button,
                0,
                wxALL,
                1);
        }

        vSizer->Add(gridSizer, 0, wxALIGN_CENTER | wxALL, 5);
    }

    vSizer->AddSpacer(10);

    // NPC control
    {
        auto * gridSizer = new wxFlexGridSizer(2, 2, 2);

        // Human panic level
        {
            mNpcHumanPanicLevelSlider = new SliderControl<float>(
                this,
                SliderWidth,
                SliderHeight,
                "Panic",
                "The panic level for all humans",
                [this](float value)
                {
                    humanNpcPanicLevelChangedEvent evt(EVT_HUMAN_NPC_PANIC_LEVEL_CHANGED, this->GetId(),
                        value);
                    ProcessWindowEvent(evt);
                },
                std::make_unique<LinearSliderCore>(
                    0.0f,
                    4.0f));

            mNpcHumanPanicLevelSlider->SetValue(0.0f);

            gridSizer->Add(
                mNpcHumanPanicLevelSlider,
                1,
                wxALL,
                1);
        }

        vSizer->Add(gridSizer, 0, wxALIGN_CENTER | wxALL, 5);
    }

    vSizer->AddSpacer(10);

    this->SetSizer(vSizer);
}

void ControlToolbar::ReconcialiteUI(
    bool isGravityEnabled,
    bool isViewGridEnabled)
{
    mSetParticleGravityButton->SetValue(isGravityEnabled);
    mViewControlGridButton->SetValue(isViewGridEnabled);
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
            mSelectParticleButton->SetValue(false);

            mAddHumanNpcButton->SetValue(false);
            mMoveNpcButton->SetValue(false);
            mRemoveNpcButton->SetValue(false);

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
            mSelectParticleButton->SetValue(false);

            mAddHumanNpcButton->SetValue(false);
            mMoveNpcButton->SetValue(false);
            mRemoveNpcButton->SetValue(false);

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
            mSelectParticleButton->SetValue(false);

            mAddHumanNpcButton->SetValue(false);
            mMoveNpcButton->SetValue(false);
            mRemoveNpcButton->SetValue(false);

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
            mSelectParticleButton->SetValue(false);

            mAddHumanNpcButton->SetValue(false);
            mMoveNpcButton->SetValue(false);
            mRemoveNpcButton->SetValue(false);

            break;
        }

        case ToolType::SelectParticle:
        {
            mMoveParticleButton->SetValue(false);
            mMoveVertexButton->SetValue(false);
            mRotateMeshByPositionButton->SetValue(false);
            mRotateMeshByParticleButton->SetValue(false);
            mSetParticleTrajectoryButton->SetValue(false);
            mSetOriginTriangleButton->SetValue(false);
            mSelectParticleButton->SetValue(true);

            mAddHumanNpcButton->SetValue(false);
            mMoveNpcButton->SetValue(false);
            mRemoveNpcButton->SetValue(false);

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
            mSelectParticleButton->SetValue(false);

            mAddHumanNpcButton->SetValue(false);
            mMoveNpcButton->SetValue(false);
            mRemoveNpcButton->SetValue(false);

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
            mSelectParticleButton->SetValue(false);

            mAddHumanNpcButton->SetValue(false);
            mMoveNpcButton->SetValue(false);
            mRemoveNpcButton->SetValue(false);

            break;
        }

        case ToolType::AddHumanNpc:
        {
            mMoveParticleButton->SetValue(false);
            mMoveVertexButton->SetValue(false);
            mRotateMeshByPositionButton->SetValue(false);
            mRotateMeshByParticleButton->SetValue(false);
            mSetParticleTrajectoryButton->SetValue(false);
            mSetOriginTriangleButton->SetValue(false);
            mSelectParticleButton->SetValue(false);

            mAddHumanNpcButton->SetValue(true);
            mMoveNpcButton->SetValue(false);
            mRemoveNpcButton->SetValue(false);

            break;
        }

        case ToolType::MoveNpc:
        {
            mMoveParticleButton->SetValue(false);
            mMoveVertexButton->SetValue(false);
            mRotateMeshByPositionButton->SetValue(false);
            mRotateMeshByParticleButton->SetValue(false);
            mSetParticleTrajectoryButton->SetValue(false);
            mSetOriginTriangleButton->SetValue(false);
            mSelectParticleButton->SetValue(false);

            mAddHumanNpcButton->SetValue(false);
            mMoveNpcButton->SetValue(true);
            mRemoveNpcButton->SetValue(false);

            break;
        }

        case ToolType::RemoveNpc:
        {
            mMoveParticleButton->SetValue(false);
            mMoveVertexButton->SetValue(false);
            mRotateMeshByPositionButton->SetValue(false);
            mRotateMeshByParticleButton->SetValue(false);
            mSetParticleTrajectoryButton->SetValue(false);
            mSetOriginTriangleButton->SetValue(false);
            mSelectParticleButton->SetValue(false);

            mAddHumanNpcButton->SetValue(false);
            mMoveNpcButton->SetValue(false);
            mRemoveNpcButton->SetValue(true);

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
        if (!mSimulationControlPauseButton->GetValue())
        {
            // Pause
            mSimulationControlPauseButton->SetFocus();
            mSimulationControlPauseButton->SetValue(true);
            OnSimulationControlButton(mSimulationControlPauseButton);
            return true;
        }
        else
        {
            // Step
            OnSimulationControlStepButton();
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