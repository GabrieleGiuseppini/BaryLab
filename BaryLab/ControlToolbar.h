/***************************************************************************************
* Original Author:      Gabriele Giuseppini
* Created:              2020-05-21
* Copyright:            Gabriele Giuseppini  (https://github.com/GabrieleGiuseppini)
***************************************************************************************/
#pragma once

#include "Tools.h"

#include "UIControls/SliderControl.h"

#include <GameCore/Vectors.h>

#include <wx/bmpbuttn.h>
#include <wx/event.h>
#include <wx/panel.h>
#include <wx/tglbtn.h>

class ControlToolbar final : public wxPanel
{
public:

    static wxEventType const wxEVT_TOOLBAR_ACTION;

    static long const ID_MOVE_PARTICLE;
    static long const ID_SET_PARTICLE_TRAJECTORY;
    static long const ID_SET_ORIGIN_TRIANGLE;
    static long const ID_SELECT_PARTICLE;
    static long const ID_MOVE_VERTEX;
    static long const ID_ROTATE_MESH_BY_POSITION;
    static long const ID_ROTATE_MESH_BY_PARTICLE;

    static long const ID_ADD_HUMAN_NPC;
    static long const ID_MOVE_NPC;
    static long const ID_REMOVE_NPC;

    static long const ID_SET_PARTICLE_GRAVITY;

    static long const ID_SIMULATION_CONTROL_PLAY;
    static long const ID_SIMULATION_CONTROL_PAUSE;
    static long const ID_SIMULATION_CONTROL_STEP;

    static long const ID_ACTION_RESET;
    static long const ID_ACTION_LOAD_MESH;
    static long const ID_ACTION_SETTINGS;

    static long const ID_VIEW_CONTROL_GRID;

    class meshTransformationChangedEvent : public wxEvent
    {
    public:

        meshTransformationChangedEvent(
            wxEventType eventType,
            int winid,
            vec2f velocity)
            : wxEvent(winid, eventType)
            , mVelocity(velocity)
        {
            m_propagationLevel = wxEVENT_PROPAGATE_MAX;
        }

        meshTransformationChangedEvent(meshTransformationChangedEvent const & other)
            : wxEvent(other)
            , mVelocity(other.mVelocity)
        {
            m_propagationLevel = wxEVENT_PROPAGATE_MAX;
        }

        virtual wxEvent * Clone() const override
        {
            return new meshTransformationChangedEvent(*this);
        }

        vec2f const & GetVelocity() const
        {
            return mVelocity;
        }

    private:

        vec2f const mVelocity;
    };

    class humanNpcPanicLevelChangedEvent : public wxEvent
    {
    public:

        humanNpcPanicLevelChangedEvent(
            wxEventType eventType,
            int winid,
            float panicLevel)
            : wxEvent(winid, eventType)
            , mPanicLevel(panicLevel)
        {
            m_propagationLevel = wxEVENT_PROPAGATE_MAX;
        }

        humanNpcPanicLevelChangedEvent(humanNpcPanicLevelChangedEvent const & other)
            : wxEvent(other)
            , mPanicLevel(other.mPanicLevel)
        {
            m_propagationLevel = wxEVENT_PROPAGATE_MAX;
        }

        virtual wxEvent * Clone() const override
        {
            return new humanNpcPanicLevelChangedEvent(*this);
        }

        float const & GetPanicLevel() const
        {
            return mPanicLevel;
        }

    private:

        float const mPanicLevel;
    };

public:

    ControlToolbar(wxWindow * parent);

    virtual ~ControlToolbar() = default;

    void ReconcialiteUI(
        SimulationControlStateType simulationControlState,
        bool isGravityEnabled,
        bool isViewGridEnabled);

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
    wxBitmapToggleButton * mSetParticleTrajectoryButton;
    wxBitmapToggleButton * mSetOriginTriangleButton;
    wxBitmapToggleButton * mSelectParticleButton;
    wxBitmapToggleButton * mMoveVertexButton;
    wxBitmapToggleButton * mRotateMeshByPositionButton;
    wxBitmapToggleButton * mRotateMeshByParticleButton;

    wxBitmapToggleButton * mAddHumanNpcButton;
    wxBitmapToggleButton * mMoveNpcButton;
    wxBitmapToggleButton * mRemoveNpcButton;

    wxBitmapToggleButton * mSetParticleGravityButton;

    wxBitmapToggleButton * mSimulationControlPlayButton;
    wxBitmapToggleButton * mSimulationControlPauseButton;
    wxBitmapButton * mSimulationControlStepButton;

    wxBitmapToggleButton * mViewControlGridButton;

    SliderControl<float> * mHorizontalMeshVelocitySlider;
    SliderControl<float> * mVerticalMeshVelocitySlider;
    SliderControl<float> * mNpcHumanPanicLevelSlider;
};

wxDECLARE_EVENT(EVT_MESH_TRANSFORMATION_CHANGED, ControlToolbar::meshTransformationChangedEvent);
wxDECLARE_EVENT(EVT_HUMAN_NPC_PANIC_LEVEL_CHANGED, ControlToolbar::humanNpcPanicLevelChangedEvent);
