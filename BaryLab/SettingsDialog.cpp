/***************************************************************************************
 * Original Author:     Gabriele Giuseppini
 * Created:             2020-05-23
 * Copyright:           Gabriele Giuseppini  (https://github.com/GabrieleGiuseppini)
 ***************************************************************************************/
#include "SettingsDialog.h"

#include "WxHelpers.h"

#include "UIControls/ExponentialSliderCore.h"
#include "UIControls/FixedTickSliderCore.h"
#include "UIControls/IntegralLinearSliderCore.h"
#include "UIControls/LinearSliderCore.h"

#include <wx/gbsizer.h>
#include <wx/intl.h>
#include <wx/notebook.h>
#include <wx/sizer.h>
#include <wx/stattext.h>
#include <wx/string.h>

#include <algorithm>
#include <stdexcept>

#ifdef _MSC_VER
 // Nothing to do here - we use RC files
#else
#include "Resources/SLabBBB.xpm"
#endif

static int constexpr SliderWidth = 80;
static int constexpr SliderHeight = 140;

static int constexpr StaticBoxTopMargin = 7;
static int constexpr StaticBoxInsetMargin = 10;
static int constexpr CellBorder = 8;

SettingsDialog::SettingsDialog(
    wxWindow * parent,
    std::shared_ptr<SettingsManager> settingsManager,
    std::shared_ptr<LabController> labsimulationController)
    : mParent(parent)
    , mSettingsManager(std::move(settingsManager))
	, mLabController(std::move(labsimulationController))
    // State
    , mLiveSettings(mSettingsManager->MakeSettings())
    , mCheckpointSettings(mSettingsManager->MakeSettings())
{
    Create(
        mParent,
        wxID_ANY,
        _("Simulation Settings"),
        wxDefaultPosition,
        wxSize(400, 200),
        wxCAPTION | wxCLOSE_BOX | wxMINIMIZE_BOX | wxFRAME_NO_TASKBAR
			| /* wxFRAME_FLOAT_ON_PARENT */ wxSTAY_ON_TOP, // See https://trac.wxwidgets.org/ticket/18535
        _T("Settings Window"));

    this->Bind(wxEVT_CLOSE_WINDOW, &SettingsDialog::OnCloseButton, this);

    SetBackgroundColour(wxSystemSettings::GetColour(wxSYS_COLOUR_BTNFACE));

    SetIcon(wxICON(BBB_SLAB_ICON));


    //
    // Lay the dialog out
    //

    wxBoxSizer * dialogVSizer = new wxBoxSizer(wxVERTICAL);


    wxNotebook * notebook = new wxNotebook(
        this,
        wxID_ANY,
        wxPoint(-1, -1),
        wxSize(-1, -1),
        wxNB_TOP);


    //
    // Simulator
    //

    wxPanel * simulatorPanel = new wxPanel(notebook);

    PopulateSimulatorPanel(simulatorPanel);

    notebook->AddPage(simulatorPanel, "Simulator");

    //
    // NPCs
    //

    wxPanel * npcsPanel = new wxPanel(notebook);

    PopulateNpcsPanel(npcsPanel);

    notebook->AddPage(npcsPanel, "NPCs");


    dialogVSizer->Add(notebook, 0, wxEXPAND);


    dialogVSizer->AddSpacer(20);



    // Buttons
    {
        wxBoxSizer * buttonsSizer = new wxBoxSizer(wxHORIZONTAL);

        buttonsSizer->AddSpacer(20);

		mRevertToDefaultsButton = new wxButton(this, wxID_ANY, "Revert to Defaults");
		mRevertToDefaultsButton->SetToolTip("Resets all settings to their default values.");
		mRevertToDefaultsButton->Bind(wxEVT_BUTTON, &SettingsDialog::OnRevertToDefaultsButton, this);
		buttonsSizer->Add(mRevertToDefaultsButton, 0, 0, 0);

		buttonsSizer->AddStretchSpacer(1);

        mOkButton = new wxButton(this, wxID_ANY, "OK");
		mOkButton->SetToolTip("Closes the window keeping all changes.");
        mOkButton->Bind(wxEVT_BUTTON, &SettingsDialog::OnOkButton, this);
        buttonsSizer->Add(mOkButton, 0, 0, 0);

        buttonsSizer->AddSpacer(20);

        mCancelButton = new wxButton(this, wxID_ANY, "Cancel");
		mCancelButton->SetToolTip("Reverts all changes effected since the window was last opened, and closes the window.");
        mCancelButton->Bind(wxEVT_BUTTON, &SettingsDialog::OnCancelButton, this);
        buttonsSizer->Add(mCancelButton, 0, 0, 0);

        buttonsSizer->AddSpacer(20);

        mUndoButton = new wxButton(this, wxID_ANY, "Undo");
		mUndoButton->SetToolTip("Reverts all changes effected since the window was last opened.");
        mUndoButton->Bind(wxEVT_BUTTON, &SettingsDialog::OnUndoButton, this);
        buttonsSizer->Add(mUndoButton, 0, 0, 0);

        buttonsSizer->AddSpacer(20);

        dialogVSizer->Add(buttonsSizer, 0, wxEXPAND, 0);
    }

    dialogVSizer->AddSpacer(20);


    //
    // Finalize dialog
    //

    SetSizerAndFit(dialogVSizer);

    Centre(wxCENTER_ON_SCREEN | wxBOTH);
}

SettingsDialog::~SettingsDialog()
{
}

void SettingsDialog::Open()
{
    if (IsShown())
        return; // Handle Ctrl^S while minimized

    assert(!!mSettingsManager);

    //
    // Initialize state
    //

    // Pull currently-enforced settings
    mSettingsManager->Pull(mLiveSettings);
    mLiveSettings.ClearAllDirty();

    // Save checkpoint for undo
    mCheckpointSettings = mLiveSettings;

    // Populate controls with live settings
	SyncControlsWithSettings(mLiveSettings);

    // Remember that the user hasn't changed anything yet in this session
    mHasBeenDirtyInCurrentSession = false;

	// Enable Revert to Defaults button only if settings are different than defaults
	mAreSettingsDirtyWrtDefaults = (mLiveSettings != mSettingsManager->GetDefaults());

	// Reconcile controls wrt dirty state
	ReconcileDirtyState();


    //
    // Open dialog
    //

    this->Raise();
    this->Show();
}

////////////////////////////////////////////////////////////

void SettingsDialog::OnRevertToDefaultsButton(wxCommandEvent& /*event*/)
{
	//
	// Enforce default settings
	//

	mLiveSettings = mSettingsManager->GetDefaults();

	// Do not update checkpoint, allow user to revert to it

	// Enforce everything as a safety net, immediately
	mLiveSettings.MarkAllAsDirty();
	mSettingsManager->EnforceDirtySettingsImmediate(mLiveSettings);

	// We are back in sync
	mLiveSettings.ClearAllDirty();

	assert(mSettingsManager->Pull() == mLiveSettings);

	// Re-populate controls with new values
	SyncControlsWithSettings(mLiveSettings);

	// Remember user has made changes wrt checkpoint
	mHasBeenDirtyInCurrentSession = true;

	// Remember we are clean now wrt defaults
	mAreSettingsDirtyWrtDefaults = false;

	ReconcileDirtyState();
}

void SettingsDialog::OnOkButton(wxCommandEvent & /*event*/)
{
    // Just close the dialog
    DoClose();
}

void SettingsDialog::OnCancelButton(wxCommandEvent & /*event*/)
{
    DoCancel();
}

void SettingsDialog::OnUndoButton(wxCommandEvent & /*event*/)
{
    assert(!!mSettingsManager);

    //
    // Undo changes done since last open, including eventual loads
    //

    mLiveSettings = mCheckpointSettings;

    // Just enforce anything in the checkpoint that is different than the current settings,
	// immediately
    mLiveSettings.SetDirtyWithDiff(mSettingsManager->Pull());
    mSettingsManager->EnforceDirtySettingsImmediate(mLiveSettings);

    mLiveSettings.ClearAllDirty();

    assert(mSettingsManager->Pull() == mCheckpointSettings);

    // Re-populate controls with new values
	SyncControlsWithSettings(mLiveSettings);

    // Remember we are clean now
    mHasBeenDirtyInCurrentSession = false;
    ReconcileDirtyState();
}

void SettingsDialog::OnCloseButton(wxCloseEvent & /*event*/)
{
    DoCancel();
}

/////////////////////////////////////////////////////////////////////////////

void SettingsDialog::DoCancel()
{
    assert(!!mSettingsManager);

    if (mHasBeenDirtyInCurrentSession)
    {
        //
        // Undo changes done since last open, including eventual loads
        //

        mLiveSettings = mCheckpointSettings;

        // Just enforce anything in the checkpoint that is different than the current settings,
		// immediately
        mLiveSettings.SetDirtyWithDiff(mSettingsManager->Pull());
        mSettingsManager->EnforceDirtySettingsImmediate(mLiveSettings);
    }

    //
    // Close the dialog
    //

    DoClose();
}

void SettingsDialog::DoClose()
{
    this->Hide();
}

void SettingsDialog::PopulateSimulatorPanel(wxPanel * panel)
{
    wxGridBagSizer * gridSizer = new wxGridBagSizer(0, 0);

    // Mechanics
    {
        wxStaticBox * mechanicsBox = new wxStaticBox(panel, wxID_ANY, _("Mechanics"));

        wxBoxSizer * mechanicsBoxSizer = new wxBoxSizer(wxVERTICAL);
        mechanicsBoxSizer->AddSpacer(StaticBoxTopMargin);

        {
            wxGridBagSizer * mechanicsSizer = new wxGridBagSizer(0, 0);

            // Elasticity Adjustment
            {
                mElasticityAdjustmentSlider = new SliderControl<float>(
                    mechanicsBox,
                    SliderWidth,
                    SliderHeight,
                    "Elasticity Adjust",
                    "Adjusts the elasticity of NPC particles' impacts.",
					[this](float value)
                    {
                        this->mLiveSettings.SetValue(SLabSettings::ElasticityAdjustment, value);
                        this->OnLiveSettingsChanged();
                    },
                    std::make_unique<LinearSliderCore>(
                        mLabController->GetMinElasticityAdjustment(),
                        mLabController->GetMaxElasticityAdjustment()));

                mechanicsSizer->Add(
                    mElasticityAdjustmentSlider,
                    wxGBPosition(0, 0),
                    wxGBSpan(1, 1),
                    wxEXPAND | wxALL,
                    CellBorder);
            }

            // Static Friction Adjustment
            {
                mStaticFrictionAdjustmentSlider = new SliderControl<float>(
                    mechanicsBox,
                    SliderWidth,
                    SliderHeight,
                    "Static Friction Adjust",
                    "The adjustment for the static friction of NPC particles against floor surfaces.",
					[this](float value)
                    {
                        this->mLiveSettings.SetValue(SLabSettings::StaticFrictionAdjustment, value);
                        this->OnLiveSettingsChanged();
                    },
                    std::make_unique<ExponentialSliderCore>(
                        mLabController->GetMinStaticFrictionAdjustment(),
                        1.0f,
                        mLabController->GetMaxStaticFrictionAdjustment()));

                mechanicsSizer->Add(
                    mStaticFrictionAdjustmentSlider,
                    wxGBPosition(0, 1),
                    wxGBSpan(1, 1),
                    wxEXPAND | wxALL,
                    CellBorder);
            }

            // Kinetic Friction Adjustment
            {
                mKineticFrictionAdjustmentSlider = new SliderControl<float>(
                    mechanicsBox,
                    SliderWidth,
                    SliderHeight,
                    "Kinetic Friction Adjust",
                    "The adjustment for the kinetic friction of NPC particles against floor surfaces.",
                    [this](float value)
                    {
                        this->mLiveSettings.SetValue(SLabSettings::KineticFrictionAdjustment, value);
                        this->OnLiveSettingsChanged();
                    },
                    std::make_unique<ExponentialSliderCore>(
                        mLabController->GetMinKineticFrictionAdjustment(),
                        1.0f,
                        mLabController->GetMaxKineticFrictionAdjustment()));

                mechanicsSizer->Add(
                    mKineticFrictionAdjustmentSlider,
                    wxGBPosition(0, 2),
                    wxGBSpan(1, 1),
                    wxEXPAND | wxALL,
                    CellBorder);
            }

            // Mass Adjustment
            {
                mMassAdjustmentSlider = new SliderControl<float>(
                    mechanicsBox,
                    SliderWidth,
                    SliderHeight,
                    "Mass Adjust",
                    "Adjusts the mass of all particles.",
					[this](float value)
                    {
                        this->mLiveSettings.SetValue(SLabSettings::MassAdjustment, value);
                        this->OnLiveSettingsChanged();
                    },
                    std::make_unique<ExponentialSliderCore>(
                        mLabController->GetMinMassAdjustment(),
                        1.0f,
                        mLabController->GetMaxMassAdjustment()));

                mechanicsSizer->Add(
                    mMassAdjustmentSlider,
                    wxGBPosition(0, 3),
                    wxGBSpan(1, 1),
                    wxEXPAND | wxALL,
                    CellBorder);
            }

            // Gravity Adjustment
            {
                mGravityAdjustmentSlider = new SliderControl<float>(
                    mechanicsBox,
                    SliderWidth,
                    SliderHeight,
                    "Gravity Adjust",
                    "Adjusts the magnitude of gravity.",
                    [this](float value)
                    {
                        this->mLiveSettings.SetValue(SLabSettings::GravityAdjustment, value);
                        this->OnLiveSettingsChanged();
                    },
                    std::make_unique<ExponentialSliderCore>(
                        mLabController->GetMinGravityAdjustment(),
                        1.0f,
                        mLabController->GetMaxGravityAdjustment()));

                mechanicsSizer->Add(
                    mGravityAdjustmentSlider,
                    wxGBPosition(0, 4),
                    wxGBSpan(1, 1),
                    wxEXPAND | wxALL,
                    CellBorder);
            }

            // Global Damping
            {
                mGlobalDampingAdjustmentSlider = new SliderControl<float>(
                    mechanicsBox,
                    SliderWidth,
                    SliderHeight,
                    "Global Damp Adjust",
                    "Adjusts the global velocity damping.",
                    [this](float value)
                    {
                        this->mLiveSettings.SetValue(SLabSettings::GlobalDampingAdjustment, value);
                        this->OnLiveSettingsChanged();
                    },
                    std::make_unique<LinearSliderCore>(
                        mLabController->GetMinGlobalDampingAdjustment(),
                        mLabController->GetMaxGlobalDampingAdjustment()));

                mechanicsSizer->Add(
                    mGlobalDampingAdjustmentSlider,
                    wxGBPosition(0, 5),
                    wxGBSpan(1, 1),
                    wxEXPAND | wxALL,
                    CellBorder);
            }

            // Sea Level
            {
                mSeaLevelSlider = new SliderControl<float>(
                    mechanicsBox,
                    SliderWidth,
                    SliderHeight,
                    "Sea Level",
                    "The level of the sea (m).",
                    [this](float value)
                    {
                        this->mLiveSettings.SetValue(SLabSettings::SeaLevel, value);
                        this->OnLiveSettingsChanged();
                    },
                    std::make_unique<LinearSliderCore>(
                        mLabController->GetMinSeaLevel(),
                        mLabController->GetMaxSeaLevel()));

                mechanicsSizer->Add(
                    mSeaLevelSlider,
                    wxGBPosition(0, 6),
                    wxGBSpan(1, 1),
                    wxEXPAND | wxALL,
                    CellBorder);
            }

            // Spring Reduction Fraction Adjustment
            {
                mSpringReductionFractionAdjustmentSlider = new SliderControl<float>(
                    mechanicsBox,
                    SliderWidth,
                    SliderHeight,
                    "Spring Stiffness Adjust",
                    "Adjusts the stiffness of springs.",
                    [this](float value)
                    {
                        this->mLiveSettings.SetValue(SLabSettings::SpringReductionFractionAdjustment, value);
                        this->OnLiveSettingsChanged();
                    },
                    std::make_unique<LinearSliderCore>(
                        mLabController->GetMinSpringReductionFractionAdjustment(),
                        mLabController->GetMaxSpringReductionFractionAdjustment()));

                mechanicsSizer->Add(
                    mSpringReductionFractionAdjustmentSlider,
                    wxGBPosition(0, 7),
                    wxGBSpan(1, 1),
                    wxEXPAND | wxALL,
                    CellBorder);
            }

            // Spring Damping Coefficient Adjustment
            {
                mSpringDampingCoefficientAdjustmentSlider = new SliderControl<float>(
                    mechanicsBox,
                    SliderWidth,
                    SliderHeight,
                    "Spring Damping Adjust",
                    "Adjusts the damping coefficient of springs.",
                    [this](float value)
                    {
                        this->mLiveSettings.SetValue(SLabSettings::SpringDampingCoefficientAdjustment, value);
                        this->OnLiveSettingsChanged();
                    },
                    std::make_unique<LinearSliderCore>(
                        mLabController->GetMinSpringDampingCoefficientAdjustment(),
                        mLabController->GetMaxSpringDampingCoefficientAdjustment()));

                mechanicsSizer->Add(
                    mSpringDampingCoefficientAdjustmentSlider,
                    wxGBPosition(0, 8),
                    wxGBSpan(1, 1),
                    wxEXPAND | wxALL,
                    CellBorder);
            }

            /////////////////////////////

            // Water Friction Drag Coefficient Adjustment
            {
                mWaterFrictionDragCoefficientAdjustmentSlider = new SliderControl<float>(
                    mechanicsBox,
                    SliderWidth,
                    SliderHeight,
                    "Water Drag Adjust",
                    "Adjusts the water drag.",
                    [this](float value)
                    {
                        this->mLiveSettings.SetValue(SLabSettings::WaterFrictionDragAdjustment, value);
                        this->OnLiveSettingsChanged();
                    },
                    std::make_unique<LinearSliderCore>(
                        mLabController->GetMinWaterFrictionDragAdjustment(),
                        mLabController->GetMaxWaterFrictionDragAdjustment()));

                mechanicsSizer->Add(
                    mWaterFrictionDragCoefficientAdjustmentSlider,
                    wxGBPosition(1, 0),
                    wxGBSpan(1, 1),
                    wxEXPAND | wxALL,
                    CellBorder);
            }

            // Buoyancy Adjustment
            {
                mBuoyancyAdjustmentSlider = new SliderControl<float>(
                    mechanicsBox,
                    SliderWidth,
                    SliderHeight,
                    "Buoyancy Adjust",
                    "Adjusts the buoyancy.",
                    [this](float value)
                    {
                        this->mLiveSettings.SetValue(SLabSettings::BuoyancyAdjustment, value);
                        this->OnLiveSettingsChanged();
                    },
                    std::make_unique<LinearSliderCore>(
                        mLabController->GetMinBuoyancyAdjustment(),
                        mLabController->GetMaxBuoyancyAdjustment()));

                mechanicsSizer->Add(
                    mBuoyancyAdjustmentSlider,
                    wxGBPosition(1, 1),
                    wxGBSpan(1, 1),
                    wxEXPAND | wxALL,
                    CellBorder);
            }


            mechanicsBoxSizer->Add(mechanicsSizer, 0, wxALL, StaticBoxInsetMargin);
        }

        mechanicsBox->SetSizerAndFit(mechanicsBoxSizer);

        gridSizer->Add(
            mechanicsBox,
            wxGBPosition(0, 0),
            wxGBSpan(1, 8),
            wxEXPAND | wxALL | wxALIGN_CENTER_HORIZONTAL,
            CellBorder);
    }

    // Finalize panel

    for (int c = 0; c < gridSizer->GetCols(); ++c)
        gridSizer->AddGrowableCol(c);

    panel->SetSizer(gridSizer);
}

void SettingsDialog::PopulateNpcsPanel(wxPanel * panel)
{
    wxGridBagSizer * gridSizer = new wxGridBagSizer(0, 0);

    // NPCs
    {
        wxStaticBox * npcsBox = new wxStaticBox(panel, wxID_ANY, _("NPCs"));

        wxBoxSizer * npcsBoxSizer = new wxBoxSizer(wxVERTICAL);
        npcsBoxSizer->AddSpacer(StaticBoxTopMargin);

        {
            wxGridBagSizer * npcsSizer = new wxGridBagSizer(0, 0);

            // Human NPC Equilibrium Torque Stiffness Coefficient
            {
                mHumanNpcEquilibriumTorqueStiffnessCoefficientSlider = new SliderControl<float>(
                    npcsBox,
                    SliderWidth,
                    SliderHeight,
                    "Torque Stiffness",
                    "The strength of the torque applied while rising.",
                    [this](float value)
                    {
                        this->mLiveSettings.SetValue(SLabSettings::HumanNpcEquilibriumTorqueStiffnessCoefficient, value);
                        this->OnLiveSettingsChanged();
                    },
                    std::make_unique<LinearSliderCore>(
                        mLabController->GetMinHumanNpcEquilibriumTorqueStiffnessCoefficient(),
                        mLabController->GetMaxHumanNpcEquilibriumTorqueStiffnessCoefficient()));

                npcsSizer->Add(
                    mHumanNpcEquilibriumTorqueStiffnessCoefficientSlider,
                    wxGBPosition(0, 0),
                    wxGBSpan(1, 1),
                    wxEXPAND | wxALL,
                    CellBorder);
            }

            // Human NPC Equilibrium Torque Damping Coefficient
            {
                mHumanNpcEquilibriumTorqueDampingCoefficientSlider = new SliderControl<float>(
                    npcsBox,
                    SliderWidth,
                    SliderHeight,
                    "Torque Damping",
                    "The damping of the torque applied while rising.",
                    [this](float value)
                    {
                        this->mLiveSettings.SetValue(SLabSettings::HumanNpcEquilibriumTorqueDampingCoefficient, value);
                        this->OnLiveSettingsChanged();
                    },
                    std::make_unique<LinearSliderCore>(
                        mLabController->GetMinHumanNpcEquilibriumTorqueDampingCoefficient(),
                        mLabController->GetMaxHumanNpcEquilibriumTorqueDampingCoefficient()));

                npcsSizer->Add(
                    mHumanNpcEquilibriumTorqueDampingCoefficientSlider,
                    wxGBPosition(0, 1),
                    wxGBSpan(1, 1),
                    wxEXPAND | wxALL,
                    CellBorder);
            }

            // Human NPC Walking Speed Adjustment
            {
                mHumanNpcWalkingSpeedAdjustmentSlider = new SliderControl<float>(
                    npcsBox,
                    SliderWidth,
                    SliderHeight,
                    "Walking Speed Adjust",
                    "Adjusts the speed at which human NPC's walk.",
                    [this](float value)
                    {
                        this->mLiveSettings.SetValue(SLabSettings::HumanNpcWalkingSpeedAdjustment, value);
                        this->OnLiveSettingsChanged();
                    },
                    std::make_unique<LinearSliderCore>(
                        mLabController->GetMinHumanNpcWalkingSpeedAdjustment(),
                        mLabController->GetMaxHumanNpcWalkingSpeedAdjustment()));

                npcsSizer->Add(
                    mHumanNpcWalkingSpeedAdjustmentSlider,
                    wxGBPosition(0, 2),
                    wxGBSpan(1, 1),
                    wxEXPAND | wxALL,
                    CellBorder);
            }

            // NPC Size Multiplier
            {
                mNpcSizeMultiplierSlider = new SliderControl<float>(
                    npcsBox,
                    SliderWidth,
                    SliderHeight,
                    "NPC Size Multiplier",
                    "Adjusts the size of NPCs.",
                    [this](float value)
                    {
                        this->mLiveSettings.SetValue(SLabSettings::NpcSizeMultiplier, value);
                        this->OnLiveSettingsChanged();
                    },
                    std::make_unique<LinearSliderCore>(
                        mLabController->GetMinNpcSizeMultiplier(),
                        mLabController->GetMaxNpcSizeMultiplier()));

                npcsSizer->Add(
                    mNpcSizeMultiplierSlider,
                    wxGBPosition(0, 3),
                    wxGBSpan(1, 1),
                    wxEXPAND | wxALL,
                    CellBorder);
            }

            npcsBoxSizer->Add(npcsSizer, 0, wxALL, StaticBoxInsetMargin);
        }

        npcsBox->SetSizerAndFit(npcsBoxSizer);

        gridSizer->Add(
            npcsBox,
            wxGBPosition(0, 0),
            wxGBSpan(1, 1),
            wxEXPAND | wxALL | wxALIGN_CENTER_HORIZONTAL,
            CellBorder);
    }

    // Finalize panel

    for (int c = 0; c < gridSizer->GetCols(); ++c)
        gridSizer->AddGrowableCol(c);

    panel->SetSizer(gridSizer);
}

void SettingsDialog::SyncControlsWithSettings(Settings<SLabSettings> const & settings)
{
    // Simulator
    mElasticityAdjustmentSlider->SetValue(settings.GetValue<float>(SLabSettings::ElasticityAdjustment));
    mStaticFrictionAdjustmentSlider->SetValue(settings.GetValue<float>(SLabSettings::StaticFrictionAdjustment));
    mKineticFrictionAdjustmentSlider->SetValue(settings.GetValue<float>(SLabSettings::KineticFrictionAdjustment));
    mMassAdjustmentSlider->SetValue(settings.GetValue<float>(SLabSettings::MassAdjustment));
    mGravityAdjustmentSlider->SetValue(settings.GetValue<float>(SLabSettings::GravityAdjustment));
    mGlobalDampingAdjustmentSlider->SetValue(settings.GetValue<float>(SLabSettings::GlobalDampingAdjustment));
    mSeaLevelSlider->SetValue(settings.GetValue<float>(SLabSettings::SeaLevel));
    mSpringReductionFractionAdjustmentSlider->SetValue(settings.GetValue<float>(SLabSettings::SpringReductionFractionAdjustment));
    mSpringDampingCoefficientAdjustmentSlider->SetValue(settings.GetValue<float>(SLabSettings::SpringDampingCoefficientAdjustment));
    mWaterFrictionDragCoefficientAdjustmentSlider->SetValue(settings.GetValue<float>(SLabSettings::WaterFrictionDragAdjustment));
    mBuoyancyAdjustmentSlider->SetValue(settings.GetValue<float>(SLabSettings::BuoyancyAdjustment));

    // NPCs
    mHumanNpcEquilibriumTorqueStiffnessCoefficientSlider->SetValue(settings.GetValue<float>(SLabSettings::HumanNpcEquilibriumTorqueStiffnessCoefficient));
    mHumanNpcEquilibriumTorqueDampingCoefficientSlider->SetValue(settings.GetValue<float>(SLabSettings::HumanNpcEquilibriumTorqueDampingCoefficient));
    mHumanNpcWalkingSpeedAdjustmentSlider->SetValue(settings.GetValue<float>(SLabSettings::HumanNpcWalkingSpeedAdjustment));
    mNpcSizeMultiplierSlider->SetValue(settings.GetValue<float>(SLabSettings::NpcSizeMultiplier));
}

void SettingsDialog::OnLiveSettingsChanged()
{
    assert(!!mSettingsManager);

    // Enforce settings that have just changed
    mSettingsManager->EnforceDirtySettings(mLiveSettings);

    // We're back in sync
    mLiveSettings.ClearAllDirty();

	// Remember that we have changed since we were opened
	mHasBeenDirtyInCurrentSession = true;
	mAreSettingsDirtyWrtDefaults = true; // Best effort, assume each change deviates from defaults
	ReconcileDirtyState();
}

void SettingsDialog::ReconcileDirtyState()
{
    //
    // Update buttons' state based on dirty state
    //

	mRevertToDefaultsButton->Enable(mAreSettingsDirtyWrtDefaults);
    mUndoButton->Enable(mHasBeenDirtyInCurrentSession);
}