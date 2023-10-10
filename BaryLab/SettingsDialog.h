/***************************************************************************************
 * Original Author:     Gabriele Giuseppini
 * Created:             2020-05-23
 * Copyright:           Gabriele Giuseppini  (https://github.com/GabrieleGiuseppini)
 ***************************************************************************************/
#pragma once

#include "SettingsManager.h"

#include "UIControls/SliderControl.h"

#include <BLabCoreLib/LabController.h>
#include <BLabCoreLib/Settings.h>

#include <wx/button.h>
#include <wx/checkbox.h>
#include <wx/radiobox.h>

#include <memory>
#include <vector>

class SettingsDialog : public wxFrame
{
public:

    SettingsDialog(
        wxWindow * parent,
        std::shared_ptr<SettingsManager> settingsManager,
		std::shared_ptr<LabController> labController);

    virtual ~SettingsDialog();

    void Open();

private:

	void OnRevertToDefaultsButton(wxCommandEvent& event);
    void OnOkButton(wxCommandEvent & event);
    void OnCancelButton(wxCommandEvent & event);
    void OnUndoButton(wxCommandEvent & event);

    void OnCloseButton(wxCloseEvent & event);

private:

    //////////////////////////////////////////////////////
    // Control tabs
    //////////////////////////////////////////////////////

    // Simulator
    SliderControl<float> * mElasticitySlider;
    SliderControl<float> * mStaticFrictionSlider;
    SliderControl<float> * mKineticFrictionSlider;
    SliderControl<float> * mMassAdjustmentSlider;
    SliderControl<float> * mGravityAdjustmentSlider;
    SliderControl<float> * mSpringReductionFractionSlider;
    SliderControl<float> * mSpringDampingCoefficientSlider;

    //////////////////////////////////////////////////////

    // Buttons
	wxButton * mRevertToDefaultsButton;
    wxButton * mOkButton;
    wxButton * mCancelButton;
    wxButton * mUndoButton;

private:

    void DoCancel();
    void DoClose();

    void PopulateSimulatorPanel(wxPanel * panel);

    void SyncControlsWithSettings(Settings<SLabSettings> const & settings);

    void OnLiveSettingsChanged();
    void ReconcileDirtyState();

private:

    wxWindow * const mParent;
    std::shared_ptr<SettingsManager> mSettingsManager;
	std::shared_ptr<LabController> mLabController;

    //
    // State
    //

    // The current settings, always enforced
    Settings<SLabSettings> mLiveSettings;

    // The settings when the dialog was last opened
    Settings<SLabSettings> mCheckpointSettings;

    // Tracks whether the user has changed any settings since the dialog
    // was last opened. When false there's a guarantee that the current live
    // settings have not been modified.
    bool mHasBeenDirtyInCurrentSession;

	// Tracks whether the current settings are (possibly) dirty wrt the defaults.
	// Best effort, we assume all changes deviate from the default.
	bool mAreSettingsDirtyWrtDefaults;
};
