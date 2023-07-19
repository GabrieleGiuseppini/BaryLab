/***************************************************************************************
* Original Author:      Gabriele Giuseppini
* Created:              2020-05-15
* Copyright:            Gabriele Giuseppini  (https://github.com/GabrieleGiuseppini)
***************************************************************************************/
#pragma once

#include "UIControls/ScalarTimeSeriesProbeControl.h"

#include <BLabCoreLib/IEventHandler.h>

#include <wx/sizer.h>
#include <wx/textctrl.h>
#include <wx/wx.h>

#include <memory>
#include <string>
#include <unordered_map>

class ProbeToolbar final
    : public wxPanel
    , public IEventHandler
{
public:

    ProbeToolbar(wxWindow* parent);

    virtual ~ProbeToolbar();

    void Update();

public:

    //
    // Simulation event handlers
    //

    void OnReset() override;

    void OnSubjectParticleBarycentricCoordinatesChanged(vec3f const & coordinates) override;

    void OnCustomProbe(
        std::string const & name,
        float value) override;

private:

    bool IsActive() const
    {
        return this->IsShown();
    }

    std::unique_ptr<ScalarTimeSeriesProbeControl> AddScalarTimeSeriesProbe(
        std::string const & name,
        int sampleCount);

private:

    //
    // UI
    //

    wxTextCtrl * mBarycentricCoordinateL1TextCtrl;
    wxTextCtrl * mBarycentricCoordinateL2TextCtrl;
    wxTextCtrl * mBarycentricCoordinateL3TextCtrl;    

    std::unordered_map<std::string, std::unique_ptr<ScalarTimeSeriesProbeControl>> mCustomProbes;

    wxBoxSizer * mProbesSizer;
};
