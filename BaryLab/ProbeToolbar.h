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

    void OnSubjectParticleBarycentricCoordinatesWrtOriginTriangleChanged(std::optional<vec3f> const & coordinates) override;

    void OnSubjectParticleUpdated(std::optional<ParticleProbe> const & particleProbe) override;

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

    wxTextCtrl * mOriginTriangleBarycentricCoordinateL1TextCtrl;
    wxTextCtrl * mOriginTriangleBarycentricCoordinateL2TextCtrl;
    wxTextCtrl * mOriginTriangleBarycentricCoordinateL3TextCtrl;

    wxTextCtrl * mProbeTriangleIndexTextCtrl;
    wxTextCtrl * mProbeBarycentricCoordinateL1TextCtrl;
    wxTextCtrl * mProbeBarycentricCoordinateL2TextCtrl;
    wxTextCtrl * mProbeBarycentricCoordinateL3TextCtrl;

    std::unordered_map<std::string, std::unique_ptr<ScalarTimeSeriesProbeControl>> mCustomProbes;

    wxBoxSizer * mProbesSizer;
};
