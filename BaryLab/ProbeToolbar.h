/***************************************************************************************
* Original Author:      Gabriele Giuseppini
* Created:              2020-05-15
* Copyright:            Gabriele Giuseppini  (https://github.com/GabrieleGiuseppini)
***************************************************************************************/
#pragma once

#include "UIControls/ScalarTimeSeriesProbeControl.h"

#include <Game/IGameEventHandlers.h>

#include <wx/sizer.h>
#include <wx/stattext.h>
#include <wx/textctrl.h>
#include <wx/wx.h>

#include <memory>
#include <string>
#include <unordered_map>

class ProbeToolbar final
    : public wxPanel
    , public IBLabEventHandler
{
public:

    ProbeToolbar(wxWindow* parent);

    virtual ~ProbeToolbar();

    void Update();

public:

    //
    // Simulation event handlers
    //

    void OnBLabReset() override;

    void OnSubjectParticleBarycentricCoordinatesWrtOriginTriangleChanged(std::optional<bcoords3f> const & coordinates) override;

    void OnSubjectParticleConstrainedRegimeUpdated(std::optional<ConstrainedRegimeParticleProbe> const & constrainedRegimeParticleProbe) override;

    void OnSubjectParticlePhysicsUpdated(std::optional<PhysicsParticleProbe> const & probe) override;

    void OnHumanNpcBehaviorChanged(std::optional<std::string> behavior) override;

    void OnHumanNpcStateQuantityChanged(std::optional<std::tuple<std::string, std::string>> nameAndValue) override;

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

    wxTextCtrl * mParticleVelocityTextCtrl;

    wxTextCtrl * mHumanBehaviorTextCtrl;
    wxStaticText * mHumanStateQuantityNameLabel;
    wxTextCtrl * mHumanStateQuantityTextCtrl;

    std::unordered_map<std::string, std::unique_ptr<ScalarTimeSeriesProbeControl>> mCustomProbes;

    wxBoxSizer * mProbesSizer;
};
