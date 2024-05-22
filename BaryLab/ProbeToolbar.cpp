/***************************************************************************************
* Original Author:      Gabriele Giuseppini
* Created:              2020-05-15
* Copyright:            Gabriele Giuseppini  (https://github.com/GabrieleGiuseppini)
***************************************************************************************/
#include "ProbeToolbar.h"

#include <wx/gbsizer.h>

#include <cassert>
#include <iomanip>
#include <sstream>

static constexpr int TopPadding = 2;
static constexpr int ProbePadding = 10;
static constexpr int ProbeHeight = 80;

ProbeToolbar::ProbeToolbar(wxWindow* parent)
    : wxPanel(
        parent,
        wxID_ANY,
        wxDefaultPosition,
        wxSize(-1, -1),
        wxBORDER_SIMPLE | wxCLIP_CHILDREN)
{
#ifdef __WXMSW__
    SetDoubleBuffered(true);
#endif

    SetBackgroundColour(wxSystemSettings::GetColour(wxSYS_COLOUR_BTNFACE));

    wxBoxSizer * vSizer = new wxBoxSizer(wxVERTICAL);

    //
    // Setup UI
    //

    int constexpr TextCtrlWidth = 100;
    int constexpr TextCtrlWidth2 = 200;

    {
        int constexpr HMargin = 10;

        wxBoxSizer * hSizer = new wxBoxSizer(wxHORIZONTAL);

        hSizer->AddSpacer(HMargin);

        // Origin triangle
        {
            wxGridBagSizer * gridSizer = new wxGridBagSizer(0, 0);

            // L1
            {
                auto label = new wxStaticText(this, wxID_ANY, _("l1:"));
                gridSizer->Add(
                    label,
                    wxGBPosition(1, 0),
                    wxGBSpan(1, 1),
                    wxALIGN_LEFT | wxALIGN_CENTER_VERTICAL,
                    0);

                mOriginTriangleBarycentricCoordinateL1TextCtrl = new wxTextCtrl(this, wxID_ANY, "", wxDefaultPosition, wxSize(TextCtrlWidth, -1), wxTE_RIGHT | wxTE_READONLY);
                gridSizer->Add(
                    mOriginTriangleBarycentricCoordinateL1TextCtrl,
                    wxGBPosition(1, 1),
                    wxGBSpan(1, 1),
                    wxEXPAND,
                    0);
            }

            // L2
            {
                auto label = new wxStaticText(this, wxID_ANY, _("l2:"));
                gridSizer->Add(
                    label,
                    wxGBPosition(2, 0),
                    wxGBSpan(1, 1),
                    wxALIGN_LEFT | wxALIGN_CENTER_VERTICAL,
                    0);

                mOriginTriangleBarycentricCoordinateL2TextCtrl = new wxTextCtrl(this, wxID_ANY, "", wxDefaultPosition, wxSize(TextCtrlWidth, -1), wxTE_RIGHT | wxTE_READONLY);
                gridSizer->Add(
                    mOriginTriangleBarycentricCoordinateL2TextCtrl,
                    wxGBPosition(2, 1),
                    wxGBSpan(1, 1),
                    wxEXPAND,
                    0);
            }

            // L3
            {
                auto label = new wxStaticText(this, wxID_ANY, _("l3:"));
                gridSizer->Add(
                    label,
                    wxGBPosition(3, 0),
                    wxGBSpan(1, 1),
                    wxALIGN_LEFT | wxALIGN_CENTER_VERTICAL,
                    0);

                mOriginTriangleBarycentricCoordinateL3TextCtrl = new wxTextCtrl(this, wxID_ANY, "", wxDefaultPosition, wxSize(TextCtrlWidth, -1), wxTE_RIGHT | wxTE_READONLY);
                gridSizer->Add(
                    mOriginTriangleBarycentricCoordinateL3TextCtrl,
                    wxGBPosition(3, 1),
                    wxGBSpan(1, 1),
                    wxEXPAND,
                    0);
            }

            {
                // Empty
                gridSizer->Add(
                    TextCtrlWidth,
                    mOriginTriangleBarycentricCoordinateL1TextCtrl->GetSize().GetHeight(),
                    wxGBPosition(0, 0),
                    wxGBSpan(1, 2),
                    wxALIGN_LEFT | wxALIGN_CENTER_VERTICAL,
                    0);
            }

            hSizer->Add(
                gridSizer,
                0,
                0,
                0);
        }

        hSizer->AddSpacer(HMargin);

        // Constrained State Probe
        {
            wxGridBagSizer * gridSizer = new wxGridBagSizer(0, 0);

            // Triangle index
            {
                auto label = new wxStaticText(this, wxID_ANY, _("t:"));
                gridSizer->Add(
                    label,
                    wxGBPosition(0, 0),
                    wxGBSpan(1, 1),
                    wxALIGN_LEFT | wxALIGN_CENTER_VERTICAL,
                    0);

                mProbeTriangleIndexTextCtrl = new wxTextCtrl(this, wxID_ANY, "", wxDefaultPosition, wxSize(TextCtrlWidth, -1), wxTE_RIGHT | wxTE_READONLY);
                gridSizer->Add(
                    mProbeTriangleIndexTextCtrl,
                    wxGBPosition(0, 1),
                    wxGBSpan(1, 1),
                    wxEXPAND,
                    0);
            }

            // L1
            {
                auto label = new wxStaticText(this, wxID_ANY, _("l1:"));
                gridSizer->Add(
                    label,
                    wxGBPosition(1, 0),
                    wxGBSpan(1, 1),
                    wxALIGN_LEFT | wxALIGN_CENTER_VERTICAL,
                    0);

                mProbeBarycentricCoordinateL1TextCtrl = new wxTextCtrl(this, wxID_ANY, "", wxDefaultPosition, wxSize(TextCtrlWidth, -1), wxTE_RIGHT | wxTE_READONLY);
                gridSizer->Add(
                    mProbeBarycentricCoordinateL1TextCtrl,
                    wxGBPosition(1, 1),
                    wxGBSpan(1, 1),
                    wxEXPAND,
                    0);
            }

            // L2
            {
                auto label = new wxStaticText(this, wxID_ANY, _("l2:"));
                gridSizer->Add(
                    label,
                    wxGBPosition(2, 0),
                    wxGBSpan(1, 1),
                    wxALIGN_LEFT | wxALIGN_CENTER_VERTICAL,
                    0);

                mProbeBarycentricCoordinateL2TextCtrl = new wxTextCtrl(this, wxID_ANY, "", wxDefaultPosition, wxSize(TextCtrlWidth, -1), wxTE_RIGHT | wxTE_READONLY);
                gridSizer->Add(
                    mProbeBarycentricCoordinateL2TextCtrl,
                    wxGBPosition(2, 1),
                    wxGBSpan(1, 1),
                    wxEXPAND,
                    0);
            }

            // L3
            {
                auto label = new wxStaticText(this, wxID_ANY, _("l3:"));
                gridSizer->Add(
                    label,
                    wxGBPosition(3, 0),
                    wxGBSpan(1, 1),
                    wxALIGN_LEFT | wxALIGN_CENTER_VERTICAL,
                    0);

                mProbeBarycentricCoordinateL3TextCtrl = new wxTextCtrl(this, wxID_ANY, "", wxDefaultPosition, wxSize(TextCtrlWidth, -1), wxTE_RIGHT | wxTE_READONLY);
                gridSizer->Add(
                    mProbeBarycentricCoordinateL3TextCtrl,
                    wxGBPosition(3, 1),
                    wxGBSpan(1, 1),
                    wxEXPAND,
                    0);
            }

            hSizer->Add(
                gridSizer,
                0,
                0,
                0);
        }

        hSizer->AddSpacer(HMargin);

        // Physics Probe
        {
            wxGridBagSizer * gridSizer = new wxGridBagSizer(0, 0);

            // Velocity
            {
                auto label = new wxStaticText(this, wxID_ANY, _("v:"));
                gridSizer->Add(
                    label,
                    wxGBPosition(0, 0),
                    wxGBSpan(1, 1),
                    wxALIGN_LEFT | wxALIGN_CENTER_VERTICAL,
                    0);

                mParticleVelocityTextCtrl = new wxTextCtrl(this, wxID_ANY, "", wxDefaultPosition, wxSize(TextCtrlWidth, -1), wxTE_RIGHT | wxTE_READONLY);
                gridSizer->Add(
                    mParticleVelocityTextCtrl,
                    wxGBPosition(0, 1),
                    wxGBSpan(1, 1),
                    wxEXPAND,
                    0);
            }

            hSizer->Add(
                gridSizer,
                0,
                0,
                0);
        }

        hSizer->AddSpacer(80);

        // Human state
        {
            wxGridBagSizer * gridSizer = new wxGridBagSizer(0, 0);

            // Behavior
            {
                auto label = new wxStaticText(this, wxID_ANY, _("Behavior:"));
                gridSizer->Add(
                    label,
                    wxGBPosition(0, 0),
                    wxGBSpan(1, 1),
                    wxALIGN_LEFT | wxALIGN_CENTER_VERTICAL,
                    0);

                mHumanBehaviorTextCtrl = new wxTextCtrl(this, wxID_ANY, "", wxDefaultPosition, wxSize(TextCtrlWidth2, -1), wxTE_RIGHT | wxTE_READONLY);
                gridSizer->Add(
                    mHumanBehaviorTextCtrl,
                    wxGBPosition(0, 1),
                    wxGBSpan(1, 1),
                    wxEXPAND,
                    0);
            }

            // State Quantity
            {
                mHumanStateQuantityNameLabel = new wxStaticText(this, wxID_ANY, "", wxDefaultPosition, wxSize(TextCtrlWidth, -1));
                gridSizer->Add(
                    mHumanStateQuantityNameLabel,
                    wxGBPosition(1, 0),
                    wxGBSpan(1, 1),
                    wxALIGN_LEFT | wxALIGN_CENTER_VERTICAL,
                    0);

                mHumanStateQuantityTextCtrl = new wxTextCtrl(this, wxID_ANY, "", wxDefaultPosition, wxSize(TextCtrlWidth2, -1), wxTE_RIGHT | wxTE_READONLY);
                gridSizer->Add(
                    mHumanStateQuantityTextCtrl,
                    wxGBPosition(1, 1),
                    wxGBSpan(1, 1),
                    wxEXPAND,
                    0);
            }

            hSizer->Add(
                gridSizer,
                0,
                0,
                0);
        }

        hSizer->AddSpacer(HMargin);

        // Hard numbers
        {
            wxGridBagSizer * gridSizer = new wxGridBagSizer(0, 0);

            // Update time
            {
                auto label = new wxStaticText(this, wxID_ANY, _("Update Time (ms):"));
                gridSizer->Add(
                    label,
                    wxGBPosition(0, 0),
                    wxGBSpan(1, 1),
                    wxALIGN_LEFT | wxALIGN_CENTER_VERTICAL,
                    0);

                mUpdateTimeTextCtrl = new wxTextCtrl(this, wxID_ANY, "0", wxDefaultPosition, wxSize(TextCtrlWidth2, -1), wxTE_RIGHT | wxTE_READONLY);
                gridSizer->Add(
                    mUpdateTimeTextCtrl,
                    wxGBPosition(0, 1),
                    wxGBSpan(1, 1),
                    wxEXPAND,
                    0);
            }

            // RenderUpload time
            {
                auto label = new wxStaticText(this, wxID_ANY, _("Render Upl. Time (ms):"));
                gridSizer->Add(
                    label,
                    wxGBPosition(1, 0),
                    wxGBSpan(1, 1),
                    wxALIGN_LEFT | wxALIGN_CENTER_VERTICAL,
                    0);

                mRenderUploadTimeTextCtrl = new wxTextCtrl(this, wxID_ANY, "0", wxDefaultPosition, wxSize(TextCtrlWidth2, -1), wxTE_RIGHT | wxTE_READONLY);
                gridSizer->Add(
                    mRenderUploadTimeTextCtrl,
                    wxGBPosition(1, 1),
                    wxGBSpan(1, 1),
                    wxEXPAND,
                    0);
            }

            // Inside
            {
                auto label = new wxStaticText(this, wxID_ANY, _("Humans Inside:"));
                gridSizer->Add(
                    label,
                    wxGBPosition(2, 0),
                    wxGBSpan(1, 1),
                    wxALIGN_LEFT | wxALIGN_CENTER_VERTICAL,
                    0);

                mHumanNpcInsideShipCountTextCtrl = new wxTextCtrl(this, wxID_ANY, "0", wxDefaultPosition, wxSize(TextCtrlWidth2, -1), wxTE_RIGHT | wxTE_READONLY);
                gridSizer->Add(
                    mHumanNpcInsideShipCountTextCtrl,
                    wxGBPosition(2, 1),
                    wxGBSpan(1, 1),
                    wxEXPAND,
                    0);
            }

            // Outside
            {
                auto label = new wxStaticText(this, wxID_ANY, _("Humans Outside:"));
                gridSizer->Add(
                    label,
                    wxGBPosition(3, 0),
                    wxGBSpan(1, 1),
                    wxALIGN_LEFT | wxALIGN_CENTER_VERTICAL,
                    0);

                mHumanNpcOutsideShipCountTextCtrl = new wxTextCtrl(this, wxID_ANY, "0", wxDefaultPosition, wxSize(TextCtrlWidth2, -1), wxTE_RIGHT | wxTE_READONLY);
                gridSizer->Add(
                    mHumanNpcOutsideShipCountTextCtrl,
                    wxGBPosition(3, 1),
                    wxGBSpan(1, 1),
                    wxEXPAND,
                    0);
            }

            hSizer->Add(
                gridSizer,
                0,
                0,
                0);
        }

        hSizer->AddSpacer(HMargin);

        // Custom probes
        {
            mProbesSizer = new wxBoxSizer(wxHORIZONTAL);

            hSizer->Add(
                mProbesSizer,
                0,
                0,
                0);
        }

        hSizer->AddSpacer(HMargin);

        vSizer->Add(
            hSizer,
            1,
            wxALL,
            5);
    }

    SetSizer(vSizer);
}

ProbeToolbar::~ProbeToolbar()
{
}

void ProbeToolbar::Update()
{
    //
    // Update all probes
    //

    if (IsActive())
    {
        for (auto const & p : mCustomProbes)
        {
            p.second->Update();
        }
    }
}

std::unique_ptr<ScalarTimeSeriesProbeControl> ProbeToolbar::AddScalarTimeSeriesProbe(
    std::string const & name,
    int sampleCount)
{
    wxBoxSizer * sizer = new wxBoxSizer(wxVERTICAL);

    sizer->AddSpacer(TopPadding);

    auto probe = std::make_unique<ScalarTimeSeriesProbeControl>(this, sampleCount, ProbeHeight);
    sizer->Add(probe.get(), 1, wxALIGN_CENTRE, 0);

    wxStaticText * label = new wxStaticText(this, wxID_ANY, name, wxDefaultPosition, wxDefaultSize, wxALIGN_CENTRE_HORIZONTAL);
    sizer->Add(label, 0, wxALIGN_CENTRE, 0);

    mProbesSizer->Add(sizer, 1, wxLEFT | wxRIGHT, ProbePadding);

    return probe;
}

///////////////////////////////////////////////////////////////////////////////////////

void ProbeToolbar::OnBLabReset()
{
    mOriginTriangleBarycentricCoordinateL1TextCtrl->SetValue("");
    mOriginTriangleBarycentricCoordinateL2TextCtrl->SetValue("");
    mOriginTriangleBarycentricCoordinateL3TextCtrl->SetValue("");

    mProbeTriangleIndexTextCtrl->SetValue("");
    mProbeBarycentricCoordinateL1TextCtrl->SetValue("");
    mProbeBarycentricCoordinateL2TextCtrl->SetValue("");
    mProbeBarycentricCoordinateL3TextCtrl->SetValue("");

    mHumanNpcInsideShipCountTextCtrl->SetValue("0");
    mHumanNpcOutsideShipCountTextCtrl->SetValue("0");

    for (auto const & p : mCustomProbes)
    {
        p.second->Reset();
    }
}

void ProbeToolbar::OnSubjectParticleBarycentricCoordinatesWrtOriginTriangleChanged(std::optional<bcoords3f> const & coordinates)
{
    std::string l1;
    std::string l2;
    std::string l3;

    if (coordinates)
    {
        {
            std::ostringstream ss;
            ss.fill('0');
            ss << std::fixed << std::setprecision(2) << (*coordinates)[0];

            l1 = ss.str();
        }

        {
            std::ostringstream ss;
            ss.fill('0');
            ss << std::fixed << std::setprecision(2) << (*coordinates)[1];

            l2 = ss.str();
        }

        {
            std::ostringstream ss;
            ss.fill('0');
            ss << std::fixed << std::setprecision(2) << (*coordinates)[2];

            l3 = ss.str();
        }
    }

    mOriginTriangleBarycentricCoordinateL1TextCtrl->SetValue(l1);
    mOriginTriangleBarycentricCoordinateL2TextCtrl->SetValue(l2);
    mOriginTriangleBarycentricCoordinateL3TextCtrl->SetValue(l3);
}

void ProbeToolbar::OnSubjectParticleConstrainedRegimeUpdated(std::optional<AbsoluteTriangleBCoords> const & constrainedRegimeParticleProbe)
{
    std::string tIndex;
    std::string l1;
    std::string l2;
    std::string l3;

    if (constrainedRegimeParticleProbe)
    {
        {
            std::ostringstream ss;
            ss << constrainedRegimeParticleProbe->TriangleElementIndex;

            tIndex = ss.str();
        }

        {
            std::ostringstream ss;
            ss.fill('0');
            ss << std::fixed << std::setprecision(2) << constrainedRegimeParticleProbe->BCoords[0];

            l1 = ss.str();
        }

        {
            std::ostringstream ss;
            ss.fill('0');
            ss << std::fixed << std::setprecision(2) << constrainedRegimeParticleProbe->BCoords[1];

            l2 = ss.str();
        }

        {
            std::ostringstream ss;
            ss.fill('0');
            ss << std::fixed << std::setprecision(2) << constrainedRegimeParticleProbe->BCoords[2];

            l3 = ss.str();
        }
    }

    mProbeTriangleIndexTextCtrl->SetValue(tIndex);
    mProbeBarycentricCoordinateL1TextCtrl->SetValue(l1);
    mProbeBarycentricCoordinateL2TextCtrl->SetValue(l2);
    mProbeBarycentricCoordinateL3TextCtrl->SetValue(l3);
}

void ProbeToolbar::OnSubjectParticlePhysicsUpdated(std::optional<PhysicsParticleProbe> const & physicsParticleProbe)
{
    std::ostringstream ss;

    if (physicsParticleProbe)
    {
        ss.fill('0');
        ss << std::fixed << std::setprecision(2) << physicsParticleProbe->Velocity.x << ", " << physicsParticleProbe->Velocity.y;
    }

    mParticleVelocityTextCtrl->SetValue(ss.str());
}

void ProbeToolbar::OnHumanNpcBehaviorChanged(std::optional<std::string> behavior)
{
    mHumanBehaviorTextCtrl->SetValue(behavior.value_or(""));
}

void ProbeToolbar::OnHumanNpcStateQuantityChanged(std::optional<std::tuple<std::string, std::string>> nameAndValue)
{
    mHumanStateQuantityNameLabel->SetLabel(nameAndValue.has_value() ? std::get<0>(*nameAndValue) : "");
    mHumanStateQuantityTextCtrl->SetValue(nameAndValue.has_value() ? std::get<1>(*nameAndValue) : "");
}

void ProbeToolbar::OnUpdateTimeMeasured(
    float updateDurationMilliseconds,
    float renderUploadDurationMilliseconds)
{
    mUpdateTimeTextCtrl->SetValue(std::to_string(updateDurationMilliseconds));
    mRenderUploadTimeTextCtrl->SetValue(std::to_string(renderUploadDurationMilliseconds));
}

void ProbeToolbar::OnCustomProbe(
    std::string const & name,
    float value)
{
    auto & probe = mCustomProbes[name];
    if (!probe)
    {
        probe = AddScalarTimeSeriesProbe(name, 100);
        Layout();
    }

    probe->RegisterSample(value);
}

void ProbeToolbar::OnHumanNpcCountsUpdated(
    unsigned int insideShipCount,
    unsigned int outsideShipCount)
{
    mHumanNpcInsideShipCountTextCtrl->SetValue(std::to_string(insideShipCount));
    mHumanNpcOutsideShipCountTextCtrl->SetValue(std::to_string(outsideShipCount));

}