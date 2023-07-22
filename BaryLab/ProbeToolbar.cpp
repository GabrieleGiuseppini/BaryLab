/***************************************************************************************
* Original Author:      Gabriele Giuseppini
* Created:              2020-05-15
* Copyright:            Gabriele Giuseppini  (https://github.com/GabrieleGiuseppini)
***************************************************************************************/
#include "ProbeToolbar.h"

#include <wx/gbsizer.h>
#include <wx/stattext.h>

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

    {
        int constexpr HMargin = 10;

        wxBoxSizer * hSizer = new wxBoxSizer(wxHORIZONTAL);

        hSizer->AddSpacer(HMargin);

        // Fixed scalars
        {
            int constexpr TextCtrlWidth = 100;

            wxGridBagSizer * gridSizer = new wxGridBagSizer(0, 0);

            // L1
            {
                auto label = new wxStaticText(this, wxID_ANY, _("l1:"));
                gridSizer->Add(
                    label,
                    wxGBPosition(0, 0),
                    wxGBSpan(1, 1),
                    wxALIGN_LEFT | wxALIGN_CENTER_VERTICAL,
                    0);

                mBarycentricCoordinateL1TextCtrl = new wxTextCtrl(this, wxID_ANY, "", wxDefaultPosition, wxSize(TextCtrlWidth, -1), wxTE_RIGHT | wxTE_READONLY);
                gridSizer->Add(
                    mBarycentricCoordinateL1TextCtrl,
                    wxGBPosition(0, 1),
                    wxGBSpan(1, 1),
                    wxEXPAND,
                    0);
            }

            // L2
            {
                auto label = new wxStaticText(this, wxID_ANY, _("l2:"));
                gridSizer->Add(
                    label,
                    wxGBPosition(1, 0),
                    wxGBSpan(1, 1),
                    wxALIGN_LEFT | wxALIGN_CENTER_VERTICAL,
                    0);

                mBarycentricCoordinateL2TextCtrl = new wxTextCtrl(this, wxID_ANY, "", wxDefaultPosition, wxSize(TextCtrlWidth, -1), wxTE_RIGHT | wxTE_READONLY);
                gridSizer->Add(
                    mBarycentricCoordinateL2TextCtrl,
                    wxGBPosition(1, 1),
                    wxGBSpan(1, 1),
                    wxEXPAND,
                    0);
            }

            // L3
            {
                auto label = new wxStaticText(this, wxID_ANY, _("l3:"));
                gridSizer->Add(
                    label,
                    wxGBPosition(2, 0),
                    wxGBSpan(1, 1),
                    wxALIGN_LEFT | wxALIGN_CENTER_VERTICAL,
                    0);

                mBarycentricCoordinateL3TextCtrl = new wxTextCtrl(this, wxID_ANY, "", wxDefaultPosition, wxSize(TextCtrlWidth, -1), wxTE_RIGHT | wxTE_READONLY);
                gridSizer->Add(
                    mBarycentricCoordinateL3TextCtrl,
                    wxGBPosition(2, 1),
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

        // Probes
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

void ProbeToolbar::OnReset()
{
    mBarycentricCoordinateL1TextCtrl->SetValue("");
    mBarycentricCoordinateL2TextCtrl->SetValue("");
    mBarycentricCoordinateL3TextCtrl->SetValue("");

    for (auto const & p : mCustomProbes)
    {
        p.second->Reset();
    }
}

void ProbeToolbar::OnSubjectParticleBarycentricCoordinatesChanged(std::optional<vec3f> const & coordinates)
{
    std::string l1;
    std::string l2;
    std::string l3;

    if (coordinates)
    {
        {
            std::ostringstream ss;
            ss.fill('0');
            ss << std::fixed << std::setprecision(2) << coordinates->x;

            l1 = ss.str();
        }

        {
            std::ostringstream ss;
            ss.fill('0');
            ss << std::fixed << std::setprecision(2) << coordinates->y;

            l2 = ss.str();
        }

        {
            std::ostringstream ss;
            ss.fill('0');
            ss << std::fixed << std::setprecision(2) << coordinates->z;

            l3 = ss.str();
        }
    }

    mBarycentricCoordinateL1TextCtrl->SetValue(l1);
    mBarycentricCoordinateL2TextCtrl->SetValue(l2);
    mBarycentricCoordinateL3TextCtrl->SetValue(l3);
}

void ProbeToolbar::OnCustomProbe(
    std::string const & name,
    float value)
{
    auto & probe = mCustomProbes[name];
    if (!probe)
    {
        probe = AddScalarTimeSeriesProbe(name, 100);
        mProbesSizer->Layout();
    }

    probe->RegisterSample(value);
}