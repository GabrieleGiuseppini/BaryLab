/***************************************************************************************
* Original Author:      Gabriele Giuseppini
* Created:              2020-05-15
* Copyright:            Gabriele Giuseppini  (https://github.com/GabrieleGiuseppini)
***************************************************************************************/
#include "ScalarTimeSeriesProbeControl.h"

#include <wx/dcbuffer.h>

#include <algorithm>
#include <cassert>
#include <iomanip>
#include <limits>
#include <numeric>
#include <sstream>

ScalarTimeSeriesProbeControl::ScalarTimeSeriesProbeControl(
    wxWindow * parent,
    int width,
    int height)
    : wxPanel(
        parent,
        wxID_ANY,
        wxDefaultPosition,
        wxDefaultSize,
        wxBORDER_SIMPLE)
    , mWidth(width)
    , mHeight(height)
    , mBufferedDCBitmap()
    , mTimeSeriesPen(wxColor("BLACK"), 2, wxPENSTYLE_SOLID)
    , mGridPen(wxColor(0xa0, 0xa0, 0xa0), 1, wxPENSTYLE_SOLID)
{
    SetMinSize(wxSize(mWidth, mHeight));
    SetMaxSize(wxSize(mWidth, mHeight));

#ifdef __WXMSW__
    SetDoubleBuffered(true);
#endif

    SetBackgroundColour(wxColour("WHITE"));

    wxFont font(wxFontInfo(wxSize(8, 8)).Family(wxFONTFAMILY_TELETYPE));
    SetFont(font);

    Connect(this->GetId(), wxEVT_LEFT_DOWN, (wxObjectEventFunction)&ScalarTimeSeriesProbeControl::OnMouseClick);
    Connect(this->GetId(), wxEVT_PAINT, (wxObjectEventFunction)&ScalarTimeSeriesProbeControl::OnPaint);
    Connect(this->GetId(), wxEVT_ERASE_BACKGROUND, (wxObjectEventFunction)&ScalarTimeSeriesProbeControl::OnEraseBackground);

    Reset();
}

ScalarTimeSeriesProbeControl::~ScalarTimeSeriesProbeControl()
{
}

void ScalarTimeSeriesProbeControl::RegisterSample(float value)
{
    mMaxValue = std::max(mMaxValue, value);
    mMinValue = std::min(mMinValue, value);

    mSamples.emplace(
        [](float) {},
        value);

    mMaxAbsoluteValue = std::max(mMaxAbsoluteValue, std::abs(value));
}

void ScalarTimeSeriesProbeControl::Update()
{
    Refresh();
}

void ScalarTimeSeriesProbeControl::Reset()
{
    mMaxValue = std::numeric_limits<float>::lowest();
    mMinValue = std::numeric_limits<float>::max();

    mGridValueSize = 0.0f;

    mSamples.clear();
    mMaxAbsoluteValue = std::numeric_limits<float>::min();
}

///////////////////////////////////////////////////////////////////////////////////////

void ScalarTimeSeriesProbeControl::OnMouseClick(wxMouseEvent & /*event*/)
{
    // Reset extent and limit
    mMaxValue = std::numeric_limits<float>::lowest();
    mMinValue = std::numeric_limits<float>::max();
    mMaxAbsoluteValue = std::numeric_limits<float>::min();
    for (auto s : mSamples)
    {
        mMaxValue = std::max(mMaxValue, s);
        mMinValue = std::min(mMinValue, s);
        mMaxAbsoluteValue = std::max(mMaxAbsoluteValue, std::abs(s));
    }

    Refresh();
}

void ScalarTimeSeriesProbeControl::OnPaint(wxPaintEvent & /*event*/)
{
    if (this->GetSize() != wxSize(0, 0))
    {
        if (!mBufferedDCBitmap || mBufferedDCBitmap->GetSize() != this->GetSize())
        {
            mBufferedDCBitmap = std::make_unique<wxBitmap>(this->GetSize());
        }

        wxBufferedPaintDC bufDc(this, *mBufferedDCBitmap);

        Render(bufDc);
    }
}

void ScalarTimeSeriesProbeControl::OnEraseBackground(wxPaintEvent & /*event*/)
{
    // Do nothing
}

int ScalarTimeSeriesProbeControl::MapValueToY(float value) const
{
    if (mMaxValue == mMinValue)
        return mHeight / 2;

    float y = static_cast<float>(mHeight - 4) * (value - mMinValue) / (mMaxValue - mMinValue);
    return mHeight - 3 - static_cast<int>(round(y));
}

void ScalarTimeSeriesProbeControl::Render(wxDC & dc)
{
    dc.Clear();

    if (!mSamples.empty())
    {
        //
        // Check if need to resize grid
        //

        // Calculate new grid step
        float numberOfGridLines = 6.0f;
        float const currentValueExtent = mMaxValue - mMinValue;
        if (currentValueExtent > 0.0f)
        {
            if (mGridValueSize == 0.0f)
                mGridValueSize = currentValueExtent / 6.0f;

            // Number of grid lines we would have with the current extent
            numberOfGridLines = currentValueExtent / mGridValueSize;
            if (numberOfGridLines > 20.0f)
            {
                // Recalc
                mGridValueSize = currentValueExtent / 6.0f;
                numberOfGridLines = 6.0f;
            }
        }

        static const int xGridStepSize = mWidth / 6;
        int yGridStepSize = std::min(mWidth, mHeight) / static_cast<int>(ceil(numberOfGridLines));


        //
        // Draw grid
        //

        dc.SetPen(mGridPen);

        for (int y = yGridStepSize; y < mHeight - 1; y += yGridStepSize)
        {
            dc.DrawLine(0, y, mWidth - 1, y);
        }

        for (int x = xGridStepSize; x < mWidth - 1; x += xGridStepSize)
        {
            dc.DrawLine(x, 0, x, mHeight - 1);
        }


        //
        // Draw chart
        //

        dc.SetPen(mTimeSeriesPen);

        auto it = mSamples.cbegin();
        int lastX = mWidth - 2;
        int lastY = MapValueToY(*it);
        ++it;

        if (it == mSamples.cend())
        {
            // Draw just a point
            dc.DrawPoint(lastX, lastY);
        }
        else
        {
            // Draw lines
            do
            {
                int newX = lastX - 1;
                if (newX == 0)
                    break;

                int newY = MapValueToY(*it);

                dc.DrawLine(newX, newY, lastX, lastY);

                lastX = newX;
                lastY = newY;

                ++it;
            }
            while (it != mSamples.cend());
        }


        //
        // Draw label
        //

        std::stringstream ss;
        ss << std::fixed << std::setprecision(3) << *mSamples.cbegin() << " (" << mMaxAbsoluteValue << ")";

        wxString labelText(ss.str());
        dc.DrawText(labelText, 0, 1);
    }
}