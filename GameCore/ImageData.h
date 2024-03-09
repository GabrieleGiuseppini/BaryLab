/***************************************************************************************
* Original Author:		Gabriele Giuseppini
* Created:				2020-05-15
* Copyright:			Gabriele Giuseppini  (https://github.com/GabrieleGiuseppini)
***************************************************************************************/
#pragma once

#include "Colors.h"
#include "Vectors.h"

#include <memory>

#pragma pack(push, 1)

template<typename TIntegralTag>
struct _IntegralSize
{
    using integral_type = int;

    integral_type width;
    integral_type height;

    constexpr _IntegralSize(
        integral_type _width,
        integral_type _height)
        : width(_width)
        , height(_height)
    {}

    static _IntegralSize<TIntegralTag> FromFloatRound(vec2f const & vec)
    {
        return _IntegralSize<TIntegralTag>(
            static_cast<integral_type>(std::round(vec.x)),
            static_cast<integral_type>(std::round(vec.y)));
    }

    static _IntegralSize<TIntegralTag> FromFloatFloor(vec2f const & vec)
    {
        return _IntegralSize<TIntegralTag>(
            static_cast<integral_type>(std::floor(vec.x)),
            static_cast<integral_type>(std::floor(vec.y)));
    }

    inline bool operator==(_IntegralSize<TIntegralTag> const & other) const
    {
        return this->width == other.width
            && this->height == other.height;
    }

    inline bool operator!=(_IntegralSize<TIntegralTag> const & other) const
    {
        return !(*this == other);
    }

    inline _IntegralSize<TIntegralTag> operator+(_IntegralSize<TIntegralTag> const & sz) const
    {
        return _IntegralSize<TIntegralTag>(
            this->width + sz.width,
            this->height + sz.height);
    }

    inline void operator+=(_IntegralSize<TIntegralTag> const & sz)
    {
        this->width += sz.width;
        this->height += sz.height;
    }

    inline _IntegralSize<TIntegralTag> operator*(integral_type factor) const
    {
        return _IntegralSize<TIntegralTag>(
            this->width * factor,
            this->height * factor);
    }

    inline size_t GetLinearSize() const
    {
        return this->width * this->height;
    }

    inline void Rotate90()
    {
        std::swap(width, height);
    }

    inline _IntegralSize<TIntegralTag> Union(_IntegralSize<TIntegralTag> const & other) const
    {
        return _IntegralSize<TIntegralTag>(
            std::max(this->width, other.width),
            std::max(this->height, other.height));
    }

    inline _IntegralSize<TIntegralTag> Intersection(_IntegralSize<TIntegralTag> const & other) const
    {
        return _IntegralSize<TIntegralTag>(
            std::min(this->width, other.width),
            std::min(this->height, other.height));
    }

    vec2f ToFloat() const
    {
        return vec2f(
            static_cast<float>(width),
            static_cast<float>(height));
    }

    template<typename TCoordsRatio>
    vec2f ToFractionalCoords(TCoordsRatio const & coordsRatio) const
    {
        assert(coordsRatio.inputUnits != 0.0f);

        return vec2f(
            static_cast<float>(width) / coordsRatio.inputUnits * coordsRatio.outputUnits,
            static_cast<float>(height) / coordsRatio.inputUnits * coordsRatio.outputUnits);
    }

    std::string ToString() const
    {
        std::stringstream ss;
        ss << "(" << width << " x " << height << ")";
        return ss.str();
    }
};

#pragma pack(pop)

template<typename TTag>
inline std::basic_ostream<char> & operator<<(std::basic_ostream<char> & os, _IntegralSize<TTag> const & is)
{
    os << is.ToString();
    return os;
}

using ImageSize = _IntegralSize<struct ImageTag>;

///////

template <typename TColor>
struct ImageData
{
public:

    using color_type = TColor;

public:

    ImageSize const Size;
    std::unique_ptr<color_type[]> Data;

    ImageData(
        int width,
        int height,
        std::unique_ptr<color_type[]> data)
        : Size(width, height)
        , Data(std::move(data))
    {
    }

    ImageData(
        ImageSize size,
        std::unique_ptr<color_type[]> data)
        : Size(size)
        , Data(std::move(data))
    {
    }

    ImageData(ImageData && other) noexcept
        : Size(other.Size)
        , Data(std::move(other.Data))
    {
    }
};

using RgbImageData = ImageData<rgbColor>;
using RgbaImageData = ImageData<rgbaColor>;
using Vec3fImageData = ImageData<vec3f>;