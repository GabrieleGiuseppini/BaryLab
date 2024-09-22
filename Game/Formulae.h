/***************************************************************************************
* Original Author:		Gabriele Giuseppini
* Created:				2021-07-30
* Copyright:			Gabriele Giuseppini  (https://github.com/GabrieleGiuseppini)
***************************************************************************************/
#pragma once

#include "GameParameters.h"

#include <GameCore/GameMath.h>

#include <cmath>
#include <optional>

namespace Physics
{

/*
 * Collection of some of the most reused formulae in the simulation.
 */
class Formulae final
{
public:

    static float CalculateAirDensity(
        float airTemperature,
        GameParameters const & gameParameters)
    {
        return GameParameters::AirMass
            / (1.0f + GameParameters::AirThermalExpansionCoefficient * (airTemperature - GameParameters::Temperature0))
            * gameParameters.AirDensityAdjustment;
    }

    static float CalculateWaterDensity(
        float waterTemperature,
        GameParameters const & gameParameters)
    {
        return GameParameters::WaterMass
            / (1.0f + GameParameters::WaterThermalExpansionCoefficient * (waterTemperature - GameParameters::Temperature0))
            * gameParameters.WaterDensityAdjustment;
    }

    // Converts a vectorial wind speed into the vectorial force it would have on a 1m2 surface
    static vec2f WindSpeedToForceDensity(
        vec2f windSpeed, // m/s
        float airDensity)
    {
        //  m/s -> Newton: F = 1/2 rho v**2 A
        return
            windSpeed.square()
            * 0.5f
            * airDensity;
    }

    static vec2f CalculateIdealFlameVector(
        vec2f const & pointVelocity,
        float pointVelocityMagnitudeThreshold) // Threshold after which interpolation is fully towards velocity vector
    {
        // Vector Q is the vector describing the ideal, final flame's
        // direction and (unscaled) length.
        //
        // At rest it's (0, 1) - simply, the flame pointing upwards.
        // When the particle has velocity V, it is the interpolation of the rest upward
        // vector (B) with the opposite of the particle's velocity:
        //      Q = (1-a) * B - a * V
        // Where 'a' depends on the magnitude of the particle's velocity.

        vec2f constexpr B = vec2f(0.0f, 1.0f);

        // The interpolation factor depends on the magnitude of the particle's velocity,
        // via a magic formula; the more the particle's velocity, the more the resultant
        // vector is aligned with the particle's velocity
        float const interpolationFactor = SmoothStep(0.0f, pointVelocityMagnitudeThreshold, pointVelocity.length());
        vec2f Q = B * (1.0f - interpolationFactor) - pointVelocity * interpolationFactor;

        // Magnitude of vector is capped
        float constexpr Qlmax = 1.8f; // Magic number
        float const Ql = Q.length();
        vec2f const Qn = Q.normalise_approx(Ql);
        Q = Qn * std::min(Ql, Qlmax);

        return Q;
    }

    static void EvolveFlameGeometry(
        vec2f & /*in out*/ flameVector,
        float & /*in out*/ flameWindRotationAngle,
        vec2f const & flamePointPosition,
        vec2f const & flamePointVelocity,
        vec2f const & windVelocity, // Km/h
        std::optional<Wind::RadialWindField> const & radialWindField)
    {
        // Vector Q is the vector describing the ideal, final flame's
        // direction and length
        vec2f const Q = Formulae::CalculateIdealFlameVector(
            flamePointVelocity,
            100.0f); // Particle's velocity has a larger impact on the final vector

        // Inertia: converge current flame vector towards target vector Q
        //
        // Convergence rate inversely depends on the magnitude of change:
        // - A big change: little rate (lots of inertia)
        // - A small change: big rate (immediately responsive)
        float constexpr MinFlameVectorConvergenceRate = 0.02f;
        float constexpr MaxFlameVectorConvergenceRate = 0.05f;
        float const flameVectorChangeMagnitude = std::abs(Q.angleCw(flameVector));
        float const flameVectorConvergenceRate =
            MinFlameVectorConvergenceRate
            + (MaxFlameVectorConvergenceRate - MinFlameVectorConvergenceRate) * (1.0f - LinearStep(0.0f, Pi<float>, flameVectorChangeMagnitude));

        flameVector +=
            (Q - flameVector)
            * flameVectorConvergenceRate;

        //
        // Calculate flame wind rotation angle
        //
        // The wind rotation angle has three components:
        //  - Global wind
        //  - Radial wind field, if any
        //  - Particle's velocity
        //
        // We simulate inertia by converging slowly to the target angle.
        //

        vec2f resultantWindSpeedVector =
            windVelocity
            - flamePointVelocity;

        if (radialWindField.has_value())
        {
            vec2f const displacement = flamePointPosition - radialWindField->SourcePos;
            float const radius = displacement.length();
            if (radius < radialWindField->PreFrontRadius)
            {
                resultantWindSpeedVector +=
                    displacement.normalise_approx(radius)
                    * radialWindField->PreFrontWindForceMagnitude
                    * 0.4f; // Magic damper
            }
        }

        // Projection of wind speed vector along flame
        vec2f const flameDir = flameVector.normalise_approx();
        float const windSpeedMagnitudeAlongFlame = resultantWindSpeedVector.dot(flameDir);

        // Our angle moves opposite to the projection of wind along the flame:
        //  - Wind aligned with flame: proj=|W|, angle = 0
        //  - Wind perpendicular to flame: proj=|0|, angle = +/-MAX
        //  - Wind against flame: proj=-|W|, angle = +/-MAX
        float constexpr MaxAngle = 0.27f;
        float const targetFlameWindRotationAngle =
            MaxAngle
            * LinearStep(0.0f, 100.0f, resultantWindSpeedVector.length() - windSpeedMagnitudeAlongFlame)
            * (resultantWindSpeedVector.cross(flameDir) > 0.0f ? -1.0f : 1.0f); // The sign of the angle is positive (CW) when the wind vector is to the right of the flame vector

        // Converge
        float constexpr FlameWindRotationAngleConvergenceRate = 0.055f;
        flameWindRotationAngle +=
            (targetFlameWindRotationAngle - flameWindRotationAngle)
            * FlameWindRotationAngleConvergenceRate;
    }
};

}
