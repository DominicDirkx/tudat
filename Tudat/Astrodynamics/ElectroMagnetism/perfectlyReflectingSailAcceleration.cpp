/*    Copyright (c) 2010-2018, Delft University of Technology
 *    All rigths reserved
 *
 *    This file is part of the Tudat. Redistribution and use in source and
 *    binary forms, with or without modification, are permitted exclusively
 *    under the terms of the Modified BSD license. You should have received
 *    a copy of the license with this file. If not, please or visit:
 *    http://tudat.tudelft.nl/LICENSE.
 *
 */

#include "Tudat/Mathematics/BasicMathematics/mathematicalConstants.h"

#include "Tudat/Astrodynamics/ElectroMagnetism/cannonBallRadiationPressureAcceleration.h"
//#include "Tudat/Astrodynamics/ElectroMagnetism/cannonBallRadiationPressureForce.h"

namespace tudat
{
namespace electro_magnetism
{

//! Compute radiation pressure acceleration using a perfectly reflecting sail model (ideal sail).
//! Returns acceleration in Heliocentric Orbital Reference Frame
Eigen::Vector3d computePerfectReflectionSailAcceleration(
        const double sailLightnessNumber,
        const double distanceToSource,
        const double sourceBodyGravitationalParameter,
		const Eigen::Vector2d sailControlAngles	)
{
	Eigen::Vector3d sailAcceleration;
	const double commonValue = sailLightnessNumber*(sourceBodyGravitationalParameter/std::pow(distanceToSource, 2.0))*std::pow(std::cos(sailControlAngles(0)),2.0);
	sailAcceleration(0) = commonValue * std::cos(sailControlAngles(0)); //Radial direction (r)
	sailAcceleration(1) = commonValue * std::sin(sailControlAngles(0)) * std::sin(sailControlAngles(1)); // In the orbital plane, normal to radial direction (v)
	sailAcceleration(2) = commonValue * std::sin(sailControlAngles(0)) * std::cos(sailControlAngles(1)); // Perpendicular to orbital plane (h)
    return sailAcceleration
}

} // namespace electro_magnetism
} // namespace tudat
