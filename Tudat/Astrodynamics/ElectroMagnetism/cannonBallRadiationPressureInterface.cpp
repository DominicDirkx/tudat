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

#include <iostream>

#include "Tudat/Astrodynamics/ElectroMagnetism/cannonBallRadiationPressureInterface.h"


namespace tudat
{

namespace electro_magnetism
{

//! Function to update the current value of the radiation pressure
void CannonBallRadiationPressureInterface::updateInterface(
        const double currentTime )
{

    if( !( currentTime_ == currentTime ) )
    {
        // Calculate current radiation pressure
        currentSolarVector_ = sourcePositionFunction_( ) - targetPositionFunction_( );
        double distanceFromSource = currentSolarVector_.norm( );
        currentRadiationPressure_ = calculateRadiationPressure(
                    sourcePower_( ), distanceFromSource );

        // Calculate total shadowing due to occulting body; note that multiple concurrent
        // occultations are not completely correctly (prints warning).

        updateShadowFunction( currentTime );

        currentRadiationPressure_ *= currentShadowFunction_;;
        radiationPressureCoefficient_ = radiationPressureCoefficientFunction_( currentTime );

        currentTime_ = currentTime;

    }
}

} // namespace electro_magnetism
} // namespace tudat
