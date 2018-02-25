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

#include "Tudat/Astrodynamics/BasicAstrodynamics/missionGeometry.h"
#include "Tudat/Astrodynamics/ElectroMagnetism/radiationPressureInterface.h"


namespace tudat
{

namespace electro_magnetism
{

//! Calculate radiation pressure at certain distance from a source.
double calculateRadiationPressure( const double sourcePower, const double distanceFromSource )
{
    return sourcePower / ( 4.0 * mathematical_constants::PI * distanceFromSource *
                           distanceFromSource * physical_constants::SPEED_OF_LIGHT );
}

} // namespace electro_magnetism
} // namespace tudat
