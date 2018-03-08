/*    Copyright (c) 2010-2017, Delft University of Technology
 *    All rigths reserved
 *
 *    This file is part of the Tudat. Redistribution and use in source and
 *    binary forms, with or without modification, are permitted exclusively
 *    under the terms of the Modified BSD license. You should have received
 *    a copy of the license with this file. If not, please or visit:
 *    http://tudat.tudelft.nl/LICENSE.
 */

#include "conservativeRadiationPressure.h"


namespace tudat
{

namespace propagators
{

namespace sst
{

namespace force_models
{

//! Update instance's members that are computed from the current auxiliary elements.
void ConservativeRadiationPressure::updateMembers( )
{
    // Update vector to source
    r3 = radiationPressureAM->getCurrentDistanceToSource() * radiationPressureAM->getCurrentVectorToSource();

    // std::cout << "SRP: " << r3.transpose() << std::endl;

    ConservativeThirdBodyPerturbed::updateMembers();

    // std::cout << "SRP: S = " << S << ", N = " << N << std::endl;
}


//! Function that updates `Ufactor`
double ConservativeRadiationPressure::computeUfactor()
{
    using namespace mathematical_constants;
    using namespace physical_constants;

    // Get current values for parameters
    const double power = radiationPressureAM->getRadiationPressureInterface()->getSourcePowerFunction()();
    const double mass  = radiationPressureAM->getCurrentMass();
    const double area  = radiationPressureAM->getCurrentArea();
    const double CR    = radiationPressureAM->getCurrentRadiationPressureCoefficient();

    // Eq. 3.5-(7)
    const double T = CR * area / mass * power / ( 4 * PI * SPEED_OF_LIGHT );

    // Eq. 3.5-(10)
    return - T / R3;
}


} // namespace force_models

} // namespace sst

} // namespace propagators

} // namespace tudat
