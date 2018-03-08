/*    Copyright (c) 2010-2017, Delft University of Technology
 *    All rigths reserved
 *
 *    This file is part of the Tudat. Redistribution and use in source and
 *    binary forms, with or without modification, are permitted exclusively
 *    under the terms of the Modified BSD license. You should have received
 *    a copy of the license with this file. If not, please or visit:
 *    http://tudat.tudelft.nl/LICENSE.
 */

#include "auxiliaryElements.h"

namespace tudat
{

namespace propagators
{

namespace sst
{

//! Update the instance's members according to the current equinoctialElements.
void AuxiliaryElements::updateMembers()
{
    using namespace tudat::orbital_element_conversions;
    using namespace tudat::mathematical_constants;

    // Properties that depend on other bodies (constant through the integration step)
    // mass = propagatedBody->getMass();
    // area = propagatedBody->getCrossSectionalArea();
    // CD   = propagatedBody->getDragCoefficient();
    // CR   = propagatedBody->getRadiationPressureCoefficient();
    mu = centralGravityAM->getCurrentGravitationalParameter();

    // Retrograde factor
    I = equinoctialElements.isRetrograde() ? -1 : 1;

    // Get equinoctial elements
    a = equinoctialElements( semiMajorAxisIndex );
    h = equinoctialElements( hIndex );
    k = equinoctialElements( kIndex );
    p = equinoctialElements( pIndex );
    q = equinoctialElements( qIndex );

    // Define repeated values
    const double h2 = std::pow( h, 2 );
    const double k2 = std::pow( k, 2 );
    const double p2 = std::pow( p, 2 );
    const double q2 = std::pow( q, 2 );

    // Eccentricity
    e = std::sqrt( h2 + k2 );

    // Check valid state
    if ( a < 0 || e < 0 || e > 1 ) {
        std::cout << equinoctialElements.getComponents( meanType ).transpose() << std::endl;
        throw std::runtime_error( "Equinoctial state is not valid. Spacecraft may be re-entering." );
    }

    // Keplerian mean motion
    meanMotion = std::sqrt( mu / a ) / a;

    // Keplerian period
    period = 2 * PI / meanMotion;

    // Get A, B and C
    A = std::sqrt( mu * a );
    B = std::sqrt( 1 - h2 - k2 );
    C = 1 + p2 + q2;

    // Get equinoctial reference frame vectors
    Eigen::Matrix< double, 3, 3 > fgw =
            getEquinoctialReferenceFrameBasisVectors( p, q, equinoctialElements.isRetrograde() );
    f = fgw.col( 0 );
    g = fgw.col( 1 );
    w = fgw.col( 2 );
}


} // namespace sst

} // namespace propagators

} // namespace tudat
