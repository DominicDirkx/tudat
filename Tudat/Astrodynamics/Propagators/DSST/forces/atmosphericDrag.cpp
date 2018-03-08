/*    Copyright (c) 2010-2017, Delft University of Technology
 *    All rigths reserved
 *
 *    This file is part of the Tudat. Redistribution and use in source and
 *    binary forms, with or without modification, are permitted exclusively
 *    under the terms of the Modified BSD license. You should have received
 *    a copy of the license with this file. If not, please or visit:
 *    http://tudat.tudelft.nl/LICENSE.
 */

#include "atmosphericDrag.h"

namespace tudat
{

namespace propagators
{

namespace sst
{

namespace force_models
{


//! Update instance's members that are computed from the current auxiliary elements.
void AtmosphericDrag::updateMembers( )
{
    // Update direction cosines
    CentralBodyPerturbed::updateMembers();

    // Update integration limits
    NonConservative::updateMembers();
}


//! Update the values of the minimum and maximum true longitude for the averaging integral.
void AtmosphericDrag::determineIntegrationLimits( ) {
    using namespace mathematical_constants;

    const double hmax = getSettings()->altitudeLimit;
    if ( hmax <= 0.0 ) {
        // No altitude limit specified, thus consider drag thorughout all the revolution
        L1 = -PI;
        L2 =  PI;
        return;
    }

    // Eq. 3.4-(2b)
    const double cosflim = ( a * ( 1 - std::pow( e, 2 ) )
                             / ( aux.centralBody->getShapeModel()->getAverageRadius() + hmax ) - 1 ) / e;
    if ( cosflim > 1 ) { // Drag can be completely neglected
        L1 = 0;
        L2 = 0;
        return;
    } else if ( cosflim < -1 ) { // Drag has to be considered throughout all the revolution
        L1 = -PI;
        L2 =  PI;
        return;
    }
    // Otherwise, get the limits:

    const double flim = std::acos( cosflim );

    // longitude of perigee = omega + I * Omega
    // Eq. 2.1.2-(1)
    const double longitudePerigeeSin = h / e;
    const double longitudePerigeeCos = k / e;
    const double longitudePerigee = std::atan2( longitudePerigeeSin, longitudePerigeeCos );

    // Eq. 3.4-(2b)
    L1 = longitudePerigee - flim;
    L2 = longitudePerigee + flim;
}


//! Get the mean element rates for the current auxiliary elements [ Eq. 3.4-(1b) ]
Eigen::Vector6d AtmosphericDrag::computeMeanElementRates( )
{
    using namespace Eigen;
    using namespace mathematical_constants;

    if ( N > 1 )  // use Gaussian Quadrature
    {
        return NonConservative::computeMeanElementRates();
    }
    else  // evaluate integrand only at perigee and use empirical polynomial fit based on perigee altitude (and F10.7)
    {
        Vector6d perigeeIntegralResult = integrand( 0.5 * ( L2 + L1 ) ) * ( L2 - L1 );
        Vector6d perigeeMeanElementRates = perigeeIntegralResult / ( 2 * PI * B );

        // Get perigee altitude
        const double hp = a * ( 1 - e ) - aux.centralBody->getShapeModel()->getAverageRadius();

        // Get F10.7 (if using NRLMSISE00 atmosphere)
        double f107 = TUDAT_NAN;
        if ( atmosphereModel != NULL ) {
            f107 = atmosphereModel->getNRLMSISE00Input().f107;
        }

        double factor;
        if ( f107 > 40 ) {
            factor = ( 2e-07 * hp + 0.115 ) * std::pow( f107 - 40, 0.1 ) + 0.0534;
        } else {
            factor = 3.72e-07 * hp + 0.208;
        }

        return factor * perigeeMeanElementRates;
    }
}


} // namespace force_models

} // namespace sst

} // namespace propagators

} // namespace tudat
