/*    Copyright (c) 2010-2017, Delft University of Technology
 *    All rigths reserved
 *
 *    This file is part of the Tudat. Redistribution and use in source and
 *    binary forms, with or without modification, are permitted exclusively
 *    under the terms of the Modified BSD license. You should have received
 *    a copy of the license with this file. If not, please or visit:
 *    http://tudat.tudelft.nl/LICENSE.
 */

#include "nonConservative.h"
#include "boost/bind.hpp"


namespace tudat
{

namespace propagators
{

namespace sst
{

namespace force_models
{


//! Update instance's members that are computed from the current auxiliary elements.
void NonConservative::updateMembers( )
{
    using namespace mathematical_constants;

    determineIntegrationLimits();

    // std::cout << "numberOfQuadratureNodes: " << getSettings()->numberOfQuadratureNodes << std::endl;
    // std::cout << "numberOfQuadratureNodesIsScalable: " << getSettings()->numberOfQuadratureNodesIsScalable << std::endl;
    // std::cout << "alwaysIncludeCentralNode: " << getSettings()->alwaysIncludeCentralNode << std::endl;

    // Determine number of nodes for the Gaussian quadrature
    N = getSettings()->numberOfQuadratureNodes;
    if ( getSettings()->numberOfQuadratureNodesIsScalable ) {
        N = std::round( N * ( L2 - L1 ) / ( 2 * PI ) );
        // Make N odd
        if ( N % 2 == 0 && getSettings()->alwaysIncludeCentralNode ) {
            N += 1;
        }
    }

    // std::cout << N << " nodes, average spacing of " << ( L2 - L1 ) / N * 180 / PI << " deg." << std::endl;
}


//! Update the acceleration model to the current epoch and get the perturbing acceleration
Eigen::Vector3d NonConservative::getPerturbingAcceleration()
{
    // Create current global Cartesian state from r and v, and central body state
    Eigen::Matrix< double, Eigen::Dynamic, 1 > cartesianState( 6 );
    cartesianState << r, v;
    cartesianState += aux.centralBody->getStateInBaseFrameFromEphemeris( epoch );

    // Update environment
    environmentUpdater->updateEnvironment( epoch, { { transational_state, cartesianState } } );
    accelerationModel->updateMembers( epoch );

    // Get acceleration
    return accelerationModel->getAcceleration();
}


//! Update r, v, X, Y, Xdot and Ydot and partial derivatvies of the equinoctial elements wrt to v.
Eigen::Vector6d NonConservative::integrand( const double trueLongitude )
{
    using namespace Eigen;
    using namespace mathematical_constants;
    using namespace orbital_element_conversions;

    // True longitude
    L = trueLongitude;

    // Get epoch
    Vector6d components = aux.equinoctialElements.getComponents( meanType );
    const double auxMeanLongitude = components( fastVariableIndex );
    components( fastVariableIndex ) = L;
    const double meanLongitude = getMeanLongitudeFromTrue( components );
    epoch = aux.epoch + ( meanLongitude - auxMeanLongitude ) / meanMotion;
    // std::cout << std::setprecision( 10 ) << epoch << std::endl;

    // Compute other things...
    const double sinL = std::sin( L );
    const double cosL = std::cos( L );

    // R  [ Eq. 2.1.4-(6) ]
    R = a * B * B / ( 1 + h * sinL + k * cosL );

    // X, Y  [ Eq. 2.1.4-(7) ]
    X = R * cosL;
    Y = R * sinL;

    // Xdot, Ydot  [ Eq. 2.1.4-(8) ]
    Xdot = -meanMotion * a / B * ( h + sinL );
    Ydot =  meanMotion * a / B * ( k + cosL );

    // Get f, g, w, µ
    Vector3d &f = aux.f;
    Vector3d &g = aux.g;
    Vector3d &w = aux.w;
    const double mu = aux.mu;

    // r, v  [ Eq. 2.1.4-(9) ]
    r = X    * f + Y    * g;
    v = Xdot * f + Ydot * g;

    // std::cout << aux.epoch << ": r = \n" << r.norm() << "\n" << std::endl;

    // V
    // V = v.norm();

    // Partial derivatives of the equinoctial elements wrt to v  [ Eq. 2.1.7-(3) ]
    const Vector3d wfirst = ( I * q * Y - p * X ) / A * w;
    const Vector3d wsecond = C / ( 2 * A * B ) * w;
    const double XdotY = Xdot * Y;
    const double XYdot = X * Ydot;
    partials.col( 0 ) =  2 / ( std::pow( meanMotion, 2 ) * a ) * v;                                          // da/dv
    partials.col( 1 ) =  ( ( 2 * XdotY - XYdot ) * f - X * Xdot * g ) / mu + k / B * wfirst;                 // dh/dv
    partials.col( 2 ) =  ( ( 2 * XYdot - XdotY ) * g - Y * Ydot * f ) / mu - h / B * wfirst;                 // dk/dv
    partials.col( 3 ) =      Y * wsecond;                                                                    // dp/dv
    partials.col( 4 ) =  I * X * wsecond;                                                                    // dq/dv
    partials.col( 5 ) = -2 / A * r + ( k * partials.col( 1 ) - h * partials.col( 2 ) ) / ( 1 + B ) + wfirst; // dλ/dv

    // Disturbing acceleration
    perturbingAcceleration = getPerturbingAcceleration();

    // Compute vector of integrands:               6x3 * 3x1                    ->  6x1  ( Vector6d )
    return std::pow( R / a, 2 ) * partials.transpose() * perturbingAcceleration;
}


//! Get the mean element rates for the current auxiliary elements [ Eq. 3.4-(1b) ]
Eigen::Vector6d NonConservative::computeMeanElementRates( )
{
    using namespace Eigen;
    using namespace numerical_quadrature;
    using namespace mathematical_constants;

    if ( L1 >= L2 ) { // Return immediately if drag can be completely neglected
        return Vector6d::Zero();
    }

    // Solve the integral by Gaussian quadrature
    aux.gaussianQuadrature.reset( boost::bind( &NonConservative::integrand, this, _1 ), L1, L2, N );
    Vector6d integralResult = aux.gaussianQuadrature.getQuadrature();
    return integralResult / ( 2 * PI * B );
}


//! Get the short period terms for the current auxiliary elements [ Section 4.4 ]
Eigen::Vector6d NonConservative::computeShortPeriodTerms( )
{
    // Numerical quadrature
    // GaussianQuadrature gaussSin( integrandSin, L1, L2, N );
    // GaussianQuadrature gaussCos( integrandCos, L1, L2, N );

    return Eigen::Vector6d::Zero();
}


} // namespace force_models

} // namespace sst

} // namespace propagators

} // namespace tudat
