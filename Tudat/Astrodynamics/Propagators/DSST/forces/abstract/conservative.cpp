/*    Copyright (c) 2010-2017, Delft University of Technology
 *    All rigths reserved
 *
 *    This file is part of the Tudat. Redistribution and use in source and
 *    binary forms, with or without modification, are permitted exclusively
 *    under the terms of the Modified BSD license. You should have received
 *    a copy of the license with this file. If not, please or visit:
 *    http://tudat.tudelft.nl/LICENSE.
 */

#include "conservative.h"

namespace tudat
{

namespace propagators
{

namespace sst
{

namespace force_models
{


//! Update instance's members that are computed from the current auxiliary elements.
void Conservative::updateMembers( )
{
    //\Chi^{-2}
    B2 = B * B;
    //\Chi^{-3}
    B3 = B * B2;

    // \Chi
    Chi = 1.0 / B;
    Chi2 = Chi * Chi;
    Chi3 = Chi * Chi2;
    // -2 * a / A
    m2aoA = -2.0 * a / A;
    // B / A
    BoA = B / A;
    // 1 / AB
    ooAB = 1.0 / (A * B);
    // -C / 2AB
    mCo2AB = -C * ooAB / 2.0;
    // B / A(1 + B)
    BoABpo = BoA / (1.0 + B);
    // h * \Chi³
    hChi3 = h * Chi3;
    // k * \Chi³
    kChi3 = k * Chi3;
}


//! Get the mean element rates for the current auxiliary elements [ Eq. 3.1-(1) ]
Eigen::Vector6d Conservative::computeMeanElementRates( )
{
    // Compute potential U derivatives
    const Eigen::Vector6d dU = computeMeanDisturbingFunctionPartialDerivatives();
    const double dUda  = dU( 0 );
    const double dUdh  = dU( 1 );
    const double dUdk  = dU( 2 );
    const double dUdAl = dU( 3 );
    const double dUdBe = dU( 4 );
    const double dUdGa = dU( 5 );

    // Compute cross derivatives [Eq. 2.2-(8)]
    // U(alpha,gamma) = alpha * dU/dgamma - gamma * dU/dalpha
    const double UAlphaGamma   = alpha * dUdGa - gamma * dUdAl;
    // U(beta,gamma) = beta * dU/dgamma - gamma * dU/dbeta
    const double UBetaGamma    =  beta * dUdGa - gamma * dUdBe;
    // Common factor
    const double pUAGmIqUBGoAB = ( p * UAlphaGamma - I * q * UBetaGamma ) * ooAB;

    // Compute mean elements rates [Eq. 3.1-(1)]
    Eigen::Vector6d meanElementRates;
    meanElementRates << 0.0,                                                    // da / dt
            BoA * dUdk + k * pUAGmIqUBGoAB,                                     // dh / dt
            -BoA * dUdh - h * pUAGmIqUBGoAB,                                    // dk / dt
            mCo2AB * UBetaGamma,                                                // dp / dt
            mCo2AB * UAlphaGamma * I,                                           // dq / dt
            m2aoA * dUda + BoABpo * ( h * dUdh + k * dUdk ) + pUAGmIqUBGoAB;    // d\lambda / dt

    return meanElementRates;
}


} // namespace force_models

} // namespace sst

} // namespace propagators

} // namespace tudat
