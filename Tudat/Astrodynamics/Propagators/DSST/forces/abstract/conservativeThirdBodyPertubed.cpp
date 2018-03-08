/*    Copyright (c) 2010-2017, Delft University of Technology
 *    All rigths reserved
 *
 *    This file is part of the Tudat. Redistribution and use in source and
 *    binary forms, with or without modification, are permitted exclusively
 *    under the terms of the Modified BSD license. You should have received
 *    a copy of the license with this file. If not, please or visit:
 *    http://tudat.tudelft.nl/LICENSE.
 */

#include "conservativeThirdBodyPertubed.h"

#include "Tudat/Astrodynamics/BasicAstrodynamics/orbitalElementConversions.h"


namespace tudat
{

namespace propagators
{

namespace sst
{

namespace force_models
{


//! Set the values of N and S for the series expansion of the disturbing potential.
void ConservativeThirdBodyPerturbed::determineTruncationValues() {
    // Get N and S specified by user.
    // If both are at least 2, use them. Otherwise, find right values based on orbit characteristics.
    N = getSettings()->maxN;
    S = getSettings()->maxS;
    if ( N >= 2 && S >= 2 ) {
        return;
    }

    // Truncation tolerance.
    const double aor = a / R3;
    const double tol = ( aor > 0.3 || ( aor > 0.15 && e > 0.25 ) ) ?
		BIG_TRUNCATION_TOLERANCE() : SMALL_TRUNCATION_TOLERANCE();

    // Utilities for truncation
    // Set a lower bound for eccentricity
    const double eo2  = std::max( 0.0025, e / 2 );
    const double x2o2 = Chi2 / 2;
    Vectord eccPwr( MAX_N() );
    Vectord chiPwr( MAX_N() );
    eccPwr( 0 ) = 1.0;
    chiPwr( 0 ) = Chi;
    for ( unsigned int i = 1; i < MAX_N(); i++ ) {
        eccPwr( i ) = eccPwr( i - 1 ) * eo2;
        chiPwr( i ) = chiPwr( i - 1 ) * x2o2;
    }

    // Auxiliary quantities.
    const double ao2rxx = aor / ( 2 * Chi2 );
    double xmuarn       = ao2rxx * ao2rxx / Chi * Ufactor;
    double term         = 0.0;

    // Compute max power for a/R3 and e.
    N                  = 2;
    S                  = 0;
    unsigned int n     = 2;
    unsigned int m     = 2;
    unsigned int nsmd2 = 0;

    do {
        // Upper bound for Tnm.
        term = xmuarn *
                ( factorial( n + m ) / ( factorial( nsmd2 ) * factorial( nsmd2 + m ) ) ) *
                ( factorial( n + m + 1 ) / ( factorial( m ) * factorial( n + 1 ) ) ) *
                ( factorial( n - m + 1 ) / factorial( n + 1 ) ) *
                eccPwr( m ) * upper_bounds::getDnl( Chi2, chiPwr( m ), n + 2, m );

        if ( term < tol ) {
            if ( m == 0 ) {
                break;
            } else if ( m < 2 ) {
                xmuarn *= ao2rxx;
                m = 0;
                n++;
                nsmd2++;
            } else {
                m -= 2;
                nsmd2++;
            }
        } else {
            N = n;
            S = std::max( m, S );
            xmuarn *= ao2rxx;
            m++;
            n++;
        }
    } while ( n < MAX_N() );

    S = std::min( N, S );

    // Ensure N and S are ≥ 2
    N = std::max( (int) N, 2 );
    S = std::max( (int) S, 2 );

    // std::cout << "N = " << N << ", S = " << S << std::endl;
}


//! Function to compute the U derivatives for the current state [ Eq. 3.2-(2) ]
Eigen::Vector6d ConservativeThirdBodyPerturbed::computeMeanDisturbingFunctionPartialDerivatives()
{
    using namespace coefficients_factories;

    // a / R3 up to power maxAR3Pow
    const double aoR3 = a / R3;
    aR3pow = Vectord( N + 1 );
    aR3pow( 0 ) = 1.0;
    for ( unsigned int i = 1; i <= N; i++ ) {
        aR3pow( i ) = aoR3 * aR3pow( i - 1 );
    }

    // Initialise mean potential function
    // U = 0.0;

    // Potential derivatives
    Eigen::Vector6d dU = Eigen::Vector6d::Zero();

    // Vns coefficients
    NestedVectord Vns = VnsFactory.getCoefficientsUpTo( N, S );

    // K0ns coefficients and derivatives with respect to \Chi
    K0nsCoefficientsFactory K0nsFactory( Chi, N, S );
    NestedVectord K0ns = K0nsFactory.getCoefficients();
    NestedVectord dK0ns = K0nsFactory.getDerivatives();

    // Qns coefficients and derivatives with respect to \gamma
    QnsCoefficientsFactory QnsFactory( gamma, N, S );
    NestedVectord Qns  = QnsFactory.getCoefficients();
    NestedVectord dQns = QnsFactory.getDerivatives();

    // Gs coefficients and partial derivatives with respect to h, k, \alpha and \beta
    GsHsCoefficientsFactory GsFactory( h, k, alpha, beta, S );
    Vectord Gs        = GsFactory.getGsCoefficients();
    Vectord Hs        = GsFactory.getHsCoefficients();
    Vectord dGsdh     = GsFactory.getGsHDerivatives();
    Vectord dGsdk     = GsFactory.getGsKDerivatives();
    Vectord dGsdalpha = GsFactory.getGsAlphaDerivatives();
    Vectord dGsdbeta  = GsFactory.getGsBetaDerivatives();

    for ( unsigned int s = 0; s <= S; s++ ) {
        // Get the current Gs coefficient
        const double G = Gs( s );

        // Kronecker symbol ( 2 - delta_{0,s} )
        const double delta0s = ( s == 0 ) ? 1.0 : 2.0;

        for ( unsigned int n = nMin( s ); n <= N; n++ ) {
            // Vns is zero for odd ( n - s ), thus only consider when ( n - s ) is even
            if ( ( n - s ) % 2 == 0 ) {
                // Coefficients for current n and s
                const double K0  =  K0ns( n, s );
                const double dK0 = dK0ns( n, s );

                // Other terms
                const double coef0 = delta0s * aR3pow( n ) * Vns( n, s );
                const double coef1 = coef0 * Qns( n, s );
                const double coef2 = coef1 * K0;

                // Compute U
                // U += coef2 * G;

                Eigen::Vector6d dUterm;
                dUterm << coef2 * n * G,                                  // dU / da
                        coef1 * ( K0 * dGsdh( s ) + hChi3 * G * dK0 ),    // dU / dh
                        coef1 * ( K0 * dGsdk( s ) + kChi3 * G * dK0 ),    // dU / dk
                        coef2 * dGsdalpha( s ),                           // dU / dalpha
                        coef2 * dGsdbeta( s ),                            // dU / dbeta
                        coef0 * K0 * dQns( n, s ) * G;                    // dU / dgamma

                dU += dUterm;
            }
        }
    }

    // Multiply by constant term outside of summation
    // U *= Ufactor;
    dU *= Ufactor;
    dU( 0 ) /= a;

    // Return U derivatives
    return dU;
}


//! Get the short period terms for the current auxiliary elements [ Sections 2.5.2 and 4.2 ]
Eigen::Vector6d ConservativeThirdBodyPerturbed::computeShortPeriodTerms( )
{
    return Eigen::Vector6d::Zero();
}


} // namespace force_models

} // namespace sst

} // namespace propagators

} // namespace tudat
