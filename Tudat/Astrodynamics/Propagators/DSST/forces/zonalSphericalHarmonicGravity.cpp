/*    Copyright (c) 2010-2017, Delft University of Technology
 *    All rigths reserved
 *
 *    This file is part of the Tudat. Redistribution and use in source and
 *    binary forms, with or without modification, are permitted exclusively
 *    under the terms of the Modified BSD license. You should have received
 *    a copy of the license with this file. If not, please or visit:
 *    http://tudat.tudelft.nl/LICENSE.
 */

#include "zonalSphericalHarmonicGravity.h"

namespace tudat
{

namespace propagators
{

namespace sst
{

namespace force_models
{


//! Transform mean elements to osculating elements using only the contribution of J2
Eigen::Vector6d transformToOsculatingUsingJ2Contribution( const Eigen::Vector6d keplerianElements,
                                                          const double gravitationalParameter,
                                                          const double J2, const double referenceRadius )
{
    using namespace Eigen;
    using namespace orbital_element_conversions;
    using namespace mathematical_constants;

    // Only J2
    // Wakker Section 23.4 (p. 632)

    // std::cout << "EQU:  " << aux.equinoctialElements.getComponents( meanType ).transpose() << std::endl;
    // std::cout << "KEP:  " << keplerianElements.transpose() << std::endl;
    const double a =    keplerianElements( semiMajorAxisIndex );                   // mean semimajor axis
    const double e =    keplerianElements( eccentricityIndex );                    // mean eccentricity
    const double i =    keplerianElements( inclinationIndex );                     // mean inclination
    const double w =    keplerianElements( argumentOfPeriapsisIndex );             // mean argument of perigee
    const double raan = keplerianElements( longitudeOfAscendingNodeIndex );        // mean RAAN
    const double f =    keplerianElements( trueAnomalyIndex );                     // mean true anomaly
    const double E =    convertTrueAnomalyToEllipticalEccentricAnomaly( f, e );    // mean eccentric anomaly
    const double M =    convertEllipticalEccentricAnomalyToMeanAnomaly( E, e );    // mean mean anomaly

    // Eq. (23.41)
    // Repeated terms
    const double ee = std::pow( e, 2 );
    const double sqrt1ee = std::sqrt( 1 - ee );
    const double p = a * ( 1 - std::pow( e, 2 ) );
    const double JR_p2 = J2 * std::pow( referenceRadius / p, 2 );
    const double JR2_p = J2 * std::pow( referenceRadius, 2 ) / p;
    const double u = w + f;
    const double uu = u + u;
    const double wwf = u + w;
    const double wwfff = wwf + f + f;
    const double sinuu    = std::sin( uu );
    const double sinwwf   = std::sin( wwf );
    const double sinwwfff = std::sin( wwfff );
    const double cosuu    = std::cos( uu );
    const double coswwf   = std::cos( wwf );
    const double coswwfff = std::cos( wwfff );
    const double ip_bracket = cosuu + e * coswwf + e / 3 * coswwfff;
    const double sin2i = std::pow( std::sin( i ), 2 );
    const double cosi = std::cos( i );
    const double cos2i = std::pow( cosi, 2 );
    const double esinf = e * std::sin( f );
    const double sqrt_mu_p = std::sqrt( gravitationalParameter / p );
    const double ecosf = e * std::cos( f );
    const double ecosf1 = ecosf + 1;
    const double tcos2i1 = 3 * cos2i - 1;
    const double sqrt1ee1 = sqrt1ee + 1;
    const double ecosf12 = std::pow( ecosf1, 2 );
    const double uterm = - tcos2i1 / sqrt1ee1;

    // Short-period terms in auxiliary element set
    const double i_sh = 0.375 * JR_p2 * std::sin( 2 * i ) * ip_bracket;
    const double p_sh = 1.5 * JR2_p * sin2i * ip_bracket;
    const double raan_sh = -1.5 * JR_p2 * cosi * ( f - M + esinf - 0.5 * sinuu - 0.5 * e * sinwwf - e / 6 * sinwwfff );
    const double r_sh = -0.25 * JR2_p * ( tcos2i1 * ( 2 * sqrt1ee / ecosf1 + ecosf / sqrt1ee1 + 1 ) - sin2i * cosuu );
    const double rdot_sh = 0.25 * sqrt_mu_p * JR_p2 *
            ( tcos2i1 * esinf * ( sqrt1ee + ecosf12 / sqrt1ee1 ) - 2 * sin2i * ecosf12 * sinuu );
    const double u_sh = -0.125 * JR_p2 * ( 6 * ( 1 - 5 * cos2i ) * ( f - M )
                                           + 4 * esinf * ( ( 1 - 6 * cos2i ) + uterm )
                                           + uterm * ee * std::sin( 2 * f )
                                           + 2 * e * ( 5 * cos2i - 2 ) * sinwwf
                                           + ( 7 * cos2i - 1 ) * sinuu
                                           + 2 * e * cos2i * sinwwfff );

    // Mean distance
    const double r = p / ecosf1;

    // Mean radial speed: http://www.physicspages.com/2015/04/24/velocity-in-an-elliptical-orbit/
    const double rdot = sqrt_mu_p * esinf;

    // Osculating elements in auxiliary element set
    const double i_osc = i + i_sh;
    const double p_osc = p + p_sh;
    const double raan_osc = raan + raan_sh;
    const double r_osc = r + r_sh;
    const double rdot_osc = rdot + rdot_sh;
    const double u_osc = u + u_sh;

    // Osculating Keplerian elements
    const double C1 = p_osc / r_osc - 1;
    const double C2 = std::sqrt( p_osc / gravitationalParameter ) * rdot_osc;
    const double e_osc = std::sqrt( std::pow( C1, 2 ) + std::pow( C2, 2 ) );
    const double f_osc = std::atan2( C2, C1 );
    const double a_osc = p_osc / ( 1 - std::pow( e_osc, 2 ) );
    const double w_osc = u_osc - f_osc;

    Vector6d oscKeplerianElements;
    oscKeplerianElements( semiMajorAxisIndex            ) = a_osc;
    oscKeplerianElements( eccentricityIndex             ) = e_osc;
    oscKeplerianElements( inclinationIndex              ) = i_osc;
    oscKeplerianElements( argumentOfPeriapsisIndex      ) = w_osc;
    oscKeplerianElements( longitudeOfAscendingNodeIndex ) = raan_osc;
    oscKeplerianElements( trueAnomalyIndex              ) = f_osc;

    return oscKeplerianElements;
}


//! Function to compute the U derivatives for the current state [ Eq. 3.1-(6) ]
Eigen::Vector6d ZonalSphericalHarmonicGravity::computeMeanDisturbingFunctionPartialDerivatives()
{
    using namespace coefficients_factories;

    // R / a up to power maxAR3Pow
    const double R_a = R / a;
    powR_a = Vectord( N + 1 );
    powR_a( 0 ) = 1.0;
    for ( unsigned int i = 1; i <= N; i++ ) {
        powR_a( i ) = R_a * powR_a( i - 1 );
    }

    // Initialise mean potential function
    // U = 0.0;

    // Potential derivatives
    Eigen::Vector6d dU = Eigen::Vector6d::Zero();

    // Vns coefficients
    NestedVectord Vns = VnsFactory.getCoefficientsUpTo( N, S );

    // K_0{-n-1,s} coefficients and derivatives with respect to \Chi
    NegativeK0nsCoefficientsFactory mK0nsFactory( Chi, N + 1, S );
    NestedVectord  mK0ns = mK0nsFactory.getCoefficients();
    NestedVectord dmK0ns = mK0nsFactory.getDerivatives();

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

        for ( unsigned int n = s + 2; n <= N; n++ ) {
            // Vns is zero for odd ( n - s ), thus only consider when ( n - s ) is even
            if ( ( n - s ) % 2 == 0 ) {
                // Coefficients for current n and s
                const double  mK0 =  mK0ns( n + 1, s );
                const double dmK0 = dmK0ns( n + 1, s );

                // Other terms
                const double R_an  = powR_a( n );
                const double coef0 = delta0s * R_an * J( n ) * Vns( n, s );
                const double coef1 = coef0 * Qns( n, s );
                const double coef2 = coef1 * mK0;
                const double coef3 = coef2 * G;

                // Compute U
                // U += coef3;

                Eigen::Vector6d dUterm;
                dUterm << coef3 * ( n + 1 ),                                // dU / da
                        coef1 * ( mK0 * dGsdh( s ) + hChi3 * G * dmK0 ),    // dU / dh
                        coef1 * ( mK0 * dGsdk( s ) + kChi3 * G * dmK0 ),    // dU / dk
                        coef2 * dGsdalpha( s ),                             // dU / dalpha
                        coef2 * dGsdbeta( s ),                              // dU / dbeta
                        coef0 * mK0 * dQns( n, s ) * G;                     // dU / dgamma

                dU += dUterm;
            }
        }
    }

    // Multiply by constant term outside of summation
    const double term = -mu / a;
    // U *= mu_a;
    dU *= term;
    dU( 0 ) /= -a;

    // Return U derivatives
    return dU;
}


//! Get the short period terms for the current auxiliary elements
Eigen::Vector6d ZonalSphericalHarmonicGravity::computeShortPeriodTerms( )
{
    // Only J2

    using namespace Eigen;
    using namespace orbital_element_conversions;

    Vector6d shortPeriodTerms = Vector6d::Zero();

    // Get osculating Keplerian elements
    Vector6d oscKeplerianElements = transformToOsculatingUsingJ2Contribution( aux.equinoctialElements.toKeplerian(),
                                                                              aux.mu, J( 2 ), R );

    // Osculating equinoctial elements
    Vector6d oscEquinoctialElements = convertKeplerianToEquinoctialElements(
                oscKeplerianElements, meanType, aux.equinoctialElements.isRetrograde() );

    // std::cout << "SP kep: " << ( oscKeplerianElements - keplerianElements ).transpose() << std::endl;

    // Short-period terms in equinoctial elements
    shortPeriodTerms = oscEquinoctialElements - aux.equinoctialElements.getComponents( meanType );

    // std::cout << "SP equ: " << shortPeriodTerms.transpose() << "\n" << std::endl;

    return shortPeriodTerms;
}


} // namespace force_models

} // namespace sst

} // namespace propagators

} // namespace tudat
