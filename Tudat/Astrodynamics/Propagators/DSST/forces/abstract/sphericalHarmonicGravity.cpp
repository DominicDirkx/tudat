/*    Copyright (c) 2010-2017, Delft University of Technology
 *    All rigths reserved
 *
 *    This file is part of the Tudat. Redistribution and use in source and
 *    binary forms, with or without modification, are permitted exclusively
 *    under the terms of the Modified BSD license. You should have received
 *    a copy of the license with this file. If not, please or visit:
 *    http://tudat.tudelft.nl/LICENSE.
 */

#include "sphericalHarmonicGravity.h"
#include "Tudat/Astrodynamics/BasicAstrodynamics/orbitalElementConversions.h"

namespace tudat
{

namespace propagators
{

namespace sst
{

namespace force_models
{


//! Update instance's members based on current auxiliary elements.
void SphericalHarmonicGravity::updateMembers( )
{
    CentralBodyPerturbed::updateMembers();
    Conservative::updateMembers();
}


//! Set the values of N and S for the series expansion of the disturbing potential.
void SphericalHarmonicGravity::determineTruncationValues() {
    N = shDegree;
    // Compute the max eccentricity power for the mean element rate expansion
    if ( N == 2 ) {
        S = 0;
    } else {
        // Utilities for truncation
        const double ax2or = 2 * a / R;
        double xmuran = mu / a;
        // Set a lower bound for eccentricity
        const double eo2  = std::max( 0.0025, e / 2.0 );
        const double x2o2 = Chi2 / 2.;
        Vectord eccPwr( N + 1 );
        Vectord chiPwr( N + 1 );
        Vectord hafPwr( N + 1 );
        eccPwr( 0 ) = 1.0;
        chiPwr( 0 ) = Chi;
        hafPwr( 0 ) = 1.0;
        for ( unsigned int i = 1; i <= N; i++ ) {
            eccPwr( i ) = eccPwr( i - 1 ) * eo2;
            chiPwr( i ) = chiPwr( i - 1 ) * x2o2;
            hafPwr( i ) = hafPwr( i - 1 ) * 0.5;
            xmuran /= ax2or;
        }

        // Set highest power of e and degree of current spherical harmonic.
        S = 0;
        int n = N;
        // Loop over n
        do {
            // Set order of current spherical harmonic.
            int m = 0;
            // Loop over m
            do {
                // Compute magnitude of current spherical harmonic coefficient.
        const double Jnm = std::fabs( J( n, m ) );
        if ( Jnm == 0 ) break;
                // Set magnitude of last spherical harmonic term.
                double lastTerm = 0;
                // Set current power of e and related indices.
                int nsld2 = ( n - S - 1 ) / 2;
                int l = n - 2 * nsld2;
                // Loop over l
                double term = 0;
                do {
                    // Compute magnitude of current spherical harmonic term.
                    if ( m < l ) {
            term = Jnm * xmuran * ( factorial( n - l )  / ( factorial( n - m ) ) ) *
                                ( factorial( n + l )  / ( factorial( nsld2 )  *  factorial( nsld2 + l ) ) ) *
                                eccPwr( l ) * upper_bounds::getDnl( Chi2, chiPwr( l ), n, l ) *
                                ( upper_bounds::getRnml(gamma, n, l, m,  1, I) +
                                  upper_bounds::getRnml(gamma, n, l, m, -1, I) );
                    } else {
            term = Jnm * xmuran * ( factorial( n + m ) / ( factorial( nsld2 ) * factorial( nsld2 + l ) ) )
                                * eccPwr( l ) * hafPwr( m - l ) * upper_bounds::getDnl( Chi2, chiPwr( l ), n, l ) *
                                ( upper_bounds::getRnml( gamma, n, m, l,  1, I ) +
                                  upper_bounds::getRnml( gamma, n, m, l, -1, I ) );
                    }

                    // Is the current spherical harmonic term bigger than the truncation tolerance ?
                    if ( term >= TRUNCATION_TOLERANCE ) {
                        S = l;
                    } else {
                        // Is the current term smaller than the last term ?
                        if ( term < lastTerm ) {
                            break;
                        }
                    }
                    // Proceed to next power of e.
                    lastTerm = term;
                    l += 2;
                    nsld2--;
                } while ( l < n );
                // Is the current spherical harmonic term bigger than the truncation tolerance ?
                if ( term >= TRUNCATION_TOLERANCE ) {
                    S = std::min( N - 2, S );
                    return;
                }
                // Proceed to next order.
                m++;
        } while ( m <= std::min( n, (int) shOrder ) );
            // Proceed to next degree.
            xmuran *= ax2or;
            n--;
        } while ( n > (int) S + 2 );

        S = std::min( N - 2, S );
    }

}


} // namespace force_models

} // namespace sst

} // namespace propagators

} // namespace tudat
