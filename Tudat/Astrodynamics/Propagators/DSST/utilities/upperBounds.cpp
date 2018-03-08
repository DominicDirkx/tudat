/*    Copyright (c) 2010-2017, Delft University of Technology
 *    All rigths reserved
 *
 *    This file is part of the Tudat. Redistribution and use in source and
 *    binary forms, with or without modification, are permitted exclusively
 *    under the terms of the Modified BSD license. You should have received
 *    a copy of the license with this file. If not, please or visit:
 *    http://tudat.tudelft.nl/LICENSE.
 */

#include "upperBounds.h"


namespace tudat
{

namespace propagators
{

namespace sst
{

namespace upper_bounds
{


//! Get the upper bound value D^n_l(\Chi)
double getDnl( const double xx, const double xpl, const int n, const int l ) {
    const int lp2 = l + 2;
    if ( n > lp2 ) {
        const int ll = l * l;
        double dM = xpl;
        double dL = dM;
        double dB = ( l + 1 ) * xx * xpl;
        for ( int j = l + 3; j <= n; j++ ) {
            const int jm1 = j - 1;
            dL = dM;
            dM = dB;
            dB = jm1 * xx * ( ( 2 * j - 3 ) * dM - ( j - 2 ) * dL) / ( jm1 * jm1 - ll );
        }
        return dB;
    } else if ( n == lp2 ) {
        return  ( l + 1 ) * xx * xpl;
    } else {
        return xpl;
    }
}


//! Get the upper bound value R^ε_{n,m,l}(γ)
double getRnml(const double gamma, const int n, const int l, const int m, const int eps, const int irf) {
    // Initialization
    const int mei = m * eps * irf;
    const double sinisq = 1. - gamma * gamma;
    // Set a lower bound for inclination
    const double sininc = std::max( 0.03, std::sqrt( sinisq ) );
    const double onepig = 1. + gamma * irf;
    const double sinincPowLmMEI = std::pow( sininc, l - mei );
    const double onepigPowLmMEI = std::pow( onepig, mei );

    // Bound for index 0
    double rBound = sinincPowLmMEI * onepigPowLmMEI;

    // If index > 0
    if ( n > l ) {
        const int lp1 = l + 1;

        double dpnml  = lp1 * eps;
        double pnml   = dpnml * gamma - m;

        // If index > 1
        if ( n > l + 1 ) {
            const int ll  = l * l;
            const int ml  = m * l;
            const int mm  = m * m;

            double pn1ml  = 1.;
            double dpn1ml = 0.;
            double pn2ml  = 1.;
            double dpn2ml = 0.;
            for ( int in = l + 2; in <= n; in++ ) {
                const int nm1   = in - 1;
                const int tnm1  = in + nm1;
                const int nnlnl = nm1 * (in * in - ll);
                const int nnmnm = in * (nm1 * nm1 - mm);
                const int c2nne = tnm1 * in * nm1 * eps;
                const int c2nml = tnm1 * ml;
                const double coef = c2nne * gamma - c2nml;

                pn2ml  = pn1ml;
                dpn2ml = dpn1ml;
                pn1ml  = pnml;
                dpn1ml = dpnml;
                pnml   = (coef * pn1ml  - nnmnm * pn2ml) / nnlnl;
                dpnml  = (coef * dpn1ml - nnmnm * dpn2ml + c2nne * pn1ml) / nnlnl;
            }
        }
        // Bound for index > 0
        rBound *= std::sqrt( pnml * pnml + dpnml * dpnml * sinisq / ((n - l) * (n + lp1)) );
    }

    return rBound;
}


} // namespace upper_bounds

} // namespace sst

} // namespace propagators

} // namespace tudat
