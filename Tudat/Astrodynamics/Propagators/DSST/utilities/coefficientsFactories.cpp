/*    Copyright (c) 2010-2017, Delft University of Technology
 *    All rigths reserved
 *
 *    This file is part of the Tudat. Redistribution and use in source and
 *    binary forms, with or without modification, are permitted exclusively
 *    under the terms of the Modified BSD license. You should have received
 *    a copy of the license with this file. If not, please or visit:
 *    http://tudat.tudelft.nl/LICENSE.
 */

#include "coefficientsFactories.h"

namespace tudat
{

namespace propagators
{

namespace sst
{

//! Factorials generator declared at a global scope.
coefficients_factories::FactorialsFactory factorial;


namespace coefficients_factories
{

///////////////////////////////////////////////////////////////////////////////
///////////////////////////// FactorialsFactory ///////////////////////////////
///////////////////////////////////////////////////////////////////////////////

//! 0! = 1
void FactorialsFactory::initializeCoefficients() {
    factorials = { 1.0 };
}

//! x! = x * ( x - 1 )!
void FactorialsFactory::computeCoefficients( const bool extended ) {
    if ( order > maximumComputedOrder ) {
        factorials.resizeForOrder( order );
        for ( unsigned int i = maximumComputedOrder + 1; i <= order; i++ ) {
            factorials( i ) = i * factorials( i - 1 );
            // std::cout << i << "! = " << factorials( i ) << std::endl;
        }
    }
}


///////////////////////////////////////////////////////////////////////////////
/////////////////////////// GsHsCoefficientsFactory ///////////////////////////
///////////////////////////////////////////////////////////////////////////////

//! From Eq. 3.1-(5)
void GsHsCoefficientsFactory::initializeCoefficients() {
    Gs = { 1.0 };
    Hs = { 0.0 };
}

//! From Orekit
void GsHsCoefficientsFactory::initializeDependentCoefficients() {
    dGs = { { "h",     { 0.0 } },
            { "k",     { 0.0 } },
            { "alpha", { 0.0 } },
            { "beta",  { 0.0 } } };
}

//! From Eq. 3.1-(5)
void GsHsCoefficientsFactory::computeCoefficients( const bool extended ) {
    const unsigned int order = getOrder( extended );
    if ( order > maximumComputedOrder ) {
        Gs.resizeForOrder( order );
        Hs.resizeForOrder( order );
        const double h     = getParameter( "h" );
        const double k     = getParameter( "k" );
        const double alpha = getParameter( "alpha" );
        const double beta  = getParameter( "beta" );
        const double hamkb = h * alpha - k * beta;
        const double kaphb = k * alpha + h * beta;
        for ( unsigned int s = maximumComputedOrder + 1; s <= order; s++ ) {
            Gs( s ) = kaphb * Gs( s - 1 ) - hamkb * Hs( s - 1 );
            Hs( s ) = hamkb * Gs( s - 1 ) + kaphb * Hs( s - 1 );
        }
    }
}

//! From Eq. 3.1-(9)
void GsHsCoefficientsFactory::computeDependentCoefficients() {
    if ( order > maximumComputedOrderDependentCoefficients ) {
        // Resize Gs derivatives
        for ( auto &ent : dGs ) {
            ent.second.resizeForOrder( order );
        }
        // Compute extra coefficients
        const double h     = getParameter( "h" );
        const double k     = getParameter( "k" );
        const double alpha = getParameter( "alpha" );
        const double beta  = getParameter( "beta" );
        for ( unsigned int s = maximumComputedOrderDependentCoefficients + 1; s <= getOrder(); s++ ) {
            dGs[ "h"     ]( s ) = s * beta  * Gs( s - 1 ) - s * alpha * Hs( s - 1 );
            dGs[ "k"     ]( s ) = s * alpha * Gs( s - 1 ) + s * beta  * Hs( s - 1 );
            dGs[ "alpha" ]( s ) = s * k     * Gs( s - 1 ) - s * h     * Hs( s - 1 );
            dGs[ "beta"  ]( s ) = s * h     * Gs( s - 1 ) + s * k     * Hs( s - 1 );
        }
    }
}


///////////////////////////////////////////////////////////////////////////////
/////////////////////////// QnsCoefficientsFactory ////////////////////////////
///////////////////////////////////////////////////////////////////////////////

//! From Eq. 2.8.3-(3)
void QnsCoefficientsFactory::initializeCoefficients() {
    Qns = { { 1.0 } };
}

void QnsCoefficientsFactory::initializeDependentCoefficients() {
    dQns = DoubleOrderCoefficients();
}

unsigned int QnsCoefficientsFactory::getSuborder( const unsigned int i ) const {
    return std::min( i, std::min( suborder + 1, order ) );
}

//! From Eq. 2.8.3-(2)
void QnsCoefficientsFactory::computeCoefficients( const bool extended ) {
    Qns.resizeForOrder( order );
    const double gamma = getParameter( "gamma" );
    for ( unsigned int n = 1; n <= order; n++ ) {
        const int oldSuborder = Qns.suborder( n );
        const unsigned int newSuborder = getSuborder( n );
        Qns.resizeForSuborder( n, newSuborder );
        for ( unsigned int s = oldSuborder + 1; s <= newSuborder; s++ ) {
            if ( n == s ) {
                Qns( n, s ) = ( 2.0 * s - 1.0 ) * Qns( s - 1, s - 1 );
            } else if ( n == s + 1 ) {
                Qns( n, s ) = ( 2.0 * s + 1.0 ) * gamma * Qns( s, s );
            } else {
                Qns( n, s ) = ( ( 2.0 * n - 1.0 ) * gamma * Qns( n - 1, s )
                                - ( n + s - 1.0 ) * Qns( n - 2, s ) ) / ( n - s );
            }
        }
    }
}

//! From Eq. 3.1-(8)
void QnsCoefficientsFactory::computeDependentCoefficients() {
    dQns.resizeForOrder( order );
    for ( int n = 0; n <= (int) order; n++ ) {
        const int oldSuborder = dQns.suborder( n );
        const int newSuborder = getSuborder( n );
        dQns.resizeForSuborder( n, newSuborder );
        for ( int s = oldSuborder + 1; s <= newSuborder; s++ ) {
            if ( n != s && s + 1 <= Qns.suborder( n ) ) {
                dQns( n, s ) = Qns( n, s + 1 );
            } else {
                dQns( n, s ) = 0.0;
            }
        }
    }
}


///////////////////////////////////////////////////////////////////////////////
/////////////////////////// K0nsCoefficientsFactory ///////////////////////////
///////////////////////////////////////////////////////////////////////////////

//! From Eq. 2.7.3-(8)
void K0nsCoefficientsFactory::initializeCoefficients() {
    K0ns = { { 1.0, -1.0 } };
}

void K0nsCoefficientsFactory::initializeDependentCoefficients() {
    dK0ns = DoubleOrderCoefficients();
}

unsigned int K0nsCoefficientsFactory::getSuborder( const unsigned int i ) const {
    return std::min( i + 1, suborder );
}

//! From Eq. 2.7.3-(7)
void K0nsCoefficientsFactory::computeCoefficients( const bool extended ) {
    K0ns.resizeForOrder( order );
    const double Chi = getParameter( "Chi" );
    for ( unsigned int n = 1; n <= order; n++ ) {
        const int oldSuborder = K0ns.suborder( n );
        const unsigned int newSuborder = getSuborder( n );
        K0ns.resizeForSuborder( n, newSuborder );
        // std::cout << oldSuborder << " -> " << newSuborder << std::endl;
        for ( unsigned int s = oldSuborder + 1; s <= newSuborder; s++ ) {
            if ( n == s && n >= 1 ) {
                K0ns( n, s ) = ( 2.0 * s + 1.0 ) / ( s + 1.0 ) * K0ns( s - 1, s );
            } else if ( n == s - 1 && n >= 1 ) {
                K0ns( n, s ) = ( 1.0 - 2.0 * s ) / s * K0ns( s - 2, s - 1 );
            } else if ( n >= s + 1 /*&& s + 1 >= 2*/ ) {
                if ( n == 1 ) {
                    K0ns( n, s ) = 0.5 * ( 3 - 1 / ( Chi * Chi ) );
                } else if ( n >= 2 ) {
                    K0ns( n, s ) = ( 2.0 * n + 1.0 ) / ( n + 1.0 ) * K0ns( n - 1, s ) -
                            double( n + s ) * double ( n - s ) / ( n * ( n + 1.0 ) * Chi * Chi ) * K0ns( n - 2, s );
                }
            }
        }
    }
}

//! From Eq. 3.2-(3)
void K0nsCoefficientsFactory::computeDependentCoefficients() {
    dK0ns.resizeForOrder( order );
    const double Chi = getParameter( "Chi" );
    for ( unsigned int n = 0; n <= order; n++ ) {
        const unsigned int oldSuborder = dK0ns.suborder( n );
        const unsigned int newSuborder = getSuborder( n );
        dK0ns.resizeForSuborder( n, newSuborder );
        for ( unsigned int s = oldSuborder + 1; s <= newSuborder; s++ ) {
            if ( n == s || n == s - 1 ) {
                dK0ns( n, s ) = 0.0;
            } else if ( n > s /*&& n >= 2*/ ) {
                if ( n == 1 ) {
                    dK0ns( n, s ) = std::pow( Chi, -3 );
                } else if ( n >= 2 ) {
                    const double term = ( n + s ) * ( n - s ) / ( n * ( n + 1.0 ) * Chi * Chi );
                    dK0ns( n, s ) = ( 2.0 * n + 1.0 ) / ( n + 1.0 ) * dK0ns( n - 1, s ) - term * dK0ns( n - 2, s )
                            + 2.0 / Chi * term * K0ns( n - 2, s );
                }
            }
        }
    }
}


///////////////////////////////////////////////////////////////////////////////
/////////////////////// NegativeK0nsCoefficientsFactory ///////////////////////
///////////////////////////////////////////////////////////////////////////////

//! From Eq. 2.7.3-(8)
void NegativeK0nsCoefficientsFactory::initializeCoefficients() {
    mK0ns = { { 1.0, -1.0 } };
}

void NegativeK0nsCoefficientsFactory::initializeDependentCoefficients() {
    dmK0ns = DoubleOrderCoefficients();
}

unsigned int NegativeK0nsCoefficientsFactory::getSuborder( const unsigned int i ) const {
    return std::min( i - 1, suborder );
}

//! From Eq. 2.7.3-(6)
void NegativeK0nsCoefficientsFactory::computeCoefficients( const bool extended ) {
    mK0ns.resizeForOrder( order );
    const double Chi = getParameter( "Chi" );
    for ( unsigned int n = 0; n < order; n++ ) {
        const int oldSuborder = mK0ns.suborder( n + 1 );
        const unsigned int newSuborder = getSuborder( n + 1 );
        mK0ns.resizeForSuborder( n + 1, newSuborder );
        for ( unsigned int s = oldSuborder + 1; s <= newSuborder; s++ ) {
            if ( n == s ) {
                mK0ns( n + 1, s ) = 0.0;
            } else if ( n == s + 1 ) {
                mK0ns( n + 1, s ) = std::pow( Chi, 1 + 2 * s ) / std::pow( 2, s );
            } else {
                mK0ns( n + 1, s ) = ( n - 1.0 ) * Chi * Chi / ( ( n + s - 1.0 ) * ( n - s - 1.0 ) ) *
                        ( ( 2.0 * n - 3.0 ) * mK0ns( n, s ) - ( n - 2.0 ) * mK0ns( n - 1, s ) );
            }
        }
    }
}

//! From Eq. 3.1-(7)
void NegativeK0nsCoefficientsFactory::computeDependentCoefficients() {
    dmK0ns.resizeForOrder( order );
    const double Chi = getParameter( "Chi" );
    for ( unsigned int n = 0; n < order; n++ ) {
        const int oldSuborder = dmK0ns.suborder( n + 1 );
        const unsigned int newSuborder = getSuborder( n + 1 );
        dmK0ns.resizeForSuborder( n + 1, newSuborder );
        for ( unsigned int s = oldSuborder + 1; s <= newSuborder; s++ ) {
            if ( n == s ) {
                dmK0ns( n + 1, s ) = 0.0;
            } else if ( n == s + 1 ) {
                dmK0ns( n + 1, s ) = ( 1 + 2 * s ) * std::pow( Chi, 2 * s ) / std::pow( 2, s );
            } else {
                dmK0ns( n + 1, s ) = ( n - 1.0 ) * Chi * Chi / ( ( n + s - 1.0 ) * ( n - s - 1.0 ) ) *
                        ( ( 2.0 * n - 3.0 ) * dmK0ns( n, s ) - ( n - 2.0 ) * dmK0ns( n - 1, s ) ) +
                        2 / Chi * mK0ns( n + 1, s );
            }
        }
    }
}


///////////////////////////////////////////////////////////////////////////////
/////////////////////////// VnsCoefficientsFactory ////////////////////////////
///////////////////////////////////////////////////////////////////////////////

//! From Eq. 2.8.2-(2)
void VnsCoefficientsFactory::initializeCoefficients() {
    Vns = { { 1.0 } };
}

void VnsCoefficientsFactory::initializeDependentCoefficients() {
    Vmns = DoubleOrderCoefficients();
}

bool VnsCoefficientsFactory::changingParameterResetsCoefficients( std::string name ) const {
    return name != "m";
}

unsigned int VnsCoefficientsFactory::getExtendedOrder() const {
    const unsigned int m = getParameter( "m" );
    return std::max( order, m );
}

unsigned int VnsCoefficientsFactory::getSuborder( const unsigned int i ) const {
    return std::min( i, suborder );
}

//! From Eq. 2.8.2-(1), 2.8.2-(2)
void VnsCoefficientsFactory::computeCoefficients( const bool extended ) {
    const int oldOrder = Vns.order();
    const int order = getOrder( extended );
    Vectori oldSuborders;
    // Resize all first, as we are setting elements n+1 and n+2
    Vns.resizeForOrder( order );
    for ( int n = 0; n <= order; n++ ) {
        oldSuborders.push_back( Vns.suborder( n ) );
        Vns.resizeForSuborder( n, getSuborder( n ) );
    }
    // Then compute
    for ( int n = std::max( 0, oldOrder - 2 ); n <= order; n++ ) {
        const int suborder = Vns.suborder( n );
        for ( int s = std::max( 0, oldSuborders( n ) - 2 ); s <= suborder; s++ ) {
            if ( ( n - s ) % 2 != 0 ) { // (n - s) odd
                Vns( n, s ) = 0.0;
            } else { // (n - s) even
                if ( n == s && n + 1 <= order && s + 1 <= Vns.suborder( n + 1 ) ) {
                    Vns( s + 1, s + 1 ) = Vns( s, s ) / ( 2.0 * s + 2.0 );
                }
                if ( n + 2 <= order ) {
                    Vns( n + 2, s ) = Vns( n, s ) * ( s - n - 1.0 ) / ( n + s + 2.0 );
                }
            }
        }
    }
}

//! From first paragraph of Section 2.8.2
void VnsCoefficientsFactory::computeDependentCoefficients() {
    const unsigned int order = getOrder( true );
    const unsigned int m = getParameter( "m" );
    Vmns.resizeForOrder( order );
    for ( unsigned int n = 0; n <= order; n++ ) {
        const int oldSuborder = Vmns.suborder( n );
        const unsigned int newSuborder = getSuborder( n );
        Vmns.resizeForSuborder( n, newSuborder );
        for ( unsigned int s = oldSuborder + 1; s <= newSuborder; s++ ) {
            // std::cout << (n+s) << "! = " << factorial( n + s ) << std::endl;
            Vmns( n, s ) = factorial( n + s ) / factorial( n - m ) * Vns( n, s );
        }
    }
}


///////////////////////////////////////////////////////////////////////////////
////////////////////////// JacobiPolynomialsFactory ///////////////////////////
///////////////////////////////////////////////////////////////////////////////

//! From Eq. 2.7.4-(3)
void JacobiPolynomialsFactory::initializeCoefficients() {
    P = { 1.0 };
}

//! From Eq. 2.7.4-(2)
void JacobiPolynomialsFactory::computeCoefficients( const bool extended ) {
    if ( order > maximumComputedOrder ) {
        P.resizeForOrder( order );
        // Get parameters
        const double gamma = getParameter( "gamma" );
        const double v     = getParameter( "v" );
        const double w     = getParameter( "w" );
        for ( unsigned int l = maximumComputedOrder + 1; l <= order; l++ ) {
            // Compute repeated terms
            const double vwl   = v + w + l;
            const double vw2l  = vwl + l;
            const double denom = 2 * l * vwl * ( vw2l - 2 );
            const double c1    = ( vw2l - 1 ) * ( vw2l * ( vw2l - 2 ) * gamma + v * v - w * w );
            const double c2    = -2 * ( l + v - 1 ) * ( l + w - 1 ) * vw2l;
            // Compute polynomial value from previous two steps
            const double P_1 = P( l - 1 );
            const double P_2 = l > 1 ? P( l - 2 ) : 0.0;
            P( l ) = c1 / denom * P_1 + c2 / denom * P_2;
        }
    }
}


///////////////////////////////////////////////////////////////////////////////
////////////////////////// CsSsCoefficientsFactory ////////////////////////////
///////////////////////////////////////////////////////////////////////////////

//! From Eq. 2.5.3-(6)
void CsSsCoefficientsFactory::initializeCoefficients() {
    Cs = { 1.0 };
    Ss = { 0.0 };
}

//! From Eq. 2.5.3-(6)
void CsSsCoefficientsFactory::computeCoefficients( const bool extended ) {
    if ( order > maximumComputedOrder ) {
        Cs.resizeForOrder( order );
        Ss.resizeForOrder( order );
        // Get parameters
        const double a = getParameter( "a" );
        const double b = getParameter( "b" );
        for ( unsigned int s = maximumComputedOrder + 1; s <= order; s++ ) {
            // Compute coefficients recursively
            Cs( s ) = a * Cs( s - 1 ) - b * Ss( s - 1 );
            Ss( s ) = b * Cs( s - 1 ) + a * Ss( s - 1 );
        }
    }
}


} // namespace coefficients_factories

} // namespace sst

} // namespace propagators

} // namespace tudat
