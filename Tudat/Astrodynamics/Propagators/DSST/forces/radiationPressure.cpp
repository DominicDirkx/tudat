/*    Copyright (c) 2010-2017, Delft University of Technology
 *    All rigths reserved
 *
 *    This file is part of the Tudat. Redistribution and use in source and
 *    binary forms, with or without modification, are permitted exclusively
 *    under the terms of the Modified BSD license. You should have received
 *    a copy of the license with this file. If not, please or visit:
 *    http://tudat.tudelft.nl/LICENSE.
 */

#include "radiationPressure.h"

namespace tudat
{

namespace propagators
{

namespace sst
{

namespace force_models
{


//! Update instance's members that are computed from the current auxiliary elements.
void RadiationPressure::updateMembers( )
{
    // Update vector to source
    r3 = radiationPressureAM->getCurrentVectorToSource();

    // Update direction cosines
    ThirdBodyPerturbed::updateMembers();

    // Update integration limits
    NonConservative::updateMembers();
}


//! Update the values of the minimum and maximum true longitude for the averaging integral.
void RadiationPressure::determineIntegrationLimits( ) {
    using namespace mathematical_constants;
    /*
    const double Ro = *radiationPressureAM->getRadiationPressureInterface()->getOccultingBodyRadii().begin();

    // Danielson, Section 3.5
    const double M_ = std::pow( Ro / ( a * B ), 2 );
    const double B_ = alpha * beta + M_ * h * k;
    const double C_ = alpha * alpha - beta * beta + M_ * ( k * k - h * h );
    const double D_ = 1 - beta * beta - M_ * ( 1 + h * h );

    // Using different notation.
    const double a = 4 * B_ * B_ + C_ * C_;
    const double b = 8 * B_ * M_ * h + 4 * C_ * M_ * k;
    const double c = -4 * B_ * B_ + 4 * M_ * M_ * h * h - 2 * D_ * C_ + 4 * M_ * M_ * k * k;
    const double d = -8 + B_ * M_ * h - 4 * D_ * M_ * k;
    const double e = -4 * M_ * M_ * h * h + D_ * D_;

    // Solving quartic equation [https://en.wikipedia.org/wiki/Quartic_function#General_formula_for_roots]
    const double p = ( 8 * a * c - 3 * b * b ) / ( 8 * a * a );
    const double q = ( b * b * b - 4 * a * b * c + 8 * a * a * d ) / ( 8 * a * a * a );
    const double d0 = c * c - 3 * b * d + 12 * a * e;
    const double d1 = 2 * c * c * c - 9 * b * c * d + 27 * b * b * e + 27 * a * d * d - 72 * a * c * e;
    const double Q = std::pow( ( d1 + std::sqrt( d1 * d1 - 4 * d0 * d0 * d0 ) ) / 2, 1/3.0 );
    const double S = 0.5 * std::sqrt( -2/3.0 * p + ( Q + d0 / Q ) / ( 3 * a ) );
    const double r = -4 * S * S - 2 * p;
    double radicand = r + q / S;
    double R;
    if ( radicand < 0 ) {
        radicand -= 2 * q / S;
        R = -b / ( 4 * a ) + S;
    } else {
        R = -b / ( 4 * a ) - S;
    }
    const double s = 0.5 * std::sqrt( radicand );
    const double cosx1 = R + s;
    const double cosx2 = R - s;
    const double x1 = std::acos( cosx1 );
    const double x2 = std::acos( cosx2 );

    if ( x1 == x1 && x2 == x2 ) {  // both x1 and x2 are numbers (i.e. not NaN)
        L1 = std::min( x1, x2 );
        L2 = std::max( x1, x2 );
        std::cout << "L1 = " << L1 << ", L2 = " << L2 << std::endl;
        return;
    }
    */
    // Could not find roots, thus consider the whole period
    L1 = -PI;
    L2 =  PI;
}


} // namespace force_models

} // namespace sst

} // namespace propagators

} // namespace tudat
