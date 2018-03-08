/*    Copyright (c) 2010-2017, Delft University of Technology
 *    All rigths reserved
 *
 *    This file is part of the Tudat. Redistribution and use in source and
 *    binary forms, with or without modification, are permitted exclusively
 *    under the terms of the Modified BSD license. You should have received
 *    a copy of the license with this file. If not, please or visit:
 *    http://tudat.tudelft.nl/LICENSE.
 */

#define BOOST_TEST_MAIN


#include <limits>
#include <boost/test/unit_test.hpp>
#include <eigen/Eigen/Eigen>

#include "Tudat/Astrodynamics/Propagators/DSST/utilities/coefficientsfactory.h"


namespace tudat
{

namespace unit_tests
{

BOOST_AUTO_TEST_SUITE( test_dsst_propagator )

BOOST_AUTO_TEST_CASE( coefficients_factories )
{
    using namespace tudat::propagators::dsst;
    using namespace coefficients_factories;

    auto tolerance = 100 * std::numeric_limits< double >::epsilon( );

    K0nsCoefficientsFactory k0nsFactory( 2.1, 3, 3 );
    DoubleVectord K0ns  = k0nsFactory.getCoefficients();
    DoubleVectord dK0ns = k0nsFactory.getDerivatives();

    // Q_{n,s}
    QnsCoefficientsFactory qnsFactory( 1.5, 4, 4 );
    DoubleVectord Qns = qnsFactory.getCoefficients();
    BOOST_CHECK_CLOSE_FRACTION( Qns( 2, 1 ), 4.5, tolerance );

    // d(Q_{n,s}/d(\gamma)
    DoubleVectord dQns = qnsFactory.getDerivatives();
    BOOST_CHECK_CLOSE_FRACTION( dQns( 3, 2 ), Qns( 3, 3 ), tolerance );

    // G_s, H_s polynomials
    GsHsCoefficientsFactory ghFactory( 1.0, 1.0, 1.5, 0.5, 2 );
    BOOST_CHECK_CLOSE_FRACTION( ghFactory.getGsCoefficients()( 2 ), 3.0, tolerance );
    BOOST_CHECK_CLOSE_FRACTION( ghFactory.getHsCoefficients()( 2 ), 4.0, tolerance );

    // G_s derivatives
    BOOST_CHECK_CLOSE_FRACTION( ghFactory.getGsHDerivatives    ()( 2 ), -1.0, tolerance );
    BOOST_CHECK_CLOSE_FRACTION( ghFactory.getGsKDerivatives    ()( 2 ),  7.0, tolerance );
    BOOST_CHECK_CLOSE_FRACTION( ghFactory.getGsAlphaDerivatives()( 2 ),  2.0, tolerance );
    BOOST_CHECK_CLOSE_FRACTION( ghFactory.getGsBetaDerivatives ()( 2 ),  6.0, tolerance );

    // V_{n,s} coefficients
    VnsCoefficientsFactory vnsFactory( 5, 5 );
    DoubleVectord Vns = vnsFactory.getCoefficients();
    BOOST_CHECK_CLOSE_FRACTION( Vns( 0, 0 ),  1,     tolerance );
    BOOST_CHECK_CLOSE_FRACTION( Vns( 3, 1 ), -0.125, tolerance );

    // V_{n,s}^m coefficient
    DoubleVectord Vmns = vnsFactory.getCoefficientsWithSuperscript( 2 );
    BOOST_CHECK_CLOSE_FRACTION( Vmns( 3, 3 ), 15, tolerance );
}

BOOST_AUTO_TEST_SUITE_END( )

} // namespace unit_tests

} // namespace tudat
