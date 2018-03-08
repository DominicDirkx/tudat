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

#include "Tudat/Mathematics/BasicMathematics/basicMathematicsFunctions.h"
#include "Tudat/Astrodynamics/BasicAstrodynamics/modifiedEquinoctialElementConversions.h"


namespace tudat
{

namespace unit_tests
{

BOOST_AUTO_TEST_SUITE( test_dsst_propagator )

BOOST_AUTO_TEST_CASE( equinoctial_elements )
{
    using namespace tudat::orbital_element_conversions;
    using namespace tudat::basic_mathematics;
    using namespace tudat::mathematical_constants;

    using namespace tudat::propagators::dsst;
    using namespace tudat::propagators::dsst::coefficients_factories;

    auto tolerance = 100 * std::numeric_limits< double >::epsilon( );

    // Convert TO equinoctial elements

    const double mu = 4e14;
    Eigen::Vector6d keplerian;
    keplerian << 7e6, 0.1, 1, 2, 3, 4;
    Eigen::Vector6d cartesian = convertKeplerianToCartesianElements( keplerian, mu );

    Eigen::Vector6d equinoctialFromKep = convertKeplerianToEquinoctialElements( keplerian, mu );
    Eigen::Vector6d equinoctialFromCar = convertCartesianToEquinoctialElements( cartesian, mu );

    for ( unsigned int i = 0; i < 6; i++ ) {
        BOOST_CHECK_CLOSE_FRACTION( equinoctialFromKep[ i ], equinoctialFromCar[ i ], tolerance );
    }
    // std::cout << equinoctialFromKep << std::endl;
    // std::cout << equinoctialFromCar << std::endl;


    // Convert FROM equinoctial elements

    Eigen::Vector6d cartesianFromEqui = convertEquinoctialToCartesianElements( equinoctialFromCar, mu );
    for ( unsigned int i = 0; i < 6; i++ ) {
        BOOST_CHECK_CLOSE_FRACTION( cartesian[ i ], cartesianFromEqui[ i ], tolerance );
    }
//    std::cout << cartesian << std::endl;
//    std::cout << cartesianFromEqui << std::endl;

    Eigen::Vector6d keplerianFromEqui = convertEquinoctialToKeplerianElements( equinoctialFromKep, mu );
    for ( unsigned int i = 0; i < 6; i++ ) {
        BOOST_CHECK_CLOSE_FRACTION( keplerian[ i ], keplerianFromEqui[ i ], tolerance );
    }
//    std::cout << keplerian << std::endl;
//    std::cout << keplerianFromEqui << std::endl;


    // Modified equinoctial elements

    Eigen::Vector6d equinoctial = equinoctialFromKep;
    Eigen::Vector6d modEquinoctialFromEqui = convertEquinoctialToModifiedEquinoctialElements( equinoctial );
    Eigen::Vector6d modEquinoctialFromKep = convertKeplerianToModifiedEquinoctialElements( keplerian, false );
    for ( unsigned int i = 0; i < 6; i++ ) {
        BOOST_CHECK_CLOSE_FRACTION( modEquinoctialFromEqui[ i ], modEquinoctialFromKep[ i ], tolerance );
    }

    Eigen::Vector6d equinoctialFromModEqui = convertModifiedEquinoctialToEquinoctialElements( equinoctial );
    for ( unsigned int i = 0; i < 6; i++ ) {
        BOOST_CHECK_CLOSE_FRACTION( equinoctial[ i ], equinoctialFromModEqui[ i ], tolerance );
    }


    // Get eccentric and true longitudes

    const double ectrLongitude = getEccentricLongitude( equinoctial );
    const double trueLongitude = getTrueLongitude( equinoctial );
    const double meanLongitude = equinoctial[ meanLongitudeIndex ];

    const double trueAnomaly = keplerian [ trueAnomalyIndex ];
    const double e = keplerian [ eccentricityIndex ];
    const double ectrAnomaly = convertTrueAnomalyToEllipticalEccentricAnomaly( trueAnomaly, e );
    const double meanAnomaly = convertEllipticalEccentricAnomalyToMeanAnomaly( ectrAnomaly, e );

    BOOST_CHECK_CLOSE_FRACTION( computeModulo( trueLongitude - trueAnomaly, 2 * PI ),
                                computeModulo( ectrLongitude - ectrAnomaly, 2 * PI ), tolerance );

    BOOST_CHECK_CLOSE_FRACTION( computeModulo( trueLongitude - trueAnomaly, 2 * PI ),
                                computeModulo( meanLongitude - meanAnomaly, 2 * PI ), tolerance );

//    std::cout << std::fmod( trueLongitude - trueAnomaly + 2 * M_PI, 2 * M_PI ) << std::endl;
//    std::cout << std::fmod( ectrLongitude - ectrAnomaly + 2 * M_PI, 2 * M_PI ) << std::endl;
//    std::cout << std::fmod( meanLongitude - meanAnomaly + 2 * M_PI, 2 * M_PI ) << std::endl;

}

BOOST_AUTO_TEST_SUITE_END( )

} // namespace unit_tests

} // namespace tudat
