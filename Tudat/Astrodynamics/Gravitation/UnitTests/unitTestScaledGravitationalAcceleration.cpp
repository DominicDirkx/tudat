#define BOOST_TEST_MAIN

#include <iostream>

#include <boost/make_shared.hpp>
#include <boost/test/unit_test.hpp>

#include "Tudat/Basics/testMacros.h"
#include <Tudat/Astrodynamics/Gravitation/sphericalHarmonicsGravityModelBase.h>
#include <Tudat/Astrodynamics/Gravitation/centralGravityModel.h>

#include "Tudat/Astrodynamics/Gravitation/scaledGravitationalAccelerationModel.h"


namespace tudat
{

namespace unit_tests
{
using namespace tudat::gravitation;

BOOST_AUTO_TEST_SUITE( test_scaled_gravitational_acceleration )

BOOST_AUTO_TEST_CASE( testScaledGravitationalAcceleration )
{
    double oldMu = 25.0;
    std::shared_ptr< CentralGravitationalAccelerationModel3d > baseAcceleration =
            std::make_shared< CentralGravitationalAccelerationModel3d >(
                boost::lambda::constant( ( Eigen::Vector3d( )<< 4.6, 3.5, -0.54 ).finished( ) ),
                boost::lambda::constant( oldMu ),
                boost::lambda::constant( ( Eigen::Vector3d( )<< -7.4, 1.4, -0.734 ).finished( ) ), 0 );

    double newMu = 20.0;
    std::shared_ptr< ScaledGravitationalAccelerationModel< CentralGravitationalAccelerationModel3d > > derivedAcceleration =
            std::make_shared< ScaledGravitationalAccelerationModel< CentralGravitationalAccelerationModel3d > >(
                baseAcceleration , boost::lambda::constant( newMu ), 1, 1 );

    derivedAcceleration->updateMembers( );


    TUDAT_CHECK_MATRIX_CLOSE_FRACTION( ( -1.0 * ( newMu / oldMu ) * baseAcceleration->getAcceleration( ) ),
                                       derivedAcceleration->getAcceleration( ), std::numeric_limits< double >::epsilon( ) );

}

BOOST_AUTO_TEST_SUITE_END( )

}

}
