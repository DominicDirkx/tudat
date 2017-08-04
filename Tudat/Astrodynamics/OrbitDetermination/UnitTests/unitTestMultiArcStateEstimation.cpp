#define BOOST_TEST_MAIN

#include <string>
#include <thread>

#include <limits>

#include <boost/make_shared.hpp>
#include <boost/test/unit_test.hpp>

#include "Tudat/Basics/testMacros.h"
#include "Tudat/SimulationSetup/tudatSimulationHeader.h"
#include "Tudat/Astrodynamics/ObservationModels/linkTypeDefs.h"
#include "Tudat/Astrodynamics/ObservationModels/simulateObservations.h"
#include "Tudat/Astrodynamics/OrbitDetermination/orbitDeterminationManager.h"
#include "Tudat/Astrodynamics/OrbitDetermination/UnitTests/orbitDeterminationTestCases.h"
#include "Tudat/SimulationSetup/EnvironmentSetup/createGroundStations.h"


namespace tudat
{
namespace unit_tests
{
BOOST_AUTO_TEST_SUITE( test_multi_arc_parameter_estimation )


BOOST_AUTO_TEST_CASE( test_MultiArcStateEstimation )
{
    std::pair< boost::shared_ptr< PodOutput< double > >, boost::shared_ptr< PodInput< double, double > > > podData;
    std::vector< double > integrationArcLimits;

    // Execute test for linked arcs and separate arcs.
    for( unsigned int testCase = 0; testCase < 2; testCase++ )
    {
        Eigen::VectorXd parameterError = executeMultiArcParameterEstimation< long double, tudat::Time, long double >(
                    podData, testCase, integrationArcLimits );
        int numberOfEstimatedArcs = ( parameterError.rows( ) - 3 ) / 6;

        std::cout<<parameterError.transpose( )<<std::endl;
        for( int i = 0; i < numberOfEstimatedArcs; i++ )
        {
            for( unsigned int j = 0; j < 3; j++ )
            {
                BOOST_CHECK_SMALL( std::fabs( parameterError( i * 6 + j ) ), 1E-4 );
                BOOST_CHECK_SMALL( std::fabs( parameterError( i * 6 + j + 3 ) ), 1.0E-10  );
            }
        }

        BOOST_CHECK_SMALL( std::fabs( parameterError( parameterError.rows( ) - 3 ) ), 1.0E-20 );
        BOOST_CHECK_SMALL( std::fabs( parameterError( parameterError.rows( ) - 2 ) ), 1.0E-12 );
        BOOST_CHECK_SMALL( std::fabs( parameterError( parameterError.rows( ) - 1 ) ), 1.0E-12 );
    }

}

BOOST_AUTO_TEST_SUITE_END( )

}

}
