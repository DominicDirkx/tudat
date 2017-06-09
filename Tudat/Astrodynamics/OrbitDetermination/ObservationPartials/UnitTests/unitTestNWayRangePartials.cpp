/*    Copyright (c) 2010-2012 Delft University of Technology.
 *
 *    This software is protected by national and international copyright.
 *    Any unauthorized use, reproduction or modification is unlawful and
 *    will be prosecuted. Commercial and non-private application of the
 *    software in any form is strictly prohibited unless otherwise granted
 *    by the authors.
 *
 *    The code is provided without any warranty; without even the implied
 *    warranty of merchantibility or fitness for a particular purpose.
 *
 *    Changelog
 *      YYMMDD    Author            Comment
 *
 *    References
 *
 */

#define BOOST_TEST_MAIN

#include <limits>
#include <string>
#include <vector>

#include <boost/test/unit_test.hpp>
#include <boost/make_shared.hpp>
#include <boost/lambda/lambda.hpp>

#include "Tudat/Basics/testMacros.h"

#include "Tudat/InputOutput/basicInputOutput.h"
#include "Tudat/External/SpiceInterface/spiceInterface.h"

#include "Tudat/SimulationSetup/EstimationSetup/createObservationModel.h"
#include "Tudat/Astrodynamics/ObservationModels/oneWayRangeObservationModel.h"
#include "Tudat/Astrodynamics/OrbitDetermination/EstimatableParameters/constantRotationRate.h"
#include "Tudat/SimulationSetup/EstimationSetup/createObservationPartials.h"
#include "Tudat/Astrodynamics/OrbitDetermination/ObservationPartials/UnitTests/numericalObservationPartial.h"
#include "Tudat/SimulationSetup/EnvironmentSetup/createGroundStations.h"
#include "Tudat/SimulationSetup/EnvironmentSetup/defaultBodies.h"

#include "Tudat/Astrodynamics/OrbitDetermination/ObservationPartials/UnitTests/observationPartialTestFunctions.h"

namespace tudat
{
namespace unit_tests
{

using namespace tudat::gravitation;
using namespace tudat::ephemerides;
using namespace tudat::observation_models;
using namespace tudat::simulation_setup;
using namespace tudat::spice_interface;
using namespace tudat::observation_partials;
using namespace tudat::estimatable_parameters;

BOOST_AUTO_TEST_SUITE( test_N_way_observation_partials)


std::vector< double > getRetransmissionDelays( const double evaluationTime, const int numberOfRetransmitters )
{
    std::vector< double > retransmissionDelays;

        for( int i = 0; i < numberOfRetransmitters; i++ )
        {
            retransmissionDelays.push_back( evaluationTime * 5.0E-17 * static_cast< double >( i + 1 ) );
        }
    return retransmissionDelays;
}

//! Test partial derivatives of one-way range observable, using general test suite of observation partials.
BOOST_AUTO_TEST_CASE( testnWayRangePartials )
{
    // Define and create ground stations.
    std::vector< std::pair< std::string, std::string > > groundStations;
    groundStations.resize( 2 );
    groundStations[ 0 ] = std::make_pair( "Earth", "Graz" );
    groundStations[ 1 ] = std::make_pair( "Mars", "MSL" );



    for( unsigned int linkNumber = 0; linkNumber < 3; linkNumber++ )
    {
        // Set link ends for observation model
        LinkEnds linkEnds;
        linkEnds[ transmitter ] = groundStations[ 0 ];
        linkEnds[ reflector1 ] = groundStations[ 1 ];

        if( linkNumber > 0 )
        {
            linkEnds[ reflector2 ] = groundStations[ 0 ];
        }
        if( linkNumber > 1 )
        {
            linkEnds[ reflector3 ] = groundStations[ 1 ];
        }

        if( linkNumber % 2 == 0 )
        {
            linkEnds[ receiver ] = groundStations[ 0 ];
        }
        else
        {
            linkEnds[ receiver ] = groundStations[ 1 ];
        }
        // Test partials with constant ephemerides (allows test of position partials)
        {
            // Create environment
            NamedBodyMap bodyMap = setupEnvironment( groundStations, 1.0E7, 1.2E7, 1.1E7, true );



            // Generate one-way range model
            std::vector< std::string > perturbingBodies;
            perturbingBodies.push_back( "Earth" );
            std::vector< boost::shared_ptr< observation_models::ObservationSettings > > legObservationModels;

            for( unsigned int i = 0; i < linkNumber + 2; i ++ )
            {
                legObservationModels.push_back(
                            boost::make_shared< observation_models::ObservationSettings >(
                                one_way_range, boost::make_shared< FirstOrderRelativisticLightTimeCorrectionSettings >(
                                    perturbingBodies ) ) );
            }

            boost::shared_ptr< ObservationModel< 1 > > nWayRangeModel =
                    observation_models::ObservationModelCreator< 1, double, double >::createObservationModel(
                        linkEnds, boost::make_shared< observation_models::NWayRangeObservationSettings >(
                            legObservationModels, boost::bind( &getRetransmissionDelays, _1, linkNumber + 1 ) ), bodyMap  );

            // Create parameter objects.
            boost::shared_ptr< EstimatableParameterSet< double > > fullEstimatableParameterSet =
                    createEstimatableParameters( bodyMap, 1.1E7 );

            testObservationPartials< 1 >(
                        nWayRangeModel, bodyMap, fullEstimatableParameterSet, linkEnds, n_way_range, 1.0E-6, true, true, 1.0,
                        ( Eigen::Vector4d( )<<10.0, 1.0, 1.0, 10.0 ).finished( ) );
        }

        // Test partials with real ephemerides (without test of position partials)
        {
            // Create environment
            NamedBodyMap bodyMap = setupEnvironment( groundStations, 1.0E7, 1.2E7, 1.1E7, false );

            // Generate one-way range model
            std::vector< std::string > perturbingBodies;
            perturbingBodies.push_back( "Earth" );
            std::vector< boost::shared_ptr< observation_models::ObservationSettings > > legObservationModels;
            for( unsigned int i = 0; i < linkNumber + 2; i ++ )
            {
                legObservationModels.push_back(
                            boost::make_shared< observation_models::ObservationSettings >(
                                one_way_range, boost::make_shared< FirstOrderRelativisticLightTimeCorrectionSettings >(
                                    perturbingBodies ) ) );
            }

            boost::shared_ptr< ObservationModel< 1 > > nWayRangeModel =
                    observation_models::ObservationModelCreator< 1, double, double >::createObservationModel(
                        linkEnds, boost::make_shared< observation_models::NWayRangeObservationSettings >(
                            legObservationModels, boost::bind( &getRetransmissionDelays, _1, linkNumber + 1 ) ), bodyMap  );

            // Create parameter objects.
            boost::shared_ptr< EstimatableParameterSet< double > > fullEstimatableParameterSet =
                    createEstimatableParameters( bodyMap, 1.1E7 );

            testObservationPartials< 1 >(
                        nWayRangeModel, bodyMap, fullEstimatableParameterSet, linkEnds, n_way_range, 1.0E-6, false, true, 1.0,
                        ( Eigen::Vector4d( )<<10.0, 1.0, 1.0, 20.0 ).finished( ) );
        }
    }
}


BOOST_AUTO_TEST_SUITE_END( )

} // namespace unit_tests

} // namespace tudat




