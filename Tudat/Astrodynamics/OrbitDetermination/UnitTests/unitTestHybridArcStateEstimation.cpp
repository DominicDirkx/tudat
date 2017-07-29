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
#include "Tudat/SimulationSetup/EnvironmentSetup/createGroundStations.h"


namespace tudat
{
namespace unit_tests
{
BOOST_AUTO_TEST_SUITE( test_hybrid_arc_estimation )

//Using declarations.
using namespace tudat::observation_models;
using namespace tudat::orbit_determination;
using namespace tudat::estimatable_parameters;
using namespace tudat::interpolators;
using namespace tudat::numerical_integrators;
using namespace tudat::spice_interface;
using namespace tudat::simulation_setup;
using namespace tudat::orbital_element_conversions;
using namespace tudat::ephemerides;
using namespace tudat::propagators;
using namespace tudat::basic_astrodynamics;
using namespace tudat::unit_conversions;

template< typename ObservationScalarType = double , typename TimeType = double , typename StateScalarType  = double >
Eigen::VectorXd  executeParameterEstimation(
        const int linkArcs )
{
    //Load spice kernels.f
    std::string kernelsPath = input_output::getSpiceKernelPath( );
    spice_interface::loadSpiceKernelInTudat( kernelsPath + "de-403-masses.tpc");
    spice_interface::loadSpiceKernelInTudat( kernelsPath + "naif0009.tls");
    spice_interface::loadSpiceKernelInTudat( kernelsPath + "pck00009.tpc");
    spice_interface::loadSpiceKernelInTudat( kernelsPath + "de421.bsp");

    //Define setting for total number of bodies and those which need to be integrated numerically.
    //The first numberOfNumericalBodies from the bodyNames vector will be integrated numerically.

    std::vector< std::string > bodyNames;
    bodyNames.push_back( "Sun" );
    bodyNames.push_back( "Mars" );
    bodyNames.push_back( "Jupiter" );
    bodyNames.push_back( "Earth" );

    // Specify initial time
    double initialEphemerisTime = 1.0E7;
    double finalEphemerisTime = initialEphemerisTime + 10.0 * 86400.0;
    double maximumTimeStep = 3600.0;
    double buffer = 5.0 * maximumTimeStep;

    // Create bodies needed in simulation
    std::map< std::string, boost::shared_ptr< BodySettings > > bodySettings =
            getDefaultBodySettings( bodyNames, initialEphemerisTime - buffer, finalEphemerisTime + buffer );

    NamedBodyMap bodyMap = createBodies( bodySettings );

    bodyMap[ "Orbiter" ] = boost::make_shared< Body >( );
    bodyMap[ "Orbiter" ]->setConstantBodyMass( 5.0E3 );
    bodyMap[ "Orbiter" ]->setEphemeris( boost::make_shared< MultiArcEphemeris >(
                                            std::map< double, boost::shared_ptr< Ephemeris > >( ),
                                            "Mars", "ECLIPJ2000" ) );

    double referenceAreaRadiation = 4.0;
    double radiationPressureCoefficient = 1.2;
    std::vector< std::string > occultingBodies;
    occultingBodies.push_back( "Earth" );
    boost::shared_ptr< RadiationPressureInterfaceSettings > orbiterRadiationPressureSettings =
            boost::make_shared< CannonBallRadiationPressureInterfaceSettings >(
                "Sun", referenceAreaRadiation, radiationPressureCoefficient, occultingBodies );

    // Create and set radiation pressure settings
    bodyMap[ "Orbiter" ]->setRadiationPressureInterface(
                "Sun", createRadiationPressureInterface(
                    orbiterRadiationPressureSettings, "Orbiter", bodyMap ) );


    // Finalize body creation.
    setGlobalFrameBodyEphemerides( bodyMap, "SSB", "ECLIPJ2000" );

    // Create ground stations
    std::pair< std::string, std::string > grazStation = std::pair< std::string, std::string >( "Earth", "" );
    std::pair< std::string, std::string > mslStation = std::pair< std::string, std::string >( "Mars", "MarsStation" );

    createGroundStation( bodyMap.at( "Mars" ), "MarsStation", ( Eigen::Vector3d( )<< 100.0, 0.5, 2.1 ).finished( ),
                         coordinate_conversions::geodetic_position );

    std::vector< std::pair< std::string, std::string > > groundStations;
    groundStations.push_back( grazStation );
    groundStations.push_back( mslStation );




    // Set accelerations between bodies that are to be taken into account.
    SelectedAccelerationMap singleArcAccelerationMap;
    std::map< std::string, std::vector< boost::shared_ptr< AccelerationSettings > > > accelerationsOfMars;
    accelerationsOfMars[ "Earth" ].push_back( boost::make_shared< AccelerationSettings >( central_gravity ) );
    accelerationsOfMars[ "Sun" ].push_back( boost::make_shared< AccelerationSettings >( central_gravity ) );
    accelerationsOfMars[ "Jupiter" ].push_back( boost::make_shared< AccelerationSettings >( central_gravity ) );
    singleArcAccelerationMap[ "Mars" ] = accelerationsOfMars;

    std::vector< std::string > singleArcBodiesToIntegrate, singleArcCentralBodies;
    singleArcBodiesToIntegrate.push_back( "Mars" );
    singleArcCentralBodies.push_back( "SSB" );

    AccelerationMap singleArcAccelerationModelMap = createAccelerationModelsMap(
                bodyMap, singleArcAccelerationMap, singleArcBodiesToIntegrate, singleArcCentralBodies );
    Eigen::VectorXd singleArcInitialStates = getInitialStatesOfBodies(
                singleArcBodiesToIntegrate, singleArcCentralBodies, bodyMap, initialEphemerisTime );


    boost::shared_ptr< TranslationalStatePropagatorSettings< > > singleArcPropagatorSettings =
            boost::make_shared< TranslationalStatePropagatorSettings< > >(
                singleArcCentralBodies, singleArcAccelerationModelMap, singleArcBodiesToIntegrate,
                singleArcInitialStates, finalEphemerisTime );


    SelectedAccelerationMap multiArcAccelerationMap;
    std::map< std::string, std::vector< boost::shared_ptr< AccelerationSettings > > > accelerationsOfOrbiter;
    accelerationsOfOrbiter[ "Mars" ].push_back( boost::make_shared< SphericalHarmonicAccelerationSettings >( 2, 2 ) );
    accelerationsOfOrbiter[ "Sun" ].push_back( boost::make_shared< AccelerationSettings >( central_gravity ) );
    accelerationsOfOrbiter[ "Sun" ].push_back( boost::make_shared< AccelerationSettings >( cannon_ball_radiation_pressure ) );
    accelerationsOfOrbiter[ "Jupiter" ].push_back( boost::make_shared< AccelerationSettings >( central_gravity ) );
    multiArcAccelerationMap[ "Orbiter" ] = accelerationsOfOrbiter;

    std::vector< std::string > multiArcBodiesToIntegrate, multiArcCentralBodies;
    multiArcBodiesToIntegrate.push_back( "Orbiter" );
    multiArcCentralBodies.push_back( "Mars" );

    AccelerationMap multiArcAccelerationModelMap = createAccelerationModelsMap(
                bodyMap, multiArcAccelerationMap, multiArcBodiesToIntegrate, multiArcCentralBodies );

    // Creater arc times

    double arcDuration = 2.1 * 86400.0;
    double arcOverlap  = 3600.0;

    std::vector< double > integrationArcStarts, integrationArcEnds;
    std::vector< double > integrationArcLimits;

    double integrationStartTime = initialEphemerisTime;
    double integrationEndTime = finalEphemerisTime - 1.0E4;
    double currentStartTime = integrationStartTime;
    double currentEndTime = integrationStartTime + arcDuration;
    do
    {
        integrationArcStarts.push_back( currentStartTime );
        integrationArcLimits.push_back( currentStartTime );
        integrationArcEnds.push_back( currentEndTime );

        currentStartTime = currentEndTime - arcOverlap;
        currentEndTime = currentStartTime + arcDuration;
    }
    while( currentEndTime < integrationEndTime );

    integrationArcLimits.push_back( currentStartTime + arcOverlap );


    // Create list of multi-arc initial states
    unsigned int numberOfIntegrationArcs = integrationArcStarts.size( );
    std::vector< Eigen::VectorXd > multiArcSystemInitialStates;
    multiArcSystemInitialStates.resize( numberOfIntegrationArcs );

    // Define (quasi-arbitrary) arc initial states
    double marsGravitationalParameter =  bodyMap.at( "Mars" )->getGravityFieldModel( )->getGravitationalParameter( );

    for( unsigned int j = 0; j < numberOfIntegrationArcs; j++ )
    {
        Eigen::Vector6d orbiterInitialStateInKeplerianElements;
        orbiterInitialStateInKeplerianElements( semiMajorAxisIndex ) = 6000.0E3;
        orbiterInitialStateInKeplerianElements( eccentricityIndex ) = 0.05;
        orbiterInitialStateInKeplerianElements( inclinationIndex ) = convertDegreesToRadians( 85.3 );
        orbiterInitialStateInKeplerianElements( argumentOfPeriapsisIndex )
                = convertDegreesToRadians( 235.7 - j );
        orbiterInitialStateInKeplerianElements( longitudeOfAscendingNodeIndex )
                = convertDegreesToRadians( 23.4 + j );
        orbiterInitialStateInKeplerianElements( trueAnomalyIndex ) = convertDegreesToRadians( 139.87 + j * 10.0 );

        // Convert state from Keplerian elements to Cartesian elements.
        multiArcSystemInitialStates[ j ]  = convertKeplerianToCartesianElements(
                    orbiterInitialStateInKeplerianElements,
                    marsGravitationalParameter );
    }


    // Create propagation settings for each arc
    std::vector< boost::shared_ptr< SingleArcPropagatorSettings< double > > > arcPropagationSettingsList;
    for( unsigned int i = 0; i < numberOfIntegrationArcs; i++ )
    {
        arcPropagationSettingsList.push_back(
                    boost::make_shared< TranslationalStatePropagatorSettings< double > >
                    ( multiArcCentralBodies, multiArcAccelerationModelMap, multiArcBodiesToIntegrate,
                      multiArcSystemInitialStates.at( i ), integrationArcEnds.at( i ) ) );
    }

    boost::shared_ptr< MultiArcPropagatorSettings< > > multiArcPropagatorSettings =
            boost::make_shared< MultiArcPropagatorSettings< > >( arcPropagationSettingsList, linkArcs );

    boost::shared_ptr< HybridArcPropagatorSettings< > > propagatorSettings =
            boost::make_shared< HybridArcPropagatorSettings< > >(
                singleArcPropagatorSettings, multiArcPropagatorSettings );

    boost::shared_ptr< IntegratorSettings< > > integratorSettings =
            boost::make_shared< IntegratorSettings< > >
            ( rungeKutta4, initialEphemerisTime, 60.0 );



    // Set parameters that are to be estimated.
    std::vector< boost::shared_ptr< EstimatableParameterSettings > > parameterNames;
    parameterNames.push_back(
                boost::make_shared< InitialTranslationalStateEstimatableParameterSettings< StateScalarType > >(
                    singleArcBodiesToIntegrate.at( 0 ), singleArcInitialStates, singleArcCentralBodies.at( 0 ) ) );
    parameterNames.push_back(
                boost::make_shared< ArcWiseInitialTranslationalStateEstimatableParameterSettings< StateScalarType > >(
                    multiArcBodiesToIntegrate.at( 0 ), multiArcPropagatorSettings->getInitialStates( ),
                    integrationArcStarts, multiArcCentralBodies.at( 0 ) ) );
    parameterNames.push_back( boost::make_shared< EstimatableParameterSettings >( "Sun", gravitational_parameter ) );
    parameterNames.push_back( boost::make_shared< EstimatableParameterSettings >( "Mars", gravitational_parameter ) );

    boost::shared_ptr< estimatable_parameters::EstimatableParameterSet< StateScalarType > > parametersToEstimate =
            createParametersToEstimate< StateScalarType >( parameterNames, bodyMap );


    // Define links in simulation.
    std::vector< LinkEnds > linkEnds2;
    linkEnds2.resize( 3 );
    linkEnds2[ 0 ][ transmitter ] = grazStation;
    linkEnds2[ 0 ][ receiver ] = mslStation;

    linkEnds2[ 1 ][ receiver ] = grazStation;
    linkEnds2[ 1 ][ transmitter ] = mslStation;

    linkEnds2[ 2 ][ transmitter ] = std::make_pair( "Orbiter", "" );
    linkEnds2[ 2 ][ receiver ] = mslStation;

    observation_models::ObservationSettingsMap observationSettingsMap;
    observationSettingsMap.insert( std::make_pair( linkEnds2[ 0 ], boost::make_shared< ObservationSettings >(
                                       one_way_range ) ) );
    observationSettingsMap.insert( std::make_pair( linkEnds2[ 1 ], boost::make_shared< ObservationSettings >(
                                       one_way_range ) ) );
    observationSettingsMap.insert( std::make_pair( linkEnds2[ 2 ], boost::make_shared< ObservationSettings >(
                                       one_way_range ) ) );


    // Create orbit determination object.
    OrbitDeterminationManager< ObservationScalarType, TimeType > orbitDeterminationManager =
            OrbitDeterminationManager< ObservationScalarType, TimeType >(
                bodyMap, parametersToEstimate,
                observationSettingsMap, integratorSettings, propagatorSettings );

    Eigen::Matrix< StateScalarType, Eigen::Dynamic, 1 > initialParameterEstimate =
            parametersToEstimate->template getFullParameterValues< StateScalarType >( );


    TimeType observationTime;
    int numberOfObservationsPerArc = 5000;
    double timeBuffer = 9000.0;


    std::vector< TimeType > initialObservationTimes;
    initialObservationTimes.resize( numberOfObservationsPerArc * integrationArcStarts.size( ) );

    for( unsigned int i = 0; i < integrationArcLimits.size( ) - 1; i++ )
    {
        double currentTimeStep = ( integrationArcLimits[ i + 1 ] - integrationArcLimits[ i ] - 2.0 * timeBuffer ) /
                static_cast< double >( numberOfObservationsPerArc - 1 );
        observationTime = integrationArcLimits[ i ] + timeBuffer;
        for( int j = 0; j < numberOfObservationsPerArc; j++ )
        {
            initialObservationTimes[ j + i * numberOfObservationsPerArc ] = observationTime;
            observationTime += currentTimeStep;
        }
    }

    std::map< ObservableType, std::map< LinkEnds, std::pair< std::vector< TimeType >, LinkEndType > > > measurementSimulationInput;
    std::map< LinkEnds, std::pair< std::vector< TimeType >, LinkEndType > > singleObservableSimulationInput;


    singleObservableSimulationInput[ linkEnds2[ 0 ] ] = std::make_pair( initialObservationTimes, receiver );
    measurementSimulationInput[ one_way_range ] = singleObservableSimulationInput;
    singleObservableSimulationInput[ linkEnds2[ 1 ] ] = std::make_pair( initialObservationTimes, receiver );
    measurementSimulationInput[ one_way_range ] = singleObservableSimulationInput;
    singleObservableSimulationInput[ linkEnds2[ 2 ] ] = std::make_pair( initialObservationTimes, receiver );
    measurementSimulationInput[ one_way_range ] = singleObservableSimulationInput;


    typedef Eigen::Matrix< ObservationScalarType, Eigen::Dynamic, 1 > ObservationVectorType;
    typedef std::map< LinkEnds, std::pair< ObservationVectorType, std::pair< std::vector< TimeType >, LinkEndType > > > SingleObservablePodInputType;
    typedef std::map< ObservableType, SingleObservablePodInputType > PodInputDataType;

    PodInputDataType observationsAndTimes = simulateObservations< ObservationScalarType, TimeType >(
                measurementSimulationInput, orbitDeterminationManager.getObservationSimulators( )  );

    Eigen::Matrix< StateScalarType, Eigen::Dynamic, 1 > truthParameters = initialParameterEstimate;

    for( unsigned int i = 0; i < ( truthParameters.rows( ) - 2 ) / 6; i++ )
    {
        if( i == 0 )
        {
            initialParameterEstimate[ 0 + 6 * i ] += 1.0E6;
            initialParameterEstimate[ 1 + 6 * i ] += 1.0E6;
            initialParameterEstimate[ 2 + 6 * i ] += 1.0E6;
            initialParameterEstimate[ 3 + 6 * i ] += 1.0E0;
            initialParameterEstimate[ 4 + 6 * i ] += 1.0E0;
            initialParameterEstimate[ 5 + 6 * i ] += 1.0E0;
        }
        else
        {
            initialParameterEstimate[ 0 + 6 * i ] += 1.0E0;
            initialParameterEstimate[ 1 + 6 * i ] += 1.0E0;
            initialParameterEstimate[ 2 + 6 * i ] += 1.0E0;
            initialParameterEstimate[ 3 + 6 * i ] += 1.0E-5;
            initialParameterEstimate[ 4 + 6 * i ] += 1.0E-5;
            initialParameterEstimate[ 5 + 6 * i ] += 1.0E-5;
        }
    }
    for( unsigned int i = truthParameters.rows( ) - 2; i < initialParameterEstimate.rows( ); i++ )
    {
        initialParameterEstimate[ i ] *= ( 1.0 + 1.0E-6 );
    }

    parametersToEstimate->resetParameterValues( initialParameterEstimate );

    boost::shared_ptr< PodInput< ObservationScalarType, TimeType > > podInput =
            boost::make_shared< PodInput< ObservationScalarType, TimeType > >(
                observationsAndTimes, ( initialParameterEstimate ).rows( ) );
    boost::shared_ptr< PodOutput< StateScalarType > > podOutput = orbitDeterminationManager.estimateParameters(
                podInput );

    return ( podOutput->parameterEstimate_ - truthParameters ).template cast< double >( );
}


BOOST_AUTO_TEST_CASE( test_HybridArcStateEstimation )
{
    // Execute test for linked arcs and separate arcs.
    for( unsigned int testCase = 0; testCase < 2; testCase++ )
    {
        Eigen::VectorXd parameterError = executeParameterEstimation< double, double, double >(
                    testCase );
        int numberOfEstimatedArcs = ( parameterError.rows( ) - 8 ) / 6;

        std::cout<<"Final parameter error: "<<std::endl<<parameterError.transpose( )<<std::endl;

        BOOST_CHECK_SMALL( std::fabs( parameterError( 0 ) ), 0.5 );
        BOOST_CHECK_SMALL( std::fabs( parameterError( 1 ) ), 0.5 );
        BOOST_CHECK_SMALL( std::fabs( parameterError( 2 ) ), 10.0 );

        BOOST_CHECK_SMALL( std::fabs( parameterError( 3 ) ), 5E-6 );
        BOOST_CHECK_SMALL( std::fabs( parameterError( 4 ) ), 5E-6 );
        BOOST_CHECK_SMALL( std::fabs( parameterError( 5 ) ), 4E-4 );


        for( int i = 1; i <= numberOfEstimatedArcs; i++ )
        {
            for( unsigned int j = 0; j < 3; j++ )
            {
                BOOST_CHECK_SMALL( std::fabs( parameterError( i * 6 + j ) ), 5E-4 );
                BOOST_CHECK_SMALL( std::fabs( parameterError( i * 6 + j + 3 ) ), 5.0E-7  );
            }
        }

        BOOST_CHECK_SMALL( std::fabs( parameterError( parameterError.rows( ) - 2 ) ), 1.0E10 );
        BOOST_CHECK_SMALL( std::fabs( parameterError( parameterError.rows( ) - 1 ) ), 1.0E2 );
    }

}

BOOST_AUTO_TEST_SUITE_END( )

}

}
