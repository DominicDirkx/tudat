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
#include <string>
#include <vector>

#include <boost/test/unit_test.hpp>
#include <boost/make_shared.hpp>
#include <boost/lambda/lambda.hpp>

#include "Tudat/Basics/testMacros.h"

#include "Tudat/InputOutput/basicInputOutput.h"
#include "Tudat/External/SpiceInterface/spiceInterface.h"

#include "Tudat/Astrodynamics/BasicAstrodynamics/sphericalBodyShapeModel.h"
#include "Tudat/Astrodynamics/Ephemerides/constantEphemeris.h"
#include "Tudat/Astrodynamics/OrbitDetermination/EstimatableParameters/constantRotationRate.h"
#include "Tudat/Astrodynamics/OrbitDetermination/EstimatableParameters/constantRotationalOrientation.h"
#include "Tudat/SimulationSetup/EstimationSetup/createCartesianStatePartials.h"
#include "Tudat/Astrodynamics/OrbitDetermination/ObservationPartials/UnitTests/numericalObservationPartial.h"
#include "Tudat/SimulationSetup/EnvironmentSetup/createGroundStations.h"
#include "Tudat/SimulationSetup/EnvironmentSetup/createBodies.h"
#include "Tudat/SimulationSetup/EnvironmentSetup/defaultBodies.h"
#include "Tudat/SimulationSetup/EstimationSetup/createEstimatableParameters.h"


namespace tudat
{
namespace unit_tests
{

using namespace tudat::basic_astrodynamics;
using namespace tudat::ground_stations;
using namespace tudat::gravitation;
using namespace tudat::ephemerides;
using namespace tudat::observation_models;
using namespace tudat::simulation_setup;
using namespace tudat::spice_interface;
using namespace tudat::orbit_determination;
using namespace tudat::estimatable_parameters;
using namespace tudat::observation_partials;

BOOST_AUTO_TEST_SUITE( test_position_partials)

//! Test partial derivatives of positions w.r.t. parameters. Partials of most observables are computed in terms of these
//! partials
BOOST_AUTO_TEST_CASE( testCartesianStatePartials )
{
    //Load spice kernels.
    std::string kernelsPath = input_output::getSpiceKernelPath( );
    loadSpiceKernelInTudat( kernelsPath + "de-403-masses.tpc");
    loadSpiceKernelInTudat( kernelsPath + "de421.bsp");
    loadSpiceKernelInTudat( kernelsPath + "naif0009.tls");
    loadSpiceKernelInTudat( kernelsPath + "pck00009.tpc");

    // Specify initial time
    double initialEphemerisTime = 1.0E7;

    // Create bodies.
    NamedBodyMap bodyMap;
    bodyMap[ "Earth" ] = boost::make_shared< Body >( );
    bodyMap[ "Moon" ] = boost::make_shared< Body >( );
    bodyMap[ "Sun" ] = boost::make_shared< Body >( );
    bodyMap[ "Mars" ] = boost::make_shared< Body >( );

    // Define properties of bodies
    bodyMap[ "Earth" ]->setShapeModel( boost::make_shared< SphericalBodyShapeModel >(
                                           spice_interface::getAverageRadius( "Earth" ) ) );
    bodyMap[ "Mars" ]->setShapeModel( boost::make_shared< SphericalBodyShapeModel >(
                                          spice_interface::getAverageRadius( "Mars" ) ) );

    bodyMap[ "Earth" ]->setEphemeris(
                boost::make_shared< ConstantEphemeris >(
                    getBodyCartesianStateAtEpoch( "Earth", "SSB", "ECLIPJ2000", "NONE", 0.0 ),
                    "SSB", "ECLIPJ2000" ) );
    bodyMap[ "Moon" ]->setEphemeris(
                boost::make_shared< ConstantEphemeris >(
                    getBodyCartesianStateAtEpoch( "Moon", "SSB", "ECLIPJ2000", "NONE", 0.0 ),
                    "SSB", "ECLIPJ2000" ) );
    bodyMap[ "Sun" ]->setEphemeris(
                boost::make_shared< ConstantEphemeris >(
                    getBodyCartesianStateAtEpoch( "Sun", "SSB", "ECLIPJ2000", "NONE", 0.0 ),
                    "SSB", "ECLIPJ2000" ) );
    bodyMap[ "Mars" ]->setEphemeris(
                boost::make_shared< ConstantEphemeris >(
                    getBodyCartesianStateAtEpoch( "Mars", "SSB", "ECLIPJ2000", "NONE", 0.0 ),
                    "SSB", "ECLIPJ2000" ) );

    bodyMap[ "Sun" ]->setGravityFieldModel(
                boost::make_shared< GravityFieldModel >( getBodyGravitationalParameter( "Sun" ) ) );
    bodyMap[ "Moon" ]->setGravityFieldModel(
                boost::make_shared< GravityFieldModel >( getBodyGravitationalParameter( "Moon" ) ) );
    bodyMap[ "Earth" ]->setGravityFieldModel(
                boost::make_shared< GravityFieldModel >( getBodyGravitationalParameter( "Earth" ) ) );



    bodyMap[ "Earth" ]->setRotationalEphemeris(
                createRotationModel(
                    boost::make_shared< SimpleRotationModelSettings >(
                        "ECLIPJ2000", "IAU_Earth",
                        spice_interface::computeRotationQuaternionBetweenFrames(
                            "ECLIPJ2000", "IAU_Earth", initialEphemerisTime ),
                        initialEphemerisTime, 2.0 * mathematical_constants::PI / physical_constants::JULIAN_DAY ), "Earth" ) );
    bodyMap[ "Mars" ]->setRotationalEphemeris(
                createRotationModel(
                    boost::make_shared< SimpleRotationModelSettings >(
                        "ECLIPJ2000", "IAU_Mars",
                        spice_interface::computeRotationQuaternionBetweenFrames(
                            "ECLIPJ2000", "IAU_Mars", initialEphemerisTime ),
                        initialEphemerisTime, 2.0 * mathematical_constants::PI / physical_constants::JULIAN_DAY ), "Mars" ) );

    // Finalize body creation.
    setGlobalFrameBodyEphemerides( bodyMap, "SSB", "ECLIPJ2000" );


    // Create ground stations
    std::map< std::pair< std::string, std::string >, Eigen::Vector3d > groundStationsToCreate;
    groundStationsToCreate[ std::make_pair( "Earth", "Graz" ) ] =
            ( Eigen::Vector3d( ) << 1.7E6, -6.2E6, 1.3E5 ).finished( );
    groundStationsToCreate[ std::make_pair( "Mars", "MSL" ) ] =
            ( Eigen::Vector3d( ) <<-2.5E5, 3.2E6, -2.65E4 ).finished( );
    createGroundStations( bodyMap, groundStationsToCreate );



    // Define list of ground station names.
    std::vector< std::pair< std::string, std::string > > groundStations;
    groundStations.push_back( std::make_pair( "Earth", "Graz" ) );
    groundStations.push_back( std::make_pair( "Mars", "MSL" ) );

    // Create link ends set.
    LinkEnds linkEnds;
    linkEnds[ observed_body ] = groundStations[ 0 ];

    LinkEnds linkEnds2;
    linkEnds2[ observed_body ] = groundStations[ 1 ];

    boost::shared_ptr< GroundStation > receivingGroundStation =
            bodyMap[ "Earth" ]->getGroundStation( "Graz" );

    //Create parameter objects.
    boost::shared_ptr< RotationRate > earthRotationRate = boost::make_shared< RotationRate >(
                boost::dynamic_pointer_cast< SimpleRotationalEphemeris >(
                    bodyMap[ "Earth" ]->getRotationalEphemeris( ) ), "Earth");
    boost::shared_ptr< ConstantRotationalOrientation > earthPolePosition =
            boost::make_shared< ConstantRotationalOrientation >(
                boost::dynamic_pointer_cast< SimpleRotationalEphemeris >(
                    bodyMap[ "Earth" ]->getRotationalEphemeris( ) ), "Earth" );




    // Create explicit position partial objects.
    boost::shared_ptr< CartesianStatePartial > partialObjectWrtReceiverPosition =
            createCartesianStatePartialsWrtBodyState( linkEnds, bodyMap, "Earth" ).begin( )->second;

    // Create explicit parameter partial objects.
    boost::shared_ptr< CartesianStatePartial > partialObjectWrtReceiverRotationRate =
            createCartesianStatePartialsWrtParameter(
                linkEnds, bodyMap, earthRotationRate ).begin( )->second;
    boost::shared_ptr< CartesianStatePartial > partialObjectWrtReceiverPolePosition =
            createCartesianStatePartialsWrtParameter(
                linkEnds, bodyMap, earthPolePosition ).begin( )->second;

    // Calculate transmission/reception times and states
    Eigen::Vector6d currentState;
    double receptionTime = 1.1E7;
    currentState = bodyMap.at( "Earth" )->getStateInBaseFrameFromEphemeris( receptionTime );

    double currentTime = receptionTime;

    // Compute partials
    Eigen::MatrixXd partialWrtReceiverPosition =
            partialObjectWrtReceiverPosition->calculatePartialOfPosition( currentState, currentTime );

    Eigen::MatrixXd partialWrtReceiverRotationRate =
            partialObjectWrtReceiverRotationRate->calculatePartialOfPosition( currentState, currentTime );
    Eigen::MatrixXd partialOfVelocityWrtReceiverRotationRate =
            partialObjectWrtReceiverRotationRate->calculatePartialOfVelocity( currentState, currentTime );

    Eigen::MatrixXd partialWrtReceiverPolePosition =
            partialObjectWrtReceiverPolePosition->calculatePartialOfPosition( currentState, currentTime );
    Eigen::MatrixXd partialOfVelocityWrtReceiverPolePosition =
            partialObjectWrtReceiverPolePosition->calculatePartialOfVelocity( currentState, currentTime );

    // Define observation function
    boost::function< Eigen::VectorXd( const double ) > observationFunctionAtReception =
            boost::bind( &Ephemeris::getCartesianState, createReferencePointEphemeris< double, double >(
                             bodyMap.at( "Earth" )->getEphemeris( ), bodyMap.at( "Earth" )->getRotationalEphemeris( ),
                             boost::bind( &GroundStation::getStateInPlanetFixedFrame< double, double >,
                                          bodyMap[ "Earth" ]->getGroundStation( "Graz" ), _1 ) ), _1 );



    // Calculate numerical partials w.r.t. Earth state.
    Eigen::Vector3d bodyPositionVariation;
    bodyPositionVariation << 10.0, 10.0, 10.0;
    boost::shared_ptr< ConstantEphemeris > earthEphemeris = boost::dynamic_pointer_cast< ConstantEphemeris >(
                bodyMap[ "Earth" ]->getEphemeris( ) );
    Eigen::Vector6d earthUnperturbedState = earthEphemeris->getCartesianState( 0.0 );
    Eigen::Vector6d perturbedEarthState;
    Eigen::Matrix< double, 3, 3 > numericalPartialWrtReceiverPosition = Eigen::Matrix< double, 3, 3 >::Zero( );
    for( int i = 0; i < 3; i++ )
    {
        perturbedEarthState = earthUnperturbedState;
        perturbedEarthState( i ) += bodyPositionVariation( i );
        earthEphemeris->updateConstantState( perturbedEarthState );
        Eigen::Vector3d upPerturbedPosition = observationFunctionAtReception( currentTime ).segment( 0, 3 );

        perturbedEarthState = earthUnperturbedState;
        perturbedEarthState( i ) -= bodyPositionVariation( i );
        earthEphemeris->updateConstantState( perturbedEarthState );
        Eigen::Vector3d downPerturbedPosition = observationFunctionAtReception( currentTime ).segment( 0, 3 );

        numericalPartialWrtReceiverPosition.block( 0, i, 3, 1 ) = ( upPerturbedPosition - downPerturbedPosition ) /
                ( 2.0 * bodyPositionVariation( i ) );
    }
    earthEphemeris->updateConstantState( earthUnperturbedState );

    // Test partial w.r.t. position
    TUDAT_CHECK_MATRIX_CLOSE_FRACTION( partialWrtReceiverPosition.block( 0, 0, 3, 3 ),
                                       numericalPartialWrtReceiverPosition, 1.0E-12 );


    // Compute numerical partial w.r.t. rotation rate.
    Eigen::Vector6d numericalPartialWrtReceiverRotationRate = calculateNumericalObservationParameterPartial(
                earthRotationRate, 1.0E-10, observationFunctionAtReception,
                receptionTime );

    // Test partial w.r.t. rotation rate
    TUDAT_CHECK_MATRIX_CLOSE_FRACTION( partialWrtReceiverRotationRate,
                                       numericalPartialWrtReceiverRotationRate.segment( 0, 3 ), 1.0E-5 );
    TUDAT_CHECK_MATRIX_CLOSE_FRACTION( partialOfVelocityWrtReceiverRotationRate,
                                       numericalPartialWrtReceiverRotationRate.segment( 3, 3 ), 1.0E-5 );


    // Compute numerical partial w.r.t. pole position
    Eigen::VectorXd polePositionPerturbation = ( Eigen::Vector2d( )<<1.0E-5, 1.0E-5 ).finished( );
    Eigen::MatrixXd numericalPartialWrtReceiverPolePosition = calculateNumericalObservationParameterPartial(
                earthPolePosition, polePositionPerturbation, observationFunctionAtReception, receptionTime );

    // Test partials w.r.t. pole position (different tolernaces due to different magnitudes of partials).
    TUDAT_CHECK_MATRIX_CLOSE_FRACTION( ( partialWrtReceiverPolePosition.block( 0, 0, 1, 2 ) ),
                                       ( numericalPartialWrtReceiverPolePosition.block( 0, 0, 1, 2 ) ), 1.0E-4 );
    TUDAT_CHECK_MATRIX_CLOSE_FRACTION( ( partialWrtReceiverPolePosition.block( 1, 0, 2, 2 ) ),
                                       ( numericalPartialWrtReceiverPolePosition.block( 1, 0, 2, 2 ) ), 1.0E-6 );

    for( int i = 0; i < 3; i++ )
    {
        for( int j = 0; j < 2; j++ )
        {
            BOOST_CHECK_SMALL( std::fabs( partialOfVelocityWrtReceiverPolePosition( i, j ) -
                                          numericalPartialWrtReceiverPolePosition( i + 3, j ) ), 1.0E-6 );
        }
    }
}

boost::shared_ptr< simulation_setup::RotationModelSettings > getMarsPerturbedRotationSettings( )
{
    std::string dataPath = input_output::getTudatRootPath( ) + "/Astrodynamics/Ephemerides/";
    Eigen::MatrixXd polynomialRotationVariations = input_output::readMatrixFromFile( dataPath + "marsPolynomialRotationVariations.dat" );
    Eigen::MatrixXd periodicRightAscensionVariations = input_output::readMatrixFromFile( dataPath + "marsRightAscensionPeriodicVariations.dat" );
    Eigen::MatrixXd periodicDeclinationVariations = input_output::readMatrixFromFile( dataPath + "marsDeclinationPeriodicVariations.dat" );
    Eigen::MatrixXd periodicPrimeMeridianVariations = input_output::readMatrixFromFile( dataPath + "marsMeridianPeriodicVariations.dat" );

    std::vector< double > rightAscensionPolynomialTerms =
    { polynomialRotationVariations( 0, 0 ), polynomialRotationVariations( 0, 1 ) };

    std::vector< double > declinationPolynomialTerms =
    { polynomialRotationVariations( 1, 0 ), polynomialRotationVariations( 1, 1 ) };

    std::vector< double > primeMeridianPolynomialTerms =
    { polynomialRotationVariations( 2, 0 ), polynomialRotationVariations( 2, 1 ) };

    std::map< double, std::pair< double, double > > rightAscensionLibrations;
    std::map< double, std::pair< double, double > > declinationLibrations;
    std::map< double, std::pair< double, double > > primeMeridianLibrations;

    for( unsigned int i = 0; i < periodicRightAscensionVariations.rows( ); i++ )
    {
        rightAscensionLibrations[ periodicRightAscensionVariations( i, 2 ) ] = std::make_pair(
                    periodicRightAscensionVariations( i, 1 ), periodicRightAscensionVariations( i, 0 ) );
    }

    for( unsigned int i = 0; i < periodicDeclinationVariations.rows( ); i++ )
    {
        declinationLibrations[ periodicDeclinationVariations( i, 2 ) ] = std::make_pair(
                    periodicDeclinationVariations( i, 1 ), periodicDeclinationVariations( i, 0 ) );
    }

    for( unsigned int i = 0; i < periodicPrimeMeridianVariations.rows( ); i++ )
    {
        primeMeridianLibrations[ periodicPrimeMeridianVariations( i, 2 ) ] = std::make_pair(
                    periodicPrimeMeridianVariations( i, 1 ), periodicPrimeMeridianVariations( i, 0 ) );
    }


    Eigen::Matrix3d fromIntermediateFrameToBaseFrame = Eigen::Matrix3d::Identity( );
    std::string originalFrame = "ECLIPJ2000";
    std::string targetFrame = "Mars_Fixed";

    return boost::make_shared< simulation_setup::DirectRotationVariationSettings >(
                rightAscensionPolynomialTerms, declinationPolynomialTerms, primeMeridianPolynomialTerms,
                rightAscensionLibrations, declinationLibrations, primeMeridianLibrations,
                fromIntermediateFrameToBaseFrame, originalFrame, targetFrame  );
}

BOOST_AUTO_TEST_CASE( test_RotationVariationCartesianStatePartials )
{
    //Load spice kernels.
    std::string kernelsPath = input_output::getSpiceKernelPath( );
    loadSpiceKernelInTudat( kernelsPath + "de-403-masses.tpc");
    loadSpiceKernelInTudat( kernelsPath + "de421.bsp");
    loadSpiceKernelInTudat( kernelsPath + "naif0009.tls");
    loadSpiceKernelInTudat( kernelsPath + "pck00009.tpc");

    // Create bodies.
    NamedBodyMap bodyMap;
    bodyMap[ "Earth" ] = boost::make_shared< Body >( );
    bodyMap[ "Sun" ] = boost::make_shared< Body >( );
    bodyMap[ "Mars" ] = boost::make_shared< Body >( );

    // Define properties of bodies
    bodyMap[ "Earth" ]->setShapeModel( boost::make_shared< SphericalBodyShapeModel >(
                                           spice_interface::getAverageRadius( "Earth" ) ) );
    bodyMap[ "Mars" ]->setShapeModel( boost::make_shared< SphericalBodyShapeModel >(
                                          spice_interface::getAverageRadius( "Mars" ) ) );

    bodyMap[ "Earth" ]->setEphemeris(
                boost::make_shared< ConstantEphemeris >(
                    getBodyCartesianStateAtEpoch( "Earth", "SSB", "ECLIPJ2000", "NONE", 0.0 ),
                    "SSB", "ECLIPJ2000" ) );
    bodyMap[ "Sun" ]->setEphemeris(
                boost::make_shared< ConstantEphemeris >(
                    getBodyCartesianStateAtEpoch( "Sun", "SSB", "ECLIPJ2000", "NONE", 0.0 ),
                    "SSB", "ECLIPJ2000" ) );
    bodyMap[ "Mars" ]->setEphemeris(
                boost::make_shared< ConstantEphemeris >(
                    getBodyCartesianStateAtEpoch( "Mars", "SSB", "ECLIPJ2000", "NONE", 0.0 ),
                    "SSB", "ECLIPJ2000" ) );

    std::cout<<"A"<<std::endl;

    bodyMap[ "Mars" ]->setRotationalEphemeris(
                createRotationModel( getMarsPerturbedRotationSettings( ), "Mars" ) );

    // Finalize body creation.
    setGlobalFrameBodyEphemerides( bodyMap, "SSB", "ECLIPJ2000" );


    // Create ground stations
    std::map< std::pair< std::string, std::string >, Eigen::Vector3d > groundStationsToCreate;
    groundStationsToCreate[ std::make_pair( "Mars", "MSL" ) ] =
            ( Eigen::Vector3d( ) <<-2.5E5, 3.2E6, -2.65E4 ).finished( );
    createGroundStations( bodyMap, groundStationsToCreate );



    // Define list of ground station names.
    std::vector< std::pair< std::string, std::string > > groundStations;
    groundStations.push_back( std::make_pair( "Mars", "MSL" ) );

    // Create link ends set.
    LinkEnds linkEnds;
    linkEnds[ observed_body ] = groundStations[ 0 ];


    //Create parameter objects.
    boost::shared_ptr< EstimatableParameterSettings > meridianPeriodicVariationSettings =
           boost::make_shared< RotationModelPeriodicVariationAmplitudeSettings >(
                "Mars" , prime_meridian_angle, boost::assign::list_of( 9040.18653077765 )( 3616.15697077033 ) );
    boost::shared_ptr< EstimatableParameterSettings > rightAscensionPeriodicVariationSettings =
           boost::make_shared< RotationModelPeriodicVariationAmplitudeSettings >(
                "Mars" , right_ascension_angle, boost::assign::list_of( 9040.23287156661 )( 686286398.465473 ) );
    boost::shared_ptr< EstimatableParameterSettings > declinationPeriodicVariationSettings =
           boost::make_shared< RotationModelPeriodicVariationAmplitudeSettings >(
                "Mars" , declination_angle, boost::assign::list_of( 9040.22715682008 )( 18080.9237229844 ) );

    boost::shared_ptr< EstimatableParameterSettings > meridianPolynomialVariationSettings =
           boost::make_shared< RotationModelPolynomialVariationSettings >(
                "Mars" , prime_meridian_angle, boost::assign::list_of( 0 )( 1 ) );
    boost::shared_ptr< EstimatableParameterSettings > rightAscensionPolynomialVariationSettings =
           boost::make_shared< RotationModelPolynomialVariationSettings >(
                "Mars" , right_ascension_angle, boost::assign::list_of( 1 )( 0 ) );
    boost::shared_ptr< EstimatableParameterSettings > declinationPolynomialVariationSettings =
           boost::make_shared< RotationModelPolynomialVariationSettings >(
                "Mars" , declination_angle, boost::assign::list_of( 0 )( 1 ) );



    boost::shared_ptr< EstimatableParameter< Eigen::VectorXd > > meridianPeriodicVariation =
            simulation_setup::createVectorParameterToEstimate( meridianPeriodicVariationSettings, bodyMap );
    boost::shared_ptr< EstimatableParameter< Eigen::VectorXd > > rightAscensionPeriodicVariation =
            simulation_setup::createVectorParameterToEstimate( rightAscensionPeriodicVariationSettings, bodyMap );
    boost::shared_ptr< EstimatableParameter< Eigen::VectorXd > > declinationPeriodicVariation =
            simulation_setup::createVectorParameterToEstimate( declinationPeriodicVariationSettings, bodyMap );

    boost::shared_ptr< EstimatableParameter< Eigen::VectorXd > > meridianPolynomialVariation =
            simulation_setup::createVectorParameterToEstimate( meridianPolynomialVariationSettings, bodyMap );
    boost::shared_ptr< EstimatableParameter< Eigen::VectorXd > > rightAscensionPolynomialVariation =
            simulation_setup::createVectorParameterToEstimate( rightAscensionPolynomialVariationSettings, bodyMap );
    boost::shared_ptr< EstimatableParameter< Eigen::VectorXd > > declinationPolynomialVariation =
            simulation_setup::createVectorParameterToEstimate( declinationPolynomialVariationSettings, bodyMap );

    // Create explicit parameter partial objects.
    boost::shared_ptr< CartesianStatePartial > partialObjectWrtMarsMeridianPeriodicVariation =
            createCartesianStatePartialsWrtParameter(
                linkEnds, bodyMap, meridianPeriodicVariation ).begin( )->second;
    boost::shared_ptr< CartesianStatePartial > partialObjectWrtMarsDeclinationPeriodicVariation =
            createCartesianStatePartialsWrtParameter(
                linkEnds, bodyMap, declinationPeriodicVariation ).begin( )->second;
    boost::shared_ptr< CartesianStatePartial > partialObjectWrtMarsRightAscensionPeriodicVariation =
            createCartesianStatePartialsWrtParameter(
                linkEnds, bodyMap, rightAscensionPeriodicVariation ).begin( )->second;

    boost::shared_ptr< CartesianStatePartial > partialObjectWrtMarsMeridianPolynomialVariation =
            createCartesianStatePartialsWrtParameter(
                linkEnds, bodyMap, meridianPolynomialVariation ).begin( )->second;
    boost::shared_ptr< CartesianStatePartial > partialObjectWrtMarsDeclinationPolynomialVariation =
            createCartesianStatePartialsWrtParameter(
                linkEnds, bodyMap, declinationPolynomialVariation ).begin( )->second;
    boost::shared_ptr< CartesianStatePartial > partialObjectWrtMarsRightAscensionPolynomialVariation =
            createCartesianStatePartialsWrtParameter(
                linkEnds, bodyMap, rightAscensionPolynomialVariation ).begin( )->second;

    // Calculate transmission/reception times and states
    Eigen::Vector6d currentState;
    double receptionTime = 1.1E7;
    currentState = bodyMap.at( "Mars" )->getStateInBaseFrameFromEphemeris( receptionTime );

    double currentTime = receptionTime;

    // Compute partials
    Eigen::MatrixXd partialWrtMarsMeridianPeriodicVariation =
            partialObjectWrtMarsMeridianPeriodicVariation->calculatePartialOfPosition( currentState, currentTime );
    Eigen::MatrixXd partialWrtMarsRightAscensionPeriodicVariation =
            partialObjectWrtMarsRightAscensionPeriodicVariation->calculatePartialOfPosition( currentState, currentTime );
    Eigen::MatrixXd partialWrtMarsDeclinationPeriodicVariation =
            partialObjectWrtMarsDeclinationPeriodicVariation->calculatePartialOfPosition( currentState, currentTime );

    Eigen::MatrixXd partialWrtMarsMeridianPolynomialVariation =
            partialObjectWrtMarsMeridianPolynomialVariation->calculatePartialOfPosition( currentState, currentTime );
    Eigen::MatrixXd partialWrtMarsRightAscensionPolynomialVariation =
            partialObjectWrtMarsRightAscensionPolynomialVariation->calculatePartialOfPosition( currentState, currentTime );
    Eigen::MatrixXd partialWrtMarsDeclinationPolynomialVariation =
            partialObjectWrtMarsDeclinationPolynomialVariation->calculatePartialOfPosition( currentState, currentTime );

    // Define observation function
    boost::function< Eigen::VectorXd( const double ) > observationFunctionAtReception =
            boost::bind( &Ephemeris::getCartesianState, createReferencePointEphemeris< double, double >(
                             bodyMap.at( "Mars" )->getEphemeris( ), bodyMap.at( "Mars" )->getRotationalEphemeris( ),
                             boost::bind( &GroundStation::getStateInPlanetFixedFrame< double, double >,
                                          bodyMap[ "Mars" ]->getGroundStation( "MSL" ), _1 ) ), _1 );



    Eigen::VectorXd amplitudePerturbation = ( Eigen::Vector4d( )<<1.0E-3, 1.0E-3, 1.0E-3, 1.0E-3 ).finished( );

    Eigen::MatrixXd numericalPartialWrtMeridianPeriodicVariationAmplitude = calculateNumericalObservationParameterPartial(
                meridianPeriodicVariation, amplitudePerturbation, observationFunctionAtReception, receptionTime );
    TUDAT_CHECK_MATRIX_CLOSE_FRACTION( ( partialWrtMarsMeridianPeriodicVariation.block( 0, 0, 3, 4 ) ),
                                           ( numericalPartialWrtMeridianPeriodicVariationAmplitude.block( 0, 0, 3, 4 ) ), 1.0E-5 );

    Eigen::MatrixXd numericalPartialWrtDeclinationPeriodicVariationAmplitude = calculateNumericalObservationParameterPartial(
                declinationPeriodicVariation, amplitudePerturbation, observationFunctionAtReception, receptionTime );
    TUDAT_CHECK_MATRIX_CLOSE_FRACTION( ( partialWrtMarsDeclinationPeriodicVariation.block( 0, 0, 3, 4 ) ),
                                           ( numericalPartialWrtDeclinationPeriodicVariationAmplitude.block( 0, 0, 3, 4 ) ), 1.0E-5 );

    Eigen::MatrixXd numericalPartialWrtRightAscensionPeriodicVariationAmplitude = calculateNumericalObservationParameterPartial(
                rightAscensionPeriodicVariation, amplitudePerturbation, observationFunctionAtReception, receptionTime );
    TUDAT_CHECK_MATRIX_CLOSE_FRACTION( ( partialWrtMarsRightAscensionPeriodicVariation.block( 0, 0, 3, 4 ) ),
                                           ( numericalPartialWrtRightAscensionPeriodicVariationAmplitude.block( 0, 0, 3, 4 ) ), 1.0E-5 );


    int parameterSize = 2;
    amplitudePerturbation = ( Eigen::Vector2d( )<<1.0E-3, 1.0E-10 ).finished( );
    Eigen::MatrixXd numericalPartialWrtMeridianPolynomialVariationAmplitude = calculateNumericalObservationParameterPartial(
                meridianPolynomialVariation, amplitudePerturbation, observationFunctionAtReception, receptionTime );
    TUDAT_CHECK_MATRIX_CLOSE_FRACTION( ( partialWrtMarsMeridianPolynomialVariation.block( 0, 0, 3, parameterSize ) ),
                                           ( numericalPartialWrtMeridianPolynomialVariationAmplitude.block( 0, 0, 3, parameterSize ) ), 1.0E-5 );

    amplitudePerturbation = ( Eigen::Vector2d( )<<1.0E-5, 1.0E-12 ).finished( );
    Eigen::MatrixXd numericalPartialWrtDeclinationPolynomialVariationAmplitude = calculateNumericalObservationParameterPartial(
                declinationPolynomialVariation, amplitudePerturbation, observationFunctionAtReception, receptionTime );
    TUDAT_CHECK_MATRIX_CLOSE_FRACTION( ( partialWrtMarsDeclinationPolynomialVariation.block( 0, 0, 3, parameterSize ) ),
                                           ( numericalPartialWrtDeclinationPolynomialVariationAmplitude.block( 0, 0, 3, parameterSize ) ), 1.0E-5 );

    amplitudePerturbation = ( Eigen::Vector2d( )<<1.0E-12, 1.0E-5 ).finished( );
    Eigen::MatrixXd numericalPartialWrtRightAscensionPolynomialVariationAmplitude = calculateNumericalObservationParameterPartial(
                rightAscensionPolynomialVariation, amplitudePerturbation, observationFunctionAtReception, receptionTime );
    TUDAT_CHECK_MATRIX_CLOSE_FRACTION( ( partialWrtMarsRightAscensionPolynomialVariation.block( 0, 0, 3, parameterSize ) ),
                                           ( numericalPartialWrtRightAscensionPolynomialVariationAmplitude.block( 0, 0, 3, parameterSize ) ), 1.0E-5 );


}


BOOST_AUTO_TEST_SUITE_END( )

} // namespace unit_tests

} // namespace tudat




