/*    Copyright (c) 2010-2018, Delft University of Technology
 *    All rigths reserved
 *
 *    This file is part of the Tudat. Redistribution and use in source and
 *    binary forms, with or without modification, are permitted exclusively
 *    under the terms of the Modified BSD license. You should have received
 *    a copy of the license with this file. If not, please or visit:
 *    http://tudat.tudelft.nl/LICENSE.
 *
 *    References
 *      Wakker, K. F. (2007), Lecture Notes Astrodynamics II (Chapter 18), TU Delft course AE4-874,
 *          Delft University of technology, Delft, The Netherlands.
 *
 */

#define BOOST_TEST_MAIN

#include <boost/test/floating_point_comparison.hpp>
#include <boost/test/unit_test.hpp>
#include <Eigen/Dense>
#include <math.h>
#include <iostream>

#include "Tudat/Astrodynamics/BasicAstrodynamics/physicalConstants.h"
#include "Tudat/Mathematics/BasicMathematics/mathematicalConstants.h"
#include "Tudat/Mathematics/BasicMathematics/basicMathematicsFunctions.h"
#include "Tudat/Mathematics/BasicMathematics/coordinateConversions.h"
#include "Tudat/Astrodynamics/Ephemerides/approximatePlanetPositions.h"
#include "Tudat/SimulationSetup/tudatSimulationHeader.h"
#include "Tudat/Astrodynamics/BasicAstrodynamics/celestialBodyConstants.h"
#include "Tudat/Astrodynamics/ShapeBasedMethods/baseFunctionsSphericalShaping.h"
#include "Tudat/Astrodynamics/ShapeBasedMethods/compositeFunctionSphericalShaping.h"
#include "Tudat/Astrodynamics/ShapeBasedMethods/sphericalShaping.h"

namespace tudat
{
namespace unit_tests
{

//! Test spherical shaping implementation.
BOOST_AUTO_TEST_SUITE( test_spherical_shaping )

//! Test.
BOOST_AUTO_TEST_CASE( test_spherical_shaping_earth_mars_transfer )
{
    spice_interface::loadStandardSpiceKernels( );

    int numberOfRevolutions = 1;
    double julianDate = 8174.5 * physical_constants::JULIAN_DAY;
    double  timeOfFlight = 580.0;

    // Create central, departure and arrival bodies.
    std::vector< std::string > bodiesToCreate;
    bodiesToCreate.push_back( "Sun" );

    std::map< std::string, std::shared_ptr< simulation_setup::BodySettings > > bodySettings =
            simulation_setup::getDefaultBodySettings( bodiesToCreate );

    std::string frameOrigin = "SSB";
    std::string frameOrientation = "ECLIPJ2000";

    // Define central body ephemeris settings.
    bodySettings[ "Sun" ]->ephemerisSettings = std::make_shared< simulation_setup::ConstantEphemerisSettings >(
                ( Eigen::Vector6d( ) << 0.0, 0.0, 0.0, 0.0, 0.0, 0.0 ).finished( ), frameOrigin, frameOrientation );

    bodySettings[ "Sun" ]->ephemerisSettings->resetFrameOrientation( frameOrientation );
    bodySettings[ "Sun" ]->rotationModelSettings->resetOriginalFrame( frameOrientation );

    // Create body map.
    simulation_setup::NamedBodyMap bodyMap = createBodies( bodySettings );

    bodyMap[ "Vehicle" ] = std::make_shared< simulation_setup::Body >( );
    bodyMap.at( "Vehicle" )->setEphemeris( std::make_shared< ephemerides::TabulatedCartesianEphemeris< > >(
                                                         std::shared_ptr< interpolators::OneDimensionalInterpolator
                                                         < double, Eigen::Vector6d > >( ), frameOrigin, frameOrientation ) );

    setGlobalFrameBodyEphemerides( bodyMap, frameOrigin, frameOrientation );

    // Set vehicle mass.
    bodyMap[ "Vehicle" ]->setConstantBodyMass( 400.0 );

    // Define integrator settings.
    std::shared_ptr< numerical_integrators::IntegratorSettings< double > > integratorSettings =
            std::make_shared< numerical_integrators::IntegratorSettings< double > > (
                numerical_integrators::rungeKutta4, 0.0, timeOfFlight * tudat::physical_constants::JULIAN_DAY / 500.0 );


    // Ephemeris departure body.
    ephemerides::EphemerisPointer pointerToDepartureBodyEphemeris = std::make_shared< ephemerides::ApproximatePlanetPositions>(
                ephemerides::ApproximatePlanetPositionsBase::BodiesWithEphemerisData::earthMoonBarycenter );

    // Ephemeris arrival body.
    ephemerides::EphemerisPointer pointerToArrivalBodyEphemeris = std::make_shared< ephemerides::ApproximatePlanetPositions >(
                ephemerides::ApproximatePlanetPositionsBase::BodiesWithEphemerisData::mars );

//    Eigen::VectorXd radialFunctionCoefficients = ( Eigen::Vector7d() << 1.0, 1.0, 1.0, 1.0, 1.0, 1.0, 1.0 ).finished();
//    Eigen::Vector4d elevationFunctionCoefficients = ( Eigen::Vector4d() << 1.0, 1.0, 1.0, 1.0 ).finished();

    Eigen::Vector6d initialState = pointerToDepartureBodyEphemeris->getCartesianState( julianDate );
    Eigen::Vector6d finalState = pointerToArrivalBodyEphemeris->getCartesianState(
                julianDate + timeOfFlight * tudat::physical_constants::JULIAN_DAY );

    // Define root finder settings (used to update the updated value of the free coefficient, so that it matches the required time of flight).
        std::shared_ptr< root_finders::RootFinderSettings > rootFinderSettings =
                std::make_shared< root_finders::RootFinderSettings >( root_finders::bisection_root_finder, 1.0e-6, 30 );

    // Compute shaped trajectory.
    shape_based_methods::SphericalShaping sphericalShaping = shape_based_methods::SphericalShaping(
                initialState, finalState, timeOfFlight * tudat::physical_constants::JULIAN_DAY,
                numberOfRevolutions, bodyMap, "Vehicle", "Sun", 0.000703,
                rootFinderSettings, 1.0e-6, 1.0e-1, integratorSettings );

    // Compute step size.
    double stepSize = ( sphericalShaping.getFinalAzimuthAngle() - sphericalShaping.getInitialAzimuthAngle() ) / 5000.0;

    // Initialise peak acceleration.
    double peakThrustAcceleration = 0.0;

    // Check that the trajectory is feasible, ie curved toward the central body.
    for ( int i = 0 ; i <= 5000 ; i++ )
    {
        double currentThetaAngle = sphericalShaping.getInitialAzimuthAngle() + i * stepSize;

        if ( sphericalShaping.computeCurrentThrustAccelerationVector( currentThetaAngle ).norm() >  peakThrustAcceleration )
        {
            peakThrustAcceleration = sphericalShaping.computeCurrentThrustAccelerationVector( currentThetaAngle ).norm();
        }
    }


    // Check results consistency w.r.t. thesis from T. Roegiers (ADD PROPER REFERENCE)

    double expectedDeltaV = 5700.0;
    double expectedPeakAcceleration = 2.4e-4;

    // DeltaV provided with a precision of 0.1 km/s
    BOOST_CHECK_SMALL( std::fabs(  sphericalShaping.computeDeltaV() - expectedDeltaV ), 100.0 );
    // Peak acceleration provided with a precision 1.0e-5 m/s^2
    BOOST_CHECK_SMALL( std::fabs(  peakThrustAcceleration - expectedPeakAcceleration ), 1e-5 );

}



//! Test.
BOOST_AUTO_TEST_CASE( test_spherical_shaping_earth_1989ML_transfer )
{
    spice_interface::loadStandardSpiceKernels( );

    int numberOfRevolutions = 1;
    double julianDate = 7799.5 * physical_constants::JULIAN_DAY;
    double timeOfFlight = 600.0;

    // Create central, departure and arrival bodies.
    std::vector< std::string > bodiesToCreate;
    bodiesToCreate.push_back( "Sun" );

    std::map< std::string, std::shared_ptr< simulation_setup::BodySettings > > bodySettings =
            simulation_setup::getDefaultBodySettings( bodiesToCreate );

    std::string frameOrigin = "SSB";
    std::string frameOrientation = "ECLIPJ2000";

    // Define central body ephemeris settings.
    bodySettings[ "Sun" ]->ephemerisSettings = std::make_shared< simulation_setup::ConstantEphemerisSettings >(
                ( Eigen::Vector6d( ) << 0.0, 0.0, 0.0, 0.0, 0.0, 0.0 ).finished( ), frameOrigin, frameOrientation );

    bodySettings[ "Sun" ]->ephemerisSettings->resetFrameOrientation( frameOrientation );
    bodySettings[ "Sun" ]->rotationModelSettings->resetOriginalFrame( frameOrientation );

    // Create body map.
    simulation_setup::NamedBodyMap bodyMap = createBodies( bodySettings );

    bodyMap[ "Vehicle" ] = std::make_shared< simulation_setup::Body >( );
    bodyMap.at( "Vehicle" )->setEphemeris( std::make_shared< ephemerides::TabulatedCartesianEphemeris< > >(
                                                         std::shared_ptr< interpolators::OneDimensionalInterpolator
                                                         < double, Eigen::Vector6d > >( ), frameOrigin, frameOrientation ) );

    setGlobalFrameBodyEphemerides( bodyMap, frameOrigin, frameOrientation );

    // Set vehicle mass.
    bodyMap[ "Vehicle" ]->setConstantBodyMass( 400.0 );

    // Define integrator settings.
    std::shared_ptr< numerical_integrators::IntegratorSettings< double > > integratorSettings =
            std::make_shared< numerical_integrators::IntegratorSettings< double > > (
                numerical_integrators::rungeKutta4, 0.0, timeOfFlight * tudat::physical_constants::JULIAN_DAY / 500.0 );

    // Ephemeris departure body.
    ephemerides::EphemerisPointer pointerToDepartureBodyEphemeris = std::make_shared< ephemerides::ApproximatePlanetPositions>(
                ephemerides::ApproximatePlanetPositionsBase::BodiesWithEphemerisData::earthMoonBarycenter );


    Eigen::Vector6d initialState = pointerToDepartureBodyEphemeris->getCartesianState( julianDate );

    // Final state derived from ML1989 ephemeris (from Spice).
    Eigen::Vector6d finalState = ( Eigen::Vector6d() <<
                                            1.197701029846094E+00 * physical_constants::ASTRONOMICAL_UNIT,
                                            1.653518856610793E-01 * physical_constants::ASTRONOMICAL_UNIT,
                                            - 9.230177854743750E-02 * physical_constants::ASTRONOMICAL_UNIT,
                                            - 4.891080912584867E-05 * physical_constants::ASTRONOMICAL_UNIT / physical_constants::JULIAN_DAY,
                                            1.588950249593135E-02 * physical_constants::ASTRONOMICAL_UNIT / physical_constants::JULIAN_DAY,
                                            - 2.980245580772588E-04 * physical_constants::ASTRONOMICAL_UNIT / physical_constants::JULIAN_DAY ).finished();


    // Define root finder settings (used to update the updated value of the free coefficient, so that it matches the required time of flight).
        std::shared_ptr< root_finders::RootFinderSettings > rootFinderSettings =
                std::make_shared< root_finders::RootFinderSettings >( root_finders::bisection_root_finder, 1.0e-6, 30 );

    // Compute shaped trajectory.
    shape_based_methods::SphericalShaping sphericalShaping = shape_based_methods::SphericalShaping(
                initialState, finalState, timeOfFlight * tudat::physical_constants::JULIAN_DAY,
                numberOfRevolutions, bodyMap, "Vehicle", "Sun", - 0.0000703,
                rootFinderSettings, - 1.0e-2, 1.0e-2, integratorSettings );


    // Compute step size.
    double numberOfSteps = 5000.0;
    double stepSize = ( sphericalShaping.getFinalAzimuthAngle() - sphericalShaping.getInitialAzimuthAngle() ) / numberOfSteps;

    // Initialise peak acceleration.
    double peakThrustAcceleration = 0.0;

    // Compute peak acceleration.
    for ( int i = 0 ; i <= numberOfSteps ; i++ )
    {
        double currentThetaAngle = sphericalShaping.getInitialAzimuthAngle() + i * stepSize;

        if ( sphericalShaping.computeCurrentThrustAccelerationVector( currentThetaAngle ).norm() >  peakThrustAcceleration )
        {
            peakThrustAcceleration = sphericalShaping.computeCurrentThrustAccelerationVector( currentThetaAngle ).norm();
        }
    }


    // Check results consistency w.r.t. thesis from T. Roegiers (ADD PROPER REFERENCE)
    // The expected differences are a bit larger than for Earth-Mars transfer due to the higher uncertainty in 1989ML's ephemeris.

    double expectedDeltaV = 4530.0;
    double expectedPeakAcceleration = 1.8e-4;

    // DeltaV provided with a precision of 0.1 km/s
    BOOST_CHECK_SMALL( std::fabs(  sphericalShaping.computeDeltaV() - expectedDeltaV ), 100.0 );
    // Peak acceleration provided with a precision 1.0e-5 m/s^2
    BOOST_CHECK_SMALL( std::fabs(  peakThrustAcceleration - expectedPeakAcceleration ), 1e-5 );

}


//! Test spherical shaping method with various number of revolutions.
BOOST_AUTO_TEST_CASE( test_spherical_shaping_earth_mars_transfer_multi_revolutions )
{
    spice_interface::loadStandardSpiceKernels( );

    double julianDate = 8174.5 * physical_constants::JULIAN_DAY;

    std::vector< int > numberOfRevolutionsVector = { 0, 1, 2 };
    std::vector< double > timeOfFlightVector = { 300.0, 580.0, 750.0 };


    // Create central, departure and arrival bodies.
    std::vector< std::string > bodiesToCreate;
    bodiesToCreate.push_back( "Sun" );

    std::map< std::string, std::shared_ptr< simulation_setup::BodySettings > > bodySettings =
            simulation_setup::getDefaultBodySettings( bodiesToCreate );

    std::string frameOrigin = "SSB";
    std::string frameOrientation = "ECLIPJ2000";

    // Define central body ephemeris settings.
    bodySettings[ "Sun" ]->ephemerisSettings = std::make_shared< simulation_setup::ConstantEphemerisSettings >(
                ( Eigen::Vector6d( ) << 0.0, 0.0, 0.0, 0.0, 0.0, 0.0 ).finished( ), frameOrigin, frameOrientation );

    bodySettings[ "Sun" ]->ephemerisSettings->resetFrameOrientation( frameOrientation );
    bodySettings[ "Sun" ]->rotationModelSettings->resetOriginalFrame( frameOrientation );

    // Create body map.
    simulation_setup::NamedBodyMap bodyMap = createBodies( bodySettings );

    bodyMap[ "Vehicle" ] = std::make_shared< simulation_setup::Body >( );
    bodyMap.at( "Vehicle" )->setEphemeris( std::make_shared< ephemerides::TabulatedCartesianEphemeris< > >(
                                                         std::shared_ptr< interpolators::OneDimensionalInterpolator
                                                         < double, Eigen::Vector6d > >( ), frameOrigin, frameOrientation ) );

    setGlobalFrameBodyEphemerides( bodyMap, frameOrigin, frameOrientation );

    // Set vehicle mass.
    bodyMap[ "Vehicle" ]->setConstantBodyMass( 400.0 );

    // Define integrator settings.
    std::shared_ptr< numerical_integrators::IntegratorSettings< double > > integratorSettings =
            std::make_shared< numerical_integrators::IntegratorSettings< double > > (
                numerical_integrators::rungeKutta4, 0.0, timeOfFlightVector[ 0 ] * tudat::physical_constants::JULIAN_DAY / 500.0 );

    // Ephemeris departure body.
    ephemerides::EphemerisPointer pointerToDepartureBodyEphemeris = std::make_shared< ephemerides::ApproximatePlanetPositions>(
                ephemerides::ApproximatePlanetPositionsBase::BodiesWithEphemerisData::earthMoonBarycenter );

    // Ephemeris arrival body.
    ephemerides::EphemerisPointer pointerToArrivalBodyEphemeris = std::make_shared< ephemerides::ApproximatePlanetPositions >(
                ephemerides::ApproximatePlanetPositionsBase::BodiesWithEphemerisData::mars );

    // Define root finder settings (used to update the updated value of the free coefficient, so that it matches the required time of flight).
        std::shared_ptr< root_finders::RootFinderSettings > rootFinderSettings =
                std::make_shared< root_finders::RootFinderSettings >( root_finders::bisection_root_finder, 1.0e-6, 30 );

    // Lower bounds for the free parameter.
    std::vector< double > freeParameterLowerBoundVector = { -1.0, 1.0e-6, -1.0e-2 };
    // Upper bounds for the free parameter.
    std::vector< double > freeParameterUpperBoundVector = { 5.0e-1, 1.0e-1, 1.0e-2 };

    // Define initial state.
    Eigen::Vector6d initialState = pointerToDepartureBodyEphemeris->getCartesianState( julianDate );

    for ( unsigned int currentTestCase = 0 ; currentTestCase < numberOfRevolutionsVector.size() ; currentTestCase++ )
    {

        // Define final state.
        Eigen::Vector6d finalState = pointerToArrivalBodyEphemeris->getCartesianState(
                    julianDate + timeOfFlightVector[ currentTestCase ] * tudat::physical_constants::JULIAN_DAY );

        // Compute shaped trajectory.
        shape_based_methods::SphericalShaping sphericalShaping = shape_based_methods::SphericalShaping(
                    initialState, finalState, timeOfFlightVector[ currentTestCase ] * tudat::physical_constants::JULIAN_DAY,
                    numberOfRevolutionsVector[ currentTestCase ], bodyMap, "Vehicle", "Sun",
                    0.000703, rootFinderSettings, freeParameterLowerBoundVector[ currentTestCase ],
                    freeParameterUpperBoundVector[ currentTestCase ], integratorSettings );

        // Check consistency of final azimuth angle value with required number of revolutions.
        double initialAzimuthAngle = sphericalShaping.getInitialAzimuthAngle();
        double finalAzimuthAngle = sphericalShaping.getFinalAzimuthAngle();

        double expectedInitialAzimuthAngle = coordinate_conversions::convertCartesianToSphericalState( initialState )[ 1 ];
        if ( expectedInitialAzimuthAngle < 0.0 )
        {
            expectedInitialAzimuthAngle += 2.0 * mathematical_constants::PI;
        }
        double expectedFinalAzimuthAngle = coordinate_conversions::convertCartesianToSphericalState( finalState )[ 1 ];
        if ( expectedFinalAzimuthAngle < 0.0 )
        {
            expectedFinalAzimuthAngle += 2.0 * mathematical_constants::PI;
        }

        if ( expectedFinalAzimuthAngle - expectedInitialAzimuthAngle < 0 )
        {
            expectedFinalAzimuthAngle += 2.0 * mathematical_constants::PI * ( numberOfRevolutionsVector[ currentTestCase ] + 1 );
        }
        else
        {
            expectedFinalAzimuthAngle += 2.0 * mathematical_constants::PI * numberOfRevolutionsVector[ currentTestCase ];
        }

        BOOST_CHECK_SMALL( std::fabs(  initialAzimuthAngle - expectedInitialAzimuthAngle ), 1.0e-15 );
        BOOST_CHECK_SMALL( std::fabs(  finalAzimuthAngle - expectedFinalAzimuthAngle ), 1.0e-15 );

        // Check consistency of expected and calculated states (both at departure and arrival).
        for ( int i = 0 ; i < 6 ; i++ )
        {
            BOOST_CHECK_SMALL( std::fabs( ( initialState[ i ] - sphericalShaping.computeCurrentStateVector( initialAzimuthAngle )[ i ] )
                                          / initialState[ i ] ), 1.0e-12 );
            BOOST_CHECK_SMALL( std::fabs( ( finalState[ i ] - sphericalShaping.computeCurrentStateVector( finalAzimuthAngle )[ i ] )
                                          / finalState[ i ] ), 1.0e-12 );
        }

    }


}


//! Test.
BOOST_AUTO_TEST_CASE( test_spherical_shaping_full_propagation )
{

    int numberOfRevolutions = 1;
    double julianDate = 8174.5 * physical_constants::JULIAN_DAY;
    double  timeOfFlight = 580.0;

    // Ephemeris departure body.
    ephemerides::EphemerisPointer pointerToDepartureBodyEphemeris = std::make_shared< ephemerides::ApproximatePlanetPositions>(
                ephemerides::ApproximatePlanetPositionsBase::BodiesWithEphemerisData::earthMoonBarycenter );

    // Ephemeris arrival body.
    ephemerides::EphemerisPointer pointerToArrivalBodyEphemeris = std::make_shared< ephemerides::ApproximatePlanetPositions >(
                ephemerides::ApproximatePlanetPositionsBase::BodiesWithEphemerisData::mars );

    // Define root finder settings (used to update the updated value of the free coefficient, so that it matches the required time of flight).
        std::shared_ptr< root_finders::RootFinderSettings > rootFinderSettings =
                std::make_shared< root_finders::RootFinderSettings >( root_finders::bisection_root_finder, 1.0e-6, 30 );

    std::map< double, Eigen::VectorXd > fullPropagationResults;
    std::map< double, Eigen::VectorXd > shapingMethodResults;
    std::map< double, Eigen::VectorXd > dependentVariablesHistory;

    spice_interface::loadStandardSpiceKernels( );

    // Create central, departure and arrival bodies.
    std::vector< std::string > bodiesToCreate;
    bodiesToCreate.push_back( "Sun" );
    bodiesToCreate.push_back( "Earth" );
    bodiesToCreate.push_back( "Mars" );
    bodiesToCreate.push_back( "Jupiter" );

    std::map< std::string, std::shared_ptr< simulation_setup::BodySettings > > bodySettings =
            simulation_setup::getDefaultBodySettings( bodiesToCreate );

    std::string frameOrigin = "SSB";
    std::string frameOrientation = "ECLIPJ2000";


    // Define central body ephemeris settings.
    bodySettings[ "Sun" ]->ephemerisSettings = std::make_shared< simulation_setup::ConstantEphemerisSettings >(
                ( Eigen::Vector6d( ) << 0.0, 0.0, 0.0, 0.0, 0.0, 0.0 ).finished( ), frameOrigin, frameOrientation );

    bodySettings[ "Sun" ]->ephemerisSettings->resetFrameOrientation( frameOrientation );
    bodySettings[ "Sun" ]->rotationModelSettings->resetOriginalFrame( frameOrientation );


    // Create body map.
    simulation_setup::NamedBodyMap bodyMap = createBodies( bodySettings );

    bodyMap[ "Vehicle" ] = std::make_shared< simulation_setup::Body >( );
    bodyMap.at( "Vehicle" )->setEphemeris( std::make_shared< ephemerides::TabulatedCartesianEphemeris< > >(
                                                         std::shared_ptr< interpolators::OneDimensionalInterpolator
                                                         < double, Eigen::Vector6d > >( ), frameOrigin, frameOrientation ) );


    setGlobalFrameBodyEphemerides( bodyMap, frameOrigin, frameOrientation );



    std::vector< std::string > bodiesToPropagate;
    bodiesToPropagate.push_back( "Vehicle" );
    std::vector< std::string > centralBodies;
    centralBodies.push_back( "Sun" );

    // Acceleration from the central body.
    std::map< std::string, std::vector< std::shared_ptr< simulation_setup::AccelerationSettings > > > bodyToPropagateAccelerations;
    bodyToPropagateAccelerations[ "Sun" ].push_back( std::make_shared< simulation_setup::AccelerationSettings >(
                                                                basic_astrodynamics::central_gravity ) );

    simulation_setup::SelectedAccelerationMap accelerationMap;
    accelerationMap[ "Vehicle" ] = bodyToPropagateAccelerations;

    // Create the acceleration map.
    basic_astrodynamics::AccelerationMap accelerationModelMap = createAccelerationModelsMap(
                bodyMap, accelerationMap, bodiesToPropagate, centralBodies );

    // Define integrator settings
    std::shared_ptr< numerical_integrators::IntegratorSettings< double > > integratorSettings =
            std::make_shared< numerical_integrators::IntegratorSettings< double > > ( numerical_integrators::rungeKutta4, 0.0,
                                                                                      timeOfFlight * physical_constants::JULIAN_DAY / ( 1000.0 ) );

    // Define mass function of the vehicle.
    std::function< double( const double ) > newMassFunction = [ = ]( const double currentTime )
    {
        return 2000.0 - 50.0 / ( timeOfFlight * physical_constants::JULIAN_DAY ) * currentTime ;
    };
    bodyMap[ "Vehicle" ]->setBodyMassFunction( newMassFunction );


    // Define specific impulse function.
    std::function< double( const double ) > specificImpulseFunction = [ = ]( const double currentTime )
    {
        return 3000.0;
    };


    // Compute shaped trajectory.
    shape_based_methods::SphericalShaping sphericalShaping = shape_based_methods::SphericalShaping(
                pointerToDepartureBodyEphemeris->getCartesianState( julianDate ),
                pointerToArrivalBodyEphemeris->getCartesianState( julianDate + timeOfFlight * tudat::physical_constants::JULIAN_DAY ),
                timeOfFlight * tudat::physical_constants::JULIAN_DAY, numberOfRevolutions, bodyMap, "Vehicle", "Sun", 0.000703,
                rootFinderSettings, 1.0e-6, 1.0e-1, integratorSettings );


    // Define list of dependent variables to save.
    std::vector< std::shared_ptr< propagators::SingleDependentVariableSaveSettings > > dependentVariablesList;
    dependentVariablesList.push_back( std::make_shared< propagators::SingleAccelerationDependentVariableSaveSettings >(
                        basic_astrodynamics::thrust_acceleration, "Vehicle", "Vehicle", 0 ) );

    // Create object with list of dependent variables
    std::shared_ptr< propagators::DependentVariableSaveSettings > dependentVariablesToSave =
            std::make_shared< propagators::DependentVariableSaveSettings >( dependentVariablesList );

    // Create termination conditions settings.
    std::pair< std::shared_ptr< propagators::PropagationTerminationSettings >, std::shared_ptr< propagators::PropagationTerminationSettings > >
            terminationConditions;

    terminationConditions.first = std::make_shared< propagators::PropagationTimeTerminationSettings >( 0.0 );
    terminationConditions.second = std::make_shared< propagators::PropagationTimeTerminationSettings >( timeOfFlight * physical_constants::JULIAN_DAY );

    // Compute halved time of flight.
    double halfOfTimeOfFlight = timeOfFlight * tudat::physical_constants::JULIAN_DAY / 2.0;

    // Compute state at half of the time of flight.
    double azimuthAngleAtHalfTimeOfFlight = sphericalShaping.convertTimeToIndependentVariable( halfOfTimeOfFlight );
    Eigen::Vector6d initialStateAtHalvedTimeOfFlight = sphericalShaping.computeCurrentStateVector( azimuthAngleAtHalfTimeOfFlight );

    // Create low thrust acceleration model.
    std::shared_ptr< propulsion::ThrustAcceleration > lowThrustAccelerationModel =
            sphericalShaping.getLowThrustAccelerationModel( specificImpulseFunction );

    accelerationModelMap[ "Vehicle" ][ "Vehicle" ].push_back( lowThrustAccelerationModel );


    // Create complete propagation settings (backward and forward propagations).
    std::pair< std::shared_ptr< propagators::PropagatorSettings< double > >,
            std::shared_ptr< propagators::PropagatorSettings< double > > > propagatorSettings;


    // Define translational state propagation settings
    std::pair< std::shared_ptr< propagators::TranslationalStatePropagatorSettings< double > >,
            std::shared_ptr< propagators::TranslationalStatePropagatorSettings< double > > > translationalStatePropagatorSettings;

    // Define backward translational state propagation settings.
    translationalStatePropagatorSettings.first = std::make_shared< propagators::TranslationalStatePropagatorSettings< double > >
                        ( centralBodies, accelerationModelMap, bodiesToPropagate, initialStateAtHalvedTimeOfFlight, terminationConditions.first,
                          propagators::cowell, dependentVariablesToSave );

    // Define forward translational state propagation settings.
    translationalStatePropagatorSettings.second = std::make_shared< propagators::TranslationalStatePropagatorSettings< double > >
                        ( centralBodies, accelerationModelMap, bodiesToPropagate,
                          initialStateAtHalvedTimeOfFlight, terminationConditions.second,
                          propagators::cowell, dependentVariablesToSave );

    // Create list of propagation settings.
    std::pair< std::vector< std::shared_ptr< propagators::SingleArcPropagatorSettings< double > > >,
            std::vector< std::shared_ptr< propagators::SingleArcPropagatorSettings< double > > > > propagatorSettingsVector;

    // Backward propagator settings vector.
    propagatorSettingsVector.first.push_back( translationalStatePropagatorSettings.first );

    // Forward propagator settings vector.
    propagatorSettingsVector.second.push_back( translationalStatePropagatorSettings.second );

    // Backward hybrid propagation settings.
    propagatorSettings.first = std::make_shared< propagators::MultiTypePropagatorSettings< double > >( propagatorSettingsVector.first,
                terminationConditions.first, dependentVariablesToSave );

    // Forward hybrid propagation settings.
    propagatorSettings.second = std::make_shared< propagators::MultiTypePropagatorSettings< double > >( propagatorSettingsVector.second,
                terminationConditions.second, dependentVariablesToSave );

    // Compute shaped trajectory and propagated trajectory.
    sphericalShaping.computeSemiAnalyticalAndFullPropagation( specificImpulseFunction, integratorSettings, propagatorSettings,
                                                              fullPropagationResults, shapingMethodResults, dependentVariablesHistory );



    // Check difference between full propagation and shaping method at arrival
    // (disregarding the very last values because of expected interpolation errors).
    int numberOfDisregardedValues = 7;
    std::map< double, Eigen::VectorXd >::iterator itr = fullPropagationResults.end();
    for( int i = 0 ; i < numberOfDisregardedValues ; i++ )
    {
        itr--;
    }

    // Check results consistency between full propagation and shaped trajectory at arrival.
    for ( int i = 0 ; i < 6 ; i++ )
    {
        BOOST_CHECK_SMALL( std::fabs( shapingMethodResults[ itr->first ][ i ] - itr->second[ i ] ) / shapingMethodResults[ itr->first ][ i ] , 1.0e-6 );
    }



    // Check difference between full propagation and shaping method at departure
    // (disregarding the very first values because of expected interpolation errors).
    itr = fullPropagationResults.begin();
    for( int i = 0 ; i < numberOfDisregardedValues ; i++ )
    {
        itr++;
    }

    // Check results consistency between full propagation and shaped trajectory at departure.
    for ( int i = 0 ; i < 6 ; i++ )
    {
        BOOST_CHECK_SMALL( std::fabs( shapingMethodResults[ itr->first ][ i ] - itr->second[ i ] ) / shapingMethodResults[ itr->first ][ i ] , 1.0e-6 );
    }


}


//! Test.
BOOST_AUTO_TEST_CASE( test_spherical_shaping_full_propagation_mass_propagation )
{

    int numberOfRevolutions = 1;
    double julianDate = 8174.5 * physical_constants::JULIAN_DAY;
    double  timeOfFlight = 580.0;


    // Ephemeris departure body.
    ephemerides::EphemerisPointer pointerToDepartureBodyEphemeris = std::make_shared< ephemerides::ApproximatePlanetPositions>(
                ephemerides::ApproximatePlanetPositionsBase::BodiesWithEphemerisData::earthMoonBarycenter );

    // Ephemeris arrival body.
    ephemerides::EphemerisPointer pointerToArrivalBodyEphemeris = std::make_shared< ephemerides::ApproximatePlanetPositions >(
                ephemerides::ApproximatePlanetPositionsBase::BodiesWithEphemerisData::mars );

    // Define root finder settings (used to update the updated value of the free coefficient, so that it matches the required time of flight).
        std::shared_ptr< root_finders::RootFinderSettings > rootFinderSettings =
                std::make_shared< root_finders::RootFinderSettings >( root_finders::bisection_root_finder, 1.0e-6, 30 );


    std::map< double, Eigen::VectorXd > fullPropagationResults;
    std::map< double, Eigen::VectorXd > shapingMethodResults;
    std::map< double, Eigen::VectorXd > dependentVariablesHistory;

    spice_interface::loadStandardSpiceKernels( );

    // Create central, departure and arrival bodies.
    std::vector< std::string > bodiesToCreate;
    bodiesToCreate.push_back( "Sun" );
    bodiesToCreate.push_back( "Earth" );
    bodiesToCreate.push_back( "Mars" );
    bodiesToCreate.push_back( "Jupiter" );

    std::map< std::string, std::shared_ptr< simulation_setup::BodySettings > > bodySettings =
            simulation_setup::getDefaultBodySettings( bodiesToCreate );

    std::string frameOrigin = "SSB";
    std::string frameOrientation = "ECLIPJ2000";


    // Define central body ephemeris settings.
    bodySettings[ "Sun" ]->ephemerisSettings = std::make_shared< simulation_setup::ConstantEphemerisSettings >(
                ( Eigen::Vector6d( ) << 0.0, 0.0, 0.0, 0.0, 0.0, 0.0 ).finished( ), frameOrigin, frameOrientation );

    bodySettings[ "Sun" ]->ephemerisSettings->resetFrameOrientation( frameOrientation );
    bodySettings[ "Sun" ]->rotationModelSettings->resetOriginalFrame( frameOrientation );


    // Create body map.
    simulation_setup::NamedBodyMap bodyMap = createBodies( bodySettings );

    bodyMap[ "Vehicle" ] = std::make_shared< simulation_setup::Body >( );
    bodyMap.at( "Vehicle" )->setEphemeris( std::make_shared< ephemerides::TabulatedCartesianEphemeris< > >(
                                                         std::shared_ptr< interpolators::OneDimensionalInterpolator
                                                         < double, Eigen::Vector6d > >( ), frameOrigin, frameOrientation ) );


    setGlobalFrameBodyEphemerides( bodyMap, frameOrigin, frameOrientation );



    std::vector< std::string > bodiesToPropagate;
    bodiesToPropagate.push_back( "Vehicle" );
    std::vector< std::string > centralBodies;
    centralBodies.push_back( "Sun" );

    // Acceleration from the central body.
    std::map< std::string, std::vector< std::shared_ptr< simulation_setup::AccelerationSettings > > > bodyToPropagateAccelerations;
    bodyToPropagateAccelerations[ "Sun" ].push_back( std::make_shared< simulation_setup::AccelerationSettings >(
                                                                basic_astrodynamics::central_gravity ) );

    simulation_setup::SelectedAccelerationMap accelerationMap;
    accelerationMap[ "Vehicle" ] = bodyToPropagateAccelerations;

    // Create the acceleration map.
    basic_astrodynamics::AccelerationMap accelerationModelMap = createAccelerationModelsMap(
                bodyMap, accelerationMap, bodiesToPropagate, centralBodies );

    // Define integrator settings
    std::shared_ptr< numerical_integrators::IntegratorSettings< double > > integratorSettings =
            std::make_shared< numerical_integrators::IntegratorSettings< double > > ( numerical_integrators::rungeKutta4, 0.0,
                                                                                      timeOfFlight * physical_constants::JULIAN_DAY / ( 1000.0 ) );

    // Set mass of the vehicle.
    bodyMap[ "Vehicle" ]->setConstantBodyMass( 2000.0 );


    // Define specific impulse function.
    std::function< double( const double ) > specificImpulseFunction = [ = ]( const double currentTime )
    {
        return 3000.0;
    };


    // Compute shaped trajectory.
    shape_based_methods::SphericalShaping sphericalShaping = shape_based_methods::SphericalShaping(
                pointerToDepartureBodyEphemeris->getCartesianState( julianDate ),
                pointerToArrivalBodyEphemeris->getCartesianState( julianDate + timeOfFlight * tudat::physical_constants::JULIAN_DAY ),
                timeOfFlight * tudat::physical_constants::JULIAN_DAY,
                numberOfRevolutions, bodyMap, "Vehicle", "Sun", 0.000703,
                rootFinderSettings, 1.0e-6, 1.0e-1, integratorSettings );

    // Define list of dependent variables to save.
    std::vector< std::shared_ptr< propagators::SingleDependentVariableSaveSettings > > dependentVariablesList;
    dependentVariablesList.push_back( std::make_shared< propagators::SingleAccelerationDependentVariableSaveSettings >(
                        basic_astrodynamics::thrust_acceleration, "Vehicle", "Vehicle", 0 ) );
    dependentVariablesList.push_back( std::make_shared< propagators::SingleDependentVariableSaveSettings >(
                    propagators::total_mass_rate_dependent_variables, "Vehicle" ) );

    // Create object with list of dependent variables
    std::shared_ptr< propagators::DependentVariableSaveSettings > dependentVariablesToSave =
            std::make_shared< propagators::DependentVariableSaveSettings >( dependentVariablesList );

    // Create termination conditions settings.
    std::pair< std::shared_ptr< propagators::PropagationTerminationSettings >, std::shared_ptr< propagators::PropagationTerminationSettings > > terminationConditions;

    terminationConditions.first = std::make_shared< propagators::PropagationTimeTerminationSettings >( 0.0 );
    terminationConditions.second = std::make_shared< propagators::PropagationTimeTerminationSettings >( timeOfFlight * physical_constants::JULIAN_DAY );

    // Compute halved time of flight.
    double halfOfTimeOfFlight = timeOfFlight * tudat::physical_constants::JULIAN_DAY / 2.0;

    // Compute state at half of the time of flight.
    double azimuthAngleAtHalfTimeOfFlight = sphericalShaping.convertTimeToIndependentVariable( halfOfTimeOfFlight );
    Eigen::Vector6d initialStateAtHalfOfTimeOfFlight = sphericalShaping.computeCurrentStateVector( azimuthAngleAtHalfTimeOfFlight );

    // Create low thrust acceleration model.
    std::shared_ptr< propulsion::ThrustAcceleration > lowThrustAccelerationModel =
            sphericalShaping.getLowThrustAccelerationModel( specificImpulseFunction );

    accelerationModelMap[ "Vehicle" ][ "Vehicle" ].push_back( lowThrustAccelerationModel );


    // Create complete propagation settings (backward and forward propagations).
    std::pair< std::shared_ptr< propagators::PropagatorSettings< double > >,
            std::shared_ptr< propagators::PropagatorSettings< double > > > propagatorSettings;



    // Define translational state propagation settings
    std::pair< std::shared_ptr< propagators::TranslationalStatePropagatorSettings< double > >,
            std::shared_ptr< propagators::TranslationalStatePropagatorSettings< double > > > translationalStatePropagatorSettings;

    // Define backward translational state propagation settings.
    translationalStatePropagatorSettings.first = std::make_shared< propagators::TranslationalStatePropagatorSettings< double > >
                        ( centralBodies, accelerationModelMap, bodiesToPropagate, initialStateAtHalfOfTimeOfFlight, terminationConditions.first,
                          propagators::cowell, dependentVariablesToSave );

    // Define forward translational state propagation settings.
    translationalStatePropagatorSettings.second = std::make_shared< propagators::TranslationalStatePropagatorSettings< double > >
                        ( centralBodies, accelerationModelMap, bodiesToPropagate, initialStateAtHalfOfTimeOfFlight, terminationConditions.second,
                          propagators::cowell, dependentVariablesToSave );

    // Create mass rate models
    std::map< std::string, std::shared_ptr< basic_astrodynamics::MassRateModel > > massRateModels;
    massRateModels[ "Vehicle" ] = simulation_setup::createMassRateModel( "Vehicle", std::make_shared< simulation_setup::FromThrustMassModelSettings >( 1 ),
                                                       bodyMap, accelerationModelMap );

    double massHalfOfTimeOfFlight = sphericalShaping.computeCurrentMass( halfOfTimeOfFlight, specificImpulseFunction, integratorSettings );

    // Create settings for propagating the mass of the vehicle.
    std::pair< std::shared_ptr< propagators::MassPropagatorSettings< double > >,
            std::shared_ptr< propagators::MassPropagatorSettings< double > > > massPropagatorSettings;

    // Define backward mass propagation settings.
    massPropagatorSettings.first = std::make_shared< propagators::MassPropagatorSettings< double > >(
                bodiesToPropagate, massRateModels, ( Eigen::Matrix< double, 1, 1 >( ) << massHalfOfTimeOfFlight ).finished( ),
                terminationConditions.first );

    // Define forward mass propagation settings.
    massPropagatorSettings.second = std::make_shared< propagators::MassPropagatorSettings< double > >(
                bodiesToPropagate, massRateModels, ( Eigen::Matrix< double, 1, 1 >( ) << massHalfOfTimeOfFlight ).finished( ),
                terminationConditions.second );

    // Create list of propagation settings.
    std::pair< std::vector< std::shared_ptr< propagators::SingleArcPropagatorSettings< double > > >,
            std::vector< std::shared_ptr< propagators::SingleArcPropagatorSettings< double > > > > propagatorSettingsVector;

    // Backward propagator settings vector.
    propagatorSettingsVector.first.push_back( translationalStatePropagatorSettings.first );
    propagatorSettingsVector.first.push_back( massPropagatorSettings.first );

    // Forward propagator settings vector.
    propagatorSettingsVector.second.push_back( translationalStatePropagatorSettings.second );
    propagatorSettingsVector.second.push_back( massPropagatorSettings.second );

    // Backward hybrid propagation settings.
    propagatorSettings.first = std::make_shared< propagators::MultiTypePropagatorSettings< double > >( propagatorSettingsVector.first,
                terminationConditions.first, dependentVariablesToSave );

    // Forward hybrid propagation settings.
    propagatorSettings.second = std::make_shared< propagators::MultiTypePropagatorSettings< double > >( propagatorSettingsVector.second,
                terminationConditions.second, dependentVariablesToSave );


    // Compute shaped trajectory and propagated trajectory.
    sphericalShaping.computeSemiAnalyticalAndFullPropagation( specificImpulseFunction, integratorSettings, propagatorSettings,
                                                              fullPropagationResults, shapingMethodResults, dependentVariablesHistory );



    // Check difference between full propagation and shaping method at arrival
    // (disregarding the very last values because of expected interpolation errors).
    int numberOfDisregardedValues = 7;
    std::map< double, Eigen::VectorXd >::iterator itr = fullPropagationResults.end();
    for( int i = 0 ; i < numberOfDisregardedValues ; i++ )
    {
        itr--;
    }

    // Check results consistency between full propagation and shaped trajectory at arrival.
    for ( int i = 0 ; i < 6 ; i++ )
    {
        BOOST_CHECK_SMALL( std::fabs( shapingMethodResults[ itr->first ][ i ] - itr->second[ i ] ) / shapingMethodResults[ itr->first ][ i ] , 1.0e-6 );
    }



    // Check difference between full propagation and shaping method at departure
    // (disregarding the very first values because of expected interpolation errors).
    itr = fullPropagationResults.begin();
    for( int i = 0 ; i < numberOfDisregardedValues ; i++ )
    {
        itr++;
    }

    // Check results consistency between full propagation and shaped trajectory at departure.
    for ( int i = 0 ; i < 6 ; i++ )
    {
        BOOST_CHECK_SMALL( std::fabs( shapingMethodResults[ itr->first ][ i ] - itr->second[ i ] ) / shapingMethodResults[ itr->first ][ i ] , 1.0e-6 );
    }

    // Check consistency between current and expected mass rates.
    for ( std::map< double, Eigen::VectorXd >::iterator itr = dependentVariablesHistory.begin() ; itr != dependentVariablesHistory.end() ; itr++ )
    {
        Eigen::Vector3d currentThrustAccelerationVector = itr->second.segment( 0, 3 );
        double currentMass = fullPropagationResults.at( itr->first )( 6 );
        double currentMassRate = - itr->second( 3 );
        double expectedMassRate = currentThrustAccelerationVector.norm() * currentMass /
                ( specificImpulseFunction( itr->first ) * physical_constants::SEA_LEVEL_GRAVITATIONAL_ACCELERATION );
        BOOST_CHECK_SMALL( std::fabs( currentMassRate - expectedMassRate ), 1.0e-15 );

    }

}

BOOST_AUTO_TEST_SUITE_END( )

} // namespace unit_tests
} // namespace tudat
