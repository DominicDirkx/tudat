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

#include <chrono>
#include <limits>
#include <boost/make_shared.hpp>
#include <boost/test/unit_test.hpp>
#include <eigen/Eigen/Eigen>


#include <Tudat/SimulationSetup/tudatSimulationHeader.h>


/*
#include "Tudat/Mathematics/BasicMathematics/basicMathematicsFunctions.h"
#include "Tudat/Astrodynamics/BasicAstrodynamics/modifiedEquinoctialElementConversions.h"

#include "Tudat/SimulationSetup/EnvironmentSetup/createAtmosphereModel.h"
#include "Tudat/Astrodynamics/ElectroMagnetism/radiationPressureInterface.h"

#include "Tudat/Astrodynamics/Propagators/DSST/forces/availableForceModels.h"
#include "Tudat/Astrodynamics/Propagators/DSST/utilities/osculatingToMean.h"
*/

namespace tudat
{

namespace unit_tests
{


BOOST_AUTO_TEST_SUITE( test_dsst_propagator )

BOOST_AUTO_TEST_CASE( dsst_propagation )
{

    //////////////////////////////////////////////////////////////////////////////////////////////////////////////////
    //////////////////            USING STATEMENTS              //////////////////////////////////////////////////////
    //////////////////////////////////////////////////////////////////////////////////////////////////////////////////

    using namespace tudat;
    using namespace tudat::simulation_setup;
    using namespace tudat::propagators;
    using namespace tudat::numerical_integrators;
    using namespace tudat::orbital_element_conversions;
    using namespace tudat::basic_mathematics;
    using namespace tudat::gravitation;
    using namespace tudat::numerical_integrators;

    using namespace std::chrono;

    //////////////////////////////////////////////////////////////////////////////////////////////////////////////////
    //////////////////     CREATE ENVIRONMENT AND VEHICLE       //////////////////////////////////////////////////////
    //////////////////////////////////////////////////////////////////////////////////////////////////////////////////

    auto t = steady_clock::now();

    // Load Spice kernels.
    spice_interface::loadSpiceKernelInTudat( input_output::getSpiceKernelPath( ) + "pck00009.tpc" );
    spice_interface::loadSpiceKernelInTudat( input_output::getSpiceKernelPath( ) + "de-403-masses.tpc" );
    spice_interface::loadSpiceKernelInTudat( input_output::getSpiceKernelPath( ) + "de421.bsp" );

    // Set simulation time settings.
    const double simulationStartEpoch = 89242986;
    const double simulationEndEpoch = simulationStartEpoch + 2 * physical_constants::SIDEREAL_YEAR;

    // Define body settings for simulation.
    std::vector< std::string > bodiesToCreate;
    bodiesToCreate.push_back( "Sun" );
    bodiesToCreate.push_back( "Earth" );
    bodiesToCreate.push_back( "Moon" );

    // Create body objects.
    std::map< std::string, boost::shared_ptr< BodySettings > > bodySettings =
            getDefaultBodySettings( bodiesToCreate/*, simulationStartEpoch - 300.0, simulationEndEpoch + 300.0*/ );
    for( unsigned int i = 0; i < bodiesToCreate.size( ); i++ )
    {
        bodySettings[ bodiesToCreate.at( i ) ]->ephemerisSettings->resetFrameOrientation( "J2000" );
        bodySettings[ bodiesToCreate.at( i ) ]->rotationModelSettings->resetOriginalFrame( "J2000" );
    }

    // EARTH
    /*bodySettings[ "Earth" ]->gravityFieldSettings =
            boost::make_shared< GraceGravityFieldSettings >( ggm02c );*/

    bodySettings[ "Earth" ]->atmosphereSettings = boost::make_shared< AtmosphereSettings >( nrlmsise00 );

    // MOON
    bodySettings[ "Moon" ]->gravityFieldSettings = boost::make_shared< GravityFieldSettings >( central_spice );

    NamedBodyMap bodyMap = createBodies( bodySettings );

    //////////////////////////////////////////////////////////////////////////////////////////////////////////////////
    //////////////////             CREATE VEHICLE            /////////////////////////////////////////////////////////
    //////////////////////////////////////////////////////////////////////////////////////////////////////////////////

    // Create spacecraft object.
    bodyMap[ "Asterix" ] = boost::make_shared< simulation_setup::Body >( );
    bodyMap[ "Asterix" ]->setConstantBodyMass( 3000.0 );

    // Create aerodynamic coefficient interface settings.
    double referenceArea = 15.0;
    double aerodynamicCoefficient = 2.2;
    boost::shared_ptr< AerodynamicCoefficientSettings > aerodynamicCoefficientSettings =
            boost::make_shared< ConstantAerodynamicCoefficientSettings >(
                referenceArea, aerodynamicCoefficient * Eigen::Vector3d::UnitX( ), 1, 1 );

    // Create and set aerodynamic coefficients object
    bodyMap[ "Asterix" ]->setAerodynamicCoefficientInterface(
                createAerodynamicCoefficientInterface( aerodynamicCoefficientSettings, "Asterix" ) );

    // Create radiation pressure settings
    double referenceAreaRadiation = 15.0;
    double radiationPressureCoefficient = 1.5;
    std::vector< std::string > occultingBodies;
    // occultingBodies.push_back( "Earth" );
    boost::shared_ptr< RadiationPressureInterfaceSettings > asterixRadiationPressureSettings =
            boost::make_shared< CannonBallRadiationPressureInterfaceSettings >(
                "Sun", referenceAreaRadiation, radiationPressureCoefficient, occultingBodies );

    // Create and set radiation pressure settings
    bodyMap[ "Asterix" ]->setRadiationPressureInterface(
                "Sun", createRadiationPressureInterface(
                    asterixRadiationPressureSettings, "Asterix", bodyMap ) );


    // Finalize body creation.
    setGlobalFrameBodyEphemerides( bodyMap, "SSB", "J2000" );

    //////////////////////////////////////////////////////////////////////////////////////////////////////////////////
    //////////////////            CREATE ACCELERATIONS          //////////////////////////////////////////////////////
    //////////////////////////////////////////////////////////////////////////////////////////////////////////////////

    // Define propagator settings variables.
    SelectedAccelerationMap accelerationMap;
    std::vector< std::string > bodiesToPropagate;
    std::vector< std::string > centralBodies;

    // Define propagation settings.
    std::map< std::string, std::vector< boost::shared_ptr< AccelerationSettings > > > accelerationsOfAsterix;

    // accelerationsOfAsterix[ "Earth" ].push_back( boost::make_shared< AccelerationSettings >( basic_astrodynamics::central_gravity ) );

    accelerationsOfAsterix[ "Earth" ].push_back( boost::make_shared< SphericalHarmonicAccelerationSettings >( 3, 0 ) );

    accelerationsOfAsterix[ "Sun" ].push_back( boost::make_shared< AccelerationSettings >( basic_astrodynamics::central_gravity ) );

    accelerationsOfAsterix[ "Moon" ].push_back( boost::make_shared< AccelerationSettings >( basic_astrodynamics::central_gravity ) );

    accelerationsOfAsterix[ "Earth" ].push_back( boost::make_shared< AccelerationSettings >( basic_astrodynamics::aerodynamic ) );

    accelerationsOfAsterix[ "Sun" ].push_back( boost::make_shared< AccelerationSettings >( basic_astrodynamics::cannon_ball_radiation_pressure ) );


    accelerationMap[  "Asterix" ] = accelerationsOfAsterix;
    bodiesToPropagate.push_back( "Asterix" );
    centralBodies.push_back( "Earth" );

    basic_astrodynamics::AccelerationMap accelerationModelMap = createAccelerationModelsMap(
                bodyMap, accelerationMap, bodiesToPropagate, centralBodies );

    //////////////////////////////////////////////////////////////////////////////////////////////////////////////////
    //////////////////             CREATE PROPAGATION SETTINGS            ////////////////////////////////////////////
    //////////////////////////////////////////////////////////////////////////////////////////////////////////////////

    // Set Keplerian elements for Asterix.
    Eigen::Vector6d asterixInitialStateInKeplerianElements;
    asterixInitialStateInKeplerianElements( semiMajorAxisIndex ) = 24361e3;
    asterixInitialStateInKeplerianElements( eccentricityIndex ) = 0.73027;
    asterixInitialStateInKeplerianElements( inclinationIndex ) = 0.1744;
    asterixInitialStateInKeplerianElements( argumentOfPeriapsisIndex ) = 0.0;
    asterixInitialStateInKeplerianElements( longitudeOfAscendingNodeIndex ) = 5.388;
    asterixInitialStateInKeplerianElements( trueAnomalyIndex ) = 1.571;

    double earthGravitationalParameter = bodyMap.at( "Earth" )->getGravityFieldModel( )->getGravitationalParameter( );
    const Eigen::Vector6d asterixInitialState = convertKeplerianToCartesianElements(
                asterixInitialStateInKeplerianElements, earthGravitationalParameter );


    TranslationalPropagatorType ptype = dsst;

    boost::shared_ptr< TranslationalStatePropagatorSettings< double > > propagatorSettings =
            boost::make_shared< TranslationalStatePropagatorSettings< double > >
            ( centralBodies, accelerationModelMap, bodiesToPropagate, asterixInitialState, simulationEndEpoch, ptype );

    const double fixedStepSize = ptype == dsst ? physical_constants::JULIAN_DAY : 30.0;
    boost::shared_ptr< IntegratorSettings< > > integratorSettings =
            boost::make_shared< IntegratorSettings< > >
            ( rungeKutta4, simulationStartEpoch, fixedStepSize );

    //////////////////////////////////////////////////////////////////////////////////////////////////////////////////
    //////////////////             PROPAGATE ORBIT            ////////////////////////////////////////////////////////
    //////////////////////////////////////////////////////////////////////////////////////////////////////////////////

    // Create simulation object and propagate dynamics.
    std::cout << "Took "
              << duration_cast< milliseconds >( steady_clock::now() - t ).count() * 1.0e-3
              << " s to set up the propagation."  << std::endl;
    t = steady_clock::now();

    SingleArcDynamicsSimulator< > dynamicsSimulator(
                bodyMap, integratorSettings, propagatorSettings, true, false, false );
    std::map< double, Eigen::VectorXd > integrationResult = dynamicsSimulator.getEquationsOfMotionNumericalSolution( );

    std::cout << "Took "
              << duration_cast< milliseconds >( steady_clock::now() - t ).count() * 1.0e-3
              << " s to perform the propagation, which was "
              << ( (--integrationResult.end())->first - simulationStartEpoch ) / physical_constants::JULIAN_DAY
              << " days long." << std::endl;


    //////////////////////////////////////////////////////////////////////////////////////////////////////////////////
    //////////////////        PROVIDE OUTPUT TO CONSOLE AND FILES           //////////////////////////////////////////
    //////////////////////////////////////////////////////////////////////////////////////////////////////////////////

    Eigen::Vector6d finalState = (--integrationResult.end( ) )->second;
    Eigen::Vector6d keplerianState = convertCartesianToKeplerianElements( finalState, earthGravitationalParameter );
    std::cout << "Final Keplerian state:\n" << keplerianState << std::endl;

    /*
    Eigen::Vector6d cowellKeplerianState;
    cowellKeplerianState << 2.42725e+07,
                            0.729089,
                            0.174239,
                            0.166021,
                            3.29148,
                            2.82826;
    std::cout << "Final Keplerian state using Cowell propagator:\n" << cowellKeplerianState << std::endl;
    */

    // FIXME: this unit test is not comparing anything... will always succeed as long as no exception is thrown.

}

BOOST_AUTO_TEST_SUITE_END( )

} // namespace unit_tests

} // namespace tudat







