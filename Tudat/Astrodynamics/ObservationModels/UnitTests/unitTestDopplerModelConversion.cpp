/*    Copyright (c) 2010-2018, Delft University of Technology
 *    All rigths reserved
 *
 *    This file is part of the Tudat. Redistribution and use in source and
 *    binary forms, with or without modification, are permitted exclusively
 *    under the terms of the Modified BSD license. You should have received
 *    a copy of the license with this file. If not, please or visit:
 *    http://tudat.tudelft.nl/LICENSE.
 */


#include <boost/test/unit_test.hpp>


#include "Tudat/Astrodynamics/Ephemerides/keplerEphemeris.h"
#include "Tudat/Astrodynamics/ObservationModels/dopplerModelConversions.h"
#include "Tudat/External/SpiceInterface/spiceEphemeris.h"
#include "Tudat/SimulationSetup/tudatEstimationHeader.h"
#include "Tudat/Mathematics/NumericalQuadrature/createNumericalQuadrature.h"


using namespace tudat;
using namespace tudat::observation_models;
using namespace tudat::spice_interface;
using namespace tudat::ephemerides;
using namespace tudat::simulation_setup;
using namespace tudat::orbital_element_conversions;
using namespace tudat::coordinate_conversions;
using namespace tudat::unit_conversions;

namespace tudat
{
namespace unit_tests
{



template< typename ScalarType, typename TimeType >
Eigen::Matrix< ScalarType, 6, 1 > getCurrentTestState(
        const TimeType currentTime )
{
    ScalarType angularVelocity = static_cast< ScalarType >( 2.0 * tudat::mathematical_constants::PI / ( 3600.0 ) );
    ScalarType radius = static_cast< ScalarType >( 3.0E6 );

    Eigen::Matrix< ScalarType, 6, 1 > currentState;
    currentState.setZero( );

    currentState( 0 ) = static_cast< ScalarType >( 1.5E8 ) + static_cast< ScalarType >( 10.0E3 ) * currentTime;
    currentState( 0 ) += radius * std::sin( angularVelocity * static_cast< ScalarType >( currentTime ) );
    currentState( 2 ) += radius * std::cos( angularVelocity * static_cast< ScalarType >( currentTime ) );

    currentState( 3 ) = static_cast< ScalarType >( 10.0E3 );
    currentState( 3 ) += angularVelocity * radius * std::cos( angularVelocity * static_cast< ScalarType >( currentTime ) );
    currentState( 5 ) -= angularVelocity * radius * std::sin( angularVelocity * static_cast< ScalarType >( currentTime ) );

    return currentState;
}

}

}

using namespace tudat::unit_tests;

template< typename ObservationScalarType, typename TimeType >
void testDopplerConversion( )
{

    std::string useDouble;
    if( sizeof( ObservationScalarType ) == 8 )
    {
        useDouble = "1";
    }
    else
    {
        useDouble = "0";
    }

    // Load Spice kernels
    std::string kernelsPath = input_output::getSpiceKernelPath( );

    spice_interface::loadStandardSpiceKernels( );
    spice_interface::loadSpiceKernelInTudat( kernelsPath + "MEX_ROB_130101_131231_001.bsp" );

    // Define bodies to use.
    std::vector< std::string > bodiesToCreate;
    bodiesToCreate.push_back( "Earth" );
    bodiesToCreate.push_back( "Sun" );
    bodiesToCreate.push_back( "Mars" );

    // Specify initial time
    double initialEphemerisTime = 1.0E6;//tudat::basic_astrodynamics::convertCalendarDateToJulianDaysSinceEpoch(
                //2013, 12, 29, 6, 0, 0.0, tudat::basic_astrodynamics::JULIAN_DAY_ON_J2000 ) * 86400.0;
    double finalEphemerisTime = initialEphemerisTime + 86400.0;
    double maximumTimeStep = 3600.0;
    double buffer = 10.0 * maximumTimeStep;

    // Create bodies settings needed in simulation
    std::map< std::string, boost::shared_ptr< BodySettings > > defaultBodySettings =
            getDefaultBodySettings(
                bodiesToCreate, initialEphemerisTime - buffer, finalEphemerisTime + buffer );


    boost::shared_ptr< InterpolatedSpiceEphemerisSettings > marsEphemerisSettings =
            boost::make_shared< InterpolatedSpiceEphemerisSettings >(
                initialEphemerisTime - buffer, finalEphemerisTime + buffer, 300.0 );


    boost::shared_ptr< InterpolatedSpiceEphemerisSettings > earthEphemerisSettings =
            boost::make_shared< InterpolatedSpiceEphemerisSettings >(
                initialEphemerisTime - buffer, finalEphemerisTime + buffer, 300.0 );

    if( sizeof( ObservationScalarType ) != 8 )
    {
        earthEphemerisSettings->setUseLongDoubleStates( true );
        earthEphemerisSettings->setUseExtendedTime( true );
        marsEphemerisSettings->setUseLongDoubleStates( true );
        marsEphemerisSettings->setUseExtendedTime( true );
    }

    defaultBodySettings[ "Mars" ]->ephemerisSettings = //marsEphemerisSettings;
                boost::make_shared< CustomEphemerisSettings< ObservationScalarType, TimeType > >(
                    &getCurrentTestState< ObservationScalarType, TimeType > );
    defaultBodySettings[ "Earth" ]->ephemerisSettings = //earthEphemerisSettings;
                boost::make_shared< ConstantEphemerisSettings >(
                    Eigen::Vector6d::Zero( ) );

    // Create bodies
    NamedBodyMap bodyMap = createBodies( defaultBodySettings );
    bodyMap[ "MEX" ] = boost::make_shared< Body >( );
    bodyMap[ "MEX" ]->setEphemeris( boost::make_shared< tudat::ephemerides::SpiceEphemeris >(
                                        "MEX", "Mars", false, true, true ) );

    // Set Keplerian elements for Asterix.
    Eigen::Vector6d initialStateInKeplerianElements;
    initialStateInKeplerianElements( semiMajorAxisIndex ) = 4000.0E3;
    initialStateInKeplerianElements( eccentricityIndex ) = 0.1;
    initialStateInKeplerianElements( inclinationIndex ) = unit_conversions::convertDegreesToRadians( 85.3 );
    initialStateInKeplerianElements( argumentOfPeriapsisIndex )
            = unit_conversions::convertDegreesToRadians( 235.7 );
    initialStateInKeplerianElements( longitudeOfAscendingNodeIndex )
            = unit_conversions::convertDegreesToRadians( 0 );
    initialStateInKeplerianElements( trueAnomalyIndex ) = unit_conversions::convertDegreesToRadians( 139.87 );

    double marsGravitationalParameter = bodyMap.at( "Mars" )->getGravityFieldModel( )->getGravitationalParameter( );

    bodyMap[ "MarsOrbiter1" ] = boost::make_shared< Body >( );
    bodyMap[ "MarsOrbiter1" ]->setEphemeris(
                boost::make_shared< KeplerEphemeris >(
                    initialStateInKeplerianElements, initialEphemerisTime, marsGravitationalParameter, "Mars" ) );

    Eigen::Vector6d initialStateInKeplerianElements2 = initialStateInKeplerianElements;
    initialStateInKeplerianElements2( longitudeOfAscendingNodeIndex )
            = unit_conversions::convertDegreesToRadians( 30.0 );

    bodyMap[ "MarsOrbiter2" ] = boost::make_shared< Body >( );
    bodyMap[ "MarsOrbiter2" ]->setEphemeris(
                boost::make_shared< KeplerEphemeris >(
                    initialStateInKeplerianElements2, initialEphemerisTime, marsGravitationalParameter, "Mars" ) );

    Eigen::Vector6d initialStateInKeplerianElements3 = initialStateInKeplerianElements;
    initialStateInKeplerianElements3( longitudeOfAscendingNodeIndex )
            = unit_conversions::convertDegreesToRadians( 60.0 );

    bodyMap[ "MarsOrbiter3" ] = boost::make_shared< Body >( );
    bodyMap[ "MarsOrbiter3" ]->setEphemeris(
                boost::make_shared< KeplerEphemeris >(
                    initialStateInKeplerianElements3, initialEphemerisTime, marsGravitationalParameter, "Mars" ) );


    Eigen::Vector6d initialStateInKeplerianElements4 = initialStateInKeplerianElements;
    initialStateInKeplerianElements4( longitudeOfAscendingNodeIndex )
            = unit_conversions::convertDegreesToRadians( 90.0 );

    bodyMap[ "MarsOrbiter4" ] = boost::make_shared< Body >( );
    bodyMap[ "MarsOrbiter4" ]->setEphemeris(
                boost::make_shared< KeplerEphemeris >(
                    initialStateInKeplerianElements4, initialEphemerisTime, marsGravitationalParameter, "Mars" ) );




    // Create ground station
    const Eigen::Vector3d stationCartesianPosition( 1917032.190, 6029782.349, -801376.113 );
    createGroundStation( bodyMap.at( "Earth" ), "Station1", stationCartesianPosition, cartesian_position );

    setGlobalFrameBodyEphemerides( bodyMap, "SSB", "ECLIPJ2000" );

    std::vector< std::string > observedBodies;
    observedBodies.push_back( "Mars" );
    //    observedBodies.push_back( "MarsOrbiter1" );
    //    observedBodies.push_back( "MarsOrbiter2" );
    //    observedBodies.push_back( "MarsOrbiter3" );
    //    observedBodies.push_back( "MarsOrbiter4" );

    for( unsigned int k = 0; k < observedBodies.size( ); k++ )
    {

        // Define link ends for observations.
        LinkEnds linkEnds;
        linkEnds[ transmitter ] = std::make_pair( "Earth" , ""  );
        linkEnds[ receiver ] = std::make_pair( observedBodies.at( k ) , ""  );

        // Create observation settings
        boost::shared_ptr< ObservationSettings > rangeSettings = boost::make_shared< ObservationSettings >
                ( one_way_range );

        boost::shared_ptr< ObservationSettings > openLoopObservableSettings = boost::make_shared< ObservationSettings >
                ( one_way_doppler );

        // Create open-loop bservation model.
        boost::shared_ptr< ObservationModel< 1, ObservationScalarType, TimeType > > rangeObservationModel =
                ObservationModelCreator< 1, ObservationScalarType, TimeType >::createObservationModel(
                    linkEnds, rangeSettings, bodyMap );

        boost::shared_ptr< ObservationModel< 1, ObservationScalarType, TimeType > > openLoopObservationModel =
                ObservationModelCreator< 1, ObservationScalarType, TimeType >::createObservationModel(
                    linkEnds, openLoopObservableSettings, bodyMap );




        boost::shared_ptr< OneWayDopplerObservationModel< ObservationScalarType, TimeType > > dopplerObservationModel =
                boost::dynamic_pointer_cast< OneWayDopplerObservationModel< ObservationScalarType, TimeType > >( openLoopObservableSettings );

        LinkEndType referenceLinkEnd = receiver;

        std::vector< TimeType > integrationTimes = { TimeType( 1.0 ), TimeType( 10.0 ),
                                                     TimeType( 100.0 ), TimeType( 1000.0 ) };
        std::vector< boost::shared_ptr< numerical_quadrature::NumericalQuadratureSettings > > quadratureSettings;
        quadratureSettings.push_back( boost::make_shared< numerical_quadrature::NumericalQuadratureSettings >( numerical_quadrature::trapezoidal_quadrature ) );
        quadratureSettings.push_back( boost::make_shared< numerical_quadrature::NumericalQuadratureSettings >( numerical_quadrature::simpsons_quadrature ) );
        quadratureSettings.push_back( boost::make_shared< numerical_quadrature::GaussianQuadratureSettings >( 3 ) );
        quadratureSettings.push_back( boost::make_shared< numerical_quadrature::GaussianQuadratureSettings >( 5 ) );
        quadratureSettings.push_back( boost::make_shared< numerical_quadrature::GaussianQuadratureSettings >( 7 ) );
        //quadratureSettings.push_back( boost::make_shared< numerical_quadrature::GaussianQuadratureSettings >( 11 ) );

        double simulationTime = 3600.0;
        for( unsigned int i = 0; i < integrationTimes.size( ); i++ )
        {
            boost::shared_ptr< ObservationSettings > closedLoopObservableSettings =
                    boost::make_shared< OneWayDifferencedRangeRateObservationSettings >( boost::lambda::constant( integrationTimes.at( i ) ) );

            // Create closed-loop observation model.
            boost::shared_ptr< ObservationModel< 1, ObservationScalarType, TimeType > > closedLoopObservationModel =
                    ObservationModelCreator< 1, ObservationScalarType, TimeType >::createObservationModel(
                        linkEnds, closedLoopObservableSettings, bodyMap );

            for( unsigned int j = 0; j < quadratureSettings.size( ); j++ )
            {
                std::cout<<i<<" "<<j<<" "<<k<<std::endl;
                std::map< TimeType, ObservationScalarType > syntheticClosedLoopData = convertOpenLoopToClosedLoopDopplerData< ObservationScalarType, TimeType >(
                            boost::dynamic_pointer_cast< OneWayDopplerObservationModel< ObservationScalarType, TimeType > >( openLoopObservationModel ),
                            quadratureSettings.at( j ), referenceLinkEnd, initialEphemerisTime,
                            initialEphemerisTime + simulationTime, integrationTimes.at( i ) );

                std::map< double, ObservationScalarType > simulatedClosedLoopData, syntheticClosedLoopErrors;
                std::map< double, ObservationScalarType > arcStartSimulatedRangeData, arcEndSimulatedRangeData, manualSimulatedClosedLoopData;
                std::map< double, ObservationScalarType > arcStartRealRangeData, arcEndRealRangeData, manualRealClosedLoopData,
                        realArcStartOpenLoopData, realArcEndOpenLoopData;

                for( auto dataIterator = syntheticClosedLoopData.begin( ); dataIterator != syntheticClosedLoopData.end( );
                     dataIterator++ )
                {
                    // Simulate closed loop data, and compare to synthetic value
                    simulatedClosedLoopData[ dataIterator->first ] = closedLoopObservationModel-> computeObservationEntry(
                                dataIterator->first, referenceLinkEnd, 0 );
                    syntheticClosedLoopErrors[ dataIterator->first ] = syntheticClosedLoopData[ dataIterator->first ] +
                            simulatedClosedLoopData[ dataIterator->first ];

                    // Simulate range data
                    arcEndSimulatedRangeData[ dataIterator->first ] = rangeObservationModel->computeObservationEntry(
                                dataIterator->first, referenceLinkEnd, 0 );
                    arcStartSimulatedRangeData[ dataIterator->first ] = rangeObservationModel->computeObservationEntry(
                                dataIterator->first - integrationTimes.at( i ), referenceLinkEnd, 0 );

                    // Compute true range data
                    arcEndRealRangeData[ dataIterator->first ] = tudat::unit_tests::getCurrentTestState< ObservationScalarType, TimeType >(
                                dataIterator->first ).segment( 0, 3 ).norm( );
                    arcStartRealRangeData[ dataIterator->first ] = tudat::unit_tests::getCurrentTestState< ObservationScalarType, TimeType >(
                                dataIterator->first - integrationTimes.at( i ) ).segment( 0, 3 ).norm( );

                    // Compute true open-loop data
                    Eigen::Matrix< ObservationScalarType, 6, 1 > currentState = tudat::unit_tests::getCurrentTestState< ObservationScalarType, TimeType >(
                                dataIterator->first );
                    realArcEndOpenLoopData[ dataIterator->first ] = currentState.segment( 0, 3 ).dot( currentState.segment( 3, 3 ) );
                    Eigen::Matrix< ObservationScalarType, 6, 1 > previousState = tudat::unit_tests::getCurrentTestState< ObservationScalarType, TimeType >(
                                dataIterator->first - integrationTimes.at( i ) );
                    realArcStartOpenLoopData[ dataIterator->first ] = previousState.segment( 0, 3 ).dot( previousState.segment( 3, 3 ) );


                    manualRealClosedLoopData[ dataIterator->first ] =
                            ( arcEndRealRangeData[ dataIterator->first ] - arcStartRealRangeData[ dataIterator->first ] ) /
                            static_cast< ObservationScalarType >( integrationTimes.at( i ) );

                    manualSimulatedClosedLoopData[ dataIterator->first ] =
                            ( arcEndSimulatedRangeData[ dataIterator->first ] - arcStartSimulatedRangeData[ dataIterator->first ] ) /
                            static_cast< ObservationScalarType >( integrationTimes.at( i ) );

                }
                tudat::input_output::writeDataMapToTextFile(
                            syntheticClosedLoopData, "syntheticClosedLoopDopplerErrors_" + std::to_string( i ) + "_" +
                            std::to_string( j ) + "_" + std::to_string( k ) + "_" +  useDouble + ".dat" );

                tudat::input_output::writeDataMapToTextFile(
                            syntheticClosedLoopErrors, "closedLoopDopplerErrors_" + std::to_string( i ) + "_" +
                            std::to_string( j ) + "_" + std::to_string( k ) + "_" +  useDouble + ".dat" );
                tudat::input_output::writeDataMapToTextFile(
                            simulatedClosedLoopData, "closedLoopDoppler_" + std::to_string( i ) + "_" +
                            std::to_string( j ) + "_" + std::to_string( k ) + "_" +  useDouble + ".dat" );

                tudat::input_output::writeDataMapToTextFile(
                            manualSimulatedClosedLoopData, "manualClosedLoopDoppler_" + std::to_string( i ) + "_" +
                            std::to_string( j ) + "_" + std::to_string( k ) + "_" +  useDouble + ".dat" );
                tudat::input_output::writeDataMapToTextFile(
                            manualRealClosedLoopData, "manualRealClosedLoopDoppler_" + std::to_string( i ) + "_" +
                            std::to_string( j ) + "_" + std::to_string( k ) + "_" +  useDouble + ".dat" );
            }

        }
    }


}

int main( )
{
    testDopplerConversion< double, double >( );
    testDopplerConversion< long double, Time >( );
}


