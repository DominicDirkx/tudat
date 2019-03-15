/*    Copyright (c) 2010-2018, Delft University of Technology
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
#include <boost/test/floating_point_comparison.hpp>
#include <boost/test/unit_test.hpp>

#include <Eigen/Core>

#include "Tudat/Basics/testMacros.h"

#include "Tudat/Astrodynamics/Ephemerides/tidallyLockedRotationalEphemeris.h"
#include "Tudat/SimulationSetup/tudatSimulationHeader.h"

namespace tudat
{
namespace unit_tests
{

using namespace ephemerides;
using namespace simulation_setup;
using namespace input_output;

BOOST_AUTO_TEST_SUITE( test_locked_rotational_ephemeris )

BOOST_AUTO_TEST_CASE( test_TidallyLockedRotationModel )
{

    tudat::spice_interface::loadStandardSpiceKernels( );
    tudat::spice_interface::loadSpiceKernelInTudat( input_output::getSpiceKernelPath( ) + "jup310_small.bsp" );
    std::vector< std::string > bodyNames;

    bodyNames.push_back( "Io" );
    bodyNames.push_back( "Europa" );
    bodyNames.push_back( "Ganymede" );
    bodyNames.push_back( "Sun" );
    bodyNames.push_back( "Jupiter" );

    std::vector< std::string > bodiesToTest;
    bodiesToTest.push_back( "Io" );
    bodiesToTest.push_back( "Europa" );
    bodiesToTest.push_back( "Ganymede" );

    // Specify initial time
    double initialTime = 1.0E7;
    double finalTime = 1.05E7;

    // Create bodies needed in simulation
    std::map< std::string, std::shared_ptr< BodySettings > > bodySettings =
            getDefaultBodySettings( bodyNames, initialTime, finalTime );
    for( int i = 0; i < bodiesToTest.size( ); i++ )
    {
       bodySettings[ bodiesToTest.at( i ) ]->rotationModelSettings = getGalileanMoonLockedRotationModel(
                   bodiesToTest.at( i ), "ECLIPJ2000" );
       bodySettings[ bodiesToTest.at( i ) ]->ephemerisSettings->resetFrameOrigin( "Jupiter" );
    }

    NamedBodyMap bodyMap = createBodies( bodySettings );

    setGlobalFrameBodyEphemerides( bodyMap, "SSB", "ECLIPJ2000" );

    double testTime = 0.0;
    for( unsigned int i = 0; i < bodiesToTest.size( ); i++ )
    {
        std::shared_ptr< TidallyLockedRotationalEphemeris > currentEphemeris =
                std::dynamic_pointer_cast< TidallyLockedRotationalEphemeris >(
                    bodyMap.at( bodiesToTest.at( i ) )->getRotationalEphemeris( ) );

        Eigen::Matrix3d toFixedPrimeMeridianFrame =
                ( reference_frames::getRotatingPlanetocentricToInertialFrameTransformationQuaternion(
                      currentEphemeris->getDeclinationFunction( )( testTime ),
                      currentEphemeris->getRightAscensionFunction( )( testTime ), 0.0 ).inverse( ) ).toRotationMatrix( );
        Eigen::Matrix3d spiceFullRotation = spice_interface::computeRotationQuaternionBetweenFrames(
                    "J2000", "IAU_" + bodiesToTest.at( i ), testTime ).toRotationMatrix( );

        Eigen::Matrix3d inferredPrimeMeridianRotation = spiceFullRotation * toFixedPrimeMeridianFrame.inverse( );
        BOOST_CHECK_SMALL( std::fabs( inferredPrimeMeridianRotation( 0, 2 ) ), 3.0 * std::numeric_limits< double >::epsilon( ) );
        BOOST_CHECK_SMALL( std::fabs( inferredPrimeMeridianRotation( 1, 2 ) ), 3.0 * std::numeric_limits< double >::epsilon( ) );
        BOOST_CHECK_SMALL( std::fabs( inferredPrimeMeridianRotation( 2, 0 ) ), 3.0 * std::numeric_limits< double >::epsilon( ) );
        BOOST_CHECK_SMALL( std::fabs( inferredPrimeMeridianRotation( 2, 1 ) ), 3.0 * std::numeric_limits< double >::epsilon( ) );
        BOOST_CHECK_SMALL( std::fabs( inferredPrimeMeridianRotation( 2, 2 ) - 1 ), 3.0 * std::numeric_limits< double >::epsilon( ) );
        BOOST_CHECK_SMALL( std::fabs( sqrt( inferredPrimeMeridianRotation( 0, 0 ) * inferredPrimeMeridianRotation( 0, 0 ) +
                                            inferredPrimeMeridianRotation( 0, 1 ) * inferredPrimeMeridianRotation( 0, 1 ) ) - 1 ),
                           3.0 * std::numeric_limits< double >::epsilon( ) );
        BOOST_CHECK_SMALL( std::fabs( sqrt( inferredPrimeMeridianRotation( 1, 0 ) * inferredPrimeMeridianRotation( 1, 0 ) +
                                            inferredPrimeMeridianRotation( 1, 1 ) * inferredPrimeMeridianRotation( 1, 1 ) ) - 1 ),
                           3.0 * std::numeric_limits< double >::epsilon( ) );
        BOOST_CHECK_SMALL( std::fabs( inferredPrimeMeridianRotation( 0, 0 ) - inferredPrimeMeridianRotation( 1, 1 ) ), 3.0 * std::numeric_limits< double >::epsilon( ) );
        BOOST_CHECK_SMALL( std::fabs( inferredPrimeMeridianRotation( 0, 1 ) + inferredPrimeMeridianRotation( 1, 0 ) ), 3.0 * std::numeric_limits< double >::epsilon( ) );

        TUDAT_CHECK_MATRIX_CLOSE_FRACTION( ( spice_interface::computeRotationQuaternionBetweenFrames(
                                                 "J2000", "ECLIPJ2000", 0.0 ).toRotationMatrix( ) ),
                                           ( currentEphemeris->getIntermediateToBaseFrameRotation( ).toRotationMatrix( ) ),
                                           ( 3.0 * std::numeric_limits< double >::epsilon( ) ) );
    }


    testTime = 1.04E7;
    setAreBodiesInPropagation( bodyMap, 0 );
    double timePerturbation = 1.0;
    for( unsigned int i = 0; i < bodiesToTest.size( ); i++ )
    {
        Eigen::Quaterniond rotationToLocalFrame = bodyMap.at( bodiesToTest.at( i ) )->getRotationalEphemeris( )->getRotationToTargetFrame( testTime );
        Eigen::Vector6d currentVectorFromMoon = -bodyMap.at( bodiesToTest.at( i ) )->getEphemeris( )->getCartesianState(
                    testTime );
        Eigen::Vector3d moonFixedVectorToPlanet = rotationToLocalFrame * ( -currentVectorFromMoon.segment( 0, 3 ) );
        BOOST_CHECK_SMALL( moonFixedVectorToPlanet( 1 ) / moonFixedVectorToPlanet.norm( ), 3.0 * std::numeric_limits< double >::epsilon( ) );

        Eigen::Matrix3d upPerturbedRotationMatrix =
                ( bodyMap.at( bodiesToTest.at( i ) )->getRotationalEphemeris( )->getRotationToTargetFrame( testTime + timePerturbation ) ).toRotationMatrix( );
        Eigen::Matrix3d downPerturbedRotationMatrix =
                ( bodyMap.at( bodiesToTest.at( i ) )->getRotationalEphemeris( )->getRotationToTargetFrame( testTime - timePerturbation ) ).toRotationMatrix( );
        Eigen::Matrix3d numericalMatrixDerivative  = ( upPerturbedRotationMatrix - downPerturbedRotationMatrix ) / ( 2.0 * timePerturbation );
        Eigen::Matrix3d analyticalMatrixDerivative  = bodyMap.at( bodiesToTest.at( i ) )->getRotationalEphemeris( )->getDerivativeOfRotationToTargetFrame(
                    testTime );

        for( unsigned int i = 0; i < 3; i++ )
        {
            for( unsigned int j = 0; j < 3; j++ )
            {
                BOOST_CHECK_SMALL( ( numericalMatrixDerivative( i, j ) - analyticalMatrixDerivative( i, j ) ) / analyticalMatrixDerivative.norm( ), 1.0E-5 );
            }
        }

        Eigen::Quaterniond directRotationToLocalFrame;
        Eigen::Matrix3d directAnalyticalMatrixDerivative;
        Eigen::Vector3d directAngularVelocityVector;
        bodyMap.at( bodiesToTest.at( i ) )->getRotationalEphemeris( )->getFullRotationalQuantitiesToTargetFrame(
                    directRotationToLocalFrame, directAnalyticalMatrixDerivative, directAngularVelocityVector, testTime );

        BOOST_CHECK_SMALL( std::fabs( rotationToLocalFrame.w( ) - directRotationToLocalFrame.w( ) ), 2.0 * std::numeric_limits< double >::epsilon( ) );
        BOOST_CHECK_SMALL( std::fabs( rotationToLocalFrame.x( ) - directRotationToLocalFrame.x( ) ), std::numeric_limits< double >::epsilon( ) );
        BOOST_CHECK_SMALL( std::fabs( rotationToLocalFrame.y( ) - directRotationToLocalFrame.y( ) ), std::numeric_limits< double >::epsilon( ) );
        BOOST_CHECK_SMALL( std::fabs( rotationToLocalFrame.z( ) - directRotationToLocalFrame.z( ) ), std::numeric_limits< double >::epsilon( ) );

        for( unsigned int i = 0; i < 3; i++ )
        {
            for( unsigned int j = 0; j < 3; j++ )
            {
                BOOST_CHECK_SMALL( ( directAnalyticalMatrixDerivative( i, j ) - analyticalMatrixDerivative( i, j ) ) / analyticalMatrixDerivative.norm( ),
                                   std::numeric_limits< double >::epsilon( ) );
            }
        }
    }


    double testTimeForEphemeris = 1.02E7;
    for( NamedBodyMap::iterator bodyIterator = bodyMap.begin( ); bodyIterator != bodyMap.end( ); bodyIterator++ )
    {
        bodyIterator->second->setStateFromEphemeris( testTimeForEphemeris );
    }

    setAreBodiesInPropagation( bodyMap, 1 );
    for( unsigned int i = 0; i < bodiesToTest.size( ); i++ )
    {
        Eigen::Quaterniond rotationToLocalFrame = bodyMap.at( bodiesToTest.at( i ) )->getRotationalEphemeris( )->getRotationToTargetFrame( testTime );
        Eigen::Vector6d currentVectorFromMoon = -bodyMap.at( bodiesToTest.at( i ) )->getEphemeris( )->getCartesianState(
                    testTimeForEphemeris );
        Eigen::Vector3d moonFixedVectorToPlanet = rotationToLocalFrame * ( -currentVectorFromMoon.segment( 0, 3 ) );
        BOOST_CHECK_SMALL( moonFixedVectorToPlanet( 1 ) / moonFixedVectorToPlanet.norm( ), 5.0E-14 );
    }

}

BOOST_AUTO_TEST_SUITE_END( )

} // namespace unit_tests
} // namespace tudat
