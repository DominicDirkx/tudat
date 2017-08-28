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

#include <boost/test/floating_point_comparison.hpp>
#include <boost/test/unit_test.hpp>

#include <Eigen/Core>

#include "Tudat/Astrodynamics/BasicAstrodynamics/unitConversions.h"
#include "Tudat/Basics/testMacros.h"

#include "Tudat/Astrodynamics/BasicAstrodynamics/timeConversions.h"
#include "Tudat/Astrodynamics/Ephemerides/simpleRotationalEphemeris.h"
//#include "Tudat/External/SpiceInterface/spiceInterface.h"
#include "Tudat/InputOutput/basicInputOutput.h"
#include "Tudat/InputOutput/matrixTextFileReader.h"
#include "Tudat/Astrodynamics/ReferenceFrames/referenceFrameTransformations.h"
#include "Tudat/SimulationSetup/EnvironmentSetup/createRotationModel.h"

namespace tudat
{
namespace unit_tests
{

BOOST_AUTO_TEST_SUITE( test_simple_rotational_ephemeris )

// Test simple rotational ephemeris class by compariing results to Venus-fixed frame with Spice.
BOOST_AUTO_TEST_CASE( testSimpleRotationalEphemeris )
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

    boost::shared_ptr< simulation_setup::RotationModelSettings > rotationModelSettings =
            boost::make_shared< simulation_setup::DirectRotationVariationSettings >(
                rightAscensionPolynomialTerms, declinationPolynomialTerms, primeMeridianPolynomialTerms,
                rightAscensionLibrations, declinationLibrations, primeMeridianLibrations,
                fromIntermediateFrameToBaseFrame, originalFrame, targetFrame  );

    boost::shared_ptr< ephemerides::RotationalEphemeris > marsRotationModel =
            simulation_setup::createRotationModel( rotationModelSettings, "Mars" );

    double currentRightAscension, currentDeclination, currentPrimeMeridian;
    Eigen::Quaterniond currentRotation, currentExpectedRotation;
    Eigen::Matrix3d currentRotationDerivative;

    currentRightAscension = rightAscensionPolynomialTerms[ 0 ];
    currentDeclination = declinationPolynomialTerms[ 0 ];
    currentPrimeMeridian = primeMeridianPolynomialTerms[ 0 ];

    for( unsigned int i = 0; i < periodicRightAscensionVariations.rows( ); i++ )
    {
        currentRightAscension += periodicRightAscensionVariations( i, 1 );
    }

    for( unsigned int i = 0; i < periodicDeclinationVariations.rows( ); i++ )
    {
        currentDeclination += periodicDeclinationVariations( i, 1 );

    }

    for( unsigned int i = 0; i < periodicPrimeMeridianVariations.rows( ); i++ )
    {
        currentPrimeMeridian += periodicPrimeMeridianVariations( i, 1 );

    }

    currentExpectedRotation = reference_frames::getInertialToPlanetocentricFrameTransformationQuaternion(
                currentDeclination, currentRightAscension, currentPrimeMeridian );
    currentRotation = marsRotationModel->getRotationToTargetFrame( 0.0 );
    currentRotationDerivative = marsRotationModel->getDerivativeOfRotationToTargetFrame( 0.0 );

    TUDAT_CHECK_MATRIX_CLOSE_FRACTION(
                ( currentExpectedRotation.toRotationMatrix( ) ),
                ( currentRotation.toRotationMatrix( ) ), std::numeric_limits< double>::epsilon( ) );

    currentRotation = marsRotationModel->getRotationToBaseFrame( 0.0 );
    TUDAT_CHECK_MATRIX_CLOSE_FRACTION(
                ( currentExpectedRotation.toRotationMatrix( ).transpose( ) ),
                ( currentRotation.toRotationMatrix( ) ), std::numeric_limits< double>::epsilon( ) );

    double timePerturbation = 1.0;
    Eigen::Matrix3d upperturbedRotationMatrix = marsRotationModel->getRotationToTargetFrame( timePerturbation ).toRotationMatrix( );
    Eigen::Matrix3d downperturbedRotationMatrix = marsRotationModel->getRotationToTargetFrame( -timePerturbation ).toRotationMatrix( );

    Eigen::Matrix3d currentExpectedRotationDerivative = ( upperturbedRotationMatrix - downperturbedRotationMatrix ) /
            ( 2.0 * timePerturbation );

    double testTime = 1.0E8;

    {
        double currentArgument = TUDAT_NAN;
        currentRightAscension = rightAscensionPolynomialTerms[ 0 ] + rightAscensionPolynomialTerms[ 1 ] * testTime;
        currentDeclination = declinationPolynomialTerms[ 0 ] + declinationPolynomialTerms[ 1 ] * testTime;
        currentPrimeMeridian = primeMeridianPolynomialTerms[ 0 ] + primeMeridianPolynomialTerms[ 1 ] * testTime;

        for( unsigned int i = 0; i < periodicRightAscensionVariations.rows( ); i++ )
        {
            currentArgument = 2.0 * mathematical_constants::PI * testTime / periodicRightAscensionVariations( i, 2 );
            currentRightAscension +=
                    periodicRightAscensionVariations( i, 1 ) * std::cos( currentArgument ) +
                    periodicRightAscensionVariations( i, 0 ) * std::sin( currentArgument );
        }

        for( unsigned int i = 0; i < periodicDeclinationVariations.rows( ); i++ )
        {
            currentArgument = 2.0 * mathematical_constants::PI * testTime / periodicDeclinationVariations( i, 2 );
            currentDeclination +=
                    periodicDeclinationVariations( i, 1 ) * std::cos( currentArgument ) +
                    periodicDeclinationVariations( i, 0 ) * std::sin( currentArgument );

        }

        for( unsigned int i = 0; i < periodicPrimeMeridianVariations.rows( ); i++ )
        {
            currentArgument = 2.0 * mathematical_constants::PI * testTime / periodicPrimeMeridianVariations( i, 2 );
            currentPrimeMeridian +=
                    periodicPrimeMeridianVariations( i, 1 ) * std::cos( currentArgument ) +
                    periodicPrimeMeridianVariations( i, 0 ) * std::sin( currentArgument );

        }

        currentExpectedRotation = reference_frames::getInertialToPlanetocentricFrameTransformationQuaternion(
                    currentDeclination, currentRightAscension, currentPrimeMeridian );
        currentRotation = marsRotationModel->getRotationToTargetFrame( testTime );

        for( unsigned int i = 0; i < 3; i++ )
        {
            for( unsigned int j = 0; j < 3; j++ )
            {
                BOOST_CHECK_SMALL( std::fabs( currentExpectedRotation.toRotationMatrix( )( i, j ) -
                                   currentRotation.toRotationMatrix( )( i, j ) ), 4.0 * std::numeric_limits< double>::epsilon( ) );
            }
        }

        currentRotation = marsRotationModel->getRotationToBaseFrame( testTime );

        for( unsigned int i = 0; i < 3; i++ )
        {
            for( unsigned int j = 0; j < 3; j++ )
            {
                BOOST_CHECK_SMALL( std::fabs( ( currentExpectedRotation.inverse( ) ).toRotationMatrix( )( i, j ) -
                                   currentRotation.toRotationMatrix( )( i, j ) ), 4.0 * std::numeric_limits< double>::epsilon( ) );
            }
        }
    }


    std::cout<<currentExpectedRotationDerivative<<std::endl<<std::endl<<
            currentRotationDerivative<<std::endl;

    throw std::runtime_error( "Error, directly perturbed rotation matrix derivative not fully implemented" );





}

BOOST_AUTO_TEST_SUITE_END( )

} // namespace unit_tests
} // namespace tudat
