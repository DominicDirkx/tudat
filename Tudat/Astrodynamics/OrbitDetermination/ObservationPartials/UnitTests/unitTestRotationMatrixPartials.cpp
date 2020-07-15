/*    Copyright (c) 2010-2019, Delft University of Technology
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
#include <boost/test/unit_test.hpp>
#include <boost/make_shared.hpp>
#include <boost/lambda/lambda.hpp>

#include "Tudat/External/SpiceInterface/spiceInterface.h"
#include "Tudat/SimulationSetup/EstimationSetup/createCartesianStatePartials.h"
#include "Tudat/Astrodynamics/OrbitDetermination/ObservationPartials/rotationMatrixPartial.h"
#include "Tudat/Astrodynamics/ReferenceFrames/referenceFrameTransformations.h"
#include "Tudat/InputOutput/basicInputOutput.h"

namespace tudat
{
namespace unit_tests
{

using namespace ephemerides;
using namespace estimatable_parameters;
using namespace observation_partials;

BOOST_AUTO_TEST_SUITE( test_rotation_matrix_partaisl )

//! Test whether partial derivatives of rotation matrix computed by SimpleRotationalEphemeris works correctly
BOOST_AUTO_TEST_CASE( testSimpleRotationalEphemerisPartials )
{
    // Load spice kernels.
    spice_interface::loadStandardSpiceKernels( );

    // Create rotation model
    double nominalRotationRate = 2.0 * mathematical_constants::PI / 86400.0;
    std::shared_ptr< SimpleRotationalEphemeris > rotationalEphemeris =
            std::make_shared< SimpleRotationalEphemeris >(
                spice_interface::computeRotationQuaternionBetweenFrames( "ECLIPJ2000", "IAU_Earth", 1.0E7 ),
                nominalRotationRate, 1.0E7, "ECLIPJ2000", "IAU_Earth" );

    {
        // Create partial object.
        std::shared_ptr< RotationMatrixPartialWrtConstantRotationRate > rotationMatrixPartialObject =
                std::make_shared< RotationMatrixPartialWrtConstantRotationRate >( rotationalEphemeris );

        // Compute partial analytically
        double testTime = 1.0E6;
        Eigen::Matrix3d rotationMatrixPartial =
                rotationMatrixPartialObject->calculatePartialOfRotationMatrixToBaseFrameWrParameter(
                    testTime ).at( 0 );

        Eigen::Matrix3d rotationMatrixDerivativePartial =
                rotationMatrixPartialObject->calculatePartialOfRotationMatrixDerivativeToBaseFrameWrParameter(
                    testTime ).at( 0 );

        // Compute partial numerically.
        double perturbation = 1.0E-12;
        rotationalEphemeris->resetRotationRate( nominalRotationRate + perturbation );
        Eigen::Matrix3d upperturbedRotationMatrix = rotationalEphemeris->getRotationToBaseFrame(
                    testTime).toRotationMatrix( );
        Eigen::Matrix3d upperturbedRotationMatrixDerivative = rotationalEphemeris->getDerivativeOfRotationToBaseFrame(
                    testTime );

        rotationalEphemeris->resetRotationRate( nominalRotationRate - perturbation );
        Eigen::Matrix3d downperturbedRotationMatrix = rotationalEphemeris->getRotationToBaseFrame(
                    testTime).toRotationMatrix( );
        Eigen::Matrix3d downperturbedRotationMatrixDerivative = rotationalEphemeris->getDerivativeOfRotationToBaseFrame(
                    testTime );

        Eigen::Matrix3d numericalRotationMatrixPartial =
                ( upperturbedRotationMatrix - downperturbedRotationMatrix ) / ( 2.0 * perturbation );
        Eigen::Matrix3d numericalRotationMatrixDerivativePartial =
                ( upperturbedRotationMatrixDerivative - downperturbedRotationMatrixDerivative ) / ( 2.0 * perturbation );

        Eigen::Matrix3d matrixDifference = rotationMatrixPartial - numericalRotationMatrixPartial;

        // Compare analytical and numerical result.
        for( unsigned int i = 0; i < 3; i++ )
        {
            for( unsigned int j = 0; j < 3; j++ )
            {
                BOOST_CHECK_SMALL( std::fabs( matrixDifference( i, j ) ), 0.1 );
            }
        }

        matrixDifference = rotationMatrixDerivativePartial - numericalRotationMatrixDerivativePartial;

        // Compare analytical and numerical result.
        for( unsigned int i = 0; i < 3; i++ )
        {
            for( unsigned int j = 0; j < 3; j++ )
            {
                BOOST_CHECK_SMALL( std::fabs( matrixDifference( i, j ) ), 1.0E-5 );
            }
        }
    }

    {
        // Create partial object.
        std::shared_ptr< RotationMatrixPartialWrtPoleOrientation > rotationMatrixPartialObject =
                std::make_shared< RotationMatrixPartialWrtPoleOrientation >( rotationalEphemeris );

        // Compute partial analytically
        double testTime = 1.0E6;
        std::vector< Eigen::Matrix3d > rotationMatrixPartials =
                rotationMatrixPartialObject->calculatePartialOfRotationMatrixToBaseFrameWrParameter(
                    testTime );

        std::vector< Eigen::Matrix3d > rotationMatrixDerivativePartials =
                rotationMatrixPartialObject->calculatePartialOfRotationMatrixDerivativeToBaseFrameWrParameter(
                    testTime );

        Eigen::Vector3d nominalEulerAngles = rotationalEphemeris->getInitialEulerAngles( );
        double perturbedAngle;

        // Compute partial numerically.
        double perturbation = 1.0E-6;
        {


            // Compute partial for right ascension numerically.
            {
                perturbedAngle = nominalEulerAngles( 0 ) + perturbation;
                rotationalEphemeris->resetInitialPoleRightAscensionAndDeclination( perturbedAngle, nominalEulerAngles( 1 ) );
                Eigen::Matrix3d upperturbedRotationMatrix = rotationalEphemeris->getRotationToBaseFrame(
                            testTime).toRotationMatrix( );
                Eigen::Matrix3d upperturbedRotationMatrixDerivative = rotationalEphemeris->getDerivativeOfRotationToBaseFrame(
                            testTime );

                perturbedAngle = nominalEulerAngles( 0 ) - perturbation;
                rotationalEphemeris->resetInitialPoleRightAscensionAndDeclination( perturbedAngle, nominalEulerAngles( 1 ) );
                Eigen::Matrix3d downperturbedRotationMatrix = rotationalEphemeris->getRotationToBaseFrame(
                            testTime).toRotationMatrix( );
                Eigen::Matrix3d downperturbedRotationMatrixDerivative =
                        rotationalEphemeris->getDerivativeOfRotationToBaseFrame(
                            testTime );

                Eigen::Matrix3d numericalRotationMatrixPartial =
                        ( upperturbedRotationMatrix - downperturbedRotationMatrix ) / ( 2.0 * perturbation );
                Eigen::Matrix3d numericalRotationMatrixDerivativePartial =
                        ( upperturbedRotationMatrixDerivative - downperturbedRotationMatrixDerivative ) /
                        ( 2.0 * perturbation );

                Eigen::Matrix3d matrixDifference = rotationMatrixPartials.at( 0 ) - numericalRotationMatrixPartial;

                // Compare analytical and numerical result.
                for( unsigned int i = 0; i < 3; i++ )
                {
                    for( unsigned int j = 0; j < 3; j++ )
                    {
                        BOOST_CHECK_SMALL( std::fabs( matrixDifference( i, j ) ), 1.0E-8 );
                    }
                }

                matrixDifference = rotationMatrixDerivativePartials.at( 0 ) - numericalRotationMatrixDerivativePartial;

                // Compare analytical and numerical result.
                for( unsigned int i = 0; i < 3; i++ )
                {
                    for( unsigned int j = 0; j < 3; j++ )
                    {
                        BOOST_CHECK_SMALL( std::fabs( matrixDifference( i, j ) ), 1.0E-13 );
                    }
                }
            }

            // Compute partial for declination numerically.
            {
                perturbedAngle = nominalEulerAngles( 1 ) + perturbation;
                rotationalEphemeris->resetInitialPoleRightAscensionAndDeclination( nominalEulerAngles( 0 ), perturbedAngle );
                Eigen::Matrix3d upperturbedRotationMatrix = rotationalEphemeris->getRotationToBaseFrame(
                            testTime).toRotationMatrix( );
                Eigen::Matrix3d upperturbedRotationMatrixDerivative = rotationalEphemeris->getDerivativeOfRotationToBaseFrame(
                            testTime );

                perturbedAngle = nominalEulerAngles( 1 ) - perturbation;
                rotationalEphemeris->resetInitialPoleRightAscensionAndDeclination( nominalEulerAngles( 0 ), perturbedAngle );
                Eigen::Matrix3d downperturbedRotationMatrix = rotationalEphemeris->getRotationToBaseFrame(
                            testTime).toRotationMatrix( );
                Eigen::Matrix3d downperturbedRotationMatrixDerivative =
                        rotationalEphemeris->getDerivativeOfRotationToBaseFrame( testTime );

                Eigen::Matrix3d numericalRotationMatrixPartial =
                        ( upperturbedRotationMatrix - downperturbedRotationMatrix ) / ( 2.0 * perturbation );
                Eigen::Matrix3d numericalRotationMatrixDerivativePartial =
                        ( upperturbedRotationMatrixDerivative -
                          downperturbedRotationMatrixDerivative ) / ( 2.0 * perturbation );

                Eigen::Matrix3d matrixDifference = rotationMatrixPartials.at( 1 ) - numericalRotationMatrixPartial;

                // Compare analytical and numerical result.
                for( unsigned int i = 0; i < 3; i++ )
                {
                    for( unsigned int j = 0; j < 3; j++ )
                    {
                        BOOST_CHECK_SMALL( std::fabs( matrixDifference( i, j ) ), 1.0E-8 );
                    }
                }

                matrixDifference = rotationMatrixDerivativePartials.at( 1 ) - numericalRotationMatrixDerivativePartial;

                // Compare analytical and numerical result.
                for( unsigned int i = 0; i < 3; i++ )
                {
                    for( unsigned int j = 0; j < 3; j++ )
                    {
                        BOOST_CHECK_SMALL( std::fabs( matrixDifference( i, j ) ), 1.0E-13 );
                    }
                }
            }
        }
    }
}


//BOOST_AUTO_TEST_CASE( testSynchronousRotationPartialsTest )
//{
//    // Define nominal state
//    Eigen::Vector6d nominalState =
//            spice_interface::getBodyCartesianStateAtEpoch(
//                "Moon", "Earth", "ECLIPJ2000", "None", 1.0E7 );

//    Eigen::MatrixXd analyticalPartial = getTestPartial( nominalState );

//    double positionPerturbation = 100.0;
//    double velocityPerturbation = 0.1;
//    Eigen::VectorXd upPerturbedTestValue, downPerturbedTestValue;
//    Eigen::Vector6d currentState;

//    int partialSize = 1;
//    Eigen::Matrix< double, Eigen::Dynamic, 6 > numericalPartial =
//            Eigen::MatrixXd::Zero( partialSize, 6 );

//    for( int i = 0; i < 3; i++ )
//    {
//        currentState = nominalState;
//        currentState( i ) += positionPerturbation;
//        upPerturbedTestValue = getTestFunction( currentState );

//        currentState = nominalState;
//        currentState( i ) -= positionPerturbation;
//        downPerturbedTestValue = getTestFunction( currentState );

//        numericalPartial.block( 0, i, partialSize, 1 ) = ( upPerturbedTestValue - downPerturbedTestValue )
//                / ( 2.0 * positionPerturbation );


//        currentState = nominalState;
//        currentState( i + 3 ) += velocityPerturbation;
//        upPerturbedTestValue = getTestFunction( currentState );

//        currentState = nominalState;
//        currentState( i + 3 ) -= velocityPerturbation;
//        downPerturbedTestValue = getTestFunction( currentState );

//        numericalPartial.block( 0, i + 3, partialSize, 1 ) = ( upPerturbedTestValue - downPerturbedTestValue )
//                / ( 2.0 * velocityPerturbation );
//    }
//    std::cout<<"ANALYTICAL: "<<std::endl<<analyticalPartial<<std::endl;
//    std::cout<<"REL. DIFF.: "<<std::endl<<( numericalPartial - analyticalPartial ).cwiseQuotient(
//                   analyticalPartial )<<std::endl;

//}

//! Test whether partial derivatives of rotation matrix computed by SynchronousRotationalEphemeris works correctly
BOOST_AUTO_TEST_CASE( testSynchronousRotationPartials )
{
    // Define nominal state
    Eigen::Vector6d nominalState =
            spice_interface::getBodyCartesianStateAtEpoch(
                "Moon", "Earth", "ECLIPJ2000", "None", 1.0E7 + 6.0E5 );

    double scaledLibationAmplitude = -3.0;
    std::shared_ptr< DirectLongitudeLibrationCalculator > librationCalculator =
            std::make_shared< DirectLongitudeLibrationCalculator >( scaledLibationAmplitude );
    Eigen::MatrixXd librationAnglePartials =
            calculatePartialOfDirectLibrationAngleWrtCartesianStates( nominalState, scaledLibationAmplitude );

    // Define nominal state function
    Eigen::Vector6d currentState = nominalState;
    std::function< Eigen::Vector6d( const double, bool ) > relativeStateFunction =
            [ & ]( const double, bool ){ return currentState; };

    // Create rotation model
    std::shared_ptr< ephemerides::SynchronousRotationalEphemeris > synchronousRotationModel =
            std::make_shared< ephemerides::SynchronousRotationalEphemeris >(
                relativeStateFunction, "SSB", "Mercury_Fixed", "ECLIPJ2000" );

    // Create rotation partial model
    std::shared_ptr< RotationMatrixPartial > rotationMatrixPartialObject =
            std::make_shared< SynchronousRotationMatrixPartialWrtTranslationalState >( synchronousRotationModel );

    // Define test settings
    double testTime = 1.0E7;

    double positionPerturbation = 1000.0;
    double velocityPerturbation = 0.01;

    // Test partials w.r.t. position and velocity components
    std::vector< Eigen::Matrix3d > rotationMatrixPartials =
            rotationMatrixPartialObject->calculatePartialOfRotationMatrixToBaseFrameWrParameter( testTime );

    Eigen::Matrix< double, 1, 6 > numericalLibrationAnglePartials =
            Eigen::Matrix< double, 1, 6 >::Zero( );
    Eigen::Matrix3d upPerturbedRotationMatrix, downPerturbedRotationMatrix;
    for( int i = 0; i < 3; i++ )
    {
        currentState = nominalState;
        currentState( i ) += positionPerturbation;
        upPerturbedRotationMatrix = synchronousRotationModel->getRotationToBaseFrame(
                    1.0E7 ).toRotationMatrix( );
        double upPerturbedLibrationAngle = librationCalculator->getLibrationAngleWrtFullySynchronousRotation(
                    currentState, testTime );

        currentState = nominalState;
        currentState( i ) -= positionPerturbation;
        downPerturbedRotationMatrix = synchronousRotationModel->getRotationToBaseFrame(
                    1.0E7 ).toRotationMatrix( );
        double downPerturbedLibrationAngle = librationCalculator->getLibrationAngleWrtFullySynchronousRotation(
                    currentState, testTime );
        numericalLibrationAnglePartials( i ) = ( upPerturbedLibrationAngle - downPerturbedLibrationAngle )
                / ( 2.0 * positionPerturbation );

        Eigen::Matrix3d relativePartialError =
                ( ( upPerturbedRotationMatrix - downPerturbedRotationMatrix ) / ( 2.0 * positionPerturbation ) -
                  rotationMatrixPartials.at( i ) ) / rotationMatrixPartials.at( i ).norm( );

        for( int j = 0; j < 3; j++ )
        {
            for( int k = 0; k < 3; k++ )
            {
                BOOST_CHECK_SMALL( std::fabs( relativePartialError( j, k ) ), 1.0E-9 );
            }
        }

        currentState = nominalState;
        currentState( i + 3 ) += velocityPerturbation;
        upPerturbedRotationMatrix = synchronousRotationModel->getRotationToBaseFrame(
                    1.0E7 ).toRotationMatrix( );
        upPerturbedLibrationAngle = librationCalculator->getLibrationAngleWrtFullySynchronousRotation(
                    currentState, testTime );

        currentState = nominalState;
        currentState( i + 3 ) -= velocityPerturbation;
        downPerturbedRotationMatrix = synchronousRotationModel->getRotationToBaseFrame(
                    1.0E7 ).toRotationMatrix( );
        downPerturbedLibrationAngle = librationCalculator->getLibrationAngleWrtFullySynchronousRotation(
                    currentState, testTime );

        relativePartialError =
                ( ( upPerturbedRotationMatrix - downPerturbedRotationMatrix ) / ( 2.0 * velocityPerturbation ) -
                  rotationMatrixPartials.at( i + 3 ) ) / rotationMatrixPartials.at( i + 3 ).norm( );

        numericalLibrationAnglePartials( i + 3 ) = ( upPerturbedLibrationAngle - downPerturbedLibrationAngle )
                / ( 2.0 * velocityPerturbation );
        for( int j = 0; j < 3; j++ )
        {
            for( int k = 0; k < 3; k++ )
            {
                BOOST_CHECK_SMALL( std::fabs( relativePartialError( j, k ) ), 1.0E-9 );
            }
        }
    }

    {
        std::shared_ptr< DirectLongitudeLibrationCalculator > librationModel =
                std::make_shared< DirectLongitudeLibrationCalculator >( 0.0 );
        synchronousRotationModel->setLibrationCalculation( librationModel );

        std::shared_ptr< RotationMatrixPartial > rotationMatrixPartialObjectWrtLibration =
                std::make_shared< RotationMatrixPartialWrtScaledLongitudeLibrationAmplitude >( synchronousRotationModel );

        double amplitudePerturbation = 0.001;

        librationModel->setScaledLibrationAmplitude( amplitudePerturbation );
        upPerturbedRotationMatrix = synchronousRotationModel->getRotationToBaseFrame(
                    1.0E7 ).toRotationMatrix( );

        librationModel->setScaledLibrationAmplitude( -amplitudePerturbation );
        downPerturbedRotationMatrix = synchronousRotationModel->getRotationToBaseFrame(
                    1.0E7 ).toRotationMatrix( );

        Eigen::Matrix3d numericalRotationMatrixWrtPartial = ( upPerturbedRotationMatrix - downPerturbedRotationMatrix ) /
                ( 2.0 * amplitudePerturbation );

        Eigen::Matrix3d analyticalRotationMatrixWrtPartial =
                rotationMatrixPartialObjectWrtLibration->calculatePartialOfRotationMatrixToBaseFrameWrParameter( testTime ).at( 0 );

        std::cout<<numericalRotationMatrixWrtPartial<<std::endl<<std::endl;
        std::cout<<analyticalRotationMatrixWrtPartial<<std::endl<<std::endl;

        std::cout<<( analyticalRotationMatrixWrtPartial - numericalRotationMatrixWrtPartial ).cwiseQuotient(
                       numericalRotationMatrixWrtPartial )<<std::endl;
    }




}

BOOST_AUTO_TEST_SUITE_END( )

} // namespace unit_tests

} // namespace tudat





