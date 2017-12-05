/*    Copyright (c) 2010-2017, Delft University of Technology
 *    All rigths reserved
 *
 *    This file is part of the Tudat. Redistribution and use in source and
 *    binary forms, with or without modification, are permitted exclusively
 *    under the terms of the Modified BSD license. You should have received
 *    a copy of the license with this file. If not, please or visit:
 *    http://tudat.tudelft.nl/LICENSE.
 *
 *    References
 *      Easy calculation. Newton's Law of Gravity Tutorial,
 *          http://easycalculation.com/physics/classical-physics/learn-newtons-law.php, last
 *          accessed: 12th February, 2012.
 *
 */

#define BOOST_TEST_MAIN

#include <limits>

#include <boost/test/floating_point_comparison.hpp>
#include <boost/test/unit_test.hpp>

#include <Eigen/Core>

#include "Tudat/SimulationSetup/tudatSimulationHeader.h"

namespace tudat
{
namespace unit_tests
{

BOOST_AUTO_TEST_SUITE( test_gravitational_torque )

BOOST_AUTO_TEST_CASE( testDegreeTwoGravitationalTorque )
{
    using namespace tudat::simulation_setup;
    using namespace tudat::basic_astrodynamics;
    using namespace tudat::gravitation;

    // Load Spice kernels.
    spice_interface::loadStandardSpiceKernels( );

    // Set simulation time settings.
    const double simulationStartEpoch = 0.0;
    const double simulationEndEpoch = tudat::physical_constants::JULIAN_DAY;

    // Define body settings for simulation.
    std::vector< std::string > bodiesToCreate;
    bodiesToCreate.push_back( "Earth" );
    bodiesToCreate.push_back( "Moon" );

    for( unsigned int testCase = 0; testCase < 7; testCase++ )
    {

        Eigen::Matrix3d cosineCoefficients = Eigen::Matrix3d::Zero( );
        Eigen::Matrix3d sineCoefficients = Eigen::Matrix3d::Zero( );
        Eigen::Matrix3d inertiaTensorDeviation = Eigen::Matrix3d::Zero( );

        // Create body objects.
        std::map< std::string, boost::shared_ptr< BodySettings > > bodySettings =
                getDefaultBodySettings( bodiesToCreate );

        if( testCase > 0 )
        {
            cosineCoefficients( 0, 0 ) = 1.0;

            if( testCase == 2 )
            {
                cosineCoefficients( 2, 0 ) = 0.01;
            }
            else if( testCase == 3 )
            {
                cosineCoefficients( 2, 2 ) = 0.01;
            }
            else if( testCase == 4 )
            {
                sineCoefficients( 2, 2 ) = 0.01;
            }
            else if( testCase == 5 )
            {
                cosineCoefficients( 2, 1 ) = 0.01;
            }
            else if( testCase == 6 )
            {
                sineCoefficients( 2, 1 ) = 0.01;
            }

            bodySettings[ "Moon" ]->gravityFieldSettings = boost::make_shared< SphericalHarmonicsGravityFieldSettings >(
                        spice_interface::getBodyGravitationalParameter( "Moon" ), spice_interface::getAverageRadius( "Moon" ),
                        cosineCoefficients, sineCoefficients, "IAU_Moon" );
        }


        if( testCase > 0 )
        {
            double c20InertiaContribution  = cosineCoefficients( 2, 0 ) *
                    tudat::basic_mathematics::calculateLegendreGeodesyNormalizationFactor( 2, 0 ) / 3.0;
            inertiaTensorDeviation( 0, 0 ) += c20InertiaContribution;
            inertiaTensorDeviation( 1, 1 ) += c20InertiaContribution;
            inertiaTensorDeviation( 2, 2 ) -= 2.0 * c20InertiaContribution;

            double c21InertiaContribution  = cosineCoefficients( 2, 1 ) *
                    tudat::basic_mathematics::calculateLegendreGeodesyNormalizationFactor( 2, 1 );
            inertiaTensorDeviation( 0, 2 ) += c21InertiaContribution;
            inertiaTensorDeviation( 2, 0 ) += c21InertiaContribution;

            double c22InertiaContribution  = 2.0 * cosineCoefficients( 2, 2 ) *
                    tudat::basic_mathematics::calculateLegendreGeodesyNormalizationFactor( 2, 2 );
            inertiaTensorDeviation( 0, 0 ) -= c22InertiaContribution;
            inertiaTensorDeviation( 1, 1 ) += c22InertiaContribution;

            double s21InertiaContribution  = sineCoefficients( 2, 1 ) *
                    tudat::basic_mathematics::calculateLegendreGeodesyNormalizationFactor( 2, 1 );
            inertiaTensorDeviation( 1, 2 ) += s21InertiaContribution;
            inertiaTensorDeviation( 2, 1 ) += s21InertiaContribution;

            double s22InertiaContribution  = 2.0 * sineCoefficients( 2, 2 ) *
                    tudat::basic_mathematics::calculateLegendreGeodesyNormalizationFactor( 2, 2 );
            inertiaTensorDeviation( 0, 1 ) += s22InertiaContribution;
            inertiaTensorDeviation( 1, 0 ) += s22InertiaContribution;

            inertiaTensorDeviation *= spice_interface::getAverageRadius( "Moon" ) * spice_interface::getAverageRadius( "Moon" ) *
                    spice_interface::getBodyGravitationalParameter( "Moon" ) / physical_constants::GRAVITATIONAL_CONSTANT;
        }

        NamedBodyMap bodyMap = createBodies( bodySettings );
        setGlobalFrameBodyEphemerides( bodyMap, "SSB", "ECLIPJ2000" );

        SelectedTorqueMap selectedTorqueModelMap;
        selectedTorqueModelMap[ "Moon" ][ "Earth" ].push_back(
                    boost::make_shared< TorqueSettings >( second_order_gravitational_torque ) );

        basic_astrodynamics::TorqueModelMap torqueModelMap = createTorqueModelsMap(
                    bodyMap, selectedTorqueModelMap );

        boost::shared_ptr< TorqueModel > secondDegreeGravitationalTorque =
                torqueModelMap.at( "Moon" ).at( "Earth" ).at( 0 );

        double evaluationTime = tudat::physical_constants::JULIAN_DAY / 2.0;

        bodyMap.at( "Moon" )->setStateFromEphemeris( evaluationTime );
        bodyMap.at( "Moon" )->setCurrentRotationalStateToLocalFrameFromEphemeris( evaluationTime );
        bodyMap.at( "Moon" )->setBodyInertiaTensorFromGravityField( 0.4 );

        bodyMap.at( "Earth" )->setStateFromEphemeris( evaluationTime );
        bodyMap.at( "Earth" )->setCurrentRotationalStateToLocalFrameFromEphemeris( evaluationTime );
        bodyMap.at( "Earth" )->setBodyInertiaTensorFromGravityField( 0.4 );

        if( testCase == 0 )
        {
            inertiaTensorDeviation = bodyMap.at( "Moon" )->getBodyInertiaTensor( );
        }

        secondDegreeGravitationalTorque->updateMembers( evaluationTime );

        Eigen::Vector3d earthRelativePosition =
                bodyMap.at( "Moon" )->getCurrentRotationToLocalFrame( ) *
                ( bodyMap.at( "Earth" )->getPosition( ) - bodyMap.at( "Moon" )->getPosition( ) );

        Eigen::Vector3d manualTorque = 3.0 * earthRelativePosition.cross( inertiaTensorDeviation * earthRelativePosition ) *
                bodyMap.at( "Earth" )->getGravityFieldModel( )->getGravitationalParameter( ) /
                std::pow( earthRelativePosition.norm( ), 5.0 );

        Eigen::Vector3d currentTorque = secondDegreeGravitationalTorque->getTorque( );

        Eigen::Vector3d torqueError = ( currentTorque - manualTorque );

        for( unsigned int i = 0; i < 3; i++ )
        {
            BOOST_CHECK_SMALL( std::fabs( torqueError( i ) ), 1.0E-14 * currentTorque.norm( ) );
        }
    }
}

BOOST_AUTO_TEST_CASE( testSphericalGravitationalTorque )
{
    using namespace tudat::simulation_setup;
    using namespace tudat::basic_astrodynamics;
    using namespace tudat::gravitation;

    // Load Spice kernels.
    spice_interface::loadStandardSpiceKernels( );

    // Set simulation time settings.
    const double simulationStartEpoch = 0.0;
    const double simulationEndEpoch = tudat::physical_constants::JULIAN_DAY;

    for( unsigned int testCase = 0; testCase < 3; testCase++ )
    {
        // Define body settings for simulation.
        std::vector< std::string > bodiesToCreate;
        bodiesToCreate.push_back( "Earth" );
        bodiesToCreate.push_back( "Moon" );

        // Create body objects.
        std::map< std::string, boost::shared_ptr< BodySettings > > bodySettings =
                getDefaultBodySettings( bodiesToCreate );
        bodySettings[ "Earth" ]->rotationModelSettings = boost::make_shared< SimpleRotationModelSettings >(
                    "ECLIPJ2000", "IAU_Earth", Eigen::Quaterniond( Eigen::Matrix3d::Identity( ) ), 0.0, 0.0 );
        bodySettings[ "Earth" ]->ephemerisSettings = boost::make_shared< ConstantEphemerisSettings >(
                    Eigen::Vector6d::Zero( ) );

        bodySettings[ "Moon" ]->rotationModelSettings = boost::make_shared< SimpleRotationModelSettings >(
                    "ECLIPJ2000", "IAU_Moon", Eigen::Quaterniond( Eigen::Matrix3d::Identity( ) ), 0.0, 0.0 );
        Eigen::Vector6d moonState = Eigen::Vector6d::Zero( );
        moonState( 0 ) = 1.0E8;
        moonState( 1 ) = 2.0E8;

        bodySettings[ "Moon" ]->ephemerisSettings = boost::make_shared< ConstantEphemerisSettings >(
                    moonState );

        Eigen::Matrix3d cosineCoefficients = Eigen::Matrix3d::Zero( );
        Eigen::Matrix3d sineCoefficients = Eigen::Matrix3d::Zero( );
        //cosineCoefficients( 0, 0 ) = 1.0;
        //cosineCoefficients( 2, 0 ) = 0.1;
        //cosineCoefficients( 2, 1 ) = 0.1;
        //cosineCoefficients( 2, 2 ) = 0.1;

        //sineCoefficients( 2, 1 ) = 0.1;
        if( testCase == 0 )
        {
            cosineCoefficients( 2, 0 ) = 0.1;
        }

        if( testCase == 1 )
        {
            cosineCoefficients( 2, 2 ) = 0.1;
        }

        if( testCase == 2 )
        {
            sineCoefficients( 2, 2 ) = 0.1;
        }

        bodySettings[ "Moon" ]->gravityFieldSettings = boost::make_shared< SphericalHarmonicsGravityFieldSettings >(
                    spice_interface::getBodyGravitationalParameter( "Moon" ), spice_interface::getAverageRadius( "Moon" ),
                    cosineCoefficients, sineCoefficients, "IAU_Moon" );

        NamedBodyMap bodyMap = createBodies( bodySettings );
        setGlobalFrameBodyEphemerides( bodyMap, "SSB", "ECLIPJ2000" );

        SelectedTorqueMap selectedTorqueModelMap;
        selectedTorqueModelMap[ "Moon" ][ "Earth" ].push_back(
                    boost::make_shared< TorqueSettings >( second_order_gravitational_torque ) );
        selectedTorqueModelMap[ "Moon" ][ "Earth" ].push_back(
                    boost::make_shared< SphericalHarmonicTorqueSettings >( 2, 2 ) );

        basic_astrodynamics::TorqueModelMap torqueModelMap = createTorqueModelsMap(
                    bodyMap, selectedTorqueModelMap );

        boost::shared_ptr< TorqueModel > secondDegreeGravitationalTorque =
                torqueModelMap.at( "Moon" ).at( "Earth" ).at( 0 );
        boost::shared_ptr< TorqueModel > sphercialHarmonicGravitationalTorque =
                torqueModelMap.at( "Moon" ).at( "Earth" ).at( 1 );

        double evaluationTime = tudat::physical_constants::JULIAN_DAY / 2.0;

        bodyMap.at( "Moon" )->setStateFromEphemeris( evaluationTime );
        bodyMap.at( "Moon" )->setCurrentRotationalStateToLocalFrameFromEphemeris( evaluationTime );
        bodyMap.at( "Moon" )->setBodyInertiaTensorFromGravityField( 0.0 );

        bodyMap.at( "Earth" )->setStateFromEphemeris( evaluationTime );
        bodyMap.at( "Earth" )->setCurrentRotationalStateToLocalFrameFromEphemeris( evaluationTime );
        bodyMap.at( "Earth" )->setBodyInertiaTensorFromGravityField( 0.0 );

        secondDegreeGravitationalTorque->updateMembers( evaluationTime );
        Eigen::Vector3d currentExplicitTorque = secondDegreeGravitationalTorque->getTorque( );


        sphercialHarmonicGravitationalTorque->updateMembers( evaluationTime );
        Eigen::Vector3d currentSphericalHarmonicTorque = sphercialHarmonicGravitationalTorque->getTorque( );

        std::cout<<"Explicit: "<<currentExplicitTorque.transpose( )<<std::endl;
        std::cout<<"SH      : "<<currentSphericalHarmonicTorque.transpose( )<<std::endl;
        std::cout<<( currentSphericalHarmonicTorque - currentExplicitTorque ).transpose( )<<std::endl;
        std::cout<<( ( currentSphericalHarmonicTorque - currentExplicitTorque ).cwiseQuotient( currentExplicitTorque ) ).transpose( )
                <<std::endl;
    }


}

BOOST_AUTO_TEST_SUITE_END( )

} // namespace unit_tests
} // namespace tudat
