#define BOOST_TEST_MAIN

#include <limits>

#include <boost/test/unit_test.hpp>

#include "Tudat/Basics/testMacros.h"
#include "Tudat/SimulationSetup/tudatSimulationHeader.h"
#include "Tudat/Astrodynamics/Gravitation/scaledGravitationalAccelerationModel.h"
#include "Tudat/Mathematics/Statistics/randomVariableGenerator.h"


namespace tudat
{

namespace unit_tests
{

using namespace tudat::basic_astrodynamics;
using namespace tudat::simulation_setup;
using namespace tudat::gravitation;
using namespace tudat::basic_astrodynamics;


BOOST_AUTO_TEST_SUITE( test_gravitational_acceleration_setup )



template< class BaseAcceleration >
std::shared_ptr< BaseAcceleration > castToBaseAcceleration(
        const std::shared_ptr< AccelerationModel< Eigen::Vector3d > > accelerationModel )
{
    return std::dynamic_pointer_cast< BaseAcceleration >( accelerationModel );
}

template< class BaseAcceleration >
std::shared_ptr< ThirdBodyAcceleration< BaseAcceleration > > castToThirdBodyAcceleration(
        const std::shared_ptr< AccelerationModel< Eigen::Vector3d > > accelerationModel )
{
    return std::dynamic_pointer_cast< ThirdBodyAcceleration< BaseAcceleration > >( accelerationModel );
}

template< class BaseAcceleration >
std::shared_ptr< ScaledGravitationalAccelerationModel< BaseAcceleration > > castToScaledAcceleration(
        const std::shared_ptr< AccelerationModel< Eigen::Vector3d > > accelerationModel )
{
    return std::dynamic_pointer_cast< ScaledGravitationalAccelerationModel< BaseAcceleration > >( accelerationModel );
}

double getBodyGravitationalParameter(
        const::std::shared_ptr< Body > body )
{
    return body->getGravityFieldModel( )->getGravitationalParameter( );
}

template< class BaseAcceleration >
double getScaledAccelerationGravitationalParameter(
        const std::shared_ptr< ScaledGravitationalAccelerationModel< BaseAcceleration > > scaledAcceleration );

template< >
double getScaledAccelerationGravitationalParameter< MutualSphericalHarmonicsGravitationalAccelerationModel >(
        const std::shared_ptr< ScaledGravitationalAccelerationModel< MutualSphericalHarmonicsGravitationalAccelerationModel > > scaledAcceleration )
{
    double gravitationalParameterOfBodyUndergoingAcceleration =
            scaledAcceleration->getAccelerationModelFromShExpansionOfBodyUndergoingAcceleration( )->getGravitationalParameterFunction( )( );
    double gravitationalParameterOfBodyExertingAcceleration =
            scaledAcceleration->getAccelerationModelFromShExpansionOfBodyExertingAcceleration( )->getGravitationalParameterFunction( )( );

    BOOST_CHECK_EQUAL( gravitationalParameterOfBodyUndergoingAcceleration, gravitationalParameterOfBodyExertingAcceleration );
    return gravitationalParameterOfBodyUndergoingAcceleration;
}

template< >
double getScaledAccelerationGravitationalParameter< SphericalHarmonicsGravitationalAccelerationModel >(
        const std::shared_ptr< ScaledGravitationalAccelerationModel< SphericalHarmonicsGravitationalAccelerationModel > > scaledAcceleration )
{
    return scaledAcceleration->getGravitationalParameterFunction( )( );
}

template< >
double getScaledAccelerationGravitationalParameter< CentralGravitationalAccelerationModel3d >(
        const std::shared_ptr< ScaledGravitationalAccelerationModel< CentralGravitationalAccelerationModel3d > > scaledAcceleration )
{
    return scaledAcceleration->getGravitationalParameterFunction( )( );
}

template< class BaseAcceleration >
double getScaledAccelerationGravitationalParameterFromBase(
        const std::shared_ptr< ScaledGravitationalAccelerationModel< BaseAcceleration > > scaledAcceleration );

template< >
double getScaledAccelerationGravitationalParameterFromBase< SphericalHarmonicsGravitationalAccelerationModel >(
        const std::shared_ptr< ScaledGravitationalAccelerationModel< SphericalHarmonicsGravitationalAccelerationModel > > scaledAcceleration )
{
    return scaledAcceleration->getOriginalAccelerationModel( )->getGravitationalParameterFunction( )( );
}

template< >
double getScaledAccelerationGravitationalParameterFromBase< CentralGravitationalAccelerationModel3d >(
        const std::shared_ptr< ScaledGravitationalAccelerationModel< CentralGravitationalAccelerationModel3d > > scaledAcceleration )
{
    return scaledAcceleration->getOriginalAccelerationModel( )->getGravitationalParameterFunction( )( );
}

template< >
double getScaledAccelerationGravitationalParameterFromBase< MutualSphericalHarmonicsGravitationalAccelerationModel >(
        const std::shared_ptr< ScaledGravitationalAccelerationModel< MutualSphericalHarmonicsGravitationalAccelerationModel > > scaledAcceleration )
{
    double gravitationalParameterOfBodyUndergoingAccelerationBase =
            getScaledAccelerationGravitationalParameterFromBase(
                std::dynamic_pointer_cast< ScaledGravitationalAccelerationModel< SphericalHarmonicsGravitationalAccelerationModel > >(
                    scaledAcceleration->getAccelerationModelFromShExpansionOfBodyUndergoingAcceleration( ) ) );
    double gravitationalParameterOfBodyExertingAccelerationBase =
            getScaledAccelerationGravitationalParameterFromBase(
                std::dynamic_pointer_cast< ScaledGravitationalAccelerationModel< SphericalHarmonicsGravitationalAccelerationModel > >(
                    scaledAcceleration->getAccelerationModelFromShExpansionOfBodyExertingAcceleration( ) ) );

    BOOST_CHECK_EQUAL( gravitationalParameterOfBodyUndergoingAccelerationBase, gravitationalParameterOfBodyExertingAccelerationBase );
    return gravitationalParameterOfBodyUndergoingAccelerationBase;
}


template< class BaseAcceleration >
void checkGravitationalAccelerationScaling(
        const std::shared_ptr< ScaledGravitationalAccelerationModel< BaseAcceleration > > scaledAcceleration,
        const NamedBodyMap& bodyMap,
        const std::string bodyUndergoingAcceleration,
        const std::string bodyExertingAcceleration,
        const bool expectedInvertPositions = 1 )
{
    double gravitationalParameterOfUndergoingBody = getBodyGravitationalParameter( bodyMap.at( bodyUndergoingAcceleration ) );
    double gravitationalParameterOfExertingBody = getBodyGravitationalParameter( bodyMap.at( bodyExertingAcceleration ) );

    if( bodyUndergoingAcceleration == "Jupiter" )
    {
        gravitationalParameterOfUndergoingBody += gravitationalParameterOfExertingBody;
    }

    if( bodyExertingAcceleration == "Jupiter" )
    {
        gravitationalParameterOfExertingBody += gravitationalParameterOfUndergoingBody;
    }

    BOOST_CHECK_CLOSE_FRACTION(
                gravitationalParameterOfExertingBody, getScaledAccelerationGravitationalParameter( scaledAcceleration ),
                std::numeric_limits< double >::epsilon( ) );

    BOOST_CHECK_CLOSE_FRACTION(
                gravitationalParameterOfUndergoingBody, getScaledAccelerationGravitationalParameterFromBase( scaledAcceleration ),
                std::numeric_limits< double >::epsilon( ) );

    BOOST_CHECK_EQUAL( scaledAcceleration->getInvertPositionVectors( ), expectedInvertPositions );
}

template< class BaseAcceleration >
void verifyAccelerationModelTypesAndScalings(
        const AccelerationMap& accelerationModelMap,
        const NamedBodyMap& bodyMap,
        const bool expectScaledAccelerations = 1 )
{
    std::vector< std::string > calculatedBodies;
    std::vector< std::string > calculatedExertingBodies;

    std::cout<<"TESTING"<<std::endl;
    bool isFirstIterationDone = 0;

    // Iterate over all bodies undergoing acceleration
    for( AccelerationMap::const_iterator fullAccelerationIterator = accelerationModelMap.begin( );
         fullAccelerationIterator != accelerationModelMap.end( ); fullAccelerationIterator++ )
    {
        std::cout<<"undergoing: "<<fullAccelerationIterator->first<<std::endl;

        // Iterate over all bodies exerting acceleration
        for( SingleBodyAccelerationMap::const_iterator accelerationIterator = fullAccelerationIterator->second.begin( );
             accelerationIterator != fullAccelerationIterator->second.end( ); accelerationIterator++ )
        {
            std::cout<<"exerting: "<<accelerationIterator->first<<std::endl;

            // Iterate over all accelerations in current list
            for( unsigned int i = 0; i < accelerationIterator->second.size( ); i++ )
            {
                // Test if body exerting acceleration is Jupiter
                if( accelerationIterator->first == "Jupiter" )
                {
                    // Acceleration is not third-body
                    BOOST_CHECK_EQUAL( ( castToBaseAcceleration< BaseAcceleration >( accelerationIterator->second.at( i ) )
                                         != nullptr ), 1 );
                    BOOST_CHECK_EQUAL( ( castToThirdBodyAcceleration< BaseAcceleration >(
                                             accelerationIterator->second.at( i ) ) == nullptr ), 1 );

                    if( !isFirstIterationDone )
                    {
                        BOOST_CHECK_EQUAL( ( castToScaledAcceleration< BaseAcceleration >(
                                                 accelerationIterator->second.at( i ) ) == nullptr ), 1 );
                    }
                    else
                    {
                        std::shared_ptr< ScaledGravitationalAccelerationModel< BaseAcceleration > > scaledAcceleration =
                                castToScaledAcceleration< BaseAcceleration >(
                                    accelerationIterator->second.at( i ) );

                        BOOST_CHECK_EQUAL( ( scaledAcceleration != nullptr ), expectScaledAccelerations );
                        if( expectScaledAccelerations && ( scaledAcceleration != nullptr ) )
                        {
                            checkGravitationalAccelerationScaling(
                                        scaledAcceleration, bodyMap, fullAccelerationIterator->first, accelerationIterator->first );
                        }

                    }
                }
                else if( std::find( calculatedBodies.begin( ), calculatedBodies.end( ), accelerationIterator->first ) ==
                         calculatedBodies.end( ) )
                {
                    // Check that acceleration is not base acceleration, or scaled acceleration
                    BOOST_CHECK_EQUAL( ( castToBaseAcceleration< BaseAcceleration >(
                                             accelerationIterator->second.at( i ) ) == nullptr ), 1 );
                    BOOST_CHECK_EQUAL( ( castToScaledAcceleration< BaseAcceleration >(
                                             accelerationIterator->second.at( i ) ) == nullptr ), 1 );

                    // Check that acceleration is third-body acceleration
                    std::shared_ptr< ThirdBodyAcceleration< BaseAcceleration > > thirdBodyAcceleration =
                            castToThirdBodyAcceleration< BaseAcceleration >(
                                accelerationIterator->second.at( i ) );
                    BOOST_CHECK_EQUAL( ( thirdBodyAcceleration != nullptr ), 1 );

                    // Check that acceleration of body undergoing acceleration is a base acceleration
                    BOOST_CHECK_EQUAL(
                                ( castToBaseAcceleration< BaseAcceleration >(
                                      thirdBodyAcceleration->getAccelerationModelForBodyUndergoingAcceleration( ) ) != nullptr ), 1 );
                    BOOST_CHECK_EQUAL(
                                ( castToThirdBodyAcceleration< BaseAcceleration >(
                                      thirdBodyAcceleration->getAccelerationModelForBodyUndergoingAcceleration( ) ) == nullptr ), 1 );
                    BOOST_CHECK_EQUAL(
                                ( castToScaledAcceleration< BaseAcceleration >(
                                      thirdBodyAcceleration->getAccelerationModelForBodyUndergoingAcceleration( ) ) == nullptr ), 1 );

                    // Check that acceleration of body undergoing acceleration is a base acceleration
                    BOOST_CHECK_EQUAL( ( castToBaseAcceleration< BaseAcceleration >(
                                             thirdBodyAcceleration->getAccelerationModelForCentralBody( ) ) != nullptr ), 1 );
                    BOOST_CHECK_EQUAL( ( castToThirdBodyAcceleration< BaseAcceleration >(
                                             thirdBodyAcceleration->getAccelerationModelForCentralBody( ) ) == nullptr ), 1 );

                    // If first iteration is done, acceleration acing on central body is scaled acceleration
                    if( !isFirstIterationDone || !expectScaledAccelerations )
                    {
                        BOOST_CHECK_EQUAL( ( castToScaledAcceleration< BaseAcceleration >(
                                                 thirdBodyAcceleration->getAccelerationModelForCentralBody( ) ) == nullptr ), 1 );
                    }
                    else
                    {
                        BOOST_CHECK_EQUAL( ( castToScaledAcceleration< BaseAcceleration >(
                                                 thirdBodyAcceleration->getAccelerationModelForCentralBody( ) ) != nullptr ), 1 );
                    }

                }
                else
                {
                    BOOST_CHECK_EQUAL( ( castToBaseAcceleration< BaseAcceleration >(
                                             accelerationIterator->second.at( i ) ) == nullptr ), 1 );
                    BOOST_CHECK_EQUAL( ( castToScaledAcceleration< BaseAcceleration >(
                                             accelerationIterator->second.at( i ) ) == nullptr ), 1 );

                    std::shared_ptr< ThirdBodyAcceleration< BaseAcceleration > > thirdBodyAcceleration =
                            castToThirdBodyAcceleration< BaseAcceleration >(
                                accelerationIterator->second.at( i ) );
                    BOOST_CHECK_EQUAL( ( thirdBodyAcceleration != nullptr ), 1 );
                    BOOST_CHECK_EQUAL( ( castToBaseAcceleration< BaseAcceleration >(
                                             thirdBodyAcceleration->getAccelerationModelForBodyUndergoingAcceleration( ) ) != nullptr ), 1 );
                    BOOST_CHECK_EQUAL( ( castToThirdBodyAcceleration< BaseAcceleration >(
                                             thirdBodyAcceleration->getAccelerationModelForBodyUndergoingAcceleration( ) ) == nullptr ), 1 );

                    std::shared_ptr< ScaledGravitationalAccelerationModel< BaseAcceleration > > scaledAcceleration =
                            castToScaledAcceleration< BaseAcceleration >(
                                thirdBodyAcceleration->getAccelerationModelForBodyUndergoingAcceleration( ) );

                    BOOST_CHECK_EQUAL( ( scaledAcceleration != nullptr ), expectScaledAccelerations );
                    if( expectScaledAccelerations && ( scaledAcceleration != nullptr ) )
                    {
                        checkGravitationalAccelerationScaling(
                                    scaledAcceleration, bodyMap, fullAccelerationIterator->first, accelerationIterator->first );
                    }

                    BOOST_CHECK_EQUAL( ( castToBaseAcceleration< BaseAcceleration >(
                                             thirdBodyAcceleration->getAccelerationModelForCentralBody( ) ) != nullptr ), 1 );
                    BOOST_CHECK_EQUAL( ( castToThirdBodyAcceleration< BaseAcceleration >(
                                             thirdBodyAcceleration->getAccelerationModelForCentralBody( ) ) == nullptr ), 1 );

                    scaledAcceleration = castToScaledAcceleration< BaseAcceleration >(
                                thirdBodyAcceleration->getAccelerationModelForCentralBody( ) );

                    bool expectThirdBodyScaledAcceleration = 0;
                    if( std::find( calculatedExertingBodies.begin( ), calculatedExertingBodies.end( ), accelerationIterator->first ) != calculatedExertingBodies.end( ) ||
                            expectScaledAccelerations )
                    {
                        if( expectScaledAccelerations )
                        {
                            expectThirdBodyScaledAcceleration = 1;
                        }
                    }

                    BOOST_CHECK_EQUAL( ( scaledAcceleration != nullptr ), expectThirdBodyScaledAcceleration );
                    if( expectThirdBodyScaledAcceleration &&  ( scaledAcceleration != nullptr ) )
                    {
                        if( accelerationIterator->first == accelerationModelMap.begin( )->first && expectScaledAccelerations )
                        {
                            checkGravitationalAccelerationScaling(
                                        scaledAcceleration, bodyMap, "Jupiter", accelerationIterator->first );
                        }
                        else
                        {
                            checkGravitationalAccelerationScaling(
                                        scaledAcceleration, bodyMap, accelerationIterator->first, accelerationIterator->first, 0 );
                        }
                    }
                }
            }
            calculatedExertingBodies.push_back( accelerationIterator->first );
        }
        calculatedBodies.push_back( fullAccelerationIterator->first );
        isFirstIterationDone = 1;
    }
}

BOOST_AUTO_TEST_CASE( testIdenticalGravitySetup )
{
    // Load Spice kernels
    tudat::spice_interface::loadStandardSpiceKernels( );
    tudat::spice_interface::loadSpiceKernelInTudat( input_output::getSpiceKernelPath( ) + "jup310_small.bsp" );

    // Create bodies needed in simulation
    std::vector< std::string > bodyNames;
    bodyNames.push_back( "Europa" );
    bodyNames.push_back( "Ganymede" );
    bodyNames.push_back( "Io" );
    bodyNames.push_back( "Jupiter" );


    std::map< std::string, std::shared_ptr< BodySettings > > bodySettings =
            getDefaultBodySettings( bodyNames );

    std::function< double( ) > randomNumberGenerator =
            statistics::createBoostContinuousRandomVariableGeneratorFunction(
                statistics::normal_boost_distribution,
    { 0.0, 1.0E-4 }, 0.0 );

    // Set spherical harmonic gravity fields (dummy coefficients) for all bodies in Jovian system

    for( unsigned int i = 0; i < bodyNames.size( ); i++ )
    {
        Eigen::MatrixXd cosineCoefficients = Eigen::MatrixXd::Zero( 5, 5 );
        Eigen::MatrixXd sineCoefficients = Eigen::MatrixXd::Zero( 5, 5 );
        cosineCoefficients( 0, 0 ) = 1.0;
        for( int l = 2; l < 5; l++ )
        {
            for( int m = 0; m <= l; m++ )
            {
                cosineCoefficients( l, m ) = randomNumberGenerator( );
                sineCoefficients( l, m ) = randomNumberGenerator( );
            }
        }

        bodySettings[ bodyNames.at( i ) ]->gravityFieldSettings =
                std::make_shared< SphericalHarmonicsGravityFieldSettings >(
                    tudat::spice_interface::getBodyGravitationalParameter( bodyNames.at( i ) ),
                    tudat::spice_interface::getAverageRadius( bodyNames.at( i ) ), cosineCoefficients,
                    sineCoefficients, "IAU_" + bodyNames.at( i ) );
    }

    // Finalize body creation
    NamedBodyMap bodyMap = createBodies( bodySettings );
    setGlobalFrameBodyEphemerides( bodyMap, "SSB", "ECLIPJ2000" );


    std::vector< std::string > acceleratedBodies;
    acceleratedBodies.push_back( "Europa" );
    acceleratedBodies.push_back( "Ganymede" );
    acceleratedBodies.push_back( "Io" );


    std::map< std::string, std::string > centralBodyMap;
    SelectedAccelerationMap centralAccelerationMap;
    SelectedAccelerationMap sphericalHarmonicAccelerationMap;
    SelectedAccelerationMap mutualSphericalHarmonicAccelerationMap;

    for( unsigned int i = 0; i < acceleratedBodies.size( ); i++ )
    {
        centralBodyMap[ acceleratedBodies.at( i ) ] = "Jupiter";
        std::map< std::string, std::vector< std::shared_ptr< AccelerationSettings > > > centralAccelerationsOfCurrentMoon;
        std::map< std::string, std::vector< std::shared_ptr< AccelerationSettings > > > sphericalHarmonicAccelerationsOfCurrentMoon;
        std::map< std::string, std::vector< std::shared_ptr< AccelerationSettings > > > mutualSphericalHarmonicAccelerationsOfCurrentMoon;

        for( unsigned int j = 0; j < bodyNames.size( ); j++ )
        {
            if( acceleratedBodies.at( i ) != bodyNames.at( j ) )
            {
                centralAccelerationsOfCurrentMoon[ bodyNames.at( j ) ].push_back(
                            std::make_shared< AccelerationSettings >( central_gravity ) );
                sphericalHarmonicAccelerationsOfCurrentMoon[ bodyNames.at( j ) ].push_back(
                            std::make_shared< SphericalHarmonicAccelerationSettings >( 2, 2 ) );

                if( bodyNames.at( j ) == "Jupiter" )
                {
                    mutualSphericalHarmonicAccelerationsOfCurrentMoon[ bodyNames.at( j ) ].push_back(
                                std::make_shared< MutualSphericalHarmonicAccelerationSettings >( 4, 4, 2, 2 ) );
                }
                else
                {
                    mutualSphericalHarmonicAccelerationsOfCurrentMoon[ bodyNames.at( j ) ].push_back(
                                std::make_shared< MutualSphericalHarmonicAccelerationSettings >( 2, 2, 2, 2, 4, 4 ) );
                }
            }

        }

        centralAccelerationMap[ acceleratedBodies.at( i ) ] = centralAccelerationsOfCurrentMoon;
        sphericalHarmonicAccelerationMap[ acceleratedBodies.at( i ) ] = sphericalHarmonicAccelerationsOfCurrentMoon;
        mutualSphericalHarmonicAccelerationMap[ acceleratedBodies.at( i ) ] = mutualSphericalHarmonicAccelerationsOfCurrentMoon;

        centralAccelerationsOfCurrentMoon.clear( );
    }

    std::vector< AccelerationMap > accelerationModelMaps;

    accelerationModelMaps.push_back( createAccelerationModelsMap(
                                         bodyMap, centralAccelerationMap, centralBodyMap, false ) );
    accelerationModelMaps.push_back( createAccelerationModelsMap(
                                         bodyMap, centralAccelerationMap, centralBodyMap, true ) );
    accelerationModelMaps.push_back( createAccelerationModelsMap(
                                         bodyMap, sphericalHarmonicAccelerationMap, centralBodyMap, false ) );
    accelerationModelMaps.push_back( createAccelerationModelsMap(
                                         bodyMap, sphericalHarmonicAccelerationMap, centralBodyMap, true ) );
    accelerationModelMaps.push_back( createAccelerationModelsMap(
                                         bodyMap, mutualSphericalHarmonicAccelerationMap, centralBodyMap, false ) );
    accelerationModelMaps.push_back( createAccelerationModelsMap(
                                         bodyMap, mutualSphericalHarmonicAccelerationMap, centralBodyMap, true ) );
    for( int i = 0; i < 6; i++ )
    {
        std::cout<<"******************************************************* "<<i<<std::endl;
        if( i < 2 )
        {
            verifyAccelerationModelTypesAndScalings< CentralGravitationalAccelerationModel3d >(
                        accelerationModelMaps.at( i ), bodyMap, i % 2 );
        }
        else if( i < 4 )
        {
            verifyAccelerationModelTypesAndScalings< SphericalHarmonicsGravitationalAccelerationModel >(
                        accelerationModelMaps.at( i ), bodyMap, i % 2 );
        }
        else
        {
            verifyAccelerationModelTypesAndScalings< MutualSphericalHarmonicsGravitationalAccelerationModel >(
                        accelerationModelMaps.at( i ), bodyMap, i % 2 );
        }
    }

    for( int i = 0; i < 100; i++ )
    {
        double currentTime = 86400.0 * static_cast< double >( i );
        for( unsigned int j = 0; j < bodyNames.size( ); j++ )
        {
            bodyMap.at( bodyNames.at( j ) )->setStateFromEphemeris( currentTime );
            bodyMap.at( bodyNames.at( j ) )->setCurrentRotationalStateToLocalFrameFromEphemeris( currentTime );
        }

        for( int j = 0; j < 3; j++ )
        {
            AccelerationMap unscaledAccelerationMap = accelerationModelMaps.at( 2 * j );
            AccelerationMap scaledAccelerationMap = accelerationModelMaps.at( 2 * j + 1 );

            auto unscaledIterator = unscaledAccelerationMap.begin( );
            auto scaledIterator = scaledAccelerationMap.begin( );

            for( unsigned int k = 0; k < scaledAccelerationMap.size( ); k++ )
            {
                auto unscaledInnerIterator = unscaledIterator->second.begin( );
                auto scaledInnerIterator = scaledIterator->second.begin( );

                for( unsigned int l = 0; l < scaledAccelerationMap.size( ); l++ )
                {
                    std::vector< std::shared_ptr< basic_astrodynamics::AccelerationModel< Eigen::Vector3d > > >
                            scaledAccelerationList = scaledInnerIterator->second;
                    std::vector< std::shared_ptr< basic_astrodynamics::AccelerationModel< Eigen::Vector3d > > >
                            unscaledAccelerationList = unscaledInnerIterator->second;

                    for( int m = 0; m < scaledAccelerationList.size( ); m++ )
                    {
                        scaledAccelerationList.at( m )->resetTime( TUDAT_NAN );
                        unscaledAccelerationList.at( m )->resetTime( TUDAT_NAN );

                        Eigen::Vector3d scaledAcceleration = updateAndGetAcceleration(
                                    scaledAccelerationList.at( m ), currentTime );
                        Eigen::Vector3d unscaledAcceleration = updateAndGetAcceleration(
                                    unscaledAccelerationList.at( m ), currentTime );

                        TUDAT_CHECK_MATRIX_CLOSE_FRACTION( scaledAcceleration, unscaledAcceleration, 1.0E-15 );

                        std::cout<<unscaledIterator->first<<" "<<unscaledInnerIterator->first<<" "<<j<<std::endl;
                        std::cout<<std::setprecision( 16 )<<scaledAcceleration.transpose( )<<std::endl;
                        std::cout<<std::setprecision( 16 )<<unscaledAcceleration.transpose( )<<std::endl<<std::endl;

                    }

                    unscaledInnerIterator++;
                    scaledInnerIterator++;
                }

                unscaledIterator++;
                scaledIterator++;
            }

        }

    }

}

BOOST_AUTO_TEST_SUITE_END( )


}

}
