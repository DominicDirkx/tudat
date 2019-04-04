#define BOOST_TEST_MAIN

#include <limits>

#include <boost/test/unit_test.hpp>

#include "Tudat/SimulationSetup/tudatSimulationHeader.h"
#include "Tudat/Astrodynamics/Gravitation/scaledGravitationalAccelerationModel.h"


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
                        if( expectScaledAccelerations )
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
                    std::cout<<"B"<<std::endl;

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
                    if( expectScaledAccelerations )
                    {
                        checkGravitationalAccelerationScaling(
                                    scaledAcceleration, bodyMap, fullAccelerationIterator->first, accelerationIterator->first );
                    }

                    std::cout<<"C"<<std::endl;

                    BOOST_CHECK_EQUAL( ( castToBaseAcceleration< BaseAcceleration >(
                                             thirdBodyAcceleration->getAccelerationModelForCentralBody( ) ) != nullptr ), 1 );
                    BOOST_CHECK_EQUAL( ( castToThirdBodyAcceleration< BaseAcceleration >(
                                             thirdBodyAcceleration->getAccelerationModelForCentralBody( ) ) == nullptr ), 1 );

                    scaledAcceleration = castToScaledAcceleration< BaseAcceleration >(
                                thirdBodyAcceleration->getAccelerationModelForCentralBody( ) );
                    std::cout<<"D"<<std::endl;

                    bool expectThirdBodyScaledAcceleration = 0;
                    if( std::find( calculatedExertingBodies.begin( ), calculatedExertingBodies.end( ), accelerationIterator->first ) != calculatedExertingBodies.end( ) ||
                            expectScaledAccelerations )
                    {
                        if( expectScaledAccelerations )
                        {
                            expectThirdBodyScaledAcceleration = 1;
                        }
                    }
                    std::cout<<"E "<<expectThirdBodyScaledAcceleration<<std::endl;

                    BOOST_CHECK_EQUAL( ( scaledAcceleration != nullptr ), expectThirdBodyScaledAcceleration );
                    if( expectThirdBodyScaledAcceleration &&  ( scaledAcceleration != nullptr ) )
                    {
                        if( accelerationIterator->first == accelerationModelMap.begin( )->first && expectScaledAccelerations )
                        {
                            std::cout<<"G"<<std::endl;
                            checkGravitationalAccelerationScaling(
                                        scaledAcceleration, bodyMap, "Jupiter", accelerationIterator->first );
                        }
                        else
                        {
                            std::cout<<"H"<<std::endl;
                            checkGravitationalAccelerationScaling(
                                        scaledAcceleration, bodyMap, accelerationIterator->first, accelerationIterator->first, 0 );
                        }
                    }
                    std::cout<<"F"<<std::endl;

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
    bodyNames.push_back( "Io" );
    bodyNames.push_back( "Jupiter" );
    //    bodyNames.push_back( "Ganymede" );
//    bodyNames.push_back( "Sun" );


    std::map< std::string, std::shared_ptr< BodySettings > > bodySettings =
            getDefaultBodySettings( bodyNames );

    // Set spherical harmonic gravity fields (dummy coefficients) for all bodies in Jovian system
    Eigen::MatrixXd cosineCoefficients = Eigen::MatrixXd::Zero( 5, 5 );
    Eigen::MatrixXd sineCoefficients = Eigen::MatrixXd::Zero( 5, 5 );
    for( unsigned int i = 1; i < bodyNames.size( ); i++ )
    {
        bodySettings[ bodyNames.at( i ) ]->gravityFieldSettings =
                std::make_shared< SphericalHarmonicsGravityFieldSettings >(
                    tudat::spice_interface::getBodyGravitationalParameter( bodyNames.at( i ) ),
                    1.0E5, cosineCoefficients, sineCoefficients, "IAU_" + bodyNames.at( i ) + "_SIMPLIFIED" );
    }

    // Finalize body creation
    NamedBodyMap bodyMap = createBodies( bodySettings );
    setGlobalFrameBodyEphemerides( bodyMap, "SSB", "ECLIPJ2000" );


    std::vector< std::string > acceleratedBodies;
    acceleratedBodies.push_back( "Europa" );
//    acceleratedBodies.push_back( "Ganymede" );
    acceleratedBodies.push_back( "Io" );

    for( int useScaledAccelerations = 1; useScaledAccelerations < 2; useScaledAccelerations++ )
    {
        std::cout<<"Central gravity; scaled "<<useScaledAccelerations<<std::endl;

        std::map< std::string, std::string > centralBodyMap;
        SelectedAccelerationMap accelerationMap;
        for( unsigned int i = 0; i < acceleratedBodies.size( ); i++ )
        {
            centralBodyMap[ acceleratedBodies.at( i ) ] = "Jupiter";
            std::map< std::string, std::vector< std::shared_ptr< AccelerationSettings > > > accelerationsOfCurrentMoon;

            for( unsigned int j = 0; j < bodyNames.size( ); j++ )
            {
                if( acceleratedBodies.at( i ) != bodyNames.at( j ) )
                {
                    accelerationsOfCurrentMoon[ bodyNames.at( j ) ].push_back(
                                std::make_shared< AccelerationSettings >( central_gravity ) );
                }
            }

            accelerationMap[ acceleratedBodies.at( i ) ] = accelerationsOfCurrentMoon;
            accelerationsOfCurrentMoon.clear( );
        }

        AccelerationMap accelerationModelMap = createAccelerationModelsMap(
                    bodyMap, accelerationMap, centralBodyMap, useScaledAccelerations );

        std::cout<<std::endl;

        verifyAccelerationModelTypesAndScalings< CentralGravitationalAccelerationModel3d >(
                    accelerationModelMap, bodyMap, useScaledAccelerations );

    }

//    {
//        std::map< std::string, std::string > centralBodyMap;

//        SelectedAccelerationMap accelerationMap;
//        for( unsigned int i = 0; i < acceleratedBodies.size( ); i++ )
//        {
//            centralBodyMap[ acceleratedBodies.at( i ) ] = "Jupiter";
//            std::map< std::string, std::vector< std::shared_ptr< AccelerationSettings > > > accelerationsOfCurrentMoon;

//            for( unsigned int j = 1; j < bodyNames.size( ); j++ )
//            {
//                if( acceleratedBodies.at( i ) != bodyNames.at( j ) )
//                {
//                    accelerationsOfCurrentMoon[ bodyNames.at( j ) ].push_back( std::make_shared< SphericalHarmonicAccelerationSettings >( 2, 2 ) );
//                }
//            }

//            accelerationMap[ acceleratedBodies.at( i ) ] = accelerationsOfCurrentMoon;
//            accelerationsOfCurrentMoon.clear( );
//        }

//        AccelerationMap accelerationModelMap = createAccelerationModelsMap(
//                    bodyMap, accelerationMap, centralBodyMap );

//        verifyAccelerationModelTypesAndScalings< SphericalHarmonicsGravitationalAccelerationModel >(
//                    accelerationModelMap, bodyMap, 0 );
//    }
//    std::cout<<"AAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAA"<<std::endl;

//    {
//        std::map< std::string, std::string > centralBodyMap;

//        SelectedAccelerationMap accelerationMap;
//        for( unsigned int i = 0; i < acceleratedBodies.size( ); i++ )
//        {
//            centralBodyMap[ acceleratedBodies.at( i ) ] = "Jupiter";
//            std::map< std::string, std::vector< std::shared_ptr< AccelerationSettings > > > accelerationsOfCurrentMoon;

//            for( unsigned int j = 1; j < bodyNames.size( ); j++ )
//            {
//                if( acceleratedBodies.at( i ) != bodyNames.at( j ) )
//                {
//                    if( bodyNames.at( j ) == "Jupiter" )
//                    {
//                        accelerationsOfCurrentMoon[ bodyNames.at( j ) ].push_back( std::make_shared< MutualSphericalHarmonicAccelerationSettings >( 4, 4, 2, 2 ) );
//                    }
//                    else
//                    {
//                        accelerationsOfCurrentMoon[ bodyNames.at( j ) ].push_back( std::make_shared< MutualSphericalHarmonicAccelerationSettings >( 2, 2, 2, 2, 4, 4 ) );
//                    }
//                }
//            }

//            accelerationMap[ acceleratedBodies.at( i ) ] = accelerationsOfCurrentMoon;
//            accelerationsOfCurrentMoon.clear( );
//        }

//        AccelerationMap accelerationModelMap = createAccelerationModelsMap(
//                    bodyMap, accelerationMap, centralBodyMap );

//        verifyAccelerationModelTypesAndScalings< MutualSphericalHarmonicsGravitationalAccelerationModel >(
//                    accelerationModelMap, bodyMap );
//    }
}

BOOST_AUTO_TEST_SUITE_END( )


}

}
