#define BOOST_TEST_MAIN

#include <string>
#include <thread>

#include <boost/test/unit_test.hpp>

#include "Tudat/Basics/testMacros.h"

#include "Tudat/SimulationSetup/tudatSimulationHeader.h"
#include "Tudat/Mathematics/Statistics/randomVariableGenerator.h"


namespace tudat
{

namespace unit_tests
{

using namespace tudat::basic_astrodynamics;
using namespace tudat::simulation_setup;
using namespace tudat::spice_interface;
using namespace tudat::statistics;
using namespace tudat::gravitation;
using namespace tudat;

std::pair< Eigen::MatrixXd, Eigen::MatrixXd > generateCosineSineCoefficients(
        boost::shared_ptr< RandomVariableGenerator< double > > randomNumberGenerator,
        const int maximumDegree, const int maximumOrder )
{
    Eigen::MatrixXd cosineCoefficients = Eigen::MatrixXd::Zero( maximumDegree + 1, maximumOrder + 1 );
    Eigen::MatrixXd sineCoefficients = Eigen::MatrixXd::Zero( maximumDegree + 1, maximumOrder + 1 );

    for( int i = 0; i <= maximumDegree; i++ )
    {
        for( int j = 0; ( ( j <= i ) && ( j <= maximumOrder ) ); j++ )
        {
            cosineCoefficients( i, j ) = randomNumberGenerator->getRandomVariableValue( );
            if( j != 0 )
            {
                sineCoefficients( i, j ) = randomNumberGenerator->getRandomVariableValue( );
            }
        }
    }

    //cosineCoefficients( 1, 0 ) = randomNumberGenerator->getRandomVariableValue( );
    //    sineCoefficients( 2, 1 ) = randomNumberGenerator->getRandomVariableValue( );

    cosineCoefficients( 0, 0 ) = 1.0;


    return std::make_pair( cosineCoefficients, sineCoefficients );
}

boost::shared_ptr< tudat::simulation_setup::GravityFieldSettings > getDummyJovianSystemGravityField(
        const std::string& bodyName,
        const int isNormalized )
{
    boost::shared_ptr< GravityFieldSettings > gravityFieldSettings;

    std::vector< double > randomNumberSettings;
    randomNumberSettings.push_back( 0.0 );
    randomNumberSettings.push_back( 1.0E-4 );

    std::pair< Eigen::MatrixXd, Eigen::MatrixXd > coefficients;

    if( bodyName == "Jupiter" )
    {
        boost::shared_ptr< RandomVariableGenerator< double > > randomCoefficientGenerator = createBoostContinuousRandomVariableGenerator(
                    normal_boost_distribution, randomNumberSettings , 0.0 );
        coefficients = generateCosineSineCoefficients( randomCoefficientGenerator, 10, 10 );
        gravityFieldSettings = boost::make_shared< SphericalHarmonicsGravityFieldSettings >
                ( getBodyGravitationalParameter( "Jupiter" ), getAverageRadius( "Jupiter" ),
                  coefficients.first, coefficients.second, "IAU_Jupiter_SIMPLIFIED" );
    }
    else if( bodyName == "Io" )
    {
        boost::shared_ptr< RandomVariableGenerator< double > > randomCoefficientGenerator = createBoostContinuousRandomVariableGenerator(
                    normal_boost_distribution, randomNumberSettings , 1.0 );
        coefficients = generateCosineSineCoefficients( randomCoefficientGenerator, 10, 10 );
        gravityFieldSettings = boost::make_shared< SphericalHarmonicsGravityFieldSettings >
                ( 5.959916033410404E012, getAverageRadius( "Io" ),
                  coefficients.first, coefficients.second, "IAU_Io_SIMPLIFIED" );
    }
    else if( bodyName == "Europa" )
    {
        boost::shared_ptr< RandomVariableGenerator< double > > randomCoefficientGenerator = createBoostContinuousRandomVariableGenerator(
                    normal_boost_distribution, randomNumberSettings , 2.0 );
        coefficients = generateCosineSineCoefficients( randomCoefficientGenerator, 10, 10 );
        gravityFieldSettings = boost::make_shared< SphericalHarmonicsGravityFieldSettings >
                ( 3.202738774922892E12, getAverageRadius( "Europa" ),
                  coefficients.first, coefficients.second, "IAU_Europa_SIMPLIFIED" );
    }

    return gravityFieldSettings;

}



//void getBodyCoefficientEulerAnglePartials(
//        const Eigen::Quaterniond nominalQuaternion,
//        boost::shared_ptr< gravitation::EffectiveMutualSphericalHarmonicsField > mutualShField,
//        std::vector< Eigen::MatrixXd >& numericalTransformedCosineCoefficientsOfBody2Partials,
//        std::vector< Eigen::MatrixXd >& numericalTransformedSineCoefficientsOfBody2Partials,
//        const double perturbation )
//{
//    numericalTransformedCosineCoefficientsOfBody2Partials.resize( 3 );
//    numericalTransformedSineCoefficientsOfBody2Partials.resize( 3 );


//    Eigen::Vector3d nominalEulerAngles = basic_mathematics::get313EulerAnglesFromQuaternion( nominalQuaternion );
//    Eigen::Vector3d perturbedEulerAngles;
//    Eigen::MatrixXd upPerturbedCosineCoefficients, upPerturbedSineCoefficients;
//    Eigen::MatrixXd downPerturbedCosineCoefficients, downPerturbedSineCoefficients;

//    for( unsigned int i = 0; i < 3; i ++ )
//    {
//        perturbedEulerAngles = nominalEulerAngles;
//        perturbedEulerAngles( i ) += perturbation;
//        mutualShField->computeCurrentEffectiveCoefficients( perturbedEulerAngles( 2 ), perturbedEulerAngles( 1 ), perturbedEulerAngles( 0 ) );
//        upPerturbedCosineCoefficients = mutualShField->getTransformedCosineCoefficientsOfBody2( );
//        upPerturbedSineCoefficients = mutualShField->getTransformedSineCoefficientsOfBody2( );


//        perturbedEulerAngles = nominalEulerAngles;
//        perturbedEulerAngles( i ) -= perturbation;
//        mutualShField->computeCurrentEffectiveCoefficients( perturbedEulerAngles( 2 ), perturbedEulerAngles( 1 ), perturbedEulerAngles( 0 ) );
//        downPerturbedCosineCoefficients = mutualShField->getTransformedCosineCoefficientsOfBody2( );
//        downPerturbedSineCoefficients = mutualShField->getTransformedSineCoefficientsOfBody2( );

//        numericalTransformedCosineCoefficientsOfBody2Partials[ i ] =
//                ( upPerturbedCosineCoefficients - downPerturbedCosineCoefficients ) / ( 2.0 * perturbation );
//        numericalTransformedSineCoefficientsOfBody2Partials[ i ] =
//                ( upPerturbedSineCoefficients - downPerturbedSineCoefficients ) / ( 2.0 * perturbation );

//    }
//}

Eigen::Matrix2d getEffectiveMutualPotentialCoefficientWrtBody2Coefficient(
        const int degreeOfBody1, const int orderOfBody1, const int degreeOfBody2, const int orderOfBody2,
        const Eigen::MatrixXd& nominalTransformedCosineCoefficients,
        const Eigen::MatrixXd& nominalTransformedSineCoefficients,
        const boost::shared_ptr< gravitation::EffectiveMutualSphericalHarmonicsField > mutualPotentialField,
        const double coefficientPerturbation )
{
    Eigen::MatrixXd perturbedTransformedCosineCoefficients = nominalTransformedCosineCoefficients;
    Eigen::MatrixXd perturbedTransformedSineCoefficients = nominalTransformedSineCoefficients;


    perturbedTransformedCosineCoefficients( degreeOfBody2, std::abs( orderOfBody2 ) ) += coefficientPerturbation;

    mutualPotentialField->computeCurrentEffectiveCoefficientsFromManualTransformedCoefficients(
                perturbedTransformedCosineCoefficients, perturbedTransformedSineCoefficients );
    double upPerturbedCosineCoefficientFromCosinePerturbation = mutualPotentialField->getEffectiveCosineCoefficient(
                degreeOfBody1, orderOfBody1, degreeOfBody2, orderOfBody2 );
    double upPerturbedSineCoefficientFromCosinePerturbation = mutualPotentialField->getEffectiveSineCoefficient(
                degreeOfBody1, orderOfBody1, degreeOfBody2, orderOfBody2 );



    perturbedTransformedCosineCoefficients = nominalTransformedCosineCoefficients;
    perturbedTransformedSineCoefficients = nominalTransformedSineCoefficients;

    perturbedTransformedCosineCoefficients( degreeOfBody2, std::abs( orderOfBody2 ) ) -= coefficientPerturbation;

    mutualPotentialField->computeCurrentEffectiveCoefficientsFromManualTransformedCoefficients(
                perturbedTransformedCosineCoefficients, perturbedTransformedSineCoefficients );
    double downPerturbedCosineCoefficientFromCosinePerturbation = mutualPotentialField->getEffectiveCosineCoefficient(
                degreeOfBody1, orderOfBody1, degreeOfBody2, orderOfBody2 );
    double downPerturbedSineCoefficientFromCosinePerturbation = mutualPotentialField->getEffectiveSineCoefficient(
                degreeOfBody1, orderOfBody1, degreeOfBody2, orderOfBody2 );



    perturbedTransformedCosineCoefficients = nominalTransformedCosineCoefficients;
    perturbedTransformedSineCoefficients = nominalTransformedSineCoefficients;

    perturbedTransformedSineCoefficients( degreeOfBody2, std::abs( orderOfBody2 ) ) += coefficientPerturbation;

    mutualPotentialField->computeCurrentEffectiveCoefficientsFromManualTransformedCoefficients(
                perturbedTransformedCosineCoefficients, perturbedTransformedSineCoefficients );
    double upPerturbedCosineCoefficientFromSinePerturbation = mutualPotentialField->getEffectiveCosineCoefficient(
                degreeOfBody1, orderOfBody1, degreeOfBody2, orderOfBody2 );
    double upPerturbedSineCoefficientFromSinePerturbation = mutualPotentialField->getEffectiveSineCoefficient(
                degreeOfBody1, orderOfBody1, degreeOfBody2, orderOfBody2 );



    perturbedTransformedCosineCoefficients = nominalTransformedCosineCoefficients;
    perturbedTransformedSineCoefficients = nominalTransformedSineCoefficients;

    perturbedTransformedSineCoefficients( degreeOfBody2, std::abs( orderOfBody2 ) ) -= coefficientPerturbation;

    mutualPotentialField->computeCurrentEffectiveCoefficientsFromManualTransformedCoefficients(
                perturbedTransformedCosineCoefficients, perturbedTransformedSineCoefficients );
    double downPerturbedCosineCoefficientFromSinePerturbation = mutualPotentialField->getEffectiveCosineCoefficient(
                degreeOfBody1, orderOfBody1, degreeOfBody2, orderOfBody2 );
    double downPerturbedSineCoefficientFromSinePerturbation = mutualPotentialField->getEffectiveSineCoefficient(
                degreeOfBody1, orderOfBody1, degreeOfBody2, orderOfBody2 );

    return ( ( Eigen::Matrix2d( )<<
               ( upPerturbedCosineCoefficientFromCosinePerturbation  - downPerturbedCosineCoefficientFromCosinePerturbation ),
               ( upPerturbedCosineCoefficientFromSinePerturbation  - downPerturbedCosineCoefficientFromSinePerturbation ),
               ( upPerturbedSineCoefficientFromCosinePerturbation  - downPerturbedSineCoefficientFromCosinePerturbation ),
               ( upPerturbedSineCoefficientFromSinePerturbation  - downPerturbedSineCoefficientFromSinePerturbation ) ).finished( ) ) /
            ( 2.0 * coefficientPerturbation );


}

BOOST_AUTO_TEST_SUITE( test_mutual_spherical_harmonic_gravity )

BOOST_AUTO_TEST_CASE( testMutualSphericalHarmonicGravity )
{
    // Load spice kernels.
    spice_interface::loadStandardSpiceKernels( { input_output::getSpiceKernelPath( ) + "de430_jup310_small.bsp" } );

    // Create list of bodies to create.
    std::vector< std::string > bodyNames;
    bodyNames.push_back( "Jupiter" );
    bodyNames.push_back( "Io" );
    bodyNames.push_back( "Europa" );
    bodyNames.push_back( "Sun" );

    // Specify initial time
    double initialTime = 1.0E7;
    double finalTime = 1.2E7;

    int expansionDegree = 5;

    //for( int isNormalized = 1; isNormalized < 1; isNormalized++ )
    int isNormalized = 1;
    {
        // Get body settings.
        std::map< std::string, boost::shared_ptr< BodySettings > > bodySettings =
                getDefaultBodySettings( bodyNames, initialTime, finalTime );
        bodySettings[ "Jupiter" ]->gravityFieldSettings = getDummyJovianSystemGravityField( "Jupiter", isNormalized );
        bodySettings[ "Io" ]->gravityFieldSettings = getDummyJovianSystemGravityField( "Io", isNormalized );
        bodySettings[ "Europa" ]->gravityFieldSettings = getDummyJovianSystemGravityField( "Europa", isNormalized );

        bodySettings[ "Jupiter" ]->rotationModelSettings = boost::make_shared<
                SimpleRotationModelSettings >( "ECLIPJ2000", "IAU_Jupiter_SIMPLIFIED", Eigen::Quaterniond(
                                                   Eigen::AngleAxisd( -1.4, Eigen::Vector3d::UnitZ( ) ) *
                                                   Eigen::AngleAxisd( 0.4, Eigen::Vector3d::UnitX( ) ) *
                                                   Eigen::AngleAxisd( 2.6, Eigen::Vector3d::UnitZ( ) ) ), 0.0, 0.0 );

        bodySettings[ "Io" ]->rotationModelSettings = boost::make_shared<
                SimpleRotationModelSettings >( "ECLIPJ2000", "IAU_Io_SIMPLIFIED", Eigen::Quaterniond(
                                                   Eigen::AngleAxisd( mathematical_constants::PI * 3.0 / 8.0, Eigen::Vector3d::UnitZ( ) ) *
                                                   Eigen::AngleAxisd( mathematical_constants::PI * 5.0 / 8.0, Eigen::Vector3d::UnitX( ) ) *
                                                   Eigen::AngleAxisd( mathematical_constants::PI * -2.5 / 8.0, Eigen::Vector3d::UnitZ( ) ) ), 0.0, 0.0 );
        bodySettings[ "Europa" ]->rotationModelSettings = boost::make_shared<
                SimpleRotationModelSettings >( "ECLIPJ2000", "IAU_Europa_SIMPLIFIED", Eigen::Quaterniond(
                                                   Eigen::AngleAxisd( mathematical_constants::PI * 3.0 / 8.0, Eigen::Vector3d::UnitZ( ) ) ), 0.0, 0.0 );

        // Create bodies needed in simulation
        NamedBodyMap bodyMap = createBodies( bodySettings );
        setGlobalFrameBodyEphemerides( bodyMap, "SSB", "ECLIPJ2000" );

        // Set current state and rotation of bodies.
        double currentTime = 1.1E7;
        bodyMap[ "Jupiter" ]->setCurrentRotationToLocalFrameFromEphemeris( currentTime );
        bodyMap[ "Jupiter" ]->setStateFromEphemeris( currentTime );
        bodyMap[ "Io" ]->setCurrentRotationToLocalFrameFromEphemeris( currentTime );
        bodyMap[ "Io" ]->setStateFromEphemeris( currentTime );
        bodyMap[ "Europa" ]->setCurrentRotationToLocalFrameFromEphemeris( currentTime );
        bodyMap[ "Europa" ]->setStateFromEphemeris( currentTime );

        // Retrieve gravity fields.
        boost::shared_ptr< SphericalHarmonicsGravityField > jupiterGravityField = boost::dynamic_pointer_cast< SphericalHarmonicsGravityField >(
                    ( bodyMap.at( "Jupiter" ) )->getGravityFieldModel( ) );
        boost::shared_ptr< SphericalHarmonicsGravityField > ioGravityField = boost::dynamic_pointer_cast< SphericalHarmonicsGravityField >(
                    ( bodyMap.at( "Io" ) )->getGravityFieldModel( ) );
        boost::shared_ptr< SphericalHarmonicsGravityField > europaGravityField = boost::dynamic_pointer_cast< SphericalHarmonicsGravityField >(
                    ( bodyMap.at( "Europa" ) )->getGravityFieldModel( ) );

        // Create central gravity acceleration (mu = Io + Jupiter)
        boost::shared_ptr< AccelerationSettings > centralGravitySettings = boost::make_shared< AccelerationSettings >(
                    central_gravity );
        boost::shared_ptr< CentralGravitationalAccelerationModel3d > centralGravity =
                boost::dynamic_pointer_cast< CentralGravitationalAccelerationModel3d >(
                    createAccelerationModel( bodyMap.at( "Io" ), bodyMap.at( "Jupiter" ), centralGravitySettings, "Io", "Jupiter",
                                             bodyMap.at( "Jupiter" ), "Jupiter" ) );

        // Calculate central gravity acceleration.
        centralGravity->updateMembers( );
        Eigen::Vector3d centralGravityAcceleration = centralGravity->getAcceleration( );

        // Create spherical harmonic gravity of Jupiter on Io, Jupiter-fixed (mu = Io + Jupiter)
        boost::shared_ptr< AccelerationSettings > sphericalHarmonicGravityOnIoFromJupiterSettings =
                boost::make_shared< SphericalHarmonicAccelerationSettings >( expansionDegree, expansionDegree );
        boost::shared_ptr< SphericalHarmonicsGravitationalAccelerationModel > sphericalHarmonicGravityOnIoFromJupiter =
                boost::dynamic_pointer_cast< SphericalHarmonicsGravitationalAccelerationModel >(
                    createAccelerationModel(  bodyMap.at( "Io" ), bodyMap.at( "Jupiter" ), sphericalHarmonicGravityOnIoFromJupiterSettings,
                                              "Io", "Jupiter" ) );

        // Calculate spherical harmonic gravity of Jupiter on Io.
        sphericalHarmonicGravityOnIoFromJupiter->updateMembers( );
        Eigen::Vector3d sphericalHarmonicGravityOnIoFromJupiterAcceleration = sphericalHarmonicGravityOnIoFromJupiter->getAcceleration( );
        double sphericalHarmonicPotentialAtIoFromJupiter = jupiterGravityField->getGravitationalPotential(
                    sphericalHarmonicGravityOnIoFromJupiter->getCurrentRelativePosition( ), expansionDegree, expansionDegree );

        // Create spherical harmonic gravity of Io on Jupiter, Io-fixed (mu = Io + Jupiter)
        boost::shared_ptr< AccelerationSettings > sphericalHarmonicGravityOnJupiterFromIoSettings =
                boost::make_shared< SphericalHarmonicAccelerationSettings >( expansionDegree, expansionDegree );
        boost::shared_ptr< SphericalHarmonicsGravitationalAccelerationModel > sphericalHarmonicGravityOnJupiterFromIo =
                boost::dynamic_pointer_cast< SphericalHarmonicsGravitationalAccelerationModel >(
                    createAccelerationModel( bodyMap.at( "Jupiter" ), bodyMap.at( "Io" ), sphericalHarmonicGravityOnJupiterFromIoSettings,
                                             "Jupiter", "Io" ) );

        // Calculate spherical harmonic gravity of Io on Jupiter.
        sphericalHarmonicGravityOnJupiterFromIo->updateMembers( );
        Eigen::Vector3d sphericalHarmonicGravityOnJupiterFromIoAcceleration = sphericalHarmonicGravityOnJupiterFromIo->getAcceleration( );
        double sphericalHarmonicPotentialAtJupiterFromIo = ioGravityField->getGravitationalPotential(
                    sphericalHarmonicGravityOnJupiterFromIo->getCurrentRelativePosition( ), expansionDegree, expansionDegree );

        // Create mutual spherical harmonic gravity between Io and Jupiter on Io, Jupiter fixed (mu = Io + Jupiter)
        boost::shared_ptr< AccelerationSettings > mutualDirectJupiterIoShGravitySettings =
                boost::make_shared< MutualSphericalHarmonicAccelerationSettings >( expansionDegree, expansionDegree, expansionDegree, expansionDegree );
        boost::shared_ptr< MutualSphericalHarmonicsGravitationalAccelerationModel > mutualDirectJupiterIoShGravity =
                boost::dynamic_pointer_cast< MutualSphericalHarmonicsGravitationalAccelerationModel >(
                    createAccelerationModel( bodyMap.at( "Io" ), bodyMap.at( "Jupiter" ), mutualDirectJupiterIoShGravitySettings,
                                             "Io", "Jupiter" ) );

        mutualDirectJupiterIoShGravity->updateMembers( );
        Eigen::Vector3d mutualSphericalHarmonicGravityOnIoFromJupiterAcceleration = mutualDirectJupiterIoShGravity->getAcceleration( );


        {
            // Create (through mutual extended body interface) central gravity acceleration (mu = Io + Jupiter)
            boost::shared_ptr< AccelerationSettings > extendedBodySettings = boost::make_shared< MutualExtendedBodySphericalHarmonicAccelerationSettings >(
                        expansionDegree, expansionDegree, 0, 0 );
            boost::shared_ptr< MutualExtendedBodySphericalHarmonicAcceleration > sh1ExtendedGravity =
                    boost::dynamic_pointer_cast< MutualExtendedBodySphericalHarmonicAcceleration >(
                        createAccelerationModel( bodyMap.at( "Io" ), bodyMap.at( "Jupiter" ), extendedBodySettings, "Io", "Jupiter" ) );
            sh1ExtendedGravity->updateMembers( );
            Eigen::Vector3d sh1ExtendedGravityAcceleration = sh1ExtendedGravity->getAcceleration( );

            Eigen::Vector3d precomputedAcceleration = ( ( -1.0 * spice_interface::getBodyGravitationalParameter( "Jupiter" ) ) /
                                                        5.959916033410404E012 * sphericalHarmonicGravityOnJupiterFromIoAcceleration );

            Eigen::Vector3d accelerationDifference = precomputedAcceleration - sh1ExtendedGravityAcceleration;

            for( unsigned int i = 0; i < 3; i++ )
            {
                BOOST_CHECK_SMALL( std::fabs( accelerationDifference( i ) ), 1.0E-15 );
            }
        }


        {
            // Create (through mutual extended body interface) central gravity acceleration (mu = Io + Jupiter)
            boost::shared_ptr< AccelerationSettings > extendedBodySettings = boost::make_shared< MutualExtendedBodySphericalHarmonicAccelerationSettings >(
                        0, 0, expansionDegree, expansionDegree );
            boost::shared_ptr< MutualExtendedBodySphericalHarmonicAcceleration > sh1ExtendedGravity =
                    boost::dynamic_pointer_cast< MutualExtendedBodySphericalHarmonicAcceleration >(
                        createAccelerationModel( bodyMap.at( "Io" ), bodyMap.at( "Jupiter" ), extendedBodySettings, "Io", "Jupiter" ) );
            sh1ExtendedGravity->updateMembers( );
            Eigen::Vector3d sh1ExtendedGravityAcceleration = sh1ExtendedGravity->getAcceleration( );;

            Eigen::Vector3d accelerationDifference = sphericalHarmonicGravityOnIoFromJupiterAcceleration - sh1ExtendedGravityAcceleration;

            for( unsigned int i = 0; i < 3; i++ )
            {
                BOOST_CHECK_SMALL( std::fabs( accelerationDifference( i ) ), 1.0E-15 );
            }
        }

        {
            // Create (through mutual extended body interface) central gravity acceleration (mu = Io + Jupiter)
            boost::shared_ptr< AccelerationSettings > extendedBodySettings = boost::make_shared< MutualExtendedBodySphericalHarmonicAccelerationSettings >(
                        0, 0, expansionDegree, expansionDegree );
            boost::shared_ptr< MutualExtendedBodySphericalHarmonicAcceleration > sh1ExtendedGravity =
                    boost::dynamic_pointer_cast< MutualExtendedBodySphericalHarmonicAcceleration >(
                        createAccelerationModel( bodyMap.at( "Jupiter" ), bodyMap.at( "Io" ), extendedBodySettings, "Jupiter", "Io" ) );
            sh1ExtendedGravity->updateMembers( );
            Eigen::Vector3d sh1ExtendedGravityAcceleration = sh1ExtendedGravity->getAcceleration( );

            Eigen::Vector3d accelerationDifference = sphericalHarmonicGravityOnJupiterFromIoAcceleration - sh1ExtendedGravityAcceleration;

            for( unsigned int i = 0; i < 3; i++ )
            {
                BOOST_CHECK_SMALL( std::fabs( accelerationDifference( i ) ), 1.0E-15 );
            }
        }

        {
            // Create (through mutual extended body interface) central gravity acceleration (mu = Io + Jupiter)
            boost::shared_ptr< AccelerationSettings > extendedBodySettings = boost::make_shared< MutualExtendedBodySphericalHarmonicAccelerationSettings >(
                        expansionDegree, expansionDegree, 0, 0 );
            boost::shared_ptr< MutualExtendedBodySphericalHarmonicAcceleration > sh1ExtendedGravity =
                    boost::dynamic_pointer_cast< MutualExtendedBodySphericalHarmonicAcceleration >(
                        createAccelerationModel( bodyMap.at( "Jupiter" ), bodyMap.at( "Io" ), extendedBodySettings, "Jupiter", "Io" ) );
            sh1ExtendedGravity->updateMembers( );
            Eigen::Vector3d sh1ExtendedGravityAcceleration = sh1ExtendedGravity->getAcceleration( );

            Eigen::Vector3d precomputedAcceleration = -1.0 * spice_interface::getBodyGravitationalParameter( "Jupiter" ) /
                    5.959916033410404E012 * sh1ExtendedGravityAcceleration;

            Eigen::Vector3d accelerationDifference = sphericalHarmonicGravityOnIoFromJupiterAcceleration - precomputedAcceleration;

            for( unsigned int i = 0; i < 3; i++ )
            {
                BOOST_CHECK_SMALL( std::fabs( accelerationDifference( i ) ), 1.0E-15 );
            }
        }

        {
            // Create (through mutual extended body interface) central gravity acceleration (mu = Io + Jupiter)
            boost::shared_ptr< AccelerationSettings > extendedBodySettings = boost::make_shared< MutualExtendedBodySphericalHarmonicAccelerationSettings >(
                        getExtendedSinglePointMassInteractions( expansionDegree, expansionDegree, expansionDegree, expansionDegree ) );
            boost::shared_ptr< MutualExtendedBodySphericalHarmonicAcceleration > sh1ExtendedGravity =
                    boost::dynamic_pointer_cast< MutualExtendedBodySphericalHarmonicAcceleration >(
                        createAccelerationModel( bodyMap.at( "Io" ), bodyMap.at( "Jupiter" ),  extendedBodySettings, "Io", "Jupiter" ) );
            sh1ExtendedGravity->updateMembers( );
            Eigen::Vector3d sh1ExtendedGravityAcceleration = sh1ExtendedGravity->getAcceleration( );

            Eigen::Vector3d accelerationDifference = mutualSphericalHarmonicGravityOnIoFromJupiterAcceleration - sh1ExtendedGravityAcceleration;

            for( unsigned int i = 0; i < 3; i++ )
            {
                BOOST_CHECK_SMALL( std::fabs( accelerationDifference( i ) ), 1.0E-15 );
            }
        }

        {
            // Create (through mutual extended body interface) central gravity acceleration (mu = Io + Jupiter)
            boost::shared_ptr< AccelerationSettings > extendedBodySettings = boost::make_shared< MutualExtendedBodySphericalHarmonicAccelerationSettings >(
                        getExtendedSinglePointMassInteractions( expansionDegree, expansionDegree, expansionDegree, expansionDegree ) );
            boost::shared_ptr< MutualExtendedBodySphericalHarmonicAcceleration > sh1ExtendedGravity =
                    boost::dynamic_pointer_cast< MutualExtendedBodySphericalHarmonicAcceleration >(
                        createAccelerationModel( bodyMap.at( "Jupiter" ), bodyMap.at( "Io" ),  extendedBodySettings, "Jupiter", "Io" ) );
            sh1ExtendedGravity->updateMembers( );
            Eigen::Vector3d sh1ExtendedGravityAcceleration = sh1ExtendedGravity->getAcceleration( );

            Eigen::Vector3d precomputedAcceleration = -1.0 * spice_interface::getBodyGravitationalParameter( "Jupiter" ) /
                    5.959916033410404E012 * sh1ExtendedGravityAcceleration;

            Eigen::Vector3d accelerationDifference = mutualSphericalHarmonicGravityOnIoFromJupiterAcceleration - precomputedAcceleration;

            for( unsigned int i = 0; i < 3; i++ )
            {
                BOOST_CHECK_SMALL( std::fabs( accelerationDifference( i ) ), 1.0E-15 );
            }

        }

        {
            std::cout<<"A"<<std::endl<<std::endl<<std::endl<<std::endl<<std::endl<<std::endl<<std::endl<<std::endl<<std::endl<<std::endl;
            // Create (through mutual extended body interface) central gravity acceleration (mu = Io + Jupiter)
            boost::shared_ptr< AccelerationSettings > extendedBodySettings = boost::make_shared< MutualExtendedBodySphericalHarmonicAccelerationSettings >(
                        getExtendedSinglePointMassInteractions( expansionDegree, expansionDegree, 0, 0 ) );
            boost::shared_ptr< TorqueSettings > extendedBodyTorqueSettings = boost::make_shared< MutualExtendedBodySphericalHarmonicTorqueSettings >(
                        2, 2, 2, 2 );
//                        getExtendedSinglePointMassInteractions( 2, 2, 2, 2 ) );
//            boost::shared_ptr< MutualExtendedBodySphericalHarmonicAcceleration > extendedGravityAcceleration =
//                    boost::dynamic_pointer_cast< MutualExtendedBodySphericalHarmonicAcceleration >(
//                        createAccelerationModel( bodyMap.at( "Jupiter" ), bodyMap.at( "Io" ),  extendedBodySettings, "Jupiter", "Io" ) );

            std::cout<<"B"<<std::endl;

            boost::shared_ptr< MutualExtendedBodySphericalHarmonicTorque > extendedGravityTorque =
                    boost::dynamic_pointer_cast< MutualExtendedBodySphericalHarmonicTorque >(
                        createTorqueModel( bodyMap.at( "Jupiter" ), bodyMap.at( "Io" ),  extendedBodyTorqueSettings, "Jupiter", "Io" ) );

            std::cout<<"C"<<std::endl;

            //extendedGravityAcceleration->updateMembers( );

            std::cout<<"D"<<std::endl;

            extendedGravityTorque->updateMembers( );

            std::cout<<"E"<<std::endl;


        }

    }
}
BOOST_AUTO_TEST_SUITE_END( )

}

}

