#define BOOST_TEST_MAIN

#include <boost/test/unit_test.hpp>

#include "Tudat/External/SpiceInterface/spiceInterface.h"
#include "Tudat/InputOutput/basicInputOutput.h"

#include "Tudat/External/SpiceInterface/spiceRotationalEphemeris.h"
#include "Tudat/Astrodynamics/Ephemerides/fullPlanetaryRotationModel.h"
#include "Tudat/SimulationSetup/EnvironmentSetup/defaultBodies.h"
#include "Tudat/SimulationSetup/EnvironmentSetup/createRotationModel.h"



namespace tudat
{
namespace unit_tests
{

using namespace tudat::ephemerides;
using namespace tudat::simulation_setup;

BOOST_AUTO_TEST_SUITE( test_planetary_rotation_model )

BOOST_AUTO_TEST_CASE( testPlanetaryRotationModel )
{
    double initialTime = 0.0;
    double finalTime = 1.0E8;

    spice_interface::loadStandardSpiceKernels( );

    NamedBodyMap bodyMap;
    bodyMap[ "Mars" ] = std::make_shared< Body >( );

    std::shared_ptr< RotationModelSettings > defaultMarsRotationSettings = getHighAccuracyMarsRotationModel(
            initialTime, finalTime );

    std::shared_ptr< RotationalEphemeris > marsRotationModel = createRotationModel(
                defaultMarsRotationSettings, "Mars" );

    std::shared_ptr< PlanetaryRotationModelSettings > modifiedMarsRotationModel = std::dynamic_pointer_cast< PlanetaryRotationModelSettings >(
                defaultMarsRotationSettings );

    modifiedMarsRotationModel->setPeriodTermsToZero( );

    std::shared_ptr< RotationalEphemeris > marsSimplifiedRotationModel = createRotationModel(
                modifiedMarsRotationModel, "Mars" );

    std::shared_ptr< RotationalEphemeris > spiceMarsRotationModel = std::make_shared< SpiceRotationalEphemeris >(
                "ECLIPJ2000", "IAU_Mars" );

    double timeStep = 3600.0;
    double currentTime = initialTime + timeStep;

    std::map< double, Eigen::MatrixXd > rotationDifferenceMap;

    while( currentTime < finalTime - timeStep )
    {
        rotationDifferenceMap[ currentTime ] =  ( marsRotationModel->getRotationToTargetFrame( currentTime ) ).toRotationMatrix( ) -
                ( spiceMarsRotationModel->getRotationToTargetFrame( currentTime ) ).toRotationMatrix( );

        std::cout<<rotationDifferenceMap[ currentTime ]<<std::endl<<std::endl;
        rotationDifferenceMap[ currentTime ] =  ( marsSimplifiedRotationModel->getRotationToTargetFrame( currentTime ) ).toRotationMatrix( ) -
                ( spiceMarsRotationModel->getRotationToTargetFrame( currentTime ) ).toRotationMatrix( );
        std::cout<<rotationDifferenceMap[ currentTime ]<<std::endl<<std::endl<<std::endl<<std::endl;

        for( unsigned int i = 0; i < 3; i++ )
        {
            for( unsigned j = 0; j < 3; j++ )
            {
                BOOST_CHECK_SMALL( rotationDifferenceMap[ currentTime ]( i, j ), 1.0E-5 );
            }
        }
        currentTime += timeStep;
    }
}

BOOST_AUTO_TEST_SUITE_END( )

}

}
