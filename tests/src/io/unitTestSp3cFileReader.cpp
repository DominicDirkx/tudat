#include <Eigen/Geometry>

#include <iostream>

#include "tudat/astro/ephemerides/itrsToGcrsRotationModel.h"
#include "tudat/astro/earth_orientation/earthOrientationCalculator.h"
#include "tudat/io/sp3cFileReader.h"

int main( )
{
    using namespace tudat;
    using namespace tudat::input_output;
    using namespace tudat::basic_astrodynamics;
    using namespace tudat::earth_orientation;
    using namespace tudat::ephemerides;

    std::string baseFolder = "/home/dominic/Software/tudat-bundle/tudat-bundle/tudat/tests/src/io/";
    boost::shared_ptr< SP3cFileContents > fileContents = readSp3cFile( baseFolder + "graA.2004-02-05.sp3" );
    std::cout<<fileContents->frameName<<" "<<fileContents->timeScale<<" "<<fileContents->vehiclesStates.size( )<<" "<<fileContents->vehiclesStates[ "L51" ].size( )<<std::endl;

    // Create rotation model
    std::shared_ptr< GcrsToItrsRotationModel > earthRotationModel =
            std::make_shared< GcrsToItrsRotationModel >(
                earth_orientation::createStandardEarthOrientationCalculator( ) );


    std::map< double, Eigen::VectorXd > itrsPositionHistory = fileContents->vehiclesStates[ "L09" ];
    std::map< double, Eigen::VectorXd > inertialPositionHistory;

    Eigen::Matrix3d rotationMatrix;
    Eigen::Matrix3d rotationRateMatrix;

    double currentTime, currentTaiTime;
    for( std::map< double, Eigen::VectorXd >::iterator it = itrsPositionHistory.begin( ); it != itrsPositionHistory.end( ); it++ )
    {
        currentTaiTime = it->first + 19;
        currentTime = earthRotationModel->getAnglesCalculator( )->getTerrestrialTimeScaleConverter( )->getCurrentTime(
                    tai_scale, tdb_scale, currentTaiTime, Eigen::Vector3d::Zero( ) );

        rotationMatrix = earthRotationModel->getRotationToBaseFrame( currentTime ).toRotationMatrix( );
        rotationRateMatrix =  earthRotationModel->getDerivativeOfRotationToBaseFrame( currentTime );

        inertialPositionHistory[ currentTime ] = Eigen::VectorXd::Zero( 3 );
        inertialPositionHistory[ currentTime ].segment( 0, 3 ) =
                rotationMatrix * it->second.segment( 0, 3 );
        std::cout<<std::setprecision( 16 );
        std::cout<<currentTime<<" "<<inertialPositionHistory[ currentTime ].transpose( )<<std::endl;
    }

    std::map< double, Eigen::VectorXd >::iterator iteratorMinus3 = inertialPositionHistory.begin( );
    std::map< double, Eigen::VectorXd >::iterator iteratorMinus2 = inertialPositionHistory.begin( );
    std::advance( iteratorMinus2, 1 );
    std::map< double, Eigen::VectorXd >::iterator iteratorMinus1 = inertialPositionHistory.begin( );
    std::advance( iteratorMinus1, 2 );
    std::map< double, Eigen::VectorXd >::iterator iterator = inertialPositionHistory.begin( );
    std::advance( iterator, 3 );
    std::map< double, Eigen::VectorXd >::iterator iteratorPlus1 = inertialPositionHistory.begin( );
    std::advance( iteratorPlus1, 4 );
    std::map< double, Eigen::VectorXd >::iterator iteratorPlus2 = inertialPositionHistory.begin( );
    std::advance( iteratorPlus2, 5 );
    std::map< double, Eigen::VectorXd >::iterator iteratorPlus3 = inertialPositionHistory.begin( );
    std::advance( iteratorPlus3, 6 );

    double deltaT = 5.0;

    while( iteratorPlus3 != inertialPositionHistory.end( ) )
    {
        Eigen::Vector3d stateMinus3 = iteratorMinus3->second;
        Eigen::Vector3d stateMinus2 = iteratorMinus2->second;
        Eigen::Vector3d stateMinus1 = iteratorMinus1->second;
        Eigen::Vector3d state = iterator->second;
        Eigen::Vector3d statePlus1 = iteratorPlus1->second;
        Eigen::Vector3d statePlus2 = iteratorPlus2->second;
        Eigen::Vector3d statePlus3 = iteratorPlus3->second;

        Eigen::Vector3d secondOrder = ( - stateMinus1 + statePlus1 ) / ( 2 * deltaT );
        Eigen::Vector3d fourthOrder =
                ( stateMinus2 / 12.0 - stateMinus1 * 2.0 / 3.0 +
                  statePlus1 * 2.0 / 3.0 - statePlus2 / 12.0 ) / ( deltaT );
        Eigen::Vector3d sixthOrder =
                ( -stateMinus3 / 60.0 +  stateMinus2 * 3.0 / 20.0 - stateMinus1 * 3.0 / 4.0 +
                  statePlus1 * 3.0 / 4.0 - statePlus2 * 3.0 / 20.0 + statePlus3 / 60.0 ) / ( deltaT );


//        coefficients[ order6 ][ -3 ] = -1.0 / 60.0;
//        coefficients[ order6 ][ -2 ] = 3.0 / 20.0;
//        coefficients[ order6 ][ -1 ] = -3.0 / 4.0;
//        coefficients[ order6 ][ 1 ] = 1.0 / 60.0;
//        coefficients[ order6 ][ 2 ] = -3.0 / 20.0;
//        coefficients[ order6 ][ 3 ] = 3.0 / 4.0;

        iteratorMinus3++;
        iteratorMinus2++;
        iteratorMinus1++;
        iterator++;
        iteratorPlus1++;
        iteratorPlus2++;
        iteratorPlus3++;
        std::cout<<iterator->first<<std::endl;
        std::cout<<secondOrder.transpose( )<<std::endl;
        std::cout<<fourthOrder.transpose( )<<std::endl;
        std::cout<<sixthOrder.transpose( )<<std::endl;
        std::cout<<( sixthOrder - fourthOrder ).transpose( )<<std::endl;

    }



}


