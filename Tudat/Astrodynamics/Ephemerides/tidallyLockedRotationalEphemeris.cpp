#include <iostream>
#include <iomanip>

#include "Tudat/Astrodynamics/BasicAstrodynamics/physicalConstants.h"
#include "Tudat/Astrodynamics/BasicAstrodynamics/orbitalElementConversions.h"
#include "Tudat/Astrodynamics/Ephemerides/tidallyLockedRotationalEphemeris.h"


namespace tudat
{

namespace ephemerides
{

Eigen::Vector6d getStateFromSelectedStateFunction(
        const double currentTime,
        const bool useFirstFunction,
        const std::function< Eigen::Vector6d( const double ) > stateFunction1,
        const std::function< Eigen::Vector6d( const double ) > stateFunction2 )
{
    return ( useFirstFunction ) ? ( stateFunction1( currentTime ) ) : ( stateFunction2( currentTime ) );
}


double evaluateRightAscensionDeclinationFunction(
        const double nominalValue,
        const double secularRate,
        const std::map< double, std::pair< double, double > >& periodicTerms,
        const double currentTime,
        const bool useCosines )
{
    double angleValue = nominalValue;
    double centuriesSinceEpoch = currentTime / ( 100.0 * 86400.0 * 365.25 );
    angleValue += secularRate * centuriesSinceEpoch;

    for( std::map< double, std::pair< double, double > >::const_iterator termIterator = periodicTerms.begin( );
         termIterator != periodicTerms.end( ); termIterator++ )
    {
        if( useCosines )
        {
            angleValue += termIterator->first * std::cos( termIterator->second.first *  centuriesSinceEpoch + termIterator->second.second );
        }
        else
        {
            angleValue += termIterator->first * std::sin( termIterator->second.first *  centuriesSinceEpoch + termIterator->second.second );
        }
    }

    return angleValue;
}

Eigen::Quaterniond TidallyLockedRotationalEphemeris::getRotationToBaseFrame(
        const double currentRightAscension,
        const double currentDeclination,
        const Eigen::Vector6d relativeState )
{
    Eigen::Quaterniond secondIntermediateRotationToTargetFrame = getSecondIntermediateRotationToTargetFrame( currentRightAscension, currentDeclination );

    throw std::runtime_error( "Need to correct line of code ");
    double currentRotationAngle;// = getFullyLockedRotationAngle(
//                secondIntermediateRotationToTargetFrame, relativeState );

    return secondIntermediateRotationToTargetFrame.inverse( ) *
            reference_frames::getRotatingPlanetocentricToInertialFrameTransformationQuaternion(
                currentRotationAngle );
}


Eigen::Quaterniond TidallyLockedRotationalEphemeris::getRotationToBaseFrame( const double currentTime )
{
    double currentRightAscension = rightAscensionFunction_( currentTime );
    double currentDeclination = declinationFunction_( currentTime );

    Eigen::Vector6d relativeState = relativeStateFunction_( currentTime, isBodyInPropagation_ );

    return getRotationToBaseFrame( currentRightAscension, currentDeclination, relativeState );
}

Eigen::Matrix3d TidallyLockedRotationalEphemeris::getDerivativeOfRotationToBaseFrame( const double currentTime )
{
    double currentRightAscension = rightAscensionFunction_( currentTime );
    double currentDeclination = declinationFunction_( currentTime );

    Eigen::Vector6d currentRelativeState = ( relativeStateFunction_( currentTime, isBodyInPropagation_ ) );

    Eigen::Quaterniond intermediateRotationToTargetFrame = getSecondIntermediateRotationToTargetFrame( currentRightAscension, currentDeclination );

    Eigen::Vector2d currentRelativeCentralBodyPosition = ( intermediateRotationToTargetFrame * ( -currentRelativeState.segment( 0, 3 ) ) ).segment( 0, 2 );
    Eigen::Vector2d currentRelativeCentralBodyVelocity  = ( intermediateRotationToTargetFrame * ( -currentRelativeState.segment( 3, 3 ) ) ).segment( 0, 2 );

    double currentRotationAngle =
            std::atan2( currentRelativeCentralBodyPosition( 1 ), currentRelativeCentralBodyPosition( 0 ) );
    double currentRotationRate = currentRelativeCentralBodyPosition( 0 ) * currentRelativeCentralBodyVelocity( 1 ) -
            currentRelativeCentralBodyPosition( 1  ) * currentRelativeCentralBodyVelocity( 0 );
    currentRotationRate /= ( currentRelativeCentralBodyPosition.norm( ) * currentRelativeCentralBodyPosition.norm( ) );

    Eigen::Matrix3d auxiliaryMatrix;
    auxiliaryMatrix<< 0.0, 1.0, 0.0, -1.0, 0.0, 0.0, 0.0, 0.0, 0.0;

    return -currentRotationRate * Eigen::Matrix3d( intermediateToBaseFrameRotation_ ) *
            Eigen::Matrix3d( reference_frames::getRotatingPlanetocentricToInertialFrameTransformationQuaternion(
                            currentDeclination, currentRightAscension, 0.0 ) ) * auxiliaryMatrix * tudat::reference_frames::
                                   getInertialToPlanetocentricFrameTransformationQuaternion( -currentRotationAngle );
}

}

}

