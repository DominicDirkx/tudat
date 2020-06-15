#include <iostream>
#include <iomanip>

#include "Tudat/Astrodynamics/BasicAstrodynamics/physicalConstants.h"
#include "Tudat/Astrodynamics/BasicAstrodynamics/orbitalElementConversions.h"
#include "Tudat/Astrodynamics/Ephemerides/synchronousRotationalEphemeris.h"


namespace tudat
{

namespace ephemerides
{

//! Calculate rotation quaternion from target frame to base frame.
Eigen::Matrix3d SynchronousRotationalEphemeris::getFullyLockedRotationToBaseFrame( const double currentTime )
{
    // Get rotation to RSW frame
    Eigen::Vector6d relativeState = relativeStateFunction_( currentTime, isBodyInPropagation_ );
    Eigen::Matrix3d rotationToBaseFrame = reference_frames::getInertialToRswSatelliteCenteredFrameRotationMatrix(
                relativeState ).transpose( );
    rotationToBaseFrame.block( 0, 0, 3, 2 ) *= -1.0;

    return rotationToBaseFrame;
}

Eigen::Matrix3d SynchronousRotationalEphemeris::getLibrationRotation( const double currentTime )
{
    if( isLibrationOn_ )
    {
        double librationAngle = librationAngleFunction_( currentTime, longitudeLibrationAmlpitude_ );
        return Eigen::Matrix3d( Eigen::AngleAxisd( -librationAngle, Eigen::Vector3d::UnitZ( ) ) );
    }
    else
    {
        return Eigen::Matrix3d::Identity( );
    }
}

//! Calculate rotation quaternion from target frame to base frame.
Eigen::Quaterniond SynchronousRotationalEphemeris::getRotationToBaseFrame( const double currentTime )
{

    Eigen::Matrix3d rotationToBaseFrame = getFullyLockedRotationToBaseFrame( currentTime );
    Eigen::Matrix3d deviationFromSynchronicity = Eigen::Matrix3d::Identity( );
    if( isLibrationOn_ )
    {
        deviationFromSynchronicity = getLibrationRotation( currentTime );
    }
    // Flip sign of x and y axes

    rotationToBaseFrame = rotationToBaseFrame * deviationFromSynchronicity;

    return Eigen::Quaterniond( rotationToBaseFrame );
}

//! Function to calculate the derivative of the rotation matrix from target frame to base frame.
Eigen::Matrix3d SynchronousRotationalEphemeris::getDerivativeOfRotationToBaseFrame( const double currentTime )
{
    if( !warningPrinted_ )
    {
        std::cerr<<"Warning, time-derivative of synchronous rotation matrix not yet implemented (using zero matrix)"<<std::endl;
        warningPrinted_ = true;
    }
    return Eigen::Matrix3d::Zero( );
}

}

}

