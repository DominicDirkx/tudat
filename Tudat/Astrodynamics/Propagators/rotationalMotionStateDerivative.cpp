/*    Copyright (c) 2010-2018, Delft University of Technology
 *    All rigths reserved
 *
 *    This file is part of the Tudat. Redistribution and use in source and
 *    binary forms, with or without modification, are permitted exclusively
 *    under the terms of the Modified BSD license. You should have received
 *    a copy of the license with this file. If not, please or visit:
 *    http://tudat.tudelft.nl/LICENSE.
 *
 */

#include "Tudat/Astrodynamics/Gravitation/sphericalHarmonicsGravityField.h"
#include "Tudat/Astrodynamics/Propagators/rotationalMotionStateDerivative.h"

namespace tudat
{

namespace propagators
{

//! Function to evaluated the classical rotational equations of motion (Euler equations)
Eigen::Vector3d evaluateRotationalEquationsOfMotion(
        const Eigen::Matrix3d& inertiaTensor, const Eigen::Vector3d& totalTorque,
        const Eigen::Vector3d& angularVelocityVector, const Eigen::Matrix3d& inertiaTensorTimeDerivative )
{
    std::cout<<"Derivative: "<<totalTorque.transpose( )<<std::endl<<inertiaTensor<<std::endl;
    return inertiaTensor.inverse( ) * ( totalTorque );
}

//! Function to obtain the matrix by which a quaternion vector is to be pre-multiplied to obtain this quaternion's time-derivative=
Eigen::Matrix4d getQuaterionToQuaternionRateMatrix( const Eigen::Vector3d& angularVelocityVectorInBodyFixedFrame )
{
    Eigen::Matrix4d conversionMatrix = Eigen::Matrix4d::Zero( );
    conversionMatrix( 1, 0 ) = angularVelocityVectorInBodyFixedFrame( 0 );
    conversionMatrix( 2, 0 ) = angularVelocityVectorInBodyFixedFrame( 1 );
    conversionMatrix( 3, 0 ) = angularVelocityVectorInBodyFixedFrame( 2 );

    conversionMatrix( 2, 1 ) = -angularVelocityVectorInBodyFixedFrame( 2 );
    conversionMatrix( 3, 1 ) = angularVelocityVectorInBodyFixedFrame( 1 );

    conversionMatrix( 3, 2 ) = -angularVelocityVectorInBodyFixedFrame( 0 );

    conversionMatrix( 0, 1 ) = -angularVelocityVectorInBodyFixedFrame( 0 );
    conversionMatrix( 0, 2 ) = -angularVelocityVectorInBodyFixedFrame( 1 );
    conversionMatrix( 0, 3 ) = -angularVelocityVectorInBodyFixedFrame( 2 );

    conversionMatrix( 1, 2 ) = angularVelocityVectorInBodyFixedFrame( 2 );
    conversionMatrix( 1, 3 ) = -angularVelocityVectorInBodyFixedFrame( 1 );

    conversionMatrix( 2, 3 ) = angularVelocityVectorInBodyFixedFrame( 0 );

    return 0.5 * ( conversionMatrix );
}

Eigen::Matrix< double, 4, 3 > getAngularVelocityToQuaternionRateMatrix( const Eigen::Vector4d& quaternionVector )
{
    Eigen::Matrix< double, 4, 3 > conversionMatrix = Eigen::Matrix< double, 4, 3 >::Zero( );

    conversionMatrix( 0, 0 ) = -quaternionVector( 1 );
    conversionMatrix( 0, 1 ) = -quaternionVector( 2 );
    conversionMatrix( 0, 2 ) = -quaternionVector( 3 );

    conversionMatrix( 1, 0 ) = quaternionVector( 0 );
    conversionMatrix( 1, 1 ) = -quaternionVector( 3 );
    conversionMatrix( 1, 2 ) = quaternionVector( 2 );

    conversionMatrix( 2, 0 ) = quaternionVector( 3 );
    conversionMatrix( 2, 1 ) = quaternionVector( 0 );
    conversionMatrix( 2, 2 ) = -quaternionVector( 1 );

    conversionMatrix( 3, 0 ) = -quaternionVector( 2 );
    conversionMatrix( 3, 1 ) = quaternionVector( 1 );
    conversionMatrix( 3, 2 ) = quaternionVector( 0 );

    return 0.5 * conversionMatrix;

}

//! Function to obtain the time derivative of a quaternion (in vector representation) of body-fixed to inertial frame
Eigen::Vector4d calculateQuaternionDerivative(
        const Eigen::Vector4d& currentQuaternionToBaseFrame, const Eigen::Vector3d& angularVelocityVectorInBodyFixedFrame )
{
    return getQuaterionToQuaternionRateMatrix( angularVelocityVectorInBodyFixedFrame ) * currentQuaternionToBaseFrame;
}

}

}
