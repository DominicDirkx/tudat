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

#include <cmath>

#include "Tudat/Astrodynamics/Gravitation/centralGravityModel.h"

namespace tudat
{
namespace gravitation
{

 //! Compute gravitational acceleration.
Eigen::Vector3d computeGravitationalAcceleration(
        const double universalGravitationalConstant,
        const Eigen::Vector3d& positionOfBodySubjectToAcceleration,
        const double massOfBodyExertingAcceleration,
        const Eigen::Vector3d& positionOfBodyExertingAcceleration )
{
    return computeGravitationalAcceleration(
                positionOfBodySubjectToAcceleration,
                universalGravitationalConstant * massOfBodyExertingAcceleration,
                positionOfBodyExertingAcceleration );
}

void computeGravitationalAccelerationByReference(
        const Eigen::Vector3d& relativePositionOfBodySubjectToAcceleration,
        const double gravitationalParameterOfBodyExertingAcceleration,
        Eigen::Vector3d& acceleration )
{

    double distance = relativePositionOfBodySubjectToAcceleration.norm( );
    acceleration = -( gravitationalParameterOfBodyExertingAcceleration / ( distance * distance * distance ) )
            * ( relativePositionOfBodySubjectToAcceleration );
}

Eigen::Vector3d computeGravitationalAcceleration(
        const Eigen::Vector3d& relativePositionOfBodySubjectToAcceleration,
        const double gravitationalParameterOfBodyExertingAcceleration )
{

    double distance = relativePositionOfBodySubjectToAcceleration.norm( );
    return -gravitationalParameterOfBodyExertingAcceleration
            * ( relativePositionOfBodySubjectToAcceleration )
            / ( distance * distance * distance );
}

//! Compute gravitational acceleration.
Eigen::Vector3d computeGravitationalAcceleration(
        const Eigen::Vector3d& positionOfBodySubjectToAcceleration,
        const double gravitationalParameterOfBodyExertingAcceleration,
        const Eigen::Vector3d& positionOfBodyExertingAcceleration )
{
    return computeGravitationalAcceleration(
                positionOfBodySubjectToAcceleration - positionOfBodyExertingAcceleration,
                gravitationalParameterOfBodyExertingAcceleration );
}

//! Compute gravitational force.
Eigen::Vector3d computeGravitationalForce(
        const double universalGravitationalParameter,
        const double massOfBodySubjectToForce,
        const Eigen::Vector3d& positionOfBodySubjectToForce,
        const double massOfBodyExertingForce,
        const Eigen::Vector3d& positionOfBodyExertingForce )
{
    return massOfBodySubjectToForce * computeGravitationalAcceleration(
                universalGravitationalParameter, positionOfBodySubjectToForce,
                massOfBodyExertingForce, positionOfBodyExertingForce );
}

//! Compute gravitational force.
Eigen::Vector3d computeGravitationalForce(
        const double massOfBodySubjectToForce,
        const Eigen::Vector3d& positionOfBodySubjectToForce,
        const double gravitationalParameterOfBodyExertingForce,
        const Eigen::Vector3d& positionOfBodyExertingForce )
{
    return massOfBodySubjectToForce * computeGravitationalAcceleration(
                positionOfBodySubjectToForce, gravitationalParameterOfBodyExertingForce,
                positionOfBodyExertingForce );
}

} // namespace gravitation
} // namespace tudat
