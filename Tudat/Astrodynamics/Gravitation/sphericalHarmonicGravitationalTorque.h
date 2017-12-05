/*    Copyright (c) 2010-2017, Delft University of Technology
 *    All rigths reserved
 *
 *    This file is part of the Tudat. Redistribution and use in source and
 *    binary forms, with or without modification, are permitted exclusively
 *    under the terms of the Modified BSD license. You should have received
 *    a copy of the license with this file. If not, please or visit:
 *    http://tudat.tudelft.nl/LICENSE.
 */

#ifndef TUDAT_SPHERICALHARMONICGRAVITATIONALTORQUE_H
#define TUDAT_SPHERICALHARMONICGRAVITATIONALTORQUE_H


#include <boost/function.hpp>
#include <boost/lambda/lambda.hpp>

#include <Eigen/Geometry>

#include "Tudat/Astrodynamics/BasicAstrodynamics/physicalConstants.h"
#include "Tudat/Astrodynamics/BasicAstrodynamics/torqueModel.h"
#include "Tudat/Astrodynamics/Gravitation/sphericalHarmonicsGravityModel.h"

namespace tudat
{

namespace gravitation
{


class SphericalHarmonicGravitationalTorqueModel: public basic_astrodynamics::TorqueModel
{
public:

    SphericalHarmonicGravitationalTorqueModel(
            const boost::shared_ptr< SphericalHarmonicsGravitationalAccelerationModel > sphericalHarmonicAcceleration,
            const boost::function< Eigen::Quaterniond( ) > rotationToBodyUndergoingTorque,
            const boost::function< double( ) > perturberMassFunction ):
        sphericalHarmonicAcceleration_( sphericalHarmonicAcceleration ),
        rotationToBodyUndergoingTorque_( rotationToBodyUndergoingTorque ),
        perturberMassFunction_( perturberMassFunction ){ }

    Eigen::Vector3d getTorque( )
    {
        return currentTorque_;
    }

    void updateMembers( const double currentTime )
    {
        sphericalHarmonicAcceleration_->updateMembers( currentTime );

        currentTorque_ = perturberMassFunction_( ) *
                ( ( sphericalHarmonicAcceleration_->getCurrentRelativePosition( ) ).cross(
                      sphericalHarmonicAcceleration_->getAccelerationInBodyFixedFrame( ) ) );
    }

protected:

private:


    boost::shared_ptr< SphericalHarmonicsGravitationalAccelerationModel > sphericalHarmonicAcceleration_;

    boost::function< Eigen::Quaterniond( ) > rotationToBodyUndergoingTorque_;

    boost::function< double( ) > perturberMassFunction_;


    Eigen::Vector3d currentTorque_;
};

}

}

#endif // TUDAT_SPHERICALHARMONICGRAVITATIONALTORQUE_H
