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
            const boost::shared_ptr< SphericalHarmonicsGravitationalAccelerationModel > sphericalHarmonicAcceleration ){ }

    //! Get gravitational torque.
    /*!
     * Returns the gravitational torque. All data required for the computation is taken
     * from member variables, which are set to their latest values by the last call of the
     * updateMembers function.
     * \return Gravitational torque.
     * \sa updateMembers().
     */
    Eigen::Vector3d getTorque( )
    {
        return currentTorque_;
    }

    //! Update member variables used by the gravitational torque model.
    /*!
     * Updates member variables used by the gravitational accfeleration model.
     * Function pointers to retrieve the current values of quantities from which the
     * torque is to be calculated are set by constructor. This function calls
     * them to update the associated variables to their current state.
     * \param currentTime Time at which torque model is to be updated.
     */
    void updateMembers( const double currentTime )
    {
        sphericalHarmonicAcceleration_->updateMembers( currentTime );

        currentTorque_ =
    }

protected:

private:


    boost::shared_ptr< SphericalHarmonicsGravitationalAccelerationModel > sphericalHarmonicAcceleration_;

    Eigen::Vector3d currentTorque_;
};

}

}

#endif // TUDAT_SPHERICALHARMONICGRAVITATIONALTORQUE_H
