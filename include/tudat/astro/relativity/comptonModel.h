//
// Created by Michael Plumaris on 24/01/2022.
//

#ifndef TUDATBUNDLE_COMPTONMODEL_H
#define TUDATBUNDLE_COMPTONMODEL_H

#include <boost/lambda/lambda.hpp>
#include <memory>
#include <boost/make_shared.hpp>

#include <Eigen/Core>

#include "tudat/astro/basic_astro/accelerationModel.h"
#include "tudat/astro/gravitation/sphericalHarmonicsGravityModelBase.h"
#include "tudat/astro/relativity/metric.h"

//using namespace gravitation;

namespace tudat
{
namespace relativity
{


// Michael
Eigen::Vector3d computeComptonRelativisticAcceleration(
        const Eigen::Vector3d& vectorToAcceleratedBody,
        const double gravitationalParameter,
        const double ComptonWavelength = relativity::comptonWavelength );

/*
Eigen::Vector3d computeGravitationalForce(
const double universalGravitationalParameter,
const double massOfBodySubjectToForce,
const Eigen::Vector3d& positionOfBodySubjectToForce,
const double massOfBodyExertingForce,
const Eigen::Vector3d& positionOfBodyExertingForce = Eigen::Vector3d::Zero( )  );

*/

class ComptonRelativisticAcceleration : public basic_astrodynamics::AccelerationModel3d
{
private:

    //typedef gravitation::SphericalHarmonicsGravitationalAccelerationModelBase< Eigen::Vector3d > Base;
    typedef std::function< Eigen::Vector3d( ) > Vector3dReturningFunction;

public:

    // Ensure that correctly aligned pointers are generated (Eigen, 2013).
    EIGEN_MAKE_ALIGNED_OPERATOR_NEW


    ComptonRelativisticAcceleration(
            const Vector3dReturningFunction sourcePositionFunction,
            const Vector3dReturningFunction acceleratedBodyPositionFunction,
            const std::function< double( ) > gravitationalParameterFunction
            //std::function< double( ) > comptonWavelengthFunction = [ ]( ){ return std::numeric_limits<double>::infinity(); }
            )
        : sourcePositionFunction_( sourcePositionFunction ),
          acceleratedBodyPositionFunction_( acceleratedBodyPositionFunction ),
          gravitationalParameterFunction_( gravitationalParameterFunction )
          //comptonWavelengthFunction_ (comptonWavelength)
    {
        this->updateMembers();
    }


    void updateMembers(const double currentTime = TUDAT_NAN) {
        if (!(this->currentTime_ == currentTime)) {
            //this->updateBaseMembers();
            currentGravitationalParameter_ = gravitationalParameterFunction_( );
            vectorToAcceleratedBody_ = ( acceleratedBodyPositionFunction_( )
                                       - sourcePositionFunction_( ) );

            this->currentAcceleration_ = computeComptonRelativisticAcceleration(
                    this->vectorToAcceleratedBody_,
                    this->currentGravitationalParameter_
                    //this->comptonWavelength_
                    );
        }
    }

//    std::function< double( ) > getcomptonWavelengthFunction_( )
//    {
//        return comptonWavelengthFunction_;
//    }


protected:

private:

    //double comptonWavelength_;

    const Vector3dReturningFunction sourcePositionFunction_;

    const Vector3dReturningFunction acceleratedBodyPositionFunction_;

    const std::function< double( ) > gravitationalParameterFunction_;

    Eigen::Vector3d vectorToAcceleratedBody_;

    double currentGravitationalParameter_;
    //std::function< double( ) > comptonWavelengthFunction_;

};


//! Typedef for shared-pointer to CentralGravitationalAccelerationModel3d.
typedef std::shared_ptr< ComptonRelativisticAcceleration > ComptonRelativisticAccelerationPointer;

}

} // namespace tudat


#endif //TUDATBUNDLE_COMPTONMODEL_H
