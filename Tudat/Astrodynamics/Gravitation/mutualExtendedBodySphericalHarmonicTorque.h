#ifndef MUTUALEXTENDEDBODYSPHERICALHARMONICTORQUE_CPP
#define MUTUALEXTENDEDBODYSPHERICALHARMONICTORQUE_CPP

#include <boost/function.hpp>
#include <boost/bind.hpp>
#include <boost/tuple/tuple.hpp>


#include <Eigen/Core>
#include <Eigen/Geometry>

#include <Tudat/Astrodynamics/BasicAstrodynamics/accelerationModel.h>
#include <Tudat/Mathematics/BasicMathematics/linearAlgebra.h>
#include <Tudat/Mathematics/BasicMathematics/legendrePolynomials.h>
#include <Tudat/Mathematics/BasicMathematics/wignerDMatrices.h>

#include "Tudat/Astrodynamics/BasicAstrodynamics/torqueModel.h"
#include "Tudat/Astrodynamics/Gravitation/mutualExtendedBodySphericalHarmonicAcceleration.h"

namespace tudat
{

namespace gravitation
{

class MutualExtendedBodySphericalHarmonicTorque: public basic_astrodynamics::TorqueModel
{

public:

    MutualExtendedBodySphericalHarmonicTorque(
            const boost::shared_ptr< MutualExtendedBodySphericalHarmonicAcceleration > accelerationBetweenBodies,
            const boost::function< double( ) > gravitationalParameterFunctionOfBodyUndergoingTorque ):
        accelerationBetweenBodies_( accelerationBetweenBodies ),
        gravitationalParameterFunctionOfBodyUndergoingTorque_( gravitationalParameterFunctionOfBodyUndergoingTorque )
    {
        effectiveCoefficientCalculator_ = accelerationBetweenBodies->getEffectiveMutualPotentialField( );
        effectiveCoefficientCalculator_->getTransformationCache( )->getWignerDMatricesCache( )->setComputeAngularMomentumOperators(
                    true );
    }

    void updateMembers( const double currentTime = TUDAT_NAN );

    Eigen::Vector3d getTorque( )
    {
        return currentTorque_;
    }

    Eigen::Vector3d getTorqueOnBodyExertingTorque( )
    {
        return - ( accelerationBetweenBodies_->getCurrentRotationFromBody2ToBody1( ).inverse( ) *
                -currentAngularMomentumOpertorOfMutualPotential_ );
    }

    boost::shared_ptr< MutualExtendedBodySphericalHarmonicAcceleration > getAccelerationBetweenBodies( )
    {
        return accelerationBetweenBodies_;
    }

    Eigen::Vector3d getCurrentAngularMomentumOpertorOfMutualPotential( )
    {
        return currentAngularMomentumOpertorOfMutualPotential_;
    }


private:

    Eigen::Vector3d currentTorque_;

    //current torque acting on body exerting acceleration in frame fixed to body undergoing acceleration_
    Eigen::Vector3d currentAngularMomentumOpertorOfMutualPotential_;

    boost::shared_ptr< MutualExtendedBodySphericalHarmonicAcceleration > accelerationBetweenBodies_;

    boost::function< double( ) > gravitationalParameterFunctionOfBodyUndergoingTorque_;

    boost::shared_ptr< gravitation::EffectiveMutualSphericalHarmonicsField > effectiveCoefficientCalculator_;
};

}

}

#endif // MUTUALEXTENDEDBODYSPHERICALHARMONICTORQUE_CPP
