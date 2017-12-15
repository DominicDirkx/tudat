
#include "Tudat/Astrodynamics/Gravitation/mutualExtendedBodySphericalHarmonicTorque.h"

namespace tudat
{

namespace gravitation
{

void MutualExtendedBodySphericalHarmonicTorque::updateMembers( const double currentTime )
{

    accelerationBetweenBodies_->updateMembers( currentTime );

    currentAngularMomentumOpertorOfMutualPotential_ =
            gravitationalParameterFunctionOfBodyUndergoingTorque_( ) / physical_constants::GRAVITATIONAL_CONSTANT *
            effectiveCoefficientCalculator_->getAngularMomentumOpertorOfGravitationalPotential(
                accelerationBetweenBodies_->getCurrentRelativePosition( ),
                accelerationBetweenBodies_->getSphericalHarmonicsCache( ) );

    std::cout<<"Torque cross-product extended: "<<std::endl<<
               accelerationBetweenBodies_->getCurrentRelativePosition( ).norm( )<<" "<<
               accelerationBetweenBodies_->getCurrentRelativePosition( ).transpose( )<<std::endl<<
               accelerationBetweenBodies_->getAccelerationInBodyFixedFrame( ).norm( )<<" "<<
               accelerationBetweenBodies_->getAccelerationInBodyFixedFrame( ).transpose( )<<std::endl;

    currentTorque_ = currentAngularMomentumOpertorOfMutualPotential_ -
            accelerationBetweenBodies_->getCurrentBodyFixedRelativePosition( ).cross(
                gravitationalParameterFunctionOfBodyUndergoingTorque_( ) / physical_constants::GRAVITATIONAL_CONSTANT *
                accelerationBetweenBodies_->getAccelerationInBodyFixedFrame( ) );

}

}

}
