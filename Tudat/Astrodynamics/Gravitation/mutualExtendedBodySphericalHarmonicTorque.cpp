
#include "Tudat/Astrodynamics/Gravitation/mutualExtendedBodySphericalHarmonicTorque.h"

namespace tudat
{

namespace gravitation
{

void MutualExtendedBodySphericalHarmonicTorque::updateMembers( const double currentTime )
{

    return computeMutualForcePotential(
                bodyFixedPosition, gravitationalParameterFunction_( ), equatorialRadiusOfBody1_, equatorialRadiusOfBody2_,
                maximumDegree1_, maximumDegree2_,
                boost::bind(
                    &EffectiveMutualSphericalHarmonicsField::getEffectiveCosineCoefficient,
                    this, _1, _2, _3, _4 ),
                boost::bind(
                    &EffectiveMutualSphericalHarmonicsField::getEffectiveSineCoefficient,
                    this, _1, _2, _3, _4 ), coefficientCombinationsToUse_,
                sphericalHarmonicsCache );


//    if( acceleratedBodyIsBody1_ )
//    {
//        currentTorque_ = -acceleratedBodyIsBody1_ - ( accelerationBetweenBodies_->getCurrentRelativePosition( ) ).cross(
//                    accelerationBetweenBodies_->getAcceleration( ) );
//    }
//    else
//    {
//        currentTorque_ = acceleratedBodyIsBody1_;
//    }

}

}

}
