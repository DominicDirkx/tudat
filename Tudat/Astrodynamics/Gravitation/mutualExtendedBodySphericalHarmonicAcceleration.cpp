#include <iostream>

#include "Tudat/Astrodynamics/Gravitation/mutualExtendedBodySphericalHarmonicAcceleration.h"
#include "Tudat/Mathematics/BasicMathematics/basicMathematicsFunctions.h"

namespace tudat
{

namespace gravitation
{

//! Class to compute the full two-extended body gravitational attraction between two bodies
MutualExtendedBodySphericalHarmonicAcceleration::MutualExtendedBodySphericalHarmonicAcceleration(
        const boost::function< Eigen::Vector3d( ) > positionOfBody1Function,
        const boost::function< Eigen::Vector3d( ) > positionOfBody2Function,
        const boost::function< double( ) > gravitationalParameterFunction,
        const double equatorialRadiusOfBody1,
        const double equatorialRadiusOfBody2,
        const boost::function< Eigen::MatrixXd( ) > cosineHarmonicCoefficientsOfBody1Function,
        const boost::function< Eigen::MatrixXd( ) > sineHarmonicCoefficientsOfBody1Function,
        const boost::function< Eigen::MatrixXd( ) > cosineHarmonicCoefficientsOfBody2Function,
        const boost::function< Eigen::MatrixXd( ) > sineHarmonicCoefficientsOfBody2Function,
        const std::vector< boost::tuple< unsigned int, unsigned int, unsigned int, unsigned int > >& coefficientCombinationsToUse,
        const boost::function< Eigen::Quaterniond( ) > toLocalFrameOfBody1Transformation,
        const boost::function< Eigen::Quaterniond( ) > toLocalFrameOfBody2Transformation,
        const bool useCentraBodyFrame,
        const bool areCoefficientsNormalized ):
    positionOfBody1Function_( positionOfBody1Function ), positionOfBody2Function_( positionOfBody2Function ),
    gravitationalParameterFunction_( gravitationalParameterFunction ),
    equatorialRadiusOfBody1_( equatorialRadiusOfBody1 ), equatorialRadiusOfBody2_( equatorialRadiusOfBody2 ),
    coefficientCombinationsToUse_( coefficientCombinationsToUse ),
    toLocalFrameOfBody1Transformation_( toLocalFrameOfBody1Transformation ),
    toLocalFrameOfBody2Transformation_( toLocalFrameOfBody2Transformation ),
    useCentraBodyFrame_( useCentraBodyFrame ),
    areCoefficientsNormalized_( areCoefficientsNormalized )
{
    // Determine maximum degrees and orders to be evaluated.
    maximumDegree_ = 0;
    maximumOrder_ = 0;
    unsigned int degreeOfBody1, orderOfBody1, degreeOfBody2, orderOfBody2;
    for( unsigned int i = 0; i < coefficientCombinationsToUse_.size( ); i++ )
    {
        degreeOfBody1 = coefficientCombinationsToUse.at( i ).get< 0 >( );
        orderOfBody1 = coefficientCombinationsToUse.at( i ).get< 1 >( );
        degreeOfBody2 = coefficientCombinationsToUse.at( i ).get< 2 >( );
        orderOfBody2 = coefficientCombinationsToUse.at( i ).get< 3 >( );

        if( degreeOfBody1 + degreeOfBody2 > maximumDegree_ )
        {
            maximumDegree_ = degreeOfBody1 + degreeOfBody2;
        }

        if( orderOfBody1 + orderOfBody2 > maximumOrder_ )
        {
            maximumOrder_ = orderOfBody1 + orderOfBody2;
        }
    }

    // Create objects used to transform spherical harmonic coefficients, and to compute effective one-body coefficients
    sphericalHarmonicsCache_ = boost::make_shared< basic_mathematics::SphericalHarmonicsCache >(
                maximumOrder_ + 1, maximumDegree_ + 1 );
    effectiveMutualPotentialField_ =  boost::make_shared< EffectiveMutualSphericalHarmonicsField >(
                coefficientCombinationsToUse_,
                cosineHarmonicCoefficientsOfBody1Function, sineHarmonicCoefficientsOfBody1Function,
                cosineHarmonicCoefficientsOfBody2Function, sineHarmonicCoefficientsOfBody2Function,
                gravitationalParameterFunction, equatorialRadiusOfBody1_, equatorialRadiusOfBody2_, areCoefficientsNormalized );
    effectiveCosineCoefficientFunction_ = boost::bind(
                &EffectiveMutualSphericalHarmonicsField::getEffectiveCosineCoefficient,
                effectiveMutualPotentialField_, _1, _2, _3, _4 );
    effectiveSineCoefficientFunction_ = boost::bind(
                &EffectiveMutualSphericalHarmonicsField::getEffectiveSineCoefficient,
                effectiveMutualPotentialField_, _1, _2, _3, _4 );

    radius1Powers_.resize( effectiveMutualPotentialField_->getMaximumDegree1( ) + 1 );
    radius2Powers_.resize( effectiveMutualPotentialField_->getMaximumDegree2( ) + 1 );

}

//! Update member variables used by the acceleration model.
void MutualExtendedBodySphericalHarmonicAcceleration::updateMembers( const double currentTime )
{
    if( !( currentTime == currentTime_ ) )
    {
        // Compute current (relative) rotation and state
        Eigen::Quaterniond currentRotationFromInertialToBody1 = toLocalFrameOfBody1Transformation_( );
        currentRotationFromBody2ToBody1_ = currentRotationFromInertialToBody1 * toLocalFrameOfBody2Transformation_( ).inverse( );
        currentRelativePosition_ = positionOfBody1Function_( ) - positionOfBody2Function_( );
        currentBodyFixedRelativePosition_ =
                currentRotationFromInertialToBody1 * ( currentRelativePosition_ );

        std::cout<<"Current rotation to body undergoing ext.: "<<std::endl<<
                   currentRotationFromInertialToBody1.toRotationMatrix( )<<std::endl;
        std::cout<<"Current rotation to body exerting ext.: "<<std::endl<<toLocalFrameOfBody2Transformation_( ).toRotationMatrix( )<<std::endl;

        // Compute effective one-body coefficients for current state
        effectiveMutualPotentialField_->computeCurrentEffectiveCoefficients( currentRotationFromBody2ToBody1_ );

        // Compute radius powers
        double currentDistance = currentRelativePosition_.norm( );
        for( int i = 0; i <= effectiveMutualPotentialField_->getMaximumDegree1( ); i++ )
        {
            radius1Powers_[ i ] = basic_mathematics::raiseToIntegerPower( equatorialRadiusOfBody1_ / currentDistance, i );
        }

        for( int i = 0; i <= effectiveMutualPotentialField_->getMaximumDegree2( ); i++ )
        {
            radius2Powers_[ i ] = basic_mathematics::raiseToIntegerPower( equatorialRadiusOfBody2_ / currentDistance, i );
        }

        // Compute acceleration in frame fixed to body 1.
        if( areCoefficientsNormalized_ )
        {
            currentAccelerationInBodyFixedFrame_ = computeGeodesyNormalizedMutualGravitationalAccelerationSum(
                        currentBodyFixedRelativePosition_, gravitationalParameterFunction_( ),
                        equatorialRadiusOfBody1_, equatorialRadiusOfBody2_,
                        effectiveCosineCoefficientFunction_, effectiveSineCoefficientFunction_,
                        coefficientCombinationsToUse_,
                        effectiveMutualPotentialField_->getMaximumDegree1( ),
                        effectiveMutualPotentialField_->getMaximumDegree2( ),
                        maximumDegree_,
                        radius1Powers_,
                        radius2Powers_,
                        sphericalHarmonicsCache_ );
        }
        else
        {            
            currentAccelerationInBodyFixedFrame_ = computeUnnormalizedMutualGravitationalAccelerationSum(
                        currentBodyFixedRelativePosition_, gravitationalParameterFunction_( ),
                        equatorialRadiusOfBody1_, equatorialRadiusOfBody2_,
                        effectiveCosineCoefficientFunction_, effectiveSineCoefficientFunction_,
                        coefficientCombinationsToUse_, sphericalHarmonicsCache_ );
        }

        currentAcceleration_ = currentRotationFromInertialToBody1.inverse( ) * currentAccelerationInBodyFixedFrame_;
        currentTime_ = currentTime;
    }
}

}

}


