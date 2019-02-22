
#include <vector>
#include <map>

#include "Tudat/Astrodynamics/Ephemerides/fullPlanetaryRotationModel.h"

namespace tudat
{

namespace ephemerides
{

void PlanetaryOrientationAngleCalculator::updateCorrections( const double ephemerisTime )
{
    currentEphemerisTime_ = ephemerisTime;
    double currentMeanMotion = bodyMeanMotion_;
    double currentMeanAnomaly = bodyMeanAnomalyAtEpoch_ + bodyMeanMotion_ * currentEphemerisTime_;

    currentAngleICorrection_ = 0.0;
    currentAnglePsiCorrection_ = 0.0;

    for( std::map< double, std::pair< double, double > >::iterator correctionIterator = meanMotionDirectNutationCorrections_.begin( );
         correctionIterator != meanMotionDirectNutationCorrections_.end( ); correctionIterator++ )
    {
        currentAngleICorrection_ += correctionIterator->second.first * std::sin(
                    correctionIterator->first * ( currentMeanMotion * currentEphemerisTime_ + bodyMeanAnomalyAtEpoch_ )  );
        currentAnglePsiCorrection_ += correctionIterator->second.second * std::cos(
                    correctionIterator->first * ( currentMeanMotion * currentEphemerisTime_ + bodyMeanAnomalyAtEpoch_ )  );
    }

    double currentPhaseAngleCorrection;

    for( unsigned int i = 0; i < meanMotionTimeDependentPhaseNutationCorrections_.size( ); i++ )
    {
        currentPhaseAngleCorrection = phaseAngleCorrectionFunctions_[ i ]( currentEphemerisTime_ );
        for( std::map< double, std::pair< double, double > >::iterator correctionIterator =
             meanMotionTimeDependentPhaseNutationCorrections_[ i ].begin( );
             correctionIterator != meanMotionTimeDependentPhaseNutationCorrections_[ i ].end( ); correctionIterator++ )
        {
            currentAngleICorrection_ += correctionIterator->second.first * std::sin(
                        correctionIterator->first * ( currentMeanMotion * currentEphemerisTime_ + bodyMeanAnomalyAtEpoch_ ) + currentPhaseAngleCorrection );
            currentAnglePsiCorrection_ += correctionIterator->second.second * std::cos(
                        correctionIterator->first * ( currentMeanMotion * currentEphemerisTime_ + bodyMeanAnomalyAtEpoch_ ) + currentPhaseAngleCorrection );
        }
    }


    currentAnglePsi_ = anglePsiAtEpoch_ + anglePsiRateAtEpoch_ * currentEphemerisTime_ + currentAnglePsiCorrection_;
    currentAngleI_ = angleIAtEpoch_ + angleIRateAtEpoch_ * currentEphemerisTime_ + currentAngleICorrection_;

    currentAnglePhiCorrection_ = -currentAnglePsiCorrection_ * std::cos( currentAngleI_ );

    for( std::map< double, std::pair< double, double > >::iterator correctionIterator = rotationRateCorrections_.begin( );
         correctionIterator != rotationRateCorrections_.end( ); correctionIterator++ )
    {
        currentAnglePhiCorrection_ += correctionIterator->second.first * std::cos(
                    correctionIterator->first * currentMeanAnomaly );
        currentAnglePhiCorrection_ += correctionIterator->second.second * std::sin(
                    correctionIterator->first * currentMeanAnomaly );
    }

    currentAnglePhi_ = anglePhiAtEpoch_ + anglePhiRateAtEpoch_ * currentEphemerisTime_ + currentAnglePhiCorrection_;

    //std::cout<<currentAngleICorrection_<<" "<<currentAnglePsiCorrection_<<" "<<currentAnglePhiCorrection_<<std::endl;

}


boost::shared_ptr< interpolators::CubicSplineInterpolator< double, Eigen::Vector3d > >
createInterpolatorForPlanetaryRotationAngles( double intervalStart,
                                               double intervalEnd,
                                               double timeStep,
                                               boost::shared_ptr< PlanetaryOrientationAngleCalculator > planetaryOrientationCalculator )
{
    using namespace interpolators;

    std::map< double, Eigen::Vector3d > orientationMap;

    double currentTime = intervalStart;
    while( currentTime < intervalEnd )
    {
        orientationMap[ currentTime ] = planetaryOrientationCalculator->updateAndGetRotationAngles(
                    currentTime  );
        currentTime += timeStep;
    }

    boost::shared_ptr< CubicSplineInterpolator< double, Eigen::Vector3d > > interpolator =
            boost::make_shared< CubicSplineInterpolator< double, Eigen::Vector3d > >( orientationMap );
    return interpolator;
}

Eigen::Quaterniond PlanetaryRotationModel::getRotationFromBodyFixedToIntermediateInertialFrame( const double ephemerisTime )
{
    Eigen::Vector3d currentAngleCorrections = planetaryOrientationAnglesCalculator_->updateAndGetRotationAngles(
                ephemerisTime );
    return Eigen::Quaterniond(
                Eigen::AngleAxisd( currentAngleCorrections.x( ), Eigen::Vector3d::UnitZ( ) ) *
                Eigen::AngleAxisd( currentAngleCorrections.y( ), Eigen::Vector3d::UnitX( ) ) *
                Eigen::AngleAxisd( currentAngleCorrections.z( ), Eigen::Vector3d::UnitZ( ) ) );
}

Eigen::Quaterniond PlanetaryRotationModel::getRotationToBaseFrame( const double ephemerisTime )
{
    return rotationFromMeanOrbitToIcrf_ * getRotationFromBodyFixedToIntermediateInertialFrame( ephemerisTime );
}


Eigen::Matrix3d PlanetaryRotationModel::getDerivativeOfRotationToBaseFrame( const double ephemerisTime )
{
    Eigen::Vector3d currentAngleCorrections = planetaryOrientationAnglesCalculator_->updateAndGetRotationAngles(
                ephemerisTime );
    double currentPhiAngle = currentAngleCorrections.z( );

    //NOTE: Check and fix this
    std::cerr<<"Warning, derivative of rotation from frame not tested for PlanetaryRotationModel"<<std::endl;


    return ( Eigen::Matrix3d( ) << -std::sin( currentPhiAngle ), -std::cos( currentPhiAngle ), 0.0,
             std::cos( currentPhiAngle ), -std::sin( currentPhiAngle ), 0.0, 0.0, 0.0, 0.0 ).finished( );

//    return meanPhiAngleDerivativeFunction_( ) * ( rotationFromMeanOrbitToIcrf_ ).toRotationMatrix( ) *
//            ( Eigen::AngleAxisd( currentAngleCorrections.x( ), Eigen::Vector3d::UnitZ( ) ) *
//                             Eigen::AngleAxisd( currentAngleCorrections.y( ), Eigen::Vector3d::UnitX( ) ) ).toRotationMatrix( );
}

}

}

