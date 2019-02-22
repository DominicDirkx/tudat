#ifndef FULLPLANETARYROTATIONMODEL_H
#define FULLPLANETARYROTATIONMODEL_H

#include <vector>
#include <map>

#include <boost/function.hpp>

#include "Tudat/External/SpiceInterface/spiceInterface.h"

#include "Tudat/Astrodynamics/Ephemerides/rotationalEphemeris.h"
#include "Tudat/Mathematics/Interpolators/lagrangeInterpolator.h"

namespace tudat
{

namespace ephemerides
{

class PlanetaryOrientationAngleCalculator
{
public:
    PlanetaryOrientationAngleCalculator( )
    { }

    PlanetaryOrientationAngleCalculator(
            const double anglePsiAtEpoch, const double anglePsiRateAtEpoch, const double angleIAtEpoch, const double angleIRateAtEpoch,
            const double anglePhiAtEpoch, const double anglePhiRateAtEpoch, const double bodyMeanMotion, const double bodyMeanAnomalyAtEpoch,
            const std::string baseFrame,
            const std::map< double, std::pair< double, double > > meanMotionDirectNutationCorrections =
            ( std::map< double, std::pair< double, double > >( ) ), // Konopliv table 5, 1st four terms
            const std::map< double, std::pair< double, double > > rotationRateCorrections =
            ( std::map< double, std::pair< double, double > >( ) ),
            const std::vector< std::map< double, std::pair< double, double > > > meanMotionTimeDependentPhaseNutationCorrections =
            ( std::vector< std::map< double, std::pair< double, double > > > ( ) ),
            const std::vector< std::function< double( const double ) > > phaseAngleCorrectionFunctions =
            ( std::vector< std::function< double( const double ) > > ( ) ) ): // Konopliv table 7, 1st four terms
        anglePsiAtEpoch_( anglePsiAtEpoch ), anglePsiRateAtEpoch_( anglePsiRateAtEpoch ), angleIAtEpoch_( angleIAtEpoch ),
        angleIRateAtEpoch_( angleIRateAtEpoch ), anglePhiAtEpoch_( anglePhiAtEpoch ), anglePhiRateAtEpoch_( anglePhiRateAtEpoch ),
        meanMotionDirectNutationCorrections_( meanMotionDirectNutationCorrections ), baseFrame_( baseFrame ),
        rotationRateCorrections_( rotationRateCorrections ),
        meanMotionTimeDependentPhaseNutationCorrections_( meanMotionTimeDependentPhaseNutationCorrections ),
        phaseAngleCorrectionFunctions_( phaseAngleCorrectionFunctions ), bodyMeanMotion_( bodyMeanMotion ),
        bodyMeanAnomalyAtEpoch_( bodyMeanAnomalyAtEpoch ){ }

    Eigen::Vector3d updateAndGetRotationAngles( const double ephemerisTime )
    {
        if( std::fabs( currentEphemerisTime_ - ephemerisTime ) > 1.0E-8 )
        {
            updateCorrections( ephemerisTime );
        }
        return ( Eigen::Vector3d( ) << currentAnglePsi_, currentAngleI_, currentAnglePhi_ ).finished( );
    }

    Eigen::Quaterniond getRotationMatrixFromMeanOfDateEquatorToInertialPlanetCenteredAtEpoch( const double ephemerisTime )
    {
        return Eigen::Quaterniond( Eigen::AngleAxisd( anglePsiAtEpoch_ + ephemerisTime * anglePsiRateAtEpoch_, Eigen::Vector3d::UnitZ( ) )  *
                                   Eigen::AngleAxisd( angleIAtEpoch_ + ephemerisTime * angleIRateAtEpoch_, Eigen::Vector3d::UnitX( ) ) );
    }

    double getCurrentAnglePsi( )
    {
        return currentAnglePsi_;
    }

    double getCurrentAngleI( )
    {
        return currentAngleI_;
    }

    double getCurrentAnglePhi( )
    {
        return currentAnglePhi_;
    }

    double getMeanPhiAngleDerivative( )
    {
        return anglePhiRateAtEpoch_;
    }

    std::string getBaseFrame( )
    {
        return baseFrame_;
    }

    double getAnglePsiRateAtEpoch( )
    {
        return anglePsiRateAtEpoch_;
    }

    void setAnglePsiRateAtEpoch( const double anglePsiRateAtEpoch )
    {
        anglePsiRateAtEpoch_ = anglePsiRateAtEpoch;
        currentEphemerisTime_ = -1.0E100;
    }


private:

    void updateCorrections( const double ephemerisTime );

    double currentEphemerisTime_;

    double bodyMeanAnomalyAtEpoch_;
    double bodyMeanMotion_;

    double currentAnglePsiCorrection_;
    double currentAngleICorrection_;
    double currentAnglePhiCorrection_;

    double currentAnglePsi_;
    double currentAngleI_;
    double currentAnglePhi_;

    double anglePsiAtEpoch_;
    double anglePsiRateAtEpoch_;
    double angleIAtEpoch_;
    double angleIRateAtEpoch_;
    double anglePhiAtEpoch_;
    double anglePhiRateAtEpoch_;

    std::string baseFrame_;

    std::map< double, std::pair< double, double > > meanMotionDirectNutationCorrections_;

    std::vector< std::map< double, std::pair< double, double > > > meanMotionTimeDependentPhaseNutationCorrections_;

    std::vector< std::function< double( const double ) > > phaseAngleCorrectionFunctions_;

    std::map< double, std::pair< double, double > > rotationRateCorrections_;
};

std::shared_ptr< interpolators::CubicSplineInterpolator< double, Eigen::Vector3d > >
createInterpolatorForPlanetaryRotationAngles( double intervalStart,
                                              double intervalEnd,
                                              double timeStep,
                                              std::shared_ptr< PlanetaryOrientationAngleCalculator > planetaryOrientationCalculator );

class PlanetaryRotationModel: public RotationalEphemeris
{
public:
    PlanetaryRotationModel( const double angleN,
                            const double angleJ,
                            const std::shared_ptr< PlanetaryOrientationAngleCalculator > planetaryOrientationAnglesCalculator,
                            const std::string& baseFrameOrientation = "",
                            const std::string& targetFrameOrientation = "" ):
        RotationalEphemeris( baseFrameOrientation, targetFrameOrientation ),
        planetaryOrientationAnglesCalculator_( planetaryOrientationAnglesCalculator )
//      orientationAnglesFunction_( orientationAnglesFunction ),
//        meanPhiAngleDerivativeFunction_( meanPhiAngleDerivativeFunction )
    {
        rotationFromMeanOrbitToIcrf_ = spice_interface::computeRotationQuaternionBetweenFrames( "J2000", "ECLIPJ2000", 0.0 ) *
                Eigen::AngleAxisd( angleN, Eigen::Vector3d::UnitZ( ) ) *
                Eigen::AngleAxisd( angleJ, Eigen::Vector3d::UnitX( ) );
    }

    Eigen::Quaterniond getRotationFromBodyFixedToIntermediateInertialFrame( const double ephemerisTime );

    Eigen::Quaterniond getRotationToBaseFrame( const double ephemerisTime );

    Eigen::Quaterniond getRotationToTargetFrame( const double ephemerisTime )
    {
        return getRotationToBaseFrame( ephemerisTime ).inverse( );
    }

    Eigen::Matrix3d getDerivativeOfRotationToBaseFrame( const double ephemerisTime );

    Eigen::Matrix3d getDerivativeOfRotationToTargetFrame(
            const double secondsSinceEpoch )
    {
        return getDerivativeOfRotationToBaseFrame( secondsSinceEpoch ).
                transpose( );
    }


    Eigen::Quaterniond getRotationFromMeanOrbitToIcrf( )
    {
        return rotationFromMeanOrbitToIcrf_;
    }

private:

    std::shared_ptr< PlanetaryOrientationAngleCalculator > planetaryOrientationAnglesCalculator_;

    Eigen::Quaterniond rotationFromMeanOrbitToIcrf_;

//    std::function< Eigen::Vector3d( const double ) > orientationAnglesFunction_;

//    std::function< double( ) > meanPhiAngleDerivativeFunction_;

};

}

}

#endif // FULLPLANETARYROTATIONMODEL_H
