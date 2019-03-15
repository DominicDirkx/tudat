#ifndef TIDALLYLOCKEDROTATIONALEPHEMERIS_H
#define TIDALLYLOCKEDROTATIONALEPHEMERIS_H

#include <functional>

#include "Tudat/Astrodynamics/Ephemerides/rotationalEphemeris.h"
#include "Tudat/Astrodynamics/ReferenceFrames/referenceFrameTransformations.h"

namespace tudat
{

namespace ephemerides
{

Eigen::Vector6d getStateFromSelectedStateFunction(
        const double currentTime,
        const bool useFirstFunction,
        const std::function< Eigen::Vector6d( const double ) > stateFunction1,
        const std::function< Eigen::Vector6d( const double ) > stateFunction2 );

double evaluateRightAscensionDeclinationFunction(
        const double nominalValue,
        const double secularRate,
        const std::map< double, std::pair< double, double > >& periodicTerms,
        const double currentTime,
        const bool useCosines );


class TidallyLockedRotationalEphemeris: public RotationalEphemeris
{
public:
    TidallyLockedRotationalEphemeris(
            const std::function< double( const double ) > rightAscensionFunction,
            const std::function< double( const double ) > declinationFunction,
            const std::function< Eigen::Vector6d( const double, bool ) > relativeStateFunction,
            const std::string& centralBodyName,
            const std::string& baseFrameOrientation,
            const std::string& targetFrameOrientation,
            const Eigen::Quaterniond& intermediateToBaseFrameRotation = Eigen::Quaterniond( Eigen::Matrix3d::Identity( ) ) ):
        RotationalEphemeris( baseFrameOrientation, targetFrameOrientation ),
        rightAscensionFunction_( rightAscensionFunction ),
        declinationFunction_( declinationFunction ),
        relativeStateFunction_( relativeStateFunction ),
        isBodyInPropagation_( 0 ),
        centralBodyName_( centralBodyName ),
        intermediateToBaseFrameRotation_( intermediateToBaseFrameRotation )
    { }

    ~TidallyLockedRotationalEphemeris( ){ }

    Eigen::Quaterniond getRotationToBaseFrame(
            const double currentRightAscension,
            const double currentDeclination,
            const Eigen::Vector6d relativeState );

    Eigen::Quaterniond getRotationToBaseFrame( const double currentTime );

    virtual Eigen::Quaterniond getRotationToTargetFrame(
            const double currentTime )
    {
        return getRotationToBaseFrame( currentTime ).inverse( );
    }

    Eigen::Matrix3d getDerivativeOfRotationToBaseFrame( const double currentTime );

    virtual Eigen::Matrix3d getDerivativeOfRotationToTargetFrame(
            const double currentTime )
    {
        return getDerivativeOfRotationToBaseFrame( currentTime ).transpose( );
    }

    void getFullRotationalQuantitiesToTargetFrame(
            Eigen::Quaterniond& currentRotationToLocalFrame,
            Eigen::Matrix3d& currentRotationToLocalFrameDerivative,
            Eigen::Vector3d& currentAngularVelocityVectorInGlobalFrame,
            const double currentTime );

    void setIsBodyInPropagation( const bool isBodyInPropagation )
    {
        isBodyInPropagation_ = isBodyInPropagation;
    }

    Eigen::Quaterniond getIntermediateToBaseFrameRotation( )
    {
        return intermediateToBaseFrameRotation_;
    }

    std::function< double( const double ) > getRightAscensionFunction( )
    {
        return rightAscensionFunction_;
    }

    std::function< double( const double ) > getDeclinationFunction( )
    {
        return declinationFunction_;
    }

    std::string getCentralBodyName( )
    {
        return centralBodyName_;
    }

    Eigen::Vector6d getCurrentRelativeState( const double time )
    {
        return relativeStateFunction_( time, isBodyInPropagation_ );
    }

    Eigen::Quaterniond getSecondIntermediateRotationToTargetFrame(
            const double currentRightAscension, const double currentDeclination )
    {
        return reference_frames::getInertialToPlanetocentricFrameTransformationQuaternion(
                    currentDeclination, currentRightAscension, 0.0 ) * intermediateToBaseFrameRotation_.inverse( );
    }

    Eigen::Quaterniond getSecondIntermediateRotationToTargetFrame( const double currentTime )
    {
        return reference_frames::getInertialToPlanetocentricFrameTransformationQuaternion(
                    declinationFunction_( currentTime ), rightAscensionFunction_( currentTime ), 0.0 ) * intermediateToBaseFrameRotation_.inverse( );
    }

    double getFullyLockedRotationAngle( const Eigen::Vector3d& relativeCentralBodyStateInSecondIntermediateFrame )
    {
        return std::atan2( relativeCentralBodyStateInSecondIntermediateFrame( 1 ),
                           relativeCentralBodyStateInSecondIntermediateFrame( 0 ) );
    }

    double getTotalRotationAngle( const Eigen::Quaterniond& secondIntermediateRotationToTargetFrame,
                                  const Eigen::Vector6d& relativeState )
    {
        double currentRotationAngle = getFullyLockedRotationAngle( secondIntermediateRotationToTargetFrame * ( -relativeState.segment( 0, 3 ) ) );
        return currentRotationAngle;
    }


private:

    std::function< double( const double ) > rightAscensionFunction_;

    std::function< double( const double ) > declinationFunction_;

    const std::function< Eigen::Vector6d( const double, bool ) > relativeStateFunction_;

    bool isBodyInPropagation_;

    std::string centralBodyName_;

    Eigen::Quaterniond intermediateToBaseFrameRotation_;

    std::function< double( ) > gravitationalParameterFunction_;



};

}

}


#endif // TIDALLYLOCKEDROTATIONALEPHEMERIS_H
