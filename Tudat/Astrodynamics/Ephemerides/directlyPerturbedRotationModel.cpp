#include <map>
#include <vector>

#include "Tudat/Mathematics/BasicMathematics/mathematicalConstants.h"
#include "Tudat/Mathematics/BasicMathematics/basicMathematicsFunctions.h"

#include "Tudat/Basics/utilities.h"
#include "Tudat/Astrodynamics/ReferenceFrames/referenceFrameTransformations.h"
#include "Tudat/Astrodynamics/Ephemerides/rotationalEphemeris.h"
#include "Tudat/Astrodynamics/Ephemerides/directlyPerturbedRotationModel.h"

namespace tudat
{

namespace ephemerides
{

//! Function to calculate the rotation quaternion from target frame (body-fixed) to original frame.
Eigen::Quaterniond DirectlyPerturbedRotationModel::getRotationToBaseFrame( const double time )
{
    // Check if update function needs to be called.
    //if( std::fabs( time - currentTime_ ) > 1.0E-7 )
    {
        update( time );
    }

    return currentRotationFromLocalFrame_;
}

//! Function to calculate the derivative of the rotation matrix from target frame (body-fixed) to original frame.
Eigen::Matrix3d DirectlyPerturbedRotationModel::getDerivativeOfRotationToBaseFrame( const double time )
{
    //if( std::fabs( time - currentTime_ ) > 1.0E-7 )
    {
        update( time );
    }

    Eigen::AngleAxisd secondRotationAroundZaxis =
            Eigen::AngleAxisd( currentPrimeMeridian_, Eigen::Vector3d::UnitZ( ) );

    Eigen::AngleAxisd rotationAroundXaxis =
            Eigen::AngleAxisd(
                ( mathematical_constants::PI / 2.0 - currentDeclination_ ),
                Eigen::Vector3d::UnitX( ) );

    Eigen::AngleAxisd firstRotationAroundZaxis = Eigen::AngleAxisd(
                ( currentRightAscension_ + mathematical_constants::PI / 2.0 ),
                Eigen::Vector3d::UnitZ( ) );

    //NOTE: Check and fix this
    //std::cerr<<"Warning, derivative of rotation from frame not tested for DirectlyPerturbedRotationModel "<<primeMeridianPolynomialTerms_[ 1 ]<<std::endl;

    return fromIntermediateFrameToBaseFrame_ *
             primeMeridianPolynomialTerms_[ 1 ] * ( reference_frames::Z_AXIS_ROTATION_MATRIX_DERIVATIVE_PREMULTIPLIER *
            ( ( currentRotationFromLocalFrame_.inverse( ) ).toRotationMatrix( ) ) ).transpose( );
}

std::pair< double, double > DirectlyPerturbedRotationModel::getPeriodicComponentAmplitudes( const PerturbedRotationModelComponents componentType,
                                                                                            const double componentPeriod )
{
    std::pair< double, double > componentAmplitudes;

    switch( componentType )
    {
    case prime_meridian_angle:
        if( primeMeridianLibrations_.count( componentPeriod ) != 0 )
        {
            componentAmplitudes = std::make_pair( primeMeridianLibrations_.at( componentPeriod ).first,
                                                  primeMeridianLibrations_.at( componentPeriod ).second );

        }
        else
        {
            std::cerr<<"Error, could not find periodic rotation variation in prime meridian of period"<<componentPeriod<<
                       " when getting components"<<std::endl;
        }
        break;
    case declination_angle:
        if( declinationLibrations_.count( componentPeriod ) != 0 )
        {
            componentAmplitudes = std::make_pair( declinationLibrations_.at( componentPeriod ).first,
                                                  declinationLibrations_.at( componentPeriod ).second );
        }
        else
        {
            std::cerr<<"Error, could not find periodic rotation variation in declination of period"<<componentPeriod<<
                       " when getting components"<<std::endl;
        }
        break;
    case right_ascension_angle:
        if( rightAscensionLibrations_.count( componentPeriod ) != 0 )
        {
            componentAmplitudes = std::make_pair( rightAscensionLibrations_.at( componentPeriod ).first,
                                                  rightAscensionLibrations_.at( componentPeriod ).second );
        }
        else
        {
            std::cerr<<"Error, could not find periodic rotation variation in right ascension of period"<<componentPeriod<<
                       " when getting components"<<std::endl;
        }
        break;
    default:
        std::cerr<<"Error when getting rotation variation amplitudes, could not find component "<<componentType<<std::endl;
    }

    return componentAmplitudes;
}

double DirectlyPerturbedRotationModel::getPolynomialComponent( const PerturbedRotationModelComponents componentType, const int power )
{
    double component = TUDAT_NAN;

    switch( componentType )
    {
    case prime_meridian_angle:
        if( static_cast< int >( primeMeridianPolynomialTerms_.size( ) ) > power )
        {
            component = primeMeridianPolynomialTerms_.at( power );

        }
        else
        {
            std::cerr<<"Error, could not find polynomial rotation variation in prime meridian of power "<<power<<
                       " when getting components"<<std::endl;
        }
        break;
    case declination_angle:
        if( static_cast< int >( declinationPolynomialTerms_.size( ) ) > power )
        {
            component = declinationPolynomialTerms_.at( power );

        }
        else
        {
            std::cerr<<"Error, could not find polynomial rotation variation in declination of power "<<power<<
                       " when getting components"<<std::endl;
        }
        break;
    case right_ascension_angle:
        if( static_cast< int >( rightAscensionPolynomialTerms_.size( ) ) > power )
        {
            component = rightAscensionPolynomialTerms_.at( power );

        }
        else
        {
            std::cerr<<"Error, could not find polynomial rotation variation in right ascension of power "<<power<<
                       " when getting components"<<std::endl;
        }
        break;
    default:
        std::cerr<<"Error when getting rotation polynomial amplitudes, could not find component "<<componentType<<std::endl;
    }

    return component;
}


void DirectlyPerturbedRotationModel::setPeriodicComponentAmplitudes( const PerturbedRotationModelComponents componentType,
                                                                     const double componentPeriod,
                                                                     std::pair< double, double > newAmplitudes )
{
    std::pair< double, double > componentAmplitudes;

    switch( componentType )
    {
    case prime_meridian_angle:
        if( primeMeridianLibrations_.count( componentPeriod ) != 0 )
        {
            primeMeridianLibrations_.at( componentPeriod ) = newAmplitudes;
        }
        else
        {
            std::cerr<<"Error, could not find periodic rotation variation in prime meridian of period"<<componentPeriod<<
                       " when setting components"<<std::endl;
        }
        break;
    case declination_angle:
        if( declinationLibrations_.count( componentPeriod ) != 0 )
        {
            declinationLibrations_.at( componentPeriod ) = newAmplitudes;
        }
        else
        {
            std::cerr<<"Error, could not find periodic rotation variation in declination of period"<<componentPeriod<<
                       " when setting components"<<std::endl;
        }
        break;
    case right_ascension_angle:
        if( rightAscensionLibrations_.count( componentPeriod ) != 0 )
        {
            rightAscensionLibrations_.at( componentPeriod ) = newAmplitudes;
        }
        else
        {
            std::cerr<<"Error, could not find periodic rotation variation in right ascension of period"<<componentPeriod<<
                       " when setting components"<<std::endl;
        }
        break;
    default:
        std::cerr<<"Error when setting rotation variation amplitudes, could not find component "<<componentType<<std::endl;
    }
}

void DirectlyPerturbedRotationModel::setPolynomialComponent(
        const PerturbedRotationModelComponents componentType, const int componentPower, double newAmplitude )
{
    switch( componentType )
    {
    case prime_meridian_angle:
        if( static_cast< int >( primeMeridianPolynomialTerms_.size( ) ) > componentPower )
        {
            primeMeridianPolynomialTerms_[ componentPower ] = newAmplitude;
        }
        else
        {
            std::cerr<<"Error, could not find polynomial rotation variation in prime meridian of period"<<componentPower<<
                       " when setting components"<<std::endl;
        }
        break;
    case declination_angle:
        if( static_cast< int >( declinationPolynomialTerms_.size( ) ) > componentPower )
        {
            declinationPolynomialTerms_[ componentPower ] = newAmplitude;
        }
        else
        {
            std::cerr<<"Error, could not find polynomial rotation variation in declination of period"<<componentPower<<
                       " when setting components"<<std::endl;
        }
        break;
    case right_ascension_angle:
        if( static_cast< int >( rightAscensionPolynomialTerms_.size( ) ) > componentPower )
        {
            rightAscensionPolynomialTerms_[ componentPower ] = newAmplitude;
        }
        else
        {
            std::cerr<<"Error, could not find polynomial rotation variation in right ascension of period"<<componentPower<<
                       " when setting components"<<std::endl;
        }
        break;
    default:
        std::cerr<<"Error when setting polynomial rotation variation amplitudes, could not find component "<<componentType<<std::endl;
    }
}

void DirectlyPerturbedRotationModel::update( const double time )
{

    double timeToPower = 1.0;
    currentRightAscension_ = 0.0;
    for( unsigned int i = 0; i < rightAscensionPolynomialTerms_.size( ); i++ )
    {
        currentRightAscension_ += rightAscensionPolynomialTerms_[ i ] * timeToPower;
        timeToPower *= time;
    }

    timeToPower = 1.0;
    currentDeclination_ = 0.0;
    for( unsigned int i = 0; i < declinationPolynomialTerms_.size( ); i++ )
    {
        currentDeclination_ += declinationPolynomialTerms_[ i ] * timeToPower;
        timeToPower *= time;
    }

    timeToPower = 1.0;
    currentPrimeMeridian_ = 0.0;
    for( unsigned int i = 0; i < primeMeridianPolynomialTerms_.size( ); i++ )
    {
        currentPrimeMeridian_ += primeMeridianPolynomialTerms_[ i ] * timeToPower;
        timeToPower *= time;
    }

    double currentLibrationPhase = 0.0;
    for( librationIterator = rightAscensionLibrations_.begin( ); librationIterator != rightAscensionLibrations_.end( );
         librationIterator++ )
    {
        currentLibrationPhase = 2.0 * mathematical_constants::PI * time / librationIterator->first;
        currentRightAscension_ += librationIterator->second.first * std::cos( currentLibrationPhase ) +
                librationIterator->second.second * std::sin( currentLibrationPhase );
    }

    for( librationIterator = declinationLibrations_.begin( ); librationIterator != declinationLibrations_.end( );
         librationIterator++ )
    {
        currentLibrationPhase = 2.0 * mathematical_constants::PI * time / librationIterator->first;
        currentDeclination_ += librationIterator->second.first * std::cos( currentLibrationPhase ) +
                librationIterator->second.second * std::sin( currentLibrationPhase );
    }

    for( librationIterator = primeMeridianLibrations_.begin( ); librationIterator != primeMeridianLibrations_.end( );
         librationIterator++ )
    {
        currentLibrationPhase = 2.0 * mathematical_constants::PI * time / librationIterator->first;
        currentPrimeMeridian_ += librationIterator->second.first * std::cos( currentLibrationPhase ) +
                librationIterator->second.second * std::sin( currentLibrationPhase );
    }

    currentRotationFromLocalFrame_ =
            fromIntermediateFrameToBaseFrame_ * (
                reference_frames::getInertialToPlanetocentricFrameTransformationQuaternion(
                    currentDeclination_, currentRightAscension_, currentPrimeMeridian_ ) ).inverse( );

    currentTime_ = time;
}

}

}


