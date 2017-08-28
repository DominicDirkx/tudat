#ifndef DIRECTPERTURBEDROTATIONMODELPARTIALS_H
#define DIRECTPERTURBEDROTATIONMODELPARTIALS_H

#include <boost/shared_ptr.hpp>

#include "Tudat/Astrodynamics/Ephemerides/directlyPerturbedRotationModel.h"
#include "Tudat/Astrodynamics/OrbitDetermination/ObservationPartials/rotationMatrixPartial.h"

namespace tudat
{

namespace observation_partials
{


class DirectPerturbedRotationModelPartialManager
{
public:
    DirectPerturbedRotationModelPartialManager( const boost::shared_ptr< ephemerides::DirectlyPerturbedRotationModel > rotationModel ):
        rotationModel_( rotationModel )
    {
        subRotationsFromLocalFrame_.resize( 3 );
        subRotationDerivativesFromLocalFrame_.resize( 3 );

        zDerivativeAuxiliaryMatrix << 0.0, 1.0, 0.0, -1.0, 0.0, 0.0, 0.0, 0.0, 0.0;
        xDerivativeAuxiliaryMatrix << 0.0, 0.0, 0.0, 0.0, 0.0, 1.0, 0.0, -1.0, 0.0;

    }

    Eigen::Matrix3d getDerivativeOfRotationMatrixWrtPeriodicComponent(
            const double time, const double componentPeriod, const PerturbedRotationModelComponents componentType, const bool isComponentCosine )
    {
        if( std::fabs( time - currentTime_ ) > 1.0E-7 )
        {
            update( time );
            currentTime_ = time;
        }

        Eigen::Matrix3d rotationMatrixPartial = Eigen::Matrix3d::Identity( );

        for( unsigned int i = 0; i < 3; i++ )
        {
            if( i == componentType )
            {
                rotationMatrixPartial = subRotationDerivativesFromLocalFrame_[ i ] * rotationMatrixPartial;
            }
            else
            {
                rotationMatrixPartial = subRotationsFromLocalFrame_[ i ] * rotationMatrixPartial;
            }
        }

        rotationMatrixPartial = getPeriodicComponent( time, componentPeriod, isComponentCosine ) * rotationMatrixPartial;

        return rotationModel_->getFromIntermediateFrameToBaseFrame( ) * rotationMatrixPartial;
    }

    Eigen::Matrix3d getDerivativeOfRotationMatrixWrtPolynomialComponent(
            const double time, const double componentPower, const PerturbedRotationModelComponents componentType )
    {
        if( std::fabs( time - currentTime_ ) > 1.0E-7 )
        {
            update( time );
            currentTime_ = time;
        }

        Eigen::Matrix3d rotationMatrixPartial = Eigen::Matrix3d::Identity( );

        for( unsigned int i = 0; i < 3; i++ )
        {
            if( i == componentType )
            {
                rotationMatrixPartial = subRotationDerivativesFromLocalFrame_[ i ] * rotationMatrixPartial;
            }
            else
            {
                rotationMatrixPartial = subRotationsFromLocalFrame_[ i ] * rotationMatrixPartial;
            }
        }

        rotationMatrixPartial = std::pow( time, componentPower ) * rotationMatrixPartial;

        return rotationModel_->getFromIntermediateFrameToBaseFrame( ) * rotationMatrixPartial;
    }

    double getPeriodicComponent( const double time, const double componentPeriod, const bool isComponentCosine )
    {
        double component = 0.0;
        if( isComponentCosine )
        {
            component = std::cos( 2.0 * mathematical_constants::PI / componentPeriod * time );
        }
        else
        {
            component = std::sin( 2.0 * mathematical_constants::PI / componentPeriod * time );
        }
        return component;
    }

    boost::shared_ptr< ephemerides::DirectlyPerturbedRotationModel > getRotationModel( )
    {
        return rotationModel_;
    }


private:
    void update( const double time )
    {
        Eigen::Matrix3d primeMeridianRotation =
                Eigen::AngleAxisd( rotationModel_->getCurrentPrimeMeridian( time ), Eigen::Vector3d::UnitZ( ) ).toRotationMatrix( );

        Eigen::Matrix3d declinationRotation =
                Eigen::AngleAxisd( ( mathematical_constants::PI / 2.0 -
                                     rotationModel_->getCurrentDeclination( time ) ), Eigen::Vector3d::UnitX( ) ).toRotationMatrix( );

        Eigen::Matrix3d rightAscensionRotation =
                    Eigen::AngleAxisd( rotationModel_->getCurrentRightAscension( time ) + mathematical_constants::PI / 2.0,
                                       Eigen::Vector3d::UnitZ( ) ).toRotationMatrix( );

        subRotationsFromLocalFrame_[ 2 ] = rightAscensionRotation;
        subRotationsFromLocalFrame_[ 1 ] = declinationRotation;
        subRotationsFromLocalFrame_[ 0 ] = primeMeridianRotation;

        subRotationDerivativesFromLocalFrame_[ 2 ] = - zDerivativeAuxiliaryMatrix * rightAscensionRotation;
        subRotationDerivativesFromLocalFrame_[ 1 ] = xDerivativeAuxiliaryMatrix * declinationRotation;
        subRotationDerivativesFromLocalFrame_[ 0 ] = -zDerivativeAuxiliaryMatrix * primeMeridianRotation;

    }

    boost::shared_ptr< ephemerides::DirectlyPerturbedRotationModel > rotationModel_;

    std::vector< Eigen::Matrix3d > subRotationsFromLocalFrame_;

    std::vector< Eigen::Matrix3d > subRotationDerivativesFromLocalFrame_;

    double currentTime_;

    Eigen::Matrix3d zDerivativeAuxiliaryMatrix;

    Eigen::Matrix3d xDerivativeAuxiliaryMatrix;
};

class RotationMatrixPartialWrtPeriodicRotationVariation: public RotationMatrixPartial
{
public:
    RotationMatrixPartialWrtPeriodicRotationVariation(
            const boost::shared_ptr< DirectPerturbedRotationModelPartialManager > rotationVariationPartialManager,
            const std::vector< double > rotationPeriods,
            const PerturbedRotationModelComponents componentType ):
        RotationMatrixPartial( rotationVariationPartialManager->getRotationModel( ) ),
        rotationVariationPartialManager_( rotationVariationPartialManager ), rotationPeriods_( rotationPeriods ),
        componentType_( componentType ){ }

    ~RotationMatrixPartialWrtPeriodicRotationVariation( ){ }


    std::vector< Eigen::Matrix3d > calculatePartialOfRotationMatrixToBaseFrameWrParameter(
            const double time )
    {
        std::vector< Eigen::Matrix3d > rotationMatrixPartials;

        for( unsigned int i = 0; i < rotationPeriods_.size( ); i++ )
        {
            rotationMatrixPartials.push_back(
                        rotationVariationPartialManager_->getDerivativeOfRotationMatrixWrtPeriodicComponent(
                            time, rotationPeriods_[ i ], componentType_, 1 ) );
            rotationMatrixPartials.push_back(
                        rotationVariationPartialManager_->getDerivativeOfRotationMatrixWrtPeriodicComponent(
                            time, rotationPeriods_[ i ], componentType_, 0 ) );
        }

        return rotationMatrixPartials;
    }

    std::vector< Eigen::Matrix3d > calculatePartialOfRotationMatrixDerivativeToBaseFrameWrParameter(
                const double time )
    {
        throw std::runtime_error( "Partial A not yet implemented" );
    }

    std::string getSecondaryIdentifier( )
    {
        return boost::lexical_cast< std::string >( static_cast< int >( componentType_ ) );
    }
private:

    boost::shared_ptr< DirectPerturbedRotationModelPartialManager > rotationVariationPartialManager_;
    std::vector< double > rotationPeriods_;
    PerturbedRotationModelComponents componentType_;


};

class RotationMatrixPartialWrtPolynimialRotationVariation: public RotationMatrixPartial
{
public:
    RotationMatrixPartialWrtPolynimialRotationVariation(
            const boost::shared_ptr< DirectPerturbedRotationModelPartialManager > rotationVariationPartialManager,
            const std::vector< int > componentPowers,
            const PerturbedRotationModelComponents componentType ):
        RotationMatrixPartial( rotationVariationPartialManager->getRotationModel( ) ),
        rotationVariationPartialManager_( rotationVariationPartialManager ), componentPowers_( componentPowers ),
        componentType_( componentType ){ }

    ~RotationMatrixPartialWrtPolynimialRotationVariation( ){ }

    std::vector< Eigen::Matrix3d > calculatePartialOfRotationMatrixToBaseFrameWrParameter(
            const double time )
    {
        std::vector< Eigen::Matrix3d > rotationMatrixPartials;

        for( unsigned int i = 0; i < componentPowers_.size( ); i++ )
        {
            rotationMatrixPartials.push_back(
                        rotationVariationPartialManager_->getDerivativeOfRotationMatrixWrtPolynomialComponent(
                            time, componentPowers_[ i ], componentType_ ) );
        }

        return rotationMatrixPartials;
    }

    std::vector< Eigen::Matrix3d > calculatePartialOfRotationMatrixDerivativeToBaseFrameWrParameter(
                const double time )
    {
        throw std::runtime_error( "Partial B not yet implemented" );
    }


    std::string getSecondaryIdentifier( )
    {
        return boost::lexical_cast< std::string >( static_cast< int >( componentType_ ) );
    }
private:

    boost::shared_ptr< DirectPerturbedRotationModelPartialManager > rotationVariationPartialManager_;
    std::vector< int > componentPowers_;
    PerturbedRotationModelComponents componentType_;


};

}

}

#endif // DIRECTPERTURBEDROTATIONMODELPARTIALS_H
