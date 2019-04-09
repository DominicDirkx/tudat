/*    Copyright (c) 2010-2014, Delft University of Technology
 *    All rights reserved.
 *
 *    Redistribution and use in source and binary forms, with or without modification, are
 *    permitted provided that the following conditions are met:
 *      - Redistributions of source code must retain the above copyright notice, this list of
 *        conditions and the following disclaimer.
 *      - Redistributions in binary form must reproduce the above copyright notice, this list of
 *        conditions and the following disclaimer in the documentation and/or other materials
 *        provided with the distribution.
 *      - Neither the name of the Delft University of Technology nor the names of its contributors
 *        may be used to endorse or promote products derived from this software without specific
 *        prior written permission.
 *
 *    THIS SOFTWARE IS PROVIDED BY THE COPYRIGHT HOLDERS AND CONTRIBUTORS "AS IS" AND ANY EXPRESS
 *    OR IMPLIED WARRANTIES, INCLUDING, BUT NOT LIMITED TO, THE IMPLIED WARRANTIES OF
 *    MERCHANTABILITY AND FITNESS FOR A PARTICULAR PURPOSE ARE DISCLAIMED. IN NO EVENT SHALL THE
 *    COPYRIGHT HOLDER OR CONTRIBUTORS BE LIABLE FOR ANY DIRECT, INDIRECT, INCIDENTAL, SPECIAL,
 *    EXEMPLARY, OR CONSEQUENTIAL DAMAGES (INCLUDING, BUT NOT LIMITED TO, PROCUREMENT OF SUBSTITUTE
 *    GOODS OR SERVICES; LOSS OF USE, DATA, OR PROFITS; OR BUSINESS INTERRUPTION) HOWEVER CAUSED
 *    AND ON ANY THEORY OF LIABILITY, WHETHER IN CONTRACT, STRICT LIABILITY, OR TORT (INCLUDING
 *    NEGLIGENCE OR OTHERWISE) ARISING IN ANY WAY OUT OF THE USE OF THIS SOFTWARE, EVEN IF ADVISED
 *    OF THE POSSIBILITY OF SUCH DAMAGE.
 *
 *    Changelog
 *      YYMMDD    Author            Comment
 *      120827    K. Kumar          File created.
 *      121105    K. Kumar          Simplified base class definition.
 *      121210    D. Dirkx          Simplified class by removing template parameters.
 *
 *    References
 *
 *    Notes
 *
 */

#ifndef TUDAT_SCALED_GRAVITATIONAL_ACCELERATION_PARTIAL_H
#define TUDAT_SCALED_GRAVITATIONAL_ACCELERATION_PARTIAL_H

#include <iostream>

#include <boost/function.hpp>
#include <boost/lambda/lambda.hpp>

#include <Eigen/Core>

#include "Tudat/Astrodynamics/OrbitDetermination/AccelerationPartials/centralGravityAccelerationPartial.h"
#include "Tudat/Astrodynamics/OrbitDetermination/AccelerationPartials/sphericalHarmonicAccelerationPartial.h"
#include "Tudat/Astrodynamics/OrbitDetermination/AccelerationPartials/mutualSphericalHarmonicGravityPartial.h"


namespace tudat
{

namespace acceleration_partials
{

template< typename BaseAccelerationPartial >
class ScaledGravitationalAccelerationPartial: public BaseAccelerationPartial
{
};

template< >
class ScaledGravitationalAccelerationPartial< CentralGravitationPartial >: public CentralGravitationPartial
{
public:
    ScaledGravitationalAccelerationPartial(
            const std::shared_ptr< CentralGravitationPartial > originalAccelerationPartial,
            const std::function< double( ) > newGravitationalParameterFunction,
            const bool newIsMutualAttractionUsed,
            const bool invertPositionVectors ):
        CentralGravitationPartial( originalAccelerationPartial ),
        originalAccelerationPartial_( originalAccelerationPartial ),
        invertPositionVectors_( invertPositionVectors )
    {
        this->gravitationalParameterFunction_ = newGravitationalParameterFunction;
        this->accelerationUsesMutualAttraction_ = newIsMutualAttractionUsed;

        currentGravitationalParameterRatio_ =
                gravitationalParameterFunction_( ) /  originalAccelerationPartial_->getGravitationalParameterFunction( )( );

        if( invertPositionVectors )
        {
            std::function< Eigen::Vector3d( ) > tempFunction = this->centralBodyState_;
            this->centralBodyState_ = this->acceleratedBodyState_;
            this->acceleratedBodyState_ = tempFunction;

            std::string tempBodyId = acceleratedBody_;
            acceleratedBody_ = acceleratingBody_;
            acceleratingBody_ = tempBodyId;
        }
    }

    ~ScaledGravitationalAccelerationPartial( ){ }

    void update( const double currentTime = TUDAT_NAN )
    {
        accelerationUpdateFunction_( currentTime );

        if( !( currentTime_ == currentTime ) )
        {
            currentAcceleratedBodyState_ = acceleratedBodyState_( );
            currentCentralBodyState_ = centralBodyState_( );
            currentGravitationalParameterRatio_ = gravitationalParameterFunction_( ) /  originalAccelerationPartial_->getGravitationalParameterFunction( )( );

            originalAccelerationPartial_->update( currentTime );
            currentPartialWrtPosition_ = originalAccelerationPartial_->getCurrentPartialWrtPosition( );
            scalePartial( currentPartialWrtPosition_ );

            currentTime_ = currentTime;
        }
    }

    std::pair< std::function< void( Eigen::MatrixXd& ) >, int > getParameterPartialFunction(
            std::shared_ptr< estimatable_parameters::EstimatableParameter< double > > parameter )

    {
        // Initialize function (automatically empty function) and number of columns to zero (no dependency)
        std::function< void( Eigen::MatrixXd& ) > partialFunction;
        int numberOfColumns = 0;

        std::pair< std::function< void( Eigen::MatrixXd& ) >, int > originalPartialFunction =
                originalAccelerationPartial_->getParameterPartialFunction( parameter );

        if( parameter->getParameterName( ).first != estimatable_parameters::gravitational_parameter )
        {
            partialFunction = std::bind( &ScaledGravitationalAccelerationPartial< CentralGravitationPartial >::scalePartialFunction,
                                         this, originalPartialFunction.first, std::placeholders::_1 );
            numberOfColumns =  originalPartialFunction.second;
        }
        else
        {
            std::pair< std::function< void( Eigen::MatrixXd& ) >, int > partialFunctionPair =
                    this->getGravitationalParameterPartialFunction( parameter->getParameterName( ) );
            partialFunction = partialFunctionPair.first;
            numberOfColumns = partialFunctionPair.second;
        }

        return std::make_pair( partialFunction, numberOfColumns );
    }

    std::pair< std::function< void( Eigen::MatrixXd& ) >, int > getParameterPartialFunction(
            std::shared_ptr< estimatable_parameters::EstimatableParameter< Eigen::VectorXd > > parameter )
    {
        std::function< void( Eigen::MatrixXd& ) > partialFunction;
        int numberOfColumns = 0;

        std::pair< std::function< void( Eigen::MatrixXd& ) >, int > originalPartialFunction =
                originalAccelerationPartial_->getParameterPartialFunction( parameter );

        partialFunction = std::bind( &ScaledGravitationalAccelerationPartial< CentralGravitationPartial >::scalePartialFunction,
                                     this, originalPartialFunction.first, std::placeholders::_1 );
        numberOfColumns =  originalPartialFunction.second;


        return std::make_pair( partialFunction, numberOfColumns );
    }

    int setParameterPartialUpdateFunction(
            std::shared_ptr< estimatable_parameters::EstimatableParameter< double > > parameter )
    {
        std::pair< std::function< void( Eigen::MatrixXd& ) >, int > parameterPartialFunction = getParameterPartialFunction( parameter );

        if( parameterPartialFunction.second > 0 && parameterDoublePartialFunctions_.count( parameter ) == 0 )
        {
            if( parameter->getParameterName( ).first == estimatable_parameters::gravitational_parameter )
            {
                if( parameterDoublePartialFunctions_.count( parameter ) == 0 )
                {
                    parameterDoublePartialFunctions_[ parameter ] = parameterPartialFunction.first;
                    isCurrentDoubleParameterPartialSet_[ parameter ] = 0;
                    currentDoubleParameterPartials_[ parameter ] = Eigen::MatrixXd( 3, 1 );
                }
            }
            else
            {
                parameterDoublePartialFunctions_[ parameter ] =
                        std::bind( static_cast< void( ScaledGravitationalAccelerationPartial< CentralGravitationPartial >::*)
                                   ( std::shared_ptr< estimatable_parameters::EstimatableParameter< double > >,
                                     Eigen::MatrixXd& ) >
                                   ( &ScaledGravitationalAccelerationPartial< CentralGravitationPartial >::getScaledParameterPartial) ,
                                   this, parameter,  std::placeholders::_1 );
                isCurrentDoubleParameterPartialSet_[ parameter ] = 0;
                currentDoubleParameterPartials_[ parameter ] = Eigen::MatrixXd( 3, 1 );

                originalAccelerationPartial_->setParameterPartialUpdateFunction( parameter );
            }
        }

        return parameterPartialFunction.second;
    }

    int setParameterPartialUpdateFunction(
            std::shared_ptr< estimatable_parameters::EstimatableParameter< Eigen::VectorXd > > parameter )
    {
        std::pair< std::function< void( Eigen::MatrixXd& ) >, int > parameterPartialFunction = this->getParameterPartialFunction( parameter );

        if( parameterPartialFunction.second > 0 && parameterVectorPartialFunctions_.count( parameter ) == 0 )
        {
            parameterVectorPartialFunctions_[ parameter ] =
                    std::bind( static_cast< void( ScaledGravitationalAccelerationPartial< CentralGravitationPartial >::*)
                               ( std::shared_ptr< estimatable_parameters::EstimatableParameter< Eigen::VectorXd > >,
                                 Eigen::MatrixXd&  ) >
                               ( &ScaledGravitationalAccelerationPartial< CentralGravitationPartial >::getScaledParameterPartial ) ,
                               this, parameter, std::placeholders::_1 );
            isCurrentVectorParameterPartialSet_[ parameter ] = 0;
            currentVectorParameterPartials_[ parameter ] = Eigen::MatrixXd( 3, parameter->getParameterSize( ) );

            originalAccelerationPartial_->setParameterPartialUpdateFunction( parameter );
        }

        return parameterPartialFunction.second;
    }

    bool getInvertPositionVectors( )
    {
        return invertPositionVectors_;
    }

    std::shared_ptr< CentralGravitationPartial > getOriginalAccelerationPartial( )
    {
        return originalAccelerationPartial_;
    }

protected:


    void resetTimeOfMemberObjects( )
    {
        originalAccelerationPartial_->resetTime( currentTime_ );
    }

    void updateParameterPartialsOfMemberObjects( )
    {
        originalAccelerationPartial_->updateParameterPartials( );
    }

private:

    double getPartialScalingValue( )
    {
        return currentGravitationalParameterRatio_;
    }

    void scalePartial( Eigen::MatrixXd& partialToScale )
    {
        partialToScale *= getPartialScalingValue( );
    }

    void getScaledParameterPartial(
            std::shared_ptr< estimatable_parameters::EstimatableParameter< double > > parameter,
            Eigen::MatrixXd& partialsMatrix )
    {
        partialsMatrix.setZero( 3, 1 );
        originalAccelerationPartial_->getCurrentParameterPartial(
                    parameter, partialsMatrix.block( 0, 0, 3, 1 ) );
        scalePartial( partialsMatrix );
    }

    void getScaledParameterPartial(
            std::shared_ptr< estimatable_parameters::EstimatableParameter< Eigen::VectorXd > > parameter,
            Eigen::MatrixXd& partialsMatrix )
    {
        partialsMatrix.setZero( 3, parameter->getParameterSize( ) );
        originalAccelerationPartial_->getCurrentParameterPartial(
                    parameter, partialsMatrix.block( 0, 0, 3, parameter->getParameterSize( ) ) );
        scalePartial( partialsMatrix );
    }

    void scalePartialFunction( const std::function< void( Eigen::MatrixXd& ) > originalPartialFunction,
                               Eigen::MatrixXd& partialsMatrix )
    {
        originalPartialFunction( partialsMatrix );
        scalePartial( partialsMatrix );
    }

    const std::shared_ptr< CentralGravitationPartial > originalAccelerationPartial_;

    bool invertPositionVectors_;

    double currentGravitationalParameterRatio_;

};

template< >
class ScaledGravitationalAccelerationPartial< SphericalHarmonicsGravityPartial >: public SphericalHarmonicsGravityPartial
{
public:
    ScaledGravitationalAccelerationPartial(
            const std::shared_ptr< SphericalHarmonicsGravityPartial > originalAccelerationPartial,
            const std::function< double( ) > newGravitationalParameterFunction,
            const bool newIsMutualAttractionUsed,
            const bool invertPositionVectors ):
        SphericalHarmonicsGravityPartial( originalAccelerationPartial ),
        originalAccelerationPartial_( originalAccelerationPartial ),
        invertPositionVectors_( invertPositionVectors )
    {
        this->gravitationalParameterFunction_ =
                newGravitationalParameterFunction;
        this->accelerationUsesMutualAttraction_ =
                newIsMutualAttractionUsed;
        this->accelerationFunction_ = std::bind(
                    &ScaledGravitationalAccelerationPartial< SphericalHarmonicsGravityPartial >::getScaledAcceleration,
                    this );

        currentGravitationalParameterRatio_ = gravitationalParameterFunction_( ) /
                originalAccelerationPartial_->getGravitationalParameterFunction( )( );

        if( invertPositionVectors )
        {
            throw std::runtime_error(
                        "Error when making scaled spherical harmonic gravity partial, cannot invert position vectors" );
        }
    }

    ~ScaledGravitationalAccelerationPartial( ){ }

    void update( const double currentTime = TUDAT_NAN )
    {
        originalAccelerationPartial_->update( currentTime );

        if( this->currentTime_ != currentTime )
        {
            currentGravitationalParameter_ = gravitationalParameterFunction_( );

            currentGravitationalParameterRatio_ =
                    currentGravitationalParameter_ /  originalAccelerationPartial_->getCurrentGravitationalParameter( );

            bodyFixedSphericalPosition_ = originalAccelerationPartial_->getBodyFixedSphericalPosition( );
            bodyFixedPosition_ = originalAccelerationPartial_->getBodyFixedPosition( );
            this->currentAcceleration_ = accelerationFunction_( );

            this->currentPartialWrtPosition_ = originalAccelerationPartial_->getCurrentPartialWrtPosition( );
            scalePartial( this->currentPartialWrtPosition_ );
            this->currentPartialWrtVelocity_ = originalAccelerationPartial_->getCurrentPartialWrtVelocity( ) ;
            scalePartial( this->currentPartialWrtVelocity_ );

            this->currentTime_ = currentTime;


        }
    }

    std::pair< std::function< void( Eigen::MatrixXd& ) >, int >
    getParameterPartialFunction( std::shared_ptr< estimatable_parameters::EstimatableParameter< double > > parameter )

    {
        // Initialize function (automatically empty function) and number of columns to zero (no dependency)
        std::function< void( Eigen::MatrixXd& ) > partialFunction;
        int numberOfColumns = 0;

        std::pair< std::function< void( Eigen::MatrixXd& ) >, int > originalPartialFunction =
                originalAccelerationPartial_->getParameterPartialFunction( parameter );

        if( parameter->getParameterName( ).first != estimatable_parameters::gravitational_parameter )
        {
            partialFunction = std::bind( &ScaledGravitationalAccelerationPartial< SphericalHarmonicsGravityPartial >::scalePartialFunction,
                                         this, originalPartialFunction.first, std::placeholders::_1 );
            numberOfColumns =  originalPartialFunction.second;
        }
        else
        {
            std::pair< std::function< void( Eigen::MatrixXd& ) >, int > partialFunctionPair =
                    this->getGravitationalParameterPartialFunction( parameter->getParameterName( ) );
            partialFunction = partialFunctionPair.first;
            numberOfColumns = partialFunctionPair.second;
        }

        return std::make_pair( partialFunction, numberOfColumns );
    }

    std::pair< std::function< void( Eigen::MatrixXd& ) >, int >
    getParameterPartialFunction( std::shared_ptr< estimatable_parameters::EstimatableParameter< Eigen::VectorXd > > parameter )

    {
        // Initialize function (automatically empty function) and number of columns to zero (no dependency)
        std::function< void( Eigen::MatrixXd& ) > partialFunction;
        int numberOfColumns = 0;

        std::pair< std::function< void( Eigen::MatrixXd& ) >, int > originalPartialFunction =
                originalAccelerationPartial_->getParameterPartialFunction( parameter );

        partialFunction = std::bind( &ScaledGravitationalAccelerationPartial< SphericalHarmonicsGravityPartial >::scalePartialFunction,
                                     this, originalPartialFunction.first, std::placeholders::_1 );
        numberOfColumns =  originalPartialFunction.second;

        return std::make_pair( partialFunction, numberOfColumns );
    }

    int setParameterPartialUpdateFunction(
            std::shared_ptr< estimatable_parameters::EstimatableParameter< double > > parameter )
    {
        std::pair< std::function< void( Eigen::MatrixXd& ) >, int > parameterPartialFunction =
                getParameterPartialFunction( parameter );

        if( parameterPartialFunction.second > 0 && parameterDoublePartialFunctions_.count( parameter ) == 0 )
        {
            if( parameter->getParameterName( ).first == estimatable_parameters::gravitational_parameter )
            {
                if( parameterDoublePartialFunctions_.count( parameter ) == 0 )
                {
                    parameterDoublePartialFunctions_[ parameter ] = parameterPartialFunction.first;
                    isCurrentDoubleParameterPartialSet_[ parameter ] = 0;
                    currentDoubleParameterPartials_[ parameter ] = Eigen::MatrixXd( 3, 1 );
                }
            }
            else
            {
                parameterDoublePartialFunctions_[ parameter ] =
                        std::bind( static_cast< void( ScaledGravitationalAccelerationPartial< SphericalHarmonicsGravityPartial >::*)
                                   ( std::shared_ptr< estimatable_parameters::EstimatableParameter< double > >,
                                     Eigen::MatrixXd& ) >
                                   ( &ScaledGravitationalAccelerationPartial< SphericalHarmonicsGravityPartial >::getScaledParameterPartial) ,
                                   this, parameter,  std::placeholders::_1 );
                isCurrentDoubleParameterPartialSet_[ parameter ] = 0;
                currentDoubleParameterPartials_[ parameter ] = Eigen::MatrixXd( 3, 1 );
                originalAccelerationPartial_->setParameterPartialUpdateFunction( parameter );
            }
        }

        return parameterPartialFunction.second;
    }

    int setParameterPartialUpdateFunction(
            std::shared_ptr< estimatable_parameters::EstimatableParameter< Eigen::VectorXd > > parameter )
    {
        std::pair< std::function< void( Eigen::MatrixXd& ) >, int > parameterPartialFunction = getParameterPartialFunction( parameter );

        if( parameterPartialFunction.second > 0 && parameterVectorPartialFunctions_.count( parameter ) == 0 )
        {
            parameterVectorPartialFunctions_[ parameter ] =
                    std::bind( static_cast< void( ScaledGravitationalAccelerationPartial< SphericalHarmonicsGravityPartial >::*)
                               ( std::shared_ptr< estimatable_parameters::EstimatableParameter< Eigen::VectorXd > >,
                                 Eigen::MatrixXd&  ) >
                               ( &ScaledGravitationalAccelerationPartial< SphericalHarmonicsGravityPartial >::getScaledParameterPartial ) ,
                               this, parameter, std::placeholders::_1 );

            originalAccelerationPartial_->setParameterPartialUpdateFunction( parameter );
            isCurrentVectorParameterPartialSet_[ parameter ] = 0;
            currentVectorParameterPartials_[ parameter ] = Eigen::MatrixXd( 3, parameter->getParameterSize( ) );

        }

        return parameterPartialFunction.second;
    }

    bool getInvertPositionVectors( )
    {
        return invertPositionVectors_;
    }

    std::shared_ptr< SphericalHarmonicsGravityPartial > getoriginalAccelerationPartial( )
    {
        return originalAccelerationPartial_;
    }

    virtual Eigen::Vector3d wrtGravitationalParameterOfCentralBody( )
    {
        return currentAcceleration_ / originalAccelerationPartial_->getCurrentGravitationalParameter( );
    }



protected:

    Eigen::Vector3d getScaledAcceleration( )
    {
        return originalAccelerationPartial_->getAccelerationFunction( )( ) *
                ( this->getCurrentGravitationalParameter( ) /
                  originalAccelerationPartial_->getCurrentGravitationalParameter( ) ) *
                ( invertPositionVectors_ ? -1.0 : 1.0 );
    }

    void updateParameterPartialsOfMemberObjects( )
    {
        //std::cout<<"Updating parameter partials (scaled)"<<std::endl;
        originalAccelerationPartial_->updateParameterPartials( );
    }


    void resetTimeOfMemberObjects( )
    {
        originalAccelerationPartial_->resetTime( currentTime_ );
    }

private:

    double getPartialScalingValue( )
    {
        return currentGravitationalParameterRatio_;
    }

    void scalePartial( Eigen::Matrix3d& originalPartial )
    {
        originalPartial *= getPartialScalingValue( );
    }

    void scalePartial( Eigen::MatrixXd& originalPartial )
    {
        originalPartial *= getPartialScalingValue( );
    }

    void getScaledParameterPartial(
            std::shared_ptr< estimatable_parameters::EstimatableParameter< double > > parameter,
            Eigen::MatrixXd& partialsMatrix )
    {
        partialsMatrix.setZero( 3, 1 );
        originalAccelerationPartial_->getCurrentParameterPartial(
                    parameter, partialsMatrix.block( 0, 0, 3, 1 ) );
//        std::cout<<"Undergoing: "<<this->acceleratedBody_<<std::endl;
//        std::cout<<"Exerting: "<<this->acceleratingBody_<<std::endl;
//        std::cout<<"Unscaled partial "<<std::endl<<partialsMatrix<<std::endl;
        scalePartial( partialsMatrix );
//        std::cout<<"Scaled partial "<<std::endl<<partialsMatrix<<std::endl<<std::endl;
    }

    void getScaledParameterPartial(
            std::shared_ptr< estimatable_parameters::EstimatableParameter< Eigen::VectorXd > > parameter,
            Eigen::MatrixXd& partialsMatrix )
    {
        partialsMatrix.setZero( 3, parameter->getParameterSize( ) );
        originalAccelerationPartial_->getCurrentParameterPartial(
                    parameter, partialsMatrix.block( 0, 0, 3, parameter->getParameterSize( ) ) );

//        std::cout<<"Undergoing: "<<this->acceleratedBody_<<std::endl;
//        std::cout<<"Exerting: "<<this->acceleratingBody_<<std::endl;
//        std::cout<<"Unscaled partial "<<std::endl<<partialsMatrix<<std::endl;
        scalePartial( partialsMatrix );
//        std::cout<<"Scaled partial "<<std::endl<<partialsMatrix<<std::endl<<std::endl;

    }


    void scalePartialFunction( const std::function< void( Eigen::MatrixXd& ) > originalPartialFunction,
                               Eigen::MatrixXd& partialsMatrix )
    {
        originalPartialFunction( partialsMatrix );
        scalePartial( partialsMatrix );
    }

    const std::shared_ptr< SphericalHarmonicsGravityPartial > originalAccelerationPartial_;

    bool invertPositionVectors_;

    double currentGravitationalParameterRatio_;

};


template< >
class ScaledGravitationalAccelerationPartial< MutualSphericalHarmonicsGravityPartial >: public MutualSphericalHarmonicsGravityPartial
{
public:
    ScaledGravitationalAccelerationPartial(
            const std::shared_ptr< MutualSphericalHarmonicsGravityPartial > originalAccelerationPartial,
            const std::function< double( ) > newGravitationalParameterFunction,
            const bool newIsMutualAttractionUsed,
            const bool invertPositionVectors ):
        MutualSphericalHarmonicsGravityPartial( originalAccelerationPartial ),
        invertPositionVectors_( invertPositionVectors )
    {
        if( invertPositionVectors )
        {
            this->accelerationPartialOfShExpansionOfBodyExertingAcceleration_ =
                    std::make_shared< ScaledGravitationalAccelerationPartial< SphericalHarmonicsGravityPartial > >(
                        originalAccelerationPartial->getAccelerationPartialOfShExpansionOfBodyUndergoingAcceleration( ),
                        newGravitationalParameterFunction, newIsMutualAttractionUsed, 0 );
            this->accelerationPartialOfShExpansionOfBodyUndergoingAcceleration_ =
                    std::make_shared< ScaledGravitationalAccelerationPartial< SphericalHarmonicsGravityPartial > >(
                        originalAccelerationPartial->getAccelerationPartialOfShExpansionOfBodyExertingAcceleration( ),
                        newGravitationalParameterFunction, newIsMutualAttractionUsed, 0 );

            std::string tempBodyId = this->acceleratedBody_;
            this->acceleratedBody_ = this->acceleratingBody_;
            this->acceleratingBody_ = tempBodyId;
        }
        else
        {
            this->accelerationPartialOfShExpansionOfBodyExertingAcceleration_ =
                    std::make_shared< ScaledGravitationalAccelerationPartial< SphericalHarmonicsGravityPartial > >(
                        originalAccelerationPartial->getAccelerationPartialOfShExpansionOfBodyExertingAcceleration( ),
                        newGravitationalParameterFunction, newIsMutualAttractionUsed, 0 );
            this->accelerationPartialOfShExpansionOfBodyUndergoingAcceleration_ =
                    std::make_shared< ScaledGravitationalAccelerationPartial< SphericalHarmonicsGravityPartial > >(
                        originalAccelerationPartial->getAccelerationPartialOfShExpansionOfBodyUndergoingAcceleration( ),
                        newGravitationalParameterFunction, newIsMutualAttractionUsed, 0 );
        }


        this->accelerationUsesMutualAttraction_ = newIsMutualAttractionUsed;

    }

    ~ScaledGravitationalAccelerationPartial( ){ }

    bool getInvertPositionVectors( )
    {
        return invertPositionVectors_;
    }

private:

    bool invertPositionVectors_;

};

}

} // namespace tudat

#endif // TUDAT_SCALED_GRAVITATIONAL_ACCELERATION_PARTIAL_H
