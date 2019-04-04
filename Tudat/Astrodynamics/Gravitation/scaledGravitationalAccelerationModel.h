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

#ifndef TUDAT_SCALED_GRAVITATIONAL_ACCELERATION_MODEL_H
#define TUDAT_SCALED_GRAVITATIONAL_ACCELERATION_MODEL_H

#include <iostream>

#include <boost/function.hpp>
#include <boost/lambda/lambda.hpp>

#include <Eigen/Core>

#include "Tudat/Astrodynamics/Gravitation/sphericalHarmonicsGravityModel.h"
#include "Tudat/Astrodynamics/Gravitation/centralGravityModel.h"

#include "Tudat/Astrodynamics/Gravitation/mutualSphericalHarmonicGravityModel.h"


namespace tudat
{
namespace gravitation
{

template< typename BaseAccelerationModel >
class ScaledGravitationalAccelerationModel: public BaseAccelerationModel
{
};

template< >
class ScaledGravitationalAccelerationModel< CentralGravitationalAccelerationModel3d >: public CentralGravitationalAccelerationModel3d
{
public:
    ScaledGravitationalAccelerationModel(
            const std::shared_ptr< CentralGravitationalAccelerationModel3d > originalAccelerationModel,
            const std::function< double( ) > newGravitationalParameterFunction,
            const bool newIsMutualAttractionUsed,
            const bool invertPositionVectors ):
        CentralGravitationalAccelerationModel3d( originalAccelerationModel ),
        originalAccelerationModel_( originalAccelerationModel ),
        invertPositionVectors_( invertPositionVectors )
    {
        this->gravitationalParameterFunction = newGravitationalParameterFunction;
        this->isMutualAttractionUsed_ = newIsMutualAttractionUsed;

        if( invertPositionVectors )
        {
            std::function< Eigen::Vector3d( ) > tempFunction = this->subjectPositionFunction;
            this->subjectPositionFunction = this->sourcePositionFunction;
            this->sourcePositionFunction = tempFunction;
        }
    }

    ~ScaledGravitationalAccelerationModel( ){ }

    //! Update members.
    /*!
     * Updates class members relevant for computing the central gravitational acceleration. In this
     * case the function simply updates the members in the base class.
     * \sa SphericalHarmonicsGravitationalAccelerationModelBase.
     */
    void updateMembers( const double currentTime = TUDAT_NAN )
    {
        originalAccelerationModel_->updateMembers( currentTime );

        if( !( this->currentTime_ == currentTime ) )
        {

            this->currentAcceleration_ = originalAccelerationModel_->getAcceleration( ) *
                    ( this->getGravitationalParameterFunction( )( ) /
                      originalAccelerationModel_->getGravitationalParameterFunction( )( ) );

            this->updateBaseMembers( );

            if( invertPositionVectors_ )
            {
                this->currentAcceleration_ *= -1.0;
            }
            this->currentTime_ = currentTime;
        }
    }

    bool getInvertPositionVectors( )
    {
        return invertPositionVectors_;
    }

    std::shared_ptr< CentralGravitationalAccelerationModel3d > getOriginalAccelerationModel( )
    {
        return originalAccelerationModel_;
    }

    void resetTime( const double currentTime = TUDAT_NAN )
    {
        if( !( currentTime_ == currentTime  ) )
        {
            originalAccelerationModel_->resetTime( currentTime );

        }
        currentTime_ = currentTime;
    }


private:
    const std::shared_ptr< CentralGravitationalAccelerationModel3d > originalAccelerationModel_;

    bool invertPositionVectors_;

    Eigen::Vector3d currentAcceleration_;
};

template< >
class ScaledGravitationalAccelerationModel< SphericalHarmonicsGravitationalAccelerationModel >: public SphericalHarmonicsGravitationalAccelerationModel
{
public:
    ScaledGravitationalAccelerationModel(
            const std::shared_ptr< SphericalHarmonicsGravitationalAccelerationModel > originalAccelerationModel,
            const std::function< double( ) > newGravitationalParameterFunction,
            const bool newIsMutualAttractionUsed,
            const bool invertPositionVectors ):
        SphericalHarmonicsGravitationalAccelerationModel( originalAccelerationModel ),
        originalAccelerationModel_( originalAccelerationModel ),
        invertPositionVectors_( invertPositionVectors )
    {
        this->gravitationalParameterFunction = newGravitationalParameterFunction;
        this->isMutualAttractionUsed_ = newIsMutualAttractionUsed;

        if( invertPositionVectors )
        {
            std::cerr<<"Error when making scaled spherical harmonic gravity acceleration, cannot invert position vectors"<<std::endl;
        }
    }

    ~ScaledGravitationalAccelerationModel( ){ }

    //! Update members.
    /*!
     * Updates class members relevant for computing the central gravitational acceleration. In this
     * case the function simply updates the members in the base class.
     * \sa SphericalHarmonicsGravitationalAccelerationModelBase.
     */
    void updateMembers( const double currentTime = TUDAT_NAN )
    {
        originalAccelerationModel_->updateMembers( currentTime );

        if( !( this->currentTime_ == currentTime ) )
        {
            this->currentAcceleration_ = originalAccelerationModel_->getAcceleration( ) *
                    ( this->getGravitationalParameterFunction( )( ) /
                      originalAccelerationModel_->getGravitationalParameterFunction( )( ) );

            this->updateBaseMembers( );

            cosineHarmonicCoefficients = originalAccelerationModel_->getCurrentCosineHarmonicCoefficients( );
            sineHarmonicCoefficients = originalAccelerationModel_->getCurrentSineHarmonicCoefficients( );

            if( invertPositionVectors_ )
            {
                this->currentAcceleration_ *= -1.0;
            }
            this->currentTime_ = currentTime;
        }
    }

    void resetTime( const double currentTime = TUDAT_NAN )
    {
        if( !( currentTime_ == currentTime  ) )
        {
            originalAccelerationModel_->resetTime( currentTime );

        }
        currentTime_ = currentTime;
    }

    bool getInvertPositionVectors( )
    {
        return invertPositionVectors_;
    }

    std::shared_ptr< SphericalHarmonicsGravitationalAccelerationModel > getOriginalAccelerationModel( )
    {
        return originalAccelerationModel_;
    }

private:
    const std::shared_ptr< SphericalHarmonicsGravitationalAccelerationModel > originalAccelerationModel_;

    bool invertPositionVectors_;
};


template< >
class ScaledGravitationalAccelerationModel< MutualSphericalHarmonicsGravitationalAccelerationModel >: public MutualSphericalHarmonicsGravitationalAccelerationModel
{
public:
    ScaledGravitationalAccelerationModel(
            const std::shared_ptr< MutualSphericalHarmonicsGravitationalAccelerationModel > originalAccelerationModel,
            const std::function< double( ) > newGravitationalParameterFunction,
            const bool newIsMutualAttractionUsed,
            const bool invertPositionVectors ):
        MutualSphericalHarmonicsGravitationalAccelerationModel( originalAccelerationModel ),
        invertPositionVectors_( invertPositionVectors ),
        isMutualAttractionUsed_( newIsMutualAttractionUsed )
    {
        if( invertPositionVectors )
        {
            accelerationModelFromShExpansionOfBodyExertingAcceleration_ =
                    std::make_shared< ScaledGravitationalAccelerationModel< SphericalHarmonicsGravitationalAccelerationModel > >(
                        originalAccelerationModel->getAccelerationModelFromShExpansionOfBodyUndergoingAcceleration( ),
                        newGravitationalParameterFunction, newIsMutualAttractionUsed, 0 );
            accelerationModelFromShExpansionOfBodyUndergoingAcceleration_ =
                    std::make_shared< ScaledGravitationalAccelerationModel< SphericalHarmonicsGravitationalAccelerationModel > >(
                        originalAccelerationModel->getAccelerationModelFromShExpansionOfBodyExertingAcceleration( ),
                        newGravitationalParameterFunction, newIsMutualAttractionUsed, 0 );
        }
        else
        {
            accelerationModelFromShExpansionOfBodyExertingAcceleration_ =
                    std::make_shared< ScaledGravitationalAccelerationModel< SphericalHarmonicsGravitationalAccelerationModel > >(
                        originalAccelerationModel->getAccelerationModelFromShExpansionOfBodyExertingAcceleration( ),
                        newGravitationalParameterFunction, newIsMutualAttractionUsed, 0 );
            accelerationModelFromShExpansionOfBodyUndergoingAcceleration_ =
                    std::make_shared< ScaledGravitationalAccelerationModel< SphericalHarmonicsGravitationalAccelerationModel > >(
                        originalAccelerationModel->getAccelerationModelFromShExpansionOfBodyUndergoingAcceleration( ),
                        newGravitationalParameterFunction, newIsMutualAttractionUsed, 0 );
        }



        this->useCentralBodyFixedFrame_ = newIsMutualAttractionUsed;
        this->gravitationalParameterFunction_ = newGravitationalParameterFunction;

    }

    ~ScaledGravitationalAccelerationModel( ){ }

    //! Update members.
    /*!
     * Updates class members relevant for computing the central gravitational acceleration. In this
     * case the function simply updates the members in the base class.
     * \sa SphericalHarmonicsGravitationalAccelerationModelBase.
     */
    void updateMembers( const double currentTime = TUDAT_NAN )
    {
        accelerationModelFromShExpansionOfBodyExertingAcceleration_->updateMembers( currentTime );
        accelerationModelFromShExpansionOfBodyUndergoingAcceleration_->updateMembers( currentTime );
    }

    bool getInvertPositionVectors( )
    {
        return invertPositionVectors_;
    }

    bool getIsMutualAttractionUsed( )
    {
        return isMutualAttractionUsed_;
    }

    void resetTime( const double currentTime = TUDAT_NAN )
    {
        accelerationModelFromShExpansionOfBodyExertingAcceleration_->resetTime( currentTime );
        accelerationModelFromShExpansionOfBodyUndergoingAcceleration_->resetTime( currentTime );

        currentTime_ = currentTime;
    }


private:

    bool invertPositionVectors_;

    bool isMutualAttractionUsed_;

    Eigen::Vector3d currentAcceleration_;

};

} // namespace gravitation
} // namespace tudat

#endif // TUDAT_SCALED_GRAVITATIONAL_ACCELERATION_MODEL_H
