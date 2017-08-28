#ifndef ROTATIONMODELPOLYNOMIALVARIATIONS_H
#define ROTATIONMODELPOLYNOMIALVARIATIONS_H

#include "Tudat/Astrodynamics/OrbitDetermination/EstimatableParameters/estimatableParameter.h"
#include "Tudat/Astrodynamics/Ephemerides/directlyPerturbedRotationModel.h"

namespace tudat
{

namespace estimatable_parameters
{


class RotationModelPolynomialVariations: public EstimatableParameter< Eigen::VectorXd >
{

public:
    RotationModelPolynomialVariations(
            const boost::shared_ptr< ephemerides::DirectlyPerturbedRotationModel > rotationModel,
            const std::vector< int > componentPowers,
            const PerturbedRotationModelComponents componentType,
            const std::string& associatedBody ):
        EstimatableParameter< Eigen::VectorXd >( rotation_model_polynomial_compoment, associatedBody ),
        rotationModel_( rotationModel ), componentPowers_( componentPowers ), componentType_( componentType )
    { }

    ~RotationModelPolynomialVariations( ) { }

    Eigen::VectorXd getParameterValue( )
    {
        Eigen::VectorXd parameterValue = Eigen::VectorXd::Zero( getParameterSize( ) );

        for( unsigned int i = 0; i < componentPowers_.size( ); i++ )
        {
            parameterValue( i ) = rotationModel_->getPolynomialComponent( componentType_, componentPowers_[ i ] );

        }
        return parameterValue;
    }

    void setParameterValue( Eigen::VectorXd parameterValue )
    {
        for( unsigned int i = 0; i < componentPowers_.size( ); i++ )
        {
            rotationModel_->setPolynomialComponent( componentType_, componentPowers_[ i ], parameterValue( i ) );
        }
    }

    int getParameterSize( )
    {
        return componentPowers_.size( ) * 2;
    }

    PerturbedRotationModelComponents getComponentType( )
    {
        return componentType_;
    }

    std::vector< int > getComponentPowers( )
    {
        return componentPowers_;
    }

    virtual std::string getSecondaryIdentifier( )
    {
        return boost::lexical_cast< std::string >( static_cast< int >( componentType_ ) );
    }

protected:

private:
    boost::shared_ptr< ephemerides::DirectlyPerturbedRotationModel > rotationModel_;

    std::vector< int > componentPowers_;

    PerturbedRotationModelComponents componentType_;

};

}

}


#endif // ROTATIONMODELPOLYNOMIALVARIATIONS_H
