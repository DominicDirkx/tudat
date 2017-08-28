#ifndef ROTATIONMODELPERIODICVARIATIONAMPLITUDES_H
#define ROTATIONMODELPERIODICVARIATIONAMPLITUDES_H

#include "Tudat/Astrodynamics/OrbitDetermination/EstimatableParameters/estimatableParameter.h"
#include "Tudat/Astrodynamics/Ephemerides/directlyPerturbedRotationModel.h"

namespace tudat
{

namespace estimatable_parameters
{


class RotationModelPeriodicVariationAmplitudes: public EstimatableParameter< Eigen::VectorXd >
{

public:
    RotationModelPeriodicVariationAmplitudes(
            const boost::shared_ptr< ephemerides::DirectlyPerturbedRotationModel > rotationModel,
            const std::vector< double > componentPeriods,
            const PerturbedRotationModelComponents componentType,
            const std::string& associatedBody ):
        EstimatableParameter< Eigen::VectorXd >( rotation_model_component_perturbation_amplitude, associatedBody ),
        rotationModel_( rotationModel ), componentPeriods_( componentPeriods ), componentType_( componentType )
    { }

    ~RotationModelPeriodicVariationAmplitudes( ) { }

    Eigen::VectorXd getParameterValue( )
    {
        Eigen::VectorXd parameterValue = Eigen::VectorXd::Zero( getParameterSize( ) );
        std::pair< double, double > currentComponents;

        for( unsigned int i = 0; i < componentPeriods_.size( ); i++ )
        {
            currentComponents = rotationModel_->getPeriodicComponentAmplitudes(
                        componentType_, componentPeriods_[ i ] );
            parameterValue( 2 * i ) = currentComponents.first;
            parameterValue( 2 * i + 1 ) = currentComponents.second;

        }
        return parameterValue;
    }

    void setParameterValue( Eigen::VectorXd parameterValue )
    {
        std::pair< double, double > currentComponents;

        for( unsigned int i = 0; i < componentPeriods_.size( ); i++ )
        {
            currentComponents = std::make_pair( parameterValue( 2 * i ), parameterValue( 2 * i + 1 ) );
            rotationModel_->setPeriodicComponentAmplitudes( componentType_, componentPeriods_[ i ], currentComponents );
        }
    }

    int getParameterSize( ){ return componentPeriods_.size( ) * 2; }

    PerturbedRotationModelComponents getComponentType( )
    {
        return componentType_;
    }

    std::vector< double > getComponentPeriods( )
    {
        return componentPeriods_;
    }

    virtual std::string getSecondaryIdentifier( )
    {
        return boost::lexical_cast< std::string >( static_cast< int >( componentType_ ) );
    }

protected:

private:
    boost::shared_ptr< ephemerides::DirectlyPerturbedRotationModel > rotationModel_;

    std::vector< double > componentPeriods_;

    PerturbedRotationModelComponents componentType_;

};

}

}

#endif // ROTATIONMODELPERIODICVARIATIONAMPLITUDES_H
