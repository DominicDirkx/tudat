/*    Copyright (c) 2010-2018, Delft University of Technology
 *    All rigths reserved
 *
 *    This file is part of the Tudat. Redistribution and use in source and
 *    binary forms, with or without modification, are permitted exclusively
 *    under the terms of the Modified BSD license. You should have received
 *    a copy of the license with this file. If not, please or visit:
 *    http://tudat.tudelft.nl/LICENSE.
 */

#include "Tudat/SimulationSetup/EstimationSetup/createAccelerationPartials.h"
#include "Tudat/Astrodynamics/Gravitation/basicSolidBodyTideGravityFieldVariations.h"
#include "Tudat/Astrodynamics/Gravitation/gravityFieldVariations.h"

namespace tudat
{

namespace simulation_setup
{

bool isAccelerationPartialScaledOriginal(
        const std::shared_ptr< acceleration_partials::AccelerationPartial > accelerationPartial )
{
    bool isAccelerationPartialScaled = 0;
    switch( accelerationPartial->getAccelerationType( ) )
    {
    case basic_astrodynamics::central_gravity:
        if( std::dynamic_pointer_cast<  acceleration_partials::ScaledGravitationalAccelerationPartial< acceleration_partials::CentralGravitationPartial > >( accelerationPartial ) != NULL )
        {
            isAccelerationPartialScaled = 1;
        }
        break;
    case basic_astrodynamics::spherical_harmonic_gravity:
        if( std::dynamic_pointer_cast<  acceleration_partials::ScaledGravitationalAccelerationPartial< acceleration_partials::SphericalHarmonicsGravityPartial > >( accelerationPartial ) != NULL )
        {
            isAccelerationPartialScaled = 1;
        }
        break;
    case basic_astrodynamics::mutual_spherical_harmonic_gravity:
        if( std::dynamic_pointer_cast<  acceleration_partials::ScaledGravitationalAccelerationPartial< acceleration_partials::MutualSphericalHarmonicsGravityPartial > >( accelerationPartial ) != NULL )
        {
            isAccelerationPartialScaled = 1;
        }
        break;
    default:
        break;
    }
    return isAccelerationPartialScaled;
}

std::shared_ptr< acceleration_partials::AccelerationPartial > findExistingAccelerationPartial(
        const std::shared_ptr< acceleration_partials::AccelerationPartial >& accelerationPartial,
        const basic_astrodynamics::AvailableAcceleration accelerationType,
        const bool checkThirdBodyCentralAcceleration,
        const std::string& nameOfBodyUndergoingAcceleration )
{
    std::shared_ptr< acceleration_partials::AccelerationPartial > existingAccelerationPartial;

    if( accelerationPartial->getAccelerationType( ) == accelerationType )
    {
        if( !checkThirdBodyCentralAcceleration )
        {
            if( !isAccelerationPartialScaledOriginal( accelerationPartial ) )
            {
                existingAccelerationPartial = accelerationPartial;
            }
        }
    }
    else if( accelerationPartial->getAccelerationType( ) == getThirdBodyEquivalentAccelerationModelType( accelerationType ) )
    {
        std::shared_ptr< acceleration_partials::AccelerationPartial > potentialExistingAccelerationPartial;

        switch( accelerationType )
        {
        case basic_astrodynamics::central_gravity:
            if( !checkThirdBodyCentralAcceleration )
            {
                potentialExistingAccelerationPartial = std::dynamic_pointer_cast< acceleration_partials::ThirdBodyGravityPartial< acceleration_partials::CentralGravitationPartial > >(
                            accelerationPartial )->getPartialOfDirectGravityOnBodyUndergoingAcceleration( );
            }
            else
            {
                std::shared_ptr< acceleration_partials::ThirdBodyGravityPartial< acceleration_partials::CentralGravitationPartial > > thirdBodyPartial =
                        std::dynamic_pointer_cast< acceleration_partials::ThirdBodyGravityPartial< acceleration_partials::CentralGravitationPartial > >(
                            accelerationPartial );
                if( thirdBodyPartial->getCentralBodyName( ) == nameOfBodyUndergoingAcceleration )
                {
                    potentialExistingAccelerationPartial = thirdBodyPartial->getPartialOfDirectGravityOnCentralBody( );
                }
            }
            break;
        case basic_astrodynamics::spherical_harmonic_gravity:
            if( !checkThirdBodyCentralAcceleration )
            {
                potentialExistingAccelerationPartial = std::dynamic_pointer_cast< acceleration_partials::ThirdBodyGravityPartial< acceleration_partials::SphericalHarmonicsGravityPartial > >(
                            accelerationPartial )->getPartialOfDirectGravityOnBodyUndergoingAcceleration( );
            }
            else
            {
                std::shared_ptr< acceleration_partials::ThirdBodyGravityPartial< acceleration_partials::SphericalHarmonicsGravityPartial > > thirdBodyPartial =
                        std::dynamic_pointer_cast< acceleration_partials::ThirdBodyGravityPartial< acceleration_partials::SphericalHarmonicsGravityPartial > >(
                            accelerationPartial );
                if( thirdBodyPartial->getCentralBodyName( ) == nameOfBodyUndergoingAcceleration )
                {
                    potentialExistingAccelerationPartial = thirdBodyPartial->getPartialOfDirectGravityOnCentralBody( );
                }
            }
            break;
        case basic_astrodynamics::mutual_spherical_harmonic_gravity:
            if( !checkThirdBodyCentralAcceleration )
            {
                potentialExistingAccelerationPartial = std::dynamic_pointer_cast< acceleration_partials::ThirdBodyGravityPartial< acceleration_partials::MutualSphericalHarmonicsGravityPartial > >(
                            accelerationPartial )->getPartialOfDirectGravityOnBodyUndergoingAcceleration( );
            }
            else
            {
                std::shared_ptr< acceleration_partials::ThirdBodyGravityPartial< acceleration_partials::MutualSphericalHarmonicsGravityPartial > > thirdBodyPartial =
                        std::dynamic_pointer_cast< acceleration_partials::ThirdBodyGravityPartial< acceleration_partials::MutualSphericalHarmonicsGravityPartial > >(
                            accelerationPartial );
                if( thirdBodyPartial->getCentralBodyName( ) == nameOfBodyUndergoingAcceleration )
                {
                    potentialExistingAccelerationPartial = thirdBodyPartial->getPartialOfDirectGravityOnCentralBody( );
                }
            }
            break;
        default:
            std::cerr<<"Error when finding scaled acceleration partial, did not recognize third body partial."<<std::endl;
        }

        if( potentialExistingAccelerationPartial != NULL )
        {
            if( !isAccelerationPartialScaledOriginal( potentialExistingAccelerationPartial ) )
            {
                existingAccelerationPartial = potentialExistingAccelerationPartial;
            }
        }

    }
    return existingAccelerationPartial;
}

std::shared_ptr< acceleration_partials::AccelerationPartial > findExistingAccelerationPartial(
        const std::map< std::string, std::map< std::string, std::vector< std::shared_ptr<
        acceleration_partials::AccelerationPartial > > > >& accelerationPartialsMap,
        const std::string& originalAcceleratedBodyName,
        const std::string& originalAcceleratingBodyName,
        const basic_astrodynamics::AvailableAcceleration accelerationType )
{
    std::shared_ptr< acceleration_partials::AccelerationPartial > existingAccelerationPartial;

    int numberOfFeasibleAccelerations = 0;

    if( accelerationPartialsMap.count( originalAcceleratedBodyName ) > 0 )
    {
        if( accelerationPartialsMap.at( originalAcceleratedBodyName ).count( originalAcceleratingBodyName ) > 0 )
        {
            for( unsigned int i = 0; i < accelerationPartialsMap.at( originalAcceleratedBodyName ).at( originalAcceleratingBodyName ).size( ); i++ )
            {
                std::shared_ptr< acceleration_partials::AccelerationPartial > currentTestedPartial = findExistingAccelerationPartial(
                            accelerationPartialsMap.at( originalAcceleratedBodyName ).at( originalAcceleratingBodyName ).at( i ), accelerationType,
                            0, "" );
                if( currentTestedPartial != NULL )
                {
                    existingAccelerationPartial = currentTestedPartial;
                    numberOfFeasibleAccelerations++;
                }
            }
        }
    }

    for( std::map< std::string, std::map< std::string, std::vector< std::shared_ptr< acceleration_partials::AccelerationPartial > > > >::const_iterator partialIterator =
         accelerationPartialsMap.begin( ); partialIterator != accelerationPartialsMap.end( ); partialIterator++ )
    {
        if( partialIterator->second.count( originalAcceleratingBodyName ) > 0 )
        {
            for( unsigned int i = 0; i < partialIterator->second.at( originalAcceleratingBodyName ).size( ); i++ )
            {
                std::shared_ptr< acceleration_partials::AccelerationPartial > currentTestedPartial = findExistingAccelerationPartial(
                            partialIterator->second.at( originalAcceleratingBodyName ).at( i ), accelerationType,
                            1, originalAcceleratedBodyName );

                if( currentTestedPartial != NULL )
                {
                    existingAccelerationPartial = currentTestedPartial;
                    numberOfFeasibleAccelerations++;
                }
            }
        }
    }

    if( numberOfFeasibleAccelerations != 1 )
    {
        std::cerr<<"Error when finding scaled acceleration partial, found "<<numberOfFeasibleAccelerations<<" feasible options for "<<
                   "type "<<accelerationType<<" of "<<originalAcceleratingBodyName<<" on "<<originalAcceleratedBodyName<<std::endl;
    }

    if( existingAccelerationPartial == NULL )
    {
        std::cerr<<"Error when finding scaled acceleration partial, found no feasible partials for "<<
                   "type "<<accelerationType<<" of "<<originalAcceleratingBodyName<<" on "<<originalAcceleratedBodyName<<std::endl;
    }

    return existingAccelerationPartial;
}


//! Function to create a list of objects that can be used to compute partials of tidal gravity field variations
std::vector< std::shared_ptr< orbit_determination::TidalLoveNumberPartialInterface > > createTidalLoveNumberInterfaces(
        const NamedBodyMap& bodyMap,
        const std::string& acceleratingBodyName )
{
    // Create return map.
    std::vector< std::shared_ptr< orbit_determination::TidalLoveNumberPartialInterface > > loveNumberInterfaces;

    // Check if any gravity field variations are present
    if( bodyMap.at( acceleratingBodyName )->getGravityFieldVariationSet( ) != nullptr )
    {
        // Get list of tidal gravity field variations.
        std::vector< std::shared_ptr< gravitation::BasicSolidBodyTideGravityFieldVariations > >  variationObjectList =
                utilities::dynamicCastSVectorToTVector< gravitation::GravityFieldVariations,
                gravitation::BasicSolidBodyTideGravityFieldVariations >(
                    bodyMap.at( acceleratingBodyName )->getGravityFieldVariationSet( )->
                    getDirectTidalGravityFieldVariations( ) );

        // Create partial for each gravity field variation objet
        if( variationObjectList.size( ) > 0 )
        {
            // Get state/rotation functions for deformed body
            std::function< Eigen::Vector3d( ) > deformedBodyPositionFunction =
                    std::bind( &Body::getPosition, bodyMap.at( acceleratingBodyName ) );
            std::function< Eigen::Quaterniond( ) > rotationToDeformedBodyFrameFrameFunction =
                    std::bind( &Body::getCurrentRotationToLocalFrame, bodyMap.at( acceleratingBodyName ) );

            for( unsigned int i = 0; i < variationObjectList.size( ); i++ )
            {
                if( variationObjectList.at( i ) != nullptr )
                {
                    // Get state/rotation functions for deforming bodyies
                    std::vector< std::function< Eigen::Vector3d( ) > > deformingBodyStateFunctions;
                    std::vector< std::string > deformingBodies = variationObjectList.at( i )->getDeformingBodies( );
                    for( unsigned int i = 0; i < deformingBodies.size( ); i++ )
                    {
                        deformingBodyStateFunctions.push_back(
                                    std::bind( &Body::getPosition, bodyMap.at( deformingBodies.at( i ) ) ) );
                    }
                    // Get state/rotation functions for deformed body

                    // Create partial object
                    loveNumberInterfaces.push_back(
                                std::make_shared< orbit_determination::TidalLoveNumberPartialInterface >(
                                    variationObjectList.at( i ),
                                    deformedBodyPositionFunction,
                                    deformingBodyStateFunctions,
                                    rotationToDeformedBodyFrameFrameFunction,
                                    acceleratingBodyName ) );
                }
            }
        }
    }
    return loveNumberInterfaces;
}

template std::shared_ptr< acceleration_partials::AccelerationPartial > createAnalyticalAccelerationPartial< double >(
        std::shared_ptr< basic_astrodynamics::AccelerationModel< Eigen::Vector3d > > accelerationModel,
        const std::pair< std::string, std::shared_ptr< simulation_setup::Body > > acceleratedBody,
        const std::pair< std::string, std::shared_ptr< simulation_setup::Body > > acceleratingBody,
        const simulation_setup::NamedBodyMap& bodyMap,
        const std::shared_ptr< estimatable_parameters::EstimatableParameterSet< double > >
        parametersToEstimate,
        const std::map< std::string, std::map< std::string,
        std::vector< std::shared_ptr< acceleration_partials::AccelerationPartial > > > >& accelerationPartialsMap  );

#if( BUILD_EXTENDED_PRECISION_PROPAGATION_TOOLS )
template std::shared_ptr< acceleration_partials::AccelerationPartial > createAnalyticalAccelerationPartial< long double >(
        std::shared_ptr< basic_astrodynamics::AccelerationModel< Eigen::Vector3d > > accelerationModel,
        const std::pair< std::string, std::shared_ptr< simulation_setup::Body > > acceleratedBody,
        const std::pair< std::string, std::shared_ptr< simulation_setup::Body > > acceleratingBody,
        const simulation_setup::NamedBodyMap& bodyMap,
        const std::shared_ptr< estimatable_parameters::EstimatableParameterSet< long double > >
        parametersToEstimate,
        const std::map< std::string, std::map< std::string,
        std::vector< std::shared_ptr< acceleration_partials::AccelerationPartial > > > >& accelerationPartialsMap  );
#endif

template orbit_determination::StateDerivativePartialsMap createAccelerationPartialsMap< double >(
        const basic_astrodynamics::AccelerationMap& accelerationMap,
        const simulation_setup::NamedBodyMap& bodyMap,
        const std::shared_ptr< estimatable_parameters::EstimatableParameterSet< double > >
        parametersToEstimate );
        template orbit_determination::StateDerivativePartialsMap createAccelerationPartialsMap< long double >(
                const basic_astrodynamics::AccelerationMap& accelerationMap,
                const simulation_setup::NamedBodyMap& bodyMap,
                const std::shared_ptr< estimatable_parameters::EstimatableParameterSet< long double > >
                parametersToEstimate );

}

}
