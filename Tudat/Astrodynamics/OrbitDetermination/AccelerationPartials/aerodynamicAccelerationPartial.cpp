/*    Copyright (c) 2010-2018, Delft University of Technology
 *    All rigths reserved
 *
 *    This file is part of the Tudat. Redistribution and use in source and
 *    binary forms, with or without modification, are permitted exclusively
 *    under the terms of the Modified BSD license. You should have received
 *    a copy of the license with this file. If not, please or visit:
 *    http://tudat.tudelft.nl/LICENSE.
 */

#include "Tudat/Astrodynamics/OrbitDetermination/AccelerationPartials/aerodynamicAccelerationPartial.h"

namespace tudat
{

namespace acceleration_partials
{


//! Function for updating partial w.r.t. the bodies' positions
void AerodynamicAccelerationPartial::update( const double currentTime )
{
    nominalState_ = vehicleStateGetFunction_( );

    // Compute state partial by numerical difference
    for( unsigned int i = 0; i < 6; i++ )
    {
        // Perturb state upwards
        perturbedState_ = nominalState_;
        perturbedState_( i ) += bodyStatePerturbations_( i );

        // Update environment/acceleration to perturbed state.
        flightConditions_->resetCurrentTime( TUDAT_NAN );
        aerodynamicAcceleration_->resetTime( TUDAT_NAN );
        vehicleStateSetFunction_( perturbedState_ );
        flightConditions_->updateConditions( currentTime );
        aerodynamicAcceleration_->updateMembers( currentTime );

        // Retrieve perturbed acceleration.
        upperturbedAcceleration_ = aerodynamicAcceleration_->getAcceleration( );

        // Perturb state downwards
        perturbedState_ = nominalState_;
        perturbedState_( i ) -= bodyStatePerturbations_( i );

        // Update environment/acceleration to perturbed state.
        flightConditions_->resetCurrentTime( TUDAT_NAN );
        aerodynamicAcceleration_->resetTime( TUDAT_NAN );
        vehicleStateSetFunction_( perturbedState_ );
        flightConditions_->updateConditions( currentTime );
        aerodynamicAcceleration_->updateMembers( currentTime );

        // Retrieve perturbed acceleration.
        downperturbedAcceleration_ = aerodynamicAcceleration_->getAcceleration( );

        // Compute partial
        currentAccelerationStatePartials_.block( 0, i, 3, 1 ) =
                ( upperturbedAcceleration_ - downperturbedAcceleration_ ) / ( 2.0 * bodyStatePerturbations_( i ) );
    }

    // Reset environment/acceleration mode to nominal conditions
    flightConditions_->resetCurrentTime( TUDAT_NAN );
    aerodynamicAcceleration_->resetTime( TUDAT_NAN );

    vehicleStateSetFunction_( nominalState_ );
    flightConditions_->updateConditions( currentTime );
    aerodynamicAcceleration_->updateMembers( currentTime );
}


} // namespace acceleration_partials

} // namespace tudat
