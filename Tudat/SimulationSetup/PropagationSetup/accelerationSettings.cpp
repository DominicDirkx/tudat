#include "Tudat/SimulationSetup/PropagationSetup/accelerationSettings.h"

namespace tudat
{

namespace simulation_setup
{

//! Function to get the full list of combinations of degrees/orders of two bodies for full two-body potential
std::vector< boost::tuple< unsigned int, unsigned int, unsigned int, unsigned int > > getExtendedSinglePointMassInteractions(
        const int maximumDegreeOfBodyUndergoingAcceleration,
        const int maximumOrderOfBodyUndergoingAcceleration,
        const int maximumDegreeOfBodyExertingAcceleration,
        const int maximumOrderOfBodyExertingAcceleration,
        const bool includePointMass )
{
    std::vector< boost::tuple< unsigned int, unsigned int, unsigned int, unsigned int > > coefficientCombinationsToUse;

    for( int i = 1; i <= maximumDegreeOfBodyUndergoingAcceleration; i++ )
    {
        for( int j = 0; ( j <= maximumOrderOfBodyUndergoingAcceleration && j <= i ); j++ )
        {
            coefficientCombinationsToUse.push_back( boost::make_tuple( i, j, 0, 0 ) );
        }
    }

    for( int k = 1; k <= maximumDegreeOfBodyExertingAcceleration; k++ )
    {
        for( int l = 0; ( l <= maximumOrderOfBodyExertingAcceleration && l <= k ); l++ )
        {
            coefficientCombinationsToUse.push_back( boost::make_tuple( 0, 0, k, l ) );
        }
    }

    if( includePointMass )
    {
        coefficientCombinationsToUse.push_back( boost::make_tuple( 0, 0, 0, 0 ) );
    }

    return coefficientCombinationsToUse;
}

//! Constructor
MutualExtendedBodySphericalHarmonicAccelerationSettings::MutualExtendedBodySphericalHarmonicAccelerationSettings(
        const int maximumDegreeOfBodyUndergoingAcceleration,
        const int maximumOrderOfBodyUndergoingAcceleration,
        const int maximumDegreeOfBodyExertingAcceleration,
        const int maximumOrderOfBodyExertingAcceleration ):
    AccelerationSettings( basic_astrodynamics::mutual_extended_body_spherical_harmonic_gravity ),
    maximumDegreeOfBody1_( maximumDegreeOfBodyUndergoingAcceleration ),
    maximumDegreeOfBody2_( maximumDegreeOfBodyExertingAcceleration )
{
    for( int i = 0; i <= maximumDegreeOfBodyUndergoingAcceleration; i++ )
    {
        for( int j = 0; ( j <= maximumOrderOfBodyUndergoingAcceleration && j <= i ); j++ )
        {
            for( int k = 0; k <= maximumDegreeOfBodyExertingAcceleration; k++ )
            {
                for( int l = 0; ( l <= maximumOrderOfBodyExertingAcceleration && l <= k ); l++ )
                {
                    coefficientCombinationsToUse_.push_back( boost::make_tuple( i, j, k, l ) );
                }
            }
        }
    }
}

//! Constructor
MutualExtendedBodySphericalHarmonicAccelerationSettings::MutualExtendedBodySphericalHarmonicAccelerationSettings(
       const std::vector< boost::tuple< unsigned int, unsigned int, unsigned int, unsigned int > >&
        coefficientCombinationsToUse ):
    AccelerationSettings( basic_astrodynamics::mutual_extended_body_spherical_harmonic_gravity ),
    coefficientCombinationsToUse_( coefficientCombinationsToUse )
{

    maximumDegreeOfBody1_ = 0;
    maximumDegreeOfBody2_ = 0;

    int degreeOfBody1, degreeOfBody2;
    for( unsigned int i = 0; i < coefficientCombinationsToUse_.size( ); i++ )
    {
        degreeOfBody1 = coefficientCombinationsToUse.at( i ).get< 0 >( );
        degreeOfBody2 = coefficientCombinationsToUse.at( i ).get< 2 >( );

        if( degreeOfBody1 > maximumDegreeOfBody1_ )
        {
            maximumDegreeOfBody1_ = degreeOfBody1;
        }

        if( degreeOfBody2 > maximumDegreeOfBody2_ )
        {
            maximumDegreeOfBody2_ = degreeOfBody2;
        }
    }
}

} // namespace simulation_setup

} // namespace tudat
