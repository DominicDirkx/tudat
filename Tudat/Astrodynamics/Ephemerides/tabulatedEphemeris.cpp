/*    Copyright (c) 2010-2017, Delft University of Technology
 *    All rigths reserved
 *
 *    This file is part of the Tudat. Redistribution and use in source and
 *    binary forms, with or without modification, are permitted exclusively
 *    under the terms of the Modified BSD license. You should have received
 *    a copy of the license with this file. If not, please or visit:
 *    http://tudat.tudelft.nl/LICENSE.
 */
#include "Tudat/Astrodynamics/Ephemerides/tabulatedEphemeris.h"

namespace tudat
{

namespace ephemerides
{

//! Get cartesian state from ephemeris (in double precision), for double StateScalarType
template< >
Eigen::Vector6d TabulatedCartesianEphemeris< double, double >::getCartesianState(
        const double ephemerisTime)
{
    return interpolator_->interpolate( ephemerisTime );
}

//! Get cartesian state from ephemeris (in long double precision), for double StateScalarType
template< >
Eigen::Matrix< long double, 6, 1 > TabulatedCartesianEphemeris< double, double >::getCartesianLongState(
        const double secondsSinceEpoch )
{
    return interpolator_->interpolate( secondsSinceEpoch ).cast< long double >( );
}

//! Get cartesian state from ephemeris (in double precision from Time input), for double StateScalarType
template< >
Eigen::Vector6d TabulatedCartesianEphemeris< double, double >::getCartesianStateFromExtendedTime(
        const Time& time )
{
    return interpolator_->interpolate( time.getSeconds< double >( ) );
}

//! Get cartesian state from ephemeris (in long double precision from Time input), for double StateScalarType
template< >
Eigen::Matrix< long double, 6, 1 > TabulatedCartesianEphemeris< double, double >::getCartesianLongStateFromExtendedTime(
        const Time& time )
{
    return interpolator_->interpolate( time.getSeconds< double >( ) ).cast< long double >( );
}





//! Get cartesian state from ephemeris (in double precision), for long double StateScalarType
template< >
Eigen::Vector6d TabulatedCartesianEphemeris< long double, double >::getCartesianState(
        const double ephemerisTime )
{
    return interpolator_->interpolate( ephemerisTime ).cast< double >( );
}

//! Get cartesian state from ephemeris (in long double precision), for long double StateScalarType
template< >
Eigen::Matrix< long double, 6, 1 > TabulatedCartesianEphemeris< long double, double >::getCartesianLongState(
        const double secondsSinceEpoch )
{
    return interpolator_->interpolate( secondsSinceEpoch );
}

//! Get cartesian state from ephemeris (in double precision from Time input), for double StateScalarType
template< >
Eigen::Vector6d TabulatedCartesianEphemeris< long double, double >::getCartesianStateFromExtendedTime(
        const Time& time )
{
    return interpolator_->interpolate( time.getSeconds< double >( ) ).cast< double >( );
}

//! Get cartesian state from ephemeris (in long double precision from Time input), for double StateScalarType
template< >
Eigen::Matrix< long double, 6, 1 > TabulatedCartesianEphemeris< long double, double >::getCartesianLongStateFromExtendedTime(
        const Time& time )
{
    return interpolator_->interpolate( time.getSeconds< double >( ) );
}






//! Get cartesian state from ephemeris (in double precision), for long double StateScalarType
template< >
Eigen::Vector6d TabulatedCartesianEphemeris< long double, Time >::getCartesianState(
        const double ephemerisTime )
{
    return interpolator_->interpolate( Time( ephemerisTime ) ).cast< double >( );
}

//! Get cartesian state from ephemeris (in long double precision), for long double StateScalarType
template< >
Eigen::Matrix< long double, 6, 1 > TabulatedCartesianEphemeris< long double, Time >::getCartesianLongState(
        const double secondsSinceEpoch )
{
    return interpolator_->interpolate( Time( secondsSinceEpoch ) );
}

//! Get cartesian state from ephemeris (in double precision from Time input).
template< >
Eigen::Vector6d TabulatedCartesianEphemeris< long double, Time >::getCartesianStateFromExtendedTime(
        const Time& time )
{
    return interpolator_->interpolate( time ).cast< double >( );
}

//! Get cartesian state from ephemeris (in long double precision from Time input).
template< >
Eigen::Matrix< long double, 6, 1 > TabulatedCartesianEphemeris< long double, Time >::getCartesianLongStateFromExtendedTime(
        const Time& time )
{
    return interpolator_->interpolate( time );
}


//! Function to check whether an ephemeris is a (type of) tabulated ephemeris
bool isTabulatedEphemeris( const boost::shared_ptr< Ephemeris > ephemeris )
{
    bool objectIsTabulated = 0;
    if( ( boost::dynamic_pointer_cast< TabulatedCartesianEphemeris< double, double > >( ephemeris ) != NULL ) ||
            ( boost::dynamic_pointer_cast< TabulatedCartesianEphemeris< long double, double > >( ephemeris ) != NULL ) ||
            ( boost::dynamic_pointer_cast< TabulatedCartesianEphemeris< long double, Time > >( ephemeris ) != NULL ) ||
            ( boost::dynamic_pointer_cast< TabulatedCartesianEphemeris< double, Time > >( ephemeris ) != NULL ) )
    {
        objectIsTabulated = 1;
    }
    return objectIsTabulated;
}

void clearTabulatedEphemerisContents(
        const boost::shared_ptr< Ephemeris > ephemerisToClear )
{
    if( isTabulatedEphemeris( ephemerisToClear ) )
    {
        if( ( boost::dynamic_pointer_cast< TabulatedCartesianEphemeris< double, double > >( ephemerisToClear ) != NULL ) )
        {
             boost::shared_ptr< interpolators::OneDimensionalInterpolator< double, Eigen::Matrix< double, 6, 1 > > >
                     stateInterpolator = boost::dynamic_pointer_cast< TabulatedCartesianEphemeris< double, double > >(
                         ephemerisToClear )->getInterpolator( );
             stateInterpolator->clearInterpolator( );
        }
        else if( boost::dynamic_pointer_cast< TabulatedCartesianEphemeris< long double, double > >( ephemerisToClear ) != NULL )
        {
            boost::shared_ptr< interpolators::OneDimensionalInterpolator< double, Eigen::Matrix< long double, 6, 1 > > >
                    stateInterpolator = boost::dynamic_pointer_cast< TabulatedCartesianEphemeris< long double, double > >(
                        ephemerisToClear )->getInterpolator( );
            stateInterpolator->clearInterpolator( );
        }
        else if( boost::dynamic_pointer_cast< TabulatedCartesianEphemeris< long double, Time > >( ephemerisToClear ) != NULL )
        {
            boost::shared_ptr< interpolators::OneDimensionalInterpolator< Time, Eigen::Matrix< long double, 6, 1 > > >
                    stateInterpolator = boost::dynamic_pointer_cast< TabulatedCartesianEphemeris< long double, Time > >(
                        ephemerisToClear )->getInterpolator( );
            stateInterpolator->clearInterpolator( );
        }
        else if( boost::dynamic_pointer_cast< TabulatedCartesianEphemeris< double, Time > >( ephemerisToClear ) != NULL )
        {
            boost::shared_ptr< interpolators::OneDimensionalInterpolator< Time, Eigen::Matrix< double, 6, 1 > > >
                    stateInterpolator = boost::dynamic_pointer_cast< TabulatedCartesianEphemeris<  double, Time > >(
                        ephemerisToClear )->getInterpolator( );
            stateInterpolator->clearInterpolator( );
        }
    }
}

} // namespace ephemerides

} // namespace tudat

