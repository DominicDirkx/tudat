/*    Copyright (c) 2010-2017, Delft University of Technology
 *    All rigths reserved
 *
 *    This file is part of the Tudat. Redistribution and use in source and
 *    binary forms, with or without modification, are permitted exclusively
 *    under the terms of the Modified BSD license. You should have received
 *    a copy of the license with this file. If not, please or visit:
 *    http://tudat.tudelft.nl/LICENSE.
 */

#include "Tudat/Astrodynamics/Ephemerides/compositeEphemeris.h"

namespace tudat
{

namespace ephemerides
{

//! Function to check whether an ephemeris is a (type of) tabulated ephemeris
bool isCompositeEphemeris( const boost::shared_ptr< Ephemeris > ephemeris )
{
    bool objectIsTabulated = 0;
    if( ( boost::dynamic_pointer_cast< CompositeEphemeris< double, double > >( ephemeris ) != NULL ) ||
            ( boost::dynamic_pointer_cast< CompositeEphemeris< long double, double > >( ephemeris ) != NULL ) ||
            ( boost::dynamic_pointer_cast< CompositeEphemeris< long double, Time > >( ephemeris ) != NULL ) ||
            ( boost::dynamic_pointer_cast< CompositeEphemeris< double, Time > >( ephemeris ) != NULL ) )
    {
        objectIsTabulated = 1;
    }
    return objectIsTabulated;
}

void clearCompositeEphemerisContents(
        const boost::shared_ptr< Ephemeris > ephemerisToClear )
{
//    if( isTabulatedEphemeris( ephemerisToClear ) )
//    {
//        if( ( boost::dynamic_pointer_cast< CompositeEphemeris< double, double > >( ephemeris ) != NULL ) )
//        {
//             boost::shared_ptr< interpolators::OneDimensionalInterpolator< double, Eigen::Matrix< double, 6, 1 > > >
//                     stateInterpolator = boost::dynamic_pointer_cast< CompositeEphemeris< double, double > >(
//                         ephemeris )->getInterpolator( );
//        }
//        else if( boost::dynamic_pointer_cast< CompositeEphemeris< long double, double > >( ephemeris ) != NULL )
//        {
//            boost::shared_ptr< interpolators::OneDimensionalInterpolator< double, Eigen::Matrix< long double, 6, 1 > > >
//                    stateInterpolator = boost::dynamic_pointer_cast< CompositeEphemeris< long double, double > >(
//                        ephemeris )->getInterpolator( );
//        }
//        else if( boost::dynamic_pointer_cast< CompositeEphemeris< long double, Time > >( ephemeris ) != NULL )
//        {
//            boost::shared_ptr< interpolators::OneDimensionalInterpolator< Time, Eigen::Matrix< long double, 6, 1 > > >
//                    stateInterpolator = boost::dynamic_pointer_cast< CompositeEphemeris< long double, Time > >(
//                        ephemeris )->getInterpolator( );
//        }
//        else if( boost::dynamic_pointer_cast< CompositeEphemeris< double, Time > >( ephemeris ) != NULL )
//        {
//            boost::shared_ptr< interpolators::OneDimensionalInterpolator< Time, Eigen::Matrix< double, 6, 1 > > >
//                    stateInterpolator = boost::dynamic_pointer_cast< CompositeEphemeris< double, Time > >(
//                        ephemeris )->getInterpolator( );
//        }
//    }
}

}

}
