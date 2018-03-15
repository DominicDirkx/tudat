/*    Copyright (c) 2010-2018, Delft University of Technology
 *    All rigths reserved
 *
 *    This file is part of the Tudat. Redistribution and use in source and
 *    binary forms, with or without modification, are permitted exclusively
 *    under the terms of the Modified BSD license. You should have received
 *    a copy of the license with this file. If not, please or visit:
 *    http://tudat.tudelft.nl/LICENSE.
 */

#ifndef TUDAT_DOPPLERMODELCONVERSIONS_H
#define TUDAT_DOPPLERMODELCONVERSIONS_H

#include "Tudat/Astrodynamics/ObservationModels/oneWayDopplerObservationModel.h"
#include "Tudat/Mathematics/NumericalQuadrature/createNumericalQuadrature.h"

namespace tudat
{

namespace observation_models
{

template< typename ObservableScalarType, typename TimeType >
std::map< TimeType, ObservableScalarType > convertOpenLoopToClosedLoopDopplerData(
        const boost::shared_ptr< OneWayDopplerObservationModel< ObservableScalarType, TimeType > > oneWayDopplerModel,
        boost::shared_ptr< numerical_quadrature::NumericalQuadratureSettings > quadratureSettings,
        const LinkEndType linkEndAssociatedWithTime,
        const TimeType startTime,
        const TimeType endTime,
        const TimeType timeStep )
{
    std::map< TimeType, ObservableScalarType > syntheticClosedLoopDopplerData;

    boost::function< ObservableScalarType( const TimeType ) > oneWayDopplerFunction =
            boost::bind( &OneWayDopplerObservationModel< ObservableScalarType, TimeType >::computeObservationEntry,
                         oneWayDopplerModel, _1, linkEndAssociatedWithTime, 0 );

    std::vector< TimeType > currentObservableTimeLimits;
    currentObservableTimeLimits.resize( 2 );
    currentObservableTimeLimits[ 0 ] = startTime - timeStep / 2.0;
    currentObservableTimeLimits[ 1 ] = startTime + timeStep / 2.0;

    boost::shared_ptr< numerical_quadrature::NumericalQuadrature< TimeType, ObservableScalarType > > quadrature =
            numerical_quadrature::createNumericalQuadrature(
                currentObservableTimeLimits, oneWayDopplerFunction, quadratureSettings );

    TimeType currentTime = startTime;
    while( currentTime <= endTime )
    {
        currentObservableTimeLimits[ 0 ] = currentTime - timeStep;
        currentObservableTimeLimits[ 1 ] = currentTime;

        quadrature->resetData(
                    currentObservableTimeLimits, oneWayDopplerFunction );
        syntheticClosedLoopDopplerData[ currentTime ] = quadrature->getQuadrature( ) / static_cast< ObservableScalarType >(
                    timeStep ) * physical_constants::getSpeedOfLight< ObservableScalarType >( );
        currentTime += timeStep;
    }

    return syntheticClosedLoopDopplerData;
}


}

}

#endif // TUDAT_DOPPLERMODELCONVERSIONS_H
