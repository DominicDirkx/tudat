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
#include "Tudat/Mathematics/NumericalQuadrature/gaussianQuadrature.h"
#include "Tudat/Mathematics/NumericalQuadrature/simpsonsQuadrature.h"
#include "Tudat/Mathematics/NumericalQuadrature/trapezoidQuadrature.h"

namespace tudat
{

namespace observation_models
{

template< typename ObservableScalarType, typename TimeType >
std::map< TimeType, ObservableScalarType > convertOpenLoopToClosedLoopDopplerData(
        const boost::shared_ptr< OneWayDopplerObservationModel< ObservableScalarType, TimeType > > oneWayDopplerModel,
        const LinkEndType linkEndAssociatedWithTime,
        const TimeType startTime,
        const TimeType endTime,
        const TimeType timeStep )
{
    std::cout<<"Link end to use "<<linkEndAssociatedWithTime<<std::endl;
    std::map< TimeType, ObservableScalarType > syntheticClosedLoopDopplerData;

    boost::function< ObservableScalarType( const TimeType ) > oneWayDopplerFunction =
            boost::bind( &OneWayDopplerObservationModel< ObservableScalarType, TimeType >::computeObservationEntry,
                         oneWayDopplerModel, _1, linkEndAssociatedWithTime, 0 );

    boost::shared_ptr< numerical_quadrature::SimpsonNumericalQuadrature< TimeType, ObservableScalarType > > quadrature =
            boost::make_shared< numerical_quadrature::SimpsonNumericalQuadrature< TimeType, ObservableScalarType > >( );

    std::vector< TimeType > currentObservableTimeLimits;
    currentObservableTimeLimits.resize( 2 );
    TimeType currentTime = startTime;
    while( currentTime <= endTime )
    {
        currentObservableTimeLimits[ 0 ] = currentTime - timeStep / 2.0;
        currentObservableTimeLimits[ 1 ] = currentTime + timeStep / 2.0;

        quadrature->resetData(
                    currentObservableTimeLimits, oneWayDopplerFunction );
        syntheticClosedLoopDopplerData[ currentTime ] = quadrature->getQuadrature( ) / timeStep;
        currentTime += timeStep;
    }

    return syntheticClosedLoopDopplerData;
}


}

}

#endif // TUDAT_DOPPLERMODELCONVERSIONS_H
