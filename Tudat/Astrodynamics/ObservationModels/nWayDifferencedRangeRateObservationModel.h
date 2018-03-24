/*    Copyright (c) 2010-2018, Delft University of Technology
 *    All rigths reserved
 *
 *    This file is part of the Tudat. Redistribution and use in source and
 *    binary forms, with or without modification, are permitted exclusively
 *    under the terms of the Modified BSD license. You should have received
 *    a copy of the license with this file. If not, please or visit:
 *    http://tudat.tudelft.nl/LICENSE.
 */

#ifndef TUDAT_NWAYDIFFERENCEDRANGERATEOBSERVATIONMODEL_H
#define TUDAT_NWAYDIFFERENCEDRANGERATEOBSERVATIONMODEL_H

#include <map>
#include <iomanip>

#include <boost/function.hpp>
#include <boost/make_shared.hpp>

#include <Eigen/Core>

#include "Tudat/Astrodynamics/BasicAstrodynamics/physicalConstants.h"
#include "Tudat/Astrodynamics/ObservationModels/observationModel.h"
#include "Tudat/Astrodynamics/ObservationModels/lightTimeSolution.h"
#include "Tudat/Astrodynamics/ObservationModels/nWayRangeObservationModel.h"

namespace tudat
{

namespace observation_models
{

//! Class for simulating n-way differenced range (e.g. closed-loop Doppler) observable
/*!
 *  Class for simulating n-way differenced range (e.g. closed-loop Doppler) observable. The observable is obtained by
 *  subtracting the range at two time intervals, and dividing by the time difference. It represents the time-averages value
 *  of the range-rate over this integration time.
 *  The user may add observation biases to model system-dependent deviations between measured and true observation.
 */
template< typename ObservationScalarType = double, typename TimeType = double >
class NWayDifferencedRangeObservationModel: public ObservationModel< 1, ObservationScalarType, TimeType >
{
public:

    typedef Eigen::Matrix< ObservationScalarType, 6, 1 > StateType;
    typedef Eigen::Matrix< ObservationScalarType, 3, 1 > PositionType;


    //! Constructor.
    /*!
     *  Constructor,
     *  \param arcStartLightTimeCalculator Object to compute the light-time (including any corrections w.r.t. Euclidean case)
     *  \param arcEndLightTimeCalculator Object to compute the light-time (including any corrections w.r.t. Euclidean case)
     *  \param observationBiasCalculator Object for calculating system-dependent errors in the
     *  observable, i.e. deviations from the physically ideal observable between reference points (default none).
     *  \param integrationTimeFunction Function returning the integration time of the observable as a function of the
     *  current observation time.
     */
    NWayDifferencedRangeObservationModel(
            const boost::shared_ptr< NWayRangeObservationModel< ObservationScalarType, TimeType > >
            nWayRangeModel,
            const boost::function< double( const double ) > integrationTimeFunction,
            const boost::shared_ptr< ObservationBias< 1 > > observationBiasCalculator = NULL ):
        ObservationModel< 1, ObservationScalarType, TimeType >( one_way_differenced_range, observationBiasCalculator ),
        nWayRangeModel_( nWayRangeModel ), integrationTimeFunction_( integrationTimeFunction )
    {
        numberOfLinkEnds_ = nWayRangeModel_->getNumberOfLinkEnds( );

        arcStartLinkEndTimes_.resize( numberOfLinkEnds_ );
        arcEndLinkEndTimes_.resize( numberOfLinkEnds_ );
        arcStartLinkEndStates_.resize( numberOfLinkEnds_ );
        arcEndLinkEndStates_.resize( numberOfLinkEnds_ );
    }

    //! Destructor
    ~NWayDifferencedRangeObservationModel( ){ }

    //! Function to compute ideal n-way differenced range observation without any corrections at given time.
    /*!
     *  This function compute ideal n-way differenced range without any corrections at a given time.
     *  The time argument can be either the reception or transmission time (defined by linkEndAssociatedWithTime input).
     *  It does not include system-dependent measurement
     *  errors, such as biases or clock errors.
     *  \param time Time at which observation is to be simulated
     *  \param linkEndAssociatedWithTime Link end at which given time is valid, i.e. link end for which associated time
     *  is kept constant (to input value)
     *  \return Calculated observed n-way differenced range
     */
    Eigen::Matrix< ObservationScalarType, 1, 1 > computeIdealObservations(
            const TimeType time,
            const LinkEndType linkEndAssociatedWithTime )

    {
        TimeType currentIntegrationTime = integrationTimeFunction_( time );

        ObservationScalarType arcEndRange = nWayRangeModel_->computeIdealObservations(
                    time, linkEndAssociatedWithTime )( 0 );
        ObservationScalarType arcStartRange = nWayRangeModel_->computeIdealObservations(
                    time - currentIntegrationTime, linkEndAssociatedWithTime )( 0 );

        return ( Eigen::Matrix< ObservationScalarType, 1, 1 >( ) << ( arcEndRange - arcStartRange ) /
                 static_cast< ObservationScalarType >( currentIntegrationTime ) ).finished( );
    }

    //! Function to compute n-way differenced range observable without any corrections.
    /*!
     *  Function to compute n-way differenced range  observable without any corrections, i.e. the true physical differenced
     *  range as computed from the defined link ends. It does not include system-dependent measurement
     *  errors, such as biases or clock errors.
     *  The times and states of the link ends are also returned in double precision. These states and times are returned by
     *  reference.
     *  \param time Time at which observable is to be evaluated.
     *  \param linkEndAssociatedWithTime Link end at which given time is valid, i.e. link end for which associated time
     *  is kept constant (to input value)
     *  \param linkEndTimes List of times at each link end during observation.
     *  \param linkEndStates List of states at each link end during observation.
     *  \return Ideal n-way differenced range observable.
     */
    Eigen::Matrix< ObservationScalarType, 1, 1 > computeIdealObservationsWithLinkEndData(
            const TimeType time,
            const LinkEndType linkEndAssociatedWithTime,
            std::vector< double >& linkEndTimes,
            std::vector< Eigen::Matrix< double, 6, 1 > >& linkEndStates )
    {
        TimeType currentIntegrationTime = integrationTimeFunction_( time );
        Eigen::Matrix< ObservationScalarType, 1, 1 > arcEndRange = nWayRangeModel_->computeIdealObservationsWithLinkEndData(
                    time, linkEndAssociatedWithTime, arcEndLinkEndTimes_, arcEndLinkEndStates_ );
        Eigen::Matrix< ObservationScalarType, 1, 1 > arcStartRange = nWayRangeModel_->computeIdealObservationsWithLinkEndData(
                    time - currentIntegrationTime, linkEndAssociatedWithTime, arcStartLinkEndTimes_, arcStartLinkEndStates_ );

        linkEndTimes.resize( 2 * numberOfLinkEnds_ );
        linkEndStates.resize( 2 * numberOfLinkEnds_ );

        for( int i = 0; i < numberOfLinkEnds_;i++ )
        {
            linkEndTimes[ i ] = arcStartLinkEndTimes_.at( i );
            linkEndTimes[ 2 * i ] = arcEndLinkEndTimes_.at( i );

            linkEndStates[ i ] = arcStartLinkEndStates_.at( i );
            linkEndStates[ 2 * i ] = arcStartLinkEndStates_.at( i );

        }

        return ( Eigen::Matrix< ObservationScalarType, 1, 1 >( ) << ( arcEndRange - arcStartRange ) /
                 static_cast< ObservationScalarType >( currentIntegrationTime ) ).finished( );
    }


private:

    boost::shared_ptr< NWayRangeObservationModel< ObservationScalarType, TimeType > > nWayRangeModel_;

    //! Function returning the integration time of the observable as a function of the current observation time.
    boost::function< double( const double ) > integrationTimeFunction_;

    int numberOfLinkEnds_;

    std::vector< double > arcStartLinkEndTimes_;

    std::vector< double > arcEndLinkEndTimes_;

    std::vector< Eigen::Matrix< double, 6, 1 > > arcStartLinkEndStates_;

    std::vector< Eigen::Matrix< double, 6, 1 > > arcEndLinkEndStates_;

};

}

}

#endif // TUDAT_NWAYDIFFERENCEDRANGERATEOBSERVATIONMODEL_H
