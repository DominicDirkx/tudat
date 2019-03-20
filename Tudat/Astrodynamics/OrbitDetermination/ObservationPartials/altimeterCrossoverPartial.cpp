/*    Copyright (c) 2010-2018, Delft University of Technology
 *    All rigths reserved
 *
 *    This file is part of the Tudat. Redistribution and use in source and
 *    binary forms, with or without modification, are permitted exclusively
 *    under the terms of the Modified BSD license. You should have received
 *    a copy of the license with this file. If not, please or visit:
 *    http://tudat.tudelft.nl/LICENSE.
 */

#include "Tudat/Astrodynamics/OrbitDetermination/ObservationPartials/altimeterCrossoverPartial.h"

namespace tudat
{

namespace observation_partials
{


////! Update the scaling object to the current times and states
//void AltimeterCrossoverScaling::update( const std::vector< Eigen::Vector6d >& linkEndStates,
//                                 const std::vector< double >& times,
//                                 const observation_models::LinkEndType fixedLinkEnd,
//                                 const Eigen::VectorXd currentObservation )
//{
//    // Compute Euclidean distance vector
//    Eigen::Vector3d rangeVector = linkEndStates[ 1 ].segment( 0, 3 ) - linkEndStates[ 0 ].segment( 0, 3 );
//    Eigen::Matrix< double, 1, 3 > rangeVectorNormalized = rangeVector.transpose( ) / rangeVector.norm( );

//    // Compute scaling for receiver reference
//    if( fixedLinkEnd == observation_models::receiver )
//    {
//        referenceLightTimeCorrectionScaling_ = 1.0 / ( 1.0 - rangeVectorNormalized.transpose( ).dot( linkEndStates[ 0 ].segment( 3, 3 ) ) /
//                physical_constants::SPEED_OF_LIGHT );
//        referenceScalingFactor_ =  rangeVectorNormalized * referenceLightTimeCorrectionScaling_;
//    }

//    // Compute scaling for transmitter reference
//    else if( fixedLinkEnd == observation_models::transmitter )
//    {
//        referenceLightTimeCorrectionScaling_ =
//                1.0 / ( 1.0 - rangeVectorNormalized.transpose( ).dot( linkEndStates[ 1 ].segment( 3, 3 ) ) /
//                physical_constants::SPEED_OF_LIGHT );
//        referenceScalingFactor_ =  rangeVectorNormalized * referenceLightTimeCorrectionScaling_;
//    }

//    currentLinkEndType_ = fixedLinkEnd;
//}

//! Function to calculate the observation partial(s) at required time and state
AltimeterCrossoverPartial::AltimeterCrossoverPartialReturnType AltimeterCrossoverPartial::calculatePartial(
        const std::vector< Eigen::Vector6d >& states,
        const std::vector< double >& times,
        const observation_models::LinkEndType linkEndOfFixedTime,
        const Eigen::Vector1d& currentObservation )
{
//    std::cout << "You've made it to altimeterCrossoverPartial.cpp!" << std::endl;
    // std::vector< std::pair< Eigen::Matrix< double, 1, Eigen::Dynamic >, double > >
    AltimeterCrossoverPartialReturnType returnPartial;

    Eigen::Matrix< double, 3, 1 > firstArcPartialWrtCurrentPosition;
    Eigen::Matrix< double, 3, 1 > secondArcPartialWrtCurrentPosition;
    Eigen::Matrix< double, 3, 1 > observationPartialWrtCurrentPosition;

    // Iterate over all link ends
    for( positionPartialIterator_ = positionPartialList_.begin( ); positionPartialIterator_ != positionPartialList_.end( );
         positionPartialIterator_++ )
    {
//        observationPartialWrtCurrentPosition = TUDAT_NAN;

        // The current partial relates to the state at arc 1.
        if( positionPartialIterator_->first == observation_models::first_arc_body )
        {
            currentState_  = states[ 0 ];
            currentTime_ = times[ 0 ];

            if(currentTime_ == 1045478283.708825 )
            {
                std::cout << std::setprecision(17) << std::endl << "for t1: " << currentTime_ << std::endl;
                std::cout << "the state vector x(t1) is (from altimeterCrossoverPartial.cpp): \n" << currentState_.transpose() << std::endl;
            }
            Eigen::Matrix< double, 3, Eigen::Dynamic > currentInertialPositionPartialWrtParameter =
                    positionPartialIterator_->second->calculatePartialOfPosition(
                                          currentState_ , currentTime_ );
            double rho = currentState_.segment( 0, 3 ).norm( );
            firstArcPartialWrtCurrentPosition << ( currentInertialPositionPartialWrtParameter *
                                                   ( (1/rho)*currentState_.segment( 0, 3 ) ) );
            // beware the MINUS!
            observationPartialWrtCurrentPosition = - ( currentInertialPositionPartialWrtParameter *
                                                    firstArcPartialWrtCurrentPosition );
//            if( currentTime_ == 1045468254.9337766 )
//            {
//                std::cout << std::setprecision(17) << "\n state first_arc_body says: \n" << currentState_ << std::endl;
//                std::cout << std::setprecision(17) << "Partials first_arc_body says: \n" << observationPartialWrtCurrentPosition << std::endl;
//            }
            returnPartial.push_back(
                        std::make_pair( observationPartialWrtCurrentPosition, currentTime_ ) );
        }
        // The current partial relates to the state at arc 2.
        else if( positionPartialIterator_->first == observation_models::second_arc_body )
        {
            currentState_  = states[ 1 ];
            currentTime_ = times[ 1 ];

            if(currentTime_ == 1045489321.9936566 )
            {
                std::cout << std::setprecision(17) << std::endl << "for t2: " << currentTime_ << std::endl;
                std::cout << "the state vector x(t2) is (from altimeterCrossoverPartial.cpp): \n" << currentState_.transpose() << std::endl;
            }
            Eigen::Matrix< double, 3, Eigen::Dynamic > currentInertialPositionPartialWrtParameter =
                    positionPartialIterator_->second->calculatePartialOfPosition(
                                          currentState_ , currentTime_ );
            double rho = currentState_.segment( 0, 3 ).norm( );
            secondArcPartialWrtCurrentPosition << ( currentInertialPositionPartialWrtParameter *
                                                    ( (1/rho)*currentState_.segment( 0, 3 ) ) );
            observationPartialWrtCurrentPosition = ( currentInertialPositionPartialWrtParameter *
                                                    secondArcPartialWrtCurrentPosition );
//            if( currentTime_ == 1045731864.6734738 )
//            {
//                std::cout << std::setprecision(17) << "\n state second_arc_body says: \n" << currentState_ << std::endl;
//                std::cout << std::setprecision(17) << "Partials second_arc_body says: \n" << observationPartialWrtCurrentPosition << std::endl;
//            }
            returnPartial.push_back(
                        std::make_pair( observationPartialWrtCurrentPosition, currentTime_ ) );

        }
    }

//        observationPartialWrtCurrentPosition << 1, 1, 1;
//        observationPartialWrtCurrentPosition << 0, 0, 0;
//    observationPartialWrtCurrentPosition << ( secondArcPartialWrtCurrentPosition -
//                                              firstArcPartialWrtCurrentPosition);
    // Set partial output
//    returnPartial.push_back(
//                std::make_pair( observationPartialWrtCurrentPosition, times[ 0 ] ) );

//    if( times[ 0 ] == 1045468254.9337766 || times[ 1 ] == 1045731864.6734738 )
//    {
//        std::cout << std::setprecision(17)<<times[ 0 ] << std::endl;
//        std::cout << std::setprecision(17)<<times[ 1 ] << std::endl;
//        std::cout << returnPartial.size() << std::endl;
//        std::cout << std::setprecision(17)<< returnPartial[0].first << std::endl;
//        std::cout << std::setprecision(17)<< returnPartial[1].first << std::endl;
//    }

    return returnPartial;
}

} // namespace observation_partials

} // namespace tudat
