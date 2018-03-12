/*    Copyright (c) 2010-2018, Delft University of Technology
 *    All rigths reserved
 *
 *    This file is part of the Tudat. Redistribution and use in source and
 *    binary forms, with or without modification, are permitted exclusively
 *    under the terms of the Modified BSD license. You should have received
 *    a copy of the license with this file. If not, please or visit:
 *    http://tudat.tudelft.nl/LICENSE.
 */

#define BOOST_TEST_MAIN


#include <boost/test/unit_test.hpp>

#include "Tudat/Basics/testMacros.h"

#include "Tudat/Astrodynamics/ObservationModels/dopplerModelConversions.h"
#include "Tudat/External/SpiceInterface/spiceEphemeris.h"
#include "Tudat/SimulationSetup/tudatEstimationHeader.h"

namespace tudat
{
namespace unit_tests
{

using namespace tudat;
using namespace tudat::observation_models;
using namespace tudat::spice_interface;
using namespace tudat::ephemerides;
using namespace tudat::simulation_setup;
using namespace tudat::orbital_element_conversions;
using namespace tudat::coordinate_conversions;
using namespace tudat::unit_conversions;


BOOST_AUTO_TEST_SUITE( test_doppler_models )


BOOST_AUTO_TEST_CASE( testOneWayDoppplerModel )
{
    // Load Spice kernels
    std::string kernelsPath = input_output::getSpiceKernelPath( );

    spice_interface::loadStandardSpiceKernels( );
    spice_interface::loadSpiceKernelInTudat( kernelsPath + "MEX_ROB_130101_131231_001.bsp" );

    // Define bodies to use.
    std::vector< std::string > bodiesToCreate;
    bodiesToCreate.push_back( "Earth" );
    bodiesToCreate.push_back( "Sun" );
    bodiesToCreate.push_back( "Mars" );

    // Specify initial time
    double initialEphemerisTime = tudat::basic_astrodynamics::convertCalendarDateToJulianDaysSinceEpoch(
                2013, 12, 29, 0, 0, 0.0, tudat::basic_astrodynamics::JULIAN_DAY_ON_J2000 ) * 86400.0;
    double finalEphemerisTime = initialEphemerisTime + 7.0 * 86400.0;
    double maximumTimeStep = 3600.0;
    double buffer = 10.0 * maximumTimeStep;

    // Create bodies settings needed in simulation
    std::map< std::string, boost::shared_ptr< BodySettings > > defaultBodySettings =
            getDefaultBodySettings(
                bodiesToCreate, initialEphemerisTime - buffer, finalEphemerisTime + buffer );

    // Create bodies
    NamedBodyMap bodyMap = createBodies( defaultBodySettings );
    bodyMap[ "MEX" ] = boost::make_shared< Body >( );
    bodyMap[ "MEX" ]->setEphemeris( boost::make_shared< tudat::ephemerides::SpiceEphemeris >(
                                        "MEX", "Mars", false, true, true ) );

    // Create ground station
    const Eigen::Vector3d stationCartesianPosition( 1917032.190, 6029782.349, -801376.113 );
    createGroundStation( bodyMap.at( "Earth" ), "Station1", stationCartesianPosition, cartesian_position );

    setGlobalFrameBodyEphemerides( bodyMap, "SSB", "ECLIPJ2000" );

    // Define link ends for observations.
    LinkEnds linkEnds;
    linkEnds[ transmitter ] = std::make_pair( "Earth" , "Station1"  );
    linkEnds[ receiver ] = std::make_pair( "MEX" , ""  );

    // Create observation settings
    boost::shared_ptr< ObservationSettings > rangeSettings = boost::make_shared< ObservationSettings >
            ( one_way_range );

    boost::shared_ptr< ObservationSettings > openLoopObservableSettings = boost::make_shared< ObservationSettings >
            ( one_way_doppler );

    boost::shared_ptr< ObservationSettings > closedLoopObservableSettings =
            boost::make_shared< OneWayDifferencedRangeRateObservationSettings >( boost::lambda::constant( 30.0 ) );

    // Create open-loop bservation model.
    boost::shared_ptr< ObservationModel< 1, double, double> > rangeObservationModel =
            ObservationModelCreator< 1, double, double>::createObservationModel(
                linkEnds, rangeSettings, bodyMap );

    boost::shared_ptr< ObservationModel< 1, double, double> > openLoopObservationModel =
            ObservationModelCreator< 1, double, double>::createObservationModel(
                linkEnds, openLoopObservableSettings, bodyMap );

    // Create closed-loop observation model.
    boost::shared_ptr< ObservationModel< 1, double, double> > closedLoopObservationModel =
            ObservationModelCreator< 1, double, double>::createObservationModel(
                linkEnds, closedLoopObservableSettings, bodyMap );


    boost::shared_ptr< OneWayDopplerObservationModel< double, double> > dopplerObservationModel =
            boost::dynamic_pointer_cast< OneWayDopplerObservationModel< double, double> >( openLoopObservableSettings );

    LinkEndType referenceLinkEnd = receiver;

//    double openLoopDopplerValue = openLoopObservationModel->computeObservationEntry( 0.0, referenceLinkEnd, 0 );
//    double closedLoopDopplerValue = closedLoopObservationModel->computeObservationEntry( 30.0 / 2.0, referenceLinkEnd, 0 );
//    double manualRangeDifference = ( rangeObservationModel->computeObservationEntry( 30.0 / 2.0, referenceLinkEnd, 0 ) -
//            rangeObservationModel->computeObservationEntry( -30.0 / 2.0, referenceLinkEnd, 0 ) ) / 30.0;

//    boost::function< double( const double ) > oneWayDopplerFunction =
//            boost::bind( &OneWayDopplerObservationModel< double, double >::computeObservationEntry,
//                         openLoopObservationModel, _1, referenceLinkEnd, 0 );
//    std::vector< double > currentObservableTimeLimits;
//    currentObservableTimeLimits.push_back( -15.0 );
//    currentObservableTimeLimits.push_back( 15.0 );

//    boost::shared_ptr< numerical_quadrature::TrapezoidNumericalQuadrature< double, double > > trapezoidalQuadrature =
//            boost::make_shared< numerical_quadrature::TrapezoidNumericalQuadrature< double, double > >( );
//    trapezoidalQuadrature->resetData(
//                currentObservableTimeLimits, oneWayDopplerFunction );
//    double integratedOpenLoopDopplerTrapezoidal = trapezoidalQuadrature->getQuadrature( ) / 30.0;

////    boost::shared_ptr< numerical_quadrature::SimpsonNumericalQuadrature< double, double > > simpsonQuadrature =
////            boost::make_shared< numerical_quadrature::SimpsonNumericalQuadrature< double, double > >( );
////    simpsonQuadrature->resetData(
////                currentObservableTimeLimits, oneWayDopplerFunction );
////    double integratedOpenLoopDopplerSimpson = simpsonQuadrature->getQuadrature( ) / 60.0;

//    std::cout<<"Difference: "<<integratedOpenLoopDopplerTrapezoidal<<" "<<//integratedOpenLoopDopplerSimpson<<" "<<
//               -closedLoopDopplerValue / physical_constants::SPEED_OF_LIGHT<<" "<<
//               integratedOpenLoopDopplerTrapezoidal + closedLoopDopplerValue / physical_constants::SPEED_OF_LIGHT<<" "<<std::endl;
              //integratedOpenLoopDopplerSimpson + closedLoopDopplerValue / physical_constants::SPEED_OF_LIGHT<<" "<<std::endl;

    double integrationTime = 60.0;

    std::map< double, double > syntheticClosedLoopData = convertOpenLoopToClosedLoopDopplerData(
                boost::dynamic_pointer_cast< OneWayDopplerObservationModel< double, double > >( openLoopObservationModel ),
                referenceLinkEnd, initialEphemerisTime, initialEphemerisTime + 600.0, integrationTime );


    for(  std::map< double, double >::const_iterator dataIterator = syntheticClosedLoopData.begin( );
          dataIterator != syntheticClosedLoopData.end( ); dataIterator++ )
    {
        double currentOpenLoopData = openLoopObservationModel->computeObservationEntry(
                    dataIterator->first, referenceLinkEnd, 0 );
        double currentClosedLoopData = closedLoopObservationModel->computeObservationEntry(
                    dataIterator->first + integrationTime / 2.0, referenceLinkEnd, 0 );
        double currentSyntheticClosedLoopData = dataIterator->second;
        double currentManualClosedLoopData =
                ( rangeObservationModel->computeObservationEntry( dataIterator->first + integrationTime / 2.0, referenceLinkEnd, 0 ) -
                  rangeObservationModel->computeObservationEntry( dataIterator->first - integrationTime / 2.0, referenceLinkEnd, 0 ) ) / integrationTime;

        double currentTime = dataIterator->first;

        std::cout<<"Difference in loop: "<<currentTime<<" "<<currentClosedLoopData<<" "<<currentSyntheticClosedLoopData<<" "<<currentManualClosedLoopData<<" "<<
                   currentSyntheticClosedLoopData + currentManualClosedLoopData / physical_constants::SPEED_OF_LIGHT<<" "<<currentOpenLoopData<<std::endl;
    }

//    // Create open-loop bservation model.
//    double testRange = rangeObservationModel->computeObservationEntry(
//                0.0, referenceLinkEnd, 0 );
//    double lightTime = testRange / physical_constants::SPEED_OF_LIGHT;

//    double testOpenLoopReceiverReference = openLoopObservationModel->computeObservationEntry(
//                0.0, referenceLinkEnd, 0 );

//    double testClosedLoopReceiverReference = closedLoopObservationModel->computeObservationEntry(
//                0.0, referenceLinkEnd, 0 );

//    double testOpenLoopTransmitterReference = openLoopObservationModel->computeObservationEntry(
//                -lightTime , transmitter, 0 );

//    double testClosedLoopTransmitterReference = closedLoopObservationModel->computeObservationEntry(
//                -lightTime, transmitter, 0 );

//    std::cout<<testOpenLoopReceiverReference - testOpenLoopTransmitterReference<<" "<<testOpenLoopTransmitterReference<<std::endl;

//    std::cout<<testClosedLoopReceiverReference - testClosedLoopTransmitterReference<<" "<<testClosedLoopTransmitterReference<<std::endl;

//    std::cout<<"Light time "<<lightTime<<std::endl;

}


BOOST_AUTO_TEST_SUITE_END( )

}

}


