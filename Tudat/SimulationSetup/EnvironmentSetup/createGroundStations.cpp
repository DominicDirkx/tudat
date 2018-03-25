/*    Copyright (c) 2010-2018, Delft University of Technology
 *    All rigths reserved
 *
 *    This file is part of the Tudat. Redistribution and use in source and
 *    binary forms, with or without modification, are permitted exclusively
 *    under the terms of the Modified BSD license. You should have received
 *    a copy of the license with this file. If not, please or visit:
 *    http://tudat.tudelft.nl/LICENSE.
 */

#include "Tudat/SimulationSetup/EnvironmentSetup/createGroundStations.h"
#include "Tudat/Astrodynamics/BasicAstrodynamics/timeConversions.h"
#include "Tudat/Astrodynamics/Ephemerides/rotationalEphemeris.h"
#include "Tudat/External/SpiceInterface/spiceInterface.h"

namespace tudat
{

namespace simulation_setup
{

//! Function to create a ground station from pre-defined station state object, and add it to a Body object
void createGroundStation(
        const boost::shared_ptr< Body >& body,
        const std::string groundStationName,
        const boost::shared_ptr< ground_stations::GroundStationState > groundStationState )
{
    boost::shared_ptr< ground_stations::PointingAnglesCalculator > pointingAnglesCalculator =
            boost::make_shared< ground_stations::PointingAnglesCalculator >(
                boost::bind( &ephemerides::RotationalEphemeris::getRotationToTargetFrame, body->getRotationalEphemeris( ), _1 ),
                boost::bind( &ground_stations::GroundStationState::getRotationFromBodyFixedToTopocentricFrame, groundStationState, _1 ) );
    body->addGroundStation( groundStationName, boost::make_shared< ground_stations::GroundStation >(
                                groundStationState, pointingAnglesCalculator, groundStationName ) );
}

//! Function to create a ground station and add it to a Body object
void createGroundStation(
        const boost::shared_ptr< Body >& body,
        const std::string groundStationName,
        const Eigen::Vector3d groundStationPosition,
        const coordinate_conversions::PositionElementTypes positionElementType )
{
    createGroundStation( body, groundStationName, boost::make_shared< ground_stations::GroundStationState >(
                             groundStationPosition, positionElementType, body->getShapeModel( ) ) );

}

void addSpiceDsnGroundStations(
        const boost::shared_ptr< Body >& body,
        const std::vector< int > groundStationIndices )
{
    double referenceTime = basic_astrodynamics::convertCalendarDateToJulianDaysSinceEpoch(
                2003, 1, 1, 0, 0, 0.0, basic_astrodynamics::JULIAN_DAY_ON_J2000 ) * physical_constants::JULIAN_DAY;

    for( unsigned int i = 0; i < groundStationIndices.size( ); i++ )
    {
        std::string currentStationName = "DSS-" + std::to_string( groundStationIndices.at( i ) );
        Eigen::Vector3d stationPosition =
                spice_interface::getBodyCartesianPositionAtEpoch(
                    currentStationName, "Earth", "ITRF93", "None", referenceTime );
        createGroundStation(
                body, currentStationName, stationPosition, coordinate_conversions::cartesian_position );
    }
}

void addAllSpiceDsnGroundStations(
        const boost::shared_ptr< Body >& body )
{
    addSpiceDsnGroundStations( body,
    { 12, 13, 14, 15, 16, 17, 23, 24, 25, 26, 27, 28, 33, 34, 42, 43, 45, 46, 49, 53, 54, 55, 61, 63, 64, 65, 66 } );
}

//! Function to create a set of ground stations and add them to the corresponding Body objects
void createGroundStations(
        const NamedBodyMap& bodyMap,
        const std::map< std::pair< std::string, std::string >, Eigen::Vector3d >& groundStationsWithPosition,
        const coordinate_conversions::PositionElementTypes positionElementType )
{
    for( std::map< std::pair< std::string, std::string >, Eigen::Vector3d >::const_iterator
         stationIterator = groundStationsWithPosition.begin( );
         stationIterator != groundStationsWithPosition.end( ); stationIterator++ )
    {
        if( bodyMap.count( stationIterator->first.first ) > 0 )
        {
            createGroundStation( bodyMap.at( stationIterator->first.first ), stationIterator->first.second,
                                 stationIterator->second, positionElementType );
        }
    }
}

void createGroundStation(
        const boost::shared_ptr< Body >& body,
        const std::string& bodyName,
        const boost::shared_ptr< GroundStationSettings > groundStationSettings )
{

    if( body->getGroundStationMap( ).count( groundStationSettings->getStationName( ) ) != 0 )
    {
        throw std::runtime_error(
                    "Error when creating ground station " + groundStationSettings->getStationName( ) +
                    " on body " + bodyName + ", station already exists." );
    }
    else
    {
        createGroundStation( body, groundStationSettings->getStationName( ),
                             groundStationSettings->getGroundStationPosition( ),
                             groundStationSettings->getPositionElementType( ) );
    }
}


}

}
