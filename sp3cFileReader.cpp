#include <fstream>
#include <sstream>
#include <string>
#include <iostream>
#include <iomanip>

#include <boost/algorithm/string.hpp>
#include <boost/algorithm/string/trim.hpp>
#include <boost/lexical_cast.hpp>
#include <boost/format.hpp>
#include <boost/make_shared.hpp>

#include "External/SofaInterface/sofaTimeConversions.h"
#include "InputOutput/sp3cFileReader.h"

namespace tudat
{

namespace input_output
{


boost::shared_ptr< SP3cFileContents > readSp3cFile(
        const std::string fileName, const double referenceJulianDay )
{
    boost::shared_ptr< SP3cFileContents > fileContents = boost::make_shared< SP3cFileContents >( );

    // Create stream from given file name.
    std::fstream stream( fileName.c_str( ), std::ios::in );

    // Check if stream was successfully created.
    if ( stream.fail( ) )
    {
        boost::throw_exception( std::runtime_error( boost::str(
                                                        boost::format( "Data file '%s' could not be opened." ) % fileName.c_str( ) ) ) );
    }

    std::vector< std::string > vectorOfIndividualStrings;
    std::string line;

    for( unsigned int i = 0; i < 22; i++ )
    {
        // Get line from file.
        std::getline( stream, line );

        // Trim input string (removes all leading and trailing whitespaces).
        boost::algorithm::trim( line );

        // Split string into multiple strings, each containing one element from a line from the
        // data file.
        boost::algorithm::split( vectorOfIndividualStrings,
                                 line,
                                 boost::algorithm::is_any_of( " " ),
                                 boost::algorithm::token_compress_on );
        switch( i )
        {
        case 0:
            if( vectorOfIndividualStrings.size( ) != 11 )
            {
                std::cerr<<"Error when reading line 1 of sp3c file"<<std::endl;
            }
            else
            {
                fileContents->frameName = vectorOfIndividualStrings.at( 8 );
            }
            break;
        case 2:
        {
            if( vectorOfIndividualStrings.size( ) != 19 )
            {
                std::cerr<<"Error when reading line 3 of sp3c file, line size is "<<vectorOfIndividualStrings.size( )<<" for line "<<line<<std::endl;
            }
            else
            {
                int numberOfSatellites = boost::lexical_cast< int >( vectorOfIndividualStrings.at( 1 ) );
                if( numberOfSatellites > 17 )
                {
                    std::cerr<<"Error when reading sp3c file, maximum of 17 satellites per file supported"<<std::endl;
                }
                else
                {
                    for( int i = 0; i < numberOfSatellites; i++ )
                    {
                        fileContents->vehiclesStates[ vectorOfIndividualStrings.at( i + 2 ) ] = std::map< double, Eigen::VectorXd >( );
                    }
                }
            }
            break;
        }
        case 12:
        {
            if( vectorOfIndividualStrings.size( ) != 13 )
            {
                std::cerr<<"Error when reading line 14 of sp3c file"<<std::endl;
            }
            else
            {
                fileContents->timeScale = vectorOfIndividualStrings.at( 3 );
            }
        }
        default:
            break;

        }
    }

    std::string currentSatelliteName;
    std::string currentStringToParse;

    double currentTime;
    std::map< std::string, Eigen::VectorXd > currentStates;
    int year, month, day;
    double hour, minute, second;

    double secondsOfDay;

    double mjd0, mjdAtEpoch;
    bool isEndOfFileReached = 0;

    while( !stream.fail( ) && !stream.eof( ) && !isEndOfFileReached )
    {
        std::getline( stream, line );

        boost::algorithm::trim( line );
        boost::algorithm::split( vectorOfIndividualStrings, line, boost::algorithm::is_any_of( " " ),
                                 boost::algorithm::token_compress_on );

        if( vectorOfIndividualStrings.at( 0 ) == "*" )
        {
            if( currentStates.size( ) != 0 )
            {
                for( std::map< std::string, Eigen::VectorXd >::iterator it = currentStates.begin( );
                     it != currentStates.end( ); it++)
                {
                    fileContents->vehiclesStates[ it->first ][ currentTime ] = it->second;
                }
                currentStates.clear( );
            }

            year = boost::lexical_cast< int >( vectorOfIndividualStrings.at( 1 ) );
            month = boost::lexical_cast< int >( vectorOfIndividualStrings.at( 2 ) );
            day = boost::lexical_cast< int >( vectorOfIndividualStrings.at( 3 ) );
            hour = boost::lexical_cast< double >( vectorOfIndividualStrings.at( 4 ) );
            minute = boost::lexical_cast< double >( vectorOfIndividualStrings.at( 5 ) );
            second = boost::lexical_cast< double >( vectorOfIndividualStrings.at( 6 ) );
            if( iauCal2jd( year, month, day, &mjd0, &mjdAtEpoch ) != 0 )
            {
                std::cerr<<"Error when parsing time stamp of 3p3c file at "<<year<<", "<<month<<", "<<day<<std::endl;
            }
            else
            {
                if( mjd0 != basic_astrodynamics::JULIAN_DAY_AT_0_MJD )
                {
                    std::cerr<<"Error, inconsistent return data at level "<<mjd0 - basic_astrodynamics::JULIAN_DAY_AT_0_MJD<<" when converting to mjd in sp3c reader."<<std::endl;
                }
                secondsOfDay = 3600.0 * hour + 60.0 * minute + second;
                currentTime = ( mjdAtEpoch - ( referenceJulianDay - mjd0 ) ) * physical_constants::JULIAN_DAY + secondsOfDay;
                //std::cout<<( mjd0 - referenceJulianDay )<<" "<<mjd0<<" "<<mjdAtEpoch<<std::endl;
            }
        }
        else if( vectorOfIndividualStrings.at( 0 ).at( 0 ) == 'P' )
        {
            currentStringToParse = vectorOfIndividualStrings.at( 0 );
            currentSatelliteName = currentStringToParse.substr( 1, currentStringToParse.size( ) - 1 );

            if( fileContents->vehiclesStates.count( currentSatelliteName ) == 0 )
            {
                std::cerr<<"Error, satellite "<<currentSatelliteName<<" not found in list of satellites "<<currentStringToParse<<std::endl;
            }
            else
            {
                if( currentStates.count( currentSatelliteName ) == 0 )
                {
                    currentStates[ currentSatelliteName ] = Eigen::VectorXd::Zero( 6 );
                }
                currentStates[ currentSatelliteName ]( 0 ) =
                        boost::lexical_cast< double >( vectorOfIndividualStrings.at( 1 ) ) * 1000.0;
                currentStates[ currentSatelliteName ]( 1 ) =
                        boost::lexical_cast< double >( vectorOfIndividualStrings.at( 2 ) ) * 1000.0;
                currentStates[ currentSatelliteName ]( 2 ) =
                        boost::lexical_cast< double >( vectorOfIndividualStrings.at( 3 ) ) * 1000.0;
            }

        }
        else if( vectorOfIndividualStrings.at( 0 ).at( 0 ) == 'V' )
        {
            currentStringToParse = vectorOfIndividualStrings.at( 0 );
            currentSatelliteName = currentStringToParse.substr( 1, currentStringToParse.size( ) - 1 );

            if( fileContents->vehiclesStates.count( currentSatelliteName ) == 0 )
            {
                std::cerr<<"Error, satellite "<<currentSatelliteName<<" not found in list of satellites"<<std::endl;
            }
            else
            {
                if( currentStates.count( currentSatelliteName ) == 0 )
                {
                    currentStates[ currentSatelliteName ] = Eigen::VectorXd::Zero( 6 );
                }
                currentStates[ currentSatelliteName ]( 3 ) =
                        boost::lexical_cast< double >( vectorOfIndividualStrings.at( 1 ) ) / 10.0;
                currentStates[ currentSatelliteName ]( 4 ) =
                        boost::lexical_cast< double >( vectorOfIndividualStrings.at( 2 ) ) / 10.0;
                currentStates[ currentSatelliteName ]( 5 ) =
                        boost::lexical_cast< double >( vectorOfIndividualStrings.at( 3 ) ) / 10.0;

            }
        }
        else if( vectorOfIndividualStrings.at( 0 ) == "EOF" )
        {
            isEndOfFileReached = 1;
        }
        else
        {
            std::cerr<<"Error, line start "<<vectorOfIndividualStrings.at( 0 )<<" not recognized when reading sp3c file."<<std::endl;
        }
    }

    if( currentStates.size( ) != 0 )
    {
        for( std::map< std::string, Eigen::VectorXd >::iterator it = currentStates.begin( );
             it != currentStates.end( ); it++)
        {
            fileContents->vehiclesStates[ it->first ][ currentTime ] = it->second;
        }
        currentStates.clear( );
    }

    return fileContents;
}

}

}
