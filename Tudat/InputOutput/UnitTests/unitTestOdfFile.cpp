#include <iostream>
#include <fstream>
#include <bitset>
#include <cmath>
#include <vector>
#include <map>

#include <boost/make_shared.hpp>

#include "Tudat/Astrodynamics/BasicAstrodynamics/timeConversions.h"
#include "Tudat/InputOutput/readOdfFile.h"
#include "Tudat/Mathematics/Interpolators/lookupScheme.h"
#include "Tudat/Astrodynamics/Propagators/parseOdfFile.h"

int main( )
{
//   boost::shared_ptr< tudat::orbit_determination::ProcessedOdfFileContents > odfContents =
//           tudat::orbit_determination::parseOdfFileContents(
//               tudat::input_output::readOdfFile( "/home/dominic/Downloads/mromagr2017_117_0745xmmmv1.odf" ) );

//   std::map< tudat::observation_models::ObservableType,
//           std::vector< boost::shared_ptr< tudat::orbit_determination::ProcessdOdfFileSingleLinkData > > > dataBlocks =
//           odfContents->dataBlocks;

//   for( auto it = dataBlocks.begin( ); it != dataBlocks.end( ); it++ )
//   {
//       int counter = 0;
//       for( unsigned int i = 0; i < it->second.size( ); i++ )
//       {
//           boost::shared_ptr< tudat::orbit_determination::ProcessdOdfFileDopplerData > currentDopplerData =
//                   boost::dynamic_pointer_cast< tudat::orbit_determination::ProcessdOdfFileDopplerData >(
//                       it->second.at( i ) );
//           std::string fileSuffix = std::to_string( it->first ) + "_" + std::to_string( counter );

//           tudat::input_output::writeDataMapToTextFile(
//                       currentDopplerData->getCompressionTimes( ),
//                           "odfTestCompressionTimes_" + fileSuffix + ".dat", "/home/dominic/Documents/" ) ;
//           tudat::input_output::writeDataMapToTextFile(
//                       currentDopplerData->getObservationData( ),
//                           "odfTestObservations_" + fileSuffix + ".dat", "/home/dominic/Documents/" ) ;
//           tudat::input_output::writeDataMapToTextFile(
//                       currentDopplerData->getReferenceFrequencies( ),
//                           "odfTestReferenceFrequencies_" + fileSuffix + ".dat", "/home/dominic/Documents/" ) ;
//           tudat::input_output::writeDataMapToTextFile(
//                       currentDopplerData->getRampFlags( ),
//                           "odfTestRampFlags_" + fileSuffix + ".dat", "/home/dominic/Documents/" ) ;
//           counter++;
//       }
//   }

   std::vector< boost::filesystem::path > files = tudat::input_output::listAllFilesInDirectory(
       "/home/dominic/Software/MercuryData/odf/", false );

   std::vector< boost::shared_ptr< tudat::orbit_determination::ProcessedOdfFileContents > > odfContentsList;

   //for( unsigned int i = 0; i < files.size( ); i++ )
   for( unsigned int i = 0; i < 1000; i++ )
   {
       std::string fileString = files.at( i ).string( );
       int stringSize = fileString.size( );

       if( fileString.substr( stringSize - 3, stringSize -1 ) == "dat" )
       {
           std::cout<<fileString<<std::endl;
           odfContentsList.push_back( tudat::orbit_determination::parseOdfFileContents(
               tudat::input_output::readOdfFile( "/home/dominic/Software/MercuryData/odf/" + fileString ) ) );
           std::cout<<std::endl;
       }
   }


   boost::shared_ptr< tudat::orbit_determination::ProcessedOdfFileContents > mergedData =
           tudat::orbit_determination::mergeOdfFileContents(
               odfContentsList );

   std::cout<<files.size( )<<std::endl;

   std::map< tudat::observation_models::ObservableType,
           std::vector< boost::shared_ptr< tudat::orbit_determination::ProcessdOdfFileSingleLinkData > > > dataBlocksMerged =
           mergedData->dataBlocks;


   for( auto it = dataBlocksMerged.begin( ); it != dataBlocksMerged.end( ); it++ )
   {
       int counter = 0;

       std::vector< double > observationTimes, obsevables, referenceFrequencies;
       std::vector< std::string > receivingStation, transmittingStation;
       std::vector< std::string > originFiles;
       std::vector< bool > rampFlags;

       for( unsigned int i = 0; i < it->second.size( ); i++ )
       {
           boost::shared_ptr< tudat::orbit_determination::ProcessdOdfFileDopplerData > currentDopplerData =
                   boost::dynamic_pointer_cast< tudat::orbit_determination::ProcessdOdfFileDopplerData >(
                       it->second.at( i ) );

           std::vector< double > currentObservationTimes = currentDopplerData->observationTimes,
                   currentObsevables = currentDopplerData->observableValues,
                   currentReferenceFrequencies = currentDopplerData->referenceFrequency;
           std::vector< std::string > currentOriginFiles = currentDopplerData->originFile;
           std::vector< bool > currentRampFlags = currentDopplerData->rampingFlag;

           observationTimes.insert( observationTimes.end( ), currentObservationTimes.begin( ), currentObservationTimes.end( ) );
           obsevables.insert( obsevables.end( ), currentObsevables.begin( ), currentObsevables.end( ) );
           originFiles.insert( originFiles.end( ), currentOriginFiles.begin( ), currentOriginFiles.end( ) );
           rampFlags.insert( rampFlags.end( ), currentRampFlags.begin( ), currentRampFlags.end( ) );
           referenceFrequencies.insert(
                       referenceFrequencies.end( ), currentReferenceFrequencies.begin( ), currentReferenceFrequencies.end( ) );

           for( unsigned int i = 0; i < currentObservationTimes.size( ); i++ )
           {
                receivingStation.push_back( currentDopplerData->receivingStation );
                transmittingStation.push_back( currentDopplerData->transmittingStation );
           }

//           std::string fileSuffix = std::to_string( it->first ) + "_" + std::to_string( counter );

//           tudat::input_output::writeDataMapToTextFile(
//                       currentDopplerData->getCompressionTimes( ),
//                           "odfTestCompressionTimes_" + fileSuffix + ".dat", "/home/dominic/Documents/" ) ;
//           tudat::input_output::writeDataMapToTextFile(
//                       currentDopplerData->getObservationData( ),
//                           "odfTestObservations_" + fileSuffix + ".dat", "/home/dominic/Documents/" ) ;
//           tudat::input_output::writeDataMapToTextFile(
//                       currentDopplerData->getReferenceFrequencies( ),
//                           "odfTestReferenceFrequencies_" + fileSuffix + ".dat", "/home/dominic/Documents/" ) ;
//           tudat::input_output::writeDataMapToTextFile(
//                       currentDopplerData->getRampFlags( ),
//                           "odfTestRampFlags_" + fileSuffix + ".dat", "/home/dominic/Documents/" ) ;
           counter++;
       }

       Eigen::MatrixXd dataMatrix = Eigen::MatrixXd( rampFlags.size( ), 6 );
       for( int j = 0; j < rampFlags.size( ); j++ )
       {
            dataMatrix( j, 0 ) = observationTimes.at( j );
            dataMatrix( j, 1 ) = obsevables.at( j );
            dataMatrix( j, 2 ) = boost::lexical_cast< double >( receivingStation.at( j ) );
            dataMatrix( j, 3 ) = boost::lexical_cast< double >( transmittingStation.at( j ) );
            dataMatrix( j, 4 ) = referenceFrequencies.at( j );
            dataMatrix( j, 5 ) = static_cast< double >( rampFlags.at( j ) );

       }

       tudat::input_output::writeMatrixToFile(
                   dataMatrix, "odfFileSummary_" + std::to_string( it->first ), 16, "/home/dominic/Documents/" );
   }


}

//uint32_t read_u32_le(std::istream& file)
//{
//    uint32_t value;
//    uint8_t  bytes[4];

//    file.read( (char*)bytes, 4);
//    value = bytes[0] << 24 | (bytes[1] << 16) | (bytes[2] << 8) | (bytes[3]);

//    //std::cout<<value<<" "<<int( bytes[0] )<<" "<<int( bytes[1] << 8 )<<" "<<int( bytes[2] << 16 )<<" "<<( bytes[3] << 0 )<<std::endl;
//    return value;
//}

//int32_t read_s32_le(std::istream& file)
//{
//    int32_t value;
//    int8_t  bytes[4];

//    file.read( (char*)bytes, 4);
//    value = bytes[0] << 24| (bytes[1] << 16) | (bytes[2] << 8) | (bytes[3]);
//    //std::cout<<value<<" "<<int( bytes[0] )<<" "<<int( bytes[1] << 8 )<<" "<<int( bytes[2] << 16 )<<" "<<( bytes[3] << 0 )<<std::endl;

//    return value;
//}

//template< int OutputBits, int InputBits >
//std::bitset< OutputBits > getBitsetSegment(
//        const std::bitset< InputBits > inputBits,
//        const int startIndex )
//{
//    std::bitset< OutputBits > outputBits;
//    for( unsigned int i = 0; i < OutputBits; i++ )
//    {
//        outputBits[ i ] = inputBits[ InputBits - OutputBits - startIndex + i  ];
//    }
//    return outputBits;
//}

//template< int NumberOfBits >
//unsigned int getUnsignedNBitInteger(
//        std::bitset< NumberOfBits > inputBits )
//{
//    int outputInteger = 0;
//    for( unsigned int i = 0; i < NumberOfBits ; i ++ )
//    {
//        outputInteger += inputBits[ i ] * std::pow( 2.0, i );
//    }
//    return outputInteger;
//}


//template< int NumberOfBits >
//int getSignedNBitInteger(
//        std::bitset< NumberOfBits > inputBits )
//{
//    int outputInteger = -inputBits[ NumberOfBits - 1 ] * std::pow( 2, NumberOfBits - 1 );

//    for( unsigned int i = 0; i < NumberOfBits - 1; i ++ )
//    {
//        outputInteger += inputBits[ i ] * std::pow( 2.0, i );
//    }
//    return outputInteger;
//}


//uint32_t convertCharactersToUnsignedInt32(
//        char data[ 4 ] )
//{
//    uint8_t bytes[4];
//    bytes[ 0 ] = data[ 0 ];
//    bytes[ 1 ] = data[ 1 ];
//    bytes[ 2 ] = data[ 2 ];
//    bytes[ 3 ] = data[ 3 ];

//    return bytes[0] << 24| (bytes[1] << 16) | (bytes[2] << 8) | (bytes[3]);
//}

//int32_t convertCharactersToSignedInt32(
//        char data[ 4 ] )
//{
//    uint8_t bytes[4];
//    bytes[ 0 ] = data[ 0 ];
//    bytes[ 1 ] = data[ 1 ];
//    bytes[ 2 ] = data[ 2 ];
//    bytes[ 3 ] = data[ 3 ];

//    int32_t returnData = bytes[0] << 24| (bytes[1] << 16) | (bytes[2] << 8) | (bytes[3]);
//    return returnData;
//}


//void readOdfFileBlock(
//        char fileBlock[ 9 ][ 4 ],
//std::istream& file )
//{
//    for( unsigned int i = 0; i < 9; i++ )
//    {
//        file.read( (char*)fileBlock[ i ], 4 );
//    }
//}

//void parseHeader( char fileBlock[ 9 ][ 4 ],
//int32_t& primaryKey,
//uint32_t& secondaryKey,
//uint32_t& logicalrecordLength,
//uint32_t& groupStartPacketNumber )
//{
//    primaryKey = convertCharactersToSignedInt32( fileBlock[ 0 ] );
//    secondaryKey = convertCharactersToUnsignedInt32( fileBlock[ 1 ] );
//    logicalrecordLength = convertCharactersToUnsignedInt32( fileBlock[ 2 ] );
//    groupStartPacketNumber = convertCharactersToUnsignedInt32( fileBlock[ 3 ] );
//    uint32_t testInteger;
//    for( int i = 0; i < 5; i++ )
//    {
//        testInteger = convertCharactersToUnsignedInt32( fileBlock[ i + 4 ] );
//        if( testInteger != 0 )
//        {
//            throw std::runtime_error( "Error when reading ODF file, header file inconsistent" );
//        }
//    }
//}

//int currentBlockIsHeader(  char fileBlock[ 9 ][ 4 ], int& secondaryKeyInt )
//{
//    int headerType = -2;

//    int32_t primaryKey;
//    uint32_t secondaryKey;
//    uint32_t logicalrecordLength;
//    uint32_t groupStartPacketNumber;

//    try
//    {
//        parseHeader( fileBlock, primaryKey, secondaryKey, logicalrecordLength, groupStartPacketNumber );
//        secondaryKeyInt = secondaryKey;
//        if( primaryKey ==  101 )
//        {
//            headerType = 1;
//        }
//        else if( primaryKey == 107 )
//        {
//            headerType = 2;
//        }
//        else if( primaryKey == 109 )
//        {
//            headerType = 3;
//        }
//        else if( primaryKey == 2030 )
//        {
//            headerType = 4;
//        }
//        else if( primaryKey == 2040 )
//        {
//            headerType = 5;
//        }
//        else if( primaryKey == -1 )
//        {
//            headerType = -1;
//        }
//    }
//    catch( std::runtime_error )
//    {
//        headerType = 0;
//    }
//    return headerType;
//}

//void parseFileLabel( char fileBlock[ 9 ][ 4 ],
//std::string& systemId, std::string& programId, std::string& fileCreationDate, std::string& fileCreationTme,
//uint32_t& spacecraftIdNumber, uint32_t& fileReferenceDate, uint32_t& fileReferenceTime )
//{
//    systemId.resize( 8 );
//    programId.resize( 8 );

//    fileCreationDate.resize( 4 );
//    fileCreationTme.resize( 4 );

//    for( unsigned int i = 0; i < 4; i++ )
//    {
//        systemId[ i ] = fileBlock[ 0 ][ i ];
//        systemId[ i + 4 ] = fileBlock[ 1 ][ i ];

//        programId[ i ] = fileBlock[ 2 ][ i ];
//        programId[ i + 4 ] = fileBlock[ 3 ][ i ];

//        fileCreationDate[ i ] = fileBlock[ 5 ][ i ];
//        fileCreationTme[ i ] = fileBlock[ 6 ][ i ];
//    }
//    spacecraftIdNumber = convertCharactersToUnsignedInt32( fileBlock[ 4 ] );
//    fileReferenceDate = convertCharactersToUnsignedInt32( fileBlock[ 7 ] );
//    fileReferenceTime = convertCharactersToUnsignedInt32( fileBlock[ 8 ] );
//}

//void parseIdentifierGroup( char fileBlock[ 9 ][ 4 ] )
//{
//    std::string testStringA, testStringB, testStringC;
//    testStringA.resize( 8 );
//    testStringB.resize( 8 );
//    testStringC.resize( 20 );

//    for( unsigned int i = 0; i < 4; i++ )
//    {
//        testStringA[ i ] = fileBlock[ 0 ][ i ];
//        testStringA[ i + 4 ] = fileBlock[ 1 ][ i ];

//        testStringB[ i ] = fileBlock[ 2 ][ i ];
//        testStringB[ i + 4 ] = fileBlock[ 3 ][ i ];

//        for( unsigned int j = 0; j < 5; j++ )
//        {
//            testStringC[ i + j * 4 ] = fileBlock[ 4 + j ][ i ];
//        }
//    }

//    std::cout<<testStringA<<" "<<testStringB<<" "<<testStringC<<std::endl;
//}

//template< int FirstInputSize, int SecondInputSize >
//std::bitset< FirstInputSize + SecondInputSize > mergeBitsets(
//        std::bitset< FirstInputSize > firstInput,
//        std::bitset< SecondInputSize > secondInput )
//{
//    std::bitset< FirstInputSize + SecondInputSize > returnBitset;
//    for( int i = 0; i < FirstInputSize; i++ )
//    {
//        returnBitset[ i + SecondInputSize ] = firstInput[ i ];
//    }

//    for( int i = 0; i < SecondInputSize; i++ )
//    {
//        returnBitset[ i ] = secondInput[ i ];
//    }
//    return returnBitset;
//}

//struct OdfDopplerDataBlock
//{

//    int receiverChannel;
//    int spacecraftId;
//    int receiverExciterFlag;
//    int referenceFrequencyHighPart;
//    int referenceFrequencyLowPart;

//    int reservedSegment;
//    int compressionTime;

//    int transmittingStationDelay;

//    double getReferenceFrequency( )
//    {
//        std::cout<<std::setprecision( 16 )<<referenceFrequencyHighPart<<" "<<referenceFrequencyLowPart<<" "<<
//                   std::pow( 2.0, 24 ) /1.0E3 * referenceFrequencyHighPart<<" "<<referenceFrequencyLowPart / 1.0E3<<std::endl;
//        return std::pow( 2.0, 24 )  /1.0E3 * referenceFrequencyHighPart + referenceFrequencyLowPart / 1.0E3;
//    }

//    void printContents( )
//    {
//        std::cout<<"Receiver: "<<receiverChannel<<" "<<receiverExciterFlag<<std::endl;
//        std::cout<<"Spacecraft: "<<spacecraftId<<std::endl;
//        std::cout<<"Reference frequency: "<<referenceFrequencyHighPart<<" "<<referenceFrequencyLowPart<<std::endl;

//        std::cout<<"Reserved: "<<reservedSegment<<std::endl;
//        std::cout<<"Compression time: "<<compressionTime<<std::endl;
//        std::cout<<"Transmission delay: "<<transmittingStationDelay<<std::endl<<std::endl;;

//    }

//};

//struct OdfCommonDataBlock
//{
//    int integerTimeTag;
//    int fractionalTimeTag;
//    int receivingStationDownlinkDelay;

//    int integerObservable;
//    int fractionalObservable;

//    int formatId;
//    int receivingStation;
//    int transmittingStation;
//    int transmittingStationNetworkId;
//    int dataType;

//    int downlinkBand;
//    int uplinkBand;
//    int referenceBand;
//    int validity;

//    void printContents( )
//    {
//        std::cout<<"Time: "<<integerTimeTag<<" "<<fractionalTimeTag<<std::endl;
//        std::cout<<"Downlink delay: "<<receivingStationDownlinkDelay<<std::endl;
//        std::cout<<"Observable: "<<integerObservable<<" "<<fractionalObservable<<std::endl;

//        std::cout<<"Format id: "<<formatId<<std::endl;
//        std::cout<<"Station data: "<<receivingStation<<" "<<transmittingStation<<" "<<transmittingStationNetworkId<<" "<<std::endl;
//        std::cout<<"Data type: "<<dataType<<std::endl;
//        std::cout<<"Bands: "<<downlinkBand<<" "<<uplinkBand<<" "<<referenceBand<<" "<<std::endl;
//        std::cout<<"Validity: "<<validity<<std::endl<<std::endl;

//    }
//};

//struct OdfRampBlock
//{
//    int integerRampStartTime;
//    int fractionalRampStartTime;

//    int integerRampEndTime;
//    int fractionalRampEndTime;

//    int integerRampRate;
//    int fractionalRampRate;

//    int integerRampStartFrequency;
//    int integerRampStartFrequencyModulo;
//    int fractionalRampStartFrequency;

//    int transmittingStationId;

//    double getRampStartFrequency( )
//    {
//        return static_cast< double >( integerRampStartFrequency ) * 1.0E9 +
//                static_cast< double >( integerRampStartFrequencyModulo ) +
//                static_cast< double >( fractionalRampStartFrequency ) * 1.0E-9;

//    }

//    double getRampRate( )
//    {
//        return static_cast< double >( integerRampRate ) +
//                static_cast< double >( fractionalRampRate ) * 1.0E-9;

//    }

//    double getRampStartTime( )
//    {
//        return static_cast< double >( integerRampStartTime ) -
//                //                ( tudat::basic_astrodynamics::JULIAN_DAY_ON_J2000 - tudat::basic_astrodynamics::convertCalendarDateToJulianDay(
//                //                        1950, 1, 1, 0, 0, 0.0 ) ) *
//                86400.0 + static_cast< double >( fractionalRampStartTime ) * 1.0E-9;
//    }

//    double getRampEndTime( )
//    {
//        return static_cast< double >( integerRampEndTime ) -
//                //                ( tudat::basic_astrodynamics::JULIAN_DAY_ON_J2000 - tudat::basic_astrodynamics::convertCalendarDateToJulianDay(
//                //                        1950, 1, 1, 0, 0, 0.0 ) ) * 86400.0 +
//                static_cast< double >( fractionalRampEndTime ) * 1.0E-9;
//    }

//    void printContents( )
//    {
//        std::cout<<"Start time "<<integerRampStartTime<<" "<<fractionalRampStartTime<<std::endl;
//        std::cout<<"End time "<<integerRampEndTime<<" "<<fractionalRampEndTime<<std::endl;
//        std::cout<<"Ramp rate "<<integerRampRate<<" "<<fractionalRampRate<<std::endl;
//        std::cout<<"Start frequency "<<integerRampStartFrequency<<" "<<integerRampStartFrequencyModulo<<" "<<
//                   fractionalRampStartFrequency<<std::endl;
//        std::cout<<"Station "<<transmittingStationId<<std::endl<<std::endl;

//    }
//};

//class RampedReferenceFrequencyInterpolator
//{
//public:
//    RampedReferenceFrequencyInterpolator(
//            std::vector< OdfRampBlock > rampBlock )
//    {
//        for( unsigned int i = 0; i < rampBlock.size( ); i++ )
//        {
//            startTimes.push_back( rampBlock.at( i ).getRampStartTime( ) );
//            endTimes.push_back( rampBlock.at( i ).getRampEndTime( ) );
//            rampRates.push_back( rampBlock.at( i ).getRampRate( ) );
//            startFrequency.push_back( rampBlock.at( i ).getRampStartFrequency( ) );
//        }

//        startTimeLookupScheme_ = boost::make_shared<
//                tudat::interpolators::HuntingAlgorithmLookupScheme< double > >(
//                    startTimes );
//    }

//    double getCurrentReferenceFrequencyIntegral(
//            const double lookupTime, const double integrationTime )
//    {
//        //        std::cout<<startTimes.size( )<<std::endl;
//        //        std::cout<<lookupTime<<" "<<lookupTime - startTimes.at( 0 )<<" "<<lookupTime - endTimes.at( endTimes.size( ) - 1 )<<std::endl;

//        int lowerNearestNeighbour = startTimeLookupScheme_->findNearestLowerNeighbour(
//                    lookupTime );
//        if( lookupTime > endTimes.at( lowerNearestNeighbour ) )
//        {
//            //std::cout<<lookupTime - endTimes.at( lowerNearestNeighbour )<<" "<<lowerNearestNeighbour<<" "<<startTimes.size( )<<std::endl;
//            //std::cout<<lookupTime - startTimes.at( lowerNearestNeighbour + 1 )<<" "<<startTimes.size( )<<std::endl<<std::endl;

//            //throw std::runtime_error( "Error when getting ramp frequency, time not inside arc range" );
//        }

//        return integrationTime * (
//                    startFrequency.at( lowerNearestNeighbour ) +
//                    rampRates.at( lowerNearestNeighbour ) * ( lookupTime - startTimes.at( lowerNearestNeighbour ) ) +
//                    0.5 * rampRates.at( lowerNearestNeighbour ) * integrationTime );
//    }

//    double getCurrentReferenceFrequency(
//            const double lookupTime, const double integrationTime )
//    {
//        int lowerNearestNeighbour = startTimeLookupScheme_->findNearestLowerNeighbour(
//                    lookupTime );
//        if( lookupTime > endTimes.at( lowerNearestNeighbour ) )
//        {
//            //std::cout<<lookupTime - endTimes.at( lowerNearestNeighbour )<<" "<<lowerNearestNeighbour<<" "<<startTimes.size( )<<std::endl;
//            //std::cout<<lookupTime - startTimes.at( lowerNearestNeighbour + 1 )<<" "<<startTimes.size( )<<std::endl<<std::endl;

//            //throw std::runtime_error( "Error when getting ramp frequency, time not inside arc range" );
//        }

//        return  startFrequency.at( lowerNearestNeighbour ) +
//                    rampRates.at( lowerNearestNeighbour ) * ( lookupTime - startTimes.at( lowerNearestNeighbour ) );
//    }



//private:

//    std::vector< double > startTimes;
//    std::vector< double > endTimes;
//    std::vector< double > rampRates;
//    std::vector< double > startFrequency;

//    boost::shared_ptr< tudat::interpolators::LookUpScheme< double > > startTimeLookupScheme_;

//};

//struct OdfDataBlock
//{
//    OdfDopplerDataBlock dopplerDataBlock;
//    OdfCommonDataBlock commonDataBlock;
//};

//OdfDopplerDataBlock parseDopplerOrbitData( char fileBlock[ 9 ][ 4 ] )
//{
//    OdfDopplerDataBlock dopplerDataBlock;

//    std::bitset< 32 > dataBits = std::bitset< 32 >( convertCharactersToUnsignedInt32( fileBlock[ 5 ] ) );

//    dopplerDataBlock.receiverChannel = getBitsetSegment< 7, 32 >( dataBits, 0 ).to_ulong( );
//    dopplerDataBlock.spacecraftId = getBitsetSegment< 10, 32 >( dataBits, 7 ).to_ulong( );
//    dopplerDataBlock.receiverExciterFlag = getBitsetSegment< 1, 32 >( dataBits, 17 ).to_ulong( );

//    std::bitset< 14 > referenceFrequencyHighPartBitsA = getBitsetSegment< 14, 32 >( dataBits, 18 );
//    dataBits = std::bitset< 32 >( convertCharactersToUnsignedInt32( fileBlock[ 6 ] ) );
//    std::bitset< 8 > referenceFrequencyHighPartBitsB = getBitsetSegment< 8, 32 >( dataBits, 0 );
//    dopplerDataBlock.referenceFrequencyHighPart = mergeBitsets< 14, 8 >(
//                referenceFrequencyHighPartBitsA, referenceFrequencyHighPartBitsB ).to_ulong();
//    dopplerDataBlock.referenceFrequencyLowPart = getBitsetSegment< 24, 32 >( dataBits, 8 ).to_ulong( );

//    //std::cout<<dopplerDataBlock.referenceFrequencyHighPart<<" "<<dopplerDataBlock.referenceFrequencyLowPart<<std::endl;

//    dataBits = std::bitset< 32 >( convertCharactersToUnsignedInt32( fileBlock[ 7 ] ) );

//    dopplerDataBlock.reservedSegment = getBitsetSegment< 20, 32 >( dataBits, 0 ).to_ulong( );
//    std::bitset< 12 > compressionTimeABits = getBitsetSegment< 12, 32 >( dataBits, 20 );
//    dataBits = std::bitset< 32 >( convertCharactersToUnsignedInt32( fileBlock[ 8 ] ) );
//    std::bitset< 10 > compressionTimeBBits = getBitsetSegment< 10, 32 >( dataBits, 0 );
//    dopplerDataBlock.compressionTime = mergeBitsets< 12, 10 >(
//                compressionTimeABits, compressionTimeBBits ).to_ulong( );
//    dopplerDataBlock.transmittingStationDelay = getBitsetSegment< 22, 32 >( dataBits, 10 ).to_ulong( );

//    //dopplerDataBlock.printContents( );

//    return dopplerDataBlock;
//}

//OdfDataBlock parseOrbitData( char fileBlock[ 9 ][ 4 ] )
//{
//    OdfCommonDataBlock commonDataBlock;

//    commonDataBlock.integerTimeTag = convertCharactersToUnsignedInt32( fileBlock[ 0 ] );

//    std::bitset< 32 > dataBits = std::bitset< 32 >( convertCharactersToUnsignedInt32( fileBlock[ 1 ] ) );
//    commonDataBlock.fractionalTimeTag = getBitsetSegment< 10, 32 >(
//                dataBits, 0 ).to_ulong( );
//    commonDataBlock.receivingStationDownlinkDelay = getBitsetSegment< 22, 32 >(
//                dataBits, 10 ).to_ulong( );

//    commonDataBlock.integerObservable = convertCharactersToSignedInt32( fileBlock[ 2 ] );
//    commonDataBlock.fractionalObservable = convertCharactersToSignedInt32( fileBlock[ 3 ] );


//    dataBits = std::bitset< 32 >( convertCharactersToUnsignedInt32( fileBlock[ 4 ] ) );


//    commonDataBlock.formatId = getBitsetSegment< 3, 32 >( dataBits, 0 ).to_ulong( );
//    commonDataBlock.receivingStation = getBitsetSegment< 7, 32 >( dataBits, 3 ).to_ulong( );
//    commonDataBlock.transmittingStation = getBitsetSegment< 7, 32 >( dataBits, 10 ).to_ulong( );

//    //    std::cout<<std::setprecision( 16 )<<static_cast< double >( commonDataBlock.integerTimeTag ) +
//    //               static_cast< double >(commonDataBlock.fractionalTimeTag ) / 1000.0<<" "<<
//    //               commonDataBlock.receivingStation<<" "<<commonDataBlock.transmittingStation<<std::endl;

//    commonDataBlock.transmittingStationNetworkId = getBitsetSegment< 2, 32 >( dataBits, 17 ).to_ulong( );
//    commonDataBlock.dataType =  getBitsetSegment< 6, 32 >( dataBits, 19 ).to_ulong( );

//    commonDataBlock.downlinkBand = getBitsetSegment< 2, 32 >( dataBits, 25 ).to_ulong( );
//    commonDataBlock.uplinkBand = getBitsetSegment< 2, 32 >( dataBits, 27 ).to_ulong( );
//    commonDataBlock.referenceBand = getBitsetSegment< 2, 32 >( dataBits, 29 ).to_ulong( );
//    commonDataBlock.validity = getBitsetSegment< 1, 32 >( dataBits, 31 ).to_ulong( );

//    //commonDataBlock.printContents( );

//    OdfDataBlock dataBlock;
//    dataBlock.commonDataBlock = commonDataBlock;
//    if( commonDataBlock.dataType == 12 || commonDataBlock.dataType == 13 )
//    {
//        dataBlock.dopplerDataBlock = parseDopplerOrbitData( fileBlock );
//    }
//    else
//    {
//        //std::cout<<commonDataBlock.dataType<<std::endl;
//        //throw std::runtime_error( "Error, data type not recognized" );
//    }
//    return dataBlock;
//}

//OdfRampBlock parseRampData( char fileBlock[ 9 ][ 4 ] )
//{

//    OdfRampBlock rampBlock;
//    rampBlock.integerRampStartTime = convertCharactersToUnsignedInt32( fileBlock[ 0 ] );
//    rampBlock.fractionalRampStartTime = convertCharactersToUnsignedInt32( fileBlock[ 1 ] );

//    rampBlock.integerRampRate = convertCharactersToSignedInt32( fileBlock[ 2 ] );
//    rampBlock.fractionalRampRate = convertCharactersToSignedInt32( fileBlock[ 3 ] );

//    std::bitset< 32 > dataBits = std::bitset< 32 >( convertCharactersToUnsignedInt32( fileBlock[ 4 ] ) );
//    rampBlock.transmittingStationId = getBitsetSegment< 10, 32 >( dataBits, 22 ).to_ulong( );
//    rampBlock.integerRampStartFrequency = getBitsetSegment< 22, 32 >( dataBits, 0 ).to_ulong( );

//    rampBlock.integerRampStartFrequencyModulo = convertCharactersToSignedInt32( fileBlock[ 5 ] );
//    rampBlock.fractionalRampStartFrequency = convertCharactersToSignedInt32( fileBlock[ 6 ] );

//    rampBlock.integerRampEndTime = convertCharactersToSignedInt32( fileBlock[ 7 ] );
//    rampBlock.fractionalRampEndTime = convertCharactersToSignedInt32( fileBlock[ 8 ] );

//    //rampBlock.printContents( );

//    return rampBlock;
//}

//std::map< std::pair< int, int >, std::map< double, OdfDataBlock > > getSortedOdfBlocks(
//        const std::vector< OdfDataBlock >& odfDataBlocks )
//{
//    std::map< std::pair< int, int >, std::map< double, OdfDataBlock > > dataMap;
//    for( unsigned int i = 0; i < odfDataBlocks.size( ); i++ )
//    {
//        if( //odfDataBlocks.at( i ).commonDataBlock.dataType == 11 ||
//                odfDataBlocks.at( i ).commonDataBlock.dataType == 12 ||
//                odfDataBlocks.at( i ).commonDataBlock.dataType == 13 )
//        {
//            double currentTime =
//                    static_cast< double >( odfDataBlocks.at( i ).commonDataBlock.integerTimeTag ) +
//                    static_cast< double >( odfDataBlocks.at( i ).commonDataBlock.fractionalTimeTag ) / 1000.0;
//            dataMap[ std::make_pair( odfDataBlocks.at( i ).commonDataBlock.transmittingStation,
//                                     odfDataBlocks.at( i ).commonDataBlock.receivingStation ) ][ currentTime ] =
//                    odfDataBlocks.at( i );
//        }
//    }
//    return dataMap;
//}


//std::map< std::pair< int, int >, std::map< double, double > > getDataMap( const std::vector< OdfDataBlock >& odfDataBlocks )
//{

//    std::map< std::pair< int, int >, std::map< double, double > > dataMap;
//    for( unsigned int i = 0; i < odfDataBlocks.size( ); i++ )
//    {
//        if( //odfDataBlocks.at( i ).commonDataBlock.dataType == 11 ||
//                odfDataBlocks.at( i ).commonDataBlock.dataType == 12 ||
//                odfDataBlocks.at( i ).commonDataBlock.dataType == 13 )
//        {
//            double currentTime =
//                    static_cast< double >( odfDataBlocks.at( i ).commonDataBlock.integerTimeTag ) +
//                    static_cast< double >( odfDataBlocks.at( i ).commonDataBlock.fractionalTimeTag ) / 1000.0;
//            dataMap[ std::make_pair( odfDataBlocks.at( i ).commonDataBlock.transmittingStation,
//                                     odfDataBlocks.at( i ).commonDataBlock.receivingStation ) ][ currentTime ] =
//                    static_cast< double >( odfDataBlocks.at( i ).commonDataBlock.integerObservable ) +
//                    static_cast< double >( odfDataBlocks.at( i ).commonDataBlock.fractionalObservable ) / 1.0E9;
//            //            std::cout<<currentTime<<" "<<dataMap[ odfDataBlocks.at( i ).commonDataBlock. receivingStation ][ currentTime ]<<" "<<
//            //                       odfDataBlocks.at( i ).commonDataBlock. receivingStation<<std::endl;
//        }
//    }
//    return dataMap;
//}


//int main( )
//{
//    std::ifstream dataFile( "/home/dominic/Downloads/mess_rs_15101_103_odf.dat", std::ios_base::binary);

//    char currentFileBlock[ 9 ][ 4 ];

//    readOdfFileBlock( currentFileBlock, dataFile );
//    {
//        int32_t primaryKey;
//        uint32_t secondaryKey, logicalrecordLength, groupStartPacketNumber;
//        parseHeader( currentFileBlock, primaryKey, secondaryKey, logicalrecordLength, groupStartPacketNumber );
//        std::cout<<primaryKey<<" "<<secondaryKey<<" "<<logicalrecordLength<<" "<<groupStartPacketNumber <<std::endl;
//    }

//    readOdfFileBlock( currentFileBlock, dataFile );
//    {
//        std::string systemId, programId, fileCreationDate, fileCreationTime;
//        uint32_t spacecraftIdNumber, fileReferenceDate, fileReferenceTime;
//        parseFileLabel( currentFileBlock, systemId, programId, fileCreationDate, fileCreationTime,
//                        spacecraftIdNumber, fileReferenceDate, fileReferenceTime );
//        std::cout<<" "<<systemId<<" "<<programId<<" "<<fileCreationDate<<" "<<fileCreationTime<<" "<<
//                   spacecraftIdNumber<<" "<<fileReferenceDate<<" "<<fileReferenceTime<<std::endl;



//    }

//    readOdfFileBlock( currentFileBlock, dataFile );
//    {
//        int32_t primaryKey;
//        uint32_t secondaryKey, logicalrecordLength, groupStartPacketNumber;
//        parseHeader( currentFileBlock, primaryKey, secondaryKey, logicalrecordLength, groupStartPacketNumber );
//        std::cout<<primaryKey<<" "<<secondaryKey<<" "<<logicalrecordLength<<" "<<groupStartPacketNumber <<std::endl;
//    }

//    readOdfFileBlock( currentFileBlock, dataFile );
//    parseIdentifierGroup( currentFileBlock );

//    readOdfFileBlock( currentFileBlock, dataFile );
//    {
//        int32_t primaryKey;
//        uint32_t secondaryKey, logicalrecordLength, groupStartPacketNumber;
//        parseHeader( currentFileBlock, primaryKey, secondaryKey, logicalrecordLength, groupStartPacketNumber );
//        std::cout<<primaryKey<<" "<<secondaryKey<<" "<<logicalrecordLength<<" "<<groupStartPacketNumber <<std::endl<<std::endl;
//    }

//    std::vector< OdfDataBlock > odfDataBlocks;
//    std::map< int, std::vector< OdfRampBlock > > odfRampBlocks;

//    bool continueFileRead = true;

//    int counter = 0;
//    int currentRampStation = -1;
//    int secondaryKey = -1;

//    int dataBlockType = 1; // 1: Orbit Data, 2: Ramp Data, 3: Clock Offset

//    while( continueFileRead )
//    {
//        readOdfFileBlock( currentFileBlock, dataFile );

//        int blockIsHeader = currentBlockIsHeader( currentFileBlock, secondaryKey );

//        if( blockIsHeader == 0 && dataBlockType == 1 )
//        {
//            odfDataBlocks.push_back( parseOrbitData( currentFileBlock ) );
//        }
//        else if( blockIsHeader == 0 && dataBlockType == 2 )
//        {
//            odfRampBlocks[ currentRampStation ].push_back( parseRampData( currentFileBlock ) );
//        }
//        else if( blockIsHeader != 0 )
//        {
//            //std::cout<<"Header "<<blockIsHeader<<" "<<counter<<" "<<secondaryKey<<std::endl;
//            if( blockIsHeader == 3 )
//            {
//                dataBlockType = 1;
//            }
//            else if( blockIsHeader == 4 )
//            {
//                currentRampStation = secondaryKey;
//                dataBlockType = 2;
//            }
//            else if( blockIsHeader == 5 )
//            {
//                dataBlockType = 3;
//            }
//            else if( blockIsHeader == -1 )
//            {
//                continueFileRead = 0;
//            }
//            else
//            {
//                throw std::runtime_error( "Error when reading ODF file, header not recognized" );
//            }
//        }
//        else
//        {
//            std::cout<<blockIsHeader<<" "<<dataBlockType<<std::endl;
//            throw std::runtime_error( "Error, did not recognized header ODF block" );
//        }
//        counter++;
//    }

//    std::map< int, boost::shared_ptr< RampedReferenceFrequencyInterpolator > > rampFrequencyInterpolators;
//    for( auto dataIterator = odfRampBlocks.begin( ); dataIterator != odfRampBlocks.end( ); dataIterator++ )
//    {
//        rampFrequencyInterpolators[ dataIterator->first ] =
//                boost::make_shared< RampedReferenceFrequencyInterpolator >( dataIterator->second );
//    }

//    std::map< std::pair< int, int >, std::map< double, OdfDataBlock > > sortedOdfBlocks = getSortedOdfBlocks( odfDataBlocks );
//    std::map< std::pair< int, int >, std::map< double, double > > dopplerData = getDataMap( odfDataBlocks );
//    for( auto dataIterator = dopplerData.begin( ); dataIterator != dopplerData.end( ); dataIterator++ )
//    {
//        std::map< double, double > currentDopplerData = dataIterator->second;
//        std::map< double, double > currentReferenceFrequency;
//        std::map< double, double > currentReferenceIntegrals;

//        tudat::input_output::writeDataMapToTextFile(
//                    dataIterator->second, "dopplerDataTest_" + std::to_string( dataIterator->first.first ) + "_" +
//                    std::to_string( dataIterator->first.second ) + ".dat", "/home/dominic/Documents/" );

//        for( auto dopplerIterator = currentDopplerData.begin( ); dopplerIterator != currentDopplerData.end( ); dopplerIterator++ )
//        {
//            currentReferenceFrequency[ dopplerIterator->first ] = sortedOdfBlocks.at( dataIterator->first ).at( dopplerIterator->first ).
//                    dopplerDataBlock.getReferenceFrequency( );
//        }

//        tudat::input_output::writeDataMapToTextFile(
//                    currentReferenceFrequency, "dopplerReferece_" + std::to_string( dataIterator->first.first ) + "_" +
//                    std::to_string( dataIterator->first.second ) + ".dat", "/home/dominic/Documents/" );

//        if( rampFrequencyInterpolators.count( dataIterator->first.first ) > 0 )
//        {
//            std::cout<<"Data points "<<currentDopplerData.size( )<<std::endl;
//            for( auto dopplerIterator = currentDopplerData.begin( ); dopplerIterator != currentDopplerData.end( ); dopplerIterator++ )
//            {
//                currentReferenceIntegrals[ dopplerIterator->first ] = rampFrequencyInterpolators[ dataIterator->first.first ]->getCurrentReferenceFrequency(
//                            dopplerIterator->first - 0.5, 1.0 );
//            }

//            tudat::input_output::writeDataMapToTextFile(
//                        currentReferenceIntegrals, "dopplerDataReference_" + std::to_string( dataIterator->first.first ) + "_" +
//                        std::to_string( dataIterator->first.second ) + ".dat", "/home/dominic/Documents/" );
//        }
//    }



//    std::cout<<"JD "<< tudat::basic_astrodynamics::convertCalendarDateToJulianDay(
//                   1950, 1, 1, 0, 0, 0.0 )<<" "<<
//               ( tudat::basic_astrodynamics::JULIAN_DAY_ON_J2000 - tudat::basic_astrodynamics::convertCalendarDateToJulianDay(
//                     1950, 1, 1, 0, 0, 0.0 ) ) * 86400.0<<std::endl;
//}
