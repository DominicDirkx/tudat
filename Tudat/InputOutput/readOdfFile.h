/*    Copyright (c) 2010-2018, Delft University of Technology
 *    All rigths reserved
 *
 *    This file is part of the Tudat. Redistribution and use in source and
 *    binary forms, with or without modification, are permitted exclusively
 *    under the terms of the Modified BSD license. You should have received
 *    a copy of the license with this file. If not, please or visit:
 *    http://tudat.tudelft.nl/LICENSE.
 *
 *    Notes
 *      The function printStandardScientificNotation() has been implemented to cope with
 *      cross-platform incompatibilities in the printed output of floating-point numbers in
 *      scientific notation.
 *
 */

#ifndef TUDAT_READ_ODF_FILE_H
#define TUDAT_READ_ODF_FILE_H

#include <iostream>
#include <fstream>
#include <bitset>
#include <cmath>
#include <vector>
#include <map>

#include <boost/make_shared.hpp>

#include "Tudat/Astrodynamics/BasicAstrodynamics/timeConversions.h"
#include "Tudat/InputOutput/basicInputOutput.h"
#include "Tudat/Mathematics/Interpolators/lookupScheme.h"

namespace tudat
{
namespace input_output
{

template< int FirstInputSize, int SecondInputSize >
std::bitset< FirstInputSize + SecondInputSize > mergeBitsets(
        std::bitset< FirstInputSize > firstInput,
        std::bitset< SecondInputSize > secondInput )
{
    std::bitset< FirstInputSize + SecondInputSize > returnBitset;
    for( int i = 0; i < FirstInputSize; i++ )
    {
        returnBitset[ i + SecondInputSize ] = firstInput[ i ];
    }

    for( int i = 0; i < SecondInputSize; i++ )
    {
        returnBitset[ i ] = secondInput[ i ];
    }
    return returnBitset;
}

template< int OutputBits, int InputBits >
std::bitset< OutputBits > getBitsetSegment(
        const std::bitset< InputBits > inputBits,
        const int startIndex )
{
    std::bitset< OutputBits > outputBits;
    for( unsigned int i = 0; i < OutputBits; i++ )
    {
        outputBits[ i ] = inputBits[ InputBits - OutputBits - startIndex + i  ];
    }
    return outputBits;
}

template< int NumberOfBits >
int getSignedNBitInteger(
        std::bitset< NumberOfBits > inputBits )
{
    int outputInteger = -inputBits[ NumberOfBits - 1 ] * std::pow( 2, NumberOfBits - 1 );

    for( unsigned int i = 0; i < NumberOfBits - 1; i ++ )
    {
        outputInteger += inputBits[ i ] * std::pow( 2.0, i );
    }
    return outputInteger;
}

uint32_t convertCharactersToUnsignedInt32(
        unsigned char data[ 4 ] );

int32_t convertCharactersToSignedInt32(
        unsigned char data[ 4 ] );

class OdfRampBlock
{
public:

    int integerRampStartTime;
    int fractionalRampStartTime;

    int integerRampEndTime;
    int fractionalRampEndTime;

    int integerRampRate;
    int fractionalRampRate;

    int integerRampStartFrequency;
    int integerRampStartFrequencyModulo;
    int fractionalRampStartFrequency;

    int transmittingStationId;

    double getRampStartFrequency( )
    {
        return static_cast< double >( integerRampStartFrequency ) * 1.0E9 +
                static_cast< double >( integerRampStartFrequencyModulo ) +
                static_cast< double >( fractionalRampStartFrequency ) * 1.0E-9;
    }

    double getRampRate( )
    {
        return static_cast< double >( integerRampRate ) +
                static_cast< double >( fractionalRampRate ) * 1.0E-9;
    }

    double getRampStartTime( )
    {
        return static_cast< double >( integerRampStartTime ) -
                86400.0 + static_cast< double >( fractionalRampStartTime ) * 1.0E-9;
    }

    double getRampEndTime( )
    {
        return static_cast< double >( integerRampEndTime ) -
                static_cast< double >( fractionalRampEndTime ) * 1.0E-9;
    }

    void printContents( )
    {
        std::cout<<"Start time "<<integerRampStartTime<<" "<<fractionalRampStartTime<<std::endl;
        std::cout<<"End time "<<integerRampEndTime<<" "<<fractionalRampEndTime<<std::endl;
        std::cout<<"Ramp rate "<<integerRampRate<<" "<<fractionalRampRate<<std::endl;
        std::cout<<"Start frequency "<<integerRampStartFrequency<<" "<<integerRampStartFrequencyModulo<<" "<<
                   fractionalRampStartFrequency<<std::endl;
        std::cout<<"Station "<<transmittingStationId<<std::endl<<std::endl;
    }
};

class OdfDataSpecificBlock
{
public:
    OdfDataSpecificBlock( int dataType_ ): dataType( dataType_ ){ }

    virtual ~OdfDataSpecificBlock( ){ }

    int dataType;
};

class OdfSequentialRangeDataBlock: public OdfDataSpecificBlock
{
public:

    OdfSequentialRangeDataBlock( ):
        OdfDataSpecificBlock( 37 ){ }

    ~OdfSequentialRangeDataBlock( ){ }

    int lowestRangingComponent;
    int spacecraftId;
    int reservedBlock;

    int referenceFrequencyHighPart;
    int referenceFrequencyLowPart;

    int coderInPhaseTimeOffset;
    int compositeTwo;
    int transmittingStationUplinkDelay;

};

class OdfDopplerDataBlock: public OdfDataSpecificBlock
{
public:
    OdfDopplerDataBlock( const int DopplerDataType ):
        OdfDataSpecificBlock( DopplerDataType ){ }

    ~OdfDopplerDataBlock( ){ }

    int receiverChannel;
    int spacecraftId;
    int receiverExciterFlag;
    int referenceFrequencyHighPart;
    int referenceFrequencyLowPart;

    int reservedSegment;
    int compressionTime;

    int transmittingStationDelay;

    double getReferenceFrequency( )
    {
        return std::pow( 2.0, 24 )  /1.0E3 * referenceFrequencyHighPart + referenceFrequencyLowPart / 1.0E3;
    }

    void printContents( )
    {
        std::cout<<"Receiver: "<<receiverChannel<<" "<<receiverExciterFlag<<std::endl;
        std::cout<<"Spacecraft: "<<spacecraftId<<std::endl;
        std::cout<<"Reference frequency: "<<referenceFrequencyHighPart<<" "<<referenceFrequencyLowPart<<std::endl;

        std::cout<<"Reserved: "<<reservedSegment<<std::endl;
        std::cout<<"Compression time: "<<compressionTime<<std::endl;
        std::cout<<"Transmission delay: "<<transmittingStationDelay<<std::endl<<std::endl;;

    }

};

class OdfCommonDataBlock
{
public:
    int integerTimeTag;
    int fractionalTimeTag;

    int integerObservable;
    int fractionalObservable;

    double getObservableValue( )
    {
        return static_cast< double >( integerObservable ) +
                static_cast< double >( fractionalObservable ) / 1.0E9;
    }

    double getObservableTime( )
    {
        return static_cast< double >( integerTimeTag ) + static_cast< double >( fractionalTimeTag ) / 1000.0;
    }

    int receivingStationDownlinkDelay;
    int formatId;
    int receivingStation;
    int transmittingStation;
    int transmittingStationNetworkId;

    int downlinkBand;
    int uplinkBand;
    int referenceBand;
    int validity;

    void printContents( )
    {
        std::cout<<"Time: "<<integerTimeTag<<" "<<fractionalTimeTag<<std::endl;
        std::cout<<"Downlink delay: "<<receivingStationDownlinkDelay<<std::endl;
        std::cout<<"Observable: "<<integerObservable<<" "<<fractionalObservable<<std::endl;

        std::cout<<"Format id: "<<formatId<<std::endl;
        std::cout<<"Station data: "<<receivingStation<<" "<<transmittingStation<<" "<<transmittingStationNetworkId<<" "<<std::endl;
        std::cout<<"Bands: "<<downlinkBand<<" "<<uplinkBand<<" "<<referenceBand<<" "<<std::endl;
        std::cout<<"Validity: "<<validity<<std::endl<<std::endl;

    }
};

class OdfDataBlock
{
public:
    boost::shared_ptr< OdfDataSpecificBlock > observableSpecificDataBlock;
    boost::shared_ptr< OdfCommonDataBlock > commonDataBlock;
};

class OdfRawFileContents
{
public:
    std::string systemId;
    std::string programId;
    uint32_t spacecraftId;

    std::string fileCreationDate;
    std::string fileCreationTime;

    uint32_t fileReferenceDate;
    uint32_t fileReferenceTime;

    std::string fileName;

    std::string identifierGroupStringA;
    std::string identifierGroupStringB;
    std::string identifierGroupStringC;

    bool eofHeaderFound;

    std::vector< boost::shared_ptr< OdfDataBlock > > dataBlocks;
    std::map< int, std::vector< OdfRampBlock > > odfRampBlocks;

};

//! Function to parse the contents of an ODF orbit data block, specific for sequenctial range data.
boost::shared_ptr< OdfSequentialRangeDataBlock > parseSequentialRangeData( unsigned char fileBlock[ 9 ][ 4 ], const int dopplerType );

//! Function to parse the contents of an ODF orbit data block, specific for Doppler data.
boost::shared_ptr< OdfDopplerDataBlock > parseDopplerOrbitData( unsigned char fileBlock[ 9 ][ 4 ], const int dopplerType );

//! Function to parse the contents of an ODF orbit data block
boost::shared_ptr< OdfDataBlock > parseOrbitData( unsigned char fileBlock[ 9 ][ 4 ] );

//! Function to parse the contents of an ODF ramp data block
OdfRampBlock parseRampData( unsigned char fileBlock[ 9 ][ 4 ] );

//! Function to parse the contents of an ODF file label block
void parseFileLabel( unsigned char fileBlock[ 9 ][ 4 ],
std::string& systemId, std::string& programId, std::string& fileCreationDate, std::string& fileCreationTme,
uint32_t& spacecraftIdNumber, uint32_t& fileReferenceDate, uint32_t& fileReferenceTime );

//! Function to parse the contents of an ODF file header block
void parseHeader( unsigned char fileBlock[ 9 ][ 4 ],
int32_t& primaryKey,
uint32_t& secondaryKey,
uint32_t& logicalrecordLength,
uint32_t& groupStartPacketNumber );

//! Function to check if the current ODF data block is a header.
int currentBlockIsHeader(  unsigned char fileBlock[ 9 ][ 4 ], unsigned int& secondaryKeyInt );

//! Function to read a single 36 byte block from ODF file
void readOdfFileBlock(
        unsigned char fileBlock[ 9 ][ 4 ], std::istream& file );

//! Function to read the contents of an ODF file into an OdfRawFileContents object
/*!
 * Function to read the contents of an ODF file into an OdfRawFileContents object. The OdfRawFileContents object contains the
 * unprocessed contents of teh ODF file, on a line-by-line basis.
 * \param odfFile File name/location of ODF file that is to be read
 * \return OdfRawFileContents object with contents of ODF file.
 */
boost::shared_ptr< OdfRawFileContents > readOdfFile(
        const std::string& odfFile );

} // namespace input_output

} // namespace tudat

#endif // TUDAT_READ_ODF_FILE_H

