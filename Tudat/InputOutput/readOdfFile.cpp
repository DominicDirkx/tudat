#include <Tudat/InputOutput/readOdfFile.h>

namespace tudat
{
namespace input_output
{

uint32_t convertCharactersToUnsignedInt32(
        char data[ 4 ] )
{
    uint8_t bytes[4];
    bytes[ 0 ] = data[ 0 ];
    bytes[ 1 ] = data[ 1 ];
    bytes[ 2 ] = data[ 2 ];
    bytes[ 3 ] = data[ 3 ];

    return bytes[0] << 24| (bytes[1] << 16) | (bytes[2] << 8) | (bytes[3]);
}

int32_t convertCharactersToSignedInt32(
        char data[ 4 ] )
{
    uint8_t bytes[4];
    bytes[ 0 ] = data[ 0 ];
    bytes[ 1 ] = data[ 1 ];
    bytes[ 2 ] = data[ 2 ];
    bytes[ 3 ] = data[ 3 ];

    int32_t returnData = bytes[0] << 24| (bytes[1] << 16) | (bytes[2] << 8) | (bytes[3]);
    return returnData;
}


boost::shared_ptr< OdfSequentialRangeDataBlock > parseSequentialRangeData( char fileBlock[ 9 ][ 4 ], const int dopplerType )
{
    boost::shared_ptr< OdfSequentialRangeDataBlock > rangeDataBlock =
            boost::make_shared< OdfSequentialRangeDataBlock >( );

    std::bitset< 32 > dataBits = std::bitset< 32 >( convertCharactersToUnsignedInt32( fileBlock[ 5 ] ) );

    rangeDataBlock->lowestRangingComponent = getBitsetSegment< 7, 32 >( dataBits, 0 ).to_ulong( );
    rangeDataBlock->spacecraftId = getBitsetSegment< 10, 32 >( dataBits, 7 ).to_ulong( );
    rangeDataBlock->reservedBlock = getBitsetSegment< 1, 32 >( dataBits, 17 ).to_ulong( );

    std::bitset< 14 > referenceFrequencyHighPartBitsA = getBitsetSegment< 14, 32 >( dataBits, 18 );
    dataBits = std::bitset< 32 >( convertCharactersToUnsignedInt32( fileBlock[ 6 ] ) );
    std::bitset< 8 > referenceFrequencyHighPartBitsB = getBitsetSegment< 8, 32 >( dataBits, 0 );
    rangeDataBlock->referenceFrequencyHighPart = mergeBitsets< 14, 8 >(
                referenceFrequencyHighPartBitsA, referenceFrequencyHighPartBitsB ).to_ulong();
    rangeDataBlock->referenceFrequencyLowPart = getBitsetSegment< 24, 32 >( dataBits, 8 ).to_ulong( );

    dataBits = std::bitset< 32 >( convertCharactersToUnsignedInt32( fileBlock[ 7 ] ) );
    rangeDataBlock->coderInPhaseTimeOffset = getSignedNBitInteger< 20 >( getBitsetSegment< 20, 32 >( dataBits, 0 ) );

    std::bitset< 12 > compositeTwoPartA = getBitsetSegment< 12, 32 >( dataBits, 20 );
    dataBits = std::bitset< 32 >( convertCharactersToUnsignedInt32( fileBlock[ 8 ] ) );
    std::bitset< 10 > compositeTwoPartB = getBitsetSegment< 10, 32 >( dataBits, 0 );
    rangeDataBlock->compositeTwo = mergeBitsets< 12, 10 >(
                compositeTwoPartA, compositeTwoPartB ).to_ulong();
    rangeDataBlock->transmittingStationUplinkDelay = getBitsetSegment< 22, 32 >( dataBits, 10 ).to_ulong( );

    return rangeDataBlock;
}

boost::shared_ptr< OdfDopplerDataBlock > parseDopplerOrbitData( char fileBlock[ 9 ][ 4 ], const int dopplerType )
{
    boost::shared_ptr< OdfDopplerDataBlock > dopplerDataBlock =
            boost::make_shared< OdfDopplerDataBlock >( dopplerType );

    std::bitset< 32 > dataBits = std::bitset< 32 >( convertCharactersToUnsignedInt32( fileBlock[ 5 ] ) );

    dopplerDataBlock->receiverChannel = getBitsetSegment< 7, 32 >( dataBits, 0 ).to_ulong( );
    dopplerDataBlock->spacecraftId = getBitsetSegment< 10, 32 >( dataBits, 7 ).to_ulong( );
    dopplerDataBlock->receiverExciterFlag = getBitsetSegment< 1, 32 >( dataBits, 17 ).to_ulong( );

    std::bitset< 14 > referenceFrequencyHighPartBitsA = getBitsetSegment< 14, 32 >( dataBits, 18 );
    dataBits = std::bitset< 32 >( convertCharactersToUnsignedInt32( fileBlock[ 6 ] ) );
    std::bitset< 8 > referenceFrequencyHighPartBitsB = getBitsetSegment< 8, 32 >( dataBits, 0 );
    dopplerDataBlock->referenceFrequencyHighPart = mergeBitsets< 14, 8 >(
                referenceFrequencyHighPartBitsA, referenceFrequencyHighPartBitsB ).to_ulong();
    dopplerDataBlock->referenceFrequencyLowPart = getBitsetSegment< 24, 32 >( dataBits, 8 ).to_ulong( );

    dataBits = std::bitset< 32 >( convertCharactersToUnsignedInt32( fileBlock[ 7 ] ) );

    dopplerDataBlock->reservedSegment = getBitsetSegment< 20, 32 >( dataBits, 0 ).to_ulong( );
    std::bitset< 12 > compressionTimeABits = getBitsetSegment< 12, 32 >( dataBits, 20 );
    dataBits = std::bitset< 32 >( convertCharactersToUnsignedInt32( fileBlock[ 8 ] ) );
    std::bitset< 10 > compressionTimeBBits = getBitsetSegment< 10, 32 >( dataBits, 0 );
    dopplerDataBlock->compressionTime = mergeBitsets< 12, 10 >(
                compressionTimeABits, compressionTimeBBits ).to_ulong( );
    dopplerDataBlock->transmittingStationDelay = getBitsetSegment< 22, 32 >( dataBits, 10 ).to_ulong( );

    return dopplerDataBlock;
}


boost::shared_ptr< OdfDataBlock > parseOrbitData( char fileBlock[ 9 ][ 4 ] )
{
    boost::shared_ptr< OdfCommonDataBlock > commonDataBlock =
            boost::make_shared< OdfCommonDataBlock >( );

    commonDataBlock->integerTimeTag = convertCharactersToUnsignedInt32( fileBlock[ 0 ] );

    std::bitset< 32 > dataBits = std::bitset< 32 >( convertCharactersToUnsignedInt32( fileBlock[ 1 ] ) );
    commonDataBlock->fractionalTimeTag = getBitsetSegment< 10, 32 >(
                dataBits, 0 ).to_ulong( );
    commonDataBlock->receivingStationDownlinkDelay = getBitsetSegment< 22, 32 >(
                dataBits, 10 ).to_ulong( );

    commonDataBlock->integerObservable = convertCharactersToSignedInt32( fileBlock[ 2 ] );
    commonDataBlock->fractionalObservable = convertCharactersToSignedInt32( fileBlock[ 3 ] );


    dataBits = std::bitset< 32 >( convertCharactersToUnsignedInt32( fileBlock[ 4 ] ) );


    commonDataBlock->formatId = getBitsetSegment< 3, 32 >( dataBits, 0 ).to_ulong( );
    commonDataBlock->receivingStation = getBitsetSegment< 7, 32 >( dataBits, 3 ).to_ulong( );
    commonDataBlock->transmittingStation = getBitsetSegment< 7, 32 >( dataBits, 10 ).to_ulong( );

    commonDataBlock->transmittingStationNetworkId = getBitsetSegment< 2, 32 >( dataBits, 17 ).to_ulong( );
    int dataType =  getBitsetSegment< 6, 32 >( dataBits, 19 ).to_ulong( );

    commonDataBlock->downlinkBand = getBitsetSegment< 2, 32 >( dataBits, 25 ).to_ulong( );
    commonDataBlock->uplinkBand = getBitsetSegment< 2, 32 >( dataBits, 27 ).to_ulong( );
    commonDataBlock->referenceBand = getBitsetSegment< 2, 32 >( dataBits, 29 ).to_ulong( );
    commonDataBlock->validity = getBitsetSegment< 1, 32 >( dataBits, 31 ).to_ulong( );

    boost::shared_ptr< OdfDataBlock > dataBlock = boost::make_shared< OdfDataBlock >( );


    dataBlock->commonDataBlock = commonDataBlock;
    if(  dataType == 11 || dataType == 12 || dataType == 13 )
    {
        if( dataType != 11 )
        {
            std::cout<<commonDataBlock->integerObservable<<" "<<std::bitset< 8 >( fileBlock[ 2 ][ 0 ] )
                    <<" "<<std::bitset< 8 >( fileBlock[ 2 ][ 1 ] )<<" "<<std::bitset< 8 >( fileBlock[ 2 ][ 2 ] )<<" "<<std::bitset< 8 >( fileBlock[ 2 ][ 3 ] )
                    <<" "<<std::bitset< 32 >( convertCharactersToUnsignedInt32( fileBlock[ 2 ] ) )<<" "<<std::endl;
        }
        dataBlock->observableSpecificDataBlock = parseDopplerOrbitData( fileBlock, dataType );
    }
    else if(  dataType == 37 )
    {
        dataBlock->observableSpecificDataBlock = parseSequentialRangeData( fileBlock, dataType );
    }
    //    else if( dataType == 14 || dataType == 7 || dataType == 2 || dataType == 40  || dataType == 50 || dataType == 20 )
    //    {
    //        dataBlock = NULL;
    //    }
    else
    {
        std::cerr<<"Data type "<<dataType<<std::endl;
        dataBlock = NULL;
        //throw std::runtime_error( "Error, ODF data type " + std::to_string( dataType ) + " not recognized" );

    }
    return dataBlock;
}

OdfRampBlock parseRampData( char fileBlock[ 9 ][ 4 ] )
{

    OdfRampBlock rampBlock;
    rampBlock.integerRampStartTime = convertCharactersToUnsignedInt32( fileBlock[ 0 ] );
    rampBlock.fractionalRampStartTime = convertCharactersToUnsignedInt32( fileBlock[ 1 ] );

    rampBlock.integerRampRate = convertCharactersToSignedInt32( fileBlock[ 2 ] );
    rampBlock.fractionalRampRate = convertCharactersToSignedInt32( fileBlock[ 3 ] );

    std::bitset< 32 > dataBits = std::bitset< 32 >( convertCharactersToUnsignedInt32( fileBlock[ 4 ] ) );
    rampBlock.transmittingStationId = getBitsetSegment< 10, 32 >( dataBits, 22 ).to_ulong( );
    rampBlock.integerRampStartFrequency = getBitsetSegment< 22, 32 >( dataBits, 0 ).to_ulong( );

    rampBlock.integerRampStartFrequencyModulo = convertCharactersToSignedInt32( fileBlock[ 5 ] );
    rampBlock.fractionalRampStartFrequency = convertCharactersToSignedInt32( fileBlock[ 6 ] );

    rampBlock.integerRampEndTime = convertCharactersToSignedInt32( fileBlock[ 7 ] );
    rampBlock.fractionalRampEndTime = convertCharactersToSignedInt32( fileBlock[ 8 ] );

    //std::cout<<rampBlock.integerRampStartFrequency<<" "<<rampBlock.integerRampStartFrequencyModulo<<" "<<rampBlock.fractionalRampStartFrequency<<std::endl;

    //rampBlock.printContents( );

    return rampBlock;
}

void parseFileLabel( char fileBlock[ 9 ][ 4 ],
std::string& systemId, std::string& programId, std::string& fileCreationDate, std::string& fileCreationTme,
uint32_t& spacecraftIdNumber, uint32_t& fileReferenceDate, uint32_t& fileReferenceTime )
{
    systemId.resize( 8 );
    programId.resize( 8 );

    fileCreationDate.resize( 4 );
    fileCreationTme.resize( 4 );

    for( unsigned int i = 0; i < 4; i++ )
    {
        systemId[ i ] = fileBlock[ 0 ][ i ];
        systemId[ i + 4 ] = fileBlock[ 1 ][ i ];

        programId[ i ] = fileBlock[ 2 ][ i ];
        programId[ i + 4 ] = fileBlock[ 3 ][ i ];

        fileCreationDate[ i ] = fileBlock[ 5 ][ i ];
        fileCreationTme[ i ] = fileBlock[ 6 ][ i ];
    }
    spacecraftIdNumber = convertCharactersToUnsignedInt32( fileBlock[ 4 ] );
    fileReferenceDate = convertCharactersToUnsignedInt32( fileBlock[ 7 ] );
    fileReferenceTime = convertCharactersToUnsignedInt32( fileBlock[ 8 ] );
}

void parseHeader( char fileBlock[ 9 ][ 4 ],
int32_t& primaryKey,
uint32_t& secondaryKey,
uint32_t& logicalrecordLength,
uint32_t& groupStartPacketNumber )
{
    primaryKey = convertCharactersToSignedInt32( fileBlock[ 0 ] );
    secondaryKey = convertCharactersToUnsignedInt32( fileBlock[ 1 ] );
    logicalrecordLength = convertCharactersToUnsignedInt32( fileBlock[ 2 ] );
    groupStartPacketNumber = convertCharactersToUnsignedInt32( fileBlock[ 3 ] );
    uint32_t testInteger;

    for( int i = 0; i < 5; i++ )
    {
        testInteger = convertCharactersToUnsignedInt32( fileBlock[ i + 4 ] );
        if( testInteger != 0 )
        {
            throw std::runtime_error( "Error when reading ODF file, header file inconsistent" );
        }
    }
}


int currentBlockIsHeader(  char fileBlock[ 9 ][ 4 ], unsigned int& secondaryKeyInt )
{
    int headerType = -2;

    int32_t primaryKey;
    uint32_t secondaryKey;
    uint32_t logicalrecordLength;
    uint32_t groupStartPacketNumber;

    try
    {
        parseHeader( fileBlock, primaryKey, secondaryKey, logicalrecordLength, groupStartPacketNumber );
        secondaryKeyInt = secondaryKey;
        if( primaryKey ==  101 )
        {
            headerType = 1;
        }
        else if( primaryKey == 107 )
        {
            headerType = 2;
        }
        else if( primaryKey == 109 )
        {
            headerType = 3;
        }
        else if( primaryKey == 2030 )
        {
            headerType = 4;
        }
        else if( primaryKey == 2040 )
        {
            headerType = 5;
        }
        else if( primaryKey == -1 )
        {
            headerType = -1;
        }
    }
    catch( std::runtime_error )
    {
        headerType = 0;
    }
    return headerType;
}

void readOdfFileBlock(
        char fileBlock[ 9 ][ 4 ],
std::istream& file )
{
    for( unsigned int i = 0; i < 9; i++ )
    {
        file.read( (char*)fileBlock[ i ], 4 );
    }
}

void parseIdentifierGroup( char fileBlock[ 9 ][ 4 ] )
{
    std::string testStringA, testStringB, testStringC;
    testStringA.resize( 8 );
    testStringB.resize( 8 );
    testStringC.resize( 20 );

    for( unsigned int i = 0; i < 4; i++ )
    {
        testStringA[ i ] = fileBlock[ 0 ][ i ];
        testStringA[ i + 4 ] = fileBlock[ 1 ][ i ];

        testStringB[ i ] = fileBlock[ 2 ][ i ];
        testStringB[ i + 4 ] = fileBlock[ 3 ][ i ];

        for( unsigned int j = 0; j < 5; j++ )
        {
            testStringC[ i + j * 4 ] = fileBlock[ 4 + j ][ i ];
        }
    }

    std::cout<<testStringA<<" "<<testStringB<<" "<<testStringC<<std::endl;
    //    if( testStringA != "TIMETAG " || testStringB != "OBSRVBL " || testStringC != "FREQ, ANCILLARY-DATA"  )
    //    {
    //        throw std::runtime_error( "Error when reading ODF file, identifier group incompatible" );
    //    }
}

boost::shared_ptr< OdfRawFileContents > readOdfFile(
        const std::string& odfFile )
{
    std::ifstream dataFile( odfFile, std::ios_base::binary);

    boost::shared_ptr< OdfRawFileContents > odfFileContents =
            boost::make_shared< OdfRawFileContents >( );

    char currentFileBlock[ 9 ][ 4 ];

    readOdfFileBlock( currentFileBlock, dataFile );

    int32_t primaryKey;
    uint32_t secondaryKey, logicalrecordLength, groupStartPacketNumber;
    parseHeader( currentFileBlock, primaryKey, secondaryKey, logicalrecordLength, groupStartPacketNumber );
    if( primaryKey != 101 || secondaryKey != 0 || logicalrecordLength != 1 || groupStartPacketNumber != 0 )
    {
        throw std::runtime_error( "Error when reding ODF file, file header incompatible" );
    }


    readOdfFileBlock( currentFileBlock, dataFile );
    parseFileLabel( currentFileBlock, odfFileContents->systemId, odfFileContents->programId,
                    odfFileContents->fileCreationDate, odfFileContents->fileCreationTime,
                    odfFileContents->spacecraftId, odfFileContents->fileReferenceDate, odfFileContents->fileReferenceTime );

    readOdfFileBlock( currentFileBlock, dataFile );
    parseHeader( currentFileBlock, primaryKey, secondaryKey, logicalrecordLength, groupStartPacketNumber );
    if( primaryKey != 107 || secondaryKey != 0 || logicalrecordLength != 1 || groupStartPacketNumber != 2 )
    {
        throw std::runtime_error( "Error when reding ODF file, file header incompatible" );
    }

    readOdfFileBlock( currentFileBlock, dataFile );
    parseIdentifierGroup( currentFileBlock );

    readOdfFileBlock( currentFileBlock, dataFile );
    parseHeader( currentFileBlock, primaryKey, secondaryKey, logicalrecordLength, groupStartPacketNumber );
    if( primaryKey != 109 || secondaryKey != 0 || logicalrecordLength != 1 || groupStartPacketNumber != 4 )
    {
        throw std::runtime_error( "Error when reding ODF file, file header incompatible" );
    }

    std::vector< boost::shared_ptr< OdfDataBlock > > unsortedOdfDataBlocks;
    std::map< int, std::vector< OdfRampBlock > > odfRampBlocks;

    bool continueFileRead = true;

    int counter = 0;
    int currentRampStation = -1;

    int dataBlockType = 3; // 1: Orbit Data, 2: Ramp Data, 3: Clock Offset

    while( continueFileRead )
    {
        readOdfFileBlock( currentFileBlock, dataFile );
        int blockIsHeader = currentBlockIsHeader( currentFileBlock, secondaryKey );

        if( blockIsHeader == 0 && dataBlockType == 3 )
        {
            boost::shared_ptr< OdfDataBlock > currentDataBlock = parseOrbitData( currentFileBlock );
            if( currentDataBlock != NULL )
            {
                unsortedOdfDataBlocks.push_back( currentDataBlock );
            }
        }
        else if( blockIsHeader == 0 && dataBlockType == 4 )
        {
            odfRampBlocks[ currentRampStation ].push_back( parseRampData( currentFileBlock ) );
        }
        else if( blockIsHeader != 0 )
        {
            dataBlockType = blockIsHeader;
            if( dataBlockType == -1 )
            {
                continueFileRead = 0;
            }
        }
        else
        {
            throw std::runtime_error( "Error, did not recognized header ODF block" );
        }
        if( dataFile.eof( ) )
        {
            continueFileRead = 0;
        }
        counter++;
    }

    odfFileContents->dataBlocks = unsortedOdfDataBlocks;
    odfFileContents->odfRampBlocks = odfRampBlocks;

    return odfFileContents;
}

} // namespace input_output

} // namespace tudat

