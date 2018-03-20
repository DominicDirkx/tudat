#include <iostream>
#include <fstream>
#include <bitset>
#include <cmath>

uint32_t read_u32_le(std::istream& file)
{
    uint32_t value;
    uint8_t  bytes[4];

    file.read( (char*)bytes, 4);
    value = bytes[0] << 24 | (bytes[1] << 16) | (bytes[2] << 8) | (bytes[3]);

    //std::cout<<value<<" "<<int( bytes[0] )<<" "<<int( bytes[1] << 8 )<<" "<<int( bytes[2] << 16 )<<" "<<( bytes[3] << 0 )<<std::endl;
    return value;
}

int32_t read_s32_le(std::istream& file)
{
    int32_t value;
    int8_t  bytes[4];

    file.read( (char*)bytes, 4);
    value = bytes[0] << 24| (bytes[1] << 16) | (bytes[2] << 8) | (bytes[3]);
    //std::cout<<value<<" "<<int( bytes[0] )<<" "<<int( bytes[1] << 8 )<<" "<<int( bytes[2] << 16 )<<" "<<( bytes[3] << 0 )<<std::endl;

    return value;
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
unsigned int getUnsignedNBitInteger(
        std::bitset< NumberOfBits > inputBits )
{
    int outputInteger = 0;
    for( unsigned int i = 0; i < NumberOfBits ; i ++ )
    {
        outputInteger += inputBits[ i ] * std::pow( 2.0, i );
    }
    return outputInteger;
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


void readOdfFileBlock(
        char fileBlock[ 9 ][ 4 ],
std::istream& file )
{
    for( unsigned int i = 0; i < 9; i++ )
    {
        file.read( (char*)fileBlock[ i ], 4 );
    }
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
}

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

void parseDopplerOrbitData( char fileBlock[ 9 ][ 4 ] )
{
    std::bitset< 32 > dataBits = std::bitset< 32 >( convertCharactersToUnsignedInt32( fileBlock[ 5 ] ) );
    std::bitset< 7 > receiverChannelBitsBits = getBitsetSegment< 7, 32 >( dataBits, 0 );
    std::bitset< 10 > spacecraftIdBits = getBitsetSegment< 10, 32 >( dataBits, 7 );
    std::bitset< 1 > receiverExciterFlagBits = getBitsetSegment< 1, 32 >( dataBits, 17 );
    std::bitset< 14 > referenceFrequencyHighPartBitsA = getBitsetSegment< 14, 32 >( dataBits, 18 );

    dataBits = std::bitset< 32 >( convertCharactersToUnsignedInt32( fileBlock[ 6 ] ) );
    std::bitset< 8 > referenceFrequencyHighPartBitsB = getBitsetSegment< 8, 32 >( dataBits, 0 );
    std::bitset< 22 > referenceFrequencyHighPartBits = mergeBitsets< 14, 8 >(
                referenceFrequencyHighPartBitsA, referenceFrequencyHighPartBitsB );

    std::bitset< 24 > referenceFrequencyLowPartBits = getBitsetSegment< 24, 32 >( dataBits, 8 );

    std::cout<<receiverChannelBitsBits.to_ulong( )<<" "<<spacecraftIdBits.to_ulong( )<<" "
               <<receiverExciterFlagBits.to_ulong( )<<" "<<referenceFrequencyHighPartBits.to_ulong( )<<" "<<
                 referenceFrequencyLowPartBits.to_ulong( )<<std::endl;

    dataBits = std::bitset< 32 >( convertCharactersToUnsignedInt32( fileBlock[ 7 ] ) );
    std::bitset< 20 > reservedSegmentBits = getBitsetSegment< 20, 32 >( dataBits, 0 );
    std::bitset< 12 > compressionTimeABits = getBitsetSegment< 12, 32 >( dataBits, 20 );

    std::cout<<"Data 7 "<<dataBits<<std::endl;

    dataBits = std::bitset< 32 >( convertCharactersToUnsignedInt32( fileBlock[ 8 ] ) );
    std::bitset< 10 > compressionTimeBBits = getBitsetSegment< 10, 32 >( dataBits, 0 );

    std::cout<<"Data 8 "<<dataBits<<std::endl;

    std::cout<<"Data comp. time "<<compressionTimeABits<<" "<<compressionTimeBBits <<std::endl;

    std::bitset< 22 > compressionTimeBits = mergeBitsets< 12, 10 >(
                compressionTimeABits, compressionTimeBBits );

    std::cout<<compressionTimeBits<<" "<<compressionTimeBits.to_ulong( )<<std::endl;

    std::bitset< 22 > transmittingStationDelayBits = getBitsetSegment< 22, 32 >( dataBits, 10 );

    std::cout<<transmittingStationDelayBits<<" "<<transmittingStationDelayBits.to_ulong( )<<std::endl;




}

struct OdfDopplerDataBlock
{

};

struct OdfCommonDataBlock
{
    int integerTimeTag;
    int fractionalTimeTag;
    int receivingStationDownlinkDelay;

    int integerObservable;
    int fractionalObservable;

    int formatId;
    int receivingStation;
    int transmittingStation;
    int transmittingStationNetworkId;
    int dataType;

    int downlinkBand;
    int uplinkBand;
    int referenceBand;
    int validity;
};

struct OdfDataBlock
{
    OdfDopplerDataBlock dopplerDataBlock;
    OdfCommonDataBlock commonDataBlock;
};


void parseOrbitData( char fileBlock[ 9 ][ 4 ] )
{
    OdfCommonDataBlock commonDataBlock;
    uint32_t integerTimeTag = convertCharactersToUnsignedInt32( fileBlock[ 0 ] );
    std::bitset< 32 > dataBits = std::bitset< 32 >( convertCharactersToUnsignedInt32( fileBlock[ 1 ] ) );

    std::cout<<integerTimeTag<<" "<<dataBits<<std::endl;

    std::bitset< 10 > fractionalTimeTagBits = getBitsetSegment< 10, 32 >(
                    dataBits, 0 );
    int fractionalTimeTag = getUnsignedNBitInteger< 10 >( fractionalTimeTagBits );

    std::bitset< 22 > receivingStationDownlinkDelayBits = getBitsetSegment< 22, 32 >(
                    dataBits, 10 );
    int receivingStationDownlinkDelay = getUnsignedNBitInteger< 22 >( receivingStationDownlinkDelayBits );

    std::cout<<fractionalTimeTagBits<<" "<<fractionalTimeTag<<std::endl;
    std::cout<<receivingStationDownlinkDelayBits<<" "<<receivingStationDownlinkDelay<<std::endl;

    int32_t integerObservable = convertCharactersToSignedInt32( fileBlock[ 2 ] );
    int32_t fractionalObservable = convertCharactersToSignedInt32( fileBlock[ 3 ] );

    std::cout<<integerObservable<<" "<<fractionalObservable<<std::endl;

    dataBits = std::bitset< 32 >( convertCharactersToUnsignedInt32( fileBlock[ 4 ] ) );
    std::bitset< 3 > formatIdBits = getBitsetSegment< 3, 32 >( dataBits, 0 );
    std::bitset< 7 > receivingStationBits = getBitsetSegment< 7, 32 >( dataBits, 3 );
    std::bitset< 7 > transmittingStationBits = getBitsetSegment< 7, 32 >( dataBits, 10 );
    std::bitset< 2 > transmittingStationNetworkIdBits = getBitsetSegment< 2, 32 >( dataBits, 17 );
    std::bitset< 6 > dataTypeBits = getBitsetSegment< 6, 32 >( dataBits, 19 );

    int formatId = getUnsignedNBitInteger< 3 >( formatIdBits );
    int receivingStation = getUnsignedNBitInteger< 7 >( receivingStationBits );
    int transmittingStation = getUnsignedNBitInteger< 7 >( transmittingStationBits );
    int transmittingStationNetworkId = getUnsignedNBitInteger< 2 >( transmittingStationNetworkIdBits );
    int dataType = getUnsignedNBitInteger< 6 >( dataTypeBits );

    std::cout<<dataBits<<std::endl;
    std::cout<<formatIdBits<<" "<<formatId<<std::endl;
    std::cout<<receivingStationBits<<" "<<receivingStation<<std::endl;
    std::cout<<transmittingStationBits<<" "<<transmittingStation<<std::endl;
    std::cout<<transmittingStationNetworkIdBits<<" "<<transmittingStationNetworkId<<std::endl;
    std::cout<<dataTypeBits<<" "<<dataType<<std::endl;

    std::bitset< 2 > downlinkBandBits = getBitsetSegment< 2, 32 >( dataBits, 25 );
    std::bitset< 2 > uplinkBandBits = getBitsetSegment< 2, 32 >( dataBits, 27 );
    std::bitset< 2 > referenceBandBits = getBitsetSegment< 2, 32 >( dataBits, 29 );
    std::bitset< 1 > validityBits = getBitsetSegment< 1, 32 >( dataBits, 31 );

    std::cout<<getUnsignedNBitInteger< 2 >( downlinkBandBits )<<" "<<
               getUnsignedNBitInteger< 2 >( uplinkBandBits )<<" "<<
               getUnsignedNBitInteger< 2 >( referenceBandBits )<<" "<<
               getUnsignedNBitInteger< 1 >( validityBits )<<std::endl;

    if( dataType == 11 || dataType == 12 || dataType == 13 )
    {
        parseDopplerOrbitData( fileBlock );
    }




}

int main( )
{
    std::ifstream dataFile( "/home/dominic/Downloads/mromagr2017_117_0745xmmmv1.odf", std::ios_base::binary);

    char currentFileBlock[ 9 ][ 4 ];

    readOdfFileBlock( currentFileBlock, dataFile );
    {
        int32_t primaryKey;
        uint32_t secondaryKey, logicalrecordLength, groupStartPacketNumber;
        parseHeader( currentFileBlock, primaryKey, secondaryKey, logicalrecordLength, groupStartPacketNumber );
        std::cout<<primaryKey<<" "<<secondaryKey<<" "<<logicalrecordLength<<" "<<groupStartPacketNumber <<std::endl;
    }

    readOdfFileBlock( currentFileBlock, dataFile );
    {
        std::string systemId, programId, fileCreationDate, fileCreationTime;
        uint32_t spacecraftIdNumber, fileReferenceDate, fileReferenceTime;
        parseFileLabel( currentFileBlock, systemId, programId, fileCreationDate, fileCreationTime,
                        spacecraftIdNumber, fileReferenceDate, fileReferenceTime );
        std::cout<<" "<<systemId<<" "<<programId<<" "<<fileCreationDate<<" "<<fileCreationTime<<" "<<
                                spacecraftIdNumber<<" "<<fileReferenceDate<<" "<<fileReferenceTime<<std::endl;



    }

    readOdfFileBlock( currentFileBlock, dataFile );
    {
        int32_t primaryKey;
        uint32_t secondaryKey, logicalrecordLength, groupStartPacketNumber;
        parseHeader( currentFileBlock, primaryKey, secondaryKey, logicalrecordLength, groupStartPacketNumber );
        std::cout<<primaryKey<<" "<<secondaryKey<<" "<<logicalrecordLength<<" "<<groupStartPacketNumber <<std::endl;
    }

    readOdfFileBlock( currentFileBlock, dataFile );
    parseIdentifierGroup( currentFileBlock );

    readOdfFileBlock( currentFileBlock, dataFile );
    {
        int32_t primaryKey;
        uint32_t secondaryKey, logicalrecordLength, groupStartPacketNumber;
        parseHeader( currentFileBlock, primaryKey, secondaryKey, logicalrecordLength, groupStartPacketNumber );
        std::cout<<primaryKey<<" "<<secondaryKey<<" "<<logicalrecordLength<<" "<<groupStartPacketNumber <<std::endl;
    }


    readOdfFileBlock( currentFileBlock, dataFile );
    {
        parseOrbitData( currentFileBlock );
    }


}

//int main( )
//{
//    std::ifstream dataFile( "/home/dominic/Downloads/mromagr2017_117_0745xmmmv1.odf", std::ios_base::binary);

//    int32_t headerType = read_s32_le( dataFile );
//    uint32_t secondaryKey = read_u32_le( dataFile );
//    uint32_t isNotEof = read_u32_le( dataFile );
//    uint32_t groupStartPacket = read_u32_le( dataFile );

//    std::cout<<headerType<<" "<<secondaryKey<<" "<<isNotEof<<" "<<groupStartPacket<<std::endl;

//    for( int i = 0; i < 20; i++ )
//    {
//       dataFile.get( );
//    }

//    char systemId[ 8 ], programId[ 8 ];

//    dataFile.read( systemId, 8 * sizeof( char ) );
//    dataFile.read( programId, 8 * sizeof( char ) );

//    uint32_t spacecraftId = read_u32_le( dataFile );
//    std::cout<<spacecraftId<<std::endl;

//    uint32_t fileCreationDate = read_u32_le( dataFile );
//    uint32_t fileCreationTime = read_u32_le( dataFile );

//    std::cout<<fileCreationDate<<" "<<fileCreationTime<<std::endl;

//    uint32_t referenceDate = read_u32_le( dataFile );
//    uint32_t referenceTime = read_u32_le( dataFile );

//    std::cout<<referenceDate<<" "<<referenceTime<<std::endl;

//     char headerGroup[ 36 ];
//     dataFile.read( headerGroup, 36 * sizeof( char ) );

//     for( int i = 0; i < 36; i++ )
//     {
//         std::cout<<int( headerGroup[ i ] )<<" ";
//     }

//     char identifierGroup[ 36 ];
//     dataFile.read( identifierGroup, 36 * sizeof( char ) );

//     std::cout<<std::endl;
//     for( int i = 0; i < 36; i++ )
//     {
//         std::cout<<( identifierGroup[ i ] );
//     }

//     std::cout<<std::endl;

//     char headerGroup2[ 36 ];
//     dataFile.read( headerGroup2, 36 * sizeof( char ) );

//     for( int i = 0; i < 36; i++ )
//     {
//         std::cout<<i<<" "<<int( headerGroup2[ i ] )<<" "<<std::endl;;
//     }

//     uint32_t recordTimeTagInt = read_u32_le( dataFile );
//     std::cout<<"Time tag: "<<recordTimeTagInt<<std::endl;

//     // Orbit Data Group, byte 4-7
//     uint32_t testBits = read_u32_le( dataFile );
//     std::cout<<testBits<<std::endl;
//     std::bitset< 32 > testBitSet( testBits );

//     std::bitset< 10 > fractionTimeTagBitset = getBitsetSegment< 10, 32 >(
//                 testBitSet, 0 );
//     std::bitset< 22 > downlinkDelayBitset = getBitsetSegment< 22, 32 >(
//                 testBitSet, 10 );

//     int fractionTimeTag = getSignedNBitInteger< 10 >( fractionTimeTagBitset );
//     int downlinkDelay = getSignedNBitInteger< 22 >( downlinkDelayBitset );

//     std::cout<<testBitSet<<std::endl;
//     std::cout<<"Time fraction: "<<fractionTimeTagBitset<<" "<<fractionTimeTag <<std::endl;
//     std::cout<<"Downlink delay: "<<downlinkDelayBitset<<" "<<downlinkDelay <<std::endl;

//     uint32_t observableIntegerPart = read_u32_le( dataFile );
//     uint32_t observableFractionalPart = read_u32_le( dataFile );
//     std::cout<<"Observable: "<<observableIntegerPart<<" "<<observableFractionalPart<<std::endl;

//     testBits = read_u32_le( dataFile );
//     testBitSet = std::bitset< 32 >( testBits );

//     std::cout<<testBitSet<<std::endl;
//     std::bitset< 3 > formatIdBitset = getBitsetSegment< 3, 32 >(
//                 testBitSet, 0 );
//     std::cout<<formatIdBitset<<" "<<formatIdBitset.to_ulong( )<<std::endl;


//}
