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
        outputBits[ i ] = inputBits[ i + startIndex ];
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


void readOdfFileBlock(
        char fileBlock[ 9 ][ 4 ],
        std::istream& file )
{
    for( unsigned int i = 0; i < 9; i++ )
    {
        file.read( (char*)fileBlock[ i ], 4 );
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
