#ifndef TUDAT_PARSEODFFILE_H
#define TUDAT_PARSEODFFILE_H

#include "Tudat/Basics/utilities.h"
#include "Tudat/InputOutput/readOdfFile.h"
#include "Tudat/Astrodynamics/ObservationModels/observableTypes.h"
#include "Tudat/Mathematics/Interpolators/lookupScheme.h"

namespace tudat
{

namespace orbit_determination
{

class ProcessdOdfFileRampData
{
public:

    ProcessdOdfFileRampData(
            const std::vector< double >& rampStartTime,
            const std::vector< double >& rampEndTime,
            const std::vector< double >& rampRate,
            const std::vector< double >& rampStartFrequency ):
        rampStartTime_( rampStartTime ), rampEndTime_( rampEndTime ), rampRate_( rampRate ),
        rampStartFrequency_( rampStartFrequency )
    {
        startTimeLookupScheme_ = boost::make_shared<
                interpolators::HuntingAlgorithmLookupScheme< double > >(
                    rampStartTime );
    }

    std::vector< double > rampStartTime_;
    std::vector< double > rampEndTime_;

    std::vector< double > rampRate_;
    std::vector< double > rampStartFrequency_;

    double getCurrentReferenceFrequencyIntegral(
            const double lookupTime, const double integrationTime )
    {

        int lowerNearestNeighbour = startTimeLookupScheme_->findNearestLowerNeighbour(
                    lookupTime );

        if( lookupTime > rampEndTime_.at( lowerNearestNeighbour ) )
        {
            std::cerr<<"Error when getting ramp frequency integral, time not inside arc range, requested "<<lookupTime<<
                       ", arc end time difference is "<<( lookupTime > rampEndTime_.at( lowerNearestNeighbour ) )<<std::endl;
        }

        return integrationTime * (
                    rampStartFrequency_.at( lowerNearestNeighbour ) +
                    rampRate_.at( lowerNearestNeighbour ) * ( lookupTime - rampStartTime_.at( lowerNearestNeighbour ) ) +
                    0.5 * rampRate_.at( lowerNearestNeighbour ) * integrationTime );
    }

    double getCurrentReferenceFrequency(
            const double lookupTime )
    {
        int lowerNearestNeighbour = startTimeLookupScheme_->findNearestLowerNeighbour(
                    lookupTime );
        if( lookupTime > rampEndTime_.at( lowerNearestNeighbour ) )
        {
            std::cerr<<"Error when getting ramp frequency, time not inside arc range, requested "<<lookupTime<<
                       ", arc end time difference is "<<( lookupTime > rampEndTime_.at( lowerNearestNeighbour ) )<<std::endl;
        }

        return  rampStartFrequency_.at( lowerNearestNeighbour ) +
                rampRate_.at( lowerNearestNeighbour ) * ( lookupTime - rampStartTime_.at( lowerNearestNeighbour ) );
    }

private:
    boost::shared_ptr< interpolators::LookUpScheme< double > > startTimeLookupScheme_;

};


class ProcessdOdfFileSingleLinkData
{
public:

    ProcessdOdfFileSingleLinkData( ){ }

    virtual ~ProcessdOdfFileSingleLinkData( ){ }

    std::vector< double > observationTimes;
    std::vector< double > observableValues;
    std::vector< double > receiverDownlinkDelay;

    std::vector< int > downlinkBand;
    std::vector< int > uplinkBand;
    std::vector< int > referenceBand;

    observation_models::ObservableType observableType;

    std::string transmittingStation;
    std::string receivingStation;
    int transmitterNetworkId;

    std::map< double, double > getObservationData( )
    {
        return utilities::createVectorFromMapValues( observationTimes, observableValues );
    }
};

class ProcessdOdfFileDopplerData: public ProcessdOdfFileSingleLinkData
{
public:

    ~ProcessdOdfFileDopplerData( ){ }

    std::vector< int > receiverChannels;
    std::vector< double > referenceFrequency;
    std::vector< double > compressionTimes;
    std::vector< double > uplinkDelays;
    std::vector< double > reservedData;
    std::vector< bool > rampingFlag;

    std::map< double, bool > getRampFlags( )
    {
        return utilities::createVectorFromMapValues( observationTimes, rampingFlag );
    }

    std::map< double, double > getReferenceFrequencies( )
    {
        return utilities::createVectorFromMapValues( observationTimes, referenceFrequency );
    }

    std::map< double, double > getCompressionTimes( )
    {
        return utilities::createVectorFromMapValues( observationTimes, compressionTimes );
    }
};


class ProcessedOdfFileContents
{
public:

    std::string spacecraftName;
    std::map< observation_models::ObservableType, std::vector< boost::shared_ptr< ProcessdOdfFileSingleLinkData > > > dataBlocks;

    std::map< std::pair< std::string, std::string >, boost::shared_ptr< ProcessdOdfFileDopplerData > > sortedDopplerDataBlocks;

    std::map< std::string, boost::shared_ptr< ProcessdOdfFileRampData > > rampDataBlocks;
};

boost::shared_ptr< ProcessedOdfFileContents > mergeOdfFileContents(
        const std::vector< boost::shared_ptr< ProcessedOdfFileContents > > odfFileContents )
{
   std::map< std::pair< std::string, std::string >, boost::shared_ptr< ProcessdOdfFileDopplerData > > sortedDopplerDataBlocks;

   std::cout<<"ODF files "<<odfFileContents.size( )<<std::endl;

   for( unsigned int i = 0; i < odfFileContents.size( ); i++ )
   {
       std::map< std::pair< std::string, std::string >, boost::shared_ptr< ProcessdOdfFileDopplerData > >
               singleFileSortedDopplerDataBlocks = odfFileContents.at( i )->sortedDopplerDataBlocks;

       for( auto it = singleFileSortedDopplerDataBlocks.begin( ); it != singleFileSortedDopplerDataBlocks.end( ); it++ )
       {
           if( sortedDopplerDataBlocks.count( it->first ) == 0 )
           {
               sortedDopplerDataBlocks[ it->first ] = it->second;
           }
           else
           {
               boost::shared_ptr< ProcessdOdfFileDopplerData > blockToAdd = it->second;

               sortedDopplerDataBlocks[ it->first ]->compressionTimes.insert(
                           sortedDopplerDataBlocks[ it->first ]->compressionTimes.end( ),
                      blockToAdd->compressionTimes.begin( ), blockToAdd->compressionTimes.end( ) );

               sortedDopplerDataBlocks[ it->first ]->receiverChannels.insert(
                           sortedDopplerDataBlocks[ it->first ]->receiverChannels.end( ),
                      blockToAdd->receiverChannels.begin( ), blockToAdd->receiverChannels.end( ) );

               sortedDopplerDataBlocks[ it->first ]->downlinkBand.insert(
                           sortedDopplerDataBlocks[ it->first ]->downlinkBand.end( ),
                      blockToAdd->downlinkBand.begin( ), blockToAdd->downlinkBand.end( ) );

               sortedDopplerDataBlocks[ it->first ]->observableValues.insert(
                           sortedDopplerDataBlocks[ it->first ]->observableValues.end( ),
                      blockToAdd->observableValues.begin( ), blockToAdd->observableValues.end( ) );

               sortedDopplerDataBlocks[ it->first ]->observationTimes.insert(
                           sortedDopplerDataBlocks[ it->first ]->observationTimes.end( ),
                      blockToAdd->observationTimes.begin( ), blockToAdd->observationTimes.end( ) );

               sortedDopplerDataBlocks[ it->first ]->referenceFrequency.insert(
                           sortedDopplerDataBlocks[ it->first ]->referenceFrequency.end( ),
                      blockToAdd->referenceFrequency.begin( ), blockToAdd->referenceFrequency.end( ) );

               sortedDopplerDataBlocks[ it->first ]->rampingFlag.insert(
                           sortedDopplerDataBlocks[ it->first ]->rampingFlag.end( ),
                      blockToAdd->rampingFlag.begin( ), blockToAdd->rampingFlag.end( ) );
           }

       }

   }

   std::cout<<"Sorting final "<<sortedDopplerDataBlocks.size( )<<std::endl;
   std::map< observation_models::ObservableType, std::vector< boost::shared_ptr< ProcessdOdfFileSingleLinkData > > > processedDataBlocks;
   for( auto it = sortedDopplerDataBlocks.begin( ); it != sortedDopplerDataBlocks.end( ); it++ )
   {
       processedDataBlocks[ it->second->observableType ].push_back( it->second );
       std::cout<<it->first.first<<" "<<it->first.second<<" "<<it->second->observableValues.size( )<<std::endl;
   }

   std::cout<<"Sorted"<<std::endl;
   boost::shared_ptr< ProcessedOdfFileContents > processedOdfFile =
           boost::make_shared< ProcessedOdfFileContents >( );

   processedOdfFile->dataBlocks = processedDataBlocks;
   processedOdfFile->spacecraftName = "AAA";
   return processedOdfFile;
}

boost::shared_ptr< ProcessedOdfFileContents > parseOdfFileContents(
        const boost::shared_ptr< input_output::OdfRawFileContents > rawOdfData )
{
    boost::shared_ptr< ProcessedOdfFileContents > processedOdfFile =
            boost::make_shared< ProcessedOdfFileContents >( );
    std::string spacecraftName = std::to_string( rawOdfData->spacecraftId );

    std::vector< boost::shared_ptr< input_output::OdfDataBlock > > rawDataBlocks = rawOdfData->dataBlocks;

    std::map< std::pair< std::string, std::string >, boost::shared_ptr< ProcessdOdfFileDopplerData > >
            processedDopplerDataBlocks;

    bool createNewObject = false;

    int currentObservableId;
    std::pair< std::string, std::string > stationIds;
    for( unsigned int i = 0; i < rawDataBlocks.size( ); i++ )
    {
        currentObservableId = rawDataBlocks.at( i )->observableSpecificDataBlock->dataType;
        stationIds = std::make_pair( std::to_string( rawDataBlocks.at( i )->commonDataBlock->transmittingStation ),
                                     std::to_string( rawDataBlocks.at( i )->commonDataBlock->receivingStation ) );

        if( rawDataBlocks.at( i )->commonDataBlock->transmittingStationNetworkId != 0 )
        {
            std::cerr<<"Warning when parsing ODF file, network ID is not DSN"<<std::endl;
        }

        createNewObject = false;

        if( rawDataBlocks.at( i )->observableSpecificDataBlock->dataType == 11 ||
                rawDataBlocks.at( i )->observableSpecificDataBlock->dataType == 12 ||
                rawDataBlocks.at( i )->observableSpecificDataBlock->dataType == 13 )
        {
            if( processedDopplerDataBlocks.count( stationIds ) == 0 )
            {
                createNewObject = true;
            }

            if( createNewObject )
            {
                observation_models::ObservableType currentObservable;
                if( currentObservableId == 11 )
                {
                    currentObservable = observation_models::one_way_differenced_range;
                }
                else
                {
                    currentObservable = observation_models::two_way_differenced_range;
                }
                processedDopplerDataBlocks[ stationIds ] = boost::make_shared< ProcessdOdfFileDopplerData >( );
                processedDopplerDataBlocks[ stationIds ]->transmittingStation =
                        std::to_string( rawDataBlocks.at( i )->commonDataBlock->transmittingStation );
                processedDopplerDataBlocks[ stationIds ]->receivingStation =
                        std::to_string( rawDataBlocks.at( i )->commonDataBlock->receivingStation );
                processedDopplerDataBlocks[ stationIds ]->observableType = currentObservable;
            }

            boost::shared_ptr< input_output::OdfDopplerDataBlock > odfDopplerDataBlock =
                    boost::dynamic_pointer_cast< input_output::OdfDopplerDataBlock >(
                        rawDataBlocks.at( i )->observableSpecificDataBlock );

            processedDopplerDataBlocks[ stationIds ]->downlinkBand.push_back( rawDataBlocks.at( i )->commonDataBlock->downlinkBand );
            processedDopplerDataBlocks[ stationIds ]->uplinkBand.push_back( rawDataBlocks.at( i )->commonDataBlock->uplinkBand );
            processedDopplerDataBlocks[ stationIds ]->referenceBand.push_back( rawDataBlocks.at( i )->commonDataBlock->referenceBand );
            processedDopplerDataBlocks[ stationIds ]->observableValues.push_back( rawDataBlocks.at( i )->commonDataBlock->getObservableValue( ) );
            processedDopplerDataBlocks[ stationIds ]->observationTimes.push_back( rawDataBlocks.at( i )->commonDataBlock->getObservableTime( ) );
            processedDopplerDataBlocks[ stationIds ]->receiverDownlinkDelay.push_back( rawDataBlocks.at( i )->commonDataBlock->receivingStationDownlinkDelay );

            processedDopplerDataBlocks[ stationIds ]->compressionTimes.push_back( odfDopplerDataBlock->compressionTime );
            processedDopplerDataBlocks[ stationIds ]->receiverChannels.push_back( odfDopplerDataBlock->receiverChannel );
            processedDopplerDataBlocks[ stationIds ]->rampingFlag.push_back( odfDopplerDataBlock->receiverExciterFlag );
            processedDopplerDataBlocks[ stationIds ]->referenceFrequency.push_back( odfDopplerDataBlock->getReferenceFrequency( ) );
            processedDopplerDataBlocks[ stationIds ]->reservedData.push_back( odfDopplerDataBlock->reservedSegment );
            processedDopplerDataBlocks[ stationIds ]->uplinkDelays.push_back( odfDopplerDataBlock->transmittingStationDelay );
        }
    }

    std::map< observation_models::ObservableType, std::vector< boost::shared_ptr< ProcessdOdfFileSingleLinkData > > > processedDataBlocks;
    for( auto it = processedDopplerDataBlocks.begin( ); it != processedDopplerDataBlocks.end( ); it++ )
    {
        processedDataBlocks[ it->second->observableType ].push_back( it->second );
        std::cout<<it->first.first<<" "<<it->first.second<<" "<<it->second->observableValues.size( )<<std::endl;
    }

    processedOdfFile->sortedDopplerDataBlocks = processedDopplerDataBlocks;
    processedOdfFile->dataBlocks = processedDataBlocks;
    processedOdfFile->spacecraftName = spacecraftName;
    return processedOdfFile;
}

} // namespace orbit_determination

} // namespace tudat

#endif // TUDAT_PARSEODFFILE_H
