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

observation_models::ObservableType getObservableTypeForOdfId(
        const int odfId )
{
    observation_models::ObservableType observableType;

    switch( odfId )
    {
    case 11:
        observableType = observation_models::one_way_differenced_range;
        break;
    case 12:
        observableType = observation_models::n_way_differenced_range;
        break;
    case 13:
        observableType = observation_models::n_way_differenced_range;
        break;
    default:
        throw std::runtime_error( "Error wen getting observable type for ODF ID, ID: " +
                                  std::to_string( odfId ) + " not recognized." );
    }
    return observableType;
}


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

    std::vector< std::string > originFile;

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

class RampedReferenceFrequencyInterpolator
{
public:
    RampedReferenceFrequencyInterpolator(
            std::vector< input_output::OdfRampBlock > rampBlock )
    {
        for( unsigned int i = 0; i < rampBlock.size( ); i++ )
        {
            startTimes.push_back( rampBlock.at( i ).getRampStartTime( ) );
            endTimes.push_back( rampBlock.at( i ).getRampEndTime( ) );
            rampRates.push_back( rampBlock.at( i ).getRampRate( ) );
            startFrequency.push_back( rampBlock.at( i ).getRampStartFrequency( ) );
        }

        startTimeLookupScheme_ = boost::make_shared<
                interpolators::HuntingAlgorithmLookupScheme< double > >(
                    startTimes );
    }

    RampedReferenceFrequencyInterpolator(
            const std::vector< double >& startTimes_,
            const std::vector< double >& endTimes_,
            const std::vector< double >& rampRates_,
            const std::vector< double >& startFrequency_ ):
        startTimes( startTimes_ ), endTimes( endTimes_ ), rampRates( rampRates_ ), startFrequency( startFrequency_ )
    {
        startTimeLookupScheme_ = boost::make_shared<
                interpolators::HuntingAlgorithmLookupScheme< double > >(
                    startTimes );
    }


    double getCurrentReferenceFrequencyIntegral(
            const double lookupTime, const double integrationTime,
            bool& issueForInterpolation )
    {
        int lowerNearestNeighbour = startTimeLookupScheme_->findNearestLowerNeighbour(
                    lookupTime );

        if( lookupTime > endTimes.at( lowerNearestNeighbour ) )
        {
            issueForInterpolation = true;
        }

        return integrationTime * (
                    startFrequency.at( lowerNearestNeighbour ) +
                    rampRates.at( lowerNearestNeighbour ) * ( lookupTime - startTimes.at( lowerNearestNeighbour ) ) +
                    0.5 * rampRates.at( lowerNearestNeighbour ) * integrationTime );
    }

    double getCurrentReferenceFrequency(
            const double lookupTime,
            bool& issueForInterpolation )
    {
        int lowerNearestNeighbour = startTimeLookupScheme_->findNearestLowerNeighbour(
                    lookupTime );


        if( lookupTime > endTimes.at( lowerNearestNeighbour ) )
        {
            issueForInterpolation = true;
        }

        return  startFrequency.at( lowerNearestNeighbour ) +
                rampRates.at( lowerNearestNeighbour ) * ( lookupTime - startTimes.at( lowerNearestNeighbour ) );
    }


    std::vector< double > startTimes;
    std::vector< double > endTimes;
    std::vector< double > rampRates;
    std::vector< double > startFrequency;

    boost::shared_ptr< interpolators::LookUpScheme< double > > startTimeLookupScheme_;

};



class ProcessedOdfFileContents
{
public:

    std::string spacecraftName;

    std::map< observation_models::ObservableType, std::map< std::pair< std::string, std::string >,
    boost::shared_ptr< ProcessdOdfFileSingleLinkData > > > processedDataBlocks;

    std::map< int, boost::shared_ptr< RampedReferenceFrequencyInterpolator > > rampInterpolators;
};


boost::shared_ptr< RampedReferenceFrequencyInterpolator > mergeRampDataInterpolators(
        const std::vector< boost::shared_ptr< RampedReferenceFrequencyInterpolator > >& interpolatorList )
{
    std::map< double, double > rampEndTimesPerStation;
    std::map< double, double > rampRatesPerStation;
    std::map< double, double > rampStartFrequenciesPerStation;

    for( unsigned int i = 0; i < interpolatorList.size( ); i++ )
    {
        for( unsigned int j = 0; j < interpolatorList.at( i )->startTimes.size( ); j++ )
        {
            rampEndTimesPerStation[ interpolatorList.at( i )->startTimes.at( j ) ] =
                    interpolatorList.at( i )->endTimes.at( j );
            rampRatesPerStation[ interpolatorList.at( i )->startTimes.at( j ) ] =
                    interpolatorList.at( i )->rampRates.at( j );
            rampStartFrequenciesPerStation[ interpolatorList.at( i )->startTimes.at( j ) ] =
                    interpolatorList.at( i )->startFrequency.at( j );
        }
    }
    return boost::make_shared< RampedReferenceFrequencyInterpolator >(
                utilities::createVectorFromMapKeys( rampEndTimesPerStation ),
                utilities::createVectorFromMapValues( rampEndTimesPerStation ),
                utilities::createVectorFromMapValues( rampRatesPerStation ),
                utilities::createVectorFromMapValues( rampStartFrequenciesPerStation ) );
}

boost::shared_ptr< ProcessedOdfFileContents > mergeOdfFileContents(
        const std::vector< boost::shared_ptr< ProcessedOdfFileContents > > odfFileContents )
{
    std::map< observation_models::ObservableType, std::map< std::pair< std::string, std::string >,
            boost::shared_ptr< ProcessdOdfFileSingleLinkData > > > mergedDataBlocks;
     std::map< int, std::vector< boost::shared_ptr< RampedReferenceFrequencyInterpolator > > > rampInterpolatorList;

    // Iterate over all ODF files
    for( unsigned int i = 0; i < odfFileContents.size( ); i++ )
    {
        // Retrieve contents of current file.
        std::map< observation_models::ObservableType, std::map< std::pair< std::string, std::string >,
                boost::shared_ptr< ProcessdOdfFileSingleLinkData > > >  dataBlocks =
                odfFileContents.at( i )->processedDataBlocks;

        for( auto it = dataBlocks.begin( ); it != dataBlocks.end( ); it++ )
        {
            for( auto linkIt = it->second.begin( ); linkIt != it->second.end( ); linkIt++ )
            {
                bool setNewObject = false;
                if( mergedDataBlocks.count( it->first ) == 0 )
                {
                    setNewObject = true;
                }
                else if( mergedDataBlocks.at( it->first ).count( linkIt->first ) == 0 )
                {
                    setNewObject = true;
                }

                boost::shared_ptr< ProcessdOdfFileSingleLinkData > blockToAdd;

                if( setNewObject )
                {
                    mergedDataBlocks[ it->first ][ linkIt->first ] = linkIt->second;
                }
                else
                {
                    blockToAdd = linkIt->second;

                    mergedDataBlocks[ it->first ][ linkIt->first ]->downlinkBand.insert(
                                mergedDataBlocks[ it->first ][ linkIt->first ]->downlinkBand.end( ),
                            blockToAdd->downlinkBand.begin( ), blockToAdd->downlinkBand.end( ) );

                    mergedDataBlocks[ it->first ][ linkIt->first ]->observableValues.insert(
                                mergedDataBlocks[ it->first ][ linkIt->first ]->observableValues.end( ),
                            blockToAdd->observableValues.begin( ), blockToAdd->observableValues.end( ) );

                    mergedDataBlocks[ it->first ][ linkIt->first ]->observationTimes.insert(
                                mergedDataBlocks[ it->first ][ linkIt->first ]->observationTimes.end( ),
                            blockToAdd->observationTimes.begin( ), blockToAdd->observationTimes.end( ) );

                    mergedDataBlocks[ it->first ][ linkIt->first ]->originFile.insert(
                                mergedDataBlocks[ it->first ][ linkIt->first ]->originFile.end( ),
                            blockToAdd->originFile.begin( ), blockToAdd->originFile.end( ) );

                    if( it->first == observation_models::one_way_differenced_range ||
                            it->first == observation_models::n_way_differenced_range )
                    {
                        boost::shared_ptr< ProcessdOdfFileDopplerData > dopplerBlockToAdd
                                = boost::dynamic_pointer_cast< ProcessdOdfFileDopplerData >(
                                    blockToAdd );
                        boost::shared_ptr< ProcessdOdfFileDopplerData > currentDopplerObservableMergedData =
                                boost::dynamic_pointer_cast< ProcessdOdfFileDopplerData >(
                                    mergedDataBlocks[ it->first ][ linkIt->first ] );


                        currentDopplerObservableMergedData->referenceFrequency.insert(
                                    currentDopplerObservableMergedData->referenceFrequency.end( ),
                                    dopplerBlockToAdd->referenceFrequency.begin( ), dopplerBlockToAdd->referenceFrequency.end( ) );

                        currentDopplerObservableMergedData->rampingFlag.insert(
                                    currentDopplerObservableMergedData->rampingFlag.end( ),
                                    dopplerBlockToAdd->rampingFlag.begin( ), dopplerBlockToAdd->rampingFlag.end( ) );

                        currentDopplerObservableMergedData->compressionTimes.insert(
                                    currentDopplerObservableMergedData->compressionTimes.end( ),
                                    dopplerBlockToAdd->compressionTimes.begin( ), dopplerBlockToAdd->compressionTimes.end( ) );

                        currentDopplerObservableMergedData->receiverChannels.insert(
                                    currentDopplerObservableMergedData->receiverChannels.end( ),
                                    dopplerBlockToAdd->receiverChannels.begin( ), dopplerBlockToAdd->receiverChannels.end( ) );
                    }
                }
            }
        }

        std::map< int, boost::shared_ptr< RampedReferenceFrequencyInterpolator > > currentRampInterpolators =
                odfFileContents.at( i )->rampInterpolators;

        for( auto it = currentRampInterpolators.begin( ); it != currentRampInterpolators.end( ); it++ )
        {
            rampInterpolatorList[ it->first ].push_back( it->second );
        }
    }

    std::map< int, boost::shared_ptr< RampedReferenceFrequencyInterpolator > > mergedRampInterpolators;
    for( auto it = rampInterpolatorList.begin( ); it != rampInterpolatorList.end( ); it++ )
    {
        mergedRampInterpolators[ it->first ] = mergeRampDataInterpolators(
                    it->second );
    }


    boost::shared_ptr< ProcessedOdfFileContents > processedOdfFile =
            boost::make_shared< ProcessedOdfFileContents >( );

    processedOdfFile->processedDataBlocks = mergedDataBlocks;
    processedOdfFile->rampInterpolators = mergedRampInterpolators;
    processedOdfFile->spacecraftName = "AAA";
    return processedOdfFile;
}


boost::shared_ptr< ProcessedOdfFileContents > parseOdfFileContents(
        const boost::shared_ptr< input_output::OdfRawFileContents > rawOdfData )
{
    // Create output object
    boost::shared_ptr< ProcessedOdfFileContents > processedOdfFile =
            boost::make_shared< ProcessedOdfFileContents >( );
    std::string spacecraftName = std::to_string( rawOdfData->spacecraftId );

    // Retrieve data blocks from ODF file raw contents
    std::vector< boost::shared_ptr< input_output::OdfDataBlock > > rawDataBlocks = rawOdfData->dataBlocks;

    // Create list of data, sorted by observable type and link ends; single object per combination of the two
    std::map< observation_models::ObservableType, std::map< std::pair< std::string, std::string >,
            boost::shared_ptr< ProcessdOdfFileSingleLinkData > > > processedDataBlocks;

    bool createNewObject = false;
    int currentObservableId;
    observation_models::ObservableType currentObservableType;
    std::pair< std::string, std::string > stationIds;

    // Iterate over all block of ODF file.
    for( unsigned int i = 0; i < rawDataBlocks.size( ); i++ )
    {
        // Retrieve observable type and link end names
        currentObservableId = rawDataBlocks.at( i )->observableSpecificDataBlock->dataType;
        currentObservableType = getObservableTypeForOdfId( currentObservableId );
        int appendedTransmittingStationId =
                rawDataBlocks.at( i )->commonDataBlock->transmittingStation + 100 *
                rawDataBlocks.at( i )->commonDataBlock->transmittingStationNetworkId;
        stationIds = std::make_pair( std::to_string( appendedTransmittingStationId ),
                                     std::to_string( rawDataBlocks.at( i )->commonDataBlock->receivingStation ) );

        // Check if data object already exists for current observable/link ends
        createNewObject = false;
        if( processedDataBlocks.count( currentObservableType ) == 0 )
        {
            createNewObject = true;
        }
        else if( processedDataBlocks.at( currentObservableType ).count( stationIds ) == 0 )
        {
            createNewObject = true;
        }

        // Create new data object, if required
        if( createNewObject )
        {
            if( currentObservableType == observation_models::one_way_differenced_range ||
                    currentObservableType == observation_models::n_way_differenced_range )
            {
                processedDataBlocks[ currentObservableType ][ stationIds ] = boost::make_shared< ProcessdOdfFileDopplerData >( );
                processedDataBlocks[ currentObservableType ][ stationIds ]->transmittingStation = appendedTransmittingStationId;
                processedDataBlocks[ currentObservableType ][ stationIds ]->receivingStation =
                        std::to_string( rawDataBlocks.at( i )->commonDataBlock->receivingStation );
                processedDataBlocks[ currentObservableType ][ stationIds ]->observableType = currentObservableType;
            }
            else
            {
                throw std::runtime_error( "Error when parsing ODF file contents, can currently only handle Doppler data" );
            }
        }

        // Add common properties to data object
        processedDataBlocks[ currentObservableType ][ stationIds ]->downlinkBand.push_back( rawDataBlocks.at( i )->commonDataBlock->downlinkBand );
        processedDataBlocks[ currentObservableType ][ stationIds ]->uplinkBand.push_back( rawDataBlocks.at( i )->commonDataBlock->uplinkBand );
        processedDataBlocks[ currentObservableType ][ stationIds ]->referenceBand.push_back( rawDataBlocks.at( i )->commonDataBlock->referenceBand );
        processedDataBlocks[ currentObservableType ][ stationIds ]->observableValues.push_back( rawDataBlocks.at( i )->commonDataBlock->getObservableValue( ) );
        processedDataBlocks[ currentObservableType ][ stationIds ]->observationTimes.push_back( rawDataBlocks.at( i )->commonDataBlock->getObservableTime( ) );
        processedDataBlocks[ currentObservableType ][ stationIds ]->receiverDownlinkDelay.push_back( rawDataBlocks.at( i )->commonDataBlock->receivingStationDownlinkDelay );

        // Add properties to data object for Doppler data
        if( currentObservableType == observation_models::one_way_differenced_range ||
                currentObservableType == observation_models::n_way_differenced_range )
        {
            boost::shared_ptr< input_output::OdfDopplerDataBlock > odfDopplerDataBlock =
                    boost::dynamic_pointer_cast< input_output::OdfDopplerDataBlock >(
                        rawDataBlocks.at( i )->observableSpecificDataBlock );
            boost::shared_ptr< ProcessdOdfFileDopplerData > odfParsedDopplerDataBlock =
                    boost::dynamic_pointer_cast< ProcessdOdfFileDopplerData >(
                        processedDataBlocks[ currentObservableType ][ stationIds ] );

            odfParsedDopplerDataBlock->compressionTimes.push_back( odfDopplerDataBlock->compressionTime );
            odfParsedDopplerDataBlock->receiverChannels.push_back( odfDopplerDataBlock->receiverChannel );
            odfParsedDopplerDataBlock->rampingFlag.push_back( odfDopplerDataBlock->receiverExciterFlag );
            odfParsedDopplerDataBlock->referenceFrequency.push_back( odfDopplerDataBlock->getReferenceFrequency( ) );
            odfParsedDopplerDataBlock->reservedData.push_back( odfDopplerDataBlock->reservedSegment );
            odfParsedDopplerDataBlock->uplinkDelays.push_back( odfDopplerDataBlock->transmittingStationDelay );
        }
    }

    // Save output and return
    processedOdfFile->processedDataBlocks = processedDataBlocks;
    processedOdfFile->spacecraftName = spacecraftName;

    std::map< int, std::vector< input_output::OdfRampBlock > > rampDataBlocks = rawOdfData->odfRampBlocks;
    std::map< int, boost::shared_ptr< RampedReferenceFrequencyInterpolator > > rampInterpolators;

    for( auto it = rampDataBlocks.begin( ); it != rampDataBlocks.end( ); it++ )
    {
        rampInterpolators[ it->first ] =
                boost::make_shared< RampedReferenceFrequencyInterpolator >( it->second );
    }
    processedOdfFile->rampInterpolators = rampInterpolators;

    return processedOdfFile;
}

} // namespace orbit_determination

} // namespace tudat

#endif // TUDAT_PARSEODFFILE_H
