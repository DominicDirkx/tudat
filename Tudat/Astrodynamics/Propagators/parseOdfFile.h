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
                    startTimes )
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

class ProcessdOdfFileDopplerData: public ProcessdOdfFileSingleLinkData
{
public:

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

class ProcessdOdfFileSingleLinkData
{
public:

    ProcessdOdfFileSingleLinkData( ){ }

    virtual ~ProcessdOdfFileSingleLinkData( ){ }

    std::vector< double > observationTimes;
    std::vector< double > observableValues;
    std::vector< double > receiverDownlinkDelay;
    std::string transmittingStation;
    std::string receivingStation;
    int transmitterNetworkId;

    int downlinkBand;
    int uplinkBand;
    int referenceBand;


    std::map< double, double > getObservationData( )
    {
        return utilities::createVectorFromMapValues( observationTimes, observableValues );
    }
};

class ProcessedOdfFileContents
{
public:

    std::string spacecraftName;
    std::map< observation_models::ObservableType, std::vector< ProcessdOdfFileSingleLinkData > > dataBlocks;
    std::map< std::string, boost::shared_ptr< ProcessdOdfFileRampData > > rampDataBlocks;
};

boost::shared_ptr< ProcessedOdfFileContents > createOdfFileContents( )
{

}

} // namespace orbit_determination

} // namespace tudat

#endif // TUDAT_VPARSEODFFILE_H
