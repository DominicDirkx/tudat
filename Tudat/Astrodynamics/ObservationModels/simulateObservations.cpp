#include "Tudat/Astrodynamics/ObservationModels/simulateObservations.h"

namespace tudat
{

namespace observation_models
{

Eigen::VectorXd getIdenticallyAndIndependentlyDistributedNoise(
        const boost::function< double( const double ) > noiseFunction,
        const int observationSize,
        const double evaluationTime )
{
    Eigen::VectorXd noiseValues = Eigen::VectorXd( observationSize );
    for( int i = 0; i < observationSize; i++ )
    {
        noiseValues( i ) = noiseFunction( evaluationTime );
    }
    return noiseValues;
}

}

}