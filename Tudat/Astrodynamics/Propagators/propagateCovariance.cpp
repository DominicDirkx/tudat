
#include<Tudat/Astrodynamics/Propagators/propagateCovariance.h>

namespace tudat
{

namespace propagators
{

//! Function to retrueve full state transition and sensitivity matrices at epochs
void getFullVariationalEquationsSolutionHistory(
        std::map< double, Eigen::MatrixXd >& fullVariationalEquationsSolutionHistory,
        const std::shared_ptr< CombinedStateTransitionAndSensitivityMatrixInterface > stateTransitionInterface,
        const double timeStep,
        const double initialTime,
        const double finalTime )
{
    double currentTime = initialTime;
    while( currentTime < finalTime )
    {
        fullVariationalEquationsSolutionHistory[ currentTime ] =
                stateTransitionInterface->getFullCombinedStateTransitionAndSensitivityMatrix( currentTime );
        currentTime += timeStep;

    }
}

//! Function to propagate full covariance at the initial time to state covariance at later times
void propagateCovariance(
        std::map< double, Eigen::MatrixXd >& propagatedCovariance,
        const Eigen::MatrixXd& initialCovariance,
        const std::map< double, Eigen::MatrixXd >& fullVariationalEquationsSolutionHistory )
{
    for( auto resultIterator : fullVariationalEquationsSolutionHistory )
    {
        propagatedCovariance[ resultIterator.first ] =
                resultIterator.second * initialCovariance *
                resultIterator.second.transpose( );
    }
}

//! Function to propagate full covariance at the initial time to state covariance at later times
void propagateCovariance(
        std::map< double, Eigen::MatrixXd >& propagatedCovariance,
        const Eigen::MatrixXd& initialCovariance,
        const std::shared_ptr< CombinedStateTransitionAndSensitivityMatrixInterface > stateTransitionInterface,
        const double timeStep,
        const double initialTime,
        const double finalTime )
{
    if( initialCovariance.rows( ) != stateTransitionInterface->getFullParameterVectorSize( ) )
    {
        throw std::runtime_error( "Error when propagating single-arc covariance, sizes are incompatible" );
    }

    std::map< double, Eigen::MatrixXd > fullVariationalEquationsSolutionHistory;
    getFullVariationalEquationsSolutionHistory(
                fullVariationalEquationsSolutionHistory, stateTransitionInterface, timeStep, initialTime, finalTime );
    propagateCovariance( propagatedCovariance, initialCovariance, fullVariationalEquationsSolutionHistory );
}

//! Function to propagate full covariance at the initial time to state formal errors at later times
void propagateFormalErrors(
        std::map< double, Eigen::VectorXd >& propagatedFormalErrors,
        const Eigen::MatrixXd& initialCovariance,
        const std::shared_ptr< CombinedStateTransitionAndSensitivityMatrixInterface > stateTransitionInterface,
        const double timeStep,
        const double initialTime,
        const double finalTime )
{
    std::map< double, Eigen::MatrixXd > propagatedCovariance;
    propagateCovariance(
                propagatedCovariance, initialCovariance, stateTransitionInterface, timeStep, initialTime, finalTime );

    for( auto covarianceIterator : propagatedCovariance )
    {
        propagatedFormalErrors[ covarianceIterator.first ] =
                Eigen::VectorXd( covarianceIterator.second.diagonal( ).array( ).sqrt( ) );
    }
}

} // namespace propagators

} // namespace tudat

