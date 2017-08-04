/*    Copyright (c) 2010-2017, Delft University of Technology
 *    All rigths reserved
 *
 *    This file is part of the Tudat. Redistribution and use in source and
 *    binary forms, with or without modification, are permitted exclusively
 *    under the terms of the Modified BSD license. You should have received
 *    a copy of the license with this file. If not, please or visit:
 *    http://tudat.tudelft.nl/LICENSE.
 *
 */

#include <cmath>
#include <iostream>
#include <iomanip>

#include <Eigen/LU>

#include "Tudat/Basics/utilities.h"
#include "Tudat/Mathematics/BasicMathematics/leastSquaresEstimation.h"

namespace tudat
{

namespace linear_algebra
{

//! Function to get condition number of matrix (using SVD decomposition)
double getConditionNumberOfInformationMatrix( const Eigen::MatrixXd informationMatrix )
{
    return getConditionNumberOfDecomposedMatrix(
                ( informationMatrix.jacobiSvd( Eigen::ComputeThinU | Eigen::ComputeFullV ) ) );
}

//! Function to get condition number of matrix from SVD decomposition
double getConditionNumberOfDecomposedMatrix( const Eigen::JacobiSVD< Eigen::MatrixXd >& singularValueDecomposition )
{
    Eigen::VectorXd singularValues = singularValueDecomposition.singularValues( );
    return singularValues( 0 ) / singularValues( singularValues.rows( ) - 1 );
}


//! Solve system of equations with SVD decomposition, checking condition number in the process
Eigen::VectorXd solveSystemOfEquationsWithSvd( const Eigen::MatrixXd matrixToInvert,
                                               const Eigen::VectorXd rightHandSideVector,
                                               const bool checkConditionNumber,
                                               const double maximumAllowedConditionNumber )
{
    Eigen::JacobiSVD< Eigen::MatrixXd > svdDecomposition = matrixToInvert.jacobiSvd(
                Eigen::ComputeThinU | Eigen::ComputeThinV );
    if( checkConditionNumber )
    {
        double conditionNumber = getConditionNumberOfDecomposedMatrix( svdDecomposition );

        if( conditionNumber > maximumAllowedConditionNumber )
        {
            std::cerr<<"Warning when performing least squares, condition number is "<<conditionNumber<<std::endl;
        }
    }
    return svdDecomposition.solve( rightHandSideVector );
}

//! Function to multiply information matrix by diagonal weights matrix
Eigen::MatrixXd multiplyInformationMatrixByDiagonalWeightMatrix(
        const Eigen::MatrixXd& informationMatrix,
        const Eigen::VectorXd& diagonalOfWeightMatrix )
{
//    std::cout<<"Pre-multiplication "<<informationMatrix.rows( )<<" "<<informationMatrix.cols( )<<std::endl;

    Eigen::MatrixXd weightedInformationMatrix = Eigen::MatrixXd::Zero( informationMatrix.rows( ), informationMatrix.cols( ) );

//    std::cout<<"Pre-multiplication"<<std::endl;
    for( int i = 0; i < informationMatrix.cols( ); i++ )
    {
        weightedInformationMatrix.block( 0, i, informationMatrix.rows( ), 1 ) =
                informationMatrix.block( 0, i, informationMatrix.rows( ), 1 ).cwiseProduct( diagonalOfWeightMatrix );
    }
//    std::cout<<"Post-multiplication"<<std::endl;

    return weightedInformationMatrix;
}

//! Function to compute inverse of covariance matrix at current iteration, including influence of a priori information
Eigen::MatrixXd calculateInverseOfUpdatedCovarianceMatrixWithoutSignatureFreeParameters(
        const Eigen::MatrixXd& informationMatrix,
        const Eigen::VectorXd& diagonalOfWeightMatrix,
        const Eigen::MatrixXd& inverseOfAPrioriCovarianceMatrix,
        std::vector< unsigned int >& removedIndices )
{

    removedIndices.clear( );
    for( unsigned int i = 0; i < informationMatrix.cols( ); i++ )
    {
        double columnMinimum = informationMatrix.block( 0, i, informationMatrix.rows( ), 1 ).minCoeff( );
        double columnMaximum = informationMatrix.block( 0, i, informationMatrix.rows( ), 1 ).maxCoeff( );

        if( std::fabs( columnMinimum ) < 2.0 * std::numeric_limits< double >::epsilon( ) &&
                std::fabs( columnMaximum ) < 2.0 * std::numeric_limits< double >::epsilon( ) )
        {
            removedIndices.push_back( i );
        }
    }


    Eigen::MatrixXd reducedInformationMatrix;
    reducedInformationMatrix.setZero( informationMatrix.rows( ), informationMatrix.cols( ) - removedIndices.size( ) );

    Eigen::MatrixXd inverseOfReducedAPrioriCovarianceMatrix;
    inverseOfReducedAPrioriCovarianceMatrix.setZero(
                informationMatrix.cols( ) - removedIndices.size( ), informationMatrix.cols( ) - removedIndices.size( ) );

    removedIndices.push_back( std::numeric_limits< int >::max( ) );

    int currentMatrixIndex = 0;
    int currentRemovalIndex = 0;



    for( unsigned int i = 0; i < informationMatrix.cols( ); i++ )
    {
        if( removedIndices.at( currentRemovalIndex ) == i )
        {
            currentRemovalIndex++;
        }
        else
        {
            reducedInformationMatrix.block( 0, currentMatrixIndex, informationMatrix.rows( ), 1 ) =
                    informationMatrix.block( 0, i, informationMatrix.rows( ), 1 );
            currentMatrixIndex++;
        }

//        int currentInnerMatrixIndex = 0;
//        int currentInnerRemovalIndex = 0;

//        for( unsigned int j = 0; j < informationMatrix.cols( ); j++ )
//        {
//            if( removedIndices.at( currentInnerRemovalIndex ) == j )
//            {
//                currentInnerRemovalIndex++;
//            }
//            else
//            {
//                std::cout<<i<<" "<<j<<" "<<currentMatrixIndex - 1<<" "<<currentInnerMatrixIndex<<std::endl;
//                inverseOfReducedAPrioriCovarianceMatrix(
//                            currentMatrixIndex - 1, currentInnerMatrixIndex ) = inverseOfAPrioriCovarianceMatrix( i, j );
//                currentInnerMatrixIndex++;
//            }
//        }
    }

//    std::cout<<"Post-it"<<std::endl;
//    std::cout<<"Original information: "<<std::endl<<std::endl<<informationMatrix<<std::endl<<std::endl;
//    std::cout<<"Reduced information: "<<std::endl<<std::endl<<reducedInformationMatrix<<std::endl<<std::endl;
//    std::cout<<"Weights: "<<std::endl<<std::endl<<diagonalOfWeightMatrix.transpose( )<<std::endl<<std::endl;

    return ( inverseOfReducedAPrioriCovarianceMatrix + reducedInformationMatrix.transpose( ) * multiplyInformationMatrixByDiagonalWeightMatrix(
                 reducedInformationMatrix, diagonalOfWeightMatrix ) );
}

//! Function to compute inverse of covariance matrix at current iteration, including influence of a priori information
Eigen::MatrixXd calculateUpdatedCovarianceMatrixWithoutSignatureFreeParameters(
        const Eigen::MatrixXd& informationMatrix,
        const Eigen::VectorXd& diagonalOfWeightMatrix,
        const Eigen::MatrixXd& inverseOfAPrioriCovarianceMatrix,
        const double valueForSignatureFreeEntries )
{
    std::vector< unsigned int > removedIndices;

//    std::cout<<"Test A "<<std::endl;

    Eigen::MatrixXd inverseOfReducedCovarianceMatrix = calculateInverseOfUpdatedCovarianceMatrixWithoutSignatureFreeParameters(
                informationMatrix, diagonalOfWeightMatrix, inverseOfAPrioriCovarianceMatrix, removedIndices );

//    std::cout<<"Indices removed: "<<removedIndices.size( ) - 1<<std::endl;

    if( inverseOfReducedCovarianceMatrix.rows( ) != inverseOfReducedCovarianceMatrix.cols( ) )
    {
        throw std::runtime_error( "Error when getting reduced covariance, inverseOfReducedCovarianceMatrix is not square" );
    }

    if( inverseOfReducedCovarianceMatrix.cols( ) + removedIndices.size( ) - 1 != informationMatrix.cols( ) )
    {
        throw std::runtime_error( "Error when getting reduced covariance, reduced matrix size is inconsistent" );

    }

//    std::cout<<"Reduced inv. cov."<<std::setprecision( 4 )<<std::endl<<inverseOfReducedCovarianceMatrix<<std::endl<<std::endl;

    Eigen::MatrixXd reducedCovarianceMatrix = inverseOfReducedCovarianceMatrix.inverse( );

//    std::cout<<"Reduced cov."<<std::endl<<reducedCovarianceMatrix<<std::endl<<std::endl;

//    std::cout<<"Weights."<<std::endl<<diagonalOfWeightMatrix.transpose( )<<std::endl<<std::endl;
//    std::cout<<"Partials ."<<std::endl<<informationMatrix<<std::endl<<std::endl;

    Eigen::MatrixXd fullCovarianceMatrix;
    fullCovarianceMatrix.setZero( informationMatrix.cols( ), informationMatrix.cols( ) );

    int removedIndex = 0;
    int reducedMatrixIndex = 0;

    for( unsigned int i = 0; i < fullCovarianceMatrix.rows( ); i++ )
    {
        if( removedIndices.at( removedIndex ) == i )
        {
            for( unsigned int j = 0; j < fullCovarianceMatrix.cols( ); j++ )
            {
                fullCovarianceMatrix( i, j ) = valueForSignatureFreeEntries;
            }
            removedIndex++;
        }
        else
        {
            int removedInnerIndex = 0;
            int reducedMatrixInnerIndex = 0;
            for( unsigned int j = 0; j < fullCovarianceMatrix.cols( ); j++ )
            {
                //std::cout<<"Index: "<<i<<" "<<j<<std::endl;
                if( removedIndices.at( removedInnerIndex ) == j )
                {
                    fullCovarianceMatrix( i, j ) = valueForSignatureFreeEntries;
                    removedInnerIndex++;;
                }
                else
                {
//                    std::cout<<"Inner Index: "<<i<<" "<<j<<" "<<reducedMatrixIndex<<" "<<reducedMatrixInnerIndex<<std::endl;
//                    std::cout<<"Inner Sizes: "<<fullCovarianceMatrix.rows( )<<" "<<fullCovarianceMatrix.cols( )<<" "<<
//                               reducedCovarianceMatrix.rows( )<<" "<<reducedCovarianceMatrix.cols( )<<std::endl;

                    fullCovarianceMatrix( i, j ) =  reducedCovarianceMatrix( reducedMatrixIndex, reducedMatrixInnerIndex );
//                    std::cout<<"Current covariance "<<fullCovarianceMatrix<<std::endl;

                    reducedMatrixInnerIndex++;
                }
            }
            reducedMatrixIndex++;
        }
    }

    return fullCovarianceMatrix;
}

//! Function to compute inverse of covariance matrix at current iteration, including influence of a priori information
Eigen::MatrixXd calculateInverseOfUpdatedCovarianceMatrix(
        const Eigen::MatrixXd& informationMatrix,
        const Eigen::VectorXd& diagonalOfWeightMatrix,
        const Eigen::MatrixXd& inverseOfAPrioriCovarianceMatrix )
{

    return inverseOfAPrioriCovarianceMatrix + informationMatrix.transpose( ) * multiplyInformationMatrixByDiagonalWeightMatrix(
                informationMatrix, diagonalOfWeightMatrix );
}

//! Function to compute inverse of covariance matrix at current iteration
Eigen::MatrixXd calculateInverseOfUpdatedCovarianceMatrix(
        const Eigen::MatrixXd& informationMatrix,
        const Eigen::VectorXd& diagonalOfWeightMatrix )
{
    return calculateInverseOfUpdatedCovarianceMatrix(
                informationMatrix, diagonalOfWeightMatrix,
                Eigen::MatrixXd::Zero( informationMatrix.cols( ), informationMatrix.cols( ) ) );
}

//! Function to perform an iteration least squares estimation from information matrix, weights and residuals and a priori
//! information
std::pair< Eigen::VectorXd, Eigen::MatrixXd > performLeastSquaresAdjustmentFromInformationMatrix(
        const Eigen::MatrixXd& informationMatrix,
        const Eigen::VectorXd& observationResiduals,
        const Eigen::VectorXd& diagonalOfWeightMatrix,
        const Eigen::MatrixXd& inverseOfAPrioriCovarianceMatrix,
        const bool checkConditionNumber,
        const double maximumAllowedConditionNumber )
{
    Eigen::VectorXd rightHandSide = informationMatrix.transpose( ) *
            ( diagonalOfWeightMatrix.cwiseProduct( observationResiduals ) );
    Eigen::MatrixXd inverseOfCovarianceMatrix = calculateInverseOfUpdatedCovarianceMatrix(
                informationMatrix, diagonalOfWeightMatrix, inverseOfAPrioriCovarianceMatrix );
    return std::make_pair( solveSystemOfEquationsWithSvd( inverseOfCovarianceMatrix, rightHandSide,
                                                          checkConditionNumber, maximumAllowedConditionNumber ),
                           inverseOfCovarianceMatrix );
}

//! Function to perform an iteration least squares estimation from information matrix, weights and residuals
std::pair< Eigen::VectorXd, Eigen::MatrixXd > performLeastSquaresAdjustmentFromInformationMatrix(
        const Eigen::MatrixXd& informationMatrix,
        const Eigen::VectorXd& observationResiduals,
        const Eigen::VectorXd& diagonalOfWeightMatrix,
        const bool checkConditionNumber,
        const double maximumAllowedConditionNumber )
{
    return performLeastSquaresAdjustmentFromInformationMatrix(
                informationMatrix, observationResiduals, diagonalOfWeightMatrix,
                Eigen::MatrixXd::Zero( informationMatrix.cols( ), informationMatrix.cols( ) ),
                checkConditionNumber, maximumAllowedConditionNumber );
}

//! Function to perform an iteration of least squares estimation from information matrix and residuals
std::pair< Eigen::VectorXd, Eigen::MatrixXd > performLeastSquaresAdjustmentFromInformationMatrix(
        const Eigen::MatrixXd& informationMatrix,
        const Eigen::VectorXd& observationResiduals,
        const bool checkConditionNumber,
        const double maximumAllowedConditionNumber )
{
    return performLeastSquaresAdjustmentFromInformationMatrix(
                informationMatrix, observationResiduals, Eigen::VectorXd::Constant( observationResiduals.size( ), 1, 1.0 ),
                checkConditionNumber, maximumAllowedConditionNumber );
}

//! Function to fit a univariate polynomial through a set of data
Eigen::VectorXd getLeastSquaresPolynomialFit(
        const Eigen::VectorXd& independentValues,
        const Eigen::VectorXd& dependentValues,
        const std::vector< double >& polynomialPowers )
{
    if( independentValues.rows( ) != dependentValues.rows( ) )
    {
        throw std::runtime_error( "Error when doing least squares polynomial fit, size of dependent and independent variable vectors is not equal" );
    }

    Eigen::MatrixXd informationMatrix = Eigen::MatrixXd::Zero( dependentValues.rows( ), polynomialPowers.size( ) );

    // Compute information matrix
    for( int i = 0; i < independentValues.rows( ); i++ )
    {
        for( unsigned int j = 0; j < polynomialPowers.size( ); j++ )
        {
            informationMatrix( i, j ) = std::pow( independentValues( i ), polynomialPowers.at( j ) );
        }
    }

    return performLeastSquaresAdjustmentFromInformationMatrix( informationMatrix, dependentValues ).first;
}

//! Function to fit a univariate polynomial through a set of data
std::vector< double > getLeastSquaresPolynomialFit(
        const std::map< double, double >& independentDependentValueMap,
        const std::vector< double >& polynomialPowers )
{
    return utilities::convertEigenVectorToStlVector(
                getLeastSquaresPolynomialFit(
                    utilities::convertStlVectorToEigenVector(
                        utilities::createVectorFromMapKeys( independentDependentValueMap ) ),
                    utilities::convertStlVectorToEigenVector(
                        utilities::createVectorFromMapValues( independentDependentValueMap ) ), polynomialPowers ) );

}

} // namespace linear_algebra

} // namespace tudat
