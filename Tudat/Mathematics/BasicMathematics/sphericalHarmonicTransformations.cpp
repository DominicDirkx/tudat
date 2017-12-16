#include <boost/math/special_functions/binomial.hpp>
#include <iostream>
#include "Tudat/Mathematics/BasicMathematics/mathematicalConstants.h"
#include "Tudat/Mathematics/BasicMathematics/basicMathematicsFunctions.h"

#include "Tudat/Mathematics/BasicMathematics/sphericalHarmonicTransformations.h"
#include "Tudat/Mathematics/BasicMathematics/cayleyKleinParameters.h"

namespace tudat
{

namespace basic_mathematics
{

//! Function to update Wigner D-matrices for current orientation, parameterized by Cayley-Klein parameters
void SphericalHarmonicTransformationCache::updateFromCayleyKleinParameters(
        const std::complex< double > cayleyKleinA,
        const std::complex< double > cayleyKleinB )
{
    wignerDMatricesCache_->updateMatrices( cayleyKleinA, cayleyKleinB );
}

//! Function to update Wigner D-matrices for current orientation, parameterized by quaternion
void SphericalHarmonicTransformationCache::updateFromQuaternion(
        const Eigen::Quaterniond& rotationQuaternion )
{
    std::complex< double > cayleyKleinA;
    std::complex< double > cayleyKleinB;

    convertQuaterionToCayleyKleinParameters( rotationQuaternion, cayleyKleinA, cayleyKleinB );
    updateFromCayleyKleinParameters( cayleyKleinA, cayleyKleinB );
}

//! Function to transform spherical harmonic coefficients using current wignerDMatricesCache_
void SphericalHarmonicTransformationCache::transformCoefficientsAtDegree(
        const Eigen::MatrixXd& originalCosineCoefficients,
        const Eigen::MatrixXd& originalSineCoefficients,
        Eigen::MatrixXd& currentCosineCoefficients,
        Eigen::MatrixXd& currentSineCoefficients,
        const bool areCoefficientsNormalized,
        const bool computeAngularMomentumOperatorProduct )
{
    // Resize transformed coefficients
    currentCosineCoefficients.setZero( originalCosineCoefficients.rows( ), originalCosineCoefficients.cols( ) );
    currentSineCoefficients.setZero( originalSineCoefficients.rows( ), originalSineCoefficients.cols( ) );

    double currentMultiplier;
    double orderMMultiplier;

    // Iterate over all degrees
    for( unsigned int l = 0; l < originalCosineCoefficients.rows( ); l++ )
    {
        // Iterate over al, orders
        for( unsigned int m = 0; ( m <= l && m < originalCosineCoefficients.cols( ) ); m++ )
        {
            // Compute coefficient multipliers for (un-)normalized coefficients
            if( !areCoefficientsNormalized )
            {
                orderMMultiplier = std::sqrt( boost::math::factorial< double >( l - m ) /
                                              boost::math::factorial< double >( l + m ) );
                currentMultiplier = orderMMultiplier;
            }
            else
            {
                orderMMultiplier = ( m == 0 ? 1.0 : 1.0 / std::sqrt( 2.0 ) );
                currentMultiplier = orderMMultiplier;
            }

            // Transform zonal coefficient to current order
            std::complex< double > orderZeroDFunction = wignerDMatricesCache_->getWignerDCoefficient( l, m, 0 );
            currentCosineCoefficients( l, m ) +=
                    ( currentMultiplier ) * orderZeroDFunction.real( ) * originalCosineCoefficients( l, 0 );
            currentSineCoefficients( l, m ) +=
                    ( currentMultiplier ) * orderZeroDFunction.imag( ) * originalCosineCoefficients( l, 0 );

            if( computeAngularMomentumOperatorProduct )
            {
                currentAngularMomentumProduceCosineCoefficients_[ l ][ m ].setZero( );
                currentAngularMomentumProduceSineCoefficients_[ l ][ m ].setZero( );

                Eigen::Vector3cd orderZeroAngularMomentumOperator =
                        wignerDMatricesCache_->getAngularMomentumOperatorOnWignerDCoefficient( l, m, 0 );

                //std::cout<<" wigner d function ang. mom "<<l<<" "<<m<<" "<<orderZeroAngularMomentumOperator.transpose( )<<std::endl;

                currentAngularMomentumProduceCosineCoefficients_[ l ][ m ] +=
                        ( currentMultiplier ) * orderZeroAngularMomentumOperator.real( ) * originalCosineCoefficients( l, 0 );
                currentAngularMomentumProduceSineCoefficients_[ l ][ m ] +=
                        ( currentMultiplier ) * orderZeroAngularMomentumOperator.imag( ) * originalCosineCoefficients( l, 0 );

            }

            // Iterate over all original orders, and transform to new coefficients
            for( int k = 1; k <= static_cast< int >( l ); k++ )
            {
                // Compute muliplier
                if( !areCoefficientsNormalized )
                {
                    currentMultiplier = std::sqrt( boost::math::factorial< double >( l + k ) /
                                                   boost::math::factorial< double >( l - k ) ) * orderMMultiplier;
                }
                else
                {
                    currentMultiplier = std::sqrt( 2.0 ) * orderMMultiplier;
                }
                double signMultiplier = ( ( ( k % 2 ) == 0 ) ? ( 1.0 ) : ( -1.0 ) );

                // Retrieve Wigner D-matrices for plus/minus current order
                std::complex< double > orderKDFunction = wignerDMatricesCache_->getWignerDCoefficient( l, m, k );
                std::complex< double > orderMinusKDFunction = wignerDMatricesCache_->getWignerDCoefficient( l, m, -k );

                // Compute addition to current order cosine coefficient
                currentCosineCoefficients( l, m ) += 0.5 * currentMultiplier * (
                            ( signMultiplier * orderKDFunction.real( ) + orderMinusKDFunction.real( ) ) *
                            originalCosineCoefficients( l, k )  +
                            ( signMultiplier * orderKDFunction.imag( ) - orderMinusKDFunction.imag( ) ) *
                            originalSineCoefficients( l, k ) );

                // Compute addition to current order sine coefficient
                currentSineCoefficients( l, m ) += 0.5 * currentMultiplier * (
                            ( signMultiplier * orderKDFunction.imag( ) + orderMinusKDFunction.imag( ) ) *
                            originalCosineCoefficients( l, k )  +
                            ( -signMultiplier * orderKDFunction.real( ) + orderMinusKDFunction.real( ) ) *
                            originalSineCoefficients( l, k ) );

                if( computeAngularMomentumOperatorProduct )
                {
                    Eigen::Vector3cd orderKAngularMomentumOperator =
                            wignerDMatricesCache_->getAngularMomentumOperatorOnWignerDCoefficient( l, m, k );
                    Eigen::Vector3cd orderMinusKAngularMomentumOperator =
                            wignerDMatricesCache_->getAngularMomentumOperatorOnWignerDCoefficient( l, m, -k );

                    currentAngularMomentumProduceCosineCoefficients_[ l ][ m ] += 0.5 * currentMultiplier * (
                                ( signMultiplier * orderKAngularMomentumOperator.real( ) +
                                  orderMinusKAngularMomentumOperator.real( ) ) *
                                originalCosineCoefficients( l, k )  +
                                ( signMultiplier * orderKAngularMomentumOperator.imag( ) -
                                  orderMinusKAngularMomentumOperator.imag( ) ) *
                                originalSineCoefficients( l, k ) );
                    currentAngularMomentumProduceSineCoefficients_[ l ][ m ] += 0.5 * currentMultiplier * (
                                ( signMultiplier * orderKAngularMomentumOperator.imag( ) +
                                  orderMinusKAngularMomentumOperator.imag( ) ) *
                                originalCosineCoefficients( l, k )  +
                                ( -signMultiplier * orderKAngularMomentumOperator.real( ) +
                                  orderMinusKAngularMomentumOperator.real( ) ) *
                                originalSineCoefficients( l, k ) );


                }
            }

            // Compute final scaling
            currentCosineCoefficients( l, m ) = currentCosineCoefficients( l, m ) *
                    ( ( ( m % 2 ) == 0 ) ? ( 1.0 ) : ( -1.0 ) );
            currentSineCoefficients( l, m ) = currentSineCoefficients( l, m ) *
                    ( ( ( ( m + 1 ) % 2 ) == 0 ) ? ( 1.0 ) : ( -1.0 ) );
            if( m > 0 )
            {
                currentCosineCoefficients( l, m ) = currentCosineCoefficients( l, m ) * 2.0;
                currentSineCoefficients( l, m ) = currentSineCoefficients( l, m ) * 2.0;
            }

            if( computeAngularMomentumOperatorProduct )
            {
                currentAngularMomentumProduceCosineCoefficients_[ l ][ m ] =
                        currentAngularMomentumProduceCosineCoefficients_[ l ][ m ] *
                        ( ( ( m % 2 ) == 0 ) ? ( 1.0 ) : ( -1.0 ) );
                currentAngularMomentumProduceSineCoefficients_[ l ][ m ] =
                        currentAngularMomentumProduceSineCoefficients_[ l ][ m ] *
                        ( ( ( ( m + 1 ) % 2 ) == 0 ) ? ( 1.0 ) : ( -1.0 ) );
                if( m > 0 )
                {
                    currentAngularMomentumProduceCosineCoefficients_[ l ][ m ] =
                            currentAngularMomentumProduceCosineCoefficients_[ l ][ m ] * 2.0;
                    currentAngularMomentumProduceSineCoefficients_[ l ][ m ] =
                            currentAngularMomentumProduceSineCoefficients_[ l ][ m ] * 2.0;
                }
            }

            //std::cout<<"Ang. mom operator of single coefficient: "<<l<<" "<<m<<std::endl<<
            //           currentAngularMomentumProduceCosineCoefficients_[ l ][ m ].transpose( )<<" "<<
            //           currentAngularMomentumProduceSineCoefficients_[ l ][ m ].transpose( )<<std::endl;
        }
    }
    //std::cout<<std::endl<<std::endl<<std::endl;
}

}

}


