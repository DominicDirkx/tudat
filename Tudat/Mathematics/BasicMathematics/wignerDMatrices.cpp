#include <iostream>

#include <Tudat/Mathematics/BasicMathematics/wignerDMatrices.h>
#include <Tudat/Mathematics/BasicMathematics/mathematicalConstants.h>

namespace tudat
{

namespace basic_mathematics
{

//! Constructor
WignerDMatricesCache::WignerDMatricesCache( const int maximumDegree ):
    maximumDegree_( maximumDegree ), computeAngularMomentumOperators_( false )
{
    wignerDMatrices_.resize( maximumDegree + 1 );

    for( int l = 0; l <= maximumDegree; l++ )
    {
        wignerDMatrices_[ l ] = Eigen::MatrixXd::Zero( 2 * l + 1, 2 * l + 1 );
    }
    wignerDMatrices_[ 0 ]( 0, 0 ) = std::complex< double >( 1.0, 0.0 );

    polarToCartesianCoordinates_ <<
            std::complex< double >( std::sqrt( 1.0 / 2.0 ), 0.0 ), std::complex< double >( 0.0, 0.0 ), std::complex< double >( -std::sqrt( 1.0 / 2.0 ), 0.0 ),
            std::complex< double >( 0.0, 1.0 / std::sqrt( 2.0 ) ), std::complex< double >( 0.0, 0.0 ), std::complex< double >( 0.0, std::sqrt( 1.0 / 2.0 ) ),
            0.0, 1.0, 0.0;

    //    transformationMatrixToCartesianBasis_ << std::complex< double >( 1.0 / std::sqrt( 2.0 ), 0.0 ),
    //            std::complex< double >( 0.0, 0.0 ), std::complex< double >( -1.0 / std::sqrt( 2.0 ), 0.0 ),
    //            std::complex< double >( 0.0, 1.0 / std::sqrt( 2.0 ) ), std::complex< double >( 0.0, 0.0 ),
    //            std::complex< double >( 0.0, 1.0 / std::sqrt( 2.0 ) ),
    //            std::complex< double >( 0.0, 0.0 ), std::complex< double >( 1.0, 0.0 ), std::complex< double >( 0.0, 0.0 ) ;


    computeCoefficients( );
}

//! Function to update contents of this object to new orientation
void WignerDMatricesCache::updateMatrices( const std::complex< double > cayleyKleinA, const std::complex< double > cayleyKleinB )
{
    // Set current orientation
    currentCayleyKleinA_ = cayleyKleinA;
    currentCayleyKleinB_ = cayleyKleinB;
    currentCayleyKleinAConjugate_ = std::conj( currentCayleyKleinA_ );
    currentCayleyKleinBConjugate_ = std::conj( currentCayleyKleinB_ );

    // Explicitly compute coefficients at degree 1
    if( maximumDegree_ > 0 )
    {
        wignerDMatrices_[ 1 ]( 2, 2 ) = currentCayleyKleinA_ * currentCayleyKleinA_;
        wignerDMatrices_[ 1 ]( 2, 1 ) = -std::sqrt( 2.0 ) * currentCayleyKleinA_ * currentCayleyKleinBConjugate_;
        wignerDMatrices_[ 1 ]( 2, 0 ) = currentCayleyKleinBConjugate_ * currentCayleyKleinBConjugate_;
        wignerDMatrices_[ 1 ]( 1, 2 ) = std::sqrt( 2.0 ) * currentCayleyKleinA_ * currentCayleyKleinB_;
        wignerDMatrices_[ 1 ]( 1, 1 ) = std::norm( currentCayleyKleinA_ ) * std::norm( currentCayleyKleinA_ ) -
                std::norm( currentCayleyKleinB_ ) * std::norm( currentCayleyKleinB_ );
        wignerDMatrices_[ 1 ]( 1, 0 ) = -std::sqrt( 2.0 ) * currentCayleyKleinAConjugate_ * currentCayleyKleinBConjugate_;
        wignerDMatrices_[ 1 ]( 0, 2 ) = currentCayleyKleinB_ * currentCayleyKleinB_;
        wignerDMatrices_[ 1 ]( 0, 1 ) = std::sqrt( 2.0 ) * currentCayleyKleinAConjugate_ * currentCayleyKleinB_;
        wignerDMatrices_[ 1 ]( 0, 0 ) = currentCayleyKleinAConjugate_ * currentCayleyKleinAConjugate_;
    }

    // Recursively compute coefficients at degree >1
    for( int l = 2; l <= maximumDegree_; l++ )
    {
        for( int i = l; i <= 2 * l; i++ )
        {
            for( int j = 0; j <= 2 * l; j++ )
            {
                if( i - 2 >= 0 )
                {
                    wignerDMatrices_[ l ]( i, j ) = 0.0;

                    // For each part in equation, check if contribution is non-zer0
                    if( j > 1 )
                    {
                        wignerDMatrices_[ l ]( i, j ) += coefficientsIndexMinusOne_[ l ]( i, j ) * wignerDMatrices_[ 1 ]( 2, 2 ) *
                                wignerDMatrices_[ l - 1 ]( i - 2, j - 2 );
                    }
                    if( ( j > 0 ) && ( j <= 2 * l - 1 ) )
                    {
                        wignerDMatrices_[ l ]( i, j ) +=
                                coefficientsIndexZero_[ l ]( i, j ) * wignerDMatrices_[ 1 ]( 2, 1 ) *
                                wignerDMatrices_[ l - 1 ]( i - 2, j - 1 );
                    }
                    if( j < 2 * l - 1 )
                    {
                        wignerDMatrices_[ l ]( i, j ) += coefficientsIndexOne_[ l ]( i, j ) * wignerDMatrices_[ 1 ]( 2, 0 ) *
                                wignerDMatrices_[ l - 1 ]( i - 2, j );
                    }
                }
                else
                {
                    wignerDMatrices_[ l ]( i, j ) = std::complex< double >( 0.0, 0.0 );
                }
            }
        }

        // Use symmetry relation to compute coefficients for negative m
        int m = 0, k = 0;
        for( int i = 0; i < l; i++ )
        {
            m = i - l;
            for( int j = 0; j <= 2 * l; j++ )
            {
                k = j - l;
                wignerDMatrices_[ l ]( i, j ) = ( ( ( ( m - k ) % 2 ) == 0 ) ? 1.0 : -1.0 ) *
                        std::conj( wignerDMatrices_[ l ]( -m + l, -k + l ) );
            }

        }
    }

    //if( computeAngularMomentumOperators_ )
    {
        computeAngularMomentumOperators( );
    }
}

void WignerDMatricesCache::computeAngularMomentumOperators( )
{
    Eigen::Vector3cd currentWignerDMatrixVector, currentAngularMomentumOperatorInCartesianCoordinates_;
    int m, k;
    for( int l = 0; l <= maximumDegree_; l++ )
    {
//        std::cout<<"Wigner: "<<l<<std::endl<<wignerDMatrices_[ l ]<<std::endl;
        for( int i = 0; i <= 2 * l; i++ )
        {
            m = i - l;
            for( int j = 0; j <= 2 * l; j++ )
            {

                currentWignerDMatrixVector.setZero( );
                k = j - l;

//                std::cout<<"Ang mom. op. "<<l<<" "<<m<<" "<<k<<" "<<std::endl;

                if( ( i + 1 ) <= 2 * l )
                {
                    currentWignerDMatrixVector( 0 ) = wignerDMatrices_[ l ]( i + 1, j );
                }
                else
                {
                    currentWignerDMatrixVector( 0 ) = 0.0;
                }

                currentWignerDMatrixVector( 1 ) = wignerDMatrices_[ l ]( i, j );

                if( i > 0 )
                {
                    currentWignerDMatrixVector( 2 ) = wignerDMatrices_[ l ]( i - 1, j );
                }
                else
                {
                    currentWignerDMatrixVector( 2 ) = 0.0;
                }

                currentAngularMomentumOperatorInCartesianCoordinates_ =
                        polarToCartesianCoordinates_ * ( angularMomentumOperatorCoefficients_.at( l ).at( i ) *
                                                         currentWignerDMatrixVector );
                angularMomentumOperator_[ l ][ m ][ k ] = currentAngularMomentumOperatorInCartesianCoordinates_;

//                std::cout<<"Wigner D vector "<<currentWignerDMatrixVector.transpose( )<<std::endl<<
//                           "Ang. mom. coefficients: "<<std::endl<<angularMomentumOperatorCoefficients_.at( l ).at( i )<<std::endl<<
//                           "Ang. mom. op. polar: "<<std::endl<<( angularMomentumOperatorCoefficients_.at( l ).at( i ) *
//                                                                 currentWignerDMatrixVector )<<std::endl<<
//                           "Conv. mat. "<<std::endl<<polarToCartesianCoordinates_<<std::endl<<
//                           "Ang. mom. op. Cart.: "<<std::endl<<currentAngularMomentumOperatorInCartesianCoordinates_<<std::endl<<std::endl;
            }
        }
    }
}

//! Function to precompute the coefficients used on the recursive formulation for Wigner D-matrices
void WignerDMatricesCache::computeCoefficients( )
{
    coefficientsIndexMinusOne_.resize( maximumDegree_ + 1 );
    coefficientsIndexZero_.resize( maximumDegree_ + 1 );
    coefficientsIndexOne_.resize( maximumDegree_ + 1 );

    Eigen::Matrix3cd angularMomentumOperatorMultiplier;
    double angularMomentumScalingEntry0, angularMomentumScalingEntry2;
    int m = 0, k = 0;
    for( int l = 0; l <= maximumDegree_; l++ )
    {
        // Allocate size for coefficients of current degree
        coefficientsIndexMinusOne_[ l ] = Eigen::MatrixXd::Zero( 2 * l + 1, 2 * l + 1 );
        coefficientsIndexZero_[ l ] = Eigen::MatrixXd::Zero( 2 * l + 1, 2 * l + 1 );
        coefficientsIndexOne_[ l ] = Eigen::MatrixXd::Zero( 2 * l + 1, 2 * l + 1 );

        // Compute coefficients at current degree.
        for( int i = 0; i <= 2 * l; i++ )
        {
            m = i - l;

            angularMomentumScalingEntry0 = std::sqrt(
                        static_cast< double >( l *( l + 1 ) - m * ( m + 1 ) ) / 2.0 );
            angularMomentumScalingEntry2 = std::sqrt(
                        static_cast< double >( l *( l + 1 ) - m * ( m - 1 ) ) / 2.0 );
            angularMomentumOperatorMultiplier.setZero( );


            angularMomentumOperatorMultiplier( 0, 0 ) = -mathematical_constants::COMPLEX_I * angularMomentumScalingEntry0;
            angularMomentumOperatorMultiplier( 2, 2 ) = mathematical_constants::COMPLEX_I * angularMomentumScalingEntry2;
            angularMomentumOperatorMultiplier( 1, 1 ) = -static_cast< double >( m );

            angularMomentumOperatorCoefficients_[ l ][ i ] = angularMomentumOperatorMultiplier;

            for( int j = 0; j <= 2 * l; j++ )
            {
                k = j - l;
                coefficientsIndexMinusOne_[ l ]( i, j ) = std::sqrt(
                            static_cast< double >( ( l + k ) * ( l + k - 1 ) ) /
                            static_cast< double >( ( l + m ) * ( l + m - 1 ) ) );
                coefficientsIndexZero_[ l ]( i, j ) = std::sqrt(
                            static_cast< double >( 2 * ( l + k ) * ( l - k ) ) /
                            static_cast< double >( ( l + m ) * ( l + m - 1 ) ) );
                coefficientsIndexOne_[ l ]( i, j ) = std::sqrt(
                            static_cast< double >( ( l - k ) * ( l - k - 1 ) ) /
                            static_cast< double >( ( l + m ) * ( l + m - 1 ) ) );
            }
        }
    }

}

} // namespace basic_mathematics

} // namespace tudat
