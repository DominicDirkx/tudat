/*    Copyright (c) 2010-2017, Delft University of Technology
 *    All rigths reserved
 *
 *    This file is part of the Tudat. Redistribution and use in source and
 *    binary forms, with or without modification, are permitted exclusively
 *    under the terms of the Modified BSD license. You should have received
 *    a copy of the license with this file. If not, please or visit:
 *    http://tudat.tudelft.nl/LICENSE.
 */

#ifndef TUDAT_MUTUALFORCEPOTENTIAL_H
#define TUDAT_MUTUALFORCEPOTENTIAL_H

#include <boost/bind.hpp>
#include <boost/make_shared.hpp>
#include <boost/function.hpp>
#include <boost/math/special_functions/factorials.hpp>

#include <Eigen/Core>
#include <Eigen/Geometry>

#include "Tudat/Basics/utilities.h"
#include "Tudat/Astrodynamics/BasicAstrodynamics/physicalConstants.h"
#include "Tudat/Mathematics/BasicMathematics/coordinateConversions.h"
#include "Tudat/Mathematics/BasicMathematics/legendrePolynomials.h"
#include "Tudat/Mathematics/BasicMathematics/sphericalHarmonics.h"
#include "Tudat/Mathematics/BasicMathematics/sphericalHarmonicTransformations.h"

namespace tudat
{

namespace gravitation
{

//! Function to get maximum degrees of used for the spherical harmonic expansions of the two bodies
/*!
 *  Function to get maximum degrees of used for the spherical harmonic expansions of the two bodies
 *  \param coefficientCombinationsToUse st of degrees/orders that are to be used for the series expansion.
 *  Each tuple contains: (degree of body 1, order of body 1, degree of body 2, order of body 2)
 *  \return Maximum degree of body 1 and body 2 (as a pair)
 */
std::pair< int, int > getMaximumDegrees(
        const std::vector< boost::tuple< unsigned int, unsigned int, unsigned int, unsigned int > >& coefficientCombinationsToUse );

//! Function to compute cross-body normalization terms for mutual two-body potential
/*!
 * Function to compute cross-body normalization term gamma for mutual two-body potential, according to formulation of Compere &
 * Lemaitre (2014)
 * \param l Parameter l used by Compere & Lemaitre (2014); degree of body 1
 * \param m Parameter m used by Compere & Lemaitre (2014); order of body 1
 * \param j Parameter j used by Compere & Lemaitre (2014); degree of body 2
 * \param k Parameter k used by Compere & Lemaitre (2014); order of body 2
 * \return Term gamma for mutual two-body potential, according to formulation of Compere & Lemaitre (2014)
 */
double getGammaCoefficientForMutualForcePotential(
        const int l, const int m, const int j, const int k );

//! Function to compute cross-body normalization terms for mutual two-body potential, for unnormalized or fully normalized
//! coefficients
/*!
 * Function to compute cross-body normalization terms for mutual two-body potential, for unnormalized or fully normalized
 * coefficients
 * \param degree1 Degree of body 1
 * \param order1 Order of body 1
 * \param degree2 Degree of body 2
 * \param order2 Order of body 2
 * \param areCoefficientsNormalized Boolean denoting whether the coefficients are fully normalized or unnormalized
 * \return Cross-body normalization terms for mutual two-body potential
 */
double getMutualPotentialEffectiveCoefficientMultiplier(
        const int degree1, const int order1, const int degree2, const int order2, const bool areCoefficientsNormalized );

//! Function to compute single-term in two-body potential, from effective one-body formulation of Dirkx et al. (2018)
/*!
 * Function to compute single-term in two-body potential, from effective one-body formulation of Dirkx et al. (2018),
 * omitting the radius power term, and the common multiplier for all terms
 * \param effectiveCosineCoefficient Effective one-body cosine coefficient
 * \param effectiveSineCoefficient Effective one-body sine coefficient
 * \param sphericalHarmonicsCache Cache object used to pre-compute Legendre polynomials and other spherical harmonic terms
 * \return Single-term in effective one-body formulation, omitting the radius power term, and the common multiplier for all terms
 */
template< typename CoefficientType = double >
CoefficientType computeSingleMutualForcePotentialTerm(
        const CoefficientType effectiveCosineCoefficient,
        const CoefficientType effectiveSineCoefficient,
        boost::shared_ptr< basic_mathematics::SphericalHarmonicsCache > sphericalHarmonicsCache,
        const int degreeOfBody1,
        const int orderOfBody1,
        const int degreeOfBody2,
        const int orderOfBody2 )
{
    return ( effectiveCosineCoefficient * sphericalHarmonicsCache->getCosineOfMultipleLongitude(
                 std::abs( orderOfBody1 + orderOfBody2 ) ) +
             effectiveSineCoefficient * sphericalHarmonicsCache->getSineOfMultipleLongitude(
                 std::abs( orderOfBody1 + orderOfBody2 ) ) ) *
            sphericalHarmonicsCache->getLegendreCache( )->getLegendrePolynomial(
                degreeOfBody1 + degreeOfBody2, std::abs( orderOfBody1 + orderOfBody2 ) );
}

template< typename CoefficientType = double >
CoefficientType computeMutualForcePotential(
        const Eigen::Vector3d& bodyFixedPosition,
        const double effectiveGravitationalParameterOfBody1,
        const double equatorialRadiusOfBody1,
        const double equatorialRadiusOfBody2,
        const int maximumDegreeOfBody1,
        const int maximumDegreeOfBody2,
        const boost::function< CoefficientType( int, int, int, int ) >& effectiveCosineCoefficientFunction,
        const boost::function< CoefficientType( int, int, int, int ) >& effectiveSineCoefficientFunction,
        const std::vector< boost::tuple< unsigned int, unsigned int, unsigned int, unsigned int > >& coefficientCombinationsToUse,
        boost::shared_ptr< basic_mathematics::SphericalHarmonicsCache > sphericalHarmonicsCache )
{

    // Determine body fixed spherical position of body udnergoing acceleration.
    Eigen::Vector3d sphericalPositon =
            coordinate_conversions::convertCartesianToSpherical( bodyFixedPosition );
    double radius = sphericalPositon.x( );
    double latitude = mathematical_constants::PI / 2.0 - sphericalPositon.y( );
    double longitude = sphericalPositon.z( );

    double sineOfLatitude = std::sin( latitude );
    sphericalHarmonicsCache->update( TUDAT_NAN, sineOfLatitude, longitude, TUDAT_NAN );

    CoefficientType potential = utilities::getZeroEntry< CoefficientType >( );

    int degreeOfBody1, degreeOfBody2, orderOfBody1, orderOfBody2;

    std::vector< double > radiusRatioOfBody1List;
    double radiusRatioOfBody1 = equatorialRadiusOfBody1 / radius;
    radiusRatioOfBody1List.push_back( 1 );
    for( int i = 1; i <= maximumDegreeOfBody1; i++ )
    {
        radiusRatioOfBody1List.push_back( radiusRatioOfBody1List.at( i - 1 ) * radiusRatioOfBody1 );
    }

    std::vector< double > radiusRatioOfBody2List;
    radiusRatioOfBody2List.push_back( 1 );
    double radiusRatioOfBody2 = equatorialRadiusOfBody2 / radius;
    for( int i = 1; i <= maximumDegreeOfBody2; i++ )
    {
        radiusRatioOfBody2List.push_back( radiusRatioOfBody2List.at( i - 1 ) * radiusRatioOfBody2 );
    }


    CoefficientType currentTerm = utilities::getZeroEntry< CoefficientType >( );
    for(  unsigned int i = 0; i < coefficientCombinationsToUse.size( ); i++ )
    {
        std::cout<<i<<std::endl;
        degreeOfBody1 = coefficientCombinationsToUse.at( i ).get< 0 >( );
        orderOfBody1 = coefficientCombinationsToUse.at( i ).get< 1 >( );
        degreeOfBody2 = coefficientCombinationsToUse.at( i ).get< 2 >( );
        orderOfBody2 = coefficientCombinationsToUse.at( i ).get< 3 >( );


        currentTerm = utilities::getZeroEntry< CoefficientType >( );
        currentTerm += computeSingleMutualForcePotentialTerm(
                    effectiveCosineCoefficientFunction( degreeOfBody1, orderOfBody1, degreeOfBody2, orderOfBody2 ),
                    effectiveSineCoefficientFunction( degreeOfBody1, orderOfBody1, degreeOfBody2, orderOfBody2 ),
                    sphericalHarmonicsCache, degreeOfBody1, orderOfBody1, degreeOfBody2, orderOfBody2 );
        std::cout<<"Current D/O: "<<degreeOfBody1<<" "<<orderOfBody1<<" "<<degreeOfBody2<<" "<<orderOfBody2<<" "<<std::endl;
        std::cout<<"Current term: "<<std::endl<<effectiveCosineCoefficientFunction( degreeOfBody1, orderOfBody1, degreeOfBody2, orderOfBody2 )<<
                   std::endl<<std::endl<<
                effectiveSineCoefficientFunction( degreeOfBody1, orderOfBody1, degreeOfBody2, orderOfBody2 )<<std::endl;

        if( orderOfBody1 != 0 )
        {
            currentTerm += computeSingleMutualForcePotentialTerm(
                        effectiveCosineCoefficientFunction( degreeOfBody1, -orderOfBody1, degreeOfBody2, orderOfBody2 ),
                        effectiveSineCoefficientFunction( degreeOfBody1, -orderOfBody1, degreeOfBody2, orderOfBody2 ),
                        sphericalHarmonicsCache, degreeOfBody1, -orderOfBody1, degreeOfBody2, orderOfBody2 );

            std::cout<<"Current term A: "<<std::endl<<effectiveCosineCoefficientFunction( degreeOfBody1, orderOfBody1, degreeOfBody2, orderOfBody2 )<<
                       std::endl<<std::endl<<
                    effectiveSineCoefficientFunction( degreeOfBody1, orderOfBody1, degreeOfBody2, orderOfBody2 )<<std::endl;
        }

        if( orderOfBody2 != 0 )
        {
            currentTerm += computeSingleMutualForcePotentialTerm(
                        effectiveCosineCoefficientFunction( degreeOfBody1, orderOfBody1, degreeOfBody2, -orderOfBody2 ),
                        effectiveSineCoefficientFunction( degreeOfBody1, orderOfBody1, degreeOfBody2, -orderOfBody2 ),
                        sphericalHarmonicsCache, degreeOfBody1, orderOfBody1, degreeOfBody2, -orderOfBody2 );
        }

        if( ( orderOfBody1 != 0 ) && ( orderOfBody2 != 0 ) )
        {
            currentTerm += computeSingleMutualForcePotentialTerm(
                        effectiveCosineCoefficientFunction( degreeOfBody1, -orderOfBody1, degreeOfBody2, -orderOfBody2 ),
                        effectiveSineCoefficientFunction( degreeOfBody1, -orderOfBody1, degreeOfBody2, -orderOfBody2 ),
                        sphericalHarmonicsCache, degreeOfBody1, -orderOfBody1, degreeOfBody2, -orderOfBody2 );
        }

        currentTerm *= radiusRatioOfBody1List.at( degreeOfBody1 );
        currentTerm *= radiusRatioOfBody2List.at( degreeOfBody2 );

        potential += currentTerm;
         std::cout<<"Current potential: "<<std::endl<<potential<<std::endl<<std::endl<<std::endl;

    }

    // Multiply by central term and return
    return potential * effectiveGravitationalParameterOfBody1 / radius;
}

Eigen::Vector3d computeGeodesyNormalizedMutualGravitationalAccelerationSum(
        const Eigen::Vector3d& positionOfBodySubjectToAcceleration,
        const double gravitationalParameterOfBody,
        const double equatorialRadiusOfBody1,
        const double equatorialRadiusOfBody2,
        const boost::function< double( int, int, int, int ) >& effectiveCosineCoefficientFunction,
        const boost::function< double( int, int, int, int ) >& effectiveSineCoefficientFunction,
        const std::vector< boost::tuple< unsigned int, unsigned int, unsigned int, unsigned int > >& coefficientCombinationsToUse,
        const int maximumEvaluationDegree,
        const std::vector< double > radius1Powers,
        const std::vector< double > radius2Powers,
        boost::shared_ptr< basic_mathematics::SphericalHarmonicsCache > sphericalHarmonicsCache,
        std::map< int, std::map< int, std::map< int, std::map< int, Eigen::Vector3d > > > >& accelerationPerTerm,
        const bool saveTermsSeparately = 0,
        const Eigen::Matrix3d& accelerationRotation = Eigen::Matrix3d::Zero( ) );

Eigen::Vector3d computeUnnormalizedMutualGravitationalAccelerationSum(
        const Eigen::Vector3d& positionOfBodySubjectToAcceleration,
        const double gravitationalParameterOfBody,
        const double equatorialRadiusOfBody1,
        const double equatorialRadiusOfBody2,
        const boost::function< double( int, int, int, int ) >& effectiveCosineCoefficientFunction,
        const boost::function< double( int, int, int, int ) >& effectiveSineCoefficientFunction,
        const std::vector< boost::tuple< unsigned int, unsigned int, unsigned int, unsigned int > >& coefficientCombinationsToUse,
        boost::shared_ptr< basic_mathematics::SphericalHarmonicsCache > sphericalHarmonicsCache );

//void computePartialDerivativesOfPotentialComponentsWrtFullCoefficients(
//        std::vector< Eigen::Matrix< double, 1, 2 > >& potentialComponentsWrtFullCoefficients,
//        const std::vector< boost::tuple< unsigned int, unsigned int, unsigned int, unsigned int > >& coefficientCombinationsToUse,
//        const double distance,
//        const std::vector< double > radius1Powers,
//        const std::vector< double > radius2Powers,
//        boost::shared_ptr< basic_mathematics::SphericalHarmonicsCache > sphericalHarmonicsCache,
//        const boost::function< int( const int, const int, const int, const int )> effectiveIndexFunction );

inline double getSigmaSignFunction( const int order )
{
    return ( ( order >= 0 ) ? ( 1.0 ) : ( ( std::abs( order ) % 2 == 0 ) ? ( 1.0 ) : ( -1.0 ) ) );
}

//! Class to compute the effective one-body coefficients used to evaluate the two-body mutual potential, according to the
//! approach of Dirkx et al. (2018)
class EffectiveMutualSphericalHarmonicsField
{
public:

    //! Constructor
    /*!
     * Constructor
     * \param coefficientCombinationsToUse List of combinations of degrees/orders of the two bodies that are to be used, given in
     * order: (degree body 1, order body 1, degree body 2, order body 2)
     * \param cosineCoefficientFunctionOfBody1 Function returning the cosine coefficients of body 1 in its body-fixed frame
     * \param sineCoefficientFunctionOfBody1 Function returning the sine coefficients of body 1 in its body-fixed frame
     * \param cosineCoefficientFunctionOfBody2 Function returning the cosine coefficients of body 2 in its body-fixed frame
     * \param sineCoefficientFunctionOfBody2 Function returning the sine coefficients of body 2 in its body-fixed frame
     * \param gravitationalParameterFunction Function returning the effective gravitational parameter
     * \param equatorialRadiusOfBody1 Reference (equatorial) radius of body 1
     * \param equatorialRadiusOfBody2 Reference (equatorial) radius of body 2
     * \param areCoefficientsNormalized Boolean denoting whether the coefficients are (4 pi) normalized
     */
    EffectiveMutualSphericalHarmonicsField(
            const std::vector< boost::tuple< unsigned int, unsigned int, unsigned int, unsigned int > >&
            coefficientCombinationsToUse,
            const boost::function< Eigen::MatrixXd( ) > cosineCoefficientFunctionOfBody1,
            const boost::function< Eigen::MatrixXd( ) > sineCoefficientFunctionOfBody1,
            const boost::function< Eigen::MatrixXd( ) > cosineCoefficientFunctionOfBody2,
            const boost::function< Eigen::MatrixXd( ) > sineCoefficientFunctionOfBody2,
            const boost::function< double( ) > gravitationalParameterFunction,
            const double equatorialRadiusOfBody1,
            const double equatorialRadiusOfBody2,
            const bool areCoefficientsNormalized = 1 ):
        coefficientCombinationsToUse_( coefficientCombinationsToUse ),
        cosineCoefficientFunctionOfBody1_( cosineCoefficientFunctionOfBody1 ),
        sineCoefficientFunctionOfBody1_( sineCoefficientFunctionOfBody1 ),
        cosineCoefficientFunctionOfBody2_( cosineCoefficientFunctionOfBody2 ),
        sineCoefficientFunctionOfBody2_( sineCoefficientFunctionOfBody2 ),
        gravitationalParameterFunction_( gravitationalParameterFunction ),
        equatorialRadiusOfBody1_( equatorialRadiusOfBody1 ),
        equatorialRadiusOfBody2_( equatorialRadiusOfBody2 ),
        areCoefficientsNormalized_( areCoefficientsNormalized ),
        calculatePartials_( false )
    {
        std::pair< int, int > maximumDegrees = getMaximumDegrees( coefficientCombinationsToUse_ );
        maximumDegree1_ = maximumDegrees.first;
        maximumDegree2_ = maximumDegrees.second;

        effectiveCosineCoefficients_.resize(
                    getTotalVectorSize( ) );
        effectiveSineCoefficients_.resize(
                    getTotalVectorSize( ) );
        multipliers_.resize(
                    getTotalVectorSize( ) );

        effectiveAngularMomentumOperatorOfCosineCoefficients_.resize( getTotalVectorSize( ) );
        effectiveAngularMomentumOperatorOfSineCoefficients_.resize( getTotalVectorSize( ) );

        transformationCache_ = boost::make_shared< basic_mathematics::SphericalHarmonicTransformationCache >(
                    cosineCoefficientFunctionOfBody1( ).rows( ) + cosineCoefficientFunctionOfBody2( ).rows( ),
                    cosineCoefficientFunctionOfBody1( ).cols( ) + cosineCoefficientFunctionOfBody2( ).cols( ) );
        initializeMultipliers( );
    }

    //! Function to update the object to the current state, with the rotation provided as a quaternion
    /*!
     * Function to update the object to the current state, with the rotation provided as a quaternion. This function
     * updates the member variables, and ultimately the effectiveCosineCoefficients_ and effectiveSineCoefficients_
     * variables containing the effective one-body coefficients
     * \param coefficientRotationQuaterion Current rotation from body-fixed frame of body 2 to body-fixed frame of body 1
     */
    void computeCurrentEffectiveCoefficients(
            const Eigen::Quaterniond coefficientRotationQuaterion );

    //! Function to update the object to the current state, with the transformed coefficient provided manually
    /*!
     * Function to update the object to the current state, with the transformed coefficient provided manually. This function
     * updates the member variables, and ultimately the effectiveCosineCoefficients_ and effectiveSineCoefficients_
     * variables containing the effective one-body coefficients. As opposed to the computeCurrentEffectiveCoefficients function,
     * the gravity field coefficients of body 2, in frame fixed to body 1, are provided directly as input to this function.
     * \param transformedCosineCoefficients Current cosine gravity field coefficients of body 2, in frame fixed to body 1
     * \param transformedSineCoefficients Current sine gravity field coefficients of body 2, in frame fixed to body 1
     */
    void computeCurrentEffectiveCoefficientsFromManualTransformedCoefficients(
            const Eigen::MatrixXd& transformedCosineCoefficients,
            const Eigen::MatrixXd& transformedSineCoefficients );

    double getGravitationalPotential(
            const Eigen::Vector3d& bodyFixedPosition,
            boost::shared_ptr< basic_mathematics::SphericalHarmonicsCache > sphericalHarmonicsCache );

    Eigen::Vector3d getAngularMomentumOpertorOfGravitationalPotential(
            const Eigen::Vector3d& bodyFixedPosition,
            boost::shared_ptr< basic_mathematics::SphericalHarmonicsCache > sphericalHarmonicsCache );

    int getMaximumDegree1( )
    {
        return maximumDegree1_;
    }

    int getMaximumDegree2( )
    {
        return maximumDegree2_;
    }

    int getTotalVectorSize( )
    {
        return ( 2 * maximumDegree1_ + 1 ) * ( maximumDegree1_ + 1 ) * ( 2 * maximumDegree2_ + 1 ) * ( maximumDegree2_ + 1 );
    }

    int getEffectiveIndex(
            const int degreeOfBody1, const int orderOfBody1, const int degreeOfBody2, const int orderOfBody2 )
    {
        return degreeOfBody1 + ( maximumDegree1_ + 1 ) * ( orderOfBody1 + maximumDegree1_ +( 2 * maximumDegree1_ + 1 ) *
                                                ( degreeOfBody2 + ( maximumDegree2_ + 1 ) * ( orderOfBody2 + maximumDegree2_ ) ) );
    }

    double getEffectiveCosineCoefficient(
            const int degreeOfBody1, const int orderOfBody1, const int degreeOfBody2, const int orderOfBody2 )
    {
        return effectiveCosineCoefficients_[ getEffectiveIndex( degreeOfBody1, orderOfBody1, degreeOfBody2, orderOfBody2 ) ];
    }

    double getEffectiveSineCoefficient(
            const int degreeOfBody1, const int orderOfBody1, const int degreeOfBody2, const int orderOfBody2 )
    {
        return effectiveSineCoefficients_[ getEffectiveIndex( degreeOfBody1, orderOfBody1, degreeOfBody2, orderOfBody2 ) ];
    }

    Eigen::Vector3d getEffectiveAngularMomentumOperatorOfCosineCoefficient(
            const int degreeOfBody1, const int orderOfBody1, const int degreeOfBody2, const int orderOfBody2 )
    {
        return effectiveAngularMomentumOperatorOfCosineCoefficients_[ getEffectiveIndex( degreeOfBody1, orderOfBody1, degreeOfBody2, orderOfBody2 ) ];
    }

    Eigen::Vector3d getEffectiveAngularMomentumOperatorOfSineCoefficient(
            const int degreeOfBody1, const int orderOfBody1, const int degreeOfBody2, const int orderOfBody2 )
    {
        return effectiveAngularMomentumOperatorOfSineCoefficients_[ getEffectiveIndex( degreeOfBody1, orderOfBody1, degreeOfBody2, orderOfBody2 ) ];
    }

    boost::shared_ptr< basic_mathematics::SphericalHarmonicTransformationCache > getTransformationCache( )
    {
        return transformationCache_;
    }

//    void setCalculatePartials( )
//    {
//        transformedCosineCoefficientsOfBody2Partials_.resize( 3 );
//        transformedSineCoefficientsOfBody2Partials_.resize( 3 );

//        transformationCache_->setUpdatePartials( );

//        calculatePartials_ = true;
//    }

    std::vector< boost::tuple< unsigned int, unsigned int, unsigned int, unsigned int > > getCoefficientCombinationsToUse( )
    {
        return coefficientCombinationsToUse_;
    }

    Eigen::MatrixXd getCosineCoefficientsOfBody1( )
    {
        return cosineCoefficientsOfBody1_;
    }

    Eigen::MatrixXd getSineCoefficientsOfBody1( )
    {
        return sineCoefficientsOfBody1_;
    }

    double getMultiplier( const int degreeOfBody1, const int orderOfBody1, const int degreeOfBody2, const int orderOfBody2 )
    {
        return multipliers_[ getEffectiveIndex(
                    degreeOfBody1, orderOfBody1, degreeOfBody2, orderOfBody2 ) ] ;
    }

//    void computePartialsOfFullCoefficientsWrtTransformedCoefficients(
//            std::vector< Eigen::Matrix2d >& fullCoefficientsWrtBody2CoefficientsList );

    Eigen::MatrixXd getTransformedCosineCoefficientsOfBody2( )
    {
        return transformedCosineCoefficientsOfBody2_;
    }

    Eigen::MatrixXd getTransformedSineCoefficientsOfBody2( )
    {
        return transformedSineCoefficientsOfBody2_;
    }

private:

    //! Function to return (by reference) single set of effective one-body coefficients
    /*!
     * Function to return (by reference) single set of effective one-body coefficients
     * \param degree1 Degree of expansion for body 1 at which coefficient is to be retrieved
     * \param order1 Order of expansion for body 1 at which coefficient is to be retrieved
     * \param degree2 Degree of expansion for body 2 at which coefficient is to be retrieved
     * \param order2 Order of expansion for body 2 at which coefficient is to be retrieved
     * \param effectiveIndex Index in multipliers_ vector corresponding to coefficient to be retrieved from getEffectiveIndex
     * function
     * \param cosineCoefficient Effective cosine coefficient (returned by reference)
     * \param sineCoefficient Effective sine coefficient (returned by reference)
     */
    void getCurrentEffectiveCoefficients(
            const int degree1, const int order1, const int degree2, const int order2, const int effectiveIndex,
            double& cosineCoefficient, double& sineCoefficient );

    void getCurrentEffectiveAngularMomentumOperatorOfCoefficients(
            const int degree1, const int order1, const int degree2, const int order2,
            const int effectiveIndex,
            Eigen::Vector3d& cosineCoefficient, Eigen::Vector3d& sineCoefficient );

    //! Function that updates the effective one-body coefficients from current member variables of object
    void updateEffectiveMutualPotential( );

    void updateEffectiveAngularMomentumOperatorOfMutualPotential( );

    //! Function to initialize state-independent terms used to compute effective one-body coefficients
    void initializeMultipliers( );

    //! List of combinations of degrees/orders of the two bodies that are to be used, given in order
    //! (degree body 1, order body 1, degree body 2, order body 2)
    std::vector< boost::tuple< unsigned int, unsigned int, unsigned int, unsigned int > > coefficientCombinationsToUse_;

    //! Function returning the cosine coefficients of body 1 in its body-fixed frame
    boost::function< Eigen::MatrixXd( ) > cosineCoefficientFunctionOfBody1_;

    //! Function returning the sine coefficients of body 1 in its body-fixed frame
    boost::function< Eigen::MatrixXd( ) > sineCoefficientFunctionOfBody1_;

    //! Function returning the cosine coefficients of body 2 in its body-fixed frame
    boost::function< Eigen::MatrixXd( ) > cosineCoefficientFunctionOfBody2_;

    //! Function returning the sine coefficients of body 2 in its body-fixed frame
    boost::function< Eigen::MatrixXd( ) > sineCoefficientFunctionOfBody2_;

    //! Current cosine coefficients of body 1 in its body-fixed frame
    Eigen::MatrixXd cosineCoefficientsOfBody1_;

    //! Current sine coefficients of body 1 in its body-fixed frame
    Eigen::MatrixXd sineCoefficientsOfBody1_;

    //! Current cosine coefficients of body 2 in its body-fixed frame
    Eigen::MatrixXd cosineCoefficientsOfBody2_;

    //! Current sine coefficients of body 2 in its body-fixed frame
    Eigen::MatrixXd sineCoefficientsOfBody2_;

    //! Current cosine coefficients of body 2, transformed to body-fixed frame of body 1
    Eigen::MatrixXd transformedCosineCoefficientsOfBody2_;

    //! Current sine coefficients of body 2, transformed to body-fixed frame of body 1
    Eigen::MatrixXd transformedSineCoefficientsOfBody2_;


    std::map< int, std::map< int, Eigen::Vector3d > > currentAngularMomentumProduceCosineCoefficients_;
    std::map< int, std::map< int, Eigen::Vector3d > > currentAngularMomentumProduceSineCoefficients_ ;

    std::vector< Eigen::MatrixXd > transformedCosineCoefficientsOfBody2Partials_;
    std::vector< Eigen::MatrixXd > transformedSineCoefficientsOfBody2Partials_;

    //! Effective one-body cosine coefficients, for given entry of coefficientCombinationsToUse_
    std::vector< double > effectiveCosineCoefficients_;

    //! Effective one-body sine coefficients, for given entry of coefficientCombinationsToUse_
    std::vector< double > effectiveSineCoefficients_;

    std::vector< Eigen::Vector3d > effectiveAngularMomentumOperatorOfCosineCoefficients_;
    std::vector< Eigen::Vector3d > effectiveAngularMomentumOperatorOfSineCoefficients_;

    //! State-independent multipliers to compute given entry of effectiveCosineCoefficients_/effectiveSineCoefficients_
    std::vector< double > multipliers_;

    //! Object used to transform spherical harmonic coefficients of body 2 to body-fixed frame of body 1
    boost::shared_ptr< basic_mathematics::SphericalHarmonicTransformationCache > transformationCache_;

    //! Function returning the effective gravitational parameter
    boost::function< double( ) > gravitationalParameterFunction_;

    //! Reference (equatorial) radius of body 1
    double equatorialRadiusOfBody1_;

    //! Reference (equatorial) radius of body 2
    double equatorialRadiusOfBody2_;

    //! Boolean denoting whether the coefficients are (4 pi) normalized
    bool areCoefficientsNormalized_;

    bool calculatePartials_;

    //! Maximum degree used for body 1
    int maximumDegree1_;

    //! Maximum degree used for body 2
    int maximumDegree2_;

    bool angularMomentumOperatorsAreSet_;

};


}

}

#endif // TUDAT_MUTUALFORCEPOTENTIAL_H
