/*    Copyright (c) 2010-2018, Delft University of Technology
 *    All rigths reserved
 *
 *    This file is part of the Tudat. Redistribution and use in source and
 *    binary forms, with or without modification, are permitted exclusively
 *    under the terms of the Modified BSD license. You should have received
 *    a copy of the license with this file. If not, please or visit:
 *    http://tudat.tudelft.nl/LICENSE.
 *
 *    References
 *      Vittaldev, V. (2010). The unified state model: Derivation and application in astrodynamics
 *          and navigation. Master thesis, Delft University of Technology.
 *      Facchinelli, M. (2018). Aerobraking Guidance, Navigation and Control.
 *          Master thesis, Delft University of Technology.
 */

#include <cmath>

#include "Tudat/Astrodynamics/BasicAstrodynamics/orbitalElementConversions.h"
#include "Tudat/Mathematics/BasicMathematics/mathematicalConstants.h"
#include "Tudat/Mathematics/BasicMathematics/basicMathematicsFunctions.h"

#include "Tudat/Astrodynamics/BasicAstrodynamics/missionGeometry.h"
#include "Tudat/Astrodynamics/BasicAstrodynamics/unifiedStateModelModifiedRodriguesParametersElementConversions.h"
#include "Tudat/Astrodynamics/BasicAstrodynamics/stateVectorIndices.h"

#include "Tudat/Astrodynamics/BasicAstrodynamics/unifiedStateModelQuaternionsElementConversions.h"

namespace tudat
{

namespace orbital_element_conversions
{

//! Convert Keplerian elements to unified state model elements with modified Rodrigues parameters.
Eigen::Vector7d convertKeplerianToUnifiedStateModelModifiedRodriguesParametersElements(
        const Eigen::Vector6d& keplerianElements,
        const double centralBodyGravitationalParameter )
{
    using mathematical_constants::PI;

    // Declaring eventual output vector.
    Eigen::Vector7d convertedUnifiedStateModelElements = Eigen::Vector7d::Zero( );

    // Define the tolerance of a singularity
    const double singularityTolerance = 20.0 * std::numeric_limits< double >::epsilon( );

    // If eccentricity is outside range [0,inf)
    if ( keplerianElements( eccentricityIndex ) < 0.0 )
    {
        //Define the error message
        std::stringstream errorMessage;
        errorMessage << "Eccentricity is expected in range [0,inf)\n"
                     << "Specified eccentricity: " << keplerianElements( eccentricityIndex ) << std::endl;

        // Throw exception
        throw std::runtime_error( std::runtime_error( errorMessage.str( ) ) );
    }

    // If inclination is outside range [0,PI]
    if ( ( keplerianElements( inclinationIndex ) < 0.0 ) || ( keplerianElements( inclinationIndex ) > PI ) )
    {
        // Define the error message.
        std::stringstream errorMessage;
        errorMessage << "Inclination is expected in range [0," << PI << "]\n"
                     << "Specified inclination: " << keplerianElements( inclinationIndex ) << " rad." << std::endl;

        // Throw exception.
        throw std::runtime_error( std::runtime_error( errorMessage.str( ) ) );
    }

    // If argument of pericenter is outside range [0,2.0 * PI]
    if ( ( keplerianElements( argumentOfPeriapsisIndex ) < 0.0 ) || ( keplerianElements( argumentOfPeriapsisIndex ) >
                                                                      2.0 * PI ) )
    {
        // Define the error message.
        std::stringstream errorMessage;
        errorMessage << "Argument of periapsis is expected in range [0," << 2.0 * PI << "]\n"
                     << "Specified inclination: " << keplerianElements( argumentOfPeriapsisIndex ) << " rad." << std::endl;

        // Throw exception.
        throw std::runtime_error( std::runtime_error( errorMessage.str( ) ) );
    }

    // If right ascension of ascending node is outside range [0,2.0 * PI]
    if ( ( keplerianElements( longitudeOfAscendingNodeIndex ) < 0.0 ) ||
         ( keplerianElements( longitudeOfAscendingNodeIndex ) > 2.0 * PI ) )
    {
        // Define the error message.
        std::stringstream errorMessage;
        errorMessage << "RAAN is expected in range [0," << 2.0 * PI << "]\n"
                     << "Specified inclination: " << keplerianElements( longitudeOfAscendingNodeIndex ) << " rad."
                     << std::endl;

        // Throw exception.
        throw std::runtime_error( std::runtime_error( errorMessage.str( ) ) );
    }

    // If true anomaly is outside range [0,2.0 * PI]
    if ( ( keplerianElements( trueAnomalyIndex ) < 0.0 ) || ( keplerianElements( trueAnomalyIndex ) > 2.0 * PI ) )
    {
        // Define the error message.
        std::stringstream errorMessage;
        errorMessage << "True anomaly is expected in range [0," << 2.0 * PI << "]\n"
                     << "Specified inclination: " << keplerianElements( trueAnomalyIndex ) << " rad." << std::endl;

        // Throw exception.
        throw std::runtime_error( std::runtime_error( errorMessage.str( ) ) );
    }

    // If inclination is zero and the right ascension of ascending node is non-zero
    if ( ( std::fabs( keplerianElements( inclinationIndex ) ) < singularityTolerance ) &&
         ( std::fabs( keplerianElements( longitudeOfAscendingNodeIndex ) ) > singularityTolerance ) )
    {
        // Define the error message.
        std::stringstream errorMessage;
        errorMessage << "When the inclination is zero, the right ascending node should be zero by definition\n"
                     << "Specified right ascension of ascending node: " <<
                        keplerianElements( longitudeOfAscendingNodeIndex ) << " rad." << std::endl;

        // Throw exception.
        throw std::runtime_error( std::runtime_error( errorMessage.str( ) ) );
    }

    // If eccentricity is zero and the argument of pericenter is non-zero
    if ( ( std::fabs( keplerianElements( eccentricityIndex ) ) < singularityTolerance ) &&
         ( std::fabs( keplerianElements( argumentOfPeriapsisIndex ) ) > singularityTolerance ) )
    {
        // Define the error message.
        std::stringstream errorMessage;
        errorMessage << "When the eccentricity is zero, the argument of pericenter should be zero by definition\n"
                     << "Specified argument of pericenter: " <<
                        keplerianElements( argumentOfPeriapsisIndex ) << " rad." << std::endl;

        // Throw exception.
        throw std::runtime_error( std::runtime_error( errorMessage.str( ) ) );
    }

    // If semi-major axis is negative and the eccentricity is smaller or equal to one
    if ( ( keplerianElements( semiMajorAxisIndex ) < 0.0 ) && ( keplerianElements( eccentricityIndex ) <= 1.0 ) )
    {
        // Define the error message.
        std::stringstream errorMessage;
        errorMessage << "When the semi-major axis is negative, the eccentricity should be larger than one\n"
                     << "Specified semi-major axis: " << keplerianElements( semiMajorAxisIndex ) << " m.\n"
                     << "Specified eccentricity: " << keplerianElements( eccentricityIndex ) << " rad." << std::endl;

        // Throw exception.
        throw std::runtime_error( std::runtime_error( errorMessage.str( ) ) );
    }

    // If semi-major axis is positive and the eccentricity is larger than one
    if ( ( keplerianElements( semiMajorAxisIndex ) > 0.0 ) && ( keplerianElements( eccentricityIndex ) > 1.0 ) )
    {
        // Define the error message.
        std::stringstream errorMessage;
        errorMessage << "When the semi-major axis is positive, the eccentricity should be smaller than or equal to one\n"
                     << "Specified semi-major axis: " << keplerianElements( semiMajorAxisIndex ) << " m.\n"
                     << "Specified eccentricity: " << keplerianElements( eccentricityIndex ) << " rad." << std::endl;

        // Throw exception.
        throw std::runtime_error( std::runtime_error( errorMessage.str( ) ) );
    }
    // Else, nothing wrong and continue

    // Compute the C hodograph element of the unified state model
    if ( std::fabs( keplerianElements( eccentricityIndex ) - 1.0) < singularityTolerance )
        // parabolic orbit -> semi-major axis is not defined
    {
        convertedUnifiedStateModelElements( CHodographModifiedRodriguesParameterIndex ) =
                std::sqrt( centralBodyGravitationalParameter / keplerianElements( semiLatusRectumIndex ) );
    }
    else
    {
        convertedUnifiedStateModelElements( CHodographModifiedRodriguesParameterIndex ) =
                std::sqrt( centralBodyGravitationalParameter / ( keplerianElements( semiMajorAxisIndex )
                                                                 * ( 1 - keplerianElements( eccentricityIndex ) *
                                                                     keplerianElements( eccentricityIndex ) ) ) );
    }

    // Calculate the additional R hodograph parameter
    double RHodographElement = keplerianElements( eccentricityIndex ) *
            convertedUnifiedStateModelElements( CHodographModifiedRodriguesParameterIndex );

    // Compute the Rf1 hodograph element of the unified state model
    convertedUnifiedStateModelElements( Rf1HodographModifiedRodriguesParameterIndex ) =
            - RHodographElement * std::sin( keplerianElements( longitudeOfAscendingNodeIndex )
                                            + keplerianElements( argumentOfPeriapsisIndex ) );

    // Compute the Rf2 hodograph element of the unified state model
    convertedUnifiedStateModelElements( Rf2HodographModifiedRodriguesParameterIndex ) =
            RHodographElement * std::cos( keplerianElements( longitudeOfAscendingNodeIndex )
                                          + keplerianElements( argumentOfPeriapsisIndex ) );

    // Calculate the additional elements
    double argumentOfLatitude = keplerianElements( argumentOfPeriapsisIndex ) + keplerianElements( trueAnomalyIndex ); // also called u
    double rightAscensionOfLatitude = argumentOfLatitude + keplerianElements( longitudeOfAscendingNodeIndex ); // also called lambda

    // Compute denominator of conversion
    double denominator = 1.0 + std::cos( 0.5 * keplerianElements( inclinationIndex ) ) *
            std::cos( 0.5 * rightAscensionOfLatitude );

    // Decide whether to switch to shadow modified Rodrigues parameters (SMRP) to avoid singularity at 2 PI
    bool shadowFlag = std::fabs( denominator ) < singularityTolerance; // select SMRP directly if denominator is null
    if ( !shadowFlag )
    {
        shadowFlag = denominator < 1.0; // check denominator
    }
    convertedUnifiedStateModelElements( shadowModifiedRodriguesParameterFlagIndex ) = shadowFlag ? 1.0 : 0.0;
    denominator = shadowFlag ? ( denominator - 2.0 ) : denominator; // redefine denominator for SMRP

    // Compute the sigma1 modified Rodrigues parameters of the unified state model
    convertedUnifiedStateModelElements( sigma1ModifiedRodriguesParameterIndex ) =
            std::sin( 0.5 * keplerianElements( inclinationIndex ) ) *
            std::cos( 0.5 * ( keplerianElements( longitudeOfAscendingNodeIndex ) - argumentOfLatitude ) ) / denominator;

    // Compute the sigma2 modified Rodrigues parameters of the unified state model
    convertedUnifiedStateModelElements( sigma2ModifiedRodriguesParameterIndex ) =
            std::sin( 0.5 * keplerianElements( inclinationIndex ) ) *
            std::sin( 0.5 * ( keplerianElements( longitudeOfAscendingNodeIndex ) - argumentOfLatitude ) ) / denominator;

    // Compute the sigma3 modified Rodrigues parameters of the unified state model
    convertedUnifiedStateModelElements( sigma3ModifiedRodriguesParameterIndex ) =
            std::cos( 0.5 * keplerianElements( inclinationIndex ) ) *
            std::sin( 0.5 * rightAscensionOfLatitude ) / denominator;

    // Give back result
    return convertedUnifiedStateModelElements;
}

//! Convert unified state model elements with modified Rodrigues parameters to Keplerian elements.
Eigen::Vector6d convertUnifiedStateModelModifiedRodriguesParametersToKeplerianElements(
        const Eigen::Vector7d& unifiedStateModelElements,
        const double centralBodyGravitationalParameter )
{
    using mathematical_constants::PI;

    // Declaring eventual output vector.
    Eigen::Vector6d convertedKeplerianElements = Eigen::Vector6d::Zero( );

    // Define the tolerance of a singularity
    const double singularityTolerance = 20.0 * std::numeric_limits< double >::epsilon( );

    // Extract modified Rodrigues parameters
    Eigen::Vector3d modifiedRodriguesParametersVector = unifiedStateModelElements.segment( sigma1ModifiedRodriguesParameterIndex, 3 );
    double modifiedRodriguesParametersMagnitude =
            modifiedRodriguesParametersVector.norm( ); // magnitude of modified Rodrigues parameters, also called sigma
    bool shadowFlag = int( unifiedStateModelElements( shadowModifiedRodriguesParameterFlagIndex ) ) == 1;

    // Precompute often used variables
    // Note the for SMRP some variable names do not match their definitions
    double modifiedRodriguesParametersMagnitudeSquared =
            modifiedRodriguesParametersMagnitude * modifiedRodriguesParametersMagnitude;
    double onePlusModifiedRodriguesParametersMagnitudeSquaredSquared =
            std::pow( 1.0 + modifiedRodriguesParametersMagnitudeSquared, 2 );
    double oneMinusModifiedRodriguesParametersMagnitudeSquared = shadowFlag ?
                ( modifiedRodriguesParametersMagnitudeSquared - 1.0 ) : // inverse definition for SMRP
                ( 1.0 - modifiedRodriguesParametersMagnitudeSquared );
    double oneMinusModifiedRodriguesParametersMagnitudeSquaredSquared =
            std::pow( oneMinusModifiedRodriguesParametersMagnitudeSquared, 2 );

    // Compute right ascension of latitude (independent on shadowing, since variables have been precomputed)
    double denominator = 4.0 * std::pow( modifiedRodriguesParametersVector( 2 ), 2 ) +
            oneMinusModifiedRodriguesParametersMagnitudeSquaredSquared; // denominator is never null
    double rightAscensionOfLatitude =
            std::atan2( 4.0 * modifiedRodriguesParametersVector( 2 ) *
                        ( 1.0 - modifiedRodriguesParametersMagnitudeSquared ) / denominator,
                        ( oneMinusModifiedRodriguesParametersMagnitudeSquaredSquared -
                          4.0 * std::pow( modifiedRodriguesParametersVector( 2 ), 2 ) ) /
                        denominator ); // also known as lambda

    // Trigonometric values of right ascension of latitude
    double cosineLambda = std::cos( rightAscensionOfLatitude );
    double sineLambda = std::sin( rightAscensionOfLatitude );

    // Compute auxiliary parameters auxiliaryParameter1 and auxiliaryParameter2
    double auxiliaryParameter1 = unifiedStateModelElements( Rf1HodographModifiedRodriguesParameterIndex ) * cosineLambda +
            unifiedStateModelElements( Rf2HodographModifiedRodriguesParameterIndex ) * sineLambda;
    double auxiliaryParameter2 = unifiedStateModelElements( CHodographModifiedRodriguesParameterIndex ) -
            unifiedStateModelElements( Rf1HodographModifiedRodriguesParameterIndex ) * sineLambda +
            unifiedStateModelElements( Rf2HodographModifiedRodriguesParameterIndex ) * cosineLambda;

    // Compute auxiliary R hodograph parameter
    double RHodographElement = std::sqrt( unifiedStateModelElements( Rf1HodographModifiedRodriguesParameterIndex ) *
                                          unifiedStateModelElements( Rf1HodographModifiedRodriguesParameterIndex ) +
                                          unifiedStateModelElements( Rf2HodographModifiedRodriguesParameterIndex ) *
                                          unifiedStateModelElements( Rf2HodographModifiedRodriguesParameterIndex ) );

    // Compute eccentricity
    convertedKeplerianElements( eccentricityIndex ) =
            RHodographElement / unifiedStateModelElements( CHodographModifiedRodriguesParameterIndex );

    // Compute semi-major axis or, in case of a parabolic orbit, the semi-latus rectum.
    if ( std::fabs( convertedKeplerianElements( eccentricityIndex ) - 1.0 ) < singularityTolerance )
        // parabolic orbit -> semi-major axis is not defined. Use semi-latus rectum instead.
    {
        convertedKeplerianElements( semiLatusRectumIndex ) = centralBodyGravitationalParameter /
                ( unifiedStateModelElements( CHodographModifiedRodriguesParameterIndex ) *
                  unifiedStateModelElements( CHodographModifiedRodriguesParameterIndex ) );
    }
    else
    {
        convertedKeplerianElements( semiMajorAxisIndex ) =
                centralBodyGravitationalParameter /
                ( std::pow( unifiedStateModelElements( CHodographModifiedRodriguesParameterIndex ), 2 ) *
                  ( 1 - std::pow( convertedKeplerianElements( eccentricityIndex ), 2 ) ) );
    }

    // Compute inclination
    double arccosineArgument = ( 4.0 * ( std::pow( modifiedRodriguesParametersVector( 2 ), 2 ) -
                                         std::pow( modifiedRodriguesParametersVector( 0 ), 2 ) -
                                         std::pow( modifiedRodriguesParametersVector( 1 ), 2 ) ) +
                                 oneMinusModifiedRodriguesParametersMagnitudeSquaredSquared ) /
            onePlusModifiedRodriguesParametersMagnitudeSquaredSquared;
    if ( ( ( std::fabs( arccosineArgument ) - 1.0 ) > 0.0 ) ||
         ( std::abs( std::fabs( arccosineArgument ) - 1.0 ) < singularityTolerance ) )
    {
        // Make sure that the cosine does not exceed 1.0 in magnitude
        arccosineArgument = ( arccosineArgument > 0.0 ) ? 1.0 : - 1.0;
    }
    convertedKeplerianElements( inclinationIndex ) = std::acos( arccosineArgument );
    // this acos is always defined correctly because the inclination is always below PI rad

    // Find sine and cosine of longitude of ascending node separately
    double signDirectionCosineMatrix = shadowFlag ? - 1.0 : 1.0; // sign depends on shadow conditions
    double sineOmega = 8.0 * modifiedRodriguesParametersVector( 0 ) *
            modifiedRodriguesParametersVector( 2 ) + signDirectionCosineMatrix *
            4.0 * modifiedRodriguesParametersVector( 1 ) *
            oneMinusModifiedRodriguesParametersMagnitudeSquared;
    double cosineOmega = - 8.0 * modifiedRodriguesParametersVector( 1 ) *
            modifiedRodriguesParametersVector( 2 ) + signDirectionCosineMatrix *
            4.0 * modifiedRodriguesParametersVector( 0 ) *
            oneMinusModifiedRodriguesParametersMagnitudeSquared;
    denominator = std::sqrt( cosineOmega * cosineOmega +
                             sineOmega * sineOmega ); // overwrite

    // Compute longitude of ascending node
    if ( std::fabs( convertedKeplerianElements( inclinationIndex ) - PI ) < singularityTolerance )
        // pure-retrograde orbit -> inclination = PI
    {
        // Define the error message
        std::stringstream errorMessage;
        errorMessage << "Pure-retrograde orbit (inclination = PI).\n"
                     << "Unified state model elements cannot be transformed to Kepler elements." << std::endl;

        // Throw exception
        throw std::runtime_error( std::runtime_error( errorMessage.str( ) ) );
    }
    else if ( ( std::fabs( convertedKeplerianElements( inclinationIndex ) ) < singularityTolerance ) ||
        // pure-prograde orbit -> inclination = 0 -> RAAN is undefined, hence give zero
              ( std::fabs( denominator ) < singularityTolerance ) )
        // denominator is null, find a workaround
    {
        convertedKeplerianElements( longitudeOfAscendingNodeIndex ) = 0.0;
    }
    else
    {
        // Compute longitude of ascending node
        convertedKeplerianElements( longitudeOfAscendingNodeIndex ) = std::atan2(
                    sineOmega / denominator, cosineOmega / denominator );

        // Round off small values of the right ascension of ascending node to zero
        if ( std::fabs( convertedKeplerianElements( longitudeOfAscendingNodeIndex ) ) < singularityTolerance )
        {
            convertedKeplerianElements( longitudeOfAscendingNodeIndex ) = 0.0;
        }

        // Ensure the longitude of ascending node is positive
        while ( convertedKeplerianElements( longitudeOfAscendingNodeIndex ) < 0.0 )
            // Because of the previous if statement, if the longitude of ascending node is smaller than 0, it will
            // always be smaller than -singularityTolerance
        {
            convertedKeplerianElements( longitudeOfAscendingNodeIndex ) += 2.0 * PI;
        }
    }

    // Compute true anomaly and argument of periapsis
    if ( std::fabs( RHodographElement ) < singularityTolerance ) // circular orbit
    {
        convertedKeplerianElements( argumentOfPeriapsisIndex ) = 0.0; // by definition
        convertedKeplerianElements( trueAnomalyIndex ) =
                rightAscensionOfLatitude - convertedKeplerianElements( longitudeOfAscendingNodeIndex );

        // Round off small theta to zero
        if ( std::fabs( convertedKeplerianElements( trueAnomalyIndex ) ) < singularityTolerance )
        {
            convertedKeplerianElements( trueAnomalyIndex ) = 0.0;
        }

        // Ensure the true anomaly is positive
        while ( convertedKeplerianElements( trueAnomalyIndex ) < 0.0 )
            // Because of the previous if statement, if the true anomaly is smaller than zero, it will always be smaller than
            // -singularityTolerance
        {
            convertedKeplerianElements( trueAnomalyIndex ) += 2.0 * PI;
        }
    }
    else
    {
        convertedKeplerianElements( trueAnomalyIndex ) =
                std::atan2( ( auxiliaryParameter1 / RHodographElement ),
                            ( ( auxiliaryParameter2 - unifiedStateModelElements( CHodographModifiedRodriguesParameterIndex ) )
                              / RHodographElement ) );

        // Round off small theta to zero
        if ( std::fabs( convertedKeplerianElements( trueAnomalyIndex ) ) < singularityTolerance )
        {
            convertedKeplerianElements( trueAnomalyIndex ) = 0.0;
        }

        // Ensure the true anomaly is positive
        while ( convertedKeplerianElements( trueAnomalyIndex ) < 0.0 )
            // Because of the previous if statement, if the true anomaly is smaller than zero, it will always
            // be smaller than -singularityTolerance
        {
            convertedKeplerianElements( trueAnomalyIndex ) += 2.0 * PI;
        }

        convertedKeplerianElements( argumentOfPeriapsisIndex ) =
                rightAscensionOfLatitude -
                convertedKeplerianElements( longitudeOfAscendingNodeIndex ) -
                convertedKeplerianElements( trueAnomalyIndex );

        // Round off small omega to zero
        if ( std::fabs( convertedKeplerianElements( argumentOfPeriapsisIndex ) ) < singularityTolerance )
        {
            convertedKeplerianElements( argumentOfPeriapsisIndex ) = 0.0;
        }

        // Ensure the argument of periapsis is positive
        while ( convertedKeplerianElements( argumentOfPeriapsisIndex ) < 0.0 )
            // Because of the previous if statement, if the argument of pericenter is smaller than zero,
            // it will be smaller than -singularityTolerance
        {
            convertedKeplerianElements( argumentOfPeriapsisIndex ) += 2.0 * PI;
        }
    }

    // Give back result
    return convertedKeplerianElements;
}

////! Convert Cartesian elements to unified state model elements with modified Rodrigues parameters.
Eigen::Vector7d convertCartesianToUnifiedStateModelModifiedRodriguesParametersElements(
        const Eigen::Vector6d& cartesianElements,
        const double centralBodyGravitationalParameter )
{
    // Declaring eventual output vector.
    Eigen::Vector7d convertedUnifiedStateModelModifiedRodriguesParametersElements = Eigen::Vector7d::Zero( );

    // Convert Cartesian to USM7
    Eigen::Vector7d unifiedStateModelQuaternionsElements =
            convertCartesianToUnifiedStateModelQuaternionsElements( cartesianElements,
                                                                    centralBodyGravitationalParameter );

    // Extract quaternions
    Eigen::Vector3d epsilonQuaternionVector =
            unifiedStateModelQuaternionsElements.segment( epsilon1QuaternionIndex, 3 );
    double etaQuaternionParameter = unifiedStateModelQuaternionsElements( etaQuaternionIndex );

    // Convert quaternions to modified Rodrigues parameters (or SMPR)
    bool shadowFlag = etaQuaternionParameter < 0;
    double conversionSign = shadowFlag ? - 1.0 : 1.0; // conversion is slightly different for SMRP and MRP
    convertedUnifiedStateModelModifiedRodriguesParametersElements.segment(
                sigma1ModifiedRodriguesParameterIndex, 3 ) = conversionSign * epsilonQuaternionVector /
            ( 1 + conversionSign * etaQuaternionParameter );
    convertedUnifiedStateModelModifiedRodriguesParametersElements(
                shadowModifiedRodriguesParameterFlagIndex ) = shadowFlag ? 1.0 : 0.0;

    // Add other elements to USMEM vector
    convertedUnifiedStateModelModifiedRodriguesParametersElements.segment(
                CHodographModifiedRodriguesParameterIndex, 3 ) =
            unifiedStateModelQuaternionsElements.segment( CHodographQuaternionIndex, 3 );

    // Give back result
    return convertedUnifiedStateModelModifiedRodriguesParametersElements;
}

//! Convert unified state model elements with modified Rodrigues parameters to Cartesian elements.
Eigen::Vector6d convertUnifiedStateModelModifiedRodriguesParametersToCartesianElements(
        const Eigen::Vector7d& unifiedStateModelModifiedRodriguesParametersElements,
        const double centralBodyGravitationalParameter )
{
    // Extract modified Rodrigues parameters
    Eigen::Vector3d modifiedRodriguesParametersVector =
            unifiedStateModelModifiedRodriguesParametersElements.segment( sigma1ModifiedRodriguesParameterIndex, 3 );
    double modifiedRodriguesParametersMagnitude =
            modifiedRodriguesParametersVector.norm( ); // magnitude of modified Rodrigues parameters, also called sigma
    bool shadowFlag =
            int( unifiedStateModelModifiedRodriguesParametersElements( shadowModifiedRodriguesParameterFlagIndex ) ) == 1;

    // Precompute often used variables
    double modifiedRodriguesParametersMagnitudeSquared = modifiedRodriguesParametersMagnitude *
            modifiedRodriguesParametersMagnitude;

    // Create USM7 vector and add velocity hodograph elements
    Eigen::Vector7d unifiedStateModelQuaternionsElements;
    unifiedStateModelQuaternionsElements.segment( CHodographQuaternionIndex, 3 ) =
            unifiedStateModelModifiedRodriguesParametersElements.segment( CHodographExponentialMapIndex, 3 );

    // Convert modified Rodrigues parameters to quaternions
    double conversionSign = shadowFlag ? - 1.0 : 1.0; // converion is slightly different for SMRP and MRP
    Eigen::Vector3d epsilonQuaternionVector = conversionSign *
            2.0 / ( 1.0 + modifiedRodriguesParametersMagnitudeSquared ) *
            modifiedRodriguesParametersVector;
    double etaQuaternionParameter = conversionSign * ( 1.0 - modifiedRodriguesParametersMagnitudeSquared ) /
            ( 1.0 + modifiedRodriguesParametersMagnitudeSquared );

    // Add quaternions to USM7 vector
    unifiedStateModelQuaternionsElements.segment( epsilon1QuaternionIndex, 3 ) = epsilonQuaternionVector;
    unifiedStateModelQuaternionsElements( etaQuaternionIndex ) = etaQuaternionParameter;

    // Give back result
    return convertUnifiedStateModelQuaternionsToCartesianElements( unifiedStateModelQuaternionsElements,
                                                                   centralBodyGravitationalParameter );
}

} // close namespace orbital_element_conversions

} // close namespace tudat