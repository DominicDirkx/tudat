/*    Copyright (c) 2010-2017, Delft University of Technology
 *    All rigths reserved
 *
 *    This file is part of the Tudat. Redistribution and use in source and
 *    binary forms, with or without modification, are permitted exclusively
 *    under the terms of the Modified BSD license. You should have received
 *    a copy of the license with this file. If not, please or visit:
 *    http://tudat.tudelft.nl/LICENSE.
 *
 *    References
 *      Chobotov, V.A. Orbital Mechanics, Third Edition, AIAA Education Series, VA, 2002.
 *      Wertz, J. R. Mission geometry; orbit and constellation design and management.
 *      Mengali, G., Quarta, A.A. Fondamenti di meccanica del volo spaziale.
 *      Wertz, J.R. Mission Geometry; Orbit and Constellation Design and Management, Spacecraft
 *          Orbit and Attitude Systems, Microcosm Press, Kluwer Academic Publishers, 2001.
 *      Advanced Concepts Team, ESA. Keplerian Toolbox, http://sourceforge.net/projects/keptoolbox,
 *          last accessed: 21st April, 2012.
 *
 *    Notes
 *      Backwards compatibility of namespaces is implemented for Tudat Core 2 in this file. The
 *      code block marked "DEPRECATED!" at the end of the file should be removed in Tudat Core 3.
 *
 */

#ifndef TUDAT_ORBITAL_ELEMENT_CONVERSIONS_H
#define TUDAT_ORBITAL_ELEMENT_CONVERSIONS_H

#include <boost/exception/all.hpp>
#include <boost/math/special_functions/atanh.hpp>

#include <iostream>
#include <cmath>
#include <limits>

#include <Eigen/Core>
#include <Eigen/Geometry>

#include "Tudat/Astrodynamics/BasicAstrodynamics/stateVectorIndices.h"
#include "Tudat/Mathematics/BasicMathematics/mathematicalConstants.h"
#include "Tudat/Mathematics/BasicMathematics/basicMathematicsFunctions.h"

namespace tudat
{

namespace orbital_element_conversions
{

//! Convert Keplerian to Cartesian orbital elements.
/*!
 * Converts Keplerian to Cartesian orbital elements (Chobotov, 2002). Use the
 * CartesianElementVectorIndices enum to access the individual orbital element components in the
 * storage vector.
 *
 * \param keplerianElements Vector containing Keplerian elements.                         \n
 *          <em>
 *                          Order of elements is important! \n
 *                          keplerianElements( 0 ) = semiMajorAxis,                   [m] \n
 *                          keplerianElements( 1 ) = eccentricity,                    [-] \n
 *                          keplerianElements( 2 ) = inclination,                   [rad] \n
 *                          keplerianElements( 3 ) = argument of periapsis,         [rad] \n
 *                          keplerianElements( 4 ) = longitude of ascending node,   [rad] \n
 *                          keplerianElements( 5 ) = true anomaly.                  [rad] \n
 *          </em>
 *        WARNING: If eccentricity is 1.0 within machine precision,
 *        keplerianElements( 0 ) = semi-latus rectum.
 *
 * \param centralBodyGravitationalParameter Gravitational parameter of central body [m^3 s^-2].
 *
 * \return Converted state in Cartesian elements.                         \n
 *         <em>
 *         Order of elements is important!                                \n
 *         cartesianElements( 0 ) = x-position coordinate,            [m] \n
 *         cartesianElements( 1 ) = y-position coordinate,            [m] \n
 *         cartesianElements( 2 ) = z-position coordinate,            [m] \n
 *         cartesianElements( 3 ) = x-velocity coordinate,          [m/s] \n
 *         cartesianElements( 4 ) = y-velocity coordinate,          [m/s] \n
 *         cartesianElements( 5 ) = z-velocity coordinate.          [m/s] \n
 *         </em>
 *
 * \sa CartesianElementVectorIndices()
 */
template< typename ScalarType = double >
Eigen::Matrix< ScalarType, 6, 1 > convertKeplerianToCartesianElements(
        const Eigen::Matrix< ScalarType, 6, 1 >& keplerianElements,
        const ScalarType centralBodyGravitationalParameter )
{
    using std::cos;
    using std::fabs;
    using std::pow;
    using std::sin;
    using std::sqrt;

    // Set tolerance.
    const ScalarType tolerance_ = std::numeric_limits< ScalarType >::epsilon( );

    // Set local keplerian elements.
    ScalarType semiMajorAxis_ = keplerianElements( semiMajorAxisIndex );
    ScalarType eccentricity_ = keplerianElements( eccentricityIndex );
    ScalarType inclination_ = keplerianElements( inclinationIndex );
    ScalarType argumentOfPeriapsis_ = keplerianElements( argumentOfPeriapsisIndex );
    ScalarType longitudeOfAscendingNode_ = keplerianElements( longitudeOfAscendingNodeIndex );
    ScalarType trueAnomaly_ = keplerianElements( trueAnomalyIndex );

    // Pre-compute sines and cosines of involved angles for efficient computation.
    ScalarType cosineOfInclination_ = cos( inclination_ );
    ScalarType sineOfInclination_ = sin( inclination_ );
    ScalarType cosineOfArgumentOfPeriapsis_ = cos( argumentOfPeriapsis_ );
    ScalarType sineOfArgumentOfPeriapsis_ = sin( argumentOfPeriapsis_ );
    ScalarType cosineOfLongitudeOfAscendingNode_ = cos( longitudeOfAscendingNode_ );
    ScalarType sineOfLongitudeOfAscendingNode_ = sin( longitudeOfAscendingNode_ );
    ScalarType cosineOfTrueAnomaly_ = cos( trueAnomaly_ );
    ScalarType sineOfTrueAnomaly_ = sin( trueAnomaly_ );

    // Declare semi-latus rectum.
    ScalarType semiLatusRectum_ = -mathematical_constants::getFloatingInteger< ScalarType >( 0 );

    // Compute semi-latus rectum in the case it is not a parabola.
    if ( fabs( eccentricity_ - mathematical_constants::getFloatingInteger< ScalarType >( 1 ) ) >
         tolerance_  )
    {
        semiLatusRectum_ = semiMajorAxis_ * (
                    mathematical_constants::getFloatingInteger< ScalarType >( 1 ) -
                    eccentricity_ * eccentricity_ );
    }

    // Else set the semi-latus rectum given for a parabola as the first element in the vector
    // of Keplerian elements..
    else
    {
        semiLatusRectum_ = semiMajorAxis_;
    }

    // Definition of position in the perifocal coordinate system.
    Eigen::Matrix< ScalarType, 2, 1 > positionPerifocal_;
    positionPerifocal_.x( ) = semiLatusRectum_ * cosineOfTrueAnomaly_
            / ( mathematical_constants::getFloatingInteger< ScalarType >( 1 ) +
                eccentricity_ * cosineOfTrueAnomaly_ );
    positionPerifocal_.y( ) = semiLatusRectum_ * sineOfTrueAnomaly_
            / ( mathematical_constants::getFloatingInteger< ScalarType >( 1 ) +
                eccentricity_ * cosineOfTrueAnomaly_ );

    // Definition of velocity in the perifocal coordinate system.
    Eigen::Matrix< ScalarType, 2, 1 > velocityPerifocal_(
                -sqrt( centralBodyGravitationalParameter / semiLatusRectum_ ) * sineOfTrueAnomaly_,
                sqrt( centralBodyGravitationalParameter / semiLatusRectum_ )
                * ( eccentricity_ + cosineOfTrueAnomaly_ ) );

    // Definition of the transformation matrix.
    Eigen::Matrix< ScalarType, 3, 2 > transformationMatrix_;

    // Compute the transformation matrix.
    transformationMatrix_( 0, 0 ) = cosineOfLongitudeOfAscendingNode_
            * cosineOfArgumentOfPeriapsis_ -sineOfLongitudeOfAscendingNode_
            * sineOfArgumentOfPeriapsis_ * cosineOfInclination_;
    transformationMatrix_( 0, 1 ) = -cosineOfLongitudeOfAscendingNode_
            * sineOfArgumentOfPeriapsis_ -sineOfLongitudeOfAscendingNode_
            * cosineOfArgumentOfPeriapsis_ * cosineOfInclination_;
    transformationMatrix_( 1, 0 ) = sineOfLongitudeOfAscendingNode_
            * cosineOfArgumentOfPeriapsis_ + cosineOfLongitudeOfAscendingNode_
            * sineOfArgumentOfPeriapsis_ * cosineOfInclination_;
    transformationMatrix_( 1, 1 ) = -sineOfLongitudeOfAscendingNode_
            * sineOfArgumentOfPeriapsis_ + cosineOfLongitudeOfAscendingNode_
            * cosineOfArgumentOfPeriapsis_ * cosineOfInclination_;
    transformationMatrix_( 2, 0 ) = sineOfArgumentOfPeriapsis_ * sineOfInclination_;
    transformationMatrix_( 2, 1 ) = cosineOfArgumentOfPeriapsis_ * sineOfInclination_;

    // Declare converted Cartesian elements.
    Eigen::Matrix< ScalarType, 6, 1 > convertedCartesianElements_;

    // Compute value of position in Cartesian coordinates.
    Eigen::Matrix< ScalarType, 3, 1 > position_ = transformationMatrix_ * positionPerifocal_;
    convertedCartesianElements_.segment( 0, 3 ) = position_;

    // Compute value of velocity in Cartesian coordinates.
    Eigen::Matrix< ScalarType, 3, 1 > velocity_ = transformationMatrix_ * velocityPerifocal_;
    convertedCartesianElements_.segment( 3, 3 ) = velocity_;

    // Return Cartesian elements.
    return convertedCartesianElements_;
}


//! Convert Cartesian to Keplerian orbital elements.
/*!
 * Converts Cartesian to Keplerian orbital elements.
 * \param cartesianElements Vector containing Cartesian elements. Order of elements is important!
 *          cartesianElements( 0 ) = x-position coordinate,                                     [m]
 *          cartesianElements( 1 ) = y-position coordinate,                                     [m]
 *          cartesianElements( 2 ) = z-position coordinate,                                     [m]
 *          cartesianElements( 3 ) = x-velocity coordinate,                                   [m/s]
 *          cartesianElements( 4 ) = y-velocity coordinate,                                   [m/s]
 *          cartesianElements( 5 ) = z-velocity coordinate.                                   [m/s]
 * \param centralBodyGravitationalParameter Gravitational parameter of central body.
 * \return Converted state in Keplerian elements. The order of elements is fixed!
 *          keplerianElements( 0 ) = semiMajorAxis,                                             [m]
 *          keplerianElements( 1 ) = eccentricity,                                              [-]
 *          keplerianElements( 2 ) = inclination,                                             [rad]
 *          keplerianElements( 3 ) = argument of periapsis,                                   [rad]
 *          keplerianElements( 4 ) = longitude of ascending node,                             [rad]
 *          keplerianElements( 5 ) = true anomaly.                                            [rad]
 *          WARNING: If eccentricity is 1.0 within 1.0e-15,
 *          keplerianElements( 0 ) = semi-latus rectum, since the orbit is parabolic.
 *          WARNING: If eccentricity is 0.0 within 1.0e-15,
 *          argument of periapsis is set to 0.0, since the orbit is circular.
 *          WARNING: If inclination is 0.0 within 1.0e-15,
 *          longitude of ascending node is set to 0.0, since the orbit is equatorial.
 *          The tolerance 1.0e-15 is hard-coded, as it should not be changed for performance
 *          reasons, unless required for specific scenarios. In those cases, users are expected
 *          to update the internal tolerance to the required value. Below these tolerance values
 *          for eccentricity and inclination, the orbit is considered to be a limit case.
 *          Essentially, special solutions are then used for parabolic, circular inclined,
 *          non-circular equatorial, and circular equatorial orbits. These special solutions are
 *          required because of singularities in the classical Keplerian elements. If high
 *          precision is required near these singularities, users are encouraged to consider using
 *          other elements, such as modified equinoctial elements. It should be noted that
 *          modified equinoctial elements also suffer from singularities, but not for zero
 *          eccentricity and inclination.
 */
template< typename ScalarType = double >
Eigen::Matrix< ScalarType, 6, 1 > convertCartesianToKeplerianElements(
        const Eigen::Matrix< ScalarType, 6, 1 >& cartesianElements,
        const ScalarType centralBodyGravitationalParameter )
{
    // Set tolerance.
    const ScalarType tolerance = 20.0 * std::numeric_limits< ScalarType >::epsilon( );

    // Declare converted Keplerian elements.
    Eigen::Matrix< ScalarType, 6, 1 > computedKeplerianElements_;

    // Set position and velocity vectors.
    const Eigen::Matrix< ScalarType, 3, 1 > position_( cartesianElements.segment( 0, 3 ) );
    const Eigen::Matrix< ScalarType, 3, 1 > velocity_( cartesianElements.segment( 3, 3 ) );

    // Compute orbital angular momentum vector.
    const Eigen::Matrix< ScalarType, 3, 1 > angularMomentum_( position_.cross( velocity_ ) );

    // Compute semi-latus rectum.
    const ScalarType semiLatusRectum_ = angularMomentum_.squaredNorm( )
            / centralBodyGravitationalParameter;

    // Compute unit vector to ascending node.
    Eigen::Matrix< ScalarType, 3, 1 > unitAscendingNodeVector_(
                ( Eigen::Matrix< ScalarType, 3, 1 >::UnitZ( ).cross(
                      angularMomentum_.normalized( ) ) ).normalized( ) );

    // Compute eccentricity vector.
    Eigen::Matrix< ScalarType, 3, 1 > eccentricityVector_(
                velocity_.cross( angularMomentum_ ) / centralBodyGravitationalParameter
                - position_.normalized( ) );

    // Store eccentricity.
    computedKeplerianElements_( eccentricityIndex ) = eccentricityVector_.norm( );

    // Compute and store semi-major axis.
    // Check if orbit is parabolic. If it is, store the semi-latus rectum instead of the
    // semi-major axis.
    if ( std::fabs( computedKeplerianElements_( eccentricityIndex ) -
                    mathematical_constants::getFloatingInteger< ScalarType >( 1 ) ) < tolerance )
    {
        computedKeplerianElements_( semiLatusRectumIndex ) = semiLatusRectum_;
    }

    // Else the orbit is either elliptical or hyperbolic, so store the semi-major axis.
    else
    {
        computedKeplerianElements_( semiMajorAxisIndex ) = semiLatusRectum_
                / ( mathematical_constants::getFloatingInteger< ScalarType >( 1 ) -
                    computedKeplerianElements_( eccentricityIndex )
                    * computedKeplerianElements_( eccentricityIndex ) );
    }

    // Compute and store inclination.
    computedKeplerianElements_( inclinationIndex ) = std::acos( angularMomentum_.z( )
                                                                / angularMomentum_.norm( ) );

    // Compute and store longitude of ascending node.
    // Define the quadrant condition for the argument of perigee.
    ScalarType argumentOfPeriapsisQuandrantCondition = eccentricityVector_.z( );

    // Check if the orbit is equatorial. If it is, set the vector to the line of nodes to the
    // x-axis.
    if ( std::fabs( computedKeplerianElements_( inclinationIndex ) ) < tolerance )
    {
        unitAscendingNodeVector_ = Eigen::Matrix< ScalarType, 3, 1 >::UnitX( );

        // If the orbit is equatorial, eccentricityVector_.z( ) is zero, therefore the quadrant
        // condition is taken to be the y-component, eccentricityVector_.y( ).
        argumentOfPeriapsisQuandrantCondition = eccentricityVector_.y( );
    }

    // Compute and store the resulting longitude of ascending node.
    computedKeplerianElements_( longitudeOfAscendingNodeIndex )
            = std::acos( unitAscendingNodeVector_.x( ) );

    // Check if the quandrant is correct.
    if ( unitAscendingNodeVector_.y( ) <
         mathematical_constants::getFloatingInteger< ScalarType >( 0 ) )
    {
        computedKeplerianElements_( longitudeOfAscendingNodeIndex ) =
                mathematical_constants::getFloatingInteger< ScalarType >( 2 ) *
                mathematical_constants::getPi< ScalarType >( ) -
                computedKeplerianElements_( longitudeOfAscendingNodeIndex );
    }

    // Compute and store argument of periapsis.
    // Define the quadrant condition for the true anomaly.
    ScalarType trueAnomalyQuandrantCondition = position_.dot( velocity_ );

    // Check if the orbit is circular. If it is, set the eccentricity vector to unit vector
    // pointing to the ascending node, i.e. set the argument of periapsis to zero.
    if ( std::fabs( computedKeplerianElements_( eccentricityIndex ) ) < tolerance )
    {
        eccentricityVector_ = unitAscendingNodeVector_;

        computedKeplerianElements_( argumentOfPeriapsisIndex ) =
                mathematical_constants::getFloatingInteger< ScalarType >( 0 );

        // Check if orbit is also equatorial and set true anomaly quandrant check condition
        // accordingly.
        if ( unitAscendingNodeVector_ == Eigen::Matrix< ScalarType, 3, 1 >::UnitX( ) )
        {
            // If the orbit is circular, position_.dot( velocity_ ) = 0, therefore this value
            // cannot be used as a quadrant condition. Moreover, if the orbit is equatorial,
            // position_.z( ) is also zero and therefore the quadrant condition is taken to be the
            // y-component, position_.y( ).
            trueAnomalyQuandrantCondition = position_.y( );
        }

        else
        {
            // If the orbit is circular, position_.dot( velocity_ ) = 0, therefore the quadrant
            // condition is taken to be the z-component of the position, position_.z( ).
            trueAnomalyQuandrantCondition = position_.z( );
        }
    }

    // Else, compute the argument of periapsis as the angle between the eccentricity vector and
    // the unit vector to the ascending node.
    else
    {
        ScalarType eccentricityAscendingNodeDotProduct = eccentricityVector_.normalized( ).dot( unitAscendingNodeVector_ );

        // Check whether dot product is in bounds (might be out of bounds due to numerical noise).
        if( eccentricityAscendingNodeDotProduct < mathematical_constants::getFloatingInteger< ScalarType >( -1 ) )
        {
            computedKeplerianElements_( argumentOfPeriapsisIndex ) = mathematical_constants::getPi< ScalarType >( );
        }
        else if( eccentricityAscendingNodeDotProduct > mathematical_constants::getFloatingInteger< ScalarType >( 1 ) )
        {
            computedKeplerianElements_( argumentOfPeriapsisIndex ) =
                    mathematical_constants::getFloatingInteger< ScalarType >( 0 );
        }
        else
        {
             computedKeplerianElements_( argumentOfPeriapsisIndex ) = std::acos( eccentricityAscendingNodeDotProduct );
        }

        // Check if the quadrant is correct.
        if ( argumentOfPeriapsisQuandrantCondition <
             mathematical_constants::getFloatingInteger< ScalarType >( 0 ) )
        {
            computedKeplerianElements_( argumentOfPeriapsisIndex ) =
                    mathematical_constants::getFloatingInteger< ScalarType >( 2 ) *
                    mathematical_constants::getPi< ScalarType >( ) -
                    computedKeplerianElements_( argumentOfPeriapsisIndex );
        }
    }

    // Compute dot-product of position and eccentricity vectors.
    ScalarType dotProductPositionAndEccentricityVectors
            = position_.normalized( ).dot( eccentricityVector_.normalized( ) );

    // Check if the dot-product is one of the limiting cases: 0.0, -1.0 or 1.0
    // (within prescribed tolerance).
    if ( std::fabs( mathematical_constants::getFloatingInteger< ScalarType >( 1 ) -
                    dotProductPositionAndEccentricityVectors ) < tolerance )
    {
        dotProductPositionAndEccentricityVectors =
                mathematical_constants::getFloatingInteger< ScalarType >( 1 );
    }

    if ( std::fabs( mathematical_constants::getFloatingInteger< ScalarType >( 1 ) +
                    dotProductPositionAndEccentricityVectors ) < tolerance )
    {
        dotProductPositionAndEccentricityVectors =
                -mathematical_constants::getFloatingInteger< ScalarType >( 1 );
    }

    if ( std::fabs( dotProductPositionAndEccentricityVectors ) < tolerance )
    {
        dotProductPositionAndEccentricityVectors  =
                mathematical_constants::getFloatingInteger< ScalarType >( 0 );
    }

    // Compute and store true anomaly.
    computedKeplerianElements_( trueAnomalyIndex )
            = std::acos( dotProductPositionAndEccentricityVectors );

    // Check if the quandrant is correct.
    if ( trueAnomalyQuandrantCondition <
         mathematical_constants::getFloatingInteger< ScalarType >( 0 ) )
    {
        computedKeplerianElements_( trueAnomalyIndex ) =
                mathematical_constants::getFloatingInteger< ScalarType >( 2 ) *
                mathematical_constants::getPi< ScalarType >( ) -
                computedKeplerianElements_( trueAnomalyIndex );
    }

    // Return converted Keplerian elements.
    return computedKeplerianElements_;
}



//! Convert true anomaly to (elliptical) eccentric anomaly.
/*!
 * Converts true anomaly to eccentric anomaly for elliptical orbits ( 0 <= eccentricity < 1.0 ).
 * The equations used can be found in (Chobotov, 2002).
 * \param trueAnomaly True anomaly.                                                           [rad]
 * \param eccentricity Eccentricity.                                                            [-]
 * \return (Elliptical) Eccentric anomaly.                                                    [rad]
 */
template< typename ScalarType = double >
ScalarType convertTrueAnomalyToEllipticalEccentricAnomaly(
        const ScalarType trueAnomaly, const ScalarType eccentricity )

{
    using std::cos;
    using std::sqrt;
    using std::atan2;

    if ( eccentricity >= mathematical_constants::getFloatingInteger< ScalarType >( 1 ) ||
         eccentricity < mathematical_constants::getFloatingInteger< ScalarType >( 0 ) )
    {
        boost::throw_exception(
                    boost::enable_error_info(
                        std::runtime_error( "Eccentricity is invalid." ) ) );
    }
    else
    {
        // Declare and compute sine and cosine of eccentric anomaly.
        ScalarType sineOfEccentricAnomaly_ =
                sqrt( mathematical_constants::getFloatingInteger< ScalarType >( 1 ) -
                      eccentricity * eccentricity ) * std::sin( trueAnomaly ) /
                ( mathematical_constants::getFloatingInteger< ScalarType >( 1 ) +
                  eccentricity * cos( trueAnomaly ) );
        ScalarType cosineOfEccentricAnomaly_ = ( eccentricity + cos( trueAnomaly ) )
                / ( mathematical_constants::getFloatingInteger< ScalarType >( 1 ) +
                    eccentricity * cos( trueAnomaly ) );

        // Return elliptic eccentric anomaly.
        return atan2( sineOfEccentricAnomaly_, cosineOfEccentricAnomaly_ );
    }
}

//! Convert true anomaly to hyperbolic eccentric anomaly.
/*!
 * Converts true anomaly to hyperbolic eccentric anomaly for hyperbolic orbits
 * ( eccentricity > 1.0 ). The equations used can be found in (Chobotov, 2002).
 * \param trueAnomaly True anomaly.                                                           [rad]
 * \param eccentricity Eccentricity.                                                            [-]
 * \return Hyperbolic eccentric anomaly.                                                      [rad]
 */
template< typename ScalarType = double >
ScalarType convertTrueAnomalyToHyperbolicEccentricAnomaly( const ScalarType trueAnomaly,
                                                           const ScalarType eccentricity )
{
    if ( eccentricity <= mathematical_constants::getFloatingInteger< ScalarType >( 1 ) )
    {
        boost::throw_exception(
                    boost::enable_error_info(
                        std::runtime_error( "Eccentricity is invalid." ) ) );
    }

    else
    {
        using std::cos;

        // Compute hyperbolic sine and hyperbolic cosine of hyperbolic eccentric anomaly.
        ScalarType hyperbolicSineOfHyperbolicEccentricAnomaly_
                = std::sqrt( eccentricity * eccentricity -
                             mathematical_constants::getFloatingInteger< ScalarType >( 1 ) )
                * std::sin( trueAnomaly ) /
                ( mathematical_constants::getFloatingInteger< ScalarType >( 1 ) +
                  cos( trueAnomaly ) );

        ScalarType hyperbolicCosineOfHyperbolicEccentricAnomaly_
                = ( cos( trueAnomaly ) + eccentricity ) /
                ( mathematical_constants::getFloatingInteger< ScalarType >( 1 ) +
                  cos( trueAnomaly ) );

        // Return hyperbolic eccentric anomaly.
        return boost::math::atanh( hyperbolicSineOfHyperbolicEccentricAnomaly_
                                   / hyperbolicCosineOfHyperbolicEccentricAnomaly_ );
    }
}

//! Convert true anomaly to eccentric anomaly.
/*!
 * Converts true anomaly to eccentric anomaly for elliptical and hyperbolic orbits
 * ( eccentricity < 1.0 && eccentricity > 1.0 ). This function is essentially a wrapper for
 * convertTrueAnomalyToEllipticalEccentricAnomaly() and
 * convertTrueAnomalyToHyperbolicEccentricAnomaly(). It should be used in cases where the
 * eccentricity of the orbit is not known a priori. Currently, this implementation performs a
 * check on the eccentricity and throws an error for eccentricity < 0.0 and parabolic orbits, which
 * have not been implemented. The equations used can be found in (Chobotov, 2002).
 * \param trueAnomaly True anomaly.                                                           [rad]
 * \param eccentricity Eccentricity.                                                            [-]
 * \return Eccentric anomaly.                                                                 [rad]
 */
template< typename ScalarType = double >
ScalarType convertTrueAnomalyToEccentricAnomaly( const ScalarType trueAnomaly,
                                                 const ScalarType eccentricity )
{
    // Declare computed eccentric anomaly.
    ScalarType eccentricAnomaly_ = 0.0;

    // Check if eccentricity is invalid and throw an error if true.
    if ( eccentricity < mathematical_constants::getFloatingInteger< ScalarType >( 0 ) )
    {
        boost::throw_exception(
                    boost::enable_error_info(
                        std::runtime_error( "Eccentricity is invalid." ) ) );
    }

    // Check if orbit is parabolic and throw an error if true.
    else if ( std::fabs( eccentricity -
                         mathematical_constants::getFloatingInteger< ScalarType >( 1 ) ) <
              std::numeric_limits< ScalarType >::epsilon( ) )
    {
        boost::throw_exception(
                    boost::enable_error_info(
                        std::runtime_error(
                            "Parabolic orbits have not yet been implemented." ) ) );
    }

    // Check if orbit is elliptical and compute eccentric anomaly.
    else if ( eccentricity >= mathematical_constants::getFloatingInteger< ScalarType >( 0 ) &&
              eccentricity < mathematical_constants::getFloatingInteger< ScalarType >( 1 ) )
    {
        eccentricAnomaly_ = convertTrueAnomalyToEllipticalEccentricAnomaly< ScalarType >(
                    trueAnomaly, eccentricity );
    }

    else if ( eccentricity > mathematical_constants::getFloatingInteger< ScalarType >( 1 ) )
    {
        eccentricAnomaly_ = convertTrueAnomalyToHyperbolicEccentricAnomaly< ScalarType >(
                    trueAnomaly, eccentricity );
    }

    // Return computed eccentric anomaly.
    return eccentricAnomaly_;
}

//! Convert (elliptical) eccentric anomaly to true anomaly.
/*!
 * Converts eccentric anomaly to true anomaly for elliptical orbits ( 0 <= eccentricity < 1.0 ).
 * The equations used can be found in (Chobotov, 2002).
 * \param ellipticEccentricAnomaly Elliptical eccentric anomaly.                              [rad]
 * \param eccentricity Eccentricity.                                                            [-]
 * \return True anomaly.                                                                      [rad]
 */
template< typename ScalarType = double >
ScalarType convertEllipticalEccentricAnomalyToTrueAnomaly(
        const ScalarType ellipticEccentricAnomaly,
        const ScalarType eccentricity )
{
    if ( eccentricity >= mathematical_constants::getFloatingInteger< ScalarType >( 1 ) ||
         eccentricity < mathematical_constants::getFloatingInteger< ScalarType >( 0 ) )
    {
        boost::throw_exception(
                    boost::enable_error_info(
                        std::runtime_error( "Eccentricity is invalid." ) ) );
    }

    else
    {
        using std::cos;
        using std::sqrt;

        // Compute sine and cosine of true anomaly.
        ScalarType sineOfTrueAnomaly_ =
                sqrt( mathematical_constants::getFloatingInteger< ScalarType >( 1 ) -
                      eccentricity * eccentricity ) *
                std::sin( ellipticEccentricAnomaly )
                / ( mathematical_constants::getFloatingInteger< ScalarType >( 1 ) -
                    eccentricity * cos( ellipticEccentricAnomaly ) );

        ScalarType cosineOfTrueAnomaly_ = ( cos( ellipticEccentricAnomaly ) - eccentricity )
                / ( mathematical_constants::getFloatingInteger< ScalarType >( 1 ) -
                    eccentricity * cos( ellipticEccentricAnomaly ) );

        // Return true anomaly.
        return std::atan2( sineOfTrueAnomaly_, cosineOfTrueAnomaly_  );
    }
}

//! Convert hyperbolic eccentric anomaly to true anomaly.
/*!
 * Converts hyperbolic eccentric anomaly to true anomaly for hyperbolic orbits
 * ( eccentricity > 1.0 ). The equations used can be found in (Chobotov, 2002).
 * \param hyperbolicEccentricAnomaly Hyperbolic eccentric anomaly.                            [rad]
 * \param eccentricity Eccentricity.                                                            [-]
 * \return True anomaly.                                                                      [rad]
 */
template< typename ScalarType = double >
ScalarType convertHyperbolicEccentricAnomalyToTrueAnomaly(
        const ScalarType hyperbolicEccentricAnomaly,
        const ScalarType eccentricity )
{
    if ( eccentricity <= mathematical_constants::getFloatingInteger< ScalarType >( 1 ) )
    {
        boost::throw_exception(
                    boost::enable_error_info(
                        std::runtime_error( "Eccentricity is invalid." ) ) );
    }

    else
    {
        using std::cosh;

        // Compute sine and cosine of true anomaly.
        ScalarType sineOfTrueAnomaly_
                = std::sqrt( eccentricity * eccentricity -
                             mathematical_constants::getFloatingInteger< ScalarType >( 1 ) )
                * std::sinh( hyperbolicEccentricAnomaly )
                / ( eccentricity * cosh( hyperbolicEccentricAnomaly ) -
                    mathematical_constants::getFloatingInteger< ScalarType >( 1 ) );

        ScalarType cosineOfTrueAnomaly_
                = ( eccentricity - cosh( hyperbolicEccentricAnomaly ) )
                / ( eccentricity * cosh( hyperbolicEccentricAnomaly ) -
                    mathematical_constants::getFloatingInteger< ScalarType >( 1 ) );

        // Return true anomaly.
        return std::atan2( sineOfTrueAnomaly_, cosineOfTrueAnomaly_ );
    }

}


//! Convert eccentric anomaly to true anomaly.
/*!
 * Converts eccentric anomaly to true anomaly for elliptical and hyperbolic orbits
 * ( eccentricity < 1.0 && eccentricity > 1.0 ). This function is essentially a wrapper for
 * convertEllipticalEccentricAnomalyToTrueAnomaly() and
 * convertHyperbolicEccentricAnomalyToTrueAnomaly(). It should be used in cases where the
 * eccentricity of the orbit is not known a priori. Currently, this implementation performs a
 * check on the eccentricity and throws an error for eccentricity < 0.0 and parabolic orbits, which
 * have not been implemented. The equations used can be found in (Chobotov, 2002).
 * \param eccentricAnomaly Eccentric anomaly.                                                 [rad]
 * \param eccentricity Eccentricity.                                                            [-]
 * \return True anomaly.                                                                      [rad]
 */
template< typename ScalarType = double >
ScalarType convertEccentricAnomalyToTrueAnomaly( const ScalarType eccentricAnomaly,
                                                 const ScalarType eccentricity )
{
    // Declare computed true anomaly.
    ScalarType trueAnomaly_ = -mathematical_constants::getFloatingInteger< ScalarType >( 0 );

    // Check if eccentricity is invalid and throw an error if true.
    if ( eccentricity < mathematical_constants::getFloatingInteger< ScalarType >( 0 ) )
    {
        boost::throw_exception(
                    boost::enable_error_info(
                        std::runtime_error( "Eccentricity is invalid." ) ) );
    }

    // Check if orbit is parabolic and throw an error if true.
    else if ( std::fabs( eccentricity -
                         mathematical_constants::getFloatingInteger< ScalarType >( 1 ) ) <
              std::numeric_limits< ScalarType >::epsilon( ) )
    {
        boost::throw_exception(
                    boost::enable_error_info(
                        std::runtime_error(
                            "Parabolic orbits have not yet been implemented." ) ) );
    }

    // Check if orbit is elliptical and compute true anomaly.
    else if ( eccentricity >= mathematical_constants::getFloatingInteger< ScalarType >( 0 ) &&
              eccentricity < mathematical_constants::getFloatingInteger< ScalarType >( 1 ) )
    {
        trueAnomaly_ = convertEllipticalEccentricAnomalyToTrueAnomaly( eccentricAnomaly,
                                                                       eccentricity );
    }

    else if ( eccentricity > mathematical_constants::getFloatingInteger< ScalarType >( 1 ) )
    {
        trueAnomaly_ = convertHyperbolicEccentricAnomalyToTrueAnomaly( eccentricAnomaly,
                                                                       eccentricity );
    }

    // Return computed true anomaly.
    return trueAnomaly_;
}


//! Convert (elliptical) eccentric anomaly to mean anomaly.
/*!
 * Converts eccentric anomaly to mean anomaly for elliptical orbits ( 0 <= eccentricity < 1.0 ).
 * The equations used can be found in (Chobotov, 2002).
 * \param eccentricity Eccentricity.                                                            [-]
 * \param ellipticalEccentricAnomaly (Elliptical) eccentric anomaly [rad].
 * \return Mean anomaly [rad].
 */
template< typename ScalarType = double >
ScalarType convertEllipticalEccentricAnomalyToMeanAnomaly(
        const ScalarType ellipticalEccentricAnomaly,
        const ScalarType eccentricity )
{
    return ellipticalEccentricAnomaly - eccentricity * std::sin( ellipticalEccentricAnomaly );
}


//! Convert hyperbolic eccentric anomaly to mean anomaly.
/*!
 * Converts hyperbolic eccentric anomaly to mean anomaly for hyperbolic orbits
 * ( eccentricity > 1.0 ). The equations used can be found in (Chobotov, 2002).
 * \param hyperbolicEccentricAnomaly Hyperbolic eccentric anomaly.                            [rad]
 * \param eccentricity Eccentricity.                                                            [-]
 * \return Mean anomaly.                                                                      [rad]
 */
template< typename ScalarType = double >
ScalarType convertHyperbolicEccentricAnomalyToMeanAnomaly(
        const ScalarType hyperbolicEccentricAnomaly,
        const ScalarType eccentricity )
{
    return eccentricity * std::sinh( hyperbolicEccentricAnomaly ) - hyperbolicEccentricAnomaly;
}

//! Convert eccentric anomaly to mean anomaly.
/*!
 * Converts eccentric anomaly to mean anomaly for elliptical and hyperbolic orbits
 * ( eccentricity < 1.0 && eccentricity > 1.0 ). This function is essentially a wrapper for
 * convertEllipticalEccentricAnomalyToMeanAnomaly() and
 * convertHyperbolicEccentricAnomalyToMeanAnomaly(). It should be used in cases where the
 * eccentricity of the orbit is not known a priori. Currently, this implementation performs a
 * check on the eccentricity and throws an error for eccentricity < 0.0 and parabolic orbits, which
 * have not been implemented. The equations used can be found in (Chobotov, 2002).
 * \param eccentricity Eccentricity.                                                            [-]
 * \param eccentricAnomaly Eccentric anomaly.                                                 [rad]
 * \return Mean anomaly.                                                                      [rad]
 */
template< typename ScalarType = double >
ScalarType convertEccentricAnomalyToMeanAnomaly(
        const ScalarType eccentricAnomaly,
        const ScalarType eccentricity )
{
    // Declare computed mean anomaly.
    ScalarType meanAnomaly_ = 0.0;

    // Check if eccentricity is invalid and throw an error if true.
    if ( eccentricity < mathematical_constants::getFloatingInteger< ScalarType >( 0 ) )
    {
        boost::throw_exception(
                    boost::enable_error_info(
                        std::runtime_error( "Eccentricity is invalid." ) ) );
    }

    // Check if orbit is parabolic and throw an error if true.
    else if ( std::fabs( eccentricity -
                         mathematical_constants::getFloatingInteger< ScalarType >( 1 ) ) <
              std::numeric_limits< ScalarType >::epsilon( ) )
    {
        boost::throw_exception(
                    boost::enable_error_info(
                        std::runtime_error(
                            "Parabolic orbits have not yet been implemented." ) ) );
    }

    // Check if orbit is elliptical and compute true anomaly.
    else if ( eccentricity >=
              mathematical_constants::getFloatingInteger< ScalarType >( 0 ) &&
              eccentricity < mathematical_constants::getFloatingInteger< ScalarType >( 1 ) )
    {
        meanAnomaly_ = convertEllipticalEccentricAnomalyToMeanAnomaly< ScalarType >(
                    eccentricAnomaly, eccentricity );
    }

    else if ( eccentricity > mathematical_constants::getFloatingInteger< ScalarType >( 1 ) )
    {
        meanAnomaly_ = convertHyperbolicEccentricAnomalyToMeanAnomaly< ScalarType >(
                    eccentricAnomaly, eccentricity );
    }

    // Return computed mean anomaly.
    return meanAnomaly_;
}

//! Convert elapsed time to (elliptical) mean anomaly change.
/*!
 * Converts elapsed time to mean anomaly change for elliptical orbits ( 0 <= eccentricity < 1.0 ).
 * The semi-major axis must be non-negative; this function will throw an error to indicate if the
 * semi-major axis is invalid. The equation used can be found in (Chobotov, 2002).
 * \param elapsedTime Elapsed time.                                                             [s]
 * \param centralBodyGravitationalParameter Gravitational parameter of central body.      [m^3/s^2]
 * \param semiMajorAxis Semi-major axis.                                                        [m]
 * \return (Elliptical) Mean anomaly change.                                                  [rad]
 */
template< typename ScalarType = double >
ScalarType convertElapsedTimeToEllipticalMeanAnomalyChange(
        const ScalarType elapsedTime, const ScalarType centralBodyGravitationalParameter,
        const ScalarType semiMajorAxis )
{
    // Check if semi-major axis is invalid and throw error if true.
    if ( semiMajorAxis < mathematical_constants::getFloatingInteger< ScalarType >( 0 ) )
    {
        boost::throw_exception(
                    boost::enable_error_info(
                        std::runtime_error( "Semi-major axis is invalid." ) ) );
    }

    // Else return elliptical mean anomaly change.
    else
    {
        return std::sqrt( centralBodyGravitationalParameter
                          / ( semiMajorAxis * semiMajorAxis * semiMajorAxis ) ) * elapsedTime;
    }
}


//! Convert elapsed time to mean anomaly change for hyperbolic orbits.
/*!
 * Converts elapsed time to mean anomaly change for hyperbolic orbits ( eccentricity > 1.0 ).
 * The semi-major axis must be non-positive; this function will throw an error to indicate if the
 * semi-major axis is invalid. The equation used can be found in (Chobotov, 2002).
 * \param elapsedTime Elapsed time.                                                             [s]
 * \param centralBodyGravitationalParameter Gravitational parameter of central body.      [m^3/s^2]
 * \param semiMajorAxis Semi-major axis.                                                        [m]
 * \return Mean anomaly change.                                                               [rad]
 */
template< typename ScalarType = double >
ScalarType convertElapsedTimeToHyperbolicMeanAnomalyChange(
        const ScalarType elapsedTime, const ScalarType centralBodyGravitationalParameter,
        const ScalarType semiMajorAxis )
{
    // Check if semi-major axis is invalid and throw error if true.
    if ( semiMajorAxis > mathematical_constants::getFloatingInteger< ScalarType >( 0 ) )
    {
        boost::throw_exception(
                    boost::enable_error_info(
                        std::runtime_error( "Semi-major axis is invalid." ) ) );
    }

    // Else return hyperbolic mean anomaly change.
    else
    {
        return std::sqrt( centralBodyGravitationalParameter
                          / ( - semiMajorAxis * semiMajorAxis * semiMajorAxis ) ) * elapsedTime;
    }
}

//! Convert elapsed time to mean anomaly change.
/*!
 * Converts elapsed time to mean anomaly change for elliptical and hyperbolic orbits
 * ( eccentricity < 1.0 && eccentricity > 1.0 ). This function is essentially a wrapper for
 * convertElapsedTimeToEllipticalMeanAnomalyChange() and
 * convertElapsedTimeToHyperbolicMeanAnomalyChange(). It should be used in cases where the
 * eccentricity of the orbit is not known a priori. The equations used can be found in
 * (Wertz, 2001).
 * \param elapsedTime Elapsed time.                                                             [s]
 * \param centralBodyGravitationalParameter Gravitational parameter of central body.      [m^3/s^2]
 * \param semiMajorAxis Semi-major axis.                                                        [m]
 * \return Mean anomaly change.                                                               [rad]
 */
template< typename ScalarType = double >
ScalarType convertElapsedTimeToMeanAnomalyChange(
        const ScalarType elapsedTime, const ScalarType centralBodyGravitationalParameter,
        const ScalarType semiMajorAxis )
{
    // Declare computed mean anomaly change.
    ScalarType meanAnomalyChange_ = -mathematical_constants::getFloatingInteger< ScalarType >( 0 );

    // Check if orbit is elliptical and compute mean anomaly change.
    if ( semiMajorAxis > mathematical_constants::getFloatingInteger< ScalarType >( 0 ) )
    {
        meanAnomalyChange_ = convertElapsedTimeToEllipticalMeanAnomalyChange(
                    elapsedTime, centralBodyGravitationalParameter, semiMajorAxis );
    }

    // Else orbit is hyperbolic; compute mean anomaly change.
    else if ( semiMajorAxis < mathematical_constants::getFloatingInteger< ScalarType >( 0 ) )
    {
        meanAnomalyChange_ = convertElapsedTimeToHyperbolicMeanAnomalyChange(
                    elapsedTime, centralBodyGravitationalParameter, semiMajorAxis );
    }

    // Return computed mean anomaly change.
    return meanAnomalyChange_;
}


//! Convert (elliptical) mean anomaly change to elapsed time.
/*!
 * Converts mean anomaly change to elapsed time for elliptical orbits ( 0 <= eccentricity < 1.0 ).
 * The equation used can be found in (Wertz, 2001). This function checks if the semi-major axis is
 * non-negative and throws an error if it not.
 * \param ellipticalMeanAnomalyChange (Elliptical) Mean anomaly change.                       [rad]
 * \param centralBodyGravitationalParameter Gravitational parameter of central body.      [m^3/s^2]
 * \param semiMajorAxis Semi-major axis.                                                        [m]
 * \return Elapsed time.                                                                        [s]
 */
template< typename ScalarType = double >
ScalarType convertEllipticalMeanAnomalyChangeToElapsedTime(
        const ScalarType ellipticalMeanAnomalyChange,
        const ScalarType centralBodyGravitationalParameter,
        const ScalarType semiMajorAxis )
{
    // Check if semi-major axis is invalid and throw error if true.
    if ( semiMajorAxis < mathematical_constants::getFloatingInteger< ScalarType >( 0 ) )
    {
        boost::throw_exception(
                    boost::enable_error_info(
                        std::runtime_error( "Semi-major axis is invalid." ) ) );
    }

    // Else return elapsed time.
    else
    {
        return ellipticalMeanAnomalyChange * std::sqrt(
                    semiMajorAxis * semiMajorAxis * semiMajorAxis
                    / centralBodyGravitationalParameter );
    }
}

//! Convert hyperbolic mean anomaly change to elapsed time.
/*!
 * Converts mean anomaly change to elapsed time for hyperbolic orbits ( eccentricity > 1.0 ).
 * The equation used can be found in (Wertz, 2001). This function checks if the semi-major axis is
 * non-positive and throws an error if it not.
 * \param hyperbolicMeanAnomalyChange Hyperbolic mean anomaly change.                         [rad]
 * \param centralBodyGravitationalParameter Gravitational parameter of central body.      [m^3/s^2]
 * \param semiMajorAxis Semi-major axis.                                                        [m]
 * \return Elapsed time.                                                                        [s]
 */
template< typename ScalarType = double >
ScalarType convertHyperbolicMeanAnomalyChangeToElapsedTime(
        const ScalarType hyperbolicMeanAnomalyChange,
        const ScalarType centralBodyGravitationalParameter,
        const ScalarType semiMajorAxis )
{
    // Check if semi-major axis is invalid and throw error if true.
    if ( semiMajorAxis > mathematical_constants::getFloatingInteger< ScalarType >( 0 ) )
    {
        boost::throw_exception(
                    boost::enable_error_info(
                        std::runtime_error( "Semi-major axis is invalid." ) ) );
    }

    // Else return elapsed time.
    else
    {
        return std::sqrt( -semiMajorAxis * semiMajorAxis * semiMajorAxis
                          / centralBodyGravitationalParameter ) * hyperbolicMeanAnomalyChange;
    }
}

//! Convert mean anomaly change to elapsed time.
/*!
 * Converts mean anomaly change to elapsed time for elliptical and hyperbolic orbits
 * ( eccentricity < 1.0 && eccentricity > 1.0 ). This function is essentially a wrapper for
 * convertEllipticalMeanAnomalyChangeToElapsedTime() and
 * convertHyperbolicMeanAnomalyChangeToElapsedTime(). It should be used in cases where the
 * eccentricity of the orbit is not known a priori. The equations used can be found in
 * (Wertz, 2001).
 * \param meanAnomalyChange Mean anomaly change.                                              [rad]
 * \param centralBodyGravitationalParameter Gravitational parameter of central body.      [m^3/s^2]
 * \param semiMajorAxis Semi-major axis.                                                        [m]
 * \return Elapsed time.                                                                        [s]
 */
template< typename ScalarType = double >
ScalarType convertMeanAnomalyChangeToElapsedTime(
        const ScalarType meanAnomalyChange, const ScalarType centralBodyGravitationalParameter,
        const ScalarType semiMajorAxis )
{
    // Declare computed elapsed time.
    ScalarType elapsedTime_ = mathematical_constants::getFloatingInteger< ScalarType >( 0 );

    // Check if orbit is elliptical and compute elapsed time.
    if ( semiMajorAxis > mathematical_constants::getFloatingInteger< ScalarType >( 0 ) )
    {
        elapsedTime_ = convertEllipticalMeanAnomalyChangeToElapsedTime< ScalarType >(
                    meanAnomalyChange, centralBodyGravitationalParameter, semiMajorAxis );
    }

    // Else orbit is hyperbolic; compute elapsed time.
    else if ( semiMajorAxis < mathematical_constants::getFloatingInteger< ScalarType >( 0 ) )
    {
        elapsedTime_ = convertHyperbolicMeanAnomalyChangeToElapsedTime< ScalarType >(
                    meanAnomalyChange, centralBodyGravitationalParameter, semiMajorAxis );
    }

    // Return computed elapsed time.
    return elapsedTime_;
}

//! Convert (elliptical) mean motion to semi-major axis.
/*!
 * Converts mean motion to semi-major axis for elliptical orbits.
 * \param ellipticalMeanMotion (Elliptical) Mean motion.                                    [rad/s]
 * \param centralBodyGravitationalParameter Gravitational parameter of central body.      [m^3/s^2]
 * \return semiMajorAxis Semi-major axis.                                                       [m]
 */
template< typename ScalarType = double >
ScalarType convertEllipticalMeanMotionToSemiMajorAxis(
        const ScalarType ellipticalMeanMotion, const ScalarType centralBodyGravitationalParameter )
{
    return std::pow( centralBodyGravitationalParameter
                     / ( ellipticalMeanMotion * ellipticalMeanMotion ),
                     mathematical_constants::getFloatingFraction< ScalarType >( 1, 3 ) );
}

//! Convert semi-major axis to elliptical mean motion.
/*!
 * Converts semi-major axis to elliptical mean motion.
 * \param semiMajorAxis Semi-major axis.                                                        [m]
 * \param centralBodyGravitationalParameter Gravitational parameter of central body.      [m^3/s^2]
 * \return ellipticalMeanMotion (Elliptical) Mean motion.                                   [rad/s]
 */
template< typename ScalarType = double >
ScalarType convertSemiMajorAxisToEllipticalMeanMotion(
        const ScalarType semiMajorAxis, const ScalarType centralBodyGravitationalParameter )
{
    // Check if semi-major axis is invalid and throw error if true.
    if ( semiMajorAxis < mathematical_constants::getFloatingInteger< ScalarType >( 0 ) )
    {
        boost::throw_exception(
                    boost::enable_error_info(
                        std::runtime_error( "Semi-major axis is invalid." ) ) );
    }

    // Else compute and return elliptical mean motion.
    {
        return std::sqrt( centralBodyGravitationalParameter /
                          ( semiMajorAxis * semiMajorAxis * semiMajorAxis ) );
    }
}




///////////////////////////////////////////////////////////////////////////////////
//////////////////////////// EQUINOCTIAL ELEMENTS /////////////////////////////////
///////////////////////////////////////////////////////////////////////////////////

//! Get eccentric longitude for set of equinoctial elements
/*!
 * Get eccentric longitude for set of equinoctial elements.
 * [Semianalytic satellite theory, Danielson (1995), Section 7.1]
 * \param equinoctialElements Set of equinoctial elements.
 * \return Eccentric longitude.                                                    [rad]
 */
template< typename ScalarType = double >
ScalarType getEccentricLongitude( const Eigen::Matrix< ScalarType, 6, 1 > equinoctialElements )
{
    using namespace tudat::mathematical_constants;
    using namespace tudat::basic_mathematics;

    /// Solve equinoctial form of Kepler's Equation [ Section 7.1 ]

    // Get necessary components
    const ScalarType h = equinoctialElements[ hIndex ];
    const ScalarType k = equinoctialElements[ kIndex ];
    const ScalarType lambda = equinoctialElements[ meanLongitudeIndex ];

    // Initialize eccentric longitude (first guess)
    ScalarType F = lambda;

    // Iterate
    ScalarType change = 1;
    const ScalarType tolerance = 1e-12;
    while ( std::abs( change ) > tolerance )
    {
        // Eq. 7.1-(2)
        const ScalarType sinF = std::sin( F );
        const ScalarType cosF = std::cos( F );
        change = ( F + h * cosF - k * sinF - lambda ) / ( 1 - h * sinF - k * cosF );
        F -= change;
    }

    return computeModulo( F, 2*PI );
}


//! Get true longitude for set of equinoctial elements
/*!
 * Get true longitude for set of equinoctial elements.
 * [Semianalytic satellite theory, Danielson (1995), Section 2.1.4]
 * \param equinoctialElements Set of equinoctial elements.
 * \return True longitude.                                                    [rad]
 */
template< typename ScalarType = double >
ScalarType getTrueLongitude( const Eigen::Matrix< ScalarType, 6, 1 > equinoctialElements )
{
    using namespace tudat::mathematical_constants;
    using namespace tudat::basic_mathematics;

    // Get eccentric longitude and sin, cos components
    const ScalarType F = getEccentricLongitude( equinoctialElements );
    const ScalarType sinF = std::sin( F );
    const ScalarType cosF = std::cos( F );

    // Get necessary components
    const ScalarType h = equinoctialElements[ hIndex ];
    const ScalarType k = equinoctialElements[ kIndex ];

    // Pre-compute repeated terms
    const ScalarType h2 = pow( h, 2 );
    const ScalarType k2 = pow( k, 2 );
    const ScalarType b = 1 / ( 1 + std::sqrt( 1 - h2 - k2 ) );
    const ScalarType hkb = h * k * b;
    const ScalarType denom = 1 - h * sinF - k * cosF;

    // Eq. 2.1.4-(5)
    const ScalarType sinL = ( (1 - k2 * b) * sinF + hkb * cosF - h ) / denom;
    const ScalarType cosL = ( (1 - h2 * b) * cosF + hkb * sinF - k ) / denom;

    return computeModulo( std::atan2( sinL, cosL ), 2*PI );
}


//! Get the basis vectors of the equinoctial reference frame
/*!
 * \param p The p component of the equinoctial elements.
 * \param p The q component of the equinoctial elements.
 * \param retrogradeElements Whether to use the retrograde set of elements (to avoid singularities for orbits
 * with inclinations of 180 degrees). Default is false, which avoids singularities for equatorial orbits.
 * \param computeF Whether to compute the f vector. If false, the first column will be [0; 0; 0].
 * \param computeG Whether to compute the g vector. If false, the second column will be [0; 0; 0].
 * \param computeW Whether to compute the w vector. If false, the third column will be [0; 0; 0].
 * \return Matrix with three columns vectors: [ f g w ] with f = [ fx; fy; fz ], etc.
 */
template< typename ScalarType = double >
Eigen::Matrix< ScalarType, 3, 3 > getEquinoctialReferenceFrameBasisVectors( const ScalarType p, const ScalarType q,
                                                                            const bool retrogradeElements = false,
                                                                            const bool computeF = true,
                                                                            const bool computeG = true,
                                                                            const bool computeW = true )
{
    Eigen::Matrix< ScalarType, 3, 3 > frameVectors;
    const ScalarType I = retrogradeElements ? -1 : 1;

    // Compute f, g and w vectors of the equinoctial reference frame [ Eq. 2.1.4-(1) ]
    const ScalarType p2 = pow( p, 2 );
    const ScalarType q2 = pow( q, 2 );
    const ScalarType oneo1pp2pq2 = 1 / ( 1 + p2 + q2 );
    const ScalarType p2mq2 = p2 - q2;
    const ScalarType twopq = 2 * p * q;

    if ( computeF ) {
        frameVectors.col( 0 ) << oneo1pp2pq2 * ( 1 - p2mq2 ),
                                 oneo1pp2pq2 * twopq,
                                 oneo1pp2pq2 * -2 * I * p;
    }

    if ( computeG ) {
        frameVectors.col( 1 ) << oneo1pp2pq2 * twopq * I,
                                 oneo1pp2pq2 * ( 1 + p2mq2 ) * I,
                                 oneo1pp2pq2 * 2 * q;
    }

    if ( computeW ) {
        frameVectors.col( 2 ) << oneo1pp2pq2 * 2 * p,
                                 oneo1pp2pq2 * -2 * q,
                                 oneo1pp2pq2 * ( 1 - p2 - q2 ) * I;
    }

    return frameVectors;
}


//! Convert Cartesian to equinoctial orbital elements.
/*!
 * Converts Cartesian to equinoctial orbital elements.
 * [Semianalytic satellite theory, Danielson (1995), Section 2.1.5]
 * \param cartesianElements Vector containing Cartesian elements. Order of elements is important!
 *          cartesianElements( 0 ) = x-position coordinate,                                     [m]
 *          cartesianElements( 1 ) = y-position coordinate,                                     [m]
 *          cartesianElements( 2 ) = z-position coordinate,                                     [m]
 *          cartesianElements( 3 ) = x-velocity coordinate,                                   [m/s]
 *          cartesianElements( 4 ) = y-velocity coordinate,                                   [m/s]
 *          cartesianElements( 5 ) = z-velocity coordinate.                                   [m/s]
 * \param centralBodyGravitationalParameter Gravitational parameter of central body.
 * \param retrogradeElements Whether to use the retrograde set of elements (to avoid singularities for orbits
 * with inclinations of 180 degrees). Default is false, which avoids singularities for equatorial orbits.
 * \return Converted state in equinoctial elements. The order of elements is fixed!
 *          equinoctialElements( 0 ) = semi-major axis,                                         [m]
 *          equinoctialElements( 1 ) = h component,                                             [-]
 *          equinoctialElements( 2 ) = k component,                                             [-]
 *          equinoctialElements( 3 ) = p component,                                             [-]
 *          equinoctialElements( 4 ) = q component,                                             [-]
 *          equinoctialElements( 5 ) = mean longitude.                                        [rad]
 */
template< typename ScalarType = double >
Eigen::Matrix< ScalarType, 6, 1 > convertCartesianToEquinoctialElements(
        const Eigen::Matrix< ScalarType, 6, 1 >& cartesianElements,
        const ScalarType centralBodyGravitationalParameter, const bool retrogradeElements = false )
{
    using std::pow;

    // Set tolerance.
    const ScalarType tolerance = 20.0 * std::numeric_limits< ScalarType >::epsilon( );

    // Determine the retrograde factor
    const ScalarType I = retrogradeElements ? -1 : 1;

    // Declare converted equinoctial elements.
    Eigen::Matrix< ScalarType, 6, 1 > equinoctialElements;

    // Set position and velocity vectors.
    const Eigen::Matrix< ScalarType, 3, 1 > r = cartesianElements.segment( 0, 3 );
    const Eigen::Matrix< ScalarType, 3, 1 > v = cartesianElements.segment( 3, 3 );

    // Compute semi-major axis [ Eq. 2.1.5-(1) ]
    const ScalarType r_norm = r.norm( );
    const ScalarType a = 1 / ( 2 / r_norm - v.squaredNorm( ) / centralBodyGravitationalParameter );

    // Compute r x v
    const Eigen::Matrix< ScalarType, 3, 1 > rxv = r.cross( v );

    // Compute w vector of the equinoctial reference frame [ Eq. 2.1.5-(2) ]
    const Eigen::Matrix< ScalarType, 3, 1 > w = rxv.normalized( );

    // Compute p and q components [ Eq. 2.1.5-(3) ]
    const ScalarType onepIwz = ( 1 + I * w[2] );
    const ScalarType p = w[0] / onepIwz;
    const ScalarType q = -w[1] / onepIwz;

    // Get f and g vectors of the equinoctial reference frame
    const Eigen::Matrix< ScalarType, 3, 3 > fg = getEquinoctialReferenceFrameBasisVectors( p, q, retrogradeElements,
                                                                                            true, true, false );
    const Eigen::Matrix< ScalarType, 3, 1 > f = fg.col( 0 );
    const Eigen::Matrix< ScalarType, 3, 1 > g = fg.col( 1 );

    // Compute eccentricity vector [ Eq. 2.1.5-(4) ]
    const Eigen::Matrix< ScalarType, 3, 1 > e = - r / r_norm + v.cross( rxv ) / centralBodyGravitationalParameter;

    // Compute h and k components [ Eq. 2.1.5-(5) ]
    const ScalarType h = e.dot( g );
    const ScalarType k = e.dot( f );

    // Compute X and Y components in the equinoctial reference frame [ Eq. 2.1.5-(6) ]
    const ScalarType X = r.dot( f );
    const ScalarType Y = r.dot( g );

    // Compute the eccentric longitude [ Eq. 2.1.5-(7) ]
    const ScalarType h2 = pow( h, 2 );
    const ScalarType k2 = pow( k, 2 );
    const ScalarType B = std::sqrt( 1 - h2 - k2 );
    const ScalarType b = 1 / ( 1 + B );
    const ScalarType hkb = h * k * b;
    const ScalarType sinF = h + ( ( 1 - h2 * b ) * Y - hkb * X ) / ( a * B );
    const ScalarType cosF = k + ( ( 1 - k2 * b ) * X - hkb * Y ) / ( a * B );
    const ScalarType F = std::atan2( sinF, cosF );

    // Compute the mean longitude [ Eq. 2.1.4-(2) ]
    const ScalarType lambda = F + h * cosF - k * sinF;

    // Store components
    equinoctialElements( semiMajorAxisIndex ) = a;
    equinoctialElements( hIndex )             = h;
    equinoctialElements( kIndex )             = k;
    equinoctialElements( pIndex )             = p;
    equinoctialElements( qIndex )             = q;
    equinoctialElements( meanLongitudeIndex ) = lambda;

    // Return converted equinoctial elements.
    return equinoctialElements;
}


//! Convert equinoctial to Cartesian orbital elements.
/*!
 * Converts equinoctial to Cartesian orbital elements.
 * [Semianalytic satellite theory, Danielson (1995), Section 2.1.4]
 * \param equinoctialElements Vector with equinoctial elements. Order of elements is important!
 *          equinoctialElements( 0 ) = semi-major axis,                                         [m]
 *          equinoctialElements( 1 ) = h component,                                             [-]
 *          equinoctialElements( 2 ) = k component,                                             [-]
 *          equinoctialElements( 3 ) = p component,                                             [-]
 *          equinoctialElements( 4 ) = q component,                                             [-]
 *          equinoctialElements( 5 ) = mean longitude.                                        [rad]
 * \param centralBodyGravitationalParameter Gravitational parameter of central body.
 * \param retrogradeElements Whether to use the retrograde set of elements (to avoid singularities for orbits
 * with inclinations of 180 degrees). Default is false, which avoids singularities for equatorial orbits.
 * \return Converted state in Cartesian elements. The order of elements is fixed!
 *          cartesianElements( 0 ) = x-position coordinate,                                     [m]
 *          cartesianElements( 1 ) = y-position coordinate,                                     [m]
 *          cartesianElements( 2 ) = z-position coordinate,                                     [m]
 *          cartesianElements( 3 ) = x-velocity coordinate,                                   [m/s]
 *          cartesianElements( 4 ) = y-velocity coordinate,                                   [m/s]
 *          cartesianElements( 5 ) = z-velocity coordinate.                                   [m/s]
 */
template< typename ScalarType = double >
Eigen::Matrix< ScalarType, 6, 1 > convertEquinoctialToCartesianElements(
        const Eigen::Matrix< ScalarType, 6, 1 >& equinoctialElements,
        const ScalarType centralBodyGravitationalParameter, const bool retrogradeElements = false )
{
    using std::pow;
    using std::sqrt;

    // Set tolerance.
    const ScalarType tolerance = 20.0 * std::numeric_limits< ScalarType >::epsilon( );

    // Declare converted Cartesian elements.
    Eigen::Matrix< ScalarType, 6, 1 > cartesianElements;

    // Read equinoctial elements
    const ScalarType a =      equinoctialElements( semiMajorAxisIndex );
    const ScalarType h =      equinoctialElements( hIndex );
    const ScalarType k =      equinoctialElements( kIndex );
    const ScalarType p =      equinoctialElements( pIndex );
    const ScalarType q =      equinoctialElements( qIndex );

    // Get f and g vectors of the equinoctial reference frame
    const Eigen::Matrix< ScalarType, 3, 3 > fg = getEquinoctialReferenceFrameBasisVectors( p, q, retrogradeElements,
                                                                                            true, true, false );
    const Eigen::Matrix< ScalarType, 3, 1 > f = fg.col( 0 );
    const Eigen::Matrix< ScalarType, 3, 1 > g = fg.col( 1 );

    // Get true longitude
    const ScalarType L = getTrueLongitude( equinoctialElements );

    // Compute repeated parameters
    const ScalarType h2 = pow( h, 2 );
    const ScalarType k2 = pow( k, 2 );
    const ScalarType B = sqrt( 1 - h2 - k2 );
    const ScalarType n = sqrt( centralBodyGravitationalParameter / pow( a, 3 ) );
    const ScalarType sinL = std::sin( L );
    const ScalarType cosL = std::cos( L );

    // Get radial distance [Eq. (6)]
    const ScalarType r = a * ( 1 - h2 - k2 ) / ( 1 + h * sinL + k * cosL );

    // Position in the equinoctial reference frame [Eq. (7)]
    const ScalarType X = r * cosL;
    const ScalarType Y = r * sinL;

    // Velocity in the equinoctial reference frame [Eq. (8)]
    const ScalarType Xdot = - n * a * ( h + sinL ) / B;
    const ScalarType Ydot = n * a * ( k + cosL ) / B;

    // Equinoctial reference frame --> Cartesian reference frame [Eq. (9)]
    const Eigen::Matrix< ScalarType, 3, 1 > position = X * f + Y * g;
    const Eigen::Matrix< ScalarType, 3, 1 > velocity = Xdot * f + Ydot * g;

    // Store Cartesian components
    cartesianElements[ xCartesianPositionIndex ] = position( 0 );
    cartesianElements[ yCartesianPositionIndex ] = position( 1 );
    cartesianElements[ zCartesianPositionIndex ] = position( 2 );
    cartesianElements[ xCartesianVelocityIndex ] = velocity( 0 );
    cartesianElements[ yCartesianVelocityIndex ] = velocity( 1 );
    cartesianElements[ zCartesianVelocityIndex ] = velocity( 2 );

    // Return converted Cartesian elements.
    return cartesianElements;
}


//! Convert Keplerian to equinoctial orbital elements.
/*!
 * Converts Keplerian to equinoctial orbital elements.
 * [Semianalytic satellite theory, Danielson (1995), Section 2.1.5]
 * \param keplerianElements Vector containing Keplerian elements. Order of elements is important!
 *          keplerianElements( 0 ) = semiMajorAxis,                                             [m]
 *          keplerianElements( 1 ) = eccentricity,                                              [-]
 *          keplerianElements( 2 ) = inclination,                                             [rad]
 *          keplerianElements( 3 ) = argument of periapsis,                                   [rad]
 *          keplerianElements( 4 ) = longitude of ascending node,                             [rad]
 *          keplerianElements( 5 ) = true anomaly.                                            [rad]
 * \param centralBodyGravitationalParameter Gravitational parameter of central body.
 * \param retrogradeElements Whether to use the retrograde set of elements (to avoid singularities for orbits
 * with inclinations of 180 degrees). Default is false, which avoids singularities for equatorial orbits.
 * \return Converted state in equinoctial elements. The order of elements is fixed!
 *          equinoctialElements( 0 ) = semi-major axis,                                         [m]
 *          equinoctialElements( 1 ) = h component,                                             [-]
 *          equinoctialElements( 2 ) = k component,                                             [-]
 *          equinoctialElements( 3 ) = p component,                                             [-]
 *          equinoctialElements( 4 ) = q component,                                             [-]
 *          equinoctialElements( 5 ) = mean longitude.                                        [rad]
 */
template< typename ScalarType = double >
Eigen::Matrix< ScalarType, 6, 1 > convertKeplerianToEquinoctialElements(
        const Eigen::Matrix< ScalarType, 6, 1 >& keplerianElements,
        const ScalarType centralBodyGravitationalParameter, const bool retrogradeElements = false )
{
    // Set tolerance.
    const ScalarType tolerance = 20.0 * std::numeric_limits< ScalarType >::epsilon( );

    // Determine the retrograde factor
    const ScalarType I = retrogradeElements ? -1 : 1;

    // Declare converted equinoctial elements.
    Eigen::Matrix< ScalarType, 6, 1 > equinoctialElements;

    // Read Keplerian elements
    const ScalarType a = keplerianElements( semiMajorAxisIndex );
    const ScalarType e = keplerianElements( eccentricityIndex );
    const ScalarType i = keplerianElements( inclinationIndex );
    const ScalarType omega = keplerianElements( argumentOfPeriapsisIndex );
    const ScalarType Omega = keplerianElements( longitudeOfAscendingNodeIndex );
    const ScalarType f = keplerianElements( trueAnomalyIndex );

    // Convert true anomaly to mean anomaly
    const ScalarType E = convertTrueAnomalyToEllipticalEccentricAnomaly( f, e );
    const ScalarType M = convertEllipticalEccentricAnomalyToMeanAnomaly( E, e );

    // Pre-compute shared terms
    const ScalarType longitudePeriapsis = omega + I * Omega;
    const ScalarType tani2 = std::pow( std::tan( i / 2 ), I );

    // Store components
    equinoctialElements( semiMajorAxisIndex ) = a;
    equinoctialElements( hIndex ) = e * std::sin( longitudePeriapsis );
    equinoctialElements( kIndex ) = e * std::cos( longitudePeriapsis );
    equinoctialElements( pIndex ) = tani2 * std::sin( Omega );
    equinoctialElements( qIndex ) = tani2 * std::cos( Omega );
    equinoctialElements( meanLongitudeIndex ) = M + longitudePeriapsis;

    // Return converted equinoctial elements.
    return equinoctialElements;
}


//! Convert equinoctial to Keplerian orbital elements.
/*!
 * Converts equinoctial to Keplerian orbital elements.
 * [Semianalytic satellite theory, Danielson (1995), Section 2.1.3]
 * \param keplerianElements Vector containing equinoctial elements. Order of elements is important!
 *          equinoctialElements( 0 ) = semi-major axis,                                         [m]
 *          equinoctialElements( 1 ) = h component,                                             [-]
 *          equinoctialElements( 2 ) = k component,                                             [-]
 *          equinoctialElements( 3 ) = p component,                                             [-]
 *          equinoctialElements( 4 ) = q component,                                             [-]
 *          equinoctialElements( 5 ) = mean longitude.                                        [rad]
 * \param centralBodyGravitationalParameter Gravitational parameter of central body.
 * \param retrogradeElements Whether to use the retrograde set of elements (to avoid singularities for orbits
 * with inclinations of 180 degrees). Default is false, which avoids singularities for equatorial orbits.
 * \return Converted state in Keplerian elements. The order of elements is fixed!
 *          keplerianElements( 0 ) = semiMajorAxis,                                             [m]
 *          keplerianElements( 1 ) = eccentricity,                                              [-]
 *          keplerianElements( 2 ) = inclination,                                             [rad]
 *          keplerianElements( 3 ) = argument of periapsis,                                   [rad]
 *          keplerianElements( 4 ) = longitude of ascending node,                             [rad]
 *          keplerianElements( 5 ) = true anomaly.                                            [rad]
 */
template< typename ScalarType = double >
Eigen::Matrix< ScalarType, 6, 1 > convertEquinoctialToKeplerianElements(
        const Eigen::Matrix< ScalarType, 6, 1 >& equinoctialElements,
        const ScalarType centralBodyGravitationalParameter, const bool retrogradeElements = false )
{
    using std::pow;
    using std::sqrt;
    using namespace tudat::mathematical_constants;
    using namespace tudat::basic_mathematics;

    // Set tolerance.
    const ScalarType tolerance = 20.0 * std::numeric_limits< ScalarType >::epsilon( );

    // Determine the retrograde factor
    const ScalarType I = retrogradeElements ? -1 : 1;

    // Declare converted equinoctial elements.
    Eigen::Matrix< ScalarType, 6, 1 > keplerianElements;

    // Read equinoctial elements
    const ScalarType a =      equinoctialElements( semiMajorAxisIndex );
    const ScalarType h =      equinoctialElements( hIndex );
    const ScalarType k =      equinoctialElements( kIndex );
    const ScalarType p =      equinoctialElements( pIndex );
    const ScalarType q =      equinoctialElements( qIndex );

    // Compute repeated parameters
    const ScalarType sqrth2k2 = sqrt( pow( h, 2 ) + pow( k, 2 ) );
    const ScalarType sqrtp2q2 = sqrt( pow( p, 2 ) + pow( q, 2 ) );

    // Compute auxiliary angle [Eq. (1)]
    const ScalarType sinAux = h / sqrth2k2;
    const ScalarType cosAux = k / sqrth2k2;
    const ScalarType aux = std::atan2( sinAux, cosAux );

    // Compute longitude of the AN [Eq. (2)]
    const ScalarType sinOmega = p / sqrtp2q2;
    const ScalarType cosOmega = q / sqrtp2q2;
    const ScalarType Omega = std::atan2( sinOmega, cosOmega );

    // Get true longitude
    const ScalarType L = getTrueLongitude( equinoctialElements );

    // Store values [Eq. (2), (3)]
    keplerianElements( semiMajorAxisIndex )             = a;
    keplerianElements( eccentricityIndex )              = sqrth2k2;
    keplerianElements( inclinationIndex )               = computeModulo( PI / 2 * ( 1 - I ) +
                                                                         2 * I * std::atan( sqrtp2q2 ), 2*PI );
    keplerianElements( argumentOfPeriapsisIndex )       = computeModulo( aux - I * Omega, 2*PI );
    keplerianElements( longitudeOfAscendingNodeIndex )  = computeModulo( Omega, 2*PI );
    keplerianElements( trueAnomalyIndex )               = computeModulo( L - aux, 2*PI );

    // Return converted equinoctial elements.
    return keplerianElements;
}



} // namespace orbital_element_conversions

} // namespace tudat

#endif // TUDAT_ORBITAL_ELEMENT_CONVERSIONS_H
