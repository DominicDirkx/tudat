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
 *      Eigen. Structures having Eigen members,
 *          http://eigen.tuxfamily.org/dox/TopicStructHavingEigenMembers.html, last accessed: 5th
 *          March, 2013.
 *		Colin R. McInnes, Solar Sailing: Technology, Dynamics and Mission Applications,
 *			Springer-Verlag Berlin Heidelberg, 2004.
 *
 */

#ifndef TUDAT_PERFECT_REFLECTION_SAIL_ACCELERATION_H
#define TUDAT_PERFECT_REFLECTION_SAIL_ACCELERATION_H 

#include <boost/function.hpp>
#include <boost/lambda/lambda.hpp>
#include <boost/shared_ptr.hpp>
#include <iostream>

#include <Eigen/Core>

#include "Tudat/Astrodynamics/BasicAstrodynamics/accelerationModel.h"
#include "Tudat/Astrodynamics/ReferenceFrames/referenceFrameTransformations.h"

//#include "Tudat/Astrodynamics/ElectroMagnetism/cannonBallRadiationPressureForce.h"

namespace tudat
{
namespace electro_magnetism
{

//! Compute radiation pressure acceleration using a perfectly reflecting sail model.
/*!
 * Computes radiation pressure acceleration using a perfectly reflecting  sail model, i.e. assuming force to be
 * pointing along the sail normal vector (given by sail control angles). This function is essentially a wrapper for the
 * function that computes the force.
 *
 * \param  sailLightnessNumber Ratio of radiation pressure acceleration to gravitational acceleration. [-]
 * \param distanceToSource Distance from target to source 										       [m]
 * \param sourceBodyGravitationalParameter Gravitational Parameter of the source body.                 [m^3/s^2]
 * \param sailControlAngles Vector containing the sail control angles - 1)cone and 2)clock angle.      [rad]
 * \return Acceleration due to radiation pressure on a perfectly reflecting sail.                      [m/s^2]
 * \sa computePerfectReflectionSailAcceleration().
 */
Eigen::Vector3d computePerfectReflectionSailAcceleration(
        const double sailLightnessNumber,
        const double distanceToSource,
        const double sourceBodyGravitationalParameter,
        const Eigen::Vector2d sailControlAngles	);

//! Perfectly Reflecting  sail acceleration model class.
/*!
 * Class that can be used to compute the radiation pressure using a perfectly reflecting sail model, i.e.,
 * assuming force to be pointing along the sail normal vector (given by cone and clock angles).
 */
class PerfectReflectionSailAcceleration: public basic_astrodynamics::AccelerationModel3d
{
private:

    //! Typedef for double-returning function.
    typedef boost::function< double( ) > DoubleReturningFunction;

    //! Typedef for Eigen::Vector3d-returning function.
    typedef boost::function< Eigen::Vector6d( ) > Vector6dReturningFunction;

    //! Typedef for Eigen::Vector2d-returning function.
    typedef boost::function< Eigen::Vector2d( ) > Vector2dReturningFunction;

public:


    //! Constructor taking function pointers for all variables.
    /*!
     * Constructor taking function pointers for all variables.
     * \param sourceBodyStateFunction Function returning state of radiation source.
     * \param acceleratedBodyStateFunction Function returning state of body undergoing
     *          radiation pressure acceleration.
     * \param controlAngleFunction Function returning current value of sail control angles.
     * \param lightnessNumberFunction Function returning current sail lightness number.
     * \param gravitationalParameterFunction Function returning gravitational parameter of source body.
     */
    PerfectReflectionSailAcceleration(
            const Vector6dReturningFunction sourceBodyStateFunction,
            const Vector6dReturningFunction acceleratedBodyStateFunction,
            const Vector2dReturningFunction controlAngleFunction,
            const DoubleReturningFunction lightnessNumberFunction,
            const DoubleReturningFunction gravitationalParameterFunction,
            const DoubleReturningFunction sailEfficiencyFunction )
        : sourceBodyStateFunction_( sourceBodyStateFunction ),
          acceleratedBodyStateFunction_( acceleratedBodyStateFunction ),
          controlAngleFunction_( controlAngleFunction ),
          lightnessNumberFunction_( lightnessNumberFunction ),
          gravitationalParameterFunction_( gravitationalParameterFunction ),
          sailEfficiencyFunction_( sailEfficiencyFunction )
    { }

    //! Constructor taking functions pointers and constant values for parameters.
    /*!
     * Constructor taking function pointers for state vectors and radiation pressure and
     * constant values for other parameters.
     * \param sourceBodyStateFunction Function returning state of radiation source.
     * \param acceleratedBodyStateFunction Function returning state of body undergoing
     *              radiation pressure acceleration.
     * \param controlAngleFunction Function returning current value of sail control angles.
     * \param sailLightnessNumber Constant Lightness number of the sail.
     * \param sourceBodyGravitationalParameter Constant Gravitational parameter of source body.
     */
    PerfectReflectionSailAcceleration(
            const Vector6dReturningFunction sourceBodyStateFunction,
            const Vector6dReturningFunction acceleratedBodyStateFunction,
            const Vector2dReturningFunction controlAngleFunction,
            const double sailLightnessNumber,
            const double sourceBodyGravitationalParameter,
            const double sailEfficiency )
        : sourceBodyStateFunction_( sourceBodyStateFunction ),
          acceleratedBodyStateFunction_( acceleratedBodyStateFunction ),
          controlAngleFunction_( controlAngleFunction ),
          lightnessNumberFunction_(
              boost::lambda::constant( sailLightnessNumber ) ),
          gravitationalParameterFunction_( boost::lambda::constant( sourceBodyGravitationalParameter ) ),
          sailEfficiencyFunction_( boost::lambda::constant( sailEfficiency ) )
    { }

    //! Get radiation pressure acceleration.
    /*!
     * Returns the radiation pressure acceleration. No arguments are passed to this function.
     * Instead, all data required for computation is to be obtained from pointers to functions/
     * classes/structs, etc., which are to be set in a derived class and evaluated by the
     * updateMembers() function below. This function is essentially a wrapper for the free
     * function that computes the radiation pressure acceleration.
     * \return Radiation pressure acceleration.
     * \sa computePerfectReflectionSailAcceleration().
     */
    Eigen::Vector3d getAcceleration( )
    {
        return currentAcceleration_;
    }

    //! Update member variables used by the radiation pressure acceleration model.
    /*!
     * Updates member variables used by the acceleration model. This function evaluates all
     * dependent variables to the 'current' values of these parameters. Only these current values,
     * not the function-pointers are then used by the getAcceleration( ) function.
     * \param currentTime Time at which acceleration model is to be updated.
     */
    void updateMembers( const double currentTime = TUDAT_NAN )
    {
        if( !( this->currentTime_ == currentTime ) )
        {
            currentSourceBodyState_ = sourceBodyStateFunction_( );
            currentAcceleratedBodyState_ = acceleratedBodyStateFunction_( );
            currentDistanceToSource_ = ( currentSourceBodyState_.segment(0,3)
                                         - currentAcceleratedBodyState_.segment(0,3) ).norm( );

            currentControlAngles_ = controlAngleFunction_( );
            currentLightnessNumber_ = lightnessNumberFunction_( );
            currentGravitationalParameter_ = gravitationalParameterFunction_( );
            currentAcceleration_ =
                    sailEfficiencyFunction_( ) *
                    reference_frames::getInertialToRswSatelliteCenteredFrameRotationMatrx(
                        currentAcceleratedBodyState_ - currentSourceBodyState_ ).transpose( ) *
                    computePerfectReflectionSailAcceleration(
                        currentLightnessNumber_, currentDistanceToSource_, currentGravitationalParameter_,
                        currentControlAngles_ );
        }
    }

    //! Function to retrieve the function pointer returning mass of accelerated body.
    /*!
     * Function to retrieve the function pointer returning mass of accelerated body.
     * \return Function pointer returning mass of accelerated body.
     */


private:

    //! Function pointer returning state of source.
    /*!
     * Function pointer returning state of source (6D vector).
     */
    const Vector6dReturningFunction sourceBodyStateFunction_;

    //! Function pointer returning state of accelerated body.
    /*!
     * Function pointer returning state of accelerated body (6D vector).
     */
    const Vector6dReturningFunction acceleratedBodyStateFunction_;

    //! Function pointer returning sail control angles.
    /*!
     * Function pointer returning sail control angles (2D vector).
        */
    const Vector2dReturningFunction controlAngleFunction_;

    //! Function pointer returning sail lightness number.
    /*!
     * Function pointer returning sail lightness number [-].
     */
    const DoubleReturningFunction lightnessNumberFunction_;

    //! Function pointer returning source body's gravitational parameter.
    /*!
     * Function pointer returning source body's gravitational parameter [m^3/s^2].
     */
    const DoubleReturningFunction gravitationalParameterFunction_;

    const DoubleReturningFunction sailEfficiencyFunction_;


    //! Current distance from accelerated body to source.
    /*!
     * Current distance from accelerated body to source (double).
     */
    double currentDistanceToSource_;

    //! Current value of sail control angles.
    /*!
     * Current value of sail control angles (3D vector).
     */

    Eigen::Vector2d currentControlAngles_;

    //! Current value of sail's lightness number.
    /*!
     * Current value of sail's lightness number (double).
     */
    double currentLightnessNumber_;

    //! Current value of source body's gravitational parameter.
    /*!
     * Current value of source body's gravitational parameter. (double).
     */
    double currentGravitationalParameter_;

    //! Current acceleration of the accelerated body
    /*!
     * Current acceleration of the accelerated body (3D vector).
     */
    Eigen::Vector3d currentAcceleration_;

    //! Current State of the source body
    /*!
     * Current state of the source body (6D vector).
     */
    Eigen::Vector6d currentSourceBodyState_;

    //! Current state of the accelerated body
    /*!
     * Current state of the accelerated body (6D vector).
     */
    Eigen::Vector6d currentAcceleratedBodyState_;

};

//! Typedef for shared-pointer to PerfectReflectionSailAcceleration.
typedef boost::shared_ptr< PerfectReflectionSailAcceleration > PerfectReflectionSailPointer;

} // namespace electro_magnetism
} // namespace tudat

#endif // TUDAT_PERFECT_REFLECTION_SAIL_ACCELERATION_H
