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
 *      Varschalovich, et al. Quantum Theory of Angular Momentum. World Scientific, February 1988
 *      Boué, (2017) The two rigid body interaction using angular momentum theory formulae, CMDA 128(2-3):261-273
 *
 *
 */

#ifndef TUDAT_CALYLEY_KLEIN_PARAMETERS_H
#define TUDAT_CALYLEY_KLEIN_PARAMETERS_H

#include <Eigen/Core>
#include <Eigen/Geometry>

#include <complex>

#include "Tudat/Mathematics/BasicMathematics/mathematicalConstants.h"

namespace tudat
{
namespace basic_mathematics
{

//! Function to convert a quaternion representing a rotation to two Cayley-Klein parameters
/*!
 * Function to convert a quaternion (q=q_{0},q_{1},q_{2},q_{3}) representing a rotation to two Cayley-Klein parameters: a and b.
 * The parameters a and b are computed as: as a = q_{0} - i*q_{3}; b = q_{1} - i*q_{2};
 * \param quaternion Quaternion from which Cayley-Klein parameters are determined
 * \param cayleyKleinA Cayley-Klein parameter a (returned by reference)
 * \param cayleyKleinB Cayley-Klein parameter b (returned by reference)
 */
void convertQuaterionToCayleyKleinParameters(
        const Eigen::Quaterniond quaternion, std::complex< double >& cayleyKleinA, std::complex< double >& cayleyKleinB );

//! Function to convert 3-2-3 Euler angles representing a rotation to two Cayley-Klein parameters
/*!
 * Function to convert 3-2-3 Euler angles representing a rotation to two Cayley-Klein parameters: a and b. Algorithm is taken
 * from Varschalovich et al. (1988)
 * \param firstZRotation Euler angle for first rotation about z-axis
 * \param yRotation Euler angle for rotation about y-axis
 * \param secondZRotation Euler angle for second rotation about z-axis
 * \param cayleyKleinA Cayley-Klein parameter a (returned by reference)
 * \param cayleyKleinB Cayley-Klein parameter b (returned by reference)
 */
void convert323EulerAnglesToCayleyKleinParameters(
        const double firstZRotation, const double yRotation, const double secondZRotation,
        std::complex< double >& cayleyKleinA, std::complex< double >& cayleyKleinB );

//! Function to convert 3-1-3 Euler angles representing a rotation to two Cayley-Klein parameters
/*!
 * Function to convert 3-1-3 Euler angles representing a rotation to two Cayley-Klein parameters: a and b. Algorithm is taken
 * from Varschalovich et al. (1988)
 * \param firstZRotation Euler angle for first rotation about z-axis
 * \param xRotation Euler angle for rotation about x-axis
 * \param secondZRotation Euler angle for second rotation about z-axis
 * \param cayleyKleinA Cayley-Klein parameter a (returned by reference)
 * \param cayleyKleinB Cayley-Klein parameter b (returned by reference)
 */
void convert313EulerAnglesToCayleyKleinParameters(
        const double firstZRotation, const double xRotation, const double secondZRotation,
        std::complex< double >& cayleyKleinA, std::complex< double >& cayleyKleinB );

} // namespace basic_mathematics

} // namespace tudat

#endif // TUDAT_CALYLEY_KLEIN_PARAMETERS_H
