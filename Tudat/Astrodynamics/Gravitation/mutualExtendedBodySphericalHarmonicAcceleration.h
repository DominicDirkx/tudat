#ifndef MUTUALEXTENDEDBODYSPHERICALHARMONICACCELERATION_H
#define MUTUALEXTENDEDBODYSPHERICALHARMONICACCELERATION_H

#include <boost/function.hpp>
#include <boost/bind.hpp>
#include <boost/tuple/tuple.hpp>


#include <Eigen/Core>
#include <Eigen/Geometry>

#include <Tudat/Astrodynamics/BasicAstrodynamics/accelerationModel.h>
#include <Tudat/Basics/basicTypedefs.h>
#include <Tudat/Mathematics/BasicMathematics/legendrePolynomials.h>

#include "Tudat/Astrodynamics/Gravitation/mutualForcePotential.h"
#include "Tudat/Astrodynamics/ReferenceFrames/referenceFrameTransformations.h"
#include "Tudat/Mathematics/BasicMathematics/sphericalHarmonicTransformations.h"

namespace tudat
{

namespace gravitation
{

//! Class to compute the full two-extended body gravitational attraction between two bodies
/*!
 *  Class to compute the full two-extended body gravitational attraction between two bodies, this involves a quadruple
 *  sum: over a set of degrees and orders of both bodies. Models for the acceleration are given by Compere and Lemaitre (2014),
 *  Boue (2017) and the details of this implementation are given in Dirkx et al. (2018).
 */
class MutualExtendedBodySphericalHarmonicAcceleration: public basic_astrodynamics::AccelerationModel3d
{

public:

    //! Constructor
    /*!
     * \brief MutualExtendedBodySphericalHarmonicAcceleration
     * \param positionOfBody1Function Function returning the position of body 1
     * \param positionOfBody2Function Function returning the position of body 2
     * \param gravitationalParameterFunction Function returning the effective gravitational parameter used as 'mu' in the
     * two-body potential
     * \param equatorialRadiusOfBody1 Equatorial radius of body 1.
     * \param equatorialRadiusOfBody2 Equatorial radius of body 2.
     * \param cosineHarmonicCoefficientsOfBody1Function Function returning the spherical harmonic cosine coefficients of body 1,
     * in a frame fixed to body 1.
     * \param sineHarmonicCoefficientsOfBody1Function Function returning the spherical harmonic sine coefficients of body 1,
     * in a frame fixed to body 1.
     * \param cosineHarmonicCoefficientsOfBody2Function Function returning the spherical harmonic cosine coefficients of body 1,
     * in a frame fixed to body 2.
     * \param sineHarmonicCoefficientsOfBody2Function Function returning the spherical harmonic sine coefficients of body 1,
     * in a frame fixed to body 2.
     * \param coefficientCombinationsToUse  List of degrees/orders that are to be used for the series expansion.
     * Each tuple contains: (degree of body 1, order of body 1, degree of body 2, order of body 2)
     * \param toLocalFrameOfBody1Transformation Function that returns the rotation from the inertial to body-fixed frame of body 1
     * \param toLocalFrameOfBody2Transformation Function that returns the rotation from the inertial to body-fixed frame of body 2
     * \param useCentraBodyFrame Variable denoting whether acceleration is expressed in body-fixed frame of body 1, or in inertial frame.
     * \param areCoefficientsNormalized Boolean denoting whether spherical harmonic coefficients are normalized
     */
    MutualExtendedBodySphericalHarmonicAcceleration(
            const boost::function< Eigen::Vector3d( ) > positionOfBody1Function,
            const boost::function< Eigen::Vector3d( ) > positionOfBody2Function,
            const boost::function< double( ) > gravitationalParameterFunction,
            const double equatorialRadiusOfBody1,
            const double equatorialRadiusOfBody2,
            const boost::function< Eigen::MatrixXd( ) > cosineHarmonicCoefficientsOfBody1Function,
            const boost::function< Eigen::MatrixXd( ) > sineHarmonicCoefficientsOfBody1Function,
            const boost::function< Eigen::MatrixXd( ) > cosineHarmonicCoefficientsOfBody2Function,
            const boost::function< Eigen::MatrixXd( ) > sineHarmonicCoefficientsOfBody2Function,
            const std::vector< boost::tuple< unsigned int, unsigned int, unsigned int, unsigned int > >& coefficientCombinationsToUse,
            const boost::function< Eigen::Quaterniond( ) > toLocalFrameOfBody1Transformation,
            const boost::function< Eigen::Quaterniond( ) > toLocalFrameOfBody2Transformation,
            const bool useCentraBodyFrame,
            const bool areCoefficientsNormalized = 1 );

    //! Update member variables used by the acceleration model.
    /*!
    *  Updates member variables used by the acceleration model.
    *  Function pointers to retrieve the current values of quantities from which the
    *  acceleration is to be calculated are set by constructor. This function calls
    *  them to update the associated variables to their current state.
    *  \param currentTime Time at which acceleration model is to be updated.
    */
    void updateMembers( const double currentTime = TUDAT_NAN );

    //! Get acceleration.
    /*!
     * Returns the acceleration, as computed by updateMembers function.
     * \return Current acceleration.
     */
    Eigen::Vector3d getAcceleration( )
    {
        return currentAcceleration_;
    }

    Eigen::Vector3d getAccelerationInBodyFixedFrame( )
    {
        return currentAccelerationInBodyFixedFrame_;
    }


    //! Function to retrieve whether acceleration is expressed in body-fixed frame of body 1, or in inertial frame.
    /*!
     *  Function to retrieve whether acceleration is expressed in body-fixed frame of body 1, or in inertial frame.
     *  \return Variable denoting whether acceleration is expressed in body-fixed frame of body 1, or in inertial frame.
     */
    bool getUseCentraBodyFrame( )
    {
        return useCentraBodyFrame_;
    }

    //! Function to retrieve position of body 1 w.r.t. body 2, expressed in inertial frame, as computed by updateMembers
    /*!
     *  Function to retrieve position of body 1 w.r.t. body 2, expressed in inertial frame, as computed by updateMembers
     *  \return Position of body 1 w.r.t. body 2, expressed in inertial frame, as computed by last call to updateMembers
     */
    Eigen::Vector3d getCurrentRelativePosition( )
    {
        return currentRelativePosition_;
    }

    Eigen::Vector3d getCurrentBodyFixedRelativePosition( )
    {
        return currentBodyFixedRelativePosition_;
    }



    //! Function to retrieve current rotation from body-fixed frame of body 2 to that of body 1, as computed by updateMembers
    /*!
     *  Function to retrieve current rotation from body-fixed frame of body 2 to that of body 1, as computed by updateMembers
     *  \return Current rotation from body-fixed frame of body 2 to that of body 1, as computed by last call to updateMembers
     */
    Eigen::Quaterniond getCurrentRotationFromBody2ToBody1( )
    {
        return currentRotationFromBody2ToBody1_;
    }

    //! Function to retrieve object that computes the effect one-body coefficients
    /*!
     *  Function to retrieve object that computes the effect one-body coefficients from the two one-body gravity fields and
     *  their relative orientation
     *  \return Object that computes the effect one-body coefficients from the two one-body gravity fields and their relative
     *   orientation
     */
    boost::shared_ptr< gravitation::EffectiveMutualSphericalHarmonicsField > getEffectiveMutualPotentialField( )
    {
        return effectiveMutualPotentialField_;
    }

    //! Function to retrieve object that computes the coefficients of body, in the body-fixed frame of body 1
    /*!
     *  Function to retrieve object that computes the coefficients of body, in the body-fixed frame of body 1
     *  \return Object that computes the coefficients of body, in the body-fixed frame of body 1
     */
    boost::shared_ptr< basic_mathematics::SphericalHarmonicsCache > getSphericalHarmonicsCache( )
    {
        return sphericalHarmonicsCache_;
    }

    //! Function to retrieve equatorial radius of body 1.
    /*!
     *  Function to retrieve equatorial radius of body 1.
     *  \return Equatorial radius of body 1.
     */
    double getEquatorialRadiusOfBody1( )
    {
        return equatorialRadiusOfBody1_;
    }

    //! Function to retrieve equatorial radius of body 2.
    /*!
     *  Function to retrieve equatorial radius of body 2.
     *  \return Equatorial radius of body 2.
     */
    double getEquatorialRadiusOfBody2( )
    {
        return equatorialRadiusOfBody2_;
    }

    //! Function to retrieve a single integer power of (distance / equatorial radius of body 1)
    /*!
     *  Function to retrieve a single integer power of (distance / equatorial radius of body 1)
     *  \param index Power to which (distance / equatorial radius of body 1) is to be computed
     *  \return Single integer power of (distance / equatorial radius of body 1)
     */
    double getRadius1Power( const int index )
    {
        return radius1Powers_.at( index );
    }

    //! Function to retrieve a single integer power of (distance / equatorial radius of body 2)
    /*!
     *  Function to retrieve a single integer power of (distance / equatorial radius of body 2)
     *  \param index Power to which (distance / equatorial radius of body 2) is to be computed
     *  \return Single integer power of (distance / equatorial radius of body 2)
     */
    double getRadius2Power( const int index )
    {
        return radius2Powers_.at( index );
    }

    //! Function to retrieve list of integer powers of (distance / equatorial radius of body 1)
    /*!
     *  Function to retrieve list of integer powers of (distance / equatorial radius of body 1)
     *  \return List of integer powers of (distance / equatorial radius of body 1)
     */
    std::vector< double > getRadius1Powers( )
    {
        return radius1Powers_;
    }

    //! Function to retrieve list of integer powers of (distance / equatorial radius of body 2)
    /*!
     *  Function to retrieve list of integer powers of (distance / equatorial radius of body 2)
     *  \return List of integer powers of (distance / equatorial radius of body 2)
     */
    std::vector< double > getRadius2Powers( )
    {
        return radius2Powers_;
    }

private:

    //! Function returning the position of body 1
    boost::function< Eigen::Vector3d( ) > positionOfBody1Function_;

    //! Function returning the position of body 2
    boost::function< Eigen::Vector3d( ) > positionOfBody2Function_;

    //! Function returning the effective gravitational parameter
    boost::function< double( ) > gravitationalParameterFunction_;

    //! Equatorial radius of body 1.
    double equatorialRadiusOfBody1_;

    //! Equatorial radius of body 2.
    double equatorialRadiusOfBody2_;

    //! List of degrees/orders that are to be used for the series expansion.
    /*!
     * List of degrees/orders that are to be used for the series expansion. Each tuple contains: (degree of body 1, order of
     * body 1, degree of body 2, order of body 2)
     */
    std::vector< boost::tuple< unsigned int, unsigned int, unsigned int, unsigned int > > coefficientCombinationsToUse_;

    //! Function that returns the effective one-body spherical harmonic cosine coefficient
    /*!
     *  Function that returns the effective one-body spherical harmonic cosine coefficient as a function of (degree of body 1,
     *  order of body 1, degree of body 2, order of body 2)
     */
    boost::function< double( int, int, int, int ) > effectiveCosineCoefficientFunction_;

    //! Function that returns the effective one-body spherical harmonic sine coefficient
    /*!
     *  Function that returns the effective one-body spherical harmonic sine coefficient as a function of (degree of body 1,
     *  order of body 1, degree of body 2, order of body 2)
     */
    boost::function< double( int, int, int, int ) > effectiveSineCoefficientFunction_;

    //! Function that returns the rotation from the inertial to body-fixed frame of body 1
    boost::function< Eigen::Quaterniond( ) > toLocalFrameOfBody1Transformation_;

    //! Function that returns the rotation from the inertial to body-fixed frame of body 2
    boost::function< Eigen::Quaterniond( ) > toLocalFrameOfBody2Transformation_;

    //! Maximum degree of effective one-body spherical harmonic expansion
    unsigned int maximumDegree_;

    //! Maximum order of effective one-body spherical harmonic expansion
    unsigned int maximumOrder_;

    //! Object that computes the effect one-body coefficients from the two one-body gravity fields and their relative orientation
    boost::shared_ptr< gravitation::EffectiveMutualSphericalHarmonicsField > effectiveMutualPotentialField_;

    //! Object that computes the coefficients of body, in the body-fixed frame of body 1
    boost::shared_ptr< basic_mathematics::SphericalHarmonicsCache > sphericalHarmonicsCache_;

    //! Acceleration, as computed by last call to updateMembers
    Eigen::Vector3d currentAcceleration_;

    Eigen::Vector3d currentAccelerationInBodyFixedFrame_;

    //! Position of body 1 w.r.t. body 2, expressed in inertial frame, as computed by last call to updateMembers
    Eigen::Vector3d currentRelativePosition_;

    //! Position of body 1 w.r.t. body 2, expressed in body-fixed frame of body 1, as computed by last call to updateMembers
    Eigen::Vector3d currentBodyFixedRelativePosition_;

    //! Current rotation from body-fixed frame of body 2 to that of body 1, as computed by last call to updateMembers
    Eigen::Quaterniond currentRotationFromBody2ToBody1_;

    //! Variable denoting whether acceleration is expressed in body-fixed frame of body 1, or in inertial frame.
    bool useCentraBodyFrame_;

    //! Boolean denoting whether spherical harmonic coefficients are normalized
    bool areCoefficientsNormalized_;

    //! List of integer powers of (distance / equatorial radius of body 1)
    std::vector< double > radius1Powers_;

    //! List of integer powers of (distance / equatorial radius of body 2)
    std::vector< double > radius2Powers_;


};

}

}

#endif // MUTUALEXTENDEDBODYSPHERICALHARMONICACCELERATION_H
