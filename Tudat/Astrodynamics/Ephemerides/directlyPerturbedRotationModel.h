#ifndef DIRECTLYPERTURBEDROTATIONMODEL_H
#define DIRECTLYPERTURBEDROTATIONMODEL_H

#include <map>
#include <vector>

#include "Tudat/Mathematics/BasicMathematics/mathematicalConstants.h"
#include "Tudat/Astrodynamics/BasicAstrodynamics/orbitalElementConversions.h"

#include "Tudat/Basics/utilities.h"
#include "Tudat/Astrodynamics/ReferenceFrames/referenceFrameTransformations.h"
#include "Tudat/Astrodynamics/Ephemerides/rotationalEphemeris.h"

namespace tudat
{

enum PerturbedRotationModelComponents
{
    prime_meridian_angle = 0,
    declination_angle = 1,
    right_ascension_angle = 2
};

namespace ephemerides
{

//! Rotation model for the orientation of a body, modelled as SimpleRotationalEphemeris with periodic and polynomial perturbations to Euler angles
/*!
 *  Rotation model for the orientation of a body, modelled as SimpleRotationalEphemeris with periodic and polynomial perturbations to Euler angles.
 *  Periods and sine and cosine amplitudes for any number of periodic perturbations on right ascension of north pole, declination of north pole
 *  and longitude of prime meridian (each separately) are provided. The model is based on the discussion model by LeMaistre et al. (2013)
 */
class DirectlyPerturbedRotationModel: public RotationalEphemeris
{
public:

    //! Constructor for rotation model.
    /*!
     *  Constructor for rotation model, takes polynomial and periodic terms for the three euler angles (right ascension, declination, longitude of
     *  prime meridian).
     *  NOTE: Put in reference time not equal to zero.
     *  \param rightAscensionPolynomialTerms Polynomial terms for right ascension, i.e. n^(th) term from vector is multiplied by currentTime^(n) to
     *  provide zeroth order term. Zeroth term is nominal, unperturbed value.
     *  \param declinationPolynomialTerms Polynomial terms for declination, i.e. n^(th) term from vector is multiplied by currentTime^(n) to
     *  provide zeroth order term. Zeroth term is nominal, unperturbed value.
     *  \param primeMeridianPolynomialTerms Polynomial terms for longitude of prime meridian, i.e. n^(th) term from vector is multiplied by
     *  currentTime^(n) to  provide zeroth order term. Zeroth term is value at reference epoch, firts term is rotation rate about z-axis
     *  (rotationRate_ in SimpleRotationalEphemeris).
     *  \param rightAscensionLibrations Periodic correction terms to right ascension of north pole, key of map is period of term, value is
     *  a pair of cosine term and sine term amplitude.
     *  \param declinationLibrations Periodic correction terms to declination of north pole, key of map is period of term, value is
     *  a pair of cosine term and sine term amplitude.
     *  \param primeMeridianLibrations Periodic correction terms to longitude of prime meridian, key of map is period of term, value is
     *  a pair of cosine term and sine term amplitude.
     *  \param fromIntermediateFrameToBaseFrame Rotation from frame wrt which the Euler angles are defined to the base frame of the rotation model.
     *  Variable is included to simplify rotation model creation (i.e. prevent use of CompositeRotation from ConstantRotationModel and this class)
     *  \param originalFrameOrientation Base frame identifier.
     *  \param targetFrameOrientation Target frame (i.e. body-fixed) identifier.
     */
    DirectlyPerturbedRotationModel(
            std::vector< double > rightAscensionPolynomialTerms,
            std::vector< double > declinationPolynomialTerms,
            std::vector< double > primeMeridianPolynomialTerms,
            std::map< double, std::pair< double, double > > rightAscensionLibrations,
            std::map< double, std::pair< double, double > > declinationLibrations,
            std::map< double, std::pair< double, double > > primeMeridianLibrations,
            const Eigen::Matrix3d fromIntermediateFrameToBaseFrame = Eigen::Matrix3d::Identity( ),
            const std::string originalFrameOrientation = "ECLIPJ2000",
            const std::string targetFrameOrientation = "" ):
        RotationalEphemeris( originalFrameOrientation, targetFrameOrientation ),
        rightAscensionPolynomialTerms_( rightAscensionPolynomialTerms ), declinationPolynomialTerms_( declinationPolynomialTerms ),
        primeMeridianPolynomialTerms_( primeMeridianPolynomialTerms) , rightAscensionLibrations_( rightAscensionLibrations ),
        declinationLibrations_( declinationLibrations ), primeMeridianLibrations_( primeMeridianLibrations ),
        fromIntermediateFrameToBaseFrame_( fromIntermediateFrameToBaseFrame )
    {
        // Initialize current time to nonsense value.
        currentTime_ = -1.0E150;

        usePowerSeriesForPrimeMeridian = 0;

    }

    //! Function to calculate the rotation quaternion from target frame (body-fixed) to original frame.
    /*!
     *  Function to calculate the rotation quaternion from target frame (body-fixed) to original frame at specified time.
     *  \param ephemerisTime Time at which rotation is to be calculated.
     *  \return Rotation from target (body-fixed) to original (typically inertial) frame at specified time.
     */
    Eigen::Quaterniond getRotationToBaseFrame( const double time );

    Eigen::Quaterniond getRotationToTargetFrame(
            const double time )
    {
        return getRotationToBaseFrame( time ).inverse( );
    }

    //! Function to calculate the derivative of the rotation matrix from target frame (body-fixed) to original frame.
    /*!
     *  Function to calculate the derivative of the rotation matrix from target frame (body-fixed) to original frame at specified time.
     *  \param ephemerisTime Time at which rotation is to be calculated.
     *  \return Derivative of the rotation matrix from target (body-fixed) to original (typically inertial) frame at specified time.
     */
    Eigen::Matrix3d getDerivativeOfRotationToBaseFrame( const double time );

    Eigen::Matrix3d getDerivativeOfRotationToTargetFrame( const double time )
    {
        return getDerivativeOfRotationToBaseFrame( time ).transpose( );
    }

    //! Function to get the right ascension of the north pole at given time.
    /*!
     *  Function to get the right ascension of the north pole at given time. Calculates all Euler angles at given time, if necessary, i.e. if currentTime_ is
     *  unequal to time.
     *  \param time Time at which reight ascension is to be returned.
     *  \return Right ascension at requested time.
     */
    double getCurrentRightAscension( const double time )
    {
        // Check if update function needs to be called.
        //if( std::fabs( time - currentTime_ ) > 1.0E-7 )
        {
            update( time );
        }

        return currentRightAscension_;
    }

    //! Function to get the declination of the north pole at given time.
    /*!
     *  Function to get the declination of the north pole at given time. Calculates all Euler angles at given time, if necessary, i.e. if currentTime_ is
     *  unequal to time.
     *  \param time Time at which declination is to be returned.
     *  \return Declination at requested time.
     */
    double getCurrentDeclination( const double time )
    {
        // Check if update function needs to be called.
        //if( std::fabs( time - currentTime_ ) > 1.0E-7 )
        {
            update( time );
        }

        return currentDeclination_;
    }

    //! Function to get the longitude of the prime meridian at given time.
    /*!
     *  Function to get the longitude of the prime meridian at given time. Calculates all Euler angles at given time, if necessary, i.e. if
     *  currentTime_ is unequal to time.
     *  \param time Time at which longitude of the prime meridian is to be returned.
     *  \return Longitude of the prime meridian at requested time.
     */
    double getCurrentPrimeMeridian( const double time )
    {
        // Check if update function needs to be called.
        //if( std::fabs( time - currentTime_ ) > 1.0E-7 )
        {
            update( time );
        }

        return currentPrimeMeridian_;
    }

    //! Function to get the cosine and sine amplitudes of the requested periodic rotation variation.
    /*!
     *  Function to get the cosine and sine amplitudes of the requested periodic rotation variation. The component of which the amplitudes
     *  are to be returned is identified by itd type (i.e. right ascension, declination, prime meridian) and the period of the component.
     *  \param componentType Euler angle for which the rotation vaiation is to be returned.
     *  \param componentPeriod Period of rotation variation of given component type which is to be returned.
     *  \return Cosine and sine amplitudes of the requested periodic rotation variation
     */
    std::pair< double, double > getPeriodicComponentAmplitudes( const PerturbedRotationModelComponents componentType,
                                                                const double componentPeriod );

    double getPolynomialComponent( const PerturbedRotationModelComponents componentType, const int power );

    //! Function to set the cosine and sine amplitudes of the requested periodic rotation variation.
    /*!
     *  Function to set the cosine and sine amplitudes of the requested periodic rotation variation. The component of which the amplitudes
     *  are to be reset is identified by itd type (i.e. right ascension, declination, prime meridian) and the period of the component.
     *  \param componentType Euler angle for which the rotation vaiation is to be reset.
     *  \param componentPeriod Period of rotation variation of given component type which is to be reset.
     */
    void setPeriodicComponentAmplitudes( const PerturbedRotationModelComponents componentType,
                                         const double componentPeriod,
                                         std::pair< double, double > newAmplitudes );

    void setPolynomialComponent( const PerturbedRotationModelComponents componentType, const int componentPower, double newAmplitude );

    std::vector< double > getRightAscensionTermPeriods( )
    {
        return utilities::createVectorFromMapKeys( rightAscensionLibrations_ );
    }

    std::vector< double > getDeclinationTermPeriods( )
    {
        return utilities::createVectorFromMapKeys( declinationLibrations_ );
    }

    std::vector< double > getPrimeMeridianTermPeriods( )
    {
        return utilities::createVectorFromMapKeys( primeMeridianLibrations_ );
    }

    Eigen::Matrix3d getFromIntermediateFrameToBaseFrame( )
    {
        return fromIntermediateFrameToBaseFrame_;
    }

    void update( const double time );

private:

    //! Current value of right ascension of north pole
    /*!
     *  Current value of right ascension of north pole at currentTime_, as set by previous call to update function.
     */
    double currentRightAscension_;

    //! Current value of declination of north pole
    /*!
     *  Current value of declination of north pole at currentTime_, as set by previous call to update function.
     */
    double currentDeclination_;

    //! Current value of longitude of prime meridian
    /*!
     *  Current value of longitude of prime meridian at currentTime_, as set by previous call to update function.
     */
    double currentPrimeMeridian_;

    //! Current rotation from local (body-fixed) to base frame.
    /*!
     *  Current rotation from local (body-fixed) to base frame at currentTime_, as set by previous call to update function.
     */
    Eigen::Quaterniond currentRotationFromLocalFrame_;

    //! Current time of euler angles
    /*!
     *  Current time of euler angles, i.e. the currentRightAscension_, currentDeclination_, currentPrimeMeridian_ and currentRotationFromLocalFrame_
     *  member variables, updated to current value by the update function.
     */
    double currentTime_;

    //! Polynomial terms for right ascension.
    /*!
     *  Polynomial terms for right ascension, i.e. n^(th) term from vector is multiplied by currentTime^(n) to
     *  provide zeroth order term. Zeroth term is nominal, unperturbed value.
     */
    std::vector< double > rightAscensionPolynomialTerms_;

    //! Polynomial terms for declination.
    /*!
     *  Polynomial terms for declination, i.e. n^(th) term from vector is multiplied by currentTime^(n) to
     *  provide zeroth order term. Zeroth term is nominal, unperturbed value.
     */
    std::vector< double > declinationPolynomialTerms_;

    //! Polynomial terms for longitude of prime meridian.
    /*!
     *  Polynomial terms for longitude of prime meridian, i.e. n^(th) term from vector is multiplied by
     *  currentTime^(n) to  provide zeroth order term. Zeroth term is value at reference epoch, firts term is rotation rate about z-axis
     *  (rotationRate_ in SimpleRotationalEphemeris).
     */
    std::vector< double > primeMeridianPolynomialTerms_;

    //! Periodic correction terms to right ascension of north pol
    /*!
     *  Periodic correction terms to right ascension of north pole, key of map is period of term, value is
     *  a pair of cosine term and sine term amplitude.
     */
    std::map< double, std::pair< double, double > > rightAscensionLibrations_;

    //! Periodic correction terms to declination of north pole
    /*!
     *  Periodic correction terms to declination of north pole, key of map is period of term, value is
     *  a pair of cosine term and sine term amplitude.
     */
    std::map< double, std::pair< double, double > > declinationLibrations_;

    //! Periodic correction terms to longitude of prime meridian
    /*!
     *  Periodic correction terms to longitude of prime meridian, key of map is period of term, value is
     *  a pair of cosine term and sine term amplitude.
     */
    std::map< double, std::pair< double, double > > primeMeridianLibrations_;

    //! Rotation from frame wrt which the Euler angles are defined to the base frame of the rotation model.
    /*!
     *  Rotation from frame wrt which the Euler angles are defined to the base frame of the rotation model.
     *  Variable is included to simplify rotation model creation (i.e. prevent use of CompositeRotation from ConstantRotationModel and this class)
     */
    Eigen::Matrix3d fromIntermediateFrameToBaseFrame_;

    //! Iterator over periodic rotation correction terms.
    /*!
     *  Iterator over periodic rotation correction terms. Used by update function. Defined here to prevent numerous construction and destructions
     *  of iterator by the update function.
     */
    std::map< double, std::pair< double, double > >::iterator librationIterator;

    bool usePowerSeriesForPrimeMeridian;

};

}

}

#endif // DIRECTLYPERTURBEDROTATIONMODEL_H
