/*    Copyright (c) 2010-2018, Delft University of Technology
 *    All rigths reserved
 *
 *    This file is part of the Tudat. Redistribution and use in source and
 *    binary forms, with or without modification, are permitted exclusively
 *    under the terms of the Modified BSD license. You should have received
 *    a copy of the license with this file. If not, please or visit:
 *    http://tudat.tudelft.nl/LICENSE.
 *
 */

#ifndef TUDAT_CANNONBALLRADIATIONPRESSUREINTERFACE_H
#define TUDAT_CANNONBALLRADIATIONPRESSUREINTERFACE_H

#include <vector>

#include <boost/function.hpp>
#include <boost/lambda/lambda.hpp>

#include <Eigen/Core>

#include "Tudat/Astrodynamics/BasicAstrodynamics/physicalConstants.h"
#include "Tudat/Astrodynamics/ElectroMagnetism/radiationPressureInterface.h"
#include "Tudat/Mathematics/BasicMathematics/mathematicalConstants.h"

namespace tudat
{

namespace electro_magnetism
{

//! Class in which the properties of a solar radiation pressure acceleration model are stored.
/*!
 *  Class in which the properties of a solar radiation pressure acceleration model are stored and
 *  the current radiation pressure is calculated based on the source power and geometry. The
 *  current implementation is limited to a cannonball model.
 */
class CannonBallRadiationPressureInterface: public RadiationPressureInterface
{
public:

    //! Constructor.
    /*!
     *  Class construtor for radiation pressure interface.
     *  \param sourcePower Function returning the current total power (in W) emitted by the source
     *  body.
     *  \param sourcePositionFunction Function returning the current position of the source body.
     *  \param targetPositionFunction Function returning the current position of the target body.
     *  \param radiationPressureCoefficient Reflectivity coefficient of the target body.
     *  \param area Reflecting area of the target body.
     *  \param occultingBodyPositions List of functions returning the positions of the bodies
     *  causing occultations (default none) NOTE: Multiple concurrent occultations may currently
     *  result in slighlty underestimted radiation pressure.
     *  \param occultingBodyRadii List of radii of the bodies causing occultations (default none).
     *  \param sourceRadius Radius of the source body (used for occultation calculations) (default 0).
     */
    CannonBallRadiationPressureInterface(
            const boost::function< double( ) > sourcePower,
            const boost::function< Eigen::Vector3d( ) > sourcePositionFunction,
            const boost::function< Eigen::Vector3d( ) > targetPositionFunction,
            const double radiationPressureCoefficient,
            const double area,
            const std::vector< std::string > occultingBodyNames = std::vector< std::string >( ),
            const std::vector< boost::function< Eigen::Vector3d( ) > > occultingBodyPositions =
            std::vector< boost::function< Eigen::Vector3d( ) > >( ),
            const std::vector< double > occultingBodyRadii = std::vector< double > ( ),
            const double sourceRadius = 0.0 ):
        RadiationPressureInterface( sourcePositionFunction, targetPositionFunction, occultingBodyNames,
                                    occultingBodyPositions, occultingBodyRadii, sourceRadius ),
        sourcePower_( sourcePower ),
        radiationPressureCoefficient_( radiationPressureCoefficient ),
        radiationPressureCoefficientFunction_( boost::lambda::constant( radiationPressureCoefficient ) ),
        area_( area ),
        currentRadiationPressure_( TUDAT_NAN ),
        currentSolarVector_( Eigen::Vector3d::Zero( ) ){ }

    //! Destructor
    virtual ~CannonBallRadiationPressureInterface( ){ }

    //! Function to update the current value of the radiation pressure
    /*!
     *  Function to update the current value of the radiation pressure, based on functions returning
     *  the positions of the bodies involved and the source power.
     * \param currentTime Time at which acceleration model is to be updated.
     */
    void updateInterface( const double currentTime = TUDAT_NAN );

    //! Function to return the current radiation pressure due to source at target (in N/m^2).
    /*!
     *  Function to return the current radiation pressure due to source at target (in N/m^2).
     *  \return Current radiation pressure due to source at target (in N/m^2).
     */
    double getCurrentRadiationPressure( ) const
    {
        return currentRadiationPressure_;
    }

    //! Function to return the current vector from the target to the source.
    /*!
     *  Function to return the current vector from the target to the source.
     *  \return Current vector from the target to the source.
     */
    Eigen::Vector3d getCurrentSolarVector( ) const
    {
        return currentSolarVector_;
    }


    //! Function to return the reflecting area of the target body.
    /*!
     *  Function to return the reflecting area of the target body.
     *  \return The reflecting area of the target body.
     */
    double getArea( ) const
    {
        return area_;
    }

    //! Function to return the radiation pressure coefficient of the target body.
    /*!
     *  Function to return the radiation pressure coefficient of the target body.
     *  \return The radiation pressure coefficient of the target body.
     */
    double getRadiationPressureCoefficient( ) const
    {
        return radiationPressureCoefficient_;
    }

    //! Function to reset a constant radiation pressure coefficient of the target body.
    /*!
     *  Function to reset a constant radiation pressure coefficient of the target body.
     *  \param radiationPressureCoefficient The new radiation pressure coefficient of the target body.
     */
    void resetRadiationPressureCoefficient( const double radiationPressureCoefficient )
    {
        radiationPressureCoefficient_ = radiationPressureCoefficient;
        radiationPressureCoefficientFunction_ = boost::lambda::constant( radiationPressureCoefficient );
    }

    //! Function to reset the function to obtain the radiation pressure coefficient of the target body.
    /*!
     *  Function to reset the function to obtain the radiation pressure coefficient of the target body.
     *  \param radiationPressureCoefficientFunction New function to obtain the radiation pressure coefficient of the target body.
     */
    void resetRadiationPressureCoefficientFunction(
            const boost::function< double( const double ) > radiationPressureCoefficientFunction )
    {
        radiationPressureCoefficientFunction_ = radiationPressureCoefficientFunction;
    }

    //! Function to return the function returning the current total power (in W) emitted by the
    //! source body.
    /*!
     *  Function to return the function returning the current total power (in W) emitted by the
     *  source body.
     *  \return  The function returning the current total power emitted by the source body.
     */
    boost::function< double( ) > getSourcePowerFunction( ) const
    {
        return sourcePower_;
    }

    //! Function to return the current time of interface (i.e. time of last updateInterface call).
    /*!
     *  Function to return the current time of interface (i.e. time of last updateInterface call).
     *  \return Current time of interface (i.e. time of last updateInterface call).
     */
    double getCurrentTime( )
    {
        return currentTime_;
    }


protected:

    //! Function returning the current total power (in W) emitted by the source body.
    boost::function< double( ) > sourcePower_;

    //! Radiation pressure coefficient of the target body.
    double radiationPressureCoefficient_;

    //! Function to reset a constant radiation pressure coefficient of the target body.
    boost::function< double( const double ) > radiationPressureCoefficientFunction_;

    //! Reflecting area of the target body.
    double area_;

    //! Current radiation pressure due to source at target (in N/m^2).
    double currentRadiationPressure_;

    //! Current vector from the target to the source.
    Eigen::Vector3d currentSolarVector_;

};

} // namespace electro_magnetism
} // namespace tudat

#endif // TUDAT_CANNONBALLRADIATIONPRESSUREINTERFACE_H
