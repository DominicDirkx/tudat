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

#ifndef TUDAT_SOLARSAILRADIATIONPRESSUREINTERFACE_H
#define TUDAT_SOLARSAILRADIATIONPRESSUREINTERFACE_H

#include <vector>

#include <boost/function.hpp>
#include <boost/lambda/lambda.hpp>

#include <Eigen/Core>

#include "Tudat/Basics/basicTypedefs.h"
#include "Tudat/Basics/utilities.h"
#include "Tudat/Astrodynamics/BasicAstrodynamics/physicalConstants.h"
#include "Tudat/Mathematics/BasicMathematics/mathematicalConstants.h"
#include "Tudat/Astrodynamics/ElectroMagnetism/radiationPressureInterface.h"

namespace tudat
{

namespace electro_magnetism
{

class SolarSailRadiationPressureInterface: public RadiationPressureInterface
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
    SolarSailRadiationPressureInterface(
            const boost::function< double( ) > lightnessNumber,
            const boost::function< Eigen::Vector3d( ) > sourcePositionFunction,
            const boost::function< Eigen::Vector3d( ) > targetPositionFunction,
            const boost::function< Eigen::Vector2d( const double ) > sailOrientationAnglesFunction,
            const double sailEfficiency,
            const std::vector< boost::function< Eigen::Vector3d( ) > > occultingBodyPositions =
            std::vector< boost::function< Eigen::Vector3d( ) > >( ),
            const std::vector< double > occultingBodyRadii = std::vector< double > ( ),
            const double sourceRadius = 0.0 ):
        RadiationPressureInterface(
            sourcePositionFunction, targetPositionFunction,
            occultingBodyPositions, occultingBodyRadii, sourceRadius ),
        lightnessNumber_( lightnessNumber ), sailOrientationAnglesFunction_( sailOrientationAnglesFunction ),
        sailEfficiency_( sailEfficiency ){ }

    //! Destructor
    virtual ~SolarSailRadiationPressureInterface( ){ }

    //! Function to update the current value of the radiation pressure
    /*!
     *  Function to update the current value of the radiation pressure, based on functions returning
     *  the positions of the bodies involved and the source power.
     * \param currentTime Time at which acceleration model is to be updated.
     */
    void updateInterface( const double currentTime = TUDAT_NAN )
    {
        if( !( currentTime == currentTime_ ) )
        {
            currentLightnessNumber_ = lightnessNumber_( );
            currentSailAngles_ = sailOrientationAnglesFunction_( currentTime );

            updateShadowFunction( currentTime );
            currentEffectiveSailEfficiency_ = sailEfficiency_ * currentShadowFunction_;

            currentTime_ = currentTime;
        }
    }

    double getCurrentLightnessNumber( )
    {
        return currentLightnessNumber_;
    }

    double getCurrentSailEfficiency( )
    {
        return currentLightnessNumber_;
    }

    Eigen::Vector2d getCurrentSailAngles( )
    {
        return currentSailAngles_;
    }

protected:

    boost::function< double( ) > lightnessNumber_;

    const boost::function< Eigen::Vector2d( const double ) > sailOrientationAnglesFunction_;

    double sailEfficiency_;

    Eigen::Vector2d currentSailAngles_;

    double currentEffectiveSailEfficiency_;

    double currentLightnessNumber_;



};

} // namespace electro_magnetism
} // namespace tudat

#endif // TUDAT_SOLARSAILRADIATIONPRESSUREINTERFACE_H
