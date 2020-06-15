/*    Copyright (c) 2010-2019, Delft University of Technology
 *    All rigths reserved
 *
 *    This file is part of the Tudat. Redistribution and use in source and
 *    binary forms, with or without modification, are permitted exclusively
 *    under the terms of the Modified BSD license. You should have received
 *    a copy of the license with this file. If not, please or visit:
 *    http://tudat.tudelft.nl/LICENSE.
 */

#ifndef TUDAT_LONGITUDELIBRATIONAMPLITUDE_H
#define TUDAT_LONGITUDELIBRATIONAMPLITUDE_H


#include "Tudat/Astrodynamics/OrbitDetermination/EstimatableParameters/estimatableParameter.h"
#include "Tudat/Astrodynamics/Ephemerides/synchronousRotationalEphemeris.h"

namespace tudat
{

namespace estimatable_parameters
{


class LongitudeLibrationAmplitude: public EstimatableParameter< double >
{

public:


    LongitudeLibrationAmplitude(
            const std::shared_ptr< ephemerides::SynchronousRotationalEphemeris > rotationModel,
            const std::string& associatedBody ):
        EstimatableParameter< double  >( longitude_libration_amplitude, associatedBody ),
        rotationModel_( rotationModel ) { }

    ~LongitudeLibrationAmplitude( ) { }


    double getParameterValue( )
    {
        return rotationModel_->getLibrationAmplitude( );
    }

    void setParameterValue( const double parameterValue )
    {
        rotationModel_->setLibrationAmplitude( parameterValue );
    }

    int getParameterSize( )
    {
        return 1;
    }

protected:

private:

    std::shared_ptr< ephemerides::SynchronousRotationalEphemeris > rotationModel_;
};

} // namespace estimatable_parameters

} // namespace tudat

#endif // TUDAT_LONGITUDELIBRATIONAMPLITUDE_H
