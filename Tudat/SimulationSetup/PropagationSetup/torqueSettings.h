/*    Copyright (c) 2010-2018, Delft University of Technology
 *    All rigths reserved
 *
 *    This file is part of the Tudat. Redistribution and use in source and
 *    binary forms, with or without modification, are permitted exclusively
 *    under the terms of the Modified BSD license. You should have received
 *    a copy of the license with this file. If not, please or visit:
 *    http://tudat.tudelft.nl/LICENSE.
 */

#ifndef TUDAT_TORQUESETTINGS_H
#define TUDAT_TORQUESETTINGS_H

#include <boost/tuple/tuple.hpp>

#include "Tudat/Astrodynamics/BasicAstrodynamics/torqueModelTypes.h"


namespace tudat
{

namespace simulation_setup
{

//! Class for providing settings for torque model.
/*!
 *  Class for providing settings for torque model. This class is a functional (base) class for
 *  settings of torque models that  require no information in addition to their type.
 *  Classes defining settings for torque models requiring additional information must be
 *  derived from this class.
 *  Bodies exerting and undergong torque are set externally from this class.
 *  This class can be used for the easy setup of torque models
 *  (see createTorqueModels.h), but users may also chose to do so manually.
 *  (Derived) Class members are all public, for ease of access and modification.
 */
class TorqueSettings
{
public:

    //! Constructor, sets type of torque.
    /*!
     *  Constructor, sets type of torque.
     *  \param torqueType Type of torque from AvailableTorque enum.
     */
    TorqueSettings( const basic_astrodynamics::AvailableTorque torqueType ):torqueType_( torqueType ){ }

    //! Destructor
    virtual ~TorqueSettings( ){ }

    //! Type of torque that is to be created.
    basic_astrodynamics::AvailableTorque torqueType_;

};

//! Class to define settings for a spherical harmonic gravitational torque exerted by a point mass.
class SphericalHarmonicTorqueSettings: public TorqueSettings
{
public:

    //! Constructor
    /*!
     * Constructor
     * \param maximumDegree Maximum degree to which gravity field of body undergoing torque is to be exerted
     * \param maximumOrder Maximum order to which gravity field of body undergoing torque is to be exerted
     */
    SphericalHarmonicTorqueSettings( const int maximumDegree,
                                     const int maximumOrder ):
        TorqueSettings( basic_astrodynamics::spherical_harmonic_gravitational_torque ),
        maximumDegree_( maximumDegree ), maximumOrder_( maximumOrder ){ }

    //! Maximum degree to which gravity field of body undergoing torque is to be exerted
    int maximumDegree_;

    //! Maximum order to which gravity field of body undergoing torque is to be exerted
    int maximumOrder_;
};

class MutualExtendedBodySphericalHarmonicTorqueSettings: public TorqueSettings
{
public:

    //! Constructor
    /*!
     * Constructor
     * \param maximumDegreeOfBodyUndergoingAcceleration Maximum degree of body undergoing acceleration
     * \param maximumOrderOfBodyUndergoingAcceleration Maximum order of body undergoing acceleration
     * \param maximumDegreeOfBodyExertingAcceleration Maximum degree of body exerting acceleration
     * \param maximumOrderOfBodyExertingAcceleration Maximum order of body exerting acceleration
     */
    MutualExtendedBodySphericalHarmonicTorqueSettings(
            const int maximumDegreeOfBodyUndergoingAcceleration,
            const int maximumOrderOfBodyUndergoingAcceleration,
            const int maximumDegreeOfBodyExertingAcceleration,
            const int maximumOrderOfBodyExertingAcceleration,
            const bool includePointMassInteractions = true ):
    TorqueSettings( basic_astrodynamics::mutual_extended_body_spherical_harmonic_gravitational_torque ),
    maximumDegreeOfBody1_( maximumDegreeOfBodyUndergoingAcceleration ),
    maximumDegreeOfBody2_( maximumDegreeOfBodyExertingAcceleration )
{
    for( int i = ( !includePointMassInteractions ); i <= maximumDegreeOfBodyUndergoingAcceleration; i++ )
    {
        for( int j = 0; ( j <= maximumOrderOfBodyUndergoingAcceleration && j <= i ); j++ )
        {
            for( int k = ( !includePointMassInteractions ); k <= maximumDegreeOfBodyExertingAcceleration; k++ )
            {
                for( int l = 0; ( l <= maximumDegreeOfBodyExertingAcceleration && l <= k ); l++ )
                {
                    coefficientCombinationsToUse_.push_back( boost::make_tuple( i, j, k, l ) );
                }
            }
        }
    }
}

    //! Constructor
    /*!
     * Constructor
     * \param coefficientCombinationsToUse List of combinations of degrees/orders of two bodies for full two-body potential, in
     * order: (degree body 1, order body 1, degree body 2, order body 2)
     */
    MutualExtendedBodySphericalHarmonicTorqueSettings(
           const std::vector< boost::tuple< unsigned int, unsigned int, unsigned int, unsigned int > >&
            coefficientCombinationsToUse ):
        TorqueSettings( basic_astrodynamics::mutual_extended_body_spherical_harmonic_gravitational_torque ),
        coefficientCombinationsToUse_( coefficientCombinationsToUse )
    {
        maximumDegreeOfBody1_ = 0;
        maximumDegreeOfBody2_ = 0;

        int degreeOfBody1, degreeOfBody2;
        for( unsigned int i = 0; i < coefficientCombinationsToUse_.size( ); i++ )
        {
            degreeOfBody1 = coefficientCombinationsToUse.at( i ).get< 0 >( );
            degreeOfBody2 = coefficientCombinationsToUse.at( i ).get< 2 >( );

            if( degreeOfBody1 > maximumDegreeOfBody1_ )
            {
                maximumDegreeOfBody1_ = degreeOfBody1;
            }

            if( degreeOfBody2 > maximumDegreeOfBody2_ )
            {
                maximumDegreeOfBody2_ = degreeOfBody2;
            }
        }
    }

    //! Full list of combinations of degrees/orders of two bodies for full two-body potential, in
    //! order: (degree body 1, order body 1, degree body 2, order body 2)
    std::vector< boost::tuple< unsigned int, unsigned int, unsigned int, unsigned int > > coefficientCombinationsToUse_;

    //! Maximum degree that is used for body 1
    int maximumDegreeOfBody1_;

    //! Maximum degree that is used for body 2
    int maximumDegreeOfBody2_;
};

typedef std::map< std::string, std::map< std::string, std::vector< boost::shared_ptr< TorqueSettings > > > > SelectedTorqueMap;

} // namespace simulation_setup

} // namespace tudat

#endif // TUDAT_TORQUESETTINGS_H
