/*    Copyright (c) 2010-2017, Delft University of Technology
 *    All rigths reserved
 *
 *    This file is part of the Tudat. Redistribution and use in source and
 *    binary forms, with or without modification, are permitted exclusively
 *    under the terms of the Modified BSD license. You should have received
 *    a copy of the license with this file. If not, please or visit:
 *    http://tudat.tudelft.nl/LICENSE.
 */

#ifndef TUDAT_CREATEROTATIONMODEL_H
#define TUDAT_CREATEROTATIONMODEL_H

#include <Eigen/Core>
#include <Eigen/Geometry>

#include <boost/shared_ptr.hpp>

#include "Tudat/SimulationSetup/EnvironmentSetup/body.h"
#include "Tudat/Astrodynamics/Ephemerides/rotationalEphemeris.h"


namespace tudat
{

namespace simulation_setup
{

//! List of rotation models available in simulations
/*!
 *  List of rotation models available in simulations. Rotation models not defined by this
 *  given enum cannot be used for automatic model setup.
 */
enum RotationModelType
{
    simple_rotation_model,
    spice_rotation_model,
    direct_rotation_variation_model
};

//! Class for providing settings for rotation model.
/*!
 *  Class for providing settings for automatic rotation model creation. This class is a
 *  functional (base) class for settings of rotation models that require no information in
 *  addition to their type. Rotation model classes defining requiring additional information
 *  must be created using an object derived from this class.
 */
class RotationModelSettings
{
public:

    //! Constructor, sets type of rotation model.
    /*!
     *  Constructor, sets type of rotation model and base and target frame identifiers.
     *  Settings for rotation models requiring additional information should be defined in a
     *  derived class.
     *  \param rotationType Type of rotation model that is to be created.
     *  \param originalFrame Base frame of rotation model.
     *  \param targetFrame Target frame of rotation model.
     */
    RotationModelSettings( const RotationModelType rotationType,
                           const std::string& originalFrame,
                           const std::string& targetFrame ):
        rotationType_( rotationType ), originalFrame_( originalFrame ),
        targetFrame_( targetFrame ){ }

    //! Destructor.
    virtual ~RotationModelSettings( ){ }

    //! Function to return the type of rotation model that is to be created.
    /*!
     *  Function to return the type of rotation model that is to be created.
     *  \return Type of rotation model that is to be created.
     */
    RotationModelType getRotationType( ){ return rotationType_; }

    //! Function to return the base frame of rotation model.
    /*!
     *  Function to return the base frame of rotation model.
     *  \return Base frame of rotation model.
     */
    std::string getOriginalFrame( ){ return originalFrame_; }

    //! Function to return the target frame of rotation model.
    /*!
     *  Function to return the target frame of rotation model.
     *  \return Target frame of rotation model.
     */
    std::string getTargetFrame( ){ return targetFrame_; }


    //! Function to rese the orientation of the base frame.
    /*!
     * Function to reset the orientation of the base frame.
     * \param originalFrame New base frame orientation
     */
    void resetOriginalFrame( const std::string& originalFrame )
    {
        originalFrame_ = originalFrame;
    }

protected:

    //! Type of rotation model that is to be created.
    RotationModelType rotationType_;

    //! Target frame of rotation model.
    std::string originalFrame_;

    //! Base frame of rotation model.
    std::string targetFrame_;

};

//! RotationModelSettings derived class for defining settings of a simple rotational ephemeris.
class SimpleRotationModelSettings: public RotationModelSettings
{
public:
    //! Constructor,
    /*!
     *  Constructor, sets simple rotational ephemeris properties.
     *  \param originalFrame Base frame of rotation model.
     *  \param targetFrame Target frame of rotation model.
     *  \param initialOrientation Rotation from base to target frame at initialTime.
     *  \param initialTime Time at which initialOrientation represents the instantaneous rotation.
     *  \param rotationRate Rotation rate of body about its local z-axis.
     */
    SimpleRotationModelSettings( const std::string& originalFrame,
                                 const std::string& targetFrame,
                                 const Eigen::Quaterniond& initialOrientation,
                                 const double initialTime,
                                 const double rotationRate ):
        RotationModelSettings( simple_rotation_model, originalFrame, targetFrame ),
        initialOrientation_( initialOrientation ),
        initialTime_( initialTime ), rotationRate_( rotationRate ){ }

    //! Function to return rotation from base to target frame at initialTime.
    /*!
     *  Function to return rotation from base to target frame at initialTime.
     *  \return Rotation from base to target frame at initialTime.
     */
    Eigen::Quaterniond getInitialOrientation( ){ return initialOrientation_; }

    //! Function to return time at which initialOrientation represents the instantaneous rotation.
    /*!
     *  Function to return time at which initialOrientation represents the instantaneous rotation.
     *  \return Time at which initialOrientation represents the instantaneous rotation.
     */
    double getInitialTime( ){ return initialTime_; }

    //! Function to return rotation rate of body about its local z-axis.
    /*!
     *  Function to return rotation rate of body about its local z-axis.
     *  \return Rotation rate of body about its local z-axis.
     */
    double getRotationRate( ){ return rotationRate_; }

private:

    //!  Rotation from base to target frame at initialTime.
    Eigen::Quaterniond initialOrientation_;

    //! Time at which initialOrientation represents the instantaneous rotation.
    double initialTime_;

    //! Rotation rate of body about its local z-axis.
    double rotationRate_;
};

class DirectRotationVariationSettings: public RotationModelSettings
{
public:
    DirectRotationVariationSettings( const std::vector< double > rightAscensionPolynomialTerms,
                                     const std::vector< double > declinationPolynomialTerms,
                                     const std::vector< double > primeMeridianPolynomialTerms,
                                     const std::map< double, std::pair< double, double > > rightAscensionLibrations,
                                     const std::map< double, std::pair< double, double > > declinationLibrations,
                                     const std::map< double, std::pair< double, double > > primeMeridianLibrations,
                                     const Eigen::Matrix3d fromIntermediateFrameToBaseFrame,
                                     std::string originalFrame,
                                     std::string targetFrame ):
        RotationModelSettings( direct_rotation_variation_model, originalFrame, targetFrame ),
        rightAscensionPolynomialTerms_( rightAscensionPolynomialTerms ), declinationPolynomialTerms_( declinationPolynomialTerms ),
        primeMeridianPolynomialTerms_( primeMeridianPolynomialTerms) , rightAscensionLibrations_( rightAscensionLibrations ),
        declinationLibrations_( declinationLibrations ), primeMeridianLibrations_( primeMeridianLibrations ),
        fromIntermediateFrameToBaseFrame_( fromIntermediateFrameToBaseFrame ){ }


    std::vector< double > getRightAscensionPolynomialTerms( ){ return rightAscensionPolynomialTerms_; }
    std::vector< double > getDeclinationPolynomialTerms_( ){ return declinationPolynomialTerms_; }
    std::vector< double > getPrimeMeridianPolynomialTerms( ){ return primeMeridianPolynomialTerms_; }
    std::map< double, std::pair< double, double > > getRightAscensionLibrations( ){ return rightAscensionLibrations_; }
    std::map< double, std::pair< double, double > > getDeclinationLibrations( ){ return declinationLibrations_; }
    std::map< double, std::pair< double, double > > getPrimeMeridianLibrations( ){ return primeMeridianLibrations_; }

    Eigen::Matrix3d getFromIntermediateFrameToBaseFrame( ){ return fromIntermediateFrameToBaseFrame_; }

private:

    std::vector< double > rightAscensionPolynomialTerms_;
    std::vector< double > declinationPolynomialTerms_;
    std::vector< double > primeMeridianPolynomialTerms_;
    std::map< double, std::pair< double, double > > rightAscensionLibrations_;
    std::map< double, std::pair< double, double > > declinationLibrations_;
    std::map< double, std::pair< double, double > > primeMeridianLibrations_;

    Eigen::Matrix3d fromIntermediateFrameToBaseFrame_;
};


//! Function to create a rotation model.
/*!
 *  Function to create a rotation model based on model-specific settings for the rotation.
 *  \param rotationModelSettings Settings for the rotation model that is to be created, defined
 *  a pointer to an object of class (derived from) RotationSettings.
 *  \param body Name of the body for which the rotation model is to be created.
 *  \return Rotation model created according to settings in rotationModelSettings.
 */
boost::shared_ptr< ephemerides::RotationalEphemeris > createRotationModel(
        const boost::shared_ptr< RotationModelSettings > rotationModelSettings,
        const std::string& body );
} // namespace simulation_setup

} // namespace tudat

#endif // TUDAT_CREATEROTATIONMODEL_H
