#ifndef TUDAT_PROPAGATORS_DSST_FORCEMODELS_NONCONSERVATIVE_H
#define TUDAT_PROPAGATORS_DSST_FORCEMODELS_NONCONSERVATIVE_H

#include "Tudat/Astrodynamics/BasicAstrodynamics/accelerationModel.h"
#include "Tudat/SimulationSetup/PropagationSetup/environmentUpdater.h"

#include "forceModel.h"

namespace tudat
{

namespace propagators
{

namespace sst
{

namespace force_models
{


typedef boost::shared_ptr< basic_astrodynamics::AccelerationModel< Eigen::Vector3d > > AccelerationModel;


//! Struct for storing settings for NonConservative ForceModel
struct NonConservativeSettings: ForceModelSettings {
    //! Default constructor
    NonConservativeSettings( unsigned int numberOfQuadratureNodes   = 40,
                             bool numberOfQuadratureNodesIsScalable = false,
                             bool alwaysIncludeCentralNode          = false )
        : ForceModelSettings( ),
          numberOfQuadratureNodes( numberOfQuadratureNodes ),
          numberOfQuadratureNodesIsScalable( numberOfQuadratureNodesIsScalable ),
          alwaysIncludeCentralNode( alwaysIncludeCentralNode ) { }

    //! Number of quadrature nodes to be used for the Gaussian quadrature in the averaging integral
    unsigned int numberOfQuadratureNodes;

    //! Whether the actual number of quadrature nodes used is scaled based on `numberOfQuadratureNodes` and the
    //! integral limits (L2 - L1) comapred to 2π
    bool numberOfQuadratureNodesIsScalable;

    //! Wether the central node (halfway between L1 and L2) should always be evaluated.
    //! When this is set to true, the actual number of nodes is always chosen odd (increasing it by 1 if necessary).
    bool alwaysIncludeCentralNode;
};


//! Abstract class for perturbations that can be expressed as a disturbing potential
class NonConservative : public virtual ForceModel {
public:

    std::string getPerturbingBodyName() {
        return perturbingBody;
    }

    AccelerationModel getAccelerationModel() {
        return accelerationModel;
    }

    void setEnvironmentUpdater( boost::shared_ptr< EnvironmentUpdater< double, double > > updater ) {
        environmentUpdater = updater;
    }

    boost::shared_ptr< NonConservativeSettings > getSettings( ) {
        boost::shared_ptr< NonConservativeSettings > castedSettings =
                boost::dynamic_pointer_cast< NonConservativeSettings >( settings );
        return castedSettings != NULL ? castedSettings : boost::make_shared< NonConservativeSettings >( );
    }


protected:

    //! Derived (non-abstract) classes must call this constructor with either
    //! maximumScalableNumberOfQuadratureNodes or fixedNumberOfQuadratureNodes ≥ 3
    NonConservative( AuxiliaryElements &auxiliaryElements, const std::string &perturbingBody,
                     AccelerationModel accelerationModel,
                     boost::shared_ptr< NonConservativeSettings > settings = NULL ) :
        ForceModel( auxiliaryElements, settings ),
        perturbingBody( perturbingBody ),
        accelerationModel( accelerationModel ) { }


    //! Name of the body exerting the acceleration
    const std::string perturbingBody;

    //! Pointer to the corresponding acceleration model
    AccelerationModel accelerationModel;

    //! Environment updater
    boost::shared_ptr< EnvironmentUpdater< double, double > > environmentUpdater;

    //! Set up the force model.
    virtual void setUp() {
        updateMembers();
    }

    //! Update instance's members that are computed from the current auxiliary elements.
    virtual void updateMembers( );


    // Common factors for integral computation
    //! b = 1 / ( 1 + B )
    // double b;


    //! Lower limit for the true longitude in the averaging integral
    double L1 = -mathematical_constants::PI;

    //! Upper limit for the true longitude in the averaging integral ( L2 - L1 ≤ 2π )
    double L2 = mathematical_constants::PI;

    //! Number of nodes for the quadrature of the averaging integral
    unsigned int N;

    //! The value of the true longitude during the current integration step
    double L;

    //! The epoch during the current integration step
    double epoch;

    //! The propagated body distance |x, y, z| during the current integration step
    double R;

    //! The propagated body speed |vx, vy, vz| during the current integration step
    // double V;

    //! The x components in the equinoctial reference frame during the current integration step
    double X;

    //! The y components in the equinoctial reference frame during the current integration step
    double Y;

    //! The vx components in the equinoctial reference frame during the current integration step
    double Xdot;

    //! The vy components in the equinoctial reference frame during the current integration step
    double Ydot;

    //! The propagated body position (x, y, z) during the current integration step
    Eigen::Vector3d r;

    //! The propagated body velocity (vx, vy, vz) during the current integration step
    Eigen::Vector3d v;

    //! Partial derivatives of the equinoctial elements wrt v (/dvx, /dvy, /dvz) during the current integration step
    Eigen::Matrix< double, 3, 6 > partials;

    //! The perturbing acceleration q (qx, qy, qz) during the current integration step
    Eigen::Vector3d perturbingAcceleration;


protected:

    //! Update the values of the minimum and maximum true longitude for the averaging integral.
    virtual void determineIntegrationLimits( ) = 0;

    //! Get the value of the integrand for a given `trueLongitude` (for each of the equinoctial elements)
    Eigen::Vector6d integrand( const double trueLongitude );

    //! Update the acceleration model to the current epoch and get the perturbing acceleration
    virtual Eigen::Vector3d getPerturbingAcceleration();

    //! Get the mean element rates for the current auxiliary elements [ Eq. 3.1-(1) ]
    Eigen::Vector6d computeMeanElementRates( );

    //! Returns the short period terms.
    Eigen::Vector6d computeShortPeriodTerms( );

};



} // namespace force_models

} // namespace sst

} // namespace propagators

} // namespace tudat

#endif // TUDAT_PROPAGATORS_DSST_FORCEMODELS_NONCONSERVATIVE_H
