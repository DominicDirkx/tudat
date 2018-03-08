#ifndef TUDAT_PROPAGATORS_DSST_FORCEMODELS_THIRDBODYPERTURBED_H
#define TUDAT_PROPAGATORS_DSST_FORCEMODELS_THIRDBODYPERTURBED_H

#include "forceModel.h"

namespace tudat
{

namespace propagators
{

namespace sst
{

namespace force_models
{


//! Abstract class for perturbations cuased by third bodies
class ThirdBodyPerturbed : public virtual ForceModel {
protected:

    //! Constructor
    ThirdBodyPerturbed( AuxiliaryElements &auxiliaryElements ) : ForceModel( auxiliaryElements ) { }


    //! Set up the force model.
    virtual void setUp() {
        updateMembers();
    }

    //! Update instance's members that are computed from the current auxiliary elements.
    virtual void updateMembers( );


    /// Derived classes are responsible for updating `r3` within their implementations of updateMembers().
    /// Then, they need to call ThirdBodyPerturbed::updateMembers(), which computes R3 and the direction cosines.

    //! Position of the third body from the centre of mass of the central body
    Eigen::Vector3d r3 = Eigen::Vector3d::Zero();

    //! Distance from center of mass of the central body to the 3rd body
    double R3;

};



} // namespace force_models

} // namespace sst

} // namespace propagators

} // namespace tudat

#endif // TUDAT_PROPAGATORS_DSST_FORCEMODELS_THIRDBODYPERTURBED_H
