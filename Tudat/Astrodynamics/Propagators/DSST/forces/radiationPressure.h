#ifndef TUDAT_PROPAGATORS_DSST_FORCEMODELS_RADIATIONPRESSURE_H
#define TUDAT_PROPAGATORS_DSST_FORCEMODELS_RADIATIONPRESSURE_H

#include "Tudat/Astrodynamics/ElectroMagnetism/cannonBallRadiationPressureAcceleration.h"

#include "Tudat/Astrodynamics/Propagators/DSST/forces/abstract/nonConservative.h"
#include "Tudat/Astrodynamics/Propagators/DSST/forces/abstract/thirdBodyPertubed.h"
#include "Tudat/Astrodynamics/Propagators/DSST/forces/conservativeRadiationPressure.h"

namespace tudat
{

namespace propagators
{

namespace sst
{

namespace force_models
{


//! Struct for storing settings for RadiationPressure ForceModel
struct RadiationPressureSettings: NonConservativeSettings {
    //! Default constructor
    RadiationPressureSettings( unsigned int numberOfQuadratureNodes   = 12,
                               bool numberOfQuadratureNodesIsScalable = false,
                               bool alwaysIncludeCentralNode          = false )
        : NonConservativeSettings( numberOfQuadratureNodes, numberOfQuadratureNodesIsScalable,
                                   alwaysIncludeCentralNode ) { }
};


//! Abstract class for perturbations that can be expressed as a disturbing potential
class RadiationPressure final : public ThirdBodyPerturbed, public NonConservative {
public:

    RadiationPressure( AuxiliaryElements &auxiliaryElements, const std::string &perturbingBody,
                       RadiationPressureAM radiationPressureAM,
                       boost::shared_ptr< RadiationPressureSettings > settings = NULL ) :
        ForceModel( auxiliaryElements, settings ),
        ThirdBodyPerturbed( auxiliaryElements ),
        NonConservative( auxiliaryElements, perturbingBody, radiationPressureAM, settings ),
        radiationPressureAM( radiationPressureAM ) { }


private:

    //! Pointer to the radiation pressure acceleration model
    RadiationPressureAM radiationPressureAM;

    //! Set up the force model.
    void setUp() {
        updateMembers();
    }

    //! Update instance's members that are computed from the current auxiliary elements.
    virtual void updateMembers( );


    //! Update the values of the minimum and maximum true longitude for the averaging integral.
    void determineIntegrationLimits( );


};



} // namespace force_models

} // namespace sst

} // namespace propagators

} // namespace tudat

#endif // TUDAT_PROPAGATORS_DSST_FORCEMODELS_RADIATIONPRESSURE_H
