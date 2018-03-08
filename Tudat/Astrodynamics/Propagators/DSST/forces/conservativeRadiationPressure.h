#ifndef TUDAT_DSST_PROPAGATORS_FORCEMODELS_CONSERVATIVESOLARRADIATIONPRESSURE_H
#define TUDAT_DSST_PROPAGATORS_FORCEMODELS_CONSERVATIVESOLARRADIATIONPRESSURE_H

#include "Tudat/Astrodynamics/BasicAstrodynamics/physicalConstants.h"

#include "Tudat/Astrodynamics/ElectroMagnetism/cannonBallRadiationPressureAcceleration.h"
#include "Tudat/Astrodynamics/ElectroMagnetism/radiationPressureInterface.h"

// #include "Tudat/Astrodynamics/Propagators/DSST/dsst.h"

#include "abstract/conservativeThirdBodyPertubed.h"


namespace tudat
{

namespace propagators
{

namespace sst
{

namespace force_models
{


typedef boost::shared_ptr< electro_magnetism::CannonBallRadiationPressureAcceleration > RadiationPressureAM;


//! Final class for the contribution of a third body's central gravity
class ConservativeRadiationPressure final : public ConservativeThirdBodyPerturbed
{
public:

    //! Constructor.
    /** Simple constructor.
     * \param auxiliaryElements Auxiliary elements used to compute the mean element rates and short period terms.
     * \param thirdBody The third body exerting the acceleration.
     */
    ConservativeRadiationPressure( AuxiliaryElements &auxiliaryElements, RadiationPressureAM radiationPressureAM,
                                   boost::shared_ptr< ConservativeSettings > settings = NULL ) :
        ForceModel( auxiliaryElements, settings ),
        ConservativeThirdBodyPerturbed( auxiliaryElements, settings ),
        radiationPressureAM( radiationPressureAM ) { }


private:

    //! Set up the force model.
    void setUp() {
        radiationPressureAM->updateMembers( aux.epoch );
        Conservative::setUp();
    }

    //! Update instance's members that are computed from the current auxiliary elements.
    void updateMembers();

    //! Pointer to the radiation pressure acceleration model
    RadiationPressureAM radiationPressureAM;

    //! Function that updates `Ufactor`
    double computeUfactor();

    //! Minimum value for n in the expansion
    unsigned int nMin( unsigned int s ) {
        return std::max( 1, (int) s );
    }

};


} // namespace force_models

} // namespace sst

} // namespace propagators

} // namespace tudat

#endif // TUDAT_DSST_PROPAGATORS_FORCEMODELS_CONSERVATIVERADIATIONPRESSURE_H
