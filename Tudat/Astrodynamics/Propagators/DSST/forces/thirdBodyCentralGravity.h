#ifndef TUDAT_PROPAGATORS_DSST_FORCEMODELS_THIRDBODYCENTRALGRAVITY_H
#define TUDAT_PROPAGATORS_DSST_FORCEMODELS_THIRDBODYCENTRALGRAVITY_H

#include "Tudat/Astrodynamics/Gravitation/thirdBodyPerturbation.h"

// #include "Tudat/Astrodynamics/Propagators/DSST/vectors.h"

#include "abstract/conservativeThirdBodyPertubed.h"


namespace tudat
{

namespace propagators
{

namespace sst
{

namespace force_models
{


typedef boost::shared_ptr< gravitation::ThirdBodyCentralGravityAcceleration > ThirdBodyAM;


//! Final class for the contribution of a third body's central gravity
class ThirdBodyCentralGravity final : public ConservativeThirdBodyPerturbed
{
public:

    //! Constructor.
    /** Simple constructor.
     * \param auxiliaryElements Auxiliary elements used to compute the mean element rates and short period terms.
     * \param thirdBody The third body exerting the acceleration.
     */
    ThirdBodyCentralGravity( AuxiliaryElements &auxiliaryElements, ThirdBodyAM thirdBodyAM,
                             boost::shared_ptr< ConservativeSettings > settings = NULL ) :
        ForceModel( auxiliaryElements, settings ),
        ConservativeThirdBodyPerturbed( auxiliaryElements, settings ),
        thirdBodyAM( thirdBodyAM ) { }


private:

    //! Set up the force model.
    void setUp() {
        thirdBodyAM->updateMembers( aux.epoch );
        Conservative::setUp();
    }

    //! Pointer to thrid body acceleration model
    ThirdBodyAM thirdBodyAM;

    //! Update instance's members that are computed from the current auxiliary elements.
    void updateMembers( );

    //! Standard gravitational parameter μ for the third body in m³/s²
    double mu3;

    //! Function that updates `Ufactor`
    double computeUfactor() {
        return mu3 / R3;
    }

    //! Minimum value for n in the expansion
    unsigned int nMin( unsigned int s ) {
        return std::max( 2, (int) s );
    }


    /*

    ///////////////////////////////////////////////////////////////////////////////
    ////////////////////////////// Short period terms /////////////////////////////
    ///////////////////////////////////////////////////////////////////////////////

    //! e / ( 1 + sqrt(1 - e²) )
     double c;

    //! Map of Jacobi polynomials factories used by getwCoefficient
    std::map< coefficients_factories::JacobiKey, coefficients_factories::JacobiPolynomialsFactory > jacobiFactories;

    //! Function to obtain the coefficients e^{-|j-s|} * w^{n,s}_j(e) and e^{-(j+s)} * w^{n,s}_{-j}(e) [Eq. 4.2-(16,10)]
    std::pair< double, double > getewjCoefficients( const int n, const int s, const int j );

    */
};


} // namespace force_models

} // namespace sst

} // namespace propagators

} // namespace tudat

#endif // TUDAT_PROPAGATORS_DSST_FORCEMODELS_THIRDBODYCENTRALGRAVITY_H
