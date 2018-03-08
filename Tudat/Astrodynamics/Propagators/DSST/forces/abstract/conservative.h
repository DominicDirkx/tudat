#ifndef TUDAT_PROPAGATORS_DSST_FORCEMODELS_CONSERVATIVE_H
#define TUDAT_PROPAGATORS_DSST_FORCEMODELS_CONSERVATIVE_H

#include "forceModel.h"
#include "Tudat/Astrodynamics/Propagators/DSST/utilities/upperBounds.h"

namespace tudat
{

namespace propagators
{

namespace sst
{

namespace force_models
{


//! Struct for storing settings for Conservative ForceModel
struct ConservativeSettings: ForceModelSettings {
    //! Default constructor
    ConservativeSettings( unsigned int maxS = 0, unsigned int maxN = 0 )
        : ForceModelSettings( ), maxS( maxS ), maxN( maxN ) { }

    //! Maximum value for the outer summatory of the serie expansion (max s)
    //! If set to 0, the force model may try to determine the right value based on orbital characteristics
    unsigned int maxS;

    //! Maximum value for the inner summatory of the serie expansion (max n)
    //! If set to 0, the force model may try to determine the right value based on orbital characteristics
    unsigned int maxN;
};


//! Abstract class for perturbations that can be expressed as a disturbing potential
class Conservative : public virtual ForceModel {
public:
    boost::shared_ptr< ConservativeSettings > getSettings( ) {
        boost::shared_ptr< ConservativeSettings > castedSettings =
                boost::dynamic_pointer_cast< ConservativeSettings >( settings );
        return castedSettings != NULL ? castedSettings : boost::make_shared< ConservativeSettings >( );
    }

protected:

    Conservative( AuxiliaryElements &auxiliaryElements, boost::shared_ptr< ConservativeSettings > settings = NULL ) :
        ForceModel( auxiliaryElements, settings  ),
        VnsFactory( auxiliaryElements.VnsFactory ) { }


    //! Reference to auxiliaryElement's Vns coefficients generator, which is shared amongst several perturbations.
    coefficients_factories::VnsCoefficientsFactory &VnsFactory;

    // Common factors for potential computation
    //! B²
    double B2;
    //! B³
    double B3;

    //! \Chi = 1 / sqrt(1 - e²) = 1 / B
    double Chi;
    //! \Chi²
    double Chi2;
    //! \Chi³
    double Chi3;
    //! -2 * a / A
    double m2aoA;
    //! B / A
    double BoA;
    //! 1 / (A * B)
    double ooAB;
    //! -C / (2 * A * B)
    double mCo2AB;
    //! B / A(1 + B)
    double BoABpo;
    //! h * \Chi³
    double hChi3;
    //! k * \Chi³
    double kChi3;

    //! Update instance's members that are computed from the current auxiliary elements.
    virtual void updateMembers( );


    //! The current value of the mean disturbing potential. Needed for the short periodic contribution.
    // double U;

    //! Maximum value of n for the series expansion of the disturbing potential.
    unsigned int N;

    //! Maximum value of s for the series expansion of the disturbing potential.
    unsigned int S;

    //! Set up the force model.
    virtual void setUp() {
        updateMembers();
        determineTruncationValues();
    }


private:

    //! Set the values of N and S for the series expansion of the disturbing potential.
    virtual void determineTruncationValues( ) = 0;

    //! Function to compute the U derivatives for the current state to be implemented by derived classes.
    //! Must update the value of the member U, used later for the computation of short period terms.
    virtual Eigen::Vector6d computeMeanDisturbingFunctionPartialDerivatives( ) = 0;

    //! Get the mean element rates for the current auxiliary elements [ Eq. 3.1-(1) ]
    Eigen::Vector6d computeMeanElementRates( );


};



} // namespace force_models

} // namespace sst

} // namespace propagators

} // namespace tudat

#endif // TUDAT_PROPAGATORS_DSST_FORCEMODELS_CONSERVATIVE_H
