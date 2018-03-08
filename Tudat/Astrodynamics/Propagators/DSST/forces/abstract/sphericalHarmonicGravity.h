#ifndef TUDAT_PROPAGATORS_DSST_FORCEMODELS_SPHERICHARMONICGRAVITY_H
#define TUDAT_PROPAGATORS_DSST_FORCEMODELS_SPHERICHARMONICGRAVITY_H

#include "Tudat/Astrodynamics/Gravitation/sphericalHarmonicsGravityModel.h"

#include "Tudat/Mathematics/BasicMathematics/legendrePolynomials.h"

#include "centralBodyPertubed.h"
#include "conservative.h"

#include "Tudat/Astrodynamics/Propagators/DSST/utilities/coefficientsFactories.h"
#include "Tudat/Astrodynamics/Propagators/DSST/utilities/upperBounds.h"


namespace tudat
{

namespace propagators
{

namespace sst
{

namespace force_models
{


typedef boost::shared_ptr< gravitation::SphericalHarmonicsGravitationalAccelerationModel > SphericalHarmonicsAM;


//! Abstract class for gravity perturbations expressed as a spheric harmonic expansion.
//! Zonal and tesseral terms are treated differently, so derived classes will implement these terms separately.
class SphericalHarmonicGravity :
        public CentralBodyPerturbed, public Conservative
{
protected:

    //! Constructor.
    /** Simple constructor.
     * \param auxiliaryElements Auxiliary elements used to compute the mean element rates and short period terms.
     */
    SphericalHarmonicGravity( AuxiliaryElements &auxiliaryElements, SphericalHarmonicsAM sphericalHarmonicsAM ) :
        ForceModel( auxiliaryElements ),
        CentralBodyPerturbed( auxiliaryElements ),
        Conservative( auxiliaryElements ),
        sphericalHarmonicsAM( sphericalHarmonicsAM ),
        R( sphericalHarmonicsAM->getReferenceRadius() )
    {
        using namespace Eigen;
        using namespace basic_mathematics;

        // Get Cterms, Sterms
        MatrixXd cosCoefficients = sphericalHarmonicsAM->getCosineHarmonicCoefficientsFunction()();
        MatrixXd sinCoefficients = sphericalHarmonicsAM->getSineHarmonicCoefficientsFunction()();

        shDegree = cosCoefficients.rows() - 1;
        shOrder  = cosCoefficients.cols() - 1;

        Cterms.resize( shDegree - 1 );
        Sterms.resize( shDegree - 1 );
        for ( unsigned int n = 2; n <= shDegree; n++ ) {
            const unsigned int currentOrder = std::min( n, shOrder );
            std::vector< double > c( currentOrder + 1 );
            std::vector< double > s( currentOrder + 1 );
            for ( unsigned int m = 0; m <= currentOrder; m++ ) {
                const double factor = calculateLegendreGeodesyNormalizationFactor( n, m );
                c[ m ] = factor * cosCoefficients( n, m );
                s[ m ] = factor * sinCoefficients( n, m );
            }
            Cterms[ n - 2 ] = c;
            Sterms[ n - 2 ] = s;
        }
    }


private:

    void setUp() {
        Conservative::setUp();
    }

    //! Update instance's members based on current auxiliary elements.
    void updateMembers( );

    //! Set the values of N and S for the series expansion of the disturbing potential.
    void determineTruncationValues( );


protected:

    //! Pointer to shperical harmonics acceleration model
    SphericalHarmonicsAM sphericalHarmonicsAM;

    //! Reference radius used for the geopotential expansion
    const double R;

    //! Cosine coefficients
    NestedVectord Cterms;

    //! Since coefficients
    NestedVectord Sterms;

    //! Returns J_{n,m}
    double J( unsigned int n, unsigned int m ) {
        return std::sqrt( std::pow( Cterms( n - 2, m ), 2 ) + std::pow( Sterms( n - 2, m ), 2 ) );
    }

    //! Returns J_n = J_{n,0}
    double J( unsigned int n ) {
        return -Cterms( n - 2, 0 );
    }

    //! Degree and order of the spherical harmonic expansion
    unsigned int shDegree;
    unsigned int shOrder;


    //! Truncation tolerance
    const double TRUNCATION_TOLERANCE = 1e-4;

    //! Number of points for interpolation
    // const int INTERPOLATION_POINTS = 3;

    //! Maximum power for eccentricity used in short periodic computation
    // const int MAX_ECCPOWER_SP = 4;

    //! Qns coefficients generator
    coefficients_factories::QnsCoefficientsFactory QnsFactory;

    //! Gs coefficients generator
    coefficients_factories::GsHsCoefficientsFactory GsFactory;

    //! R/a up to the Nth power
    Vectord powR_a;

};


} // namespace force_models

} // namespace sst

} // namespace propagators

} // namespace tudat

#endif // TUDAT_PROPAGATORS_DSST_FORCEMODELS_SPHERICHARMONICGRAVITY_H
