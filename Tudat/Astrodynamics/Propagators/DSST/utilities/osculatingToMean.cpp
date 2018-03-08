/*    Copyright (c) 2010-2017, Delft University of Technology
 *    All rigths reserved
 *
 *    This file is part of the Tudat. Redistribution and use in source and
 *    binary forms, with or without modification, are permitted exclusively
 *    under the terms of the Modified BSD license. You should have received
 *    a copy of the license with this file. If not, please or visit:
 *    http://tudat.tudelft.nl/LICENSE.
 */

#include "osculatingToMean.h"
#include "Tudat/Astrodynamics/Propagators/DSST/forces/zonalSphericalHarmonicGravity.h"

namespace tudat
{

namespace propagators
{

namespace sst
{

namespace element_conversions
{

//! Transform osculating to mean elements.
void transformOsculatingToMeanElements( AuxiliaryElements& auxiliaryElements,
                                        force_models::ForceModelMap forceModels )
{
    using namespace Eigen;
    using namespace mathematical_constants;
    using namespace orbital_element_conversions;

    // Reference to equinoctial elements
    EquinoctialElements& equinoctialElements = auxiliaryElements.equinoctialElements;

    // Store initial osculating elements
    Vector6d initialOsc = equinoctialElements.getComponents( meanType );
    // std::cout << "Initial osculating: " << initialOsc.transpose() << std::endl;

    // threshold for each parameter
    Vector6d keplerian = equinoctialElements.toKeplerian();
    const double epsilon    = 1e-10;
    const double thresholdA = epsilon * ( 1 + std::fabs( keplerian( semiMajorAxisIndex ) ) );
    const double thresholdE = epsilon * ( 1 + keplerian( eccentricityIndex ) );
    const double thresholdI = epsilon * ( 1 + keplerian( inclinationIndex ) );
    const double thresholdL = epsilon * PI;
    Vector6d thresholds;
    thresholds << thresholdA, thresholdE, thresholdE, thresholdI, thresholdI, thresholdL;

    unsigned int i = 0;
    const unsigned int maxIterations = 200;
    while ( i++ < maxIterations ) {
        // Initialize rebuilt to be equal to the current mean elements
        Vector6d rebuilt = equinoctialElements.getComponents( meanType );

        // Add short period terms to rebuilt to make it osculating
        for ( auto ent: forceModels ) {
            boost::shared_ptr< sst::force_models::ForceModel > forceModel = ent.second;
            // FIXME: only J2, calling this on ConservativeThirdBodyPerturbed breaks acceleration model
            /*
            if ( boost::dynamic_pointer_cast< sst::force_models::ZonalSphericalHarmonicGravity >( forceModel ) != NULL )
            {
                rebuilt += forceModel->getShortPeriodTerms();
            }
            */
            rebuilt += forceModel->getShortPeriodTerms();
        }

        // Get difference
        Vector6d delta = initialOsc - rebuilt;

        // Convert last element to range [-π, π]
        delta( fastVariableIndex ) -= 2 * PI * std::floor( 0.5 * ( delta( fastVariableIndex ) / PI + 1 ) );
        // std::cout << i << ": " << delta.transpose() << std::endl;

        // Check convergence
        if ( ( delta.array().abs() < thresholds.array() ).all() ) {
            return;
        }

        // Update auxiliaryElements' members

        // std::cout << "equ before: " << auxiliaryElements.equinoctialElements.getComponents( meanType ).transpose() << std::endl;
        // std::cout << "before: " << auxiliaryElements.equinoctialElements.toKeplerian().transpose() << std::endl;

        auxiliaryElements.updateComponentsByAdding( delta );

        // std::cout << "after:  " << auxiliaryElements.equinoctialElements.toKeplerian().transpose() << std::endl;
        // std::cout << "equ after: " << auxiliaryElements.equinoctialElements.getComponents( meanType ).transpose() << std::endl;

        return;  // FIXME: currently only one iteration; more iterations -> diverges
    }

    throw std::runtime_error( "Impossible to transform osculating to mean elements: "
                              "maximum number of iterations (200) exceeded." );
}


} // namespace element_conversions

} // namespace sst

} // namespace propagators

} // namespace tudat
