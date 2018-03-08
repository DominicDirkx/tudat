#ifndef TUDAT_PROPAGATORS_DSST_OSCULATINGTOMEAN_H
#define TUDAT_PROPAGATORS_DSST_OSCULATINGTOMEAN_H

#include "Tudat/Astrodynamics/Propagators/DSST/utilities/vectors.h"
#include "Tudat/Astrodynamics/Propagators/DSST/forces/abstract/forceModel.h"


namespace tudat
{

namespace propagators
{

namespace sst
{

namespace element_conversions
{

//! Transform osculating to mean elements.
/*!
 * Transforms the equinoctial elements contained by osculatingAuxiliaryElements from osculating to mean.
 * Based on a fixed-point interative process from Section 7.4-2.
 * \param auxiliaryElements Auxiliary elements containing the osculating equinoctial elements.
 * Passed by reference and modified by this function so that it ends up containing the mean equioctial elements.
 * \param forceModels The force models generating the short period terms needed to do the transformation.
 */
void transformOsculatingToMeanElements( AuxiliaryElements &auxiliaryElements,
                                        force_models::ForceModelMap forceModels );


} // namespace element_conversions

} // namespace sst

} // namespace propagators

} // namespace tudat

#endif // TUDAT_PROPAGATORS_DSST_OSCULATINGTOMEAN_H
