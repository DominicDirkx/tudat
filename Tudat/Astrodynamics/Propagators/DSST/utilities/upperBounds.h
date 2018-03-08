#ifndef TUDAT_PROPAGATORS_DSST_UPPERBOUNDS_H
#define TUDAT_PROPAGATORS_DSST_UPPERBOUNDS_H

#include "Tudat/Astrodynamics/Propagators/DSST/utilities/vectors.h"


namespace tudat
{

namespace propagators
{

namespace sst
{

namespace upper_bounds
{


//! Get the upper bound value D^n_l(\Chi).
/** Get the upper bound value D^n_l(\Chi).
 *
 * \param xx value of \Chi²
 * \param xpl value of \Chi * (\Chi² / 2)<sup>l</sup>
 * \param n index n (power of a/R)
 * \param l index l (power of eccentricity)
 * \return the upper bound D^n_l(\Chi)
 */
double getDnl( const double xx, const double xpl, const int n, const int l );


//! Get the upper bound value R^ε_{n,m,l}(γ).
/** Get the upper bound value R^ε_{n,m,l}(γ).
 *
 * \param gamma value of γ
 * \param n index n
 * \param l index l
 * \param m index m
 * \param eps ε value (+1/-1)
 * \param irf retrograde factor I (+1/-1)
 * \return the upper bound R^ε_{n,m,l}(γ)
 */
double getRnml(const double gamma, const int n, const int l, const int m, const int eps, const int irf);


} // namespace upper_bounds

} // namespace sst

} // namespace propagators

} // namespace tudat

#endif // TUDAT_PROPAGATORS_DSST_UPPERBOUNDS_H
