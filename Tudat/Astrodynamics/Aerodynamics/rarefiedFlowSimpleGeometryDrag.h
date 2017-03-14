/*    Copyright (c) 2010-2017, Delft University of Technology
 *    All rigths reserved
 *
 *    This file is part of the Tudat. Redistribution and use in source and
 *    binary forms, with or without modification, are permitted exclusively
 *    under the terms of the Modified BSD license. You should have received
 *    a copy of the license with this file. If not, please or visit:
 *    http://tudat.tudelft.nl/LICENSE.
 *
 *    References
 *      Anderson Jr., J.D. , Fundamentals of Aerodynamics, 3rd edition, McGraw Hill, 2001.
 *      Gentry, A., Smyth, D., and Oliver, W. . The Mark IV Supersonic-Hypersonic Arbitrary Body
 *          Program, Volume II - Program Formulation, Douglas Aircraft Company, 1973.
 *      Anderson Jr. , J.D, Hypersonic and High-Temperature Gas Dynamics, 2nd edition,
 *          AIAA Education Series, 2006.
 *
 */

#ifndef TUDAT_RAREFIEDFLOWSIMPLEGEOMETRYDRAG_H
#define TUDAT_RAREFIEDFLOWSIMPLEGEOMETRYDRAG_H

#include <vector>

namespace tudat
{
namespace aerodynamics
{

double computeDragCoefficientOfCubeInRarefiedFlow( const std::vector< double >& input );

double computeDragCoefficientOfSphereInRarefiedFlow( const std::vector< double >& input );

} // namespace aerodynamics

} // namespace tudat

#endif // TUDAT_RAREFIEDFLOWSIMPLEGEOMETRYDRAG_H
