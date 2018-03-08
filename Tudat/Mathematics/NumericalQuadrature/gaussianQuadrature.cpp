/*    Copyright (c) 2010-2018, Delft University of Technology
 *
 *    All rigths reserved
 *
 *    This file is part of the Tudat. Redistribution and use in source and
 *    binary forms, with or without modification, are permitted exclusively
 *    under the terms of the Modified BSD license. You should have received
 *    a copy of the license with this file. If not, please or visit:
 *    http://tudat.tudelft.nl/LICENSE.
 */

#include "Tudat/Mathematics/NumericalQuadrature/gaussianQuadrature.h"

namespace tudat
{

namespace numerical_quadrature
{


//! Function to create Gauss quadrature node/weight container with long double precision.
/*!
 *  Function to create Gauss quadrature node/weight container with long double precision.
 *  \return Gauss quadrature node/weight container
 */
template< >
boost::shared_ptr< GaussQuadratureNodesAndWeights< long double > >
getGaussQuadratureNodesAndWeights( )
{
    return longDoubleGaussQuadratureNodesAndWeights;
}

//! Function to create Gauss quadrature node/weight container with double precision.
/*!
 *  Function to create Gauss quadrature node/weight container with double precision.
 *  \return Gauss quadrature node/weight container
 */
template< >
boost::shared_ptr< GaussQuadratureNodesAndWeights< double > >
getGaussQuadratureNodesAndWeights( )
{
    return doubleGaussQuadratureNodesAndWeights;
}

//! Function to create Gauss quadrature node/weight container with float precision.
/*!
 *  Function to create Gauss quadrature node/weight container with float precision.
 *  \return Gauss quadrature node/weight container
 */
template< >
boost::shared_ptr< GaussQuadratureNodesAndWeights< float > >
getGaussQuadratureNodesAndWeights( )
{
    return floatGaussQuadratureNodesAndWeights;
}

} // namespace numerical_quadrature

} // namespace tudat

