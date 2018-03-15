/*    Copyright (c) 2010-2018, Delft University of Technology
 *    All rigths reserved
 *
 *    This file is part of the Tudat. Redistribution and use in source and
 *    binary forms, with or without modification, are permitted exclusively
 *    under the terms of the Modified BSD license. You should have received
 *    a copy of the license with this file. If not, please or visit:
 *    http://tudat.tudelft.nl/LICENSE.
 */

#ifndef TUDAT_CREATE_NUMERICAL_QUADRATURE_H
#define TUDAT_CREATE_NUMERICAL_QUADRATURE_H

#include "Tudat/Mathematics/NumericalQuadrature/numericalQuadrature.h"
#include "Tudat/Mathematics/NumericalQuadrature/trapezoidQuadrature.h"
#include "Tudat/Mathematics/NumericalQuadrature/simpsonsQuadrature.h"
#include "Tudat/Mathematics/NumericalQuadrature/gaussianQuadrature.h"

namespace tudat
{

namespace numerical_quadrature
{

enum NumericalQuadratureTypes
{
    trapezoidal_quadrature,
    simpsons_quadrature,
    gaussian_quadrature
};

class NumericalQuadratureSettings
{
public:

    NumericalQuadratureSettings( const NumericalQuadratureTypes numericalQuadratureType ):
        numericalQuadratureType_( numericalQuadratureType ){ }

    virtual ~NumericalQuadratureSettings( ){ }

    NumericalQuadratureTypes numericalQuadratureType_;

};

class GaussianQuadratureSettings: public NumericalQuadratureSettings
{
public:

    GaussianQuadratureSettings( const int numberOfNodes ):
        NumericalQuadratureSettings( gaussian_quadrature ),
        numberOfNodes_( numberOfNodes ){ }

    ~GaussianQuadratureSettings( ){ }

    int numberOfNodes_;

};

template< typename IndependentVariableType, typename DependentVariableType >
boost::shared_ptr< NumericalQuadrature< IndependentVariableType, DependentVariableType > > createNumericalQuadrature(
        const std::vector< IndependentVariableType >& independentVariables,
        const boost::function< DependentVariableType( const IndependentVariableType ) > integrand,
        const boost::shared_ptr< NumericalQuadratureSettings > quadratureSettings )
{
    boost::shared_ptr< NumericalQuadrature< IndependentVariableType, DependentVariableType > > numericalQuadrature;

    switch( quadratureSettings->numericalQuadratureType_ )
    {
    case trapezoidal_quadrature:
    {
        numericalQuadrature = boost::make_shared< TrapezoidNumericalQuadrature< IndependentVariableType, DependentVariableType > >(
                   independentVariables, integrand );
        break;
    }
    case simpsons_quadrature:
    {
        numericalQuadrature = boost::make_shared< SimpsonNumericalQuadrature< IndependentVariableType, DependentVariableType > >(
                   independentVariables, integrand );
        break;
    }
    case gaussian_quadrature:
    {
        if( independentVariables.size( ) != 2 )
        {
            throw std::runtime_error( "Error, Gaussian quadrature supports single step only" );
        }
        boost::shared_ptr< GaussianQuadratureSettings > gaussQuadratureSettings =
                boost::dynamic_pointer_cast< GaussianQuadratureSettings >( quadratureSettings );

        if( sizeof( IndependentVariableType ) == 8 )
        {
            numericalQuadrature = boost::make_shared< GaussianQuadrature< IndependentVariableType, DependentVariableType, double > >(
                        integrand, independentVariables[ 0 ], independentVariables[ 1 ],
                    gaussQuadratureSettings->numberOfNodes_ );
        }
        else
        {
            numericalQuadrature = boost::make_shared< GaussianQuadrature< IndependentVariableType, DependentVariableType, long double > >(
                        integrand, independentVariables[ 0 ], independentVariables[ 1 ],
                    gaussQuadratureSettings->numberOfNodes_ );
        }
        break;
    }
    default:
        throw std::runtime_error( "Error when making numerical quadrature, did not recognize quadrature type" );
    }

    return numericalQuadrature;
}


} // namespace numerical_quadrature

} // namespace tudat

#endif // TUDAT_CREATE_NUMERICAL_QUADRATURE_H
