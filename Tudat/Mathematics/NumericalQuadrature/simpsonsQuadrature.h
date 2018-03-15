/*    Copyright (c) 2010-2018, Delft University of Technology
 *    All rigths reserved
 *
 *    This file is part of the Tudat. Redistribution and use in source and
 *    binary forms, with or without modification, are permitted exclusively
 *    under the terms of the Modified BSD license. You should have received
 *    a copy of the license with this file. If not, please or visit:
 *    http://tudat.tudelft.nl/LICENSE.
 */

#ifndef TUDAT_SIMPSONS_QUADRATURE_H
#define TUDAT_SIMPSONS_QUADRATURE_H

#include <vector>

#include <boost/shared_ptr.hpp>
#include <boost/function.hpp>

#include <Eigen/Core>

#include "Tudat/Mathematics/NumericalQuadrature/numericalQuadrature.h"

namespace tudat
{

namespace numerical_quadrature
{


//! Function to perform numerical quadrature using the Simpson's method.
/*!
 * Function to perform numerical quadrature using the Simpson's method.
 * \param independentVariables Values of independent variables at which dependentVariables are given
 * \param dependentVariables Values of function for which the numerical quadrature is to be computed, given at
 * independentVariables.
 * \return Numerical quadrature (integral) of the data provided as input
 */
template< typename IndependentVariableType, typename DependentVariableType >
DependentVariableType performSimpsonQuadrature(
        const std::vector< IndependentVariableType >& independentVariables,
        const  boost::function< DependentVariableType( const IndependentVariableType ) > dependentVariableFunction )
{
    DependentVariableType integral;
    IndependentVariableType timeStep;
    for( unsigned int i = 0 ; i < independentVariables.size( ) - 1 ; i++ )
    {
        timeStep = independentVariables[ i + 1 ] - independentVariables[ i ];
        if( i == 0 )
        {
            integral = static_cast< DependentVariableType >( timeStep ) / 6.0 * (
                        dependentVariableFunction( independentVariables[ i ] ) +
                        4.0 * dependentVariableFunction( independentVariables[ i ] + timeStep / 2.0 ) +
                        dependentVariableFunction( independentVariables[ i + 1 ] ) );
        }
        else
        {
            integral += static_cast< DependentVariableType >( timeStep ) / 6.0 * (
                        dependentVariableFunction( independentVariables[ i ] ) +
                        4.0 * dependentVariableFunction( independentVariables[ i ] + timeStep / 2.0 ) +
                        dependentVariableFunction( independentVariables[ i + 1 ] ) );
        }

    }
    return integral;
}

//! Simpson's method numerical quadrature wrapper class.
/*!
 *  Numerical method that uses the Simpson's method to compute definite integrals of a dataset.
 */
template< typename IndependentVariableType, typename DependentVariableType >
class SimpsonNumericalQuadrature : public NumericalQuadrature< IndependentVariableType , DependentVariableType >
{
public:

    //! Constructor.
    /*!
     * Constructor
     * \param independentVariables Values of independent variables at which dependentVariables are given
     * \param dependentVariables Values of function for which the numerical quadrature is to be computed, given at
     * independentVariables.
     */
    SimpsonNumericalQuadrature(
            const std::vector< IndependentVariableType >& independentVariables,
                                const boost::function< DependentVariableType( const IndependentVariableType ) > dependentVariableFunction ):
        NumericalQuadrature< IndependentVariableType , DependentVariableType >(
            independentVariables, dependentVariableFunction )
    {
        performQuadrature( );
    }

    SimpsonNumericalQuadrature( ){ }

protected:

    //! Function that is called to perform the numerical quadrature
    /*!
     * Function that is called to perform the numerical quadrature. Sets the result in the quadratureResult_ local
     * variable.
     */
    void performQuadrature( )
    {
        this->quadratureResult_ = performSimpsonQuadrature( this->independentVariables_ , this->dependentVariableFunction_ );
    }

private:

};

//! Typede for trapezoidal quadrature with double (in)dependent variables.
typedef boost::shared_ptr< SimpsonNumericalQuadrature< double, double > > SimpsonNumericalIntegratorPointerd;

} // namespace numerical_quadrature

} // namespace tudat

#endif // TUDAT_SIMPSONS_QUADRATURE_H
