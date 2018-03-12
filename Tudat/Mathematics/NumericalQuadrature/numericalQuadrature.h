/*    Copyright (c) 2010-2018, Delft University of Technology
 *    All rigths reserved
 *
 *    This file is part of the Tudat. Redistribution and use in source and
 *    binary forms, with or without modification, are permitted exclusively
 *    under the terms of the Modified BSD license. You should have received
 *    a copy of the license with this file. If not, please or visit:
 *    http://tudat.tudelft.nl/LICENSE.
 */

#ifndef TUDAT_NUMERICAL_QUADRATURE_H
#define TUDAT_NUMERICAL_QUADRATURE_H

#include <vector>
#include <boost/function.hpp>

#include <Eigen/Core>

namespace tudat
{
namespace numerical_quadrature
{

//! Base class for numerical quadrature.
/*!
 * Base class for the numerical quadrature methods included in Tudat, the dependent and independent variable
 * types are specified as template parameters.
 */
template< typename IndependentVariableType, typename DependentVariableType >
class NumericalQuadrature
{
public:

    NumericalQuadrature( const std::vector< IndependentVariableType >& independentVariables,
                         const boost::function< DependentVariableType( const IndependentVariableType ) > dependentVariableFunction )
    {
        independentVariables_ = independentVariables;
        dependentVariableFunction_ = dependentVariableFunction;
    }

    NumericalQuadrature( ){ }

    //! Destructor.
    virtual ~NumericalQuadrature( ) { }


    //! Function to reset the (in)dependent variable values.
    /*!
     * Function to reset the (in)dependent variable values.
     * \param independentVariables Values of independent variables at which dependentVariables are given
     * \param dependentVariables Values of function for which the numerical quadrature is to be computed, given at
     * independentVariables.
     */
    virtual void resetData( const std::vector< IndependentVariableType >& independentVariables,
                    const boost::function< DependentVariableType( const IndependentVariableType ) > dependentVariableFunction )
    {
        independentVariables_ = independentVariables;
        dependentVariableFunction_ = dependentVariableFunction;
        performQuadrature( );
    }

    //! Function to return computed value of the quadrature.
    /*!
     *  Function to return computed value of the quadrature, as computed by last call to performQuadrature.
     *  \return Function to return computed value of the quadrature, as computed by last call to performQuadrature.
     */
    virtual DependentVariableType getQuadrature( )
    {
        return quadratureResult_;
    }

protected:

    //! Function that is called to perform the numerical quadrature
    /*!
     * Function that is called to perform the numerical quadrature. It must be implemented in derived classes. Any
     * implementation of performQuadrature must be made consistent with that of getQuadrature.
     */
    virtual void performQuadrature( ) = 0;

    //! Independent variables.
    std::vector< IndependentVariableType > independentVariables_;

    //! Dependent variables.
    boost::function< DependentVariableType( const IndependentVariableType ) > dependentVariableFunction_;

    //! Computed value of the quadrature, as computed by last call to performQuadrature.
    DependentVariableType quadratureResult_;

};

} // namespace numerical_quadrature

} // namespace tudat

#endif // TUDAT_NUMERICAL_QUADRATURE_H
