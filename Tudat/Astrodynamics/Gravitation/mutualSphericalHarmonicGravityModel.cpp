/*    Copyright (c) 2010-2018, Delft University of Technology
 *    All rigths reserved
 *
 *    This file is part of the Tudat. Redistribution and use in source and
 *    binary forms, with or without modification, are permitted exclusively
 *    under the terms of the Modified BSD license. You should have received
 *    a copy of the license with this file. If not, please or visit:
 *    http://tudat.tudelft.nl/LICENSE.
 */
#include "Tudat/Astrodynamics/Gravitation/mutualSphericalHarmonicGravityModel.h"

namespace tudat
{

namespace gravitation
{

//! Function to manually remove the C(0,0) term from cosine coefficients.
void setDegreeAndOrderCoefficientToZero(
        const std::function< void( Eigen::MatrixXd& ) > originalCosineCoefficientFunction,
        Eigen::MatrixXd& cosineCoefficients )
{
    originalCosineCoefficientFunction( cosineCoefficients );
    cosineCoefficients( 0, 0 ) = 0.0;
}



}

}
