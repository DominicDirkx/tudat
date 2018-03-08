/*    Copyright (c) 2010-2017, Delft University of Technology
 *    All rigths reserved
 *
 *    This file is part of the Tudat. Redistribution and use in source and
 *    binary forms, with or without modification, are permitted exclusively
 *    under the terms of the Modified BSD license. You should have received
 *    a copy of the license with this file. If not, please or visit:
 *    http://tudat.tudelft.nl/LICENSE.
 */

#include "thirdBodyPertubed.h"

namespace tudat
{

namespace propagators
{

namespace sst
{

namespace force_models
{


//! Update instance's members that are computed from the current auxiliary elements.
void ThirdBodyPerturbed::updateMembers( )
{
    // Determine distance to third body
    R3 = r3.norm( );

    // Direction cosines
    const Eigen::Vector3d bodyDir = r3 / R3;
    alpha = bodyDir.dot( aux.f );
    beta  = bodyDir.dot( aux.g );
    gamma = bodyDir.dot( aux.w );
}


} // namespace force_models

} // namespace sst

} // namespace propagators

} // namespace tudat
