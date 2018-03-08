/*    Copyright (c) 2010-2017, Delft University of Technology
 *    All rigths reserved
 *
 *    This file is part of the Tudat. Redistribution and use in source and
 *    binary forms, with or without modification, are permitted exclusively
 *    under the terms of the Modified BSD license. You should have received
 *    a copy of the license with this file. If not, please or visit:
 *    http://tudat.tudelft.nl/LICENSE.
 */

#include "centralBodyPertubed.h"

namespace tudat
{

namespace propagators
{

namespace sst
{

namespace force_models
{


//! Update instance's members that are computed from the current auxiliary elements.
void CentralBodyPerturbed::updateMembers( )
{
    // Direction cosines
    alpha = aux.f( 2 );
    beta  = aux.g( 2 );
    gamma = aux.w( 2 );
}


} // namespace force_models

} // namespace sst

} // namespace propagators

} // namespace tudat
