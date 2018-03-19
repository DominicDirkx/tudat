
#include "Tudat/Astrodynamics/OrbitDetermination/orbitDeterminationManager.h"

namespace tudat
{

namespace simulation_setup
{

template class OrbitDeterminationManager< double, double >;
template class OrbitDeterminationManager< double, Time >;
template class OrbitDeterminationManager< long double, double >;
template class OrbitDeterminationManager< long double, Time >;

}

}
