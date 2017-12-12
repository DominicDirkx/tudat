#include "Tudat/Mathematics/BasicMathematics/cayleyKleinParameters.h"

namespace tudat
{
namespace basic_mathematics
{

//! Function to convert a quaternion representing a rotation to two Cayley-Klein parameters
void convertQuaterionToCayleyKleinParameters(
        const Eigen::Quaterniond quaternion, std::complex< double >& cayleyKleinA, std::complex< double >& cayleyKleinB )
{
    cayleyKleinA = std::complex< double >( quaternion.w( ), -quaternion.z( ) );
    cayleyKleinB = std::complex< double >( quaternion.y( ), -quaternion.x( ) );
}

//! Function to convert 3-2-3 Euler angles representing a rotation to two Cayley-Klein parameters
void convert323EulerAnglesToCayleyKleinParameters(
        const double firstZRotation, const double yRotation, const double secondZRotation,
        std::complex< double >& cayleyKleinA, std::complex< double >& cayleyKleinB )
{
    cayleyKleinA = std::exp( mathematical_constants::COMPLEX_I * ( firstZRotation + secondZRotation ) / 2.0 ) *
            std::cos( yRotation/ 2.0 );

    cayleyKleinB = -std::exp( mathematical_constants::COMPLEX_I * ( firstZRotation - secondZRotation ) / 2.0 ) *
            std::sin( yRotation/ 2.0 );
}

//! Function to convert 3-1-3 Euler angles representing a rotation to two Cayley-Klein parameters
void convert313EulerAnglesToCayleyKleinParameters(
        const double firstZRotation, const double xRotation, const double secondZRotation,
        std::complex< double >& cayleyKleinA, std::complex< double >& cayleyKleinB )
{
    cayleyKleinA = std::exp( mathematical_constants::COMPLEX_I * ( firstZRotation + secondZRotation ) / 2.0 ) *
            std::cos( xRotation / 2.0 );
    cayleyKleinB = mathematical_constants::COMPLEX_I * std::exp( mathematical_constants::COMPLEX_I *
                                                                 ( firstZRotation - secondZRotation ) / 2.0 ) *
            std::sin( xRotation / 2.0 );
}


} // namespace basic_mathematics

} // namespace tudat
