#ifndef SP3CFILEREADER_H
#define SP3CFILEREADER_H

#include <string>
#include <map>
#include <vector>

#include <boost/shared_ptr.hpp>

#include <Eigen/Core>

#include "tudat/astro/basic_astro/timeConversions.h"

namespace tudat
{

namespace input_output
{

struct SP3cFileContents
{
    std::string frameName;
    std::string timeScale;
    std::map< std::string, std::map< double, Eigen::VectorXd > > vehiclesStates;
};

boost::shared_ptr< SP3cFileContents > readSp3cFile(
        const std::string fileName, const double referenceJulianDay = basic_astrodynamics::JULIAN_DAY_ON_J2000 );

}

}

#endif // SP3CFILEREADER_H
