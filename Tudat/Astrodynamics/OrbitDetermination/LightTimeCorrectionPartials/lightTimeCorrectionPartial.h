#ifndef TUDAT_LIGHTTIMECORRECTIONPARTIAL_H
#define TUDAT_LIGHTTIMECORRECTIONPARTIAL_H

#include <boost/smart_ptr/enable_shared_from_this.hpp>

#include <vector>

#include <Eigen/Core>

#include "Tudat/Mathematics/BasicMathematics/linearAlgebraTypes.h"

#include "Tudat/Astrodynamics/ObservationModels/linkTypeDefs.h"
#include "Tudat/Astrodynamics/ObservationModels/lightTimeSolution.h"
#include "Tudat/Astrodynamics/ObservationModels/ObservableCorrections/lightTimeCorrection.h"
#include "Tudat/Astrodynamics/OrbitDetermination/EstimatableParameters/estimatableParameter.h"

namespace tudat
{

namespace observation_partials
{

//! Base class for computing the partial derivatives of light-time corrections w.r.t. estimated parameters.
/*!
 *  Base class for computing the partial derivatives of light-time corrections w.r.t. estimated parameters. A single derived
 *  class is implemented per light time correction type. The same derived class is used for computing partials of single
 *  correction type w.r.t. any number of physical parameters
 */
class LightTimeCorrectionPartial
{
public:

    //! Typedef for return argument for light time partial (partial with vector and associated time)
    typedef std::pair< Eigen::Matrix< double, 1, Eigen::Dynamic >, double > SingleOneWayRangePartialReturnType;

    //! Constructor
    /*!
     * Constructor
     * \param correctionType Type of light-time correction for which partial is to be computed
     */
    LightTimeCorrectionPartial( const observation_models::LightTimeCorrectionType correctionType ):
        correctionType_( correctionType ){ }

    //! Destructor
    virtual ~LightTimeCorrectionPartial( ){ }

    //! Function to return type of light-time correction for which partial is to be computed
    /*!
     * Function to return type of light-time correction for which partial is to be computed
     * \return Type of light-time correction for which partial is to be computed
     */
    observation_models::LightTimeCorrectionType getCorrectionType( )
    {
        return correctionType_;
    }

protected:

    //! Type of light-time correction for which partial is to be computed
    observation_models::LightTimeCorrectionType correctionType_;

};

//! Function to get the function returning the light-time correction partial for given correction partial and parameter.
/*!
 * Function to get the function returning the light-time correction partial for given correction partial and parameter.
 * \param parameterId Parameter for which partial derivative is to be computed.
 * \param lightTimeCorrectionPartial Partial object from which partial function is to be retrieved.
 * \return Pair (with second entry boolean which is true if the partial is non-zero). First entry of pair gives
 * partial and associated time as a function of link-end states and times.
 */
std::pair< boost::function< LightTimeCorrectionPartial::SingleOneWayRangePartialReturnType(
        const std::vector< basic_mathematics::Vector6d >&, const std::vector< double >& ) >, bool >
getLightTimeParameterPartialFunction(
        const estimatable_parameters::EstimatebleParameterIdentifier parameterId,
        const boost::shared_ptr< LightTimeCorrectionPartial > lightTimeCorrectionPartial );

}

}

#endif // TUDAT_LIGHTTIMECORRECTIONPARTIAL_H