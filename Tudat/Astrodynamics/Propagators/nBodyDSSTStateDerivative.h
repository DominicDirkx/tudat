/*    Copyright (c) 2010-2017, Delft University of Technology
 *    All rigths reserved
 *
 *    This file is part of the Tudat. Redistribution and use in source and
 *    binary forms, with or without modification, are permitted exclusively
 *    under the terms of the Modified BSD license. You should have received
 *    a copy of the license with this file. If not, please or visit:
 *    http://tudat.tudelft.nl/LICENSE.
 */

#ifndef TUDAT_NBODYDSSTSTATEDERIVATIVE_H
#define TUDAT_NBODYDSSTSTATEDERIVATIVE_H

#include <chrono>

#include "Tudat/Astrodynamics/Propagators/nBodyStateDerivative.h"

#include "Tudat/SimulationSetup/EnvironmentSetup/body.h"
#include "Tudat/SimulationSetup/PropagationSetup/propagationSettings.h"
#include "Tudat/SimulationSetup/PropagationSetup/createEnvironmentUpdater.h"

#include "Tudat/Astrodynamics/BasicAstrodynamics/orbitalElementConversions.h"
#include "Tudat/Astrodynamics/BasicAstrodynamics/accelerationModelTypes.h"

#include "Tudat/Astrodynamics/Propagators/DSST/forces/availableForceModels.h"
#include "Tudat/Astrodynamics/Propagators/DSST/utilities/osculatingToMean.h"

namespace tudat
{

namespace propagators
{

//! Class for computing the state derivative of translational motion of N bodies, using a DSST propagator.
template< typename StateScalarType = double, typename TimeType = double >
class NBodyDSSTStateDerivative: public NBodyStateDerivative< StateScalarType, TimeType >
{
public:

    //! Constructor
    /**
     * Constructor.
     * @param propagatorSettings
     * @param centralBodyData
     * @param bodyMap
     * @param initialEpoch
     */
    NBodyDSSTStateDerivative(
            const boost::shared_ptr< TranslationalStatePropagatorSettings< StateScalarType > > propagatorSettings,
            const boost::shared_ptr< CentralBodyData< StateScalarType, TimeType > > centralBodyData,
            const simulation_setup::NamedBodyMap &bodyMap, const TimeType &initialEpoch ) :
        NBodyStateDerivative< StateScalarType, TimeType >( propagatorSettings->getAccelerationsMap( ), centralBodyData,
                                                           dsst, propagatorSettings->bodiesToIntegrate_ )
    {
        using namespace basic_astrodynamics;
        using namespace orbital_element_conversions;

        using namespace aerodynamics;
        using namespace electro_magnetism;
        using namespace gravitation;

        using namespace sst;
        using namespace sst::force_models;

        // Check only one propagated body
        if ( propagatorSettings->bodiesToIntegrate_.size() > 1 )
        {
            throw std::runtime_error("DSST propagator does not support simultaneous propagation of multiple bodies.");
        }

        // Check if user provided settings for Force Models (by using a DSSTTranslationalStatePropagatorSettings)
        boost::shared_ptr< DSSTTranslationalStatePropagatorSettings< StateScalarType > > dsstPropagatorSettings =
                boost::dynamic_pointer_cast< DSSTTranslationalStatePropagatorSettings< StateScalarType > >(
                    propagatorSettings );
        const bool dsstPropagatorSettingsProvided = dsstPropagatorSettings != NULL;

        // Get initial Cartesian elements
        const Eigen::Vector6d initialCartesianElements =
                propagatorSettings->getInitialStates().segment( 0, 6 ).template cast< double >();

        // Get central body
        const std::string centralBodyName = *centralBodyData->getCentralBodies().begin();
        BodyPtr centralBody = bodyMap.at( centralBodyName );

        // Get acceleration models for propagated body
        const std::string propagatedBodyName = *propagatorSettings->bodiesToIntegrate_.begin();
        BodyPtr propagatedBody = bodyMap.at( propagatedBodyName );
        const SingleBodyAccelerationMap accelerationModelsMap =
                propagatorSettings->getAccelerationsMap( ).at( propagatedBodyName );

        // First iteration through acceleration models, only to find central gravity model
        CentralGravityAM centralGravityAM;
        for ( auto ent: accelerationModelsMap )
        {
            const std::string bodyExertingAccelerations = ent.first;
            const std::vector< AccelerationModel3dPointer > accelerationModels = ent.second;
            for ( auto accelerationModel: accelerationModels )
            {
                centralGravityAM = boost::dynamic_pointer_cast<
                        SphericalHarmonicsGravitationalAccelerationModelBase< Eigen::Vector3d > >( accelerationModel );
                if ( centralGravityAM != NULL ) { break; }
            }
            if ( centralGravityAM != NULL ) { break; }
        }

        if ( centralGravityAM == NULL ) {
            throw std::runtime_error( "Could not find a central or spherical harmonics gravity model required "
                                      "for setting up the DSST propagator." );
        }

        // Convert from Cartesian to (eccentric) equinoctial elements
        EquinoctialElements equinoctialElements = EquinoctialElements::fromCartesian(
                    initialCartesianElements, centralGravityAM->getGravitationalParameterFunction()(), eccentricType );

        // Create auxiliary elements
        auxiliaryElements =
                AuxiliaryElements( initialEpoch, equinoctialElements, propagatedBody, centralBody, centralGravityAM );

        // Iterate through all acceleration models to create the DSST force models
        for ( auto ent: accelerationModelsMap )
        {
            const std::string perturbingBody = ent.first;
            const std::vector< AccelerationModel3dPointer > accelerationModels = ent.second;
            for ( AccelerationModel3dPointer accelerationModel: accelerationModels )
            {
                AvailableAcceleration accelerationType = getAccelerationModelType( accelerationModel );
                ForceIdentifier forceID( perturbingBody, accelerationType );
                boost::shared_ptr< ForceModelSettings > settings = dsstPropagatorSettingsProvided ?
                            dsstPropagatorSettings->getForceModelSettings( forceID ) : NULL;
                switch ( accelerationType ) {
                case central_gravity:
                {
                    // No DSST force model required
                    break;
                }
                case spherical_harmonic_gravity:
                {
                    SphericalHarmonicsAM sphericalHarmonicsAM = boost::dynamic_pointer_cast<
                            SphericalHarmonicsGravitationalAccelerationModel >( accelerationModel );
                    forceModels[ forceID ] = boost::make_shared< ZonalSphericalHarmonicGravity >(
                                               auxiliaryElements, sphericalHarmonicsAM );
                    break;
                }
                case third_body_central_gravity:
                {
                    ThirdBodyAM thirdBodyAM = boost::dynamic_pointer_cast<
                            ThirdBodyCentralGravityAcceleration >( accelerationModel );
                    forceModels[ forceID ] = boost::make_shared< ThirdBodyCentralGravity >(
                                auxiliaryElements, thirdBodyAM,
                                boost::dynamic_pointer_cast< ConservativeSettings >( settings ) );
                    break;
                }
                case aerodynamic:
                {
                    AerodynamicAM aerodynamicAM = boost::dynamic_pointer_cast<
                            AerodynamicAcceleration >( accelerationModel );
                    forceModels[ forceID ] = boost::make_shared< AtmosphericDrag >(
                                auxiliaryElements, perturbingBody, aerodynamicAM,
                                boost::dynamic_pointer_cast< AtmosphericDragSettings >( settings ) );
                    break;
                }
                case cannon_ball_radiation_pressure:
                {
                    RadiationPressureAM radiationPressureAM = boost::dynamic_pointer_cast<
                            CannonBallRadiationPressureAcceleration >( accelerationModel );
                    const unsigned int numberOfOcultingBodies =
                            radiationPressureAM->getRadiationPressureInterface()->getOccultingBodyRadii().size();
                    if ( numberOfOcultingBodies == 0 ) {
                        forceModels[ forceID ] = boost::make_shared< ConservativeRadiationPressure >(
                                    auxiliaryElements, radiationPressureAM,
                                    boost::dynamic_pointer_cast< ConservativeSettings >( settings ) );
                    } else if ( numberOfOcultingBodies == 1 ) {
                        forceModels[ forceID ] = boost::make_shared< RadiationPressure >(
                                    auxiliaryElements, perturbingBody, radiationPressureAM,
                                    boost::dynamic_pointer_cast< RadiationPressureSettings >( settings ) );
                    } else {
                        throw std::runtime_error( "Multiple occultations is not supported by the DSST propagator." );
                    }
                    break;
                }
                default:
                {
                    throw std::runtime_error( "Error, DSST propagator does not support accelerations of type: " +
                                              boost::lexical_cast< std::string >( accelerationType ) );
                    break;
                }
                }
            }
        }

        // Create environment updaters for non-conservative force models
        for ( auto ent: forceModels )
        {
            boost::shared_ptr< NonConservative > nonConservativeForceModel =
                    boost::dynamic_pointer_cast< NonConservative >( ent.second );
            if ( nonConservativeForceModel != NULL )
            {
                // Create partial acceleration map only with relevant acceleration model
                AccelerationMap partialAccelerationMap =
                { { propagatedBodyName, { { nonConservativeForceModel->getPerturbingBodyName(),
                                        { nonConservativeForceModel->getAccelerationModel() } } } } };

                // Create environment updater settings
                auto updateSettings = createTranslationalEquationsOfMotionEnvironmentUpdaterSettings(
                            partialAccelerationMap, bodyMap );
                std::map< IntegratedStateType, std::vector< std::pair< std::string, std::string > > >
                        integratedStates = { { transational_state, {
                                                   std::pair< std::string, std::string >( propagatedBodyName,"" ) } } };
                EnvironmentUpdater< double, double > updater( bodyMap, updateSettings, integratedStates );

                // Set updater
                nonConservativeForceModel->setEnvironmentUpdater(
                            boost::make_shared< EnvironmentUpdater< double, double > >( updater ) );
            }
        }

        // FIXME
        // Convert initial state from osculating to mean elements
        element_conversions::transformOsculatingToMeanElements( auxiliaryElements, forceModels );

        // std::cout << propagatorSettings->getInitialStates().transpose() << std::endl;

        const Eigen::Vector6d meanInitialCartesianElements =
                auxiliaryElements.equinoctialElements.toCartesian( auxiliaryElements.mu );
        const Eigen::Matrix< StateScalarType, Eigen::Dynamic, 1 > templatedMeanInitialCartesianElements =
                meanInitialCartesianElements.template cast< StateScalarType >();
        propagatorSettings->resetInitialStates( templatedMeanInitialCartesianElements );

        // std::cout << propagatorSettings->getInitialStates().transpose() << std::endl;

        // Store initial (total) short-period terms
        shortPeriodTermsMap[ ForceIdentifier() ][ initialEpoch ] =
                auxiliaryElements.equinoctialElements.getComponents( meanType ) -
                equinoctialElements.getComponents( meanType );
    }


    //! Destructor
    ~NBodyDSSTStateDerivative( ){ }

    //! Calculates the state derivative of the translational motion of the system.
    /*!
     * Calculates the state derivative (velocity+acceleration of each body) of the translational motion of the system
     * at the given time and position/velocity of bodies.
     *  \param time Time (TDB seconds since J2000) at which the system is to be updated.
     *  \param stateOfSystemToBeIntegrated List of 6 * bodiesToBeIntegratedNumerically_.size( ), containing Caartesian
     *  position/velocity of the bodies being integrated. The order of the values is defined by the order of bodies in
     *  bodiesToBeIntegratedNumerically_
     *  \param stateDerivative Current state derivative (velocity+acceleration) of system of bodies integrated numerically
     *  (returned by reference).
     */
    void calculateSystemStateDerivative(
            const TimeType time, const Eigen::Matrix< StateScalarType, Eigen::Dynamic, 1 >&stateOfSystemToBeIntegrated,
            Eigen::Block< Eigen::Matrix< StateScalarType, Eigen::Dynamic, Eigen::Dynamic > > stateDerivative )
    {
        using namespace Eigen;
        using namespace orbital_element_conversions;
        using namespace sst;
        using namespace sst::force_models;

        using namespace std::chrono;

        Vector6d equinoctialComponents = stateOfSystemToBeIntegrated.template cast< double >();
        auxiliaryElements.updateState( time, equinoctialComponents, meanType );

        // std::cout << time << ": " << auxiliaryElements.equinoctialElements.getComponents().transpose() << std::endl;

        Vector6d meanElementRates = Vector6d::Zero();
        Vector6d shortPeriodTerms = Vector6d::Zero();
        double totalComputationTime = 0.0;

        for ( auto ent: forceModels )
        {
            const ForceIdentifier forceID = ent.first;
            boost::shared_ptr< ForceModel > forceModel = ent.second;
            auto t = steady_clock::now();

            // Compute mean element rates and short period terms for the current force model and epoch
            Vector6d forceModelMeanElementRates = forceModel->getMeanElementRates();
            // FIXME: Vector6d forceModelShortPeriodTerms = forceModel->getShortPeriodTerms();
            Vector6d forceModelShortPeriodTerms = Vector6d::Zero();

            // Determine and store computation time
            double computationTime = duration_cast< microseconds >( steady_clock::now() - t ).count() / 1e6;  // s
            computationTimes[ forceID ][ time ] = computationTime;
            totalComputationTime += computationTime;

            // Store mean element rates and short period terms for the current force model and epoch
            meanElementRatesMap[ forceID ][ time ] = forceModelMeanElementRates;
            shortPeriodTermsMap[ forceID ][ time ] = forceModelShortPeriodTerms;

            // Add mean element rates and short period terms of the current force model to the total
            meanElementRates += forceModelMeanElementRates;
            shortPeriodTerms += forceModelShortPeriodTerms;

            // std::cout << forceID << ": " << forceModelMeanElementRates.transpose() << std::endl;
            // std::cout << forceID << ": " << forceModelShortPeriodTerms.transpose() << std::endl;
        }
/*
        for ( auto ent: computationTimes )
        {
            if ( ent.first != "ALL" ) {
                std::cout << ent.first << ": " << ent.second << " (" << ent.second / computationTimes[ "ALL" ] * 100 << "%), ";
            } else {
                std::cout << std::endl;
            }
        }
*/
        /*
        Vector6d normalizedMeanElementRates = meanElementRates;
        normalizedMeanElementRates( semiMajorAxisIndex ) /= auxiliaryElements.a;
        normalizedMeanElementRates( fastVariableIndex  ) /= 2 * mathematical_constants::PI;
        std::cout << 1 / normalizedMeanElementRates.norm() / physical_constants::JULIAN_DAY << " days" << std::endl;
        */

        // Add mean motion to the mean longitude rate
        meanElementRates( fastVariableIndex ) += auxiliaryElements.meanMotion;
        // meanElementRates( fastVariableIndex ) = 0;

        // Store mean element rates and short period terms for later reconstruction of osculating terms and/or output
        meanElementRatesMap[ ForceIdentifier() ][ time ] = meanElementRates;
        shortPeriodTermsMap[ ForceIdentifier() ][ time ] = shortPeriodTerms;
        computationTimes   [ ForceIdentifier() ][ time ] = totalComputationTime;

        // std::cout << "Total mean element rates: " << meanElementRates.transpose() << std::endl;
        /*
        // Rate of change of the perigee altitude
        const double a = auxiliaryElements.a;
        const double e = auxiliaryElements.e;
        const double h = auxiliaryElements.h;
        const double k = auxiliaryElements.k;
        const double da = meanElementRates( semiMajorAxisIndex );
        const double dh = meanElementRates( hIndex );
        const double dk = meanElementRates( kIndex );
        const double dhp = da * ( 1 - e ) - a / e * ( h * dh + k * dk );
        // std::cout << dhp / 1e3 * physical_constants::JULIAN_DAY << " km/day" << std::endl;
        */
        // std::cout << "hp = " << a * ( 1 - e ) / 1e3 << " km" << std::endl;

        // std::cout << "\n" << std::endl;

        // Update state derivative (only one body being propagated, thus only first column needs to be updated)
        stateDerivative.col( 0 ) = meanElementRates.template cast< StateScalarType >();
    }

    //! Function to convert the state in the conventional form to the DSST propagator-specific form.
    /*!
     * Function to convert the state in the conventional form to the propagator-specific form. For the Cowell propagator,
     * the two are equivalent, and this function returns the input state.
     * \param cartesianSolution State in 'conventional form'
     * \param time Current time at which the state is valid (not used in this class).
     * \return State (outputSolution), converted to the 'propagator-specific form' (which is equal to outputSolution).
     */
    Eigen::Matrix< StateScalarType, Eigen::Dynamic, Eigen::Dynamic >
    convertFromOutputSolution(
            const Eigen::Matrix< StateScalarType, Eigen::Dynamic, Eigen::Dynamic >& cartesianSolution,
            const TimeType& time )
    {
        using namespace Eigen;
        using namespace orbital_element_conversions;
        using namespace sst;

        Matrix< StateScalarType, Dynamic, Dynamic > dsstSolution;
        dsstSolution.resizeLike( cartesianSolution );

        // std::cout << "cartesianSolution: " << cartesianSolution << std::endl;

        for ( unsigned int i = 0; i < cartesianSolution.cols(); i++ )
        {
            Vector6d cartesianComponents = cartesianSolution.col( i ).template cast< double >();

            EquinoctialElements equinoctialElements =
                    EquinoctialElements::fromCartesian( cartesianComponents, auxiliaryElements.mu, meanType );

            dsstSolution.col( i ) = equinoctialElements.getComponents().template cast< StateScalarType >();
        }

        // std::cout << "Solution: " << dsstSolution.transpose() << "\n" << std::endl;

        return dsstSolution;
    }

    //! Function to convert the DSST propagator-specific form of the state to the conventional form.
    /*!
     * Function to convert the propagator-specific form of the state to the conventional form. For the Cowell propagator,
     * the two are equivalent, and this function returns the input state.
     * In contrast to the convertCurrentStateToGlobalRepresentation function, this
     * function does not provide the state in the inertial frame, but instead provides it in the
     * frame in which it is propagated.
     * \param internalSolution State in propagator-specific form (i.e. form that is used in
     * numerical integration, equal to conventional form for this class).
     * \param time Current time at which the state is valid (not used in this class).
     * \param currentCartesianLocalSoluton State (internalSolution), converted to the 'conventional form',
     *  which is equal to outputSolution for this class (returned by reference).
     */
    void convertToOutputSolution(
            const Eigen::Matrix< StateScalarType, Eigen::Dynamic, Eigen::Dynamic >& internalSolution,
            const TimeType& time,
            Eigen::Block< Eigen::Matrix< StateScalarType, Eigen::Dynamic, 1 > > currentCartesianLocalSoluton )
    {
        using namespace Eigen;
        using namespace orbital_element_conversions;
        using namespace sst;

        // std::cout << "internalSolution: " << internalSolution << std::endl;

        // Get mean elements
        Vector6d equinoctialComponents = internalSolution.template cast< double >();

        // FIXME: this is not only used to generate the output, but also during computation of k1, k2... of RK4. Why?
        // But the RK4 needs only the mean elements (in the propagator-specific representation), not the osculating.
        // Add short period terms -> osculating elements
        // equinoctialComponents += shortPeriodTermsMap.at( time );

        EquinoctialElements equinoctialElements( equinoctialComponents, meanType,
                                                 auxiliaryElements.equinoctialElements.isRetrograde() );

        Vector6d cartesian = equinoctialElements.toCartesian( auxiliaryElements.mu );

        currentCartesianLocalSoluton = cartesian.template cast< StateScalarType >();
        // std::cout << "currentCartesianLocalSoluton: \n" << cartesian.head(3).norm() << "\n\n" << std::endl;

    }


    //! Get mean element rates for the current (i.e. last) epoch for a specific ForceIdentifier
    Eigen::Vector6d getCurrentMeanElementRates( const sst::force_models::ForceIdentifier& forceID ) {
        return ( --meanElementRatesMap.at( forceID ).end( ) )->second;
    }

    //! Get short period terms for the current (i.e. last) epoch for a specific ForceIdentifier
    Eigen::Vector6d getCurrentShortPeriodTerms( const sst::force_models::ForceIdentifier& forceID ) {
        return ( --shortPeriodTermsMap.at( forceID ).end( ) )->second;
    }

    //! Get computation time for the current (i.e. last) epoch for a specific ForceIdentifier
    double getCurrentComputationTimes( const sst::force_models::ForceIdentifier& forceID ) {
        return ( --computationTimes.at( forceID ).end( ) )->second;
    }


private:

    //! Auxiliary elements shared amongst all the force models
    sst::AuxiliaryElements auxiliaryElements;

    //! Map of DSST force models
    //! Keys are the name of the forces, e.g. Earth-DRAG, Moon-3RD, etc.
    //! Values are pointers to the force models.
    sst::force_models::ForceModelMap forceModels;


    //! Computed mean element rates for each force model at each integrated epoch
    sst::force_models::HistoryForceMap< Eigen::Vector6d > meanElementRatesMap;

    //! Computed short perido terms for each force model at each integrated epoch
    sst::force_models::HistoryForceMap< Eigen::Vector6d > shortPeriodTermsMap;

    //! Cummulative time spent on computing the effects of each force model at each integrated epoch
    sst::force_models::HistoryForceMap< double > computationTimes;

};

} // namespace propagators

} // namespace tudat

#endif // TUDAT_NBODYDSSTSTATEDERIVATIVE_H
