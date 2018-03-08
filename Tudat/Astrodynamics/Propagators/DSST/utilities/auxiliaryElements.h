#ifndef TUDAT_PROPAGATORS_DSST_AUXILIARYELEMENTS_H
#define TUDAT_PROPAGATORS_DSST_AUXILIARYELEMENTS_H

#include <boost/make_shared.hpp>

#include "Tudat/SimulationSetup/EnvironmentSetup/body.h"
#include "Tudat/Astrodynamics/Gravitation/centralGravityModel.h"
#include "Tudat/SimulationSetup/PropagationSetup/environmentUpdater.h"
#include "Tudat/Mathematics/NumericalQuadrature/gaussianQuadrature.h"

// #include "Tudat/Mathematics/BasicMathematics/basicMathematicsFunctions.h"

#include "Tudat/Astrodynamics/Propagators/DSST/utilities/vectors.h"
#include "Tudat/Astrodynamics/Propagators/DSST/utilities/coefficientsFactories.h"


namespace tudat
{

namespace propagators
{

namespace sst
{


enum ElementSet
{
    automatic  = 0,  // direct for inclination ≤ 135 deg, retrograde for inclination > 135 deg
    direct     = 1,  // suitable for any inclination not close to 180 deg
    retrograde = 2   // suitable for any inclination not close to 0 deg
};


//! Class containing equinoctial elements and functionality to efficiently transform between mean, true and eccentric
class EquinoctialElements
{
public:

    //! Empty constructor
    EquinoctialElements( ) { }

    //! Basic constructor
    EquinoctialElements( const Eigen::Vector6d &components, const orbital_element_conversions::FastVariableType type,
                         const bool retrograde )
    {
        setComponents( components, type, retrograde );
    }

    //! Construct from Keplerian elements
    static EquinoctialElements fromKeplerian( const Eigen::Vector6d &components,
                                              const orbital_element_conversions::FastVariableType type,
                                              const ElementSet set = automatic ) {
        using namespace orbital_element_conversions;
        using namespace mathematical_constants;
        bool retro = set == retrograde;
        if ( set == automatic ) {
            retro = components( inclinationIndex ) > 0.75 * PI;
        }
        EquinoctialElements elements( convertKeplerianToEquinoctialElements( components, type, retro ), type, retro );
        elements.keplerianElements = components;
        elements.keplerianNeedUpdate = false;
        return elements;
    }

    //! Construct from Cartesian elements
    static EquinoctialElements fromCartesian( const Eigen::Vector6d &components, const double gravitationalParameter,
                                              const orbital_element_conversions::FastVariableType type,
                                              const ElementSet set = automatic ) {
        using namespace orbital_element_conversions;
        using namespace mathematical_constants;
        bool retro = set == retrograde;
        Eigen::Vector6d keplerianElements;
        if ( set == automatic ) {
            keplerianElements = convertCartesianToKeplerianElements( components, gravitationalParameter );
            retro = keplerianElements( inclinationIndex ) > 0.75 * PI;
        }
        EquinoctialElements elements( convertCartesianToEquinoctialElements(
                                          components, gravitationalParameter, type, retro ), type, retro );
        elements.cartesianElements = components;
        elements.cartesianNeedUpdate = false;
        if ( set == automatic ) {
            elements.keplerianElements = keplerianElements;
            elements.keplerianNeedUpdate = false;
        }
        return elements;
    }

    //! Computes Cartesian elements from the current set of equinoctial elements
    Eigen::Vector6d toCartesian( const double gravitationalParameter ) {
        using namespace orbital_element_conversions;
        if ( cartesianNeedUpdate ) {
            cartesianElements = convertEquinoctialToCartesianElements(
                        getComponents(), gravitationalParameter, fastVariableType, retrogradeSet );
            cartesianNeedUpdate = false;
        }
        return cartesianElements;
    }

    //! Computes Keplerian elements from the current set of equinoctial elements
    Eigen::Vector6d toKeplerian( ) {
        using namespace orbital_element_conversions;
        if ( keplerianNeedUpdate ) {
            keplerianElements =
                    convertEquinoctialToKeplerianElements( getComponents(), fastVariableType, retrogradeSet );
            keplerianNeedUpdate = false;
        }
        return keplerianElements;
    }


    //! Set components, with the type of the provided fast variable specified by `type`
    void setComponents( const Eigen::Vector6d &components, const orbital_element_conversions::FastVariableType type,
                        const bool retrograde = false )
    {
        using namespace orbital_element_conversions;
        using namespace basic_mathematics;
        using namespace mathematical_constants;

        a = components( semiMajorAxisIndex );
        h = components( hIndex );
        k = components( kIndex );
        p = components( pIndex );
        q = components( qIndex );

        lmean = TUDAT_NAN;
        lecc  = TUDAT_NAN;
        ltrue = TUDAT_NAN;

        const double fastVariable = computeModulo( components( fastVariableIndex ), 2 * PI );
        switch ( type ) {
        case meanType:
            lmean = fastVariable;
            break;
        case eccentricType:
            lecc  = fastVariable;
            break;
        case trueType:
            ltrue = fastVariable;
            break;
        default:
            throw std::runtime_error( "Unrecognized fast variable type." );
            break;
        }

        fastVariableType = type;
        retrogradeSet = retrograde;

        cartesianNeedUpdate = true;
        keplerianNeedUpdate = true;
    }

    //! Get components, with the type of the returned fast variable specified by `type`
    Eigen::Vector6d getComponents( const orbital_element_conversions::FastVariableType type )
    {
        using namespace orbital_element_conversions;

        double fastVariable;
        switch ( type ) {
        case meanType:
            fastVariable = getMeanLongitude();
            break;
        case eccentricType:
            fastVariable = getEccentricLongitude();
            break;
        case trueType:
            fastVariable = getTrueLongitude();
            break;
        default:
            throw std::runtime_error( "Unrecognized fast variable type." );
            break;
        }

        Eigen::Vector6d components;
        components << a, h, k, p, q, fastVariable;
        return components;
    }

    //! Get components defined by the last call to setComponents()
    Eigen::Vector6d getComponents( )
    {
        return getComponents( fastVariableType );
    }

    //! Get component using paranthesis notation
    double operator() ( int component ) {
        using namespace orbital_element_conversions;

        if ( component == fastVariableIndex ) {
            component = fastVariableType;
        }
        switch ( component ) {
        case semiMajorAxisIndex: return a; break;
        case hIndex:             return h; break;
        case kIndex:             return k; break;
        case pIndex:             return p; break;
        case qIndex:             return q; break;
        case meanType:           return getMeanLongitude();      break;
        case eccentricType:      return getEccentricLongitude(); break;
        case trueType:           return getTrueLongitude();      break;
        default:                 throw  std::runtime_error( "Unrecognized equinoctial elements component."); break;
        }
    }

    bool isRetrograde() {
        return retrogradeSet;
    }

    orbital_element_conversions::FastVariableType getFastVariableType() {
        return fastVariableType;
    }

    //! Default assignment operator
    EquinoctialElements & operator = ( const EquinoctialElements & ) = default;

    //! Comparison operator
    bool operator == ( EquinoctialElements &other ) {
        return getComponents() == other.getComponents() && isRetrograde() == other.isRetrograde();
    }

    //! Negative comparison operator
    bool operator != ( EquinoctialElements &other ) {
        return ! ( *this == other );
    }


private:

    //! Fast variable type
    orbital_element_conversions::FastVariableType fastVariableType;

    //! Set of orbital elements used
    bool retrogradeSet;


    //! Semi-major axis
    double a = TUDAT_NAN;

    //! Sine component of the eccentricity vector
    double h = TUDAT_NAN;

    //! Cosine component of the eccentricity vector
    double k = TUDAT_NAN;

    //! Sine component of the ascending node vector
    double p = TUDAT_NAN;

    //! Cosine component of the ascending node vector
    double q = TUDAT_NAN;


    //! Mean longitude
    double lmean = TUDAT_NAN;

    //! Eccentric longitude
    double lecc = TUDAT_NAN;

    //! True longitude
    double ltrue = TUDAT_NAN;

    //! Compute the mean longitude, if needed, and return its value
    double getMeanLongitude() {
        using namespace orbital_element_conversions;
        if ( lmean != lmean ) {  // lmean is not a number
            if ( lecc != lecc ) {  // lecc is not a number
                if ( ltrue != ltrue ) {  // ltrue is not a number either!
                    throw std::runtime_error( "Mean, eccentric and true longitudes are all undefined. "
                                              "Spacecraft may be re-entering." );
                }
                lecc = getEccentricLongitudeFromTrue( getComponents( trueType) );
            }
            lmean = getMeanLongitudeFromEccentric( getComponents( eccentricType ) );
        }
        return lmean;
    }

    //! Compute the eccentric longitude, if needed, and return its value
    double getEccentricLongitude() {
        using namespace orbital_element_conversions;
        if ( lecc != lecc ) {  // lecc is not a number
            if ( lmean == lmean ) {  // lmean is a number
                lecc = getEccentricLongitudeFromMean( getComponents( meanType) );
            } else {
                if ( ltrue != ltrue ) {  // ltrue is not a number either!
                    throw std::runtime_error( "Mean, eccentric and true longitudes are all undefined. "
                                              "Spacecraft may be re-entering." );
                }
                lecc = getEccentricLongitudeFromTrue( getComponents( trueType) );
            }
        }
        return lecc;
    }

    //! Compute the true longitude, if needed, and return its value
    double getTrueLongitude() {
        using namespace orbital_element_conversions;
        if ( ltrue != ltrue ) {  // ltrue is not a number
            if ( lecc != lecc ) {  // lecc is not a number
                if ( lmean != lmean ) {  // lmean is not a number either!
                    throw std::runtime_error( "Mean, eccentric and true longitudes are all undefined. "
                                              "Spacecraft may be re-entering." );
                }
                lecc = getEccentricLongitudeFromMean( getComponents( meanType) );
            }
            ltrue = getTrueLongitudeFromEccentric( getComponents( eccentricType ) );
        }
        return ltrue;
    }


    bool cartesianNeedUpdate = true;
    Eigen::Vector6d cartesianElements;

    bool keplerianNeedUpdate = true;
    Eigen::Vector6d keplerianElements;

};



typedef boost::shared_ptr< gravitation::SphericalHarmonicsGravitationalAccelerationModelBase< Eigen::Vector3d > >
CentralGravityAM;

typedef boost::shared_ptr< simulation_setup::Body > BodyPtr;

//! Class to be created at the beginning of the DSST propagation, used by all perturbations
class AuxiliaryElements
{
public:

    //! Empty constructor
    AuxiliaryElements( ) { }

    //! Constructor.
    AuxiliaryElements( const double epoch, const EquinoctialElements equinoctialElements,
                       BodyPtr propagatedBody, BodyPtr centralBody, CentralGravityAM centralGravityAM )
        : propagatedBody( propagatedBody), centralBody( centralBody ), centralGravityAM( centralGravityAM )
    {
        updateState( epoch, equinoctialElements );
    }


    //! Pointer to the propagated body
    BodyPtr propagatedBody;

    //! Pointer to the central body
    BodyPtr centralBody;

    //! Pointer to the central gravity model
    CentralGravityAM centralGravityAM;

    /*
    //! Reference to the central body
    boost::shared_ptr< Body > propagatedBody;

    //! Reference to the central body
    boost::shared_ptr< CelestialBody > centralBody;
    */

    //! Date, in seconds since J2000
    double epoch;

    //! Mean equinoctial elements
    EquinoctialElements equinoctialElements;


    //! Update to a new state.
    void updateState( const double epoch, const EquinoctialElements elements )
    {
        this->epoch = epoch;
        equinoctialElements = elements;
        updateMembers();
    }

    //! Update elements.
    void updateElements( EquinoctialElements elements ) {
        updateState( epoch, elements );
    }


    //! Update to a new state optionally changing fast variable type and/or elements set type.
    void updateState( const double epoch, const Eigen::Vector6d components,
                      const orbital_element_conversions::FastVariableType type, const bool retro )
    {
        updateState( epoch, EquinoctialElements( components, type, retro ) );
    }

    //! Update to a new state without changing fast variable type.
    void updateState( const double epoch, const Eigen::Vector6d components,
                      const orbital_element_conversions::FastVariableType type )
    {
        updateState( epoch, components, type, equinoctialElements.isRetrograde() );
    }

    //! Update to a new state without changing fast variable type or elements set type.
    void updateState( const double epoch, const Eigen::Vector6d components )
    {
        updateState( epoch, components, equinoctialElements.getFastVariableType(),
                     equinoctialElements.isRetrograde() );
    }


    //! Update elements's components optionally changing fast variable type and/or elements set type.
    void updateComponents( const Eigen::Vector6d components,
                                   const orbital_element_conversions::FastVariableType type, const bool retro ) {
        updateState( epoch, components, type, retro );
    }

    //! Update elements's components without changing fast variable type or elements set type.
    void updateComponents( const Eigen::Vector6d components ) {
        updateComponents( components, equinoctialElements.getFastVariableType(),
                                  equinoctialElements.isRetrograde() );
    }


    //! Add vector to the current element's components and update epoch
    void updateStateByAdding( const double epoch, const Eigen::Vector6d vector,
              const orbital_element_conversions::FastVariableType type = orbital_element_conversions::meanType ) {
        updateState( epoch, equinoctialElements.getComponents( type ) + vector, type );
    }

    //! Add vector to the current element's components, but don't unpdate epoch
    void updateComponentsByAdding( const Eigen::Vector6d vector,
              const orbital_element_conversions::FastVariableType type = orbital_element_conversions::meanType ) {
        updateStateByAdding( epoch, vector, type );
    }



private:

    //! Update the instance's members according to the current equinoctialElements.
    void updateMembers();


public:

    // Properties that depend on a referenced object (i.e. bodies) may change with time. However, they need to be
    // assumed constant over an integration step in order to be able to average the equations of motion.
    // Thus, they are stored (computed only once per step, used many times) and updated by updateMembers().

    /*
    //! Propagated body's mass
    double mass;

    //! Propagated body's cross-sectional area
    double area;

    //! Propagated body's drag coefficient
    double CD;

    //! Propagated body's radiation pressure coefficient
    // double CR;
    */

    //! Central body's gravitational parameter
    double mu;


    // Properties that depend on a stored object (i.e. equinoctialElements) are stored and updated by updateMembers().
    // In this way, they are only computed once for each integration step, even though they are used many times.

    //! Retrograde factor I
    int I;

    //! Eccentricity
    double e;

    //! Keplerian mean motion
    double meanMotion;

    //! Keplerian period
    double period;

    //! Semi-major axis
    double a;

    //! y component of eccentricity vector
    double h;

    //! x component of eccentricity vector
    double k;

    //! y component of inclination vector
    double p;

    //! x component of inclination vector
    double q;

    /*
    //! Mean longitude
    double meanLongitude() {
        return equinoctialElements( orbital_element_conversions::meanType );
    }

    //! Eccentric longitude
    double eccentricLongitude() {
        return equinoctialElements( orbital_element_conversions::eccentricType );
    }

    //! True longitude
    double trueLongitude() {
        return equinoctialElements( orbital_element_conversions::trueType );
    }
    */

    //! A = sqrt(μ * a)
    double A;

    //! B = sqrt(1 - h² - k²)
    double B;

    //! C = 1 + p² + q²
    double C;

    //! Equinoctial frame f vector
    Eigen::Vector3d f;

    //! Equinoctial frame g vector
    Eigen::Vector3d g;

    //! Equinoctial frame w vector
    Eigen::Vector3d w;

    //! Vns coefficients generator, used by all conservative force models
    coefficients_factories::VnsCoefficientsFactory VnsFactory;

    //! Gaussian quadrature, used by all non-conservative force models
    numerical_quadrature::GaussianQuadrature< double, Eigen::Vector6d > gaussianQuadrature;


    //! Default assignment operator
    AuxiliaryElements & operator= ( const AuxiliaryElements & ) = default;

};



} // namespace sst

} // namespace propagators

} // namespace tudat

#endif // TUDAT_PROPAGATORS_DSST_AUXILIARYELEMENTS_H
