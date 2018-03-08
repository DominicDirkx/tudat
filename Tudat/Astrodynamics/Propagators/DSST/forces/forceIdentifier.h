#ifndef TUDAT_PROPAGATORS_DSST_FORCEMODELS_FORCEIDENTIFIER_H
#define TUDAT_PROPAGATORS_DSST_FORCEMODELS_FORCEIDENTIFIER_H

#include "Tudat/Astrodynamics/BasicAstrodynamics/accelerationModelTypes.h"

namespace tudat
{

namespace propagators
{

namespace sst
{

namespace force_models
{

//! Struct for storing settings for ForceModel
//! To be extended by derived classes
struct ForceModelSettings {
    //! Empty constructor
    ForceModelSettings( ) { }

    //! Destructor (makes struct polymorphic)
    virtual ~ForceModelSettings() { }
};


//! Force identifier defined by the name of the body exerting the acceleration and the type of acceleration
struct ForceIdentifier
{
public:
    //! Empty constructor.
    //! By default, body is "" and type is undefined_acceleration, which is used by DSST propagator to represent all
    //! accelerations caused by all perturbing bodies, storing the total values in the key ForceIdentifier() of maps.
    ForceIdentifier( ) { }

    //! Constructor with a body name and acceleration type.
    ForceIdentifier( const std::string &bodyExertingAcceleration,
                     const basic_astrodynamics::AvailableAcceleration accelerationType ) :
        body( bodyExertingAcceleration ), type( accelerationType ) { }

    //! Name of the body exerting the acceleration
    std::string body = "";

    //! Type of the force model
    basic_astrodynamics::AvailableAcceleration type =
            basic_astrodynamics::AvailableAcceleration::undefined_acceleration;

    //! Description of the identifier, a string "body-type"
    std::string description( ) const
    {
        return body + "-" + boost::lexical_cast< std::string >( type );
    }

    //! Comparison operator
    bool operator == ( const ForceIdentifier& other ) const
    {
        return body == other.body && type == other.type;
    }

    //! Negative comparison operator
    bool operator != ( const ForceIdentifier& other ) const
    {
        return !( *this == other );
    }
};

//! Force identifier representing all force models
extern ForceIdentifier ALL_FORCES;

} // namespace force_models

} // namespace sst

} // namespace propagators

} // namespace tudat


namespace std
{

//! Overload << in order to be able to use a ForceIdentifier in std::cout
ostream& operator<< ( ostream& stream, const tudat::propagators::sst::force_models::ForceIdentifier& forceID );

//! Define hash for ForceIdentifier (necessary in order to be able use ForceIdentifier as an std::unordered_map's key)
template < >
struct hash< tudat::propagators::sst::force_models::ForceIdentifier >
{
    std::size_t operator() ( const tudat::propagators::sst::force_models::ForceIdentifier& forceID ) const
    {
        // Start with a hash value of 0    .
        std::size_t seed = 0;

        // Modify 'seed' by XORing and bit-shifting in
        // one member of 'ForceIdentifier' after the other:
        boost::hash_combine( seed, boost::hash_value( forceID.body ) );
        boost::hash_combine( seed, boost::hash_value( forceID.type ) );

        // Return the result.
        return seed;
    }
};

}


#endif // TUDAT_PROPAGATORS_DSST_FORCEMODELS_FORCEIDENTIFIER_H
