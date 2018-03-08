#ifndef TUDAT_PROPAGATORS_DSST_COEFFICIENTSFACTORIES_H
#define TUDAT_PROPAGATORS_DSST_COEFFICIENTSFACTORIES_H

#include "Tudat/Astrodynamics/Propagators/DSST/utilities/vectors.h"

namespace tudat
{

namespace propagators
{

namespace sst
{

namespace coefficients_factories
{


/*
bool operator == ( const CoefficientsFactory &one, const CoefficientsFactory &another ) {
    return one.isEqual( another );
}

bool operator != ( const CoefficientsFactory &one, const CoefficientsFactory &another ) {
    return !( one == another );
}
*/


class SingleOrderCoefficients : public Vectord {
public:
    //! Inherit constructors.
    using Vectord::Vector;

    /*
    //! Empty constructor
    SingleOrderCoefficients( ) { }

    //! Constructor with list of values
    SingleOrderCoefficients( const std::initializer_list< double > &list ) : Vectord( list ) { }

    // OrderCoefficients( const double initialValue ) : Vectord( 1, initialValue ) { }

    // OrderCoefficients & operator= ( const OrderCoefficients & ) = default;
    */

    int order() const {
        return size() - 1;
    }

    void resizeForOrder( const int unsigned order ) {
        if ( int( order ) > this->order() ) {
            resize( order + 1 );
        }
    }

    Vectord getUpTo( const unsigned int order ) const {
        return std::vector< double >( begin(), begin() + order + 1 );
    }

    //! Returns `this` as a Vectord
    Vectord getAll() const {
        return *this;
    }

    //! Default assignment operator
    SingleOrderCoefficients & operator= ( const SingleOrderCoefficients & ) = default;

};


class DoubleOrderCoefficients : public NestedVectord {
public:
    //! Inherit constructors.
    using NestedVectord::NestedVector;

    /*
    //! Empty constructor
    DoubleOrderCoefficients( ) { }

    //! Constructor with list of lists
    DoubleOrderCoefficients( const std::initializer_list< std::initializer_list< double > > &list )
        : NestedVectord( list ) { }

    // OrderCoefficients( const double initialValue ) : Vectord( 1, initialValue ) { }

    // OrderCoefficients & operator= ( const OrderCoefficients & ) = default;
    */

    int order() const {
        return size() - 1;
    }

    int suborder( const unsigned int i ) const {
        return ( *this ).at( i ).size() - 1;
    }

    void resizeForOrder( const unsigned int order ) {
        if ( int( order ) > this->order() ) {
            resize( order + 1 );
        }
    }

    void resizeForSuborder( const unsigned int i, const unsigned int suborder ) {
        if ( int( suborder ) > this->suborder( i ) ) {
            ( *this ).at( i ).resize( suborder + 1 );
        }
    }

    NestedVectord getUpTo( const unsigned int order, const unsigned int suborder ) const {
        NestedVectord coefficients;
        for ( unsigned int i = 0; i < size(); i++ ) {
            if ( i > order ) {
                break;
            }
            Vectord subcoefficients;
            Vectord v = ( *this ).at( i );
            for ( unsigned int j = 0; j < v.size(); j++ ) {
                if ( j > suborder ) {
                    break;
                }
                subcoefficients.push_back( v[ j ] );
            }
            coefficients.push_back( subcoefficients );
        }
        return coefficients;
    }

    //! Returns `this` as a NestedVectord
    NestedVectord getAll() const {
        return *this;
    }

    //! Default assignment operator
    DoubleOrderCoefficients & operator= ( const DoubleOrderCoefficients & ) = default;

};


class CoefficientsFactory
{
public:

    //! Default constructor
    CoefficientsFactory( ) { }

    //! Constructor with parameters
    CoefficientsFactory( const std::map< std::string, double > &parameters ) : parameters( parameters ) { }

    //! Add a parameter or change the value of an existing parameter
    virtual void setParameter( const std::string name, const double value ) {
        if ( parameters.count( name ) > 0 ) {
            if ( parameters[ name ] == value ) {
                return; // if parameter already exists and has the same value, do nothing
            }
        }
        parameters[ name ] = value;
        resetAfterParameterChange( name );
    }

    /*
    virtual bool isEquivalent( const CoefficientsFactory &other ) const {
        return parameters == other.parameters;
    }
    */

    //! Default assignment operator
    CoefficientsFactory & operator= ( const CoefficientsFactory & ) = default;


protected:
    //! Set initial values of the coefficients necessary for recursion formulas and resetting
    virtual void initializeCoefficients() { }

    //! Set initial values of the dependent coefficients necessary for recursion formulas and resetting
    virtual void initializeDependentCoefficients() { }


    //! Compute the coefficients
    /**
     * Compute the coefficients
     * @param extended Whether to compute the coefficients up to an extended order necessary for the computation of the
     * dependent coefficients. Default is false, i.e. the coefficients are only computed up to `order`.
     */
    virtual void computeCoefficients( const bool extended = false ) { }

    //! Compute the dependent coefficients
    virtual void computeDependentCoefficients() { }


    //! Whether the dependent coefficients depend on the coefficients
    virtual bool dependentCoefficientsDependOnCoefficients() {
        return true;
    }

    //! Updates the coefficients
    virtual void updateCoefficientsIfNeeded( const bool extended = false ) {
        setInitialValuesForRecursions();
        computeCoefficients( extended );
    }

    //! Updates the dependent coefficients.
    /// Updates the coefficients first if dependentCoefficientsDependOnCoefficients() returns true.
    virtual void updateDependentCoefficientsIfNeeded() {
        setInitialValuesForRecursions();
        if ( dependentCoefficientsDependOnCoefficients() ) {
            updateCoefficientsIfNeeded( true );
        }
        computeDependentCoefficients();
    }


    //! Reset the coefficients
    virtual void resetCoefficients() {
        initializeCoefficients();
    }

    //! Reset the dependent coefficients
    virtual void resetDependentCoefficients() {
        initializeDependentCoefficients();
    }


    //! Whether the coefficients are reset when the value of a parameter named `name` is changed
    virtual bool changingParameterResetsCoefficients( std::string name ) const {
        return true;
    }

    //! Whether the dependent coefficients are reset when the value of a parameter named `name` is changed
    virtual bool parameterChangeResetsDependentCoefficients( std::string name ) const {
        return true;
    }

    //! Optionally reset the (dependent) coefficients after changing the value of a parameter named `name`
    void resetAfterParameterChange( std::string parameterName ) {
        if ( changingParameterResetsCoefficients( parameterName ) ) {
            resetCoefficients();
        }
        if ( parameterChangeResetsDependentCoefficients( parameterName ) ) {
            resetDependentCoefficients();
        }
    }


    //! Get the value of a parameter named `name`
    double getParameter( const std::string name ) const {
        if ( parameters.count( name ) > 0 ) {
            return parameters.at( name );
        } else {
            std::cout << "Trying to access parameter named \"" << name << "\"." << std::endl;
            throw std::runtime_error( "CoefficientsFactory::getParameter: the requested parameter does not exist." );
        }
    }


private:

    std::map< std::string, double > parameters;

    bool initialValuesForRecursionSet = false;

    void setInitialValuesForRecursions() {
        if ( ! initialValuesForRecursionSet ) {
            initializeCoefficients();
            initializeDependentCoefficients();
            initialValuesForRecursionSet = true;
        }
    }

};


class SingleOrderCoefficientsFactory : public CoefficientsFactory
{
public:

    //! Empty constructor
    SingleOrderCoefficientsFactory( ) : CoefficientsFactory() { }

    //! Constructor with order and no parameters
    SingleOrderCoefficientsFactory( const unsigned int order ) : CoefficientsFactory(), order( order ) { }

    //! Constructor with parameters and order
    SingleOrderCoefficientsFactory( const std::map< std::string, double > &parameters, const unsigned int order ) :
        CoefficientsFactory( parameters ), order( order ) { }

    /*
    using CoefficientsFactory::isEquivalent;
    virtual bool isEquivalent( const SingleOrderCoefficientsFactory &other ) const {
        return CoefficientsFactory::isEquivalent( other) && order == other.order;
    }
    */

    //! Change the current order
    void setOrder( const unsigned int order ) {
        this->order = order;
    }

    //! Default assignment operator
    SingleOrderCoefficientsFactory & operator= ( const SingleOrderCoefficientsFactory & ) = default;


protected:

    //! Current order up to which the (dependent) coefficients have to be computed
    unsigned int order = 0;

    //! Maximum order up to which the coefficients have been computed
    unsigned int maximumComputedOrder = 0;

    //! Maximum order up to which the dependent coefficients have been computed
    unsigned int maximumComputedOrderDependentCoefficients = 0;

    //! Order for the coefficients necessary to compute the dependent coefficients
    virtual unsigned int getExtendedOrder() const {
        return order;
    }

    //! Get the current (extended) order
    /**
     * Get the current (extended) order
     * @param extended Whether the order given by getExtendedOrder() should be returned.
     * Default is false, i.e. `order` is returned.
     * @return the current (extended) order
     */
    unsigned int getOrder( const bool extended = false ) const {
        return extended ? this->getExtendedOrder() : order;
    }


    virtual void updateCoefficientsIfNeeded( const bool extended = false ) {
        CoefficientsFactory::updateCoefficientsIfNeeded( extended );
        maximumComputedOrder = std::max( maximumComputedOrder, getOrder( extended ) );
    }

    virtual void updateDependentCoefficientsIfNeeded() {
        CoefficientsFactory::updateDependentCoefficientsIfNeeded();
        maximumComputedOrderDependentCoefficients = std::max( maximumComputedOrderDependentCoefficients, order );
    }


    virtual void resetCoefficients() {
        CoefficientsFactory::resetCoefficients();
        maximumComputedOrder = 0;
    }

    virtual void resetDependentCoefficients() {
        CoefficientsFactory::resetDependentCoefficients();
        maximumComputedOrderDependentCoefficients = 0;
    }

};


class DoubleOrderCoefficientsFactory : public SingleOrderCoefficientsFactory
{
public:

    //! Empty constructor
    DoubleOrderCoefficientsFactory( ) : SingleOrderCoefficientsFactory() { }

    //! Constructor with order and suborder and no parameters
    DoubleOrderCoefficientsFactory( const unsigned int order, const unsigned int suborder ) :
        SingleOrderCoefficientsFactory( order ), suborder( suborder ) { }

    //! Constructor with parameters, order and suborder
    DoubleOrderCoefficientsFactory( const std::map< std::string, double > &parameters,
                                    const unsigned int order, const unsigned int suborder ) :
        SingleOrderCoefficientsFactory( parameters, order ), suborder( suborder ) { }

    /*
    using SingleOrderCoefficientsFactory::isEquivalent;
    virtual bool isEquivalent( const DoubleOrderCoefficientsFactory &other ) const {
        return SingleOrderCoefficientsFactory::isEquivalent( other ) && suborder == other.suborder;
    }
    */

    //! Change the current suborder
    void setSuborder( const unsigned int suborder ) {
        this->suborder = suborder;
    }

    //! Default assignment operator
    DoubleOrderCoefficientsFactory & operator= ( const DoubleOrderCoefficientsFactory & ) = default;


protected:

    //! Maximum current suborder up to which the (dependent) coefficients have to be computed
    unsigned int suborder = 0;

    // unsigned int maximumComputedSuborder = 0;

    // unsigned int suborderForMaximumComputedOrder = 0;
    // unsigned int orderForMaximumComputedSuborder = 0;


    //! Get the maximum suborder for order `i`
    virtual unsigned int getSuborder( const unsigned int i ) const {
        return suborder;
    }

    //! Get the extended maximum suborder for order `i`
    virtual unsigned int getExtendedSuborder( const unsigned int i ) const {
        return getSuborder( i );
    }

    //! Get the (extended) maximum suborder for order `i`
    unsigned int getSuborderExtended( const unsigned int i, const bool extended ) const {
        return extended ? getExtendedSuborder( i ) : getSuborder( i );
    }


    //! Inherit update methods
    using SingleOrderCoefficientsFactory::updateCoefficientsIfNeeded;
    using SingleOrderCoefficientsFactory::updateDependentCoefficientsIfNeeded;

    /*
    virtual void updateCoefficientsIfNeeded( const bool extended = false ) {
        SingleOrderCoefficientsFactory::updateCoefficientsIfNeeded( extended );
    }

    virtual void updateDependentCoefficientsIfNeeded() {
        SingleOrderCoefficientsFactory::updateDependentCoefficientsIfNeeded();
    }
    */


};


//! Factorials generator
class FactorialsFactory : public SingleOrderCoefficientsFactory
{
public:

    //! Empty constructor
    FactorialsFactory( ) : SingleOrderCoefficientsFactory() { }

    //! Constructor with order
    FactorialsFactory( const unsigned int order ) : SingleOrderCoefficientsFactory( order ) { }


    //! Returns x!
    double operator() ( const unsigned int x ) {
        if ( x > maximumComputedOrder ) {
            setOrder( x );
            updateCoefficientsIfNeeded();
        }
        return factorials( x );
    }

    //! Default assignment operator
    FactorialsFactory & operator= ( const FactorialsFactory & ) = default;


private:

    SingleOrderCoefficients factorials;

    void initializeCoefficients();
    void computeCoefficients( const bool extended = false );

};


class GsHsCoefficientsFactory : public SingleOrderCoefficientsFactory
{
public:

    //! Empty constructor
    GsHsCoefficientsFactory( ) : SingleOrderCoefficientsFactory() { }

    //! Constructor with parameters and order
    /**
     * Constructor with parameters and order
     * @param h Second component of equinoctial elements
     * @param k Third component of equinoctial elements
     * @param alpha First cosine director
     * @param beta Second cosine director
     * @param order Order up to which the coefficients should be computed
     */
    GsHsCoefficientsFactory( const double h, const double k, const double alpha, const double beta,
                             const unsigned int order ) : SingleOrderCoefficientsFactory(
    { { "h", h }, { "k", k }, { "alpha", alpha }, { "beta", beta } }, order ) { }

    //! G_s coefficients
    Vectord getGsCoefficients() {
        updateCoefficientsIfNeeded();
        return Gs.getUpTo( order );
    }

    //! H_s coefficients
    Vectord getHsCoefficients() {
        updateCoefficientsIfNeeded();
        return Hs.getUpTo( order );
    }

    //! Partial derivatives of the G_s coefficients with respecto to h
    Vectord getGsHDerivatives() {
        updateDependentCoefficientsIfNeeded();
        return dGs[ "h" ].getUpTo( order );
    }

    //! Partial derivatives of the G_s coefficients with respecto to k
    Vectord getGsKDerivatives() {
        updateDependentCoefficientsIfNeeded();
        return dGs[ "k" ].getUpTo( order );
    }

    //! Partial derivatives of the G_s coefficients with respecto to \alpha
    Vectord getGsAlphaDerivatives() {
        updateDependentCoefficientsIfNeeded();
        return dGs[ "alpha" ].getUpTo( order );
    }

    //! Partial derivatives of the G_s coefficients with respecto to \beta
    Vectord getGsBetaDerivatives() {
        updateDependentCoefficientsIfNeeded();
        return dGs[ "beta" ].getUpTo( order );
    }

    //! Default assignment operator
    GsHsCoefficientsFactory & operator= ( const GsHsCoefficientsFactory & ) = default;


private:

    SingleOrderCoefficients Gs;
    SingleOrderCoefficients Hs;
    std::map< std::string, SingleOrderCoefficients > dGs;

    void initializeCoefficients();
    void initializeDependentCoefficients();
    void computeCoefficients( const bool extended = false );
    void computeDependentCoefficients();

};


class QnsCoefficientsFactory : public DoubleOrderCoefficientsFactory
{
public:

    //! Empty constructor
    QnsCoefficientsFactory( ) : DoubleOrderCoefficientsFactory() { }

    //! Constructor with parameters, order and suborder
    /**
     * Constructor with parameters, order and suborder
     * @param gamma Third cosine director
     * @param order Order up to which the coefficients should be computed
     * @param suborder Suborder up to which the coefficients should be computed
     */
    QnsCoefficientsFactory( const double gamma, const unsigned int order, const unsigned int suborder )
        : DoubleOrderCoefficientsFactory( { { "gamma", gamma } }, order, suborder ) { }

    //! Q_{n,s} coefficients
    NestedVectord getCoefficients() {
        updateCoefficientsIfNeeded();
        return Qns.getUpTo( order, suborder );
    }

    //! Derivatives of the Q_{n,s} coefficients with respect to \gamma
    NestedVectord getDerivatives() {
        updateDependentCoefficientsIfNeeded();
        return dQns.getUpTo( order, suborder );
    }

    //! Default assignment operator
    QnsCoefficientsFactory & operator= ( const QnsCoefficientsFactory & ) = default;


private:

    DoubleOrderCoefficients Qns;
    DoubleOrderCoefficients dQns;

    void initializeCoefficients();
    void initializeDependentCoefficients();
    unsigned int getSuborder( const unsigned int i ) const;
    void computeCoefficients( const bool extended = false );
    void computeDependentCoefficients();

};


class K0nsCoefficientsFactory : public DoubleOrderCoefficientsFactory
{
public:

    //! Empty constructor
    K0nsCoefficientsFactory( ) : DoubleOrderCoefficientsFactory() { }

    //! Constructor with parameters, order and suborder
    /**
     * Constructor with parameters, order and suborder
     * @param Chi A paremeter depending on the eccentricity by 1 / sqrt(1 - e²)
     * @param order Order up to which the coefficients should be computed
     * @param suborder Suborder up to which the coefficients should be computed
     */
    K0nsCoefficientsFactory( const double Chi, const unsigned int order, const unsigned int suborder )
        : DoubleOrderCoefficientsFactory( { { "Chi", Chi } }, order, suborder ) { }

    //! K^0_{n,s} coefficients, also known as kernels of Hansen coefficients
    NestedVectord getCoefficients() {
        updateCoefficientsIfNeeded();
        return K0ns.getUpTo( order, suborder );
    }

    //! Derivatives of the K^0_{n,s} coefficients with respect to \Chi
    NestedVectord getDerivatives() {
        updateDependentCoefficientsIfNeeded();
        return dK0ns.getUpTo( order, suborder );
    }

    //! Default assignment operator
    K0nsCoefficientsFactory & operator= ( const K0nsCoefficientsFactory & ) = default;


private:

    DoubleOrderCoefficients K0ns;
    DoubleOrderCoefficients dK0ns;

    void initializeCoefficients();
    void initializeDependentCoefficients();
    unsigned int getSuborder( const unsigned int i ) const;
    void computeCoefficients( const bool extended = false );
    void computeDependentCoefficients();

};


class NegativeK0nsCoefficientsFactory : public DoubleOrderCoefficientsFactory
{
public:

    //! Empty constructor
    NegativeK0nsCoefficientsFactory( ) : DoubleOrderCoefficientsFactory() { }

    //! Constructor with parameters, order and suborder
    /**
     * Constructor with parameters, order and suborder
     * @param Chi A paremeter depending on the eccentricity by 1 / sqrt(1 - e²)
     * @param order Order up to which the coefficients should be computed
     * @param suborder Suborder up to which the coefficients should be computed
     */
    NegativeK0nsCoefficientsFactory( const double Chi, const unsigned int order, const unsigned int suborder )
        : DoubleOrderCoefficientsFactory( { { "Chi", Chi } }, order, suborder ) { }

    //! K^0_{-n-1,s} coefficients, also known as kernels of Hansen coefficients
    NestedVectord getCoefficients() {
        updateCoefficientsIfNeeded();
        return mK0ns.getUpTo( order, suborder );
    }

    //! Derivatives of the K^0_{-n-1,s} coefficients with respect to \Chi
    NestedVectord getDerivatives() {
        updateDependentCoefficientsIfNeeded();
        return dmK0ns.getUpTo( order, suborder );
    }

    //! Default assignment operator
    NegativeK0nsCoefficientsFactory & operator= ( const NegativeK0nsCoefficientsFactory & ) = default;


private:

    DoubleOrderCoefficients mK0ns;
    DoubleOrderCoefficients dmK0ns;

    void initializeCoefficients();
    void initializeDependentCoefficients();
    unsigned int getSuborder( const unsigned int i ) const;
    void computeCoefficients( const bool extended = false );
    void computeDependentCoefficients();

};


class VnsCoefficientsFactory : public DoubleOrderCoefficientsFactory
{
public:

    //! Empty constructor
    VnsCoefficientsFactory( ) : DoubleOrderCoefficientsFactory() { }

    //! Constructor with order and suborder
    /**
     * Constructor with parameters, order and suborder
     * @param order Order up to which the coefficients should be computed
     * @param suborder Suborder up to which the coefficients should be computed
     */
    VnsCoefficientsFactory( const unsigned int order, const unsigned int suborder )
        : DoubleOrderCoefficientsFactory( order, suborder ) { }

    //! Coefficients V_{n,s}
    NestedVectord getCoefficients() {
        updateCoefficientsIfNeeded();
        return Vns.getUpTo( order, suborder );
    }

    //! Coefficients V_{n,s}
    NestedVectord getCoefficientsUpTo( const unsigned int order, const unsigned int suborder ) {
        setOrder( order );
        setSuborder( suborder );
        return getCoefficients();
    }

    //! Coefficients V^m_{n,s}
    NestedVectord getCoefficientsWithSuperscript( const unsigned int m ) {
        setParameter( "m", m );
        updateDependentCoefficientsIfNeeded();
        return Vmns.getUpTo( getOrder( true ), suborder );
    }

    //! Default assignment operator
    VnsCoefficientsFactory & operator= ( const VnsCoefficientsFactory & ) = default;


private:

    DoubleOrderCoefficients Vns;
    DoubleOrderCoefficients Vmns;

    void initializeCoefficients();
    void initializeDependentCoefficients();
    bool changingParameterResetsCoefficients( std::string name ) const;
    unsigned int getExtendedOrder() const;
    unsigned int getSuborder( const unsigned int i ) const;
    void computeCoefficients( const bool extended = false );
    void computeDependentCoefficients();

};


//! Class for Jacobi polynomials keys [v,w].
class JacobiKey
{
public:

    /** Simple constructor.
     * \param v first exponent
     * \param w second exponent
     */
    JacobiKey( const int v, const int w ) : v(v), w(w) { }


    //! First exponent
    const int v;

    //! Second exponent
    const int w;


    /** Check if the instance represent the same key as another instance.
     * \param key other key
     * \return true if the instance and the other key refer to the same polynomial
     */
    bool operator == (const JacobiKey key) const {
        return v == key.v && w == key.w;
    }

    /** Comparison operator necessary to use JacobiKey as the key of std::map.
     * \param key other key
     * \return true if the instance is "smaller" than another key
     */
    bool operator < ( const JacobiKey& key ) const {
        return v == key.v ? w < key.w : v < key.v;
    }

};


//! Provider of the Jacobi polynomials P_l^{v,w}(\gamma)
class JacobiPolynomialsFactory : public SingleOrderCoefficientsFactory
{
public:

    //! Empty constructor
    JacobiPolynomialsFactory( ) : SingleOrderCoefficientsFactory() { }

    //! Constructor with parameters and order
    /**
     * Constructor with parameters, order and suborder
     * @param gamma Third cosine director
     * @param v First exponent
     * @param w Second exponent
     * @param order Order up to which the coefficients should be computed
     */
    JacobiPolynomialsFactory( const double gamma, const int v, const int w, const unsigned int l ) :
        SingleOrderCoefficientsFactory( { { "gamma", gamma }, { "v", v }, { "w", w } }, l ) { }

    JacobiPolynomialsFactory( const double gamma, const JacobiKey &key, const unsigned int l ) :
        JacobiPolynomialsFactory( gamma, key.v, key.w, l ) { }

    //! P_l^{v,w}(\gamma) polynomial values
    Vectord getPolynomials() {
        updateCoefficientsIfNeeded();
        return P.getUpTo( order );
    }

    //! P_l^{v,w}(\gamma) polynomial value
    double getPolynomial( const unsigned int l ) {
        if ( l > maximumComputedOrder ) {
            setOrder( l );
            updateCoefficientsIfNeeded();
        }
        return P( l );
    }

    //! Default assignment operator
    JacobiPolynomialsFactory & operator= ( const JacobiPolynomialsFactory & ) = default;


private:

    SingleOrderCoefficients P;

    void initializeCoefficients();
    void computeCoefficients( const bool extended = false );

};



//! Provider of the coefficients C_s(a,b) and S_s(a,b)
class CsSsCoefficientsFactory : public SingleOrderCoefficientsFactory
{
public:

    //! Empty constructor
    CsSsCoefficientsFactory( ) : SingleOrderCoefficientsFactory() { }

    //! Constructor with parameters and order
    /**
     * Constructor with parameters, order and suborder
     * @param a First argument
     * @param b Second argument
     * @param order Order up to which the coefficients should be computed
     */
    CsSsCoefficientsFactory( const double a, const double b, const unsigned int order ) :
        SingleOrderCoefficientsFactory( { { "a", a }, { "b", b } }, order ) { }

    //! C_s coefficients
    Vectord getCs() {
        updateCoefficientsIfNeeded();
        return Cs.getUpTo( order );
    }

    //! S_s coefficients
    Vectord getSs() {
        updateCoefficientsIfNeeded();
        return Ss.getUpTo( order );
    }

    //! Default assignment operator
    CsSsCoefficientsFactory & operator= ( const CsSsCoefficientsFactory & ) = default;


private:

    SingleOrderCoefficients Cs;
    SingleOrderCoefficients Ss;

    void initializeCoefficients();
    void computeCoefficients( const bool extended = false );

};



} // namespace coefficients_factories


//! Factorials generator declared at a global scope.
//! Can be accessed from anywhere within tudat::propagators::sst after including "coefficientsFactory.h".
//! To get e.g. 5! simply call factorial(5).
//! If then factorial(4) is called, no computations will be done because it already has been computed and stored.
//! If then factorial(8) is called, the computations will start from 6!, saving time.
//! Specially useful during iterative procedures in which large number factorials (eg. 40!) are needed many times.
extern coefficients_factories::FactorialsFactory factorial;


} // namespace sst

} // namespace propagators

} // namespace tudat

#endif // TUDAT_PROPAGATORS_DSST_COEFFICIENTSFACTORIES_H
