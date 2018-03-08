#ifndef TUDAT_PROPAGATORS_DSST_VECTORS_H
#define TUDAT_PROPAGATORS_DSST_VECTORS_H


#include <cmath>
#include <map>
#include <vector>
#include <set>

#include "Tudat/Astrodynamics/BasicAstrodynamics/orbitalElementConversions.h"
#include "Tudat/Basics/basicTypedefs.h"
#include "Tudat/Mathematics/BasicMathematics/mathematicalConstants.h"

#include "boost/boost/function.hpp"


namespace tudat
{

namespace propagators
{

namespace sst  // Semianalytic Satellite Theory
{


//! Support for vectors (Vector), vectors of vectors (NestedVector) and vectors of vectors of vectors (TripleVector).
//! Used throughout the DSST theory, thus defined at tudat::propagators::sst scope.

template< class T >
class Vector : public std::vector< T > {
public:
    //! Inherit std::vector< T > constructors.
    using std::vector< T >::vector;

    //! Default constructor.
    Vector( ) { }

    //! Constructor from vector.
    Vector( const std::vector< T > &vector ) : std::vector< T >( vector ) { }

    /*
    //! Constructor with list of values
    Vector( const std::initializer_list< T > &list ) : std::vector< T >( list ) { }

    //! Constructor with size
    Vector( const unsigned int m ) : std::vector< T >( m ) { }

    //! Constructor with repeated value
    Vector( const unsigned int m, const T &repeated ) : std::vector< T >( m, repeated ) { }

    Vector & operator= ( const Vector & ) = default;
    */

    T& operator() ( const unsigned int i ) {
        return ( *this ).at( i );
    }

    const T& operator() ( const unsigned int i ) const {
        return ( *this )[ i ];
    }

};
typedef Vector< double > Vectord;
typedef Vector< int >    Vectori;


template< class T >
class NestedVector : public std::vector< std::vector< T > > {
public:
    //! Inherit std::vector< std::vector< T > > constructors.
    using std::vector< std::vector< T > >::vector;

    //! Default constructor
    NestedVector( ) { }

    //! Constructor with vector of vectors
    NestedVector( const std::vector< std::vector< T > > &vectorOfVectors ) :
        std::vector< std::vector< T > >( vectorOfVectors ) { }


    //! Constructor with list of lists
    NestedVector( const std::initializer_list< std::initializer_list< T > > &listOfLists ) {
        for ( auto list : listOfLists ) {
            this->push_back( std::vector< T >( list ) );
        }
    }

    //! Constructor with size m for the first level, size n for the second level
    NestedVector( const unsigned int m, const unsigned int n ) :
        std::vector< std::vector< T > >( m, std::vector< T >( n ) ) { }

    //! Constructor from Vector: { x, y, z } -> { {x}, {y}, {z} }
    NestedVector( const std::vector< T > &vector ) {
        for ( T elem : vector ) {
            this->push_back( { elem } );
        }
    }


    /*
    //! Constructor with list of vectors
    NestedVector( const std::initializer_list< std::vector< T > > &listOfVectors ) {
        for ( auto vector : listOfVectors ) {
            this->push_back( vector );
        }
    }

    //! Constructor with size m for the first level, empty vectors for the second level
    NestedVector( const unsigned int m ) : std::vector< std::vector< T > >( m, std::vector< T >() ) { }

    //! Constructor with repeated vector
    NestedVector( const unsigned int m, const std::vector< T > &repeated ) :
        std::vector< std::vector< T > >( m, repeated ) { }

    NestedVector & operator= ( const NestedVector & ) = default;
    */

    T& operator() ( const unsigned int i, const unsigned int j ) {
        return ( *this ).at( i ).at( j );
    }

    const T& operator() ( const unsigned int i, const unsigned int j ) const {
        return ( *this )[ i ][ j ];
    }

};
typedef NestedVector< double > NestedVectord;
typedef NestedVector< int >    NestedVectori;


template< class T >
class TripleVector : public std::vector< std::vector< std::vector< T > > > {
public:
    //! Inherit std::vector< std::vector< std::vector< T > > > constructors.
    using std::vector< std::vector< std::vector< T > > >::vector;

    //! Default constructor
    TripleVector( ) { }

    //! Constructor with vector of vectors of vectors
    TripleVector( const std::vector< std::vector< std::vector< T > > > &vectorOfVectorsOfVectors ) :
        std::vector< std::vector< std::vector< T > > >( vectorOfVectorsOfVectors ) { }

    //! Constructor with list of lists of lists
    TripleVector( const std::initializer_list< std::initializer_list< std::initializer_list< T > > >
                  &listOfListsOfLists ) {
        for ( auto listOfLists : listOfListsOfLists ) {
            std::vector< std::vector< T > > vectorOfVectors;
            for ( auto list : listOfLists ) {
                vectorOfVectors.push_back( std::vector< T >( list ) );
            }
            this->push_back( vectorOfVectors );
        }
    }

    //! Constructor with list of lists of vectors
    TripleVector( const std::initializer_list< std::initializer_list< std::vector< T > > > &listOfListsOfVectors ) {
        for ( auto listOfVectors : listOfListsOfVectors ) {
            for ( auto vector : listOfVectors ) {
                this->push_back( vector );
            }
        }
    }

    //! Constructor with size m for the first level, size n for the second level, empty vectors for the thrid level
    TripleVector( const unsigned int m, const unsigned int n ) :
        std::vector< std::vector< std::vector< T > > >
        ( m, std::vector< std::vector< T > >( n, std::vector< T >() ) ) { }

    //! Constructor with size m for the first level, size n for the second level, size p for the third level
    TripleVector( const unsigned int m, const unsigned int n, const unsigned int p ) :
        std::vector< std::vector< std::vector< T > > >
        ( m, std::vector< std::vector< T > >( n, std::vector< T >( p ) ) ) { }

    /*
    //! Constructor with list of vectors of vectors
    TripleVector( const std::initializer_list< std::vector< std::vector< T > > > &listOfVectorsOfVectors ) {
        for ( auto listOfVectors : listOfVectorsOfVectors ) {
            for ( auto vector : listOfVectors ) {
                this->push_back( vector );
            }
        }
    }

    //! Constructor with size m for the first level, empty vectors of vectors for the second and third levels
    TripleVector( const unsigned int m ) :
        std::vector< std::vector< std::vector< T > > >( m, std::vector< std::vector< T > >() ) { }

    //! Constructor with repeated vector of vectors
    TripleVector( const unsigned int m, const std::vector< std::vector< T > > &repeated ) :
        std::vector< std::vector< std::vector< T > > >( m, repeated ) { }

    TripleVector & operator= ( const TripleVector & ) = default;
    */

    T& operator() ( const unsigned int i, const unsigned int j, const unsigned int k ) {
        return ( *this ).at( i ).at( j ).at( k );
    }

    const T& operator() ( const unsigned int i, const unsigned int j, const unsigned int k ) const {
        return ( *this )[ i ][ j ][ k ];
    }

};
typedef TripleVector< double > TripleVectord;
typedef TripleVector< int >    TripleVectori;


/*
class Equatable
{
public:
    virtual bool operator == ( const Equatable &other ) const;

    virtual bool operator != ( const Equatable &other ) const {
        return ! ( *this == other );
    }
};
*/


} // namespace sst

} // namespace propagators

} // namespace tudat


#endif // TUDAT_PROPAGATORS_DSST_VECTORS_H
