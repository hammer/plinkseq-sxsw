#include "matrix.h"
#include "statistics.h"

template<class T> Data::Vector<T> Data::Vector<T>::operator*( const Data::Matrix<T> & rhs ) const
{
  return Statistics::matrix_multiply( *this , rhs );
}

template<class T> Data::Vector<T> Data::Vector<T>::operator-( const Data::Vector<T> & rhs ) const
{
  Data::Vector<T> r( rhs.dim1() );
  for (int i=0; i<rhs.dim1(); i++) r[i] = (*this)[i] - rhs[i];
  return r;
}

template<class T> Data::Vector<T> Data::Vector<T>::operator+( const Data::Vector<T> & rhs ) const
{
  Data::Vector<T> r( rhs.dim1() );
  for (int i=0; i<rhs.dim1(); i++) r[i] = (*this)[i] + rhs[i];
  return r;
}

template<class T> Data::Matrix<T> Data::Matrix<T>::operator*( const Data::Matrix<T> & rhs ) const
{
  return Statistics::matrix_multiply( *this , rhs );
}

template<class T> Data::Vector<T> Data::Matrix<T>::operator*( const Data::Vector<T> & rhs ) const
{
  return Statistics::matrix_multiply( *this , rhs );
}

template<class T> Data::Matrix<T> Data::Matrix<T>::operator-( const Data::Matrix<T> & rhs ) const
{
  Data::Matrix<T> r( rhs.dim1() , rhs.dim2() );
  for (int i=0; i<rhs.dim1(); i++) 
    for (int j=0; j<rhs.dim2(); j++) 
      r(i,j) = (*this)(i,j) - rhs(i,j);
  return r;
}

template<class T> Data::Matrix<T> Data::Matrix<T>::operator+( const Data::Matrix<T> & rhs ) const
{
  Data::Matrix<T> r( rhs.dim1() , rhs.dim2() );
  for (int i=0; i<rhs.dim1(); i++) 
    for (int j=0; j<rhs.dim2(); j++) 
      r(i,j) = (*this)(i,j) + rhs(i,j);
  return r;
}


// added to avoid linker errors, given we've defined these templated functions in the .cpp file.
template class Data::Matrix<double>;
template class Data::Vector<double>;
 
