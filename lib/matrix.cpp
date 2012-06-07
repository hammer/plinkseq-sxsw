#include "plinkseq/matrix.h"
#include "plinkseq/statistics.h"
#include <sstream>

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


// pretty-printers
template<class T> std::string Data::Vector<T>::print( const std::string & label , const int nelem ) const
{
  int aelem =  nelem == 0 || nelem  > size() ? size() : nelem ;

  std::stringstream ss;
  if ( label != "" ) ss << label << "\n";
  for (int r=0;r<aelem;r++)
    {
      ss << " [" << data[r] << " ]\n";
    }
  return ss.str();
}

template<class T> std::string Data::Matrix<T>::print( const std::string & label , const int nrow , const int ncol) const
{
  int arow =  nrow == 0 || nrow > dim1() ? dim1() : nrow ; 
  int acol =  ncol == 0 || ncol > dim2() ? dim2() : ncol ;

  std::stringstream ss;
  if ( label != "" ) ss << label << "\n";

  for (int r=0;r<arow;r++)
    {
      ss << " [" ;
      for (int c=0;c<acol;c++)
	ss << " " << (*this)(r,c) ;
      ss << " ]\n";
    }
  return ss.str();
}



// added to avoid linker errors, given we've defined these templated functions in the .cpp file.
template class Data::Matrix<double>;
template class Data::Vector<double>;
 
