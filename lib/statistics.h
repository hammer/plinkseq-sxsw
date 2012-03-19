#ifndef __PSEQ_STATISTICS_H__
#define __PSEQ_STATISTICS_H__

#include "matrix.h"
#include "fisher.h"

#include <vector>

namespace Statistics { 
  
  template<class T> inline const T SQR(const T a) {return a*a;}
  template<class T> inline const T FNMAX(const T &a, const T &b) {return b > a ? (b) : (a);}
  template<class T> inline const T FNMIN(const T &a, const T &b) {return b < a ? (b) : (a);}
  template<class T> inline const T SIGN(const T &a, const T &b) {return b >= 0 ? (a >= 0 ? a : -a) : (a >= 0 ? -a : a);}
  template<class T> inline void SWAP(T &a, T &b) {T dum=a; a=b; b=dum;}

  std::vector<double> as_vector( const Data::Vector<double> & );
  
  // Singular Value Decomposition
  bool svdcmp( Data::Matrix<double> & , Data::Vector<double> & , Data::Matrix<double> & );
  void svbksb( Data::Matrix<double> & , Data::Vector<double> & , Data::Matrix<double> & , Data::Vector<double> & , Data::Vector<double> & );


  long unsigned int factorial(int n);
  long unsigned int combin(int n, int k);
  
  double dbinom( const int k , const int n , double p );
  double chi2_prob( double x, double df );
  double t_prob( double x, double df );
  double ltqnorm( double p );
  
  // mainly helper functions for the above
  double factln(int);
  double factrl(int);
  double gamln(double);
  double gammln(double);
  double factln(int n);
  double bico(int n, int k);
  double dbinom_raw( const double , const double , const double);
  double pythag( const double , const double );

  // matrix operations

  Data::Matrix<double> matrix_sqrt( const Data::Matrix<double> & );
  Data::Matrix<double> transpose( const Data::Matrix<double> & );
  Data::Matrix<double> inverse( const Data::Matrix<double> & , bool * flag = NULL );

  Data::Matrix<double> matrix_multiply( const Data::Matrix<double> & , const Data::Matrix<double> & );
  Data::Vector<double> matrix_multiply( const Data::Matrix<double> & , const Data::Vector<double> & );  
  Data::Vector<double> matrix_multiply( const Data::Vector<double> & , const Data::Matrix<double> & );
  double matrix_inner_product( const Data::Vector<double> & , const Data::Vector<double> & );
  Data::Matrix<double> matrix_outer_product( const Data::Vector<double> & , const Data::Vector<double> & );  

  // Mean, variance and covariance
  
  Data::Vector<double> mean( const Data::Matrix<double> & );
  
  Data::Vector<double> variance( const Data::Matrix<double> & );
  Data::Vector<double> variance( const Data::Matrix<double> & , const Data::Vector<double> & );
  
  Data::Matrix<double> covariance_matrix( const Data::Matrix<double> & );
  Data::Matrix<double> covariance_matrix( const Data::Matrix<double> & , const Data::Vector<double> & );

  Data::Matrix<double> covariance_matrix( const Data::Matrix<double> & , const Data::Matrix<double> & );
  Data::Matrix<double> covariance_matrix( const Data::Matrix<double> & , const Data::Vector<double> & , 
					    const Data::Matrix<double> & , const Data::Vector<double> & );
  
  std::vector<double> canonical_correlation( const Data::Matrix<double> & , const Data::Matrix<double> & , double * pv = NULL );     
  
  double bartlett(const int N, 
		  const int p, 
		  const int q, 
		  const std::vector<double> & eigen );
  

  struct Eigen
  {
    Eigen(const int n) : d(n), z(n,n) { }     
    Data::Vector<double> d; // eigenvalues
    Data::Matrix<double> z; // eigenvectors
  };
  
  Statistics::Eigen eigenvectors( Data::Matrix<double> & );
  
  void EV_tqli( Data::Vector<double> & d , 
		Data::Vector<double> & e , 
		Data::Matrix<double> & z );

  void EV_tred2( Data::Matrix<double> & a , 
		 Data::Vector<double> & d , 
		 Data::Vector<double> & e );
  

  Data::Vector<double> eigenvalues( Data::Matrix<double> & );  
  void tqli( Data::Vector<double> & d , 
	     Data::Vector<double> & e );
	     
  
  void tred2( Data::Matrix<double> & a , 
	      Data::Vector<double> & d , 
	      Data::Vector<double> & e );
  

}

  

#endif
