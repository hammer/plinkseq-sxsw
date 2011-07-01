#include "statistics.h"
#include "helper.h"
#include "matrix.h"
#include "dcdflib.h"
#include "ipmpar.h"

#include <iostream>
#include <cmath>
#include <algorithm>

extern Log plog;


#ifndef M_2PI
#define M_2PI 6.283185307179586476925286766559/* 2*pi */
#endif

#ifndef M_LN_SQRT_2PI
#define M_LN_SQRT_2PI 0.918938533204672741780329736406/* log(sqrt(2*pi)) */
#endif

#ifndef M_LN_SQRT_PId2
#define M_LN_SQRT_PId2 0.225791352644727432363097614947/* log(sqrt(pi/2)) */
#endif



Data::Vector<double> Statistics::mean( const Data::Matrix<double> & d )
{
  Data::Vector<double> m( d.dim2() );

  for (int j=0; j<d.dim2(); j++)
    {
      for (int i=0; i<d.dim1(); i++)
	m[j] += d(i,j);
      m[j] /= d.dim1();
    }
  return m;
}

Data::Vector<double> Statistics::variance( const Data::Matrix<double> & d )
{
  return variance( d , mean(d) );
}

Data::Vector<double> Statistics::variance( const Data::Matrix<double> & d , const Data::Vector<double> & u )
{
  Data::Vector<double> v( d.dim2() ); 
  Data::Matrix<double> s = covariance_matrix( d , u , d , u );
  for (int i=0; i<d.dim2(); i++) v(i) = s(i,i);
  return v;
}

Data::Matrix<double> Statistics::covariance_matrix( const Data::Matrix<double> & d )
{
  return covariance_matrix( d , mean(d) );
}

Data::Matrix<double> Statistics::covariance_matrix( const Data::Matrix<double> & d , const Data::Vector<double> & u )
{
  return covariance_matrix( d , u , d , u );  
}  



std::vector<double> Statistics::as_vector( const Data::Vector<double> & d )
{
  std::vector<double> v( d.size() ) ;
  for (int i=0; i<d.size(); i++) v[i] = d[i];
  return v;
}

Data::Matrix<double> Statistics::covariance_matrix( const Data::Matrix<double> & x , 
						      const Data::Vector<double> & u ,
						      const Data::Matrix<double> & y , 
						      const Data::Vector<double> & v )
{
  
  // calculate Sxy e.g. lower quadrant of partitioned covariance matrix
  // note -- ignores symmetry -- easy speedup to add 
  if ( x.dim1() != y.dim1() ) Helper::halt("internal error, unequal row numbers in covariance_matrix()"); 
  const int n = x.dim1();
  Data::Matrix<double> s( x.dim2() , y.dim2() ) ;
  for (int i=0; i<x.dim2(); i++)
    for (int j=0; j<y.dim2(); j++)
      {
	const double mx = u[i];
	const double my = v[j];
	for (int k=0;k<n;k++)
	  s(i,j) += ( x(k,i) - mx ) * ( y(k,j) - my ); 
	s(i,j) /= n-1;
	//	s(j,i) = s(i,j);
      }
  return s;
}

Data::Matrix<double> Statistics::covariance_matrix( const Data::Matrix<double> & x , const Data::Matrix<double> & y )
{
  return covariance_matrix( x , mean(x) , y , mean(y) );
}

Data::Matrix<double> Statistics::transpose( const Data::Matrix<double> & d )
{
  const int row = d.dim1();
  const int col = d.dim2();
  Data::Matrix<double> r( col, row );
  for (int i = 0; i < row; i++)
    for (int j = 0; j < col; j++)
      r(j,i) = d(i,j);
  return r;
}

Data::Matrix<double> Statistics::inverse( const Data::Matrix<double> & u_orig, bool * flag )
{
  
  const double eps = 1e-24; 

  Data::Matrix<double> u = u_orig;
  
  if ( u.dim1() == 0 || u.dim1() != u.dim2() ) 
    Helper::halt("cannot inverted non-square matrix");
  int n = u.dim1();
  
  Data::Vector<double> w(n);  
  Data::Matrix<double> v(n,n);

  if ( flag ) 
    *flag = Statistics::svdcmp( u, w, v ); 
  else 
    Statistics::svdcmp( u, w, v ); 
  
  // Look for singular values

  double wmax = 0;
  for (int i=0; i<n; i++)
    wmax = w[i] > wmax ? w[i] : wmax;
  double wmin = wmax * eps;
  for (int i=0; i<n; i++)
    w[i] = w[i] < wmin ? 0 : 1/w[i];
  
  // u w t(v)
  // row U * 1/w
  
  // results matrix
  Data::Matrix<double> r(n,n);

  for (int i=0; i<n; i++)
    for (int j=0; j<n; j++)
      u(i,j) *= w[j];
  
  // [nxn].[t(v)] 
  for (int i=0; i<n; i++)
    for (int j=0; j<n; j++)
      for (int k=0; k<n; k++)
	r(i,j) += u(i,k) * v(j,k);
    
  return r;
}



Data::Matrix<double> Statistics::matrix_sqrt( const Data::Matrix<double> & u_orig )
{
  
  Data::Matrix<double> u = u_orig;

  // Using SVD, square root is U . sqrt(D) . V_T
  //  msqrt <- function(m) { m <- svd(m); m$u %*% sqrt(diag(m$d)) %*% t(m$v) }  

  const double eps = 1e-12;
  
  int n = u.dim1();

  Data::Vector<double> d(n);
  Data::Matrix<double> v(n,n);

  Statistics::svdcmp(u,d,v);

  // Take square root of diagonal values                                                                                                                                                                                                                                                                          
  for (int i=0; i<n; i++)
    d[i] = sqrt(d[i]);
  
  
  // Multiplication to reconstruct original                                                                                                                                                                                                                                                                             

  Data::Matrix<double> r(n,n);
  Data::Matrix<double> r2(n,n);

  for (int i=0; i<n; i++)
    for (int j=0; j<n; j++)
      r(i,j) = u(i,j) * d[j];
  
  for (int i=0; i<n; i++)
    for (int j=0; j<n; j++)
      for (int k=0; k<n; k++)
 	r2(i,j) += r(i,k) * v(j,k);
  
  return r2;
  
}

bool Statistics::svdcmp( Data::Matrix<double> & a, Data::Vector<double> & w , Data::Matrix<double> & v )
{
  
  bool flag;
  
  int i,its,j,jj,k,l,nm;
  
  double anorm,c,f,g,h,s,scale,x,y,z;
  
  double volatile temp;
  
  int m=a.dim1();
  if ( m == 0 ) Helper::halt("Internal problem in SVD function (no observations left?)");

  int n = a.dim2();

  std::vector<double> rv1( n );
  
  g=scale=anorm=0.0;
  for (i=0;i<n;i++) {
    l=i+2;
    rv1[i]=scale*g;
    g=s=scale=0.0;
    if (i < m) {
      for (k=i;k<m;k++) scale += fabs(a(k,i));
      if (scale != 0.0) {
	for (k=i;k<m;k++) {
	  a[k][i] /= scale;
	  s += a(k,i)*a(k,i);
	}
	f=a(i,i);
	g = -Statistics::SIGN(sqrt(s),f);
	h=f*g-s;
	a(i,i)=f-g;
	for (j=l-1;j<n;j++) {
	  for (s=0.0,k=i;k<m;k++) s += a(k,i)*a(k,j);
	  f=s/h;
	  for (k=i;k<m;k++) a(k,j) += f*a(k,i);
	}
	for (k=i;k<m;k++) a[k][i] *= scale;
      }
    }
    w[i]=scale *g;
    g=s=scale=0.0;
    if (i+1 <= m && i+1 != n) {
      for (k=l-1;k<n;k++) scale += fabs(a(i,k));
      if (scale != 0.0) {
	for (k=l-1;k<n;k++) {
	  a(i,k) /= scale;
	  s += a(i,k)*a(i,k);
	}
	f=a(i,l-1);
	g = -Statistics::SIGN(sqrt(s),f);
	h=f*g-s;
	a(i,l-1)=f-g;
	for (k=l-1;k<n;k++) rv1[k]=a(i,k)/h;
	for (j=l-1;j<m;j++) {
	  for (s=0.0,k=l-1;k<n;k++) s += a(j,k)*a(i,k);
	  for (k=l-1;k<n;k++) a(j,k) += s*rv1[k];
	}
	for (k=l-1;k<n;k++) a(i,k) *= scale;
      }
    }
    anorm=FNMAX(anorm,(fabs(w[i])+fabs(rv1[i])));
  }
  for (i=n-1;i>=0;i--) {
    if (i < n-1) {
      if (g != 0.0) {
	for (j=l;j<n;j++)
	  v[j][i]=(a[i][j]/a[i][l])/g;
	for (j=l;j<n;j++) {
	  for (s=0.0,k=l;k<n;k++) s += a[i][k]*v[k][j];
	  for (k=l;k<n;k++) v[k][j] += s*v[k][i];
	}
      }
      for (j=l;j<n;j++) v[i][j]=v[j][i]=0.0;
    }
    v[i][i]=1.0;
    g=rv1[i];
    l=i;
  }
  for (i=FNMIN(m,n)-1;i>=0;i--) {
    l=i+1;
    g=w[i];
    for (j=l;j<n;j++) a[i][j]=0.0;
    if (g != 0.0) {
      g=1.0/g;
      for (j=l;j<n;j++) {
	for (s=0.0,k=l;k<m;k++) s += a[k][i]*a[k][j];
	f=(s/a[i][i])*g;
	for (k=i;k<m;k++) a[k][j] += f*a[k][i];
      }
      for (j=i;j<m;j++) a[j][i] *= g;
    } else for (j=i;j<m;j++) a[j][i]=0.0;
    ++a[i][i];
  }
  for (k=n-1;k>=0;k--) {
    for (its=0;its<30;its++) {
      flag=true;
      for (l=k;l>=0;l--) {
	nm=l-1;
	temp=fabs(rv1[l])+anorm;
	if (temp == anorm) {
	  flag=false;
	  break;
	}
	temp=fabs(w[nm])+anorm;
	if (temp == anorm) break;
      }
      if (flag) {
	c=0.0;
	s=1.0;
	for (i=l;i<k+1;i++) {
	  f=s*rv1[i];
	  rv1[i]=c*rv1[i];
	  temp = fabs(f)+anorm;
	  if (temp == anorm) break;
	  g=w[i];
	  h=pythag(f,g);
	  w[i]=h;
	  h=1.0/h;
	  c=g*h;
	  s = -f*h;
	  for (j=0;j<m;j++) {
	    y=a[j][nm];
	    z=a[j][i];
	    a[j][nm]=y*c+z*s;
	    a[j][i]=z*c-y*s;
	  }
	}
      }
      z=w[k];
      if (l == k) {
	if (z < 0.0) {
	  w[k] = -z;
	  for (j=0;j<n;j++) v[j][k] = -v[j][k];
	}
	break;
      }

      if (its == 29) 
	{
	  plog.warn("cannot converge SVD, perhaps due to multi-colinearity"); 
	  return false;
	}

      x=w[l];
      nm=k-1;
      y=w[nm];
      g=rv1[nm];
      h=rv1[k];
      f=((y-z)*(y+z)+(g-h)*(g+h))/(2.0*h*y);
      g=Statistics::pythag(f,1.0);
      f=((x-z)*(x+z)+h*((y/(f+Statistics::SIGN(g,f)))-h))/x;
      c=s=1.0;
      for (j=l;j<=nm;j++) {
	i=j+1;
	g=rv1[i];
	y=w[i];
	h=s*g;
	g=c*g;
	z=pythag(f,h);
	rv1[j]=z;
	c=f/z;
	s=h/z;
	f=x*c+g*s;
	g=g*c-x*s;
	h=y*s;
	y *= c;
	for (jj=0;jj<n;jj++) {
	  x=v[jj][j];
	  z=v[jj][i];
	  v[jj][j]=x*c+z*s;
	  v[jj][i]=z*c-x*s;
	}
	z=pythag(f,h);
	w[j]=z;
	if (z) {
	  z=1.0/z;
	  c=f*z;
	  s=h*z;
	}
	f=c*g+s*y;
	x=c*y-s*g;
	for (jj=0;jj<m;jj++) {
	  y=a[jj][j];
	  z=a[jj][i];
	  a[jj][j]=y*c+z*s;
	  a[jj][i]=z*c-y*s;
	}
      }
      rv1[l]=0.0;
      rv1[k]=f;
      w[k]=x;
    }
  }
  return true;
}

void Statistics::svbksb( Data::Matrix<double> & u , Data::Vector<double> & w , Data::Matrix<double> & v, Data::Vector<double> & b, Data::Vector<double> & x)
{
  
  int m = u.dim1();
  int n = u.dim2();

  Data::Vector<double> tmp(n);
  
  for (int j=0; j<n; j++) 
    {
      double s = 0.0;
      if ( w[j] != 0.0) 
	{
	  for (int i=0; i<m; i++) s += u(i,j)*b[i];
	  s /= w[j];
	}
      tmp[j] = s;
    }
  
  for (int j=0; j<n; j++) 
    {
      double s = 0.0;
      for (int jj=0; jj<n; jj++) 
	s += v(j,jj) * tmp[jj]; 
      x[j]=s;
    }
}


double Statistics::pythag(const double a, const double b)
{
  double absa,absb;
 
  absa=fabs(a);
  absb=fabs(b);
  if (absa > absb) return absa*sqrt(1.0+SQR(absb/absa));
  else return (absb == 0.0 ? 0.0 : absb*sqrt(1.0+SQR(absa/absb)));
}

inline double SQR(double a)
{
  return a*a;
}





std::vector<double> Statistics::canonical_correlation( const Data::Matrix<double> & x , const Data::Matrix<double> & y , double * pv )
{
  
  // 1. Partitioned covariance matrix  
  //    S_XX  S_YX
  //    S_XY  S_YY
              
  const int nx = x.dim2();
  const int ny = y.dim2();
  const int ne = nx < ny ? nx : ny ;             
  
  // total covariance matrix
  
  if ( x.dim1() != y.dim1() ) Helper::halt("different number of individuals on left and right hand of canonical correlation");
  int nind = x.dim1(); // total N individuals

  Data::Matrix<double> I11 = covariance_matrix( x, x );
  Data::Matrix<double> I12 = covariance_matrix( x, y );
  Data::Matrix<double> I21 = covariance_matrix( y, x );
  Data::Matrix<double> I22 = covariance_matrix( y, y );
  
  Data::Matrix<double> I11b( nx,nx );
  Data::Matrix<double> I22b( ny,ny );

  // 2. Calculate the p x p matrix M1 = inv(sqrt(sig11)) %*% sig12 %*% inv(sig22) %*% sig21 %*% inv(sqrt(sig11))

  bool flag = true;
  I11 = Statistics::matrix_sqrt( I11 );
  I11 = Statistics::inverse(I11,&flag);
  if ( ! flag ) plog.warn( "could not invert matrix in canonical_correlation()" );
  I22 = Statistics::inverse(I22,&flag);
  if ( ! flag ) plog.warn( "could not invert matrix in canonical_correlation()" );
  I22b = Statistics::matrix_sqrt( I22b ); // For Step 4b
  I22b = Statistics::inverse( I22b , &flag );
  if ( ! flag ) plog.warn( "could not invert matrix in canonical_correlation()" );
  I11b = Statistics::inverse( I11b, &flag );
  if ( ! flag ) plog.warn( "could not invert matrix in canonical_correlation()" );
  
  
  Data::Matrix<double> M1 = Statistics::matrix_multiply( 
                              Statistics::matrix_multiply( 
                              Statistics::matrix_multiply( 
                              Statistics::matrix_multiply( I11 , I12 ) , I22 ) , I21 ) , I11 );

  // Compute and sort eigen values only 

  std::vector<double> sorted_eigenvalues = Statistics::as_vector( Statistics::eigenvalues( M1 ) );
  std::sort( sorted_eigenvalues.begin(), sorted_eigenvalues.end(), std::greater<double>() );

  // Display largest canonical correlation and its position

  // Use Bartlett's test to get p-values for each canonical correlation

  if ( pv ) *pv = Statistics::bartlett( nind, nx , ny , sorted_eigenvalues ); 
  
  return sorted_eigenvalues;
  
}


double Statistics::bartlett(const int N, 
			    const int p, 
			    const int q, 
			    const std::vector<double> & eigen )
{
  int p2 = p < q ? p : q; // Number of canonical correlations  
  double prod_eigen=1.0;
  for (int j=0; j<p2; j++)
    prod_eigen *= (1-eigen[j]);      
  double chisq = -1*(N - 1 - 0.5*(p+q+1)) * log(prod_eigen);      
  return chi2_prob( chisq, p*q );
}


double Statistics::chi2_prob(double x, double df)
{

  if ( ! Helper::realnum(x) ) return -9;

  double p, q;
  int st = 0;      // error variable
  int w = 1;      // function variable
  double bnd = 1; // boundary function

  // NCP is set to 0
  cdfchi(&w,&p,&q,&x,&df,&st,&bnd);

  // Check status
  if (st != 0 ) return -9;

  // Return p-value
  return q;
  
}

// Inverse normal distribution

/*
 * Lower tail quantile for standard normal distribution function.
 *
 * This function returns an approximation of the inverse cumulative
 * standard normal distribution function.  I.e., given P, it returns
 * an approximation to the X satisfying P = Pr{Z <= X} where Z is a
 * random variable from the standard normal distribution.
 *
 * The algorithm uses a minimax approximation by rational functions
 * and the result has a relative error whose absolute value is less
 * than 1.15e-9.
 *
 * Author:      Peter J. Acklam
 * Time-stamp:  2002-06-09 18:45:44 +0200
 * E-mail:      jacklam@math.uio.no
 * WWW URL:     http://www.math.uio.no/~jacklam
 *
 * C implementation adapted from Peter's Perl version
 */


/* Coefficients in rational approximations. */

static const double a[] =
  {
    -3.969683028665376e+01,
    2.209460984245205e+02,
    -2.759285104469687e+02,
    1.383577518672690e+02,
    -3.066479806614716e+01,
     2.506628277459239e+00
  };

static const double b[] =
  {
    -5.447609879822406e+01,
    1.615858368580409e+02,
    -1.556989798598866e+02,
    6.680131188771972e+01,
    -1.328068155288572e+01
  };

static const double c[] =
  {
    -7.784894002430293e-03,
    -3.223964580411365e-01,
    -2.400758277161838e+00,
    -2.549732539343734e+00,
    4.374664141464968e+00,
     2.938163982698783e+00
  };

static const double d[] =
  {
    7.784695709041462e-03,
    3.224671290700398e-01,
    2.445134137142996e+00,
    3.754408661907416e+00
  };

#define LOW 0.02425
#define HIGH 0.97575

double Statistics::ltqnorm( double p )
{
  
  double q, r;
  
  if (p < 0 || p > 1)
    {
      return 0.0;
    }
  else if (p == 0)
    {
      return -HUGE_VAL /* minus "infinity" */;
    }
  else if (p == 1)
    {
      return HUGE_VAL /* "infinity" */;
    }
  else if (p < LOW)
    {
      /* Rational approximation for lower region */
      q = sqrt(-2*log(p));
      return (((((c[0]*q+c[1])*q+c[2])*q+c[3])*q+c[4])*q+c[5]) /
        ((((d[0]*q+d[1])*q+d[2])*q+d[3])*q+1);
    }
  else if (p > HIGH)
    {
      /* Rational approximation for upper region */
      q  = sqrt(-2*log(1-p));
      return -(((((c[0]*q+c[1])*q+c[2])*q+c[3])*q+c[4])*q+c[5]) /
        ((((d[0]*q+d[1])*q+d[2])*q+d[3])*q+1);
    }
  else
    {
      /* Rational approximation for central region */
      q = p - 0.5;
      r = q*q;
      return (((((a[0]*r+a[1])*r+a[2])*r+a[3])*r+a[4])*r+a[5])*q /
        (((((b[0]*r+b[1])*r+b[2])*r+b[3])*r+b[4])*r+1);
    }
}


double Statistics::t_prob(double T, double df)
{

  if ( ! Helper::realnum(T) ) return -9; 

  T = fabs(T);
  
  double p, q;
  int st = 0;      // error variable
  int w = 1;       // function variable
  double bnd = 1;  // boundary function
  
  // NCP is set to 0
  cdft(&w,&p,&q,&T,&df,&st,&bnd);
  
  // Check status
  if (st != 0 ) return -9;
  
  // Return two-sided p-value
  return 2*q;
  
}

  
//   // Sort evectors. Rows must be ordered according to cancor value (highest first)
  
//   Data::Matrix<double> sorted_eigenvectors = eigen.z;
  
//   std::vector<int> order_eigenvalues(nx);

//   for (int i=0; i<nx; i++)
//     {
 
//      // Determine position of the vector associated with the ith cancor

//       for (int j=0; j<n1; j++)
// 	{

// 	  if ( eigen.d[j] == sorted_eigenvalues_gene1[i] )               
// 	    {
// 	      if (i==0)        
// 		{
// 		  order_eigenvalues_gene1[i]=j;
// 		  break;
// 		}
// 	      else
// 		{
// 		  if (j!=order_eigenvalues_gene1[i-1])
//                     {
//                       order_eigenvalues_gene1[i]=j;
//                       break;
//                     }
// 		}
// 	    }
// 	}
//     }

//   for (int i=0; i<n1; i++)
//     {
//       sorted_eigenvectors_gene1[i] = eigen.z[order_eigenvalues_gene1[i]];
//     }

//   //   cout << "Eigenvector matrix - unsorted:\n";
//   // display(gene1_eigen.z);
//   //cout << "Eigenvector matrix - sorted:\n";
//   //display(sorted_eigenvectors_gene1);


//   ////////////////////////////////////////////////////////
//   // Step 4b. Calculate the q x q eigenvectors of M2 (f). These are
//   // required to compute the coefficients used to build the p
//   // canonical variates b[k] for gene2 (see below). The first p are
//   // given by: f[k] = (1/sqrt(eigen[k])) * inv_sqrt_I22 %*% I21 %*%
//   // inv_sqrt_sig11 %*% e[k] for (k in 1:p) { e.vectors.gene2[,k] =
//   // (1/sqrt(e.values[k])) * inv.sqrt.sig22 %*% sig21 %*%
//   // inv.sqrt.sig11 %*% e.vectors.gene1[,k] }
           
//   matrix_t M2;

//   multMatrix(I22b, I21, tmp);
//   multMatrix(tmp, I11b, M2);
//   multMatrix(M2, I12, tmp);
//   multMatrix(tmp, I22b, M2);
//   Eigen gene2_eigen = eigenvectors(M2);

//   //cout << "Eigenvalues Gene 2 - unsorted:\n";

//   //display(gene2_eigen.d);
 
//   // Sort evalues for gene2
//   vector<double> sorted_eigenvalues_gene2 = gene2_eigen.d;
//   sort(sorted_eigenvalues_gene2.begin(),sorted_eigenvalues_gene2.end(),greater<double>());

//   // Sort eigenvectors for gene2
//   matrix_t sorted_eigenvectors_gene2 = gene2_eigen.z;
//   vector<int> order_eigenvalues_gene2(ny);

//   for (int i=0; i<ny; i++)
//     {
//       // Determine position of the vector associated with the ith cancor
//       for (int j=0; j<n2; j++)
// 	{
// 	  if (gene2_eigen.d[j]==sorted_eigenvalues_gene2[i])
// 	    {
// 	      if (i==0)
// 		{
// 		  order_eigenvalues_gene2[i]=j;
// 		  break;
// 		}
// 	      else
// 		{
// 		  if (j!=order_eigenvalues_gene2[i-1])
// 		    {
// 		      order_eigenvalues_gene2[i]=j;
// 		      break;
// 		    }
// 		}
// 	    }
// 	}
//     }

//   for (int i=0; i<n2; i++)
//     {
//       sorted_eigenvectors_gene2[i] = gene2_eigen.z[order_eigenvalues_gene2[i]];
//     }

//   //cout << "Eigenvector matrix Gene 2 - unsorted:\n";
//   //display(gene2_eigen.z);

//   //cout << "Eigenvector matrix Gene 2 - sorted:\n";
//   //display(sorted_eigenvectors_gene2);

//   //exit(0);


//   //////////////////////////////////////////////////////////////////////////////////
//   // Step 5 - Calculate the gene1 (pxp) and gene2 (pxq) coefficients
//   // used to create the canonical variates associated with the p
//   // canonical correlations

//   transposeMatrix(gene1_eigen.z);
//   transposeMatrix(gene2_eigen.z);

//   matrix_t coeff_gene1;
//   matrix_t coeff_gene2;

//   multMatrix(gene1_eigen.z, I11, coeff_gene1);
//   multMatrix(gene2_eigen.z, I22b, coeff_gene2);

//   //cout << "Coefficients for Gene 1:\n";
//   //display(coeff_gene1);

//   //cout << "Coefficients for Gene 2:\n";
//   //display(coeff_gene2);

//   //exit(0);

//   ///////////////////////////////////////////////////////////////////////
//   // Step 6 - Compute the gene1 and gene2 canonical variates
//   // associated with the highest canonical correlation NOTE: the
//   // original variables of data need to have the mean subtracted  first!
//   // Otherwise, the resulting correlation between variate.gene1 and
//   // variate.gene1 != estimated cancor.

//   // For each individual, eg compos.gene1 =
//   // evector.gene1[1]*SNP1.gene1 + evector.gene1[2]*SNP2.gene1 + ...

//   /////////////////////////////////
//   // Consider each SNP in gene1

//   vector<double> gene1(nind);
            
//   for (int j=0; j<n1; j++)
//     {

//       CSNP * ps = pSNP[j];


//       ///////////////////////////
//       // Iterate over individuals

//       for (int i=0; i< P.n ; i++)
// 	{

// 	  // Only need to look at one perm set
// 	  bool a1 = ps->one[i];
// 	  bool a2 = ps->two[i];

// 	  if ( a1 )
// 	    {
// 	      if ( a2 ) // 11 homozygote
// 		{
// 		  gene1[i] += (1 - mean[j]) * coeff_gene1[order_eigenvalues_gene1[cancor1_pos]][j];
// 		}
// 	      else      // 12 
// 		{
// 		  gene1[i] += (0 - mean[j]) * coeff_gene1[order_eigenvalues_gene1[cancor1_pos]][j];
// 		}
// 	    }
// 	  else
// 	    {
// 	      if ( a2 )      // 21
// 		{
// 		  gene1[i] += (0 - mean[j]) * coeff_gene1[order_eigenvalues_gene1[cancor1_pos]][j];
// 		}
// 	      else           // 22 homozygote
// 		{
// 		  gene1[i] += (-1 - mean[j]) * coeff_gene1[order_eigenvalues_gene1[cancor1_pos]][j];
// 		}
// 	    }

// 	} // Next individual

//     } // Next SNP in gene1

//   /////////////////////////////////
//   // Consider each SNP in gene2
//   vector<double> gene2(P.n);
//   int cur_snp = -1;            
//   for (int j=n1; j<n1+n2; j++)
//     {

//       cur_snp++;
//       CSNP * ps = pSNP[j];


                
//       // Iterate over individuals

//       for (int i=0; i<P.n; i++)
// 	{
                    
// 	  // Only need to look at one perm set
// 	  bool a1 = ps->one[i];
// 	  bool a2 = ps->two[i];

// 	  if ( a1 )
// 	    {
// 	      if ( a2 ) // 11 homozygote
// 		{
// 		  gene2[i] += (1 - mean[j]) * coeff_gene2[order_eigenvalues_gene2[cancor1_pos]][cur_snp];
// 		}
// 	      else      // 12
// 		{
// 		  gene2[i] += (0 - mean[j]) * coeff_gene2[order_eigenvalues_gene2[cancor1_pos]][cur_snp];
// 		}
// 	    }
// 	  else
// 	    {
// 	      if ( a2 )      // 21
// 		{
// 		  gene2[i] += (0 - mean[j]) * coeff_gene2[order_eigenvalues_gene2[cancor1_pos]][cur_snp];
// 		}
// 	      else           // 22 homozygote
// 		{
// 		  gene2[i] += (-1 - mean[j]) * coeff_gene2[order_eigenvalues_gene2[cancor1_pos]][cur_snp];
// 		}
// 	    }
                    
// 	} // Next individual
                
//     } // Next SNP in gene2


//   // Store gene1.variate and gene2.variate in the multiple_covariates field of P.sample
//   // TO DO: NEED TO CHECK IF FIELDS ARE EMPTY FIRST!

//   for (int i=0; i<P.n; i++)
//     {
//       P.sample[i]->clist.resize(2);
//       P.sample[i]->clist[0] = gene1[i];
//       P.sample[i]->clist[1] = gene2[i];
//     }



  
//  return cc;  
//}


Data::Vector<double> Statistics::eigenvalues( Data::Matrix<double> & a )
{
  // 'a' should be a square, symmetric matrix
  int n=a.dim1();
  Data::Vector<double> e(n);
  Data::Vector<double> d(n);
  tred2(a,d,e);
  tqli(d,e);
  return d;
}

// Householder method to reduce real, symmetric matrix
// to tridiagonal form
// Modified to return only eigenvalues.
void Statistics::tred2( Data::Matrix<double> & a , 
			Data::Vector<double> & d,
			Data::Vector<double> & e)
{
  int l,k,j,i;
  double scale,hh,h,g,f;

  int n=d.dim1();
  for (i=n-1;i>0;i--) {
    l=i-1;
    h=scale=0.0;
    if (l > 0) {
      for (k=0;k<l+1;k++)
	scale += fabs(a[i][k]);
      if (scale == 0.0)
	e[i]=a[i][l];
      else {
	for (k=0;k<l+1;k++) {
	  a[i][k] /= scale;
	  h += a[i][k]*a[i][k];
	}
	f=a[i][l];
	g=(f >= 0.0 ? -sqrt(h) : sqrt(h));
	e[i]=scale*g;
	h -= f*g;
	a[i][l]=f-g;
	f=0.0;
	for (j=0;j<l+1;j++) {
	  // Next statement can be omitted if eigenvectors not wanted
// 	  a[j][i]=a[i][j]/h;
	  g=0.0;
	  for (k=0;k<j+1;k++)
	    g += a[j][k]*a[i][k];
	  for (k=j+1;k<l+1;k++)
	    g += a[k][j]*a[i][k];
	  e[j]=g/h;
	  f += e[j]*a[i][j];
	}
	hh=f/(h+h);
	for (j=0;j<l+1;j++) {
	  f=a[i][j];
	  e[j]=g=e[j]-hh*f;
	  for (k=0;k<j+1;k++)
	    a[j][k] -= (f*e[k]+g*a[i][k]);
	}
      }
    } else
      e[i]=a[i][l];
    d[i]=h;
  }
  // Next statement can be omitted if eigenvectors not wanted
//   d[0]=0.0;
  e[0]=0.0;
  // Contents of this loop can be omitted if eigenvectors not
  //	wanted except for statement d[i]=a[i][i];
  for (i=0;i<n;i++) {
//     l=i;
//     if (d[i] != 0.0) {
//       for (j=0;j<l;j++) {
// 	g=0.0;
// 	for (k=0;k<l;k++)
// 	  g += a[i][k]*a[k][j];
// 	for (k=0;k<l;k++)
// 	  a[k][j] -= g*a[k][i];
//       }
//     }
    d[i]=a[i][i];
//     a[i][i]=1.0;
//     for (j=0;j<l;j++) a[j][i]=a[i][j]=0.0;
  }
}

// Modified to return only eigenvalues.
void Statistics::tqli( Data::Vector<double> & d, Data::Vector<double> & e )
{
  int m,l,iter,i,k;
  double s,r,p,g,f,dd,c,b;
  double volatile temp;
  int n=d.dim1();
  for (i=1;i<n;i++) e[i-1]=e[i];
  e[n-1]=0.0;
  for (l=0;l<n;l++) {
    iter=0;
    do {
      for (m=l;m<n-1;m++) {
	dd=fabs(d[m])+fabs(d[m+1]);
	temp=fabs(e[m])+dd;
	if (temp == dd) break;
      }
      if (m != l) {
	if (iter++ == 30) Helper::halt("Internal problem in tqli routine");
	g=(d[l+1]-d[l])/(2.0*e[l]);
	r=pythag(g,1.0);
	g=d[m]-d[l]+e[l]/(g+SIGN(r,g));
	s=c=1.0;
	p=0.0;
	for (i=m-1;i>=l;i--) {
	  f=s*e[i];
	  b=c*e[i];
	  e[i+1]=(r=pythag(f,g));
	  if (r == 0.0) {
	    d[i+1] -= p;
	    e[m]=0.0;
	    break;
	  }
	  s=f/r;
	  c=g/r;
	  g=d[i+1]-p;
	  r=(d[i]-g)*s+2.0*c*b;
	  d[i+1]=g+(p=s*r);
	  g=c*r-b;
	  // Next loop can be omitted if eigenvectors not wanted
	  /* for (k=0;k<n;k++) {
	     f=z[k][i+1];
	     z[k][i+1]=s*z[k][i]+c*f;
	     z[k][i]=c*z[k][i]-s*f;
	     } */
	}
	if (r == 0.0 && i >= l) continue;
	d[l] -= p;
	e[l]=g;
	e[m]=0.0;
      }
    } while (m != l);
  }
}



Statistics::Eigen Statistics::eigenvectors( Data::Matrix<double> & a )
{
  // 'a' should be a square, symmetric matrix
  int n=a.dim1();
  
  Statistics::Eigen E(n);
  
  Data::Vector<double> e(n);
  Statistics::EV_tred2( a, E.d, e);
  Statistics::EV_tqli( E.d, e, a);
  E.z = a;
  return E;
}


// Householder method to reduce real, symmetric matrix
// to tridiagonal form
// Modified to return both eigenvalues and eigenvectors

void Statistics::EV_tred2( Data::Matrix<double> & a ,
			   Data::Vector<double> & d , 
			   Data::Vector<double> & e )
{
  int l,k,j,i;
 
  double scale,hh,h,g,f;
  
  int n=d.dim1();
  for (i=n-1;i>0;i--) {
    l=i-1;
    h=scale=0.0;
    if (l > 0) {
      for (k=0;k<l+1;k++)
        scale += fabs(a[i][k]);
      if (scale == 0.0)
        e[i]=a[i][l];
      else {
        for (k=0;k<l+1;k++) {
          a[i][k] /= scale;
          h += a[i][k]*a[i][k];
        }
        f=a[i][l];
        g=(f >= 0.0 ? -sqrt(h) : sqrt(h));
        e[i]=scale*g;
        h -= f*g;
        a[i][l]=f-g;
        f=0.0;
        for (j=0;j<l+1;j++) {
          a[j][i]=a[i][j]/h;
          g=0.0;
          for (k=0;k<j+1;k++)
            g += a[j][k]*a[i][k];
          for (k=j+1;k<l+1;k++)
            g += a[k][j]*a[i][k];
          e[j]=g/h;
          f += e[j]*a[i][j];
        }
        hh=f/(h+h);
        for (j=0;j<l+1;j++) {
          f=a[i][j];
          e[j]=g=e[j]-hh*f;
          for (k=0;k<j+1;k++)
            a[j][k] -= (f*e[k]+g*a[i][k]);
        }
      }
    } else
      e[i]=a[i][l];
    d[i]=h;
  }

  d[0]=0.0;
  e[0]=0.0;

  for (i=0;i<n;i++) {
    l=i;
    if (d[i] != 0.0) {
      for (j=0;j<l;j++) {
        g=0.0;
        for (k=0;k<l;k++)
          g += a[i][k]*a[k][j];
        for (k=0;k<l;k++)
          a[k][j] -= g*a[k][i];
      }
    }
    d[i]=a[i][i];
    a[i][i]=1.0;
    for (j=0;j<l;j++) a[j][i]=a[i][j]=0.0;
  }
}


void Statistics::EV_tqli( Data::Vector<double> & d , 
			  Data::Vector<double> & e , 
			  Data::Matrix<double> & z )
{
  int m,l,iter,i,k;
  double s,r,p,g,f,dd,c,b;
  
  int n=d.dim1();
  for (i=1;i<n;i++) e[i-1]=e[i];
  e[n-1]=0.0;
  for (l=0;l<n;l++) {
    iter=0;
    do {
      for (m=l;m<n-1;m++) {
        dd=fabs(d[m])+fabs(d[m+1]);
        if (fabs(e[m])+dd == dd) break;
      }
      if (m != l) {
        if (iter++ == 30) Helper::halt("Internal problem in tqli routine");
        g=(d[l+1]-d[l])/(2.0*e[l]);
        r=Statistics::pythag(g,1.0);
        g=d[m]-d[l]+e[l]/(g+Statistics::SIGN(r,g));
        s=c=1.0;
        p=0.0;
        for (i=m-1;i>=l;i--) {
          f=s*e[i];
          b=c*e[i];
          e[i+1]=(r=Statistics::pythag(f,g));
          if (r == 0.0) {
            d[i+1] -= p;
            e[m]=0.0;
            break;
          }
          s=f/r;
          c=g/r;
          g=d[i+1]-p;
          r=(d[i]-g)*s+2.0*c*b;
          d[i+1]=g+(p=s*r);
          g=c*r-b;

          for (k=0;k<n;k++) {
            f=z[k][i+1];
            z[k][i+1]=s*z[k][i]+c*f;
            z[k][i]=c*z[k][i]-s*f;
          }
        }
        if (r == 0.0 && i >= l) continue;
        d[l] -= p;
        e[l]=g;
        e[m]=0.0;
      }
    } while (m != l);
  }
}


Data::Matrix<double> Statistics::matrix_multiply( const Data::Matrix<double> & a, const Data::Matrix<double> & b )
{
  //  int ar ac x br bc
  if ( a.dim2() != b.dim1() ) Helper::halt("non-conformable matrix multiplication requested");     

  const int nrow = a.dim1();
  const int ncol = b.dim2();
  const int nk = a.dim2();
  Data::Matrix<double> r(nrow,ncol);
  for (int i=0;i<nrow;i++)
    for(int j=0; j<ncol;j++)
      for (int k=0; k<nk; k++)
	r(i,j) += a(i,k) * b(k,j); 	
  return r;
}

Data::Vector<double> Statistics::matrix_multiply( const Data::Matrix<double> & a , const Data::Vector<double> & b)
{
  if ( a.dim2() != b.dim1() ) Helper::halt("non-conformable matrix multiplication requested");     
  Data::Vector<double> r( a.dim1() );
  const int nrow = a.dim1();
  const int nk = a.dim2();
  for (int i=0;i<nrow;i++)
    for (int k=0; k<nk; k++)
      r(i) += a(i,k) * b(k); 	
  return r;
}

Data::Vector<double> Statistics::matrix_multiply( const Data::Vector<double> & a , const Data::Matrix<double> & b )
{
  if ( a.dim1() != b.dim1() ) Helper::halt("non-conformable matrix multiplication requested");     
  Data::Vector<double> r( b.dim2() );
  const int nrow = b.dim2();
  const int nk = a.dim1();
  for (int i=0;i<nrow;i++)
    for (int k=0; k<nk; k++)
      r(i) += a(k) * b(k,i);
  return r;
}

double Statistics::matrix_inner_product( const Data::Vector<double> & a , const Data::Vector<double> & b )
{
  if ( a.dim1() != b.dim1() ) 
    {
      plog.warn("internal error: non-comformable inner-product" );
      return 0;
    }
  double r;
  for (int i = 0 ; i < a.dim1() ; i++ )
    r += a[i] * b[i];
  return r;
}

Data::Matrix<double> Statistics::matrix_outer_product( const Data::Vector<double> & a , const Data::Vector<double> & b )
{
  Data::Matrix<double> r(a.dim1() , b.dim1() );
  for (int i = 0 ; i < r.dim1() ; i++ )
    for (int j = 0 ; j < r.dim2() ; j++ )
      r[i][j] = a[i] * b[j];
  return r;
}



long unsigned int Statistics::factorial(int n)
{
  long unsigned int z = 1;
  for (int i = 1 ; i <= n ; i++)
    z *= i; 
  return z;
}

long unsigned int Statistics::combin(int n, int k)
{
  if (k>n) return 0;
  long double z = 1;
  int r = k;
  if ( k > (n-k) ) r = n-k;
  for (int i = 0 ; i <= r-1 ; i++)
    z *= (long double)(n-i) / (long double)(r-i);
  return (long unsigned int)z;
}

long double factorial(int x) {
  int i;
  long double result = 1;
  for (i = 2; i <= x; i++)
    result *= i;
  return result;
}


double Statistics::dbinom( int k , int n , double p )
{
  return bico( n , k ) * pow( p , k ) * pow( 1 - p , n - k ) ;
}


double Statistics::gammln(double xx)
{
  
  //  Returns the value ln[Γ(xx)] for xx > 0. 
  
  // Internal arithmetic will be done in double precision, a nicety
  // that you can omit if five-figure accuracy is good enough.
  
  static double cof[6]={76.18009172947146,-86.50532032941677,
			24.01409824083091,-1.231739572450155,
			0.1208650973866179e-2,-0.5395239384953e-5}; 
  
  int j;
  double y=xx, x=xx; 
  double tmp=x+5.5; 
  tmp -= (x+0.5)*log(tmp); 
  double ser=1.000000000190015; 
  for (j=0;j<=5;j++) ser += cof[j]/++y; 
  return -tmp+log(2.5066282746310005*ser/x);
}

double Statistics::factrl(int n)
{  
  static int ntop=4;
  static double a[33]={1.0,1.0,2.0,6.0,24.0}; // Fill in table only as required.
  int j;
  if (n < 0) { std::cerr << "exit1\n"; }
  if (n > 32 ) return exp(gammln(n+1));
  while (ntop<n) { 
    j=ntop++;
    a[ntop]=a[j]*ntop;
  }
  return a[n];
}

double Statistics::bico(int n, int k)
{  
  // binomial coeff
  return floor( 0.5 + exp( factln(n) -factln(k) -factln(n-k) ) );
}

double Statistics::factln(int n)
{
  static double a[101];
  if (n <= 1) return 0.0;
  if (n <= 100) return a[n] ? a[n] : (a[n]=gammln(n+1.0));
  else return gammln(n+1.0);
}
