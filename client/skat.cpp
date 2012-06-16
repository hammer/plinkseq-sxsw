#include "plinkseq.h"
#include "genic.h"

#include "plinkseq/prob.h"
#include "davies.h"

#include <iostream>
#include <functional>
#include <cmath>
#include <cstdlib>
#include <cstddef>
#include <iomanip>
#include <list>


extern GStore g;

using namespace std;

bool Pseq::Assoc::Aux_skat::precalculated = false;

Data::Vector<double> Pseq::Assoc::Aux_skat::y; // phenotype --> modified (y_i - u_i)
Data::Matrix<double> Pseq::Assoc::Aux_skat::Y; // phenotype --> modified (y_i - u_i)
Data::Vector<double> Pseq::Assoc::Aux_skat::u; // phenotype --> modified (u_i)

// original, non-permuted values
Data::Vector<double> Pseq::Assoc::Aux_skat::y_orig; // phenotype --> modified (y_i - u_i)
Data::Matrix<double> Pseq::Assoc::Aux_skat::Y_orig; // phenotype --> modified (y_i - u_i)
Data::Vector<double> Pseq::Assoc::Aux_skat::u_orig; // phenotype --> modified (u_i)

std::vector<double> Pseq::Assoc::Aux_skat::rho; // for SKAT-O -- weight between burden and original SKAT
int Pseq::Assoc::Aux_skat::rho_est;
int Pseq::Assoc::Aux_skat::nr;

Data::Matrix<double> Pseq::Assoc::Aux_skat::X; // covariates
std::vector<bool> Pseq::Assoc::Aux_skat::mask; // inclusion mask
int Pseq::Assoc::Aux_skat::n_actual = 0;       // actual # of individuals

bool Pseq::Assoc::Aux_skat::logistic_model = true; // otherwise linear model for phenotype

bool Pseq::Assoc::Aux_skat::has_covar = false;
std::vector<std::string> Pseq::Assoc::Aux_skat::covars;
bool Pseq::Assoc::Aux_skat::has_weights = false;
std::string Pseq::Assoc::Aux_skat::weights = "";

bool Pseq::Assoc::Aux_skat::use_freq_weights = false;
double Pseq::Assoc::Aux_skat::a1 = 1;
double Pseq::Assoc::Aux_skat::a2 = 25;

// helper functions/structs
double SKAT_sub_integrate_func_davies(double x, void * data, bool * okay);
double SKAT_sub_integrate_func_liu(double x, void * data, bool * okay);

struct SKAT_param { 	
  SKAT_param( double muQ, double varQ, double kerQ, double df , double varRemain , 
	      const Data::Vector<double> & tau , 
	      const Data::Vector<double> & lambda ) 
    : muQ(muQ), varQ(varQ) , kerQ(kerQ), df(df), varRemain(varRemain), tau(tau), lambda(lambda) { } 
  
  double muQ, varQ, kerQ, varRemain, df;
  Data::Vector<double> tau;
  Data::Vector<double> lambda;  
};


struct SKAT_aux_bundle { 
  SKAT_aux_bundle(const Data::Vector<double> * pminq, const SKAT_param * param ) : pminq(pminq) , param(param) { } 
  const Data::Vector<double> * pminq;
  const SKAT_param * param;
};



void Pseq::Assoc::Aux_skat::fit_null()
{

  const int n = g.indmap.size();
  
  // Get # of people with non-missing phenotypes or covariates
  
  y.clear();
  X.clear();
  
  mask.resize( n , false );
  
  // Covariates?
  Data::Matrix<double> C;
  if ( has_covar ) 
    C = g.phmap.covariates( covars , g.indmap );
  
  if      ( g.phmap.type() == PHE_QT ) logistic_model = false;
  else if ( g.phmap.type() == PHE_DICHOT ) logistic_model = true;
  else Helper::halt("SKAT assumes a binary or quantitative phenotype");
  

  n_actual = 0 ;
  
  for (int i=0; i< n; i++)
    {      

      Individual * person = g.indmap.ind(i);
      
      if ( person->missing() )
	mask[i] = true; //mask out      
      else if ( C.masked(i) ) 
	mask[i] = true;
      else
	{
	  ++n_actual;
	  
	  // phenotype
	  if ( logistic_model ) 
	    y.push_back( person->affected() == CASE ? 1 : 0 );
	  else
	    y.push_back( person->qt() );
	  
	}            
    }
  
  const int ncov = C.dim2();

  X.resize( n_actual , 1 + ncov );
  int j=0;
  for (int i=0;i<n; i++)
    {
      if ( ! mask[i] ) 
	{
	  X(j,0) = 1; // intercept
	  for (int c=0;c<ncov;c++) X(j,c+1) = C[i][c];
	  ++j;
	}
    }

  
  // Fit null model

  GLM glm( logistic_model ? GLM::LOGISTIC : GLM::LINEAR );
  
  glm.set( y , X ); 
  
  glm.fit();
  
  if ( ! glm.valid() ) 
    Helper::halt( "problem attempting to fit the null model for SKAT" );
  
  // Get vector of coefficients
  
  Data::Vector<double> beta;    
  if ( ! glm.display(&beta) ) 
    Helper::halt( "problem fitting null model for SKAT" );
   
   // by default, glm() returns OR scale; so return...
   if ( logistic_model ) 
     for (int j=0;j<beta.size();j++) beta[j] = log(beta[j]);
  
  // Calculate residuals (now in a reduce format of only the 
  Data::Vector<double> o = y;
  y.resize( n_actual );
  u.resize( n_actual );
  Y.resize( 1 , n_actual );

  j=0;

  for (int i=0;i<n;i++)
    {
      if ( ! mask[i] ) 
	{
	  y[j] = beta[0];
	  
	  for (int c=0;c<ncov;c++) 
	    y[j] += C(i,c) * beta[c+1];
	  
	  // Y=logit^-1(Y)
	  if ( logistic_model ) 
	    {	      
	      double e = exp( y[j] );
	      y[j] = e / ( 1 + e ) ;	      
	    }

	  // make y-u and store
	  
	  if ( logistic_model ) 
	    u[j] = y[j] * (1-y[j]);
	  // otherwise sigma^2 calculated below to populate u[i]
	  
	  y[j] = o[j] - y[j];
	  Y(0,j) = y[j]; // fudge : duplicate to make matrix mult work... 

	  ++j;
	}
    }
  
  // we need to calculate the variance of the residuals and populate
  // u[i] accordingly:   sigma^2 = 1/(n-p) Sum(w[i] R[i]^2),              

  
  if ( ! logistic_model ) 
    {
      double mean = 0;      
      for ( int i=0;i<n_actual;i++) 
	mean += y[i];
      mean /= (double)n_actual;

      double sigma2 = 0;
      for ( int i=0;i<n_actual;i++) 
	sigma2 += ( y[i] - mean ) * ( y[i] - mean );
      sigma2 /= (double)( n_actual - ( ncov + 1 ) );
      for ( int i=0;i<n_actual;i++) u[i] = sigma2;
    }
 
  // Store originals
  y_orig = y;
  Y_orig = Y;
  u_orig = u;
  
  // Indicate that we have now fit the null model (and so do not need to do it again)
  
  precalculated = true;

}


void Pseq::Assoc::Aux_skat::populate_G( const VariantGroup & vars , Aux_prelim * aux )
{

  const int ni = vars.n_individuals();
  const int nv = vars.size();
  
  if ( mask.size() != ni ) Helper::halt( "internal error in SKAT -- mask != N" );
  
  // genotype matrix
  G.resize( n_actual , nv );

  int j = 0;
  for (int i=0;i<ni;i++)
    {
      if ( ! mask[i] )
	{
	  for (int v=0;v<nv;v++)
	    {
	      const Genotype & genotype = vars.geno(v,i);

	      // Use mean-imputation for missing data (assume diploid counts)

	      if ( genotype.null() ) 
		G(j,v) = has_weights ? w[v] * 2 * aux->maf[v] : 2 * aux->maf[v];
	      else 
		{
		  if ( has_weights )
		    G(j,v) = w[v] * genotype.minor_allele_count( aux->altmin[v] ) ;
		  else
		    G(j,v) = genotype.minor_allele_count( aux->altmin[v] ) ;
		}
	    }
	  ++j;
	}
    }
}



void Pseq::Assoc::Aux_skat::populate_K()
{

  // NOTE: never explicitly called.

  // K = G.W.G'  
  // K_{i,i'} = sum_v w_v G_iv G_i'v 
  // K is NxN matrix
  
  const int ni = G.dim1();
  const int nv = G.dim2();

  K.resize( ni , ni );
  for (int i=0;i<ni;i++)
    for (int j=i;j<ni;j++)
      {
	double s = 0;
	for (int v=0;v<nv;v++)
	  s += w[v] * G(i,v) * G(j,v);
	K(i,j) = s;
	K(j,i) = s;
      }
}
  

double Pseq::Assoc::Aux_skat::calculate_Q( Data::Matrix<double> * pW1 )
{

  const int nv = G.dim2();
    
  //
  // Q = (y-u)'K(y-u)
  //

  // Calculated as: ( (Y-U)'.Z ) . (  (Y-U)'.Z ) )
 
  Data::Matrix<double> Qtemp = Y * G;

  double Q = 0;
  for (int j=0;j<nv;j++) 
    Q += Qtemp(0,j) * Qtemp(0,j) * 0.5;
  
  // X1     n x p matrix of covariates, with column of 1's
  // pi_1   vector of U  (u)
  
  //  W.1 = t(Z) %*% (Z * pi_1) - (t(Z * pi_1) %*%X1) %*% solve(t(X1)%*%(X1 * pi_1)) %*% (t(X1) %*% (Z * pi_1)) # t(Z) P0 Z
  //        -------------------   -------------------     --------------------------     ----------------------
  //        a                     b                       c                              d


  Data::Matrix<double> W1a( nv , nv );
  for (int i=0;i<nv;i++)
    for (int j=0;j<nv;j++)
      for (int k=0;k<n_actual;k++)
	W1a(i,j) += G(k,i) * G(k,j) * u[k];
  
  const int nx = X.dim2();
 
  Data::Matrix<double> W1b( nv , nx );
  for (int i=0;i<nv;i++)
    for (int j=0;j<nx;j++)
      for (int k=0;k<n_actual;k++)
	W1b(i,j) += G(k,i) * u[k] * X(k,j);

  Data::Matrix<double> W1c( nx , nx );  
  for (int i=0;i<nx;i++)
    for (int j=0;j<nx;j++)
      for (int k=0;k<n_actual;k++)
	W1c(i,j) += X(k,i) * X(k,j) * u[k];
  W1c = Statistics::inverse( W1c );

  Data::Matrix<double> W1d( nx , nv );
  for (int i=0;i<nx;i++)
    for (int j=0;j<nv;j++)
      for (int k=0;k<n_actual;k++)
	W1d(i,j) += X(k,i) * G(k,j) * u[k];


  // keep track of W1 matrix

  *pW1 = W1a - W1b * W1c * W1d;
  
  return Q;

}




//
// Primary driver for both the basic SKAT test and SKAT-O
//

double Pseq::Assoc::stat_skat( const VariantGroup & vars , 
			       Aux_prelim * aux , 
			       Aux_skat * aux_skat , 
			       std::map<std::string,std::string> * output , 
			       bool original )
{
  
  const int nv = vars.size();  

  
  if ( original ) 
    {
      
      //
      // Fit the null model (only done once per sample/phenotype)
      //
      
      if ( ! aux_skat->precalculated ) 
	aux_skat->fit_null();
      
      //
      // Create weights?
      //

      aux_skat->w.resize( nv , 1 );
    
      Aux_skat::use_freq_weights = true;

      if ( Aux_skat::use_freq_weights ) 
	{

	  // MAF (from whole sample) should be in aux_prelim
	  
	  for (int v=0;v<nv;v++)
	    aux_skat->w[v] = Helper::PROB::beta_pdf( aux->maf[v] , Aux_skat::a1 , Aux_skat::a2 );
	  
	  Aux_skat::has_weights = true;
	}
      else if ( Aux_skat::has_weights ) 
	{
	  
	  meta_index_t midx = MetaInformation<VarMeta>::field( Aux_skat::weights );
	  
	  if ( midx.mt != META_FLOAT ) Helper::halt( "can only use Float variant tags as weights" );
	  
	  for (int v=0;v<nv;v++)
	    {
	      if ( vars(v).meta.has_field( Aux_skat::weights ) )
		aux_skat->w[v] = vars(v).meta.get1_double( Aux_skat::weights );
	      else if ( vars(v).consensus.meta.has_field( Aux_skat::weights ) )
		aux_skat->w[v] = vars(v).consensus.meta.get1_double( Aux_skat::weights );
	      else
		aux_skat->w[v] = 0;
	    }
	}
    }



  //
  // Permute residuals not the original  ( PERMUTATION NOT IMPLEMENTED YET )
  //

  
//   if ( original ) 
//     {
//       u = u_orig;
//       Y = Y_orig;
//       y = y_orig;      
//     }
//   else
//     {    
//       for (int i=0;i< Aux_skat::mask.size();i++)
// 	{
// 	  // hmm.. how to handle missing / masked instances here??
// 	  if ( Aux_skat::mask[i] ) 
// 	    {
// 	      u[ g.perm.pos(i) ] = u[ i ];
// 	      u[ g.perm.pos(i) ] = u[ i ];
// 	      u[ g.perm.pos(i) ] = u[ i ];
// 	    }
// 	}
//     }




  //
  // Populate G matrix
  //
  
  aux_skat->populate_G( vars , aux );
  
  
  //
  // Create kernel K (not done explicitly)
  //
  
  //  aux_skat->populate_K();


  
  //
  // Calculate Q 
  //

  Data::Matrix<double> W;

  double Q = aux_skat->optimal_mode() ? 
    aux_skat->calculate_optimal_Q( &W ) : 
    aux_skat->calculate_Q( &W );
  
   
  //
  // Asympotic P-value
  //
  
  double pvalue = aux_skat->optimal_mode() ? 
    Q :   // skat-O above actually returns p-value already
    aux_skat->calculate_pvalue( Q , W );
  
 

  //
  // Output
  //

  if ( original ) 
    {
      
      const std::string label = aux_skat->optimal_mode() ? "SKAT-O" : "SKAT" ; 

      if ( pvalue < 0 ) 
	(*output)[ label ] = "P=NA;";
//       else
// 	(*output)[ label ] = "P=" + Helper::dbl2str( pvalue ) + ";";
      
      if ( aux_skat->optimal_mode() )
	{
	  (*output)[ label ] += "RHO=" + Helper::dbl2str( Aux_skat::rho[ Aux_skat::rho_est ] ) + ";";
	}
      if ( aux_skat->logistic_model )
	{
	  std::map<std::string,int>::iterator i = aux->mc_a.begin();
	  while ( i != aux->mc_a.end() )
	    {
	      if ( i != aux->mc_a.begin() ) (*output)[ label ] += ";";
	      (*output)[ label ] += i->first + "(" + Helper::int2str( i->second ) + ")";
	      ++i;
	    }	  
	}      
    }


  // save is auxillary structure, so we can use in output
  
  aux_skat->returned_pvalue = pvalue < 0 ? 1.00 : pvalue ; 
  
  return aux_skat->returned_pvalue;
  
}


double Pseq::Assoc::Aux_skat::get_liu_pval( double Q , Data::Matrix<double> & A1 )
{
   
  Data::Matrix<double> A2 = A1 * A1; 
  
  double c1 = 0, c2 = 0 , c3 = 0 , c4 = 0;
  
  for (int r = 0; r < A1.dim1(); r++ )
    {
      c1 += A1(r,r);
      c2 += A2(r,r);
    }
  
  for (int r = 0; r < A1.dim1(); r++ )
    for (int c = 0; c < A1.dim2(); c++ )
      {
	c3 += A1(r,c) * A2(c,r);
	c4 += A2(r,c) * A2(c,r);
      }

  double l;  // df
  double d;  // ncp
  double muQ, sigmaQ, muX, sigmaX;

  get_liu_param( false , // do not use modified form
		 c1, c2, c3, c4 , 
		 &muX, &sigmaX, 
		 &muQ, &sigmaQ, 
		 &l, &d );
  
  double qnorm  = ( Q - muQ ) / sigmaQ;
  double qnorm1 = qnorm * sigmaX + muX;
  double pvalue = Statistics::noncentral_chi2_prob(  qnorm1 , l , d );

  return pvalue; 
}


void Pseq::Assoc::Aux_skat::get_liu_param( bool mod , 
					   double c1 , double c2 , double c3 , double c4 , 
					   double * muX , double * sigmaX , 
					   double * muQ , double * sigmaQ , 
					   double * l , double * d ) 
{
  
  *muQ = c1;
  *sigmaQ = sqrt( 2 * c2 );

  double s1 = c3 / pow( c2 , 1.5 );
  double s2 = c4 / (c2*c2) ;
  
  // sqrt(8) = 2.828427
  double beta1 = 2.828427 * s1;
  double beta2 = 12       * s2;

  double a;

  if( s1*s1 > 2 ) 
    {
      a = 1.0 / (s1 - sqrt( s1*s1 - s2));
      *d = s1 * a*a*a - a*a;
      *l = a*a - 2* (*d);      
    } 
  else 
    {
      if ( mod ) // Modified Liu 
	{
	  *l = 1.0/s2;
	  a = sqrt(*l);
	  *d = 0;
	}
      else
	{
	  a = 1.0/s1;
	  *d = 0;
	  *l = 1.0/(s1*s1);
	}
    }
  
  *muX = (*l) + (*d);
  
  //sqrt(2) = 1.414214
  *sigmaX = 1.414214 * a;
  
}
    



double Pseq::Assoc::Aux_skat::calculate_pvalue( double Q , Data::Matrix<double> & pW )
{


  // Use K = W / 2 for both Liu and Davies methods 

  Data::Matrix<double> K = pW;
  
  for (int r=0;r<K.dim1();r++)
    for (int c=0;c<K.dim2();c++)
      K(r,c) *= 0.5;
  

  // Liu original method
  
  double pvalue_liu = get_liu_pval( Q , K );

  
  // Davies method p-value

  bool okay = true;
  Data::Vector<double> eigen = Statistics::eigenvalues( K , & okay );

  // if convergence issues, use liu p-value
  if ( ! okay ) return pvalue_liu;

  std::vector<bool> idx1( eigen.size() , false );
  double mean_n = 0 , mean_s = 0;
  for (int i=0; i<idx1.size(); i++) 
    { 
      if ( eigen[i] >= 0 ) 
	{ 
	  idx1[i] = true;
	  mean_n++;
	  mean_s += eigen[i];
	}
    }
  mean_s /= mean_n;
  mean_s /= 100000.0;
  
  std::vector<double> lambda;  
  for (int i=0; i<idx1.size(); i++) 
    if ( eigen[i] > mean_s ) lambda.push_back( eigen[i] ); 

  // e.g. all monomorphic variants
  if ( lambda.size() == 0 ) 
    {
      return -1.0;
    }
  
  Davies d;
  okay = true;
  std::vector<int> df( lambda.size() , 1 );
  
  double pvalue_davies = d.pvalue( Q , lambda , df , &okay );

  // ignore davies ifault for now...
  okay = true;
  
  if ( lambda.size() == 1 || pvalue_davies <= 0 || pvalue_davies > 1 ) 
    okay = false;

  if ( okay ) return pvalue_davies; 
  return pvalue_liu;

}


double Pseq::Assoc::Aux_skat::calculate_optimal_Q( Data::Matrix<double> * pW1 )
{

  // number of variants in gene/set
  const int nv = G.dim2();
  
  // number of values of 'rho' to evaluate
  const int nr = rho.size();
  

  // assume 'linear-weighted' kernel is the default (and only) option for now.
  
  // pi_1 == 'u' (N matrix of weights/scores)
  
  // Z = t(t(Z) * (weights))
  //   our G matrix should already be weighted, if there are weights; 

  
  // Z1 = (Z * sqrt(pi_1))  -  ( X1 * sqrt(pi_1))   %*%   solve( t(X1) %*% (X1 * pi_1) )  %*%  ( t(X1) %*% (Z * pi_1) )
  //      ----------------     ------------------         -----------------------------        -------------------------
  //      a                    b                          c                                    d

  // Z     --> G  ( n_actual x nv )
  // pi_1  --> u  ( n_actual )
  // X1    --> X  ( n_actual x nx )

  const int nx = X.dim2();
  Data::Vector<double> sqrt_pi( n_actual );
  for (int i=0;i<n_actual;i++) sqrt_pi[i] = sqrt( u[i] );
  
  Data::Matrix<double> Z1a( n_actual , nv );
  for (int i=0;i<n_actual;i++) 
    for (int j=0;j<nv;j++) 
      Z1a(i,j) = G(i,j) * sqrt_pi[i];

  Data::Matrix<double> Z1b( n_actual , nx );
  for (int i=0;i<n_actual;i++) 
    for (int j=0;j<nx;j++) 
      Z1b(i,j) = X(i,j) * sqrt_pi[i];
  
  Data::Matrix<double> Z1c( nx , nx  );
  for (int i=0;i<nx;i++)
    for (int j=0;j<nx;j++)
      for (int k=0;k<n_actual;k++)
	Z1c(i,j) += X(k,i) * X(k,j) * u[k] ;

  Z1c = Statistics::inverse( Z1c );  

  // Z1 = (Z * sqrt(pi_1))  -  ( X1 * sqrt(pi_1))   %*%   solve( t(X1) %*% (X1 * pi_1) )  %*%  ( t(X1) %*% (Z * pi_1) )  
  //      ----------------     ------------------         -----------------------------        -------------------------  
  //      a                    b                          c                                    d                                               
  
  Data::Matrix<double> Z1d( nx , nv );
  for (int i=0;i<nx;i++)
    for (int j=0;j<nv;j++)
      for (int k=0;k<n_actual;k++)
	Z1d(i,j) += X(k,i) * G(k,j) * u[k] ;
  
  Data::Matrix<double> Z1 = Z1a - Z1b * Z1c * Z1d;

//   double Z1_SSQ = 0;
//   for(int i=0;i<n_actual;i++)
//     for(int j=0;j<nv;j++)
//       Z1_SSQ += Z1(i,j) * Z1(i,j);


  //
  // Calculate vector of Q's, for each value of 'rho'
  //

  Data::Vector<double> Q = sub_optimal_get_Q( G );


  //
  // Get p-values and find optimal 
  //

  double pvalue = sub_optimal_get_P( Q , Z1 );

  return pvalue;

}



double Pseq::Assoc::Aux_skat::calculate_optimal_pvalue( double , Data::Matrix<double> & )
{
  return 1.00;
}

Data::Vector<double> Pseq::Assoc::Aux_skat::sub_optimal_get_Q( const Data::Matrix<double> & Z1 )
{

  set_optimal_rcorr();

  const int nr = rho.size();

  const int nv = Z1.dim2();

  Data::Vector<double> Qr( nr );

  Data::Vector<double> temp( nv ); 
  for (int i=0; i<nv; i++) 
    for (int j=0;j<n_actual;j++) 
      temp[i] += y[j] * Z1(j,i) ; 

  double sum_of_sqr = 0;
  double s = 0;

  for (int i=0; i<nv; i++) 
    {
      sum_of_sqr += ( temp[i] * temp[i] );
      s += temp[i];
    }
  
  s /= (double)nv;

  double mean_sqr = nv * nv * s * s;

  for (int i=0; i<nr; i++) 
    {
      Qr[i] = ( 1 - rho[i] ) * sum_of_sqr  +  rho[i] * mean_sqr ;
      Qr[i] *= 0.5;
    }

  return Qr;

}


double Pseq::Assoc::Aux_skat::sub_optimal_get_P( const Data::Vector<double> & Q , 
						 const Data::Matrix<double> & Z0 )
{
  
  const int nr = rho.size();
  const int nq = Q.size();   // should always equal nr, as no resampling implemented
  
  const int ni = Z0.dim1();
  const int nv = Z0.dim2();
  
  // scale Z1 by 1/sqrt(2) (see if we can do this elsewhere/once)
  Data::Matrix<double> Z1(ni,nv);
  
  const double fac = 1.0 / sqrt(2.0) ;  
  for (int i=0;i<Z1.dim1();i++)
    for (int j=0;j<Z1.dim2();j++)
      Z1(i,j) = Z0(i,j) * fac;

  
  std::vector<std::vector<double> > lambda; 

  for (int i=0;i<nr;i++)
    {
      
      // diagonal correlation matrix
      
      Data::Matrix<double> RM( nv , nv , rho[i] );
      for (int j=0;j<nv;j++) RM(j,j) = 1.0;
      
      // cholesky decomp (returns in lower. tri )       
      Data::Matrix<double> L = Statistics::cholesky( RM );       
      Data::Matrix<double> Z2 = Z1 * L;
      Data::Matrix<double> K1 = Statistics::transpose( Z2 ) * Z2;
      lambda.push_back( * get_lambda( K1 ).data_pointer() );
    }

  
//   for (int i=0;i<lambda.size(); i++)
//     {
//       std::cout << "LAMBDA " << i << "\n";
//       for (int j=0;j<lambda[i].size();j++)
// 	std::cout << " " << lambda[i][j] ;
//       std::cout <<  "\n";
//     }


   //
   // Get mixture parameters
   //
  
   // param.m<-SKAT_Optimal_Param(Z1,r.all)
  
   double muQ, varQ, kerQ, varRemain, df;
   Data::Vector<double> tau;
   Data::Vector<double> lambda2;
   
   get_optimal_param( Z1 , &muQ, &varQ, &kerQ, &varRemain, &df, &tau, &lambda2 );
   SKAT_param param( muQ, varQ, kerQ, df, varRemain, tau, lambda2 );

   if ( lambda2.size() == 0 ) 
     {
       return 1.0; 
     }


   //
   // Find optimal rho for each Q
   //

   // Each_Info<-SKAT_Optiaml_Each_Q(param.m, Q.all, r.all, lambda.all)
   // pmin.q<-Each_Info$pmin.q

   Data::Vector<double> pminq = sub_optimal_P_foreach_Q( Q , lambda );	
   
   //   std::cout << pminq.print( "pminq" ) << "\n";
   
   
   // pval<-rep(0,n.q)  
   // (note, nq == nr here as no resampling)


   
   // 
   // Calculate both Davies/Optimal and Liu approx.
   //
      
   //  if(method == "davies" || method=="optimal"){
   
   double pvalue = sub_optimal_P_Davies( pminq , param );
   
   //   std::cout << "skat-o p-value is " << pvalue << "\n";
   
   return pvalue;
}


Data::Vector<double> Pseq::Assoc::Aux_skat::get_lambda( const Data::Matrix<double> & M )
 {
   
   Data::Matrix<double> cpy = M;
   bool okay = true; 
   Data::Vector<double> e = Statistics::eigenvalues( cpy , &okay );
   
   std::vector<double> okaye;
   
   double meane = 0;
   int cnte = 0;
   for (int i=0;i<e.size();i++) 
     {
       if ( e[i] >= 0 ) { meane += e[i]; ++cnte; } 
     }
   meane /= (double)cnte;
   meane /= 100000.0;
   
   for (int i=0;i<e.size();i++)      
     if ( e[i] > meane ) okaye.push_back(e[i]);

   if ( ( ! okay ) || okaye.size() == 0 ) plog.warn("no eigenvalue greater than 0 / convergence issues");
   
   return Data::Vector<double>( okaye );

 }



void Pseq::Assoc::Aux_skat::get_optimal_param( const Data::Matrix<double> & Z1 , 
					       double * muQ, double * varQ, double * kerQ, 
					       double * varRemain, double * df, 
					       Data::Vector<double> * tau, 
					       Data::Vector<double> * lambda )
{
  
  const int ni = Z1.dim1();
  const int nv = Z1.dim2();
  const int nr = rho.size();
  
  // row means (score for each individual)

  Data::Vector<double> z_mean( ni ) ;
  for (int j=0;j<nv;j++)
    for (int i=0;i<ni;i++)
      z_mean[i] += Z1(i,j);
  for (int i=0;i<ni;i++) z_mean[i] /= (double)nv;


  //  Z_mean<-matrix(rep(z_mean,p.m),ncol=p.m,byrow=FALSE)
  // cof1<-(t(z_mean) %*% Z1)[1,] / sum(z_mean^2)

  Data::Vector<double> cof1( nv );
  double z2 = Statistics::sum_squares( z_mean );
  for (int j=0; j<nv; j++)
    {
      for (int i=0; i<ni; i++)
	cof1[j] += z_mean[i] * Z1(i,j);
      cof1[j] /= z2;
    }
  
  // Z.item1<-Z_mean %*% diag(cof1)
  // Z.item2<-Z1 - Z.item1

  Data::Matrix<double> Z_item1( ni , nv );
  Data::Matrix<double> Z_item2( ni , nv );
  
  for (int i=0; i<ni; i++)
    for (int j=0; j<nv; j++)
      {
	Z_item1(i,j) = z_mean[i] * cof1[j];
	Z_item2(i,j) = Z1(i,j) - Z_item1(i,j);
      }
  
//   std::cout << Z_item1.print( "Z1_item1" , 10 , 10 )  << "\n";  
//   std::cout << Z_item2.print( "Z1_item2" , 10 , 10 ) ;

//   double S1 = 0;
//   double S2 = 0;
//   double S0 = 0;

//   double SS0 = 0;
//   double SS1 = 0;
//   double SS2 = 0;

//   for (int i=0;i<ni;i++)
//     for (int j=0;j<nv;j++)
//       {
// 	S0 += 1e10 * Z1(i,j);
// 	S1 += 1e10 * Z_item1(i,j);
// 	S2 += 1e10 * Z_item2(i,j);

// 	SS0 += Z1(i,j) * Z1(i,j);
// 	SS1 += Z_item1(i,j) * Z_item1(i,j);
// 	SS2 += Z_item2(i,j) * Z_item2(i,j);
	
//       }
  
//   std::cout <<" summs = " << S0 << " " << S1 << " " << S2 << "\n";
//   std::cout <<" summs^2 = " << SS0 << " " << SS1 << " " << SS2 << "\n";



  // # W3.2 Term : mixture chisq
  // W3.2.t <- t(Z.item2) %*% Z.item2
  // lambda<-Get_Lambda(W3.2.t)

  
  //  std::cout << Z_item2.print( "Z_item2" , 10 , 10 ) << "\n";

  Data::Matrix<double> W3_2_t = Statistics::transpose( Z_item2 ) * Z_item2; 

  *lambda = get_lambda( W3_2_t );

 
  // # W3.3 Term : variance of remaining ...
  // W3.3.item<- sum( (t(Z.item1) %*% Z.item1) * (t(Z.item2) %*% Z.item2) ) * 4

  Data::Matrix<double> tmp = Statistics::transpose( Z_item1 ) * Z_item1;

  double W3_3_item = 0;
  for (int i=0;i<nv;i++)
    for (int j=0;j<nv;j++)
      W3_3_item += W3_2_t(i,j) * tmp(i,j);  
  W3_3_item *= 4.0;
  
  // # Mixture Parameters
  
  //     MuQ  <- sum(lambda)
  //     VarQ <- sum(lambda^2) *2 + W3.3.item
  //     KerQ <- sum(lambda^4)/(sum(lambda^2))^2 * 12
  //     Df   <- 12/KerQ
  
  for (int i=0;i< lambda->size();i++)
    {      
      *muQ  += (*lambda)[i];
      double sqr = (*lambda)[i] * (*lambda)[i];
      *varQ += sqr;
      *kerQ += sqr * sqr;
    }
  
  *kerQ /= *varQ * *varQ ; 
  *kerQ *= 12;
  
  *varQ *= 2;
  *varQ += W3_3_item;
  
  *df = 12.0 / *kerQ; 
  

  // # W3.1 Term : tau1 * chisq_1
  // tau<-rep(0,r.n)

  tau->clear();
  tau->resize( nr );

  double sum_cof_sqr = Statistics::sum_squares( cof1 );
  double sum_z_mean_sqr = Statistics::sum_squares( z_mean ); 

  for (int i=0; i<nr; i++)
    {
      // term1<-p.m^2*r.corr + sum(cof1^2) * (1-r.corr)
      // tau[i]<-sum(term1) *  sum(z_mean^2)
      
      (*tau)[i] = ( nv * nv * rho[i] + sum_cof_sqr * ( 1 - rho[i] ) ) * sum_z_mean_sqr;
    }
  
  *varRemain = W3_3_item;  
}

Data::Vector<double> Pseq::Assoc::Aux_skat::sub_optimal_P_foreach_Q( const Data::Vector<double> & Q , 
								     const std::vector<std::vector<double> > & lambda )
{
  
  // SKAT_Optiaml_Each_Q<-function(param.m, Q.all, r.all, lambda.all){

  // n.r<-length(r.all)
  // c1<-rep(0,4)
  // n.q<-dim(Q.all)[1]

  const int nr = rho.size();
  const int nq = 1; // only ever 1, as no resampling (otherwise, # rows of Q)

  // pval<-matrix(rep(0,n.r*n.q),ncol=n.r)
  // pmin.q<-matrix(rep(0,n.r*n.q),ncol=n.r)
  // param.mat<-NULL

  // likewise, if no resampling, these are vectors, not matrices

  Data::Vector<double> pval( nr );
  Data::Vector<double> pminq( nr );
  Data::Matrix<double> param( nr , 3 ); // stores, mu, var, Df per rho/p
  
  for (int i = 0; i< nr; i++)
    {

      // Q<-Q.all[,i]

      double q = Q[i];

      // r.corr<-r.all[i]
      
      // c1[1]<-sum(lambda.temp)
      // c1[2]<-sum(lambda.temp^2)
      // c1[3]<-sum(lambda.temp^3)
      // c1[4]<-sum(lambda.temp^4)
      
      double c1 = 0 , c2 = 0 , c3 = 0 , c4 = 0;
      
      const std::vector<double> & ll = lambda[i];
      const int nl = ll.size();

      for (int j=0; j<nl; j++)
	{
	  double t = ll[j] * ll[j];
	  c1 += ll[j];
	  c2 += t;
	  c3 += t * ll[j];
	  c4 += t * t;
	}
      
      
      // param.temp<-Get_Liu_Params_Mod(c1)
	
      double l;  // df
      double d;  // ncp
      double muQ, sigmaQ, muX, sigmaX;
      
      get_liu_param( true  , // use modified form
		     c1, c2, c3, c4 , 
		     &muX, &sigmaX, 
		     &muQ, &sigmaQ, 
		     &l, &d );
      
      // muQ<-param.temp$muQ
      // varQ<-param.temp$sigmaQ^2
      // df<-param.temp$l
      
      // 
      // Get p-value
      //

      double df = l;
      double varQ = sigmaQ * sigmaQ;
      double sdQ = sqrt(varQ);

      // Q.Norm<-(Q - muQ)/sqrt(varQ) * sqrt(2*df) + df
      // pval[,i]<- 1-pchisq(Q.Norm,  df = df)
      // param.mat<-rbind(param.mat,c(muQ,varQ,df))
    
      double qnorm = ( Q[i] - muQ ) / sdQ * sqrt( 2 * df ) + df ;
      pval[i] = Statistics::chi2_prob( qnorm , df );      
      param(i,0) = muQ;
      param(i,1) = varQ;
      param(i,2) = df;
      
    }

  //
  // get minumum p-value
  //

  // pmin<-apply(pval,1,min)
  
  double pmin = pval[0];
  rho_est = 0;
  for (int i=1;i<nr;i++) 
    { 
      if ( pval[i] < pmin ) 
	{
	  pmin = pval[i]; 
	  rho_est = i;
	}
    }
  
  for (int i=0;i<nr;i++)
    {
      
      double muQ  = param(i,0);
      double varQ = param(i,1);
      double df   = param(i,2);

      double q_org = 0;
      bool okay = Statistics::qchisq( pmin , df , &q_org );

      // q.org <-qchisq(1-pmin,df=df)
      // q.q<-(q.org - df)/sqrt(2*df) *sqrt(varQ) + muQ
      // pmin.q[,i]<-q.q

      pminq[i] = ( q_org - df ) / sqrt( 2 * df ) * sqrt( varQ ) + muQ ;

    }
  
  //out<-list(pmin=pmin,pval=pval,pmin.q=pmin.q)
  
  return pminq;
  
}



double Pseq::Assoc::Aux_skat::sub_optimal_P_Davies( const Data::Vector<double> & pminq , const SKAT_param & param )
{

  //  re<-try(integrate(SKAT_Optimal_Integrate_Func_Davies, lower=0, upper=30, subdivisions=500,pmin.q=pmin.q,param.m=param.m,r.all=r.all,abs.tol \
  // 		    = 10^-15), silent = TRUE)

  SKAT_aux_bundle a(&pminq,&param);  
  bool okay = true;
  /// double pvalue = 1 - Statistics::integrate( 0 , 30 , SKAT_sub_integrate_func_davies , &okay, (void*)&a , 10e-15 );  
  double pvalue = 1 - Statistics::qromo( 0 , 30 , SKAT_sub_integrate_func_davies , &okay, (void*)&a , 3, 10e-8 );  
      
  if ( ! okay ) return sub_optimal_P_Liu( pminq , param ); 
  if ( pvalue <= 0.0 || pvalue > 1.0 ) return sub_optimal_P_Liu( pminq , param ); 
  return pvalue;
}

double Pseq::Assoc::Aux_skat::sub_optimal_P_Liu( const Data::Vector<double> & pminq , const SKAT_param & param )
{
  
  //   re<-integrate(SKAT_Optimal_Integrate_Func_Liu, lower=0, upper=30, subdivisions=500
  // 		,pmin.q=pmin.q,param.m=param.m,r.all=r.all,abs.tol = 10^-15)

  SKAT_aux_bundle a(&pminq,&param);
  bool okay = true;
  // double pvalue = 1 - Statistics::integrate( 0 , 30 , SKAT_sub_integrate_func_liu, &okay, (void*)&a , 10e-15 );
  double pvalue = 1 - Statistics::qromo( 0 , 30 , SKAT_sub_integrate_func_liu , &okay, (void*)&a , 3, 10e-8 );  
  if ( ! okay ) plog.warn("internal problem in SKAT integrate_fun_liu()");
  return pvalue; 
}


// Integration functions for Davies and Liu

double SKAT_sub_integrate_func_liu(double x, void * data, bool * okay )
{

  SKAT_aux_bundle * aux = (SKAT_aux_bundle*)data;

  const int nr = Pseq::Assoc::Aux_skat::rho.size();

  if ( x == 0 ) x = 1e-20;

//     temp1<-param.m$tau %x% t(x)
//     temp<-(pmin.q - temp1)/(1-r.all)
//     temp.min<-apply(temp,2,min)
  
  Data::Vector<double> temp = aux->param->tau;

  for (int i=0;i<nr; i++) temp[i] = ( (*aux->pminq)[i] - temp[i] * x ) / ( 1 - Pseq::Assoc::Aux_skat::rho[i] );

  double min1 = temp[0];

  for (int i=1;i<nr;i++) if ( temp[i] < min1 ) min1 = temp[i];

  //     temp.q<-(temp.min - param.m$MuQ)/sqrt(param.m$VarQ)*sqrt(2*param.m$Df) + param.m$Df  

  double tempq = ( min1 - aux->param->muQ ) / sqrt( aux->param->varQ) * sqrt( 2 * aux->param->df ) + aux->param->df;

  //     re<-pchisq(temp.q ,df=param.m$Df) * dchisq(x,df=1)
  
  double t1 = tempq > 0 ? 1 - Statistics::chi2_prob( tempq , aux->param->df ) : 0 ;
  double t2 = Helper::PROB::chi_square_pdf ( x , 1 ) ;
    
  //  std::cout << "ret value = " << x << "\t" << t1 * t2 << "\t" << t1 << "\t" << t2 << "\t" << tempq << "\n";

  return t1 * t2 ;

}

double SKAT_sub_integrate_func_davies(double x, void * data, bool * okay )
{

  SKAT_aux_bundle * aux = (SKAT_aux_bundle*)data;
  
  const int nr = Pseq::Assoc::Aux_skat::rho.size();
  
  // temp1<-param.m$tau %x% t(x)
  // temp<-(pmin.q - temp1)/(1-r.all)
  // temp.min<-apply(temp,2,min)

  Data::Vector<double> temp = aux->param->tau;

  for (int i=0;i<nr; i++) temp[i] = ( (*aux->pminq)[i] - temp[i] * x ) / ( 1 - Pseq::Assoc::Aux_skat::rho[i] );
  double min1 = temp[0];
  for (int i=1;i<nr;i++) if ( temp[i] < min1 ) min1 = temp[i];
  
  // if(min1 > sum(param.m$lambda) * 10^4){
  // temp<-0
  // } else {

  double t;

  if( min1 > Statistics::sum( aux->param->lambda ) * 100000 ) 
    {
      t = 0; 
    } 
  else 
    {

      // min1.temp <- min1 - param.m$MuQ
      
      double min1temp = min1 - aux->param->muQ;

      // sd1<-sqrt(param.m$VarQ - param.m$VarRemain)/sqrt(param.m$VarQ)

      double sd1 = sqrt( aux->param->varQ - aux->param->varRemain ) / sqrt( aux->param->varQ );
      
      // min1.st<-min1.temp *sd1 + param.m$MuQ
	
      double min1st = min1temp * sd1 + aux->param->muQ;

      // dav.re =  SKAT_davies(min1.st,param.m$lambda,acc=10^(-6))

      Davies d;
      bool okay2 = true;
      std::vector<int> df( nr , 1 );      

      std::vector<double> tl = *aux->param->lambda.data_pointer();
      std::vector<bool> ifaults;
      t = d.pvalue( min1st , tl , df , &okay2 , &ifaults );

      if ( ! okay2 ) 
	{
	  //	  std::cout << "bad value\n";
	  *okay = false;
	  return 0;
	}

      if ( t > 1 ) t = 1 ;       

//       if ( ! okay ) 
// 	{
// 	  plog.warn("problem in Davies()");
// 	  for (int z=0;z<6;z++) std::cout << "ifault " << z << " " << ifaults[z] << "\n";	  
// 	  // ifault = 1       required accuracy NOT achieved
// 	  //          2       round-off error possibly significant
// 	  //          3       invalid parameters
// 	  //          4       unable to locate integration parameters
// 	  //          5       out of memory
// 	}

    }
  

  // re[i]<-(1-temp) * dchisq( x[i] ,df=1)

  return (1-t) * Helper::PROB::chi_square_pdf ( x , 1 ) ;
  
}
