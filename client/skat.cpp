#include "pseq.h"
#include "genic.h"

#include "../lib/prob.h"
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

double Pseq::Assoc::stat_skat( const VariantGroup & vars , 
			       Aux_prelim * aux , 
			       Aux_skat * aux_skat , 
			       std::map<std::string,std::string> * output , 
			       bool original )
{
  
  const int nv = vars.size();  

  
  if ( original ) 
    {
      // create ind x variant matrix in aux_skat;


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
	    {

// 	      std::cout << "maf , w = " 
// 			<< aux->maf[v] << " " 
// 			<< Helper::PROB::beta_pdf( aux->maf[v] , Aux_skat::a1 , Aux_skat::a2 ) << "\n";
	      
	      aux_skat->w[v] = Helper::PROB::beta_pdf( aux->maf[v] , Aux_skat::a1 , Aux_skat::a2 );
	    }
	    
	    Aux_skat::has_weights = true;
	}
      else if ( Aux_skat::has_weights ) 
	{
	  
	  meta_index_t midx = MetaInformation<VarMeta>::field( Aux_skat::weights );
	  
	  if ( midx.mt != META_FLOAT ) Helper::halt( "can only use Float variant tags as weights" );
	  
	  for (int v=0;v<nv;v++)
	    {
	      if ( vars(v).meta.hasField( Aux_skat::weights ) )
		aux_skat->w[v] = vars(v).meta.get1_double( Aux_skat::weights );
	      else if ( vars(v).consensus.meta.hasField( Aux_skat::weights ) )
		aux_skat->w[v] = vars(v).consensus.meta.get1_double( Aux_skat::weights );
	      else
		aux_skat->w[v] = 0;
	    }
	}
    }


  //
  // Permute residuals not the original 
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
  // Create kernel K
  //
  
  //  aux_skat->populate_K();

  
  //
  // Calculate Q 
  //

  Data::Matrix<double> W;

  double Q = aux_skat->calculate_Q( &W );
  
   
  //
  // Asympotic P-value
  //

  double pvalue = aux_skat->calculate_pvalue( Q , W );

 

  //
  // Output
  //

  if ( original ) 
    {
      
      if ( pvalue < 0 ) 
	(*output)["SKAT"] = "P=NA;";
      else
	(*output)["SKAT"] = "P=" + Helper::dbl2str( pvalue ) + ";";

      if ( aux_skat->logistic_model )
	{
	  std::map<std::string,int>::iterator i = aux->mc_a.begin();
	  while ( i != aux->mc_a.end() )
	    {
	      if ( i != aux->mc_a.begin() ) (*output)["SKAT"] += ";";
	      (*output)["SKAT"] += i->first + "(" + Helper::int2str( i->second ) + ")";
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

  get_liu_param( c1, c2, c3, c4 , 
		 &muX, &sigmaX, 
		 &muQ, &sigmaQ, 
		 &l, &d );
  
  double qnorm  = ( Q - muQ ) / sigmaQ;
  double qnorm1 = qnorm * sigmaX + muX;
  double pvalue = Statistics::noncentral_chi2_prob(  qnorm1 , l , d );

  return pvalue; 
}


void Pseq::Assoc::Aux_skat::get_liu_param( double c1 , double c2 , double c3 , double c4 , 
					   double * muX , double * sigmaX , 
					   double * muQ , double * sigmaQ , 
					   double * l , double * d ) 
{
  
  *muQ = c1;
  *sigmaQ = sqrt( 2 * c2 );

  double s1 = c3 / pow( c2 , 1.5 );
  double s2 = c4 * c2;
  
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
      a = 1.0/s1;
      *d = 0;
      *l = 1.0/(s1*s1);
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
 
  Data::Vector<double> eigen = Statistics::eigenvalues( K );

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
  bool okay = true;
  std::vector<int> df( lambda.size() , 1 );
  
  double pvalue_davies = d.pvalue( Q , lambda , df , &okay );

  // ignore davies ifault for now...
  okay = true;
  
  if ( lambda.size() == 1 || pvalue_davies <= 0 || pvalue_davies > 1 ) 
    okay = false;

  if ( okay ) return pvalue_davies; 
  return pvalue_liu;

}

