#include "pseq.h"
#include "genic.h"

#include "../lib/prob.h"

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
Data::Vector<double> Pseq::Assoc::Aux_skat::u; // phenotype --> modified (u_i)
Data::Matrix<double> Pseq::Assoc::Aux_skat::X; // covariates
std::vector<bool> Pseq::Assoc::Aux_skat::mask; // inclusion mask
int Pseq::Assoc::Aux_skat::n_actual = 0;       // actual # of individuals

bool Pseq::Assoc::Aux_skat::logistic_model = true; // otherwise linear model for phenotype

bool Pseq::Assoc::Aux_skat::has_covar = false;
std::vector<std::string> Pseq::Assoc::Aux_skat::covars;
bool Pseq::Assoc::Aux_skat::has_weights = false;
std::string Pseq::Assoc::Aux_skat::weights = "";

bool Pseq::Assoc::Aux_skat::use_freq_weights = false;
int Pseq::Assoc::Aux_skat::a1 = 1;
int Pseq::Assoc::Aux_skat::a2 = 25;

void Pseq::Assoc::Aux_skat::fit_null()
{

  std::cout << "fitting null model...\n";

  const int n = g.indmap.size();
  
  // Get # of people with non-missing phenotypes or covariates
  
  y.clear();
  X.clear();
  
  mask.resize( n , false );
  
  // Covariates?
  Data::Matrix<double> C;
  if ( has_covar ) C = g.phmap.covariates( covars , g.indmap );
  
  
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
  
  if ( ! glm.valid() ) Helper::halt( "problem attempting to fit the null model for SKAT" );

  // Get vector of coefficients
  
  Data::Vector<double> beta;
  if ( ! glm.display(&beta) ) 
    Helper::halt( "problem fitting null model for SKAT" );


  // Calculate residuals (now in a reduce format of only the 
  Data::Vector<double> o = y;
  y.resize( n_actual );
  j=0;
  for (int i=0;i<n;i++)
    {
      if ( ! mask[i] ) 
	{
	  y[j] = beta[0];
	  for (int c=0;c<ncov;c++) y[j] += C(i,c) * beta[c+1];
	  
	  // Y=logit^-1(Y)
	  if ( logistic_model ) 
	    {
	      double e = exp( y[j] );
	      y[j] = e / ( 1 + e ) ;
	    }

	  // make y-u and store
	  u[i] = y[i];    
	  y[i] = o[i] - y[i];
	  ++j;
	}
    }

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
		G(j,v) = 2 * aux->maf[v];
	      else 
		G(j,v) = genotype.minor_allele_count( aux->altmin[v] ) ;
	      
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
  
double Pseq::Assoc::Aux_skat::calculate_Q()
{
  double Q = 0;
  const int nv = G.dim2();
  
  // Q = sum_v w_v S_v^2
  // where S_v = g' ( y - u ) 
  
  // for each variant:
  for (int v=0;v<nv;v++)
    {
      double t = 0;
      for (int i=0;i<n_actual;i++)
	t += G(i,v) * y[i] ;  
      Q += w[v] * t*t;
    }

  // Establish p-value for Q
  // P0  = V - VX(X'VX)^-1X'V

  //   where X = [ 1 , X ] matrix (X)
  //   V = sigma^2 I 
  //   where sigma is an estimator of under the null model
  
  // for dichot: 
  //   V = diag( l_i ) 
  //  where l_i 

  double sigmaSQR = 1; // how is this calculated?

  Data::Matrix<double> V(n_actual, n_actual);
  if ( logistic_model ) 
    for (int i=0;i<n_actual;i++) V(i,i) = u[i] * ( 1 - u[i] );
  else
    for (int i=0;i<n_actual;i++) V(i,i) = sigmaSQR;

  Data::Matrix<double> P0 = V - V * X * Statistics::inverse( Statistics::transpose(X) * V * X ) * Statistics::transpose(X) * V;  

//   Data::Vector<double> Statistics::eigenvalues( Data::Matrix<double> & a );
//   Data::Matrix<double> Statistics::inverse( const Data::Matrix<double> & u_orig, bool * flag );
//   Data::Matrix<double> Statistics::transpose( const Data::Matrix<double> & d );
  
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
      
      std::cout << "DONE\n";

      //
      // Create weights?
      //

      aux_skat->w.resize( nv , 1 );
      
      if ( Aux_skat::use_freq_weights ) 
	{

	  // MAF (from whole sample) should be in aux_prelim
	  std::cout << "setting freq weights\n";
	    for (int v=0;v<nv;v++)
	      {
		std::cout << "maf , w = " 
			  << aux->maf[v] << " " 
			  << Helper::PROB::beta_pdf( aux->maf[v] , Aux_skat::a1 , Aux_skat::a2 ) << "\n";

		aux_skat->w[v] = Helper::PROB::beta_pdf( aux->maf[v] , Aux_skat::a1 , Aux_skat::a2 );
	      }
	}
      else if ( Aux_skat::has_weights ) 
	{
	  Helper::halt("user weights not yet implemented");
	}

    }
  
  std::cout << "DONE WEIGHTS\n";

  
  //
  // Populate G matrix
  //
  
  aux_skat->populate_G( vars , aux );

  std::cout << "DONE G\n";
  
  //
  // Create kernel K
  //
  
  aux_skat->populate_K();

  std::cout << "DONE K\n";  
  
  //
  // Calculat Q 
  //

  double Q = aux_skat->calculate_Q();

  std::cout << "returning original " << Pseq::Assoc::Aux_skat::n_actual << " " << nv << " ==> " << Q << "\t" << original << "\n";

  return 0;

}



// # n subjects
// # p variants 
// # y_i phenotype
// # X_ij covariates
// # G_ij genotypes (allelic coding)

// # classical linear/logistic regression

// # SKAT is like fitting all betas for the G set of varinats, except it
// # assumes that they follow an arbitrary distribution with mean 0
// # and variance w_j . tau ; where w_j  weight for variant_j, tau variance component

// # thus testing H_0 : B = 0  equivalent to testing tau = 0
// # --> frame as variance-component score test in mixed model

// # Only requires fitting Y_i = u + alpha.X + e   (once)

// # Variance component score statistic

// #  Q = ( y - mean(y) )' K ( y - mean(y) )  

// # where K = GWG' 

// # mean(y) is the predicted mean of y under H0 (i.e. will be indiv spec.) 
// # mean(y) = u-hat + X * alpha-hat
// #         = logit^{-1}( .. )

// # is weight matrix for 'p' variants (diagonal matrix)
// # K is n-by-n matrix (kernel) 
// #  K(i,i') = sum_j=1/j=p w_j . G_ij G_i'j . 

// # this means K() is 'weighted linear kernel function'

// # propse w_j = Beta( MAF_j , a1, a2 ) 
// #   propose a1, a2 fixed;  suggets a1=1, a2=25 for rare variants
// # (this increases weight of rare variants, but allows non-zero weight 
// #  on moderate 1-5% MAF variants

// # smaller a1 means more strongly increasing weight of rare variants

// # a1=a2=1 means w_j =1  and so all variants are weighed equally

// # a1=a2 = 0.5 means  sqrt(w_j) = 1 / sqrt( MAF_j ( 1 - MAF_j ) ) 
// # i.e. inverse of variance ;  this puts almost 0 weight on all MAF>1%
// # and so is good for truely rare variant models

// # Q follows a mixture of chi-sqare
// # Approximated via efficient Davies method (appendex A)

// # Note -- is w_j =1 , outcom dichot, and no covariates, should equal C-alpha

// # Q is a weighted sum of individual variant tests
// # 
   
//    // Davies method to approximate distribution of Q
 

   
