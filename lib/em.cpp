
#include "em.h"
#include "helper.h"
#include "variant.h"

#include <iostream>

#include <cmath>

std::vector<double> EM::lik_to_probs( std::vector<double> & g , bool phred_scaled ) const
{
  
  std::vector<double> p(3,-1);
  
  // for now, only handle biallelic GLs/PLs
  if ( g.size() != 3 ) return p;

  double g0, g1, g2;

  if ( phred_scaled ) 
    {
      g0 = g[0] == 0 ? 1 : pow( 10 , -g[0] / 10.0 );
      g1 = g[1] == 0 ? 1 : pow( 10 , -g[1] / 10.0 );
      g2 = g[2] == 0 ? 1 : pow( 10 , -g[2] / 10.0 );
    }
  else
    {            
      // assumes GL are in order: RR,AR,AA
      // for A=alternate, R=reference (as per current GATK)
      // Assumes GL are log10(p)
      
      g0 = g[0] == -0 ? 1 : pow( 10 , g[0] );
      g1 = g[1] == -0 ? 1 : pow( 10 , g[1] );
      g2 = g[2] == -0 ? 1 : pow( 10 , g[2] );
      
    }
  
  if ( ! Helper::realnum( g0 ) ) g0 = 0;
  if ( ! Helper::realnum( g1 ) ) g1 = 0;
  if ( ! Helper::realnum( g2 ) ) g2 = 0;
  
  // scale probabilities to sum to 1.0
  
  double s = 1.0 / ( g0 + g1 + g2 );
  g0 *= s;  g1 *= s; g2 *= s;

  if ( ! ( Helper::realnum( g0 ) 
	   && Helper::realnum( g1 ) 
	   && Helper::realnum( g2 ) ) )
    {
      return p;
    }
  else
    {
      p[0] = g0;
      p[1] = g1;
      p[2] = g2;
    }

  return p;
}


/// Load EM with genotype likelihoods (AA, AB, BB)

void EM::load( Variant & v )
{
  
  var = &v;
  
  n = v.size();
  
  gl.resize( n );
  post.resize( n );
  
  for (int i=0; i<n; i++)
    {
      
      std::vector<double> g;

      bool phred = true;
      
      if ( v(i).meta.hasField( PLINKSeq::META_GENO_PHRED() ) )
	{
	  std::vector<int> tmp = v(i).meta.get_int( PLINKSeq::META_GENO_PHRED() );
	  g.resize( tmp.size() );
	  for (int t=0; t<g.size(); t++) g[t] = tmp[t];
	}
      else if ( v(i).meta.hasField( PLINKSeq::META_GENO_LIK() ) )
	{	  
	  phred = false;
	  g = v(i).meta.get_double( PLINKSeq::META_GENO_LIK() );
	}
      
      std::vector<double> p = lik_to_probs( g , phred );
      
      gl[i] = p;      
      post[i] = p; // also set P(G|reads) to 3-element vector
    }
}


/// Specify group membership for individual i
void EM::group( const int i , const int g )
{
  //
}


/// Run EM, return iteration code (-1 failed)
int EM::estimate()
{


  
  // starting value for p 

  f = 1.0 / ( (double)2*n ) ;
  
  int iter = 0;

  while ( 1 ) 
    {
      
      //   P(G) as function of 'p', assuming HWE
  
      double fhom = f*f;
      double fhet = 2*f*(1-f);
      double fref = 1-fhom-fhet;
  
      // Estimate P(G|reads)
      
      for (int i=0; i<n; i++)
	{
	  
	  std::vector<double> & g = gl[i];

	  double pref = g[0] * fref;
	  double phet = g[1] * fhet;
	  double phom = g[2] * fhom;
	  
	  double denom = 1.0 / ( phom + phet + pref ) ;
	  
	  std::vector<double> & p = post[i];      
	  p[0] = pref * denom;
	  p[1] = phet * denom;
	  p[2] = phom * denom;
	  
	}
      

      // Now estimate allele frequency
      
      double of = f;
      
      f = 0;
      
      for (int i=0; i<n; i++)
	{            
	  std::vector<double> & p = post[i];      	  
	  f += p[2] * 2 + p[1];
	}
      
      f /= (double)(2*n);

      if ( f <= 0 ) 
	{
	  f = 0;
	  break;
	}
      if ( f >= 1 ) 
	{
	  f = 1;
	  break;
	}
      
      if ( ++iter == maxiter ) break;

      if ( fabs( f - of ) < EPS ) break;

    }

  std::cout << "iter = " << iter << "\n";

}


/// For individual i, return posteriors

std::vector<double> EM::posteriors(const int i) const
{
  return post[i];
}

double EM::frequency() const
{
  return f;
}


/// For group g, return allele frequency
double EM::frequency(const int g) const
{
  return 0;
}


void EM::call( const double t ) const
{

  const int n = var->size();
  
  for (int i=0; i<n; i++)
    {
      
      Genotype & g = (*var)(i);
      
      const std::vector<double> & p = post[i];

      int m = 0;

      if ( p[1] > p[0] ) 
	{
	  m = p[2] > p[1] ? 2 : 1 ;
	}
      else
	{
	  m = p[2] > p[0] ? 2 : 0 ;
	}
      

      //
      // Does the maximum PP meet threshold?
      //
      
      if ( p[m] >= t ) 
	{
	  g.set_alternate_allele_count(m);
	}
      else // set to missing
	{
	  g.null( true );
	}

      
      //
      // Insert PP as genotype posterior probabilities, and REF/ALT allele dosages
      // TODO: assume biallelic for now... need to fix
      //

      g.meta.set( PLINKSeq::META_GENO_POSTPROB() , p );
      
      g.meta.set( PLINKSeq::META_GENO_ALT_DOSAGE() , p[1] + 2 * p[2] );
      
    }


}


void EM::entropy( double & e , double & a ) const
{
  
  // e is entropy for all individuals
  // a is entropy only considering individuals where an alternate allele is more likely
  //   than the reference

  e = a = 0;

  int na = 0; // count for 'a'

  for (int i=0; i<post.size(); i++)
    {
      
      const std::vector<double> & p = post[i];
      
      double t = 0;      
      if ( p[0] > 0 ) t -= p[0] * log( p[0] );
      if ( p[1] > 0 ) t -= p[1] * log( p[1] );
      if ( p[2] > 0 ) t -= p[2] * log( p[2] );      
      
      // Sum over all individuals
      e +=t;
      
      // Is non-reference more likely?
      
      if ( p[1] > p[0] || p[2] > p[0] )
	{
	  a += t;
	  ++na;
	}
    }

  // Return mean entropies
  e /= (double)post.size();
  a /= (double)na;
}

double EM::mean_max_posterior() const
{

  double mx = 0;

  for (int i=0; i<post.size(); i++)
    {

      const std::vector<double> & p = post[i];

      int m = 0;
      if ( p[1] > p[0] ) 
	{
	  if ( p[2] > p[1] ) m=2; else m=1;
	}
      else
	{
	  if ( p[2] > p[0] ) m=2; else m=0;
	}
      
      mx += p[m];
    }
  
  return mx/(double)n;
}


