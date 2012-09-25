#include "plinkseq/permute.h"
#include "plinkseq/indmap.h"
#include "plinkseq/phmap.h"
#include "plinkseq/individual.h"

#include <iostream>
#include <cmath>
#include <algorithm>

using namespace std;
using namespace Helper;

void Permute::initiate(int n, int s)
{

  if ( ! indmap ) 
    Helper::halt("no individual-map set in Permute::initiate()");
  
  clear();

  // the command --perm -1 implies adaptive permutation

  if ( n < 0 ) 
    adaptive();
         
  if ( n < 0 ) n = 0;
  if ( s < 0 ) s = 0;
  
  // number of permutations 

  nrep = adaptive_perm ? adaptive_max : n ;
  
  // number of test statistics per variant

  nstats = s;

  scores.resize( nstats );
  invalid.resize( nstats );
  mintie.resize( nstats );
  best_perm_score.resize( nstats );

  // Track max score per statistic per permuation [ 0..(nrep-1) ]
  // (i.e. for tie-count; for adaptive perm, this will mean we
  // have to allocate a large vector, but should be okay (can revisit
  // to build incrementally if needed)
  // note; no need to track for original

  max_score.resize( nstats );
  for (int i=0; i<nstats; i++) 
    max_score[i].resize( nrep , 0 );
  
  //
  // set any strata specified by phmap; also generally set up the
  // permutation list
  //

  set_strata();

}


void Permute::set_strata( const std::vector<bool> * fixed )
{
  
  //
  // Get total number of unique individuals in the map
  //

  int nind = indmap->size();
  
  // Set up initial permpos[] vector [ 0, 1, 2, ... , n-1 ]
  
  permpos.resize( nind );

  unpermute();
  
  strata.clear();
  
  //
  // Note: individuals designated as 'missing' do not get permuted. 
  //
  

  //
  // If no clustering variable specified, handle separately here
  //

  if ( ( ! phmap ) || (!phmap->strata_set() ) )
    {
      std::vector<int> t;
      for (int i=0; i<nind; i++)
	{	  
	  if ( indmap->ind( i )->included() )
	    {
	      if ( fixed ) 
		{
		  if (  ! (*fixed)[i] ) t.push_back( i ); 	  
		}
	      else 
		t.push_back(i);
	    }
	}
      
      strata.push_back( t );
      n_strata = 1;
      
      return;
    }
  

  //
  // Otherwise, first obtain a list of stratifying variables
  //
  
  std::map<int, std::vector<int> > strata_map;
  
  for (int i=0; i<nind; i++)
    {	  
      int k = -1;
      
      // if the person is missing, k will equal -1, which means
      // do not permute.  otherwise, classes numbered 0, 1, ...
      // (but not necessarily in order)

      if ( indmap->ind( i )->included() )
	k = indmap->ind( i )->group();
      
      if ( k != -1 ) 
	{
	  if ( fixed ) 
	    {
	      if ( ! (*fixed)[i] ) 
		strata_map[ k ].push_back( i );
	    }
	  else
	    strata_map[ k ].push_back( i );
	}
    }

  n_strata = strata_map.size();
  
  std::map<int , std::vector<int> >::iterator i = strata_map.begin();
  while ( i != strata_map.end() )
    {
      strata.push_back( i->second );
      ++i;
    }
    
}

// Reset original ordering

void Permute::unpermute()
{  
  for (int i=0; i<permpos.size(); i++) 
    {
      permpos[i] = i;
    }
}


void Permute::permute()
{

  //
  // Store remapped IDs
  //
  
  std::vector<std::vector<int> > i( n_strata );
  
  //
  // Permute phenotypes, within cluster
  //
  
  for (int k = 0; k < n_strata; k++)
    {
      std::vector<int> p( strata[k].size() );
      random_draw(p);
      i[k] = p ;
    }

  //
  // Post-permutation:
  // Iterate over clusters { s[][] }
  //

  // i[][] holds the permuted codes
  // s[][] points to individuals (non-missing)
  
  // Genotype =           sample[ s[j][k] ];        
  // Matching phenotype = sample[ s[ j ][ i[ j ][ k ] ] ];          

  // Create list of label-swapped codes

  for (int j=0; j < strata.size(); j++)
    {
      const std::vector<int> & y = strata[j];
      const std::vector<int> & z = i[j];
      for (int k=0; k < y.size(); k++)
	{
	  permpos[ y[k] ] = y[ z[k] ];	
	}
    }
}


void Permute::random_draw( vector<int> & a )
{
  
  // Generate a random permutation of 0 to n-1 where n is a.size(),
  // using Fisher-Yates shuffle.
  
  const int n = a.size( ) ;  
  for( int i = 0; i < n; i++ )
    a[ i ] = i;

  int tmp;
  for( int i = n; i > 1;  i-- )
    {
      int j = CRandom::rand(i);
      tmp = a[i-1];
      a[i-1] = a[j];
      a[j] = tmp;
    }
}


bool Permute::score( double s )
{
  vector<double> t(1);
  t[0] = s;
  return score( t );
}

bool Permute::score( const std::vector<double> & v )
{


  for (int s = 0; s < v.size(); s++)
    {
      
      bool valid = Helper::realnum( v[s] );
      
      //
      // An original observation?
      //
      
      if ( performed == 0 ) 
	{

	  //	  std::cout << "orig = " << v[s] << "\n";

	  int test = scores[s].size();
	  
	  if ( valid )
	    {
	      original_score[s] = v[s];
	      original_valid[s] = true;
	      
	      // Record original score
	      scores[s].push_back( v[s] );
	    }
	  else
	    {
	      invalid[s].insert( test );
	    }
	  
	}
      else // ... we are dealing with a score from a permuted dataset
	{
	  
	  if ( original_valid[s] ) 
	    {
	      if ( valid ) 
		{
		  

		  if ( v[s] > original_score[s] ) ++r[s];
		  else if ( v[s] == original_score[s] ) 
		    {
		      
		      // break exact ties; use whatever the last RND
		      // was, to avoid breaking the sequence of RNDs;
		      // also, this ensures that the same decision
		      // of how to break the tie is made similarly
		      // for different tests, if they give identical 
		      // results

		      // use different (default) RNG for the coin-flip
		      // to break ties -- i.e. so that sequence of 
		      // primary permutations remains preserved

		      if ( rand() / (double)RAND_MAX < 0.5 ) 
			{
			  ++r[s];
			}		      
		    }

		  
		  // track maximum test statistic across experiment
		  
		  if ( v[s] >= max_score[s][ performed - 1 ] )
		    {
		      max_score[s][ performed - 1 ] = v[s];
		    }
		  
		  //		  std::cout << "perm = " << v[s] << "\t" << mintie[s] << "\t" << best_perm_score[s] << "\n";

		  // track minimum p-value within site
		  if ( performed == 1 || v[s] > best_perm_score[s] )
		    {
		      mintie[s] = 1;
		      best_perm_score[s] = v[s];
		    }
		  else if ( v[s] == best_perm_score[s] )
		    {
		      mintie[s]++;
		    }	  
		  

		}
	      else
		{
		  // Count an invalid permuted score as greater than original
		  ++r[s];
		  ++n_invalid[s];
		}
	    }
	  
	}
    }
  
  
  // are we all done here? 
  
  // are we done yet? ( return T if more to do) 

  if ( adaptive_perm ) 
    {
      if ( performed % interval == 0 && adaptively_finished() ) return false;
    }
  else
    {
      if ( performed == nrep ) return false;
    }
  

  // If still here, more to do...

  ++performed;

  // permute before next round
  
  permute();
  
  return true;
}


bool Permute::finished() const
{
  return performed >= nrep;
}


//
// Functions to obtain/calculate empirical p-values, post-permutation
//

double Permute::pvalue(const int i) const
{
  return (double)(r[i]+1) / (double)( performed + 1 );
}

std::vector<double> Permute::pvalue() const
{
  vector<double> x( nstats, 0 );
  for (int i=0; i<nstats; i++) 
    x[i] = (double)(r[i]+1) / (double)(performed+1);
  return x;
}


// Minimum pseudo-p obtained under null
double Permute::min_pvalue(const int i) const
{
  return (double)(mintie[i]) / (double)(performed);
}

std::vector<double> Permute::min_pvalue() const
{
  vector<double> x( nstats, 0 );
  for (int i=0; i<nstats; i++) 
    x[i] = (double)(mintie[i]) / (double)(performed);
  return x;
}



bool Permute::valid_perm(const int t, const int j) const
{
  return invalid[t].find(j) != invalid[t].end();
}

vector<bool> Permute::valid_perm( const int j ) const
{
  vector<bool> b(nstats);
  for (int i=0; i<nstats; i++)
    b[i] = invalid[i].find(j) != invalid[i].end();
  return b;  
}

void Permute::calculate_max()
{
  for (int s=0; s<scores.size(); s++)
    for (int t=0; t<scores[s].size(); t++)
      { 
	double score = scores[s][t];
	int u = 1; // includes self;
	for (int p=0;p<performed; p++)
	  if ( max_score[s][p] >= score ) ++u;
	
	// replace original score with max(p)
	scores[s][t] = (double)u / (double)(performed+1);
      }  
  max_calculated = true;
}


double Permute::max_pvalue(const int s, const int t)
{
  if ( ! max_calculated ) calculate_max();
  return scores[s][t];
}


void Permute::fix( bool b )
{
  fix_samples = b;  
}


void Permute::fix( const std::vector<bool> & f )
{
  set_strata( &f );
}


bool Permute::adaptively_finished( )
{

  if ( performed < adaptive_min ) return false;
  if ( performed > adaptive_max ) return true;
  
  // update interval at which we check
  interval = (int)( adaptive_interval + performed * adaptive_interval2 );

  for (int s = 0 ; s < nstats ; s++)
    {
      // require at least one success
      if ( r[s] == 0 ) return false;
      
      // get current estimate of P , r+1 / n+1
      double pv = (double)( r[s] + 1 )/(double)( performed + 1 );
      
      double sd = sqrt( pv * (1-pv) / performed );
      double lower = pv - adaptive_zt * sd;
      double upper = pv + adaptive_zt * sd;
      if (lower<0) lower = 0;
      if (lower>1) upper = 1;
      
      // Is lower bound greater than threshold, or 
      // upper bound smaller than threshold?
      
      if ( ! ( lower > adaptive_alpha || upper < adaptive_alpha ) )
	return false;  // not done 
      
    }

  // is done
  return true;
}

  
