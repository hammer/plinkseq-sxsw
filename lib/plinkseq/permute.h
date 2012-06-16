#ifndef __PERM_H__
#define __PERM_H__

#include "crandom.h"
#include "statistics.h"

#include <vector>
#include <set>
#include <string>

class IndividualMap;
class PhenotypeMap;

class Permute { 
  
 private:
  
  // Seed, basic number of individuals, tests and replicates

  long int rseed;
  
  int performed;

  // Number of replicates
  int nrep;

  // Use adaptive permutation?
  bool adaptive_perm;
  int adaptive_min;
  int adaptive_max;
  int interval;
  int adaptive_interval;
  double adaptive_interval2;
  double adaptive_ci;
  double adaptive_zt;
  double adaptive_alpha;

  // per variant, # of statistics calculated
  int nstats; 

  // fix additional samples (above the original strat)
  bool fix_samples;
  
  std::vector<int> permpos;
  
  void unpermute();

  /////////////////////////////////////////////////////////////////////////
  
  // For a single test, number of success, number of replicates

  std::vector<double> original_score;
  
  // Was original score valid? 

  std::vector<bool> original_valid;
  
  // For single score, number of times permuted score >= original score
  
  std::vector<int> r;

  // Number of times the permuted data yielded an invalid result

  std::vector<int> n_invalid;

  // Number of times the minimum null statistic is tied

  std::vector<int> mintie;
  std::vector<double> best_perm_score;


  /////////////////////////////////////////////////////////////////////////

  
  // For experiment-wide corrected p-values, we need to track:
  //  a) the original score for all tests
  //  b) the maximum score across all tests for each replicate
  //  c) also track the tests with invalid results

  std::vector<std::vector<double> > scores;
  std::vector<std::vector<double> > max_score;
  std::vector<std::set<int>       > invalid;
  bool           max_calculated;



  /////////////////////////////////////////////////////////////////////////



  // Number of strata, and specification of strata if any exist

  int n_strata;

  std::vector< std::vector<int> > strata;
  
  
  // Helper functions
  
  void calculate_max();
  void random_draw( std::vector<int> & );
  void set_strata( const std::vector<bool> * fixed = NULL );
  
  // Maps

  IndividualMap * indmap;
  
  PhenotypeMap * phmap;
  
 public:
    
  Permute() 
    { 
      indmap = NULL;
      phmap = NULL;
      init();
    }
  
  Permute( IndividualMap & im , PhenotypeMap & pm )
    {
      set(im); set(pm);
      init();
    }
  
  void init()
  {
    rseed = time(0);
    nstats = 0;
    rseed = 0;
    nrep = 0;
    nstats = 0;
    n_strata = 1;
    fix_samples = false;
    max_calculated = false;
    reset();
    adaptive_perm = false;
    adaptive_interval  = 1;
    adaptive_interval2 = 0.001;
    adaptive_min = 5;
    adaptive_max = 1000000;
    adaptive_ci = 0.0001;
    adaptive_alpha = 0;
  }


  void seed(const long int s) { rseed = s; }

  void adaptive( const std::vector<int> & minmax )
  {
    if ( minmax.size() != 2 ) 
      {
	plog.warn("incorrect values for adaptive perm min/max");
	return;
      }
    adaptive_perm = true;
    adaptive_min = minmax[0];
    adaptive_max = minmax[1];
  }

  void adaptive( const int myadaptive_min = 5 , 
		 const int myadaptive_max = 1000000 , 
		 const int myadaptive_interval = 1 , 
		 const double myadaptive_interval2 = 0.001 ,
		 const double myadaptive_ci = 5 , 
		 const double myadaptive_alpha = 0 )
  {
    adaptive_perm = true;    
    adaptive_min = myadaptive_min;
    adaptive_max = myadaptive_max;
    adaptive_interval = myadaptive_interval;
    adaptive_interval2 = myadaptive_interval2;
    adaptive_ci = myadaptive_ci;
    adaptive_alpha = myadaptive_alpha ;
    interval = 1;
    // 20000 is fixed value for # of tests;     
    adaptive_zt = Statistics::ltqnorm( 1 - adaptive_ci / ( 2 * 20000 ) ) ; 

  }
  
  bool adaptively_finished( );
    
  // set to not permute, or unset, particular people (e.g. if missing
  // genotypes, not included in test statistic)

  void fix( bool b = true );
  void fix( const std::vector<bool> &); 
  
  void reset()
    {
      // reset class, ie. for new variant

      CRandom::srand( rseed );

      r.clear();
      r.resize( nstats , 0 );

      original_valid.clear();
      original_valid.resize(nstats,true);

      original_score.clear();
      original_score.resize( nstats, 0);

      n_invalid.clear();
      n_invalid.resize( nstats , 0 );
      
      mintie.clear();
      mintie.resize( nstats , 0 );

      best_perm_score.clear();
      best_perm_score.resize( nstats , 0 );

      performed = 0;
      unpermute();
      
      interval = 1; // adaptive perm
    }
 
  void clear()
    {
      // completely wipe
      rseed = 0;
      permpos.clear();
      nrep = 0;
      nstats = 0;
      
      r.clear();
      original_valid.clear();
      original_score.clear();
      n_invalid.clear();
      
      invalid.clear();
      mintie.clear();
      best_perm_score.clear();

      n_strata = 0;
      
      strata.clear();      
      scores.clear();
      max_score.clear();
      max_calculated = false;
      
      interval = 0; // adaptive perm
      adaptive_perm = false;
    }


  //
  // Some basic getter functions
  //

  long int seed() const { return rseed; }
  int replicates() const { return adaptive_perm ? adaptive_max : nrep; }  
  int n_stats() const { return nstats; }

  bool valid_original(const int i) const { return original_valid[i]; } 
  int n_tests(const int i) const { return scores[i].size(); }
  int n_valid_tests(const int i) const { return scores[i].size() - invalid[i].size(); } 
  int n_invalid_tests(const int i) const { return invalid[i].size(); } 
  int n_invalid_perms(const int i) const { return n_invalid[i]; } 

  // Convenience wrappers for a single statistic

  bool valid_orignal() const { return original_valid[0]; } 
  int n_tests() const { return scores[0].size(); }
  int n_valid_tests() const { return scores[0].size() - invalid[0].size(); } 
  int n_invalid_tests() const { return invalid[0].size(); } 
  int n_invalid_perms() const { return n_invalid[0]; } 



  //
  // Main setter function
  //

  void initiate(const int, const int s = 1);
  
  void set( IndividualMap & m ) { indmap = &m; }
  
  void set( PhenotypeMap & m ) { phmap = &m; }
  
  
  //
  // What is permuted index for obs 'i'?
  //

  int  pos(const int i) const { return permpos[i]; }


  //
  // Perform next permutation; return false if done
  //
  
  bool score( const std::vector<double> & );

  void permute();
  
  bool finished() const;
  
  // Convenience wrappers if only a single test specified

  bool score( double );
  
  bool valid_perm( const int , const int ) const;

  std::vector<bool> valid_perm( const int ) const;



  //
  // Obtain empirical p-values
  //
  
  std::vector<double> pvalue() const;
  std::vector<double> min_pvalue() const;
  std::vector<double> max_pvalue(const int);

  // Convenience wrappers if only a single test specified
  double pvalue( const int ) const;
  double min_pvalue( const int ) const;
  double max_pvalue(const int, const int);


  
};

#endif
