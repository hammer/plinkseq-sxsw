#ifndef __PSEQ_EM_H__
#define __PSEQ_EM_H__

#include <vector>
#include "meta.h"

class Variant;

class EM {
  
 public:

  EM() : EPS(1e-5) , maxiter(10) 
    {
      // Ensure the we have a genotype meta-tag for posteriors
      MetaInformation<GenMeta>::field( PLINKSeq::META_GENO_POSTPROB(), 
				       META_FLOAT ,
				       -1 , 
				       "P(genotype|read data)");
    }
  
  /// Load EM with genotype likelihoods (AA, AB, BB)
  void load( Variant & );

  /// Specify group membership for individual i
  void group( const int i , const int g );
  
  /// Run EM, return iteration code (-1 failed)
  int estimate();
  
  /// Call genotypes and place in a new SampleVariant
  void call( bool, const double t) const;

  /// For individual i, return posteriors
  std::vector<double> posteriors(const int i) const;

  /// For group g, return allele frequency
  double frequency(const int g) const;

  /// Return frequency estimate
  double frequency() const;

  /// Return mean entropy of posteriors (all,alternate-genotype)
  void entropy( double & h , double & halt ) const;
  
  /// Return mean max
  double mean_max_posterior() const;


 private:
  
  double EPS;
  double maxiter;

  
  Variant * var;

  // number of individuals 
  int n ;

  // Helper functions
  std::vector<double> lik_to_probs( std::vector<double> & , bool ) const;
  
  // Genotype likelihoods, P(reads|G)
  std::vector< std::vector<double> > gl;

  // Allele frequency, -> P(G)
  double f;

  // Posteriors, P(G|reads)
  std::vector< std::vector<double> > post;
  
};

#endif
