#ifndef __PSEQ_POPS_H__
#define __PSEQ_POPS_H__

class Variant;
class SeqDBase;

#include <vector>
#include <string>
#include "plinkseq/regions.h"

void f_ipop( Variant & v , void * p );


namespace Pseq {
  
  struct pop_t 
  {
    std::string label;
    std::vector<double> p; 
    double prior;
  };
  
  struct var_t 
  {
  var_t( const std::string & a , const std::string & b ) : a1(a) , a2(b) {} 
    std::string a1;
    std::string a2;
  };
  
  struct IPop 
  {
    
    IPop( const std::string & filename , SeqDBase * seqdb );

    void calculate();

    bool has( const Region & r )
    {
      return positions.find( r ) != positions.end() ; 
    }

    // population freqs/labels
    std::vector<pop_t> populations;
    
    // allele codes
    std::vector<var_t> alleles;
    
    // positions
    std::map<Region,int> positions;
    
    // accumulate genotype data ( indiv / variant )
    std::vector<std::vector<int8_t> > genotypes;
    
    // skip these 'bad' markers from calculation
    std::vector<bool> bad;
    
    // Needs a SEQDB to be attaced
    SeqDBase * seqdb;
    
  };
  
}

#endif
