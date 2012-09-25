#ifndef __ALLELE_H__
#define __ALLELE_H__

#include <string>

#include "plinkseq/meta.h"


enum aType { ALLELE_UNKNOWN      = 0 , 
	     ALLELE_REFERENCE    = 1 ,
	     ALLELE_SUBSTITUTION = 2 ,
	     ALLELE_INSERTION    = 3 ,
	     ALLELE_DELETION     = 4 ,
	     ALLELE_MULTINUCLEOTIDE_SUBSTITUTION = 5 }; 

// From VCF specification:


// REF reference base. One of A,C,G,T,N

// ALT comma separated list of alternate non-reference alleles called on at
// least one of the samples. Options are A,C,G,T,Dn (for delete n bases
// starting with the base at POS), I<seq> where <seq> is a list of ACGT bases
// to be inserted just after the base at POS, '.' (period character) if there
// are not alternate alleles.


class Allele {
  
 public:
  
  Allele() 
    { 
      atype = ALLELE_UNKNOWN;
      aname = "";
    } 
  
 Allele( const std::string & aname, aType atype) 
   : aname(aname) , atype(atype)
  {    
  } 
  
 Allele( const std::string & aname , const std::string & reference  )
   : aname(aname) 
  {    

    // allele type defined w.r.t reference 
    
    if ( aname == reference ) 
      atype = ALLELE_REFERENCE;

    // Null alleles    
    else if ( aname == "" || aname == "." || aname == "N" ) 
      atype = ALLELE_UNKNOWN; 
    
    // basic SNP
    else if ( reference.size() == 1 && aname.size() == 1 ) 
      atype = ALLELE_SUBSTITUTION;
    
    // Check for indels, MNPs     
    else if ( aname.substr(0,1) == "I" || aname == "<INS>" )
      atype = ALLELE_INSERTION;
    
    else if ( aname.substr(0,1) == "D" || aname == "<DEL>" )
      atype = ALLELE_DELETION;
        
    else if ( reference.size() < aname.size() )
      atype = ALLELE_INSERTION;

    else if ( reference.size() > aname.size() ) 
      atype = ALLELE_DELETION;

    else 
      atype = ALLELE_MULTINUCLEOTIDE_SUBSTITUTION;
  }
    
  
  // 
  // Simple type queries
  //
  
  bool indel() const { return atype == ALLELE_INSERTION || atype == ALLELE_DELETION; }

  bool ins() const { return atype == ALLELE_INSERTION; }

  bool del() const { return atype == ALLELE_DELETION; }
  
  bool snp() const { return atype == ALLELE_SUBSTITUTION; } 
  
  bool mnp() const { return atype == ALLELE_MULTINUCLEOTIDE_SUBSTITUTION; } 

  std::string name() const { return aname; }
  
  aType type() const { return atype; }
  
  //
  // Allow extensible annotation to be added per-allele
  //
  
  MetaInformation<AlleleMeta> meta;
  
  
 private:
  
  std::string     aname;
  
  aType           atype;
  
  
};


#endif
