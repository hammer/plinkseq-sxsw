#ifndef __ALLELE_H__
#define __ALLELE_H__

#include <string>

#include "meta.h"


enum aType { ALLELE_UNKNOWN      = 0 , 
	     ALLELE_NOTSEEN      = 1 ,
	     ALLELE_SUBSTITUTION = 2 , 
	     ALLELE_INSERTION    = 3 , 
	     ALLELE_DELETION     = 4 }; 

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

    Allele(std::string aname, aType atype) 
	: aname(aname) , atype(atype)
	{    
	} 
    
    Allele(std::string aname)
	: aname(aname) 
	{    
	    // Null alleles
	    if ( aname == "." || aname == "N" ) 
		atype = ALLELE_NOTSEEN;

	    // Check for indels

	    else if ( aname.size() > 1 ) 
	    {
		if ( aname.substr(0,1) == "I" )
		    atype = ALLELE_INSERTION;
		else if ( aname.substr(0,1) == "D" )
		    atype = ALLELE_DELETION;
	    }
	    
	    else 
		atype = ALLELE_SUBSTITUTION;
	}
    

    // 
    // Simple type queries
    //
    
    bool indel() const 
	{ return atype == ALLELE_INSERTION || 
	      atype == ALLELE_DELETION; }
    
    bool snp() const 
	{ return atype == ALLELE_SUBSTITUTION; } 
    
    std::string name() const 
	{ return aname; }
    
    aType type() const 
	{ return atype; }

    //
    // Allow extensible annotation to be added per-allele
    //

    MetaInformation<AlleleMeta> meta;


 private:
    
    std::string     aname;
    
    aType           atype;
    
    

};


#endif
