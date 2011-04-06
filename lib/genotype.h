
#ifndef __GENOTYPE_H__
#define __GENOTYPE_H__

#include <string>
#include <stdint.h>
#include <iostream>
#include <bitset>
#include <set>

#include "defs.h"
#include "individual.h"
#include "meta.h"


typedef std::bitset<8> bitset8;

enum gCodes { G_UNKNOWN  =  0, 
	      G_FAILED   =  1, 
	      G_PAT      =  2, 
	      G_MAT      =  3, 
	      G_PHASED   =  4, 
	      G_SWITCH   =  5,
	      G_HAPLOID  =  6,
	      G_MORE     =  7  };

// If HAPLOID, implies known original of allele
//   -> G_PAT/G_MAT indicates the original of allele
//   -> if 00 or 11, means not known 

// If PHASED, then PAT/MAT indicate 

/*!
  @class Genotype
  @brief Representation of an individual genotype
  
  A Variant object can contain 1 or more Genotype objects, 
  grouped together in a GenotypeSet object (the calls member
  of Variant).

*/


class Genotype { 
  
 private:
  
  // Which variant does this genotype belong to? 
  
  const Variant * parent;
  
  // Header byte 
  
  bitset8 bs;
  
  // G_UNKNOWN Data not available or not genotyped (T/F)
  // G_FAILED  Genotyping performed but failed (T/F)
  // G_PAT     Paternal allele (ref/alternate)
  // G_MAT     Maternal allele (ref/alternate)
  // G_PHASED  Phase-known genotype (T/F)
  // G_SWITCH  Likely switch error (T/F)  
  // G_HAPLOID Haploid (T/F)
  // G_MORE    Has integer genotype code following (T/F)
  
  
  // Core genotype store, a single integer; for a basic, non-missing
  // SNP, this will simply be 0/1/2, i.e. the number of non-reference
  // alleles
  
  int gcode;
  

  
  //////////////////////////////////////////////////////////////////////
  
  void initialise()
    {
      bs.reset();
      bs.set( G_UNKNOWN , true);  // is missing unless otherwise stated
      gcode = -1;
    }
  
 
 public:
  
  Genotype( Variant * parent = NULL ) 
    : parent(parent) 
    { 
      initialise(); 
    }
    
    Genotype( const Genotype & g )
      {
	*this = g;
      }
    
    Genotype & operator=( const Genotype & g ) 
      {
	this->bs = g.bs;
	this->parent = g.parent;
	gcode = g.gcode;
	meta = g.meta;
	return *this;
      }
    
    
    void variant( const Variant * var )
    {
      parent = var;
    }
    
    const Variant * variant()
    {
      return parent;
    }

  /// Extensible meta-information, e.g. read-depth, soft-calls, etc
  
  MetaInformation<GenMeta> meta;
  
  /// Set missing status of genotype
  void missing(bool b)  { bs.set(G_UNKNOWN,b); }

  /// Get missing statis of genotype
  bool missing() const { return bs.test( G_UNKNOWN ); } 
  
  /// Set failed status of genotype
  void failed(bool b)  { bs.set(G_FAILED,b); }

  /// Get failed status of genotype
  bool failed() const { return bs.test( G_FAILED ); } 

  /// Test for either a missing (not attempted) or failed (but attempted) genotype
  bool null() const  { return bs.test( G_UNKNOWN) || bs.test(G_FAILED); }
  bool notnull() const  { return ! ( bs.test( G_UNKNOWN ) || bs.test(G_FAILED) ) ; }

  bool null(const bool b)  { missing(b); failed(b); }

  /// Convenience function to set paternal allele, if a simple SNP
  void pat(bool b)  { bs.set(G_PAT,b); }
  
  /// Convenience function to get paternal allele, assuming a simple SNP
  bool pat() const { return bs.test( G_PAT ); }


  void mat(bool b)  { bs.set(G_MAT,b); }
  
  bool mat() const { return bs.test( G_MAT ); }

  /// Set genotype (assuming a simple SNP)
  void set_alternate_allele_count(const int g)
  {
    null(false);
    if ( g == 0 ) { pat(false); mat(false); }
    else if ( g == 1 ) { pat(true); mat(false); phased(false); }
    else if ( g == 2 ) { pat(true); mat(true); }
    else failed(true); // invalid code
  }

  void haploid(bool b)  {  bs.set( G_HAPLOID, b); }
  bool haploid() const  { return bs.test( G_HAPLOID); }

  void diploid(bool b)  {  bs.set( G_HAPLOID, !b); }
  bool diploid() const  { return ! bs.test( G_HAPLOID); }

  void phased(bool b)  { bs.set(G_PHASED,b); }
  bool phased() const  { return bs.test( G_PHASED ); }

  void pswitch(bool b)  { bs.set(G_SWITCH,b); }
  bool pswitch() const { return bs.test( G_SWITCH ); }

  void more(bool b)  { bs.set(G_MORE,b); }
  bool more() const { return bs.test( G_MORE ); }


  // BCF encoding 
  
  uint8_t bcf() const;
  void bcf( uint8_t );

  /*!
    @brief Genotype value lookups
    @return 
   */


  //
  // Primary genotype value lookups
  //

  // number of nonreference alleles 
  int allele_count( ) const;
  
  // number of a particular non-refernece allele (number 1,2,...)
  // 0 implies reference
  int allele_count( const int altcode ) const;
  int allele_count( const std::string & ) const;
  
  
  int copy_number() const 
    {
      if ( bs.test( G_HAPLOID ) ) return 1;
      else return 2;
    }
  
  bool nonreference() const
  {
    return (!more()) && (!null()) && ( pat() || mat() ) ;
  }

  bool heterozygote() const
  {
    return (!more()) && (!null()) && ( pat() || mat() ) && ! ( pat() && mat() ) ;
  }
  
  bool alternate_homozyote() const
  {
    return (!more()) && (!null()) && ( pat() && mat() ) ;
  }
  
  bool reference() const
  {
    return (!more()) && (!null()) && ! ( pat() || mat() ) ;
  }
  
  bool minor_allele(const bool altmin) const
    {
      return altmin ? nonreference() : minor_allele_count(altmin) > 0 ;
    }

  int minor_allele_count(const bool altmin) const
  {
    // ignore multi-allelic, missing, CN!=2, etc, for now...
    return altmin ? (int)pat() + (int)mat() : 2 - (int)pat() - (int)mat();
  }

  std::vector<int> allele_list() const;


  //
  // Conversion function
  //

  char pack() const { return bs.to_ulong(); }

  void unpack(char c) { bs = c; } 

  /*! 
    @brief Get genotype code 
    @return Integer that indicates genotype class for individual    
  */
  
  int code() const { return gcode; }
  
  /*!
    @brief Set non-simple genotype code 
  */

  void code(int g) { gcode = g; }

  static bool equivalent( const Genotype &, const Genotype & );
		 
  bool operator==(const Genotype & other) const 
    {
      return bs == other.bs && gcode == other.gcode;
    }
  
  bool operator!=(const Genotype & other) const 
    {
      return !(*this == other);
    }



  // Display functions


  std::string print_vcf() const
    {
      if ( more() ) return "-1/-1";       
      std::string s1, s2;
      if ( null() ) 
	s1 = s2 = ".";
      else
	{
	  s1 = pat() ? "1" : "0";
	  s2 = mat() ? "1" : "0";
	}      
      if ( haploid() ) return s1;      
      if ( pswitch() ) return s1 + "\\" + s2;
      else if ( phased() ) s1 + "|" + s2;
      else return s1 + "/" + s2;
    }

  friend std::ostream & operator<<( std::ostream & out, const Genotype & g)
  { 
    out << g.missing() << g.failed() << " - "; 
    
    out << g.pat();
    if ( !g.haploid() )
      {
	if ( g.pswitch() ) out << "\\";
	else if ( g.phased() ) out << "|";
	else out << "/";	
	out << g.mat();
      }
    
    out << " - "
	<< g.phased() 
	<< g.pswitch() 
	<< g.haploid()
	<< g.more();
    
    if ( g.more() ) 
      out << " " << g.gcode;
    
    return out;
  }
  
};


/*!
  @class GenotypeSet
  
  A collection of Genotype objects; each Variant has a GenotypeSet
  object that contains the calls for the individuals in the VCF
*/

struct SampleVariant;

class GenotypeSet{

 private:
  
    std::vector<Genotype> calls;

    // If this is non-null, it means that the alignment was flat and
    // so the genotypes are stored elsewhere NULL means either this is
    // the consensus, or it contains calls (unflat alignment) If not
    // null, we assume calls will be empty ( by definition) We also
    // assume we do not add/remove variants to a mirroring set --
    // i.e. the mirroring is designed for read-operations only

    SampleVariant * svar;

    std::vector<int> * incon;

 public:
    
    GenotypeSet( SampleVariant * p = NULL ) { svar = p; } 
     
    /// Return number of individual genotype calls in set
    int size() const;

    /*!
      @brief Allocate space for calls
      @param n Number of individuals      
    */
    
    void size(int n) { if ( ! svar ) calls.resize(n); }
    
    void set_consensus_slotmap( SampleVariant * ps , std::vector<int> * pm )
    { 
      svar = ps;
      incon = pm; 
    } 

    /// Add a Genotype to the GenotypeSet
    void add(Genotype & g) { if ( ! svar ) calls.push_back(g); }
    
    /// Add a Genotype at a particular individual index i
    void add(Genotype & g, int i) { if ( ! svar ) calls[i] = g; }

    /// Return Genotype for individual i
    Genotype & genotype(int i);

    /// Return const Genotype for individual i
    const Genotype & genotype(int i) const;
    
    /// Return genotype code an individual    
    int code(int i);
    
    /// Wipe all genotype calls (only for actual calls)
    void clear() { if ( ! svar ) calls.clear(); }

    void summarise_meta( std::map<meta_typed_key_t,std::pair<int,int> > & , 
		         std::map<meta_typed_key_t,std::string>  & ,
		         std::set<meta_typed_key_t> & ) const ;

    

  
};

#endif
