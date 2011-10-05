
#ifndef __GENOTYPE_H__
#define __GENOTYPE_H__

#include <string>
#include <stdint.h>
#include <iostream>
#include <set>

#include "defs.h"
#include "individual.h"
#include "meta.h"


/*!
  @class Genotype
  @brief Representation of an individual genotype
  
  A Variant object can contain 1 or more Genotype objects, 
  grouped together in a GenotypeSet object (the calls member
  of Variant).

*/


class Genotype { 
  
 private:

  // Core entities
  
  int       ploidy;  // 0, 1 or 2

  uint8_t   allele1; // 255 == missing
  
  uint8_t   allele2; // 255 == missing
  
  bool      is_null;
  
  bool      known_phase;
  


  //
  // Helpers: caching of known genotypes (numeric encoding)
  //

  static std::map<std::string,Genotype> gcache;

  static void clear_genotype_cache();
  
  static const Genotype * search_genotype_cache( const std::string & tok );


  //
  // Private member functions
  //

  /// Set to null and wipe all values
  void set_null() 
  {
    ploidy = 0;
    allele1 = allele2 = 0;      
    is_null = true;
    known_phase = false;     
  }


 public:
  
  static genotype_model_t model;
  
  MetaInformation<GenMeta> meta;

  
  Genotype( ) 
    { 
      set_null();
    }

  Genotype( const std::string & str , const int );
  
  Genotype( const std::string & str , 
	    const int gt_field , 
	    const std::vector<meta_index_t*> & formats , 
	    const int n_alleles ) ;
  

  /// Convert a numeric encoding of a genotype to a Genotype (w/ caching)  
  void set_from_string( const std::string & str , const int n_alleles );
  
  /// Set to null (but preserve other values)
  void null(const bool b)  { is_null = b; } 

  /// Is this genotype null?
  bool null() const { return is_null; } 

  // e.g. for PED printing, is this genotype consistent with a biallelic variant 
  //      (restricted to 0,1 ref,alt alleles only)

  /// Is this genotype 'complex' or restricted to 0/1 alleles (e.g. for PED file printing)
  bool more() const 
  { 
    if ( is_null ) return false;
    if ( allele1 > 1 ) return true;
    if ( ploidy != 2 ) return true;
    return allele2 > 1;
  }

  /// Set a diploid genotype given two allele codes
  void genotype( const uint8_t a1 , const uint8_t a2 ) 
  {
    ploidy = 2;
    allele1 = a1;
    allele2 = a2;    
    is_null = false;
    known_phase = false;
  }

  /// Set a phased diploid genotype given two allele codes
  void phased_genotype( const uint8_t a1 , const uint8_t a2 ) 
  {
    genotype(a1,a2);
    known_phase = true;
  }

  /// Set a haploid genotype given one allele code
  void genotype( const uint8_t a1 ) 
  {
    ploidy = 1 ;
    allele1 = a1;
    is_null = false;
    known_phase = false;
  }

  /// Convert diploid to haploid call, return F is problem
  bool make_haploid()
  {
    if ( ploidy == 1 ) return true;
    ploidy = 1 ;
    known_phase = false;    
    if ( is_null ) return true;
    if ( allele1 != allele2 ) 
      {
	is_null = true;
	return false;
      }
    return true;
  }
  
  /// Set to diploid genotype w/ 0-2 alt-alleles (assumes basic SNP)
  void set_alternate_allele_count(const int g)
  {    
    if      ( g == 0 ) genotype(0,0);
    else if ( g == 1 ) genotype(0,1);
    else if ( g == 2 ) genotype(1,1);
    else is_null = true;
  }

  void haploid(const bool b)  { ploidy = 1; }
  bool haploid() const  { return ploidy == 1; } 
  
  void diploid(const bool b)  { ploidy = 2; }
  bool diploid() const  { return ploidy != 2; }
  
  void phased(const bool b)  { known_phase = b; }
  bool phased() const  { return known_phase; }
  

  //
  // Primary genotype value lookups
  //
  
  // number of alternate alleles
  int allele_count( ) const 
  {
    if ( is_null || ploidy == 0 ) return 0;
    if ( ploidy == 1 ) return allele1 != 0;    
    int ac = 0;
    if ( allele1 ) ++ac;
    if ( allele2 ) ++ac;
    return ac;
  }
  
  // number of observed alleles
  int copy_number( ) const
  { 
    return is_null ? 0 : ploidy;
  }
  
  // number of a particular allele, numeric coding (0=ref, 1,2,3,...)
  int allele_count( const int a ) const
  {
    if ( is_null || ploidy == 0 ) return 0;
    if ( ploidy == 1 ) return allele1 == a;
    return allele1 == a + allele2 == a;    
  }

  // genotype scoring function
  double score( genotype_model_t model = GENOTYPE_MODEL_UNSPEC );

  // Used when recalling a genotype (i.e. merging SampleVariants with 
  // different alt-alleles

  int acode1() const { return allele1; } 
  void acode1(uint8_t a ) { allele1=a; } 

  int acode2() const { return allele2; } 
  void acode2(uint8_t a ) { allele2=a; } 


  // number of a particular allele, string encoding
  int allele_count( const std::string & , const Variant * ) const;  
  
  
  bool nonreference() const
  {
    return is_null ? false : allele1 != 0 || allele2 != 0 ;
  }
  
  bool heterozygote() const
  {
    if ( ploidy != 2 || is_null ) return false;
    return allele1 != allele2;
  }
  
  bool alternate_homozyote() const
  {
    if ( ploidy !=2 || is_null ) return false;
    return allele1 && allele2;
  }
  
  bool reference() const
  {
    if ( is_null || ploidy == 0 ) return false;
    if ( ploidy == 1 ) return allele1 == 0;
    return allele1 == 0 && allele2 == 0;
  }

  // does this genotype carry at least one copy of a minor allele-class?
  bool minor_allele( const bool reference_is_major ) const
  {
    if ( is_null || ploidy == 0 ) return false;
    if ( ploidy == 1 ) return reference_is_major ? allele1 : allele1 == 0;
    return reference_is_major ? allele1 || allele2 : allele1 == 0 || allele2 == 0;
  }
  
  // number of minor alleles
  int minor_allele_count( const bool reference_is_major ) const 
  {
    if ( is_null || ploidy == 0 ) return 0;
    if ( ploidy ==1 ) return reference_is_major ? allele1 : allele1 == 0;
    return reference_is_major ? allele1 + allele2 : ( allele1 == 0 ) + ( allele2 == 0 ); 
  }
  
  // given we have na total alleles in the population, return a vector of counts for each allele
  std::vector<int> allele_list( const int na ) const;
  
  
  //
  // Conversion function
  //
  
  
  // BCF encoding 
  
  uint8_t bcf() const;  
  void bcf( uint8_t );
  

  // VARDB encoding

  uint32_t pack() const; 
  bool unpack( uint32_t c);
  

   
  // Comparison between Genotypes
  
  static bool equivalent( const Genotype &, const Genotype & );
  
  
  bool operator==(const Genotype & rhs) const 
  {
    if ( is_null     !=  rhs.is_null     ) return false;
    if ( ploidy      !=  rhs.ploidy      ) return false;
    if ( allele1     !=  rhs.allele1     ) return false;
    if ( allele2     !=  rhs.allele2     ) return false;
    if ( known_phase !=  rhs.known_phase ) return false;
    return true;
  }
  

  bool operator!=(const Genotype & rhs ) const 
  {
    return !(*this == rhs);
  }

  
  friend std::ostream & operator<<( std::ostream & out, const Genotype & g)  
  {
    if ( g.is_null )
      {
	if ( g.ploidy == 2 ) out << ( g.known_phase ? ".|." : "./." ); 
	else out << ".";
      }
    else if ( g.ploidy == 1 ) out << (int)g.allele1; 
    else out << (int)g.allele1 << ( g.known_phase ? "|" : "/" ) << (int)g.allele2 ;
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

  friend class SampleVariant;
  
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
    
    void size(int n) { calls.resize(n); }
    
    //    void set_consensus_slotmap( SampleVariant * ps , std::vector<int> * pm );
    
    /// Add a Genotype to the GenotypeSet
    void add( const Genotype & g) { calls.push_back(g); }
    
    /// Add a Genotype at a particular individual index i
    void add( const Genotype & g, int i) { calls[i] = g; }
    
    /// Return Genotype for individual i
    Genotype & genotype(int i);
    
    /// Return const Genotype for individual i
    const Genotype & genotype(int i) const;
    
    /// Wipe all genotype calls (only for actual calls)
    void clear() { if ( ! svar ) calls.clear(); }

    void summarise_meta( std::map<meta_typed_key_t,std::pair<int,int> > & , 
		         std::map<meta_typed_key_t,std::string>  & ,
		         std::set<meta_typed_key_t> & ) const ;

    

  
};

#endif
