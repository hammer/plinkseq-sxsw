

#ifndef __VARIANT_H__
#define __VARIANT_H__

#include <vector>
#include <string>
#include <sstream>
#include <iostream>
#include <set>
#include <iterator>


#include "plinkseq/genotype.h"
#include "allele.h"
#include "plinkseq/helper.h"
#include "plinkseq/defs.h"
#include "plinkseq/meta.h"
#include "filemap.h"
#include "em.h"
#include "indmap.h"
#include "svar.h"

typedef unsigned int uint;

class Individual;
class Genotype;
class SampleVariant;


class Variant {
  
  
 public:

  //
  // Major embedded structures that are publicly-visible
  //

  /// A consensus/merged sample-variant
  SampleVariant consensus;
  
  /// Population-level meta-information
  MetaInformation<VarMeta>  meta;
    

  
  //
  // Primary interface for a Variant
  //

  
  /// Optionally specify whether Variant is 'valid' or not
  Variant( bool b = true );
  
  /// Construct variant given ID and chromosomal position
  Variant( const std::string & n, int c, int b );

  /// Attach an alignment (maps individuals to slots/IDs/phenotypes, etc) 
  void attach( IndividualMap * a) { align = a; } 
    
  /// Set name
  void name(std::string n) { vname = n; }

  /// Get name (if missing, return ".")
  const std::string name() const { return vname; }

  /// Set chromosome (integer)
  void chromosome(int c) { chr = c; }

  /// Set chromosome code (string)
  void chromosome(std::string c) { chr = Helper::chrCode(c); }

  /// Get chromosome code (integer)
  int chromosome() const { return chr; }
  
  /// Set position (leftmost position, bp1) 
  void position(int p) { bp = p; }

  /// Get position (leftmost position, bp1)
  int position() const { return bp; }

  /// Set end position (rightmost, bp2)
  void stop(int l) { bp2 = l; }

  /// Get end position (rightmost, bp2)
  int stop() const { return bp2 == 0 ? bp : bp2; }

  /// Set length span in bases (assumes position() specified)
  void length(int l) { bp2 = bp + l - 1; }

  /// Get length span in bases (bp2-bp1+1) 
  int length() const { return bp2-bp + 1 ; }

  /// Return a 'chr1:12345' type string
  std::string coordinate() const { return Helper::coordinate( chr, bp, bp2 ); } 

  
  //
  // Primary workhorse functions
  //
  
  /// Create a consensus variant given 1 or more SampleVariants
  bool  make_consensus( IndividualMap * );


  //
  // Genotype access functions (that might re-direct to consensus)
  //
  
  Genotype * genotype( const int i ) 
    {
      return & consensus.calls.genotype(i);
    }
  
  const Genotype * genotype( const int i ) const
    {
      return & consensus.calls.genotype(i);
    }

  Genotype * genotype( const int svar_id , const int i ) 
    {
      
      // We are being asked for the Genotype from a specific SampleVariant. 
      
      // The 'sv' index is 0-based and counts the total number of SVs for the Variant
      // (which may be larger than the number of samples, if there is infile-overlap)

      // The index 'i' is assumed to represent the position of the individual in the consensus
      // i.e. the single 0..N-1 mapping of all individuals
      
      // A number of scenarios:
      
      //  1) Reading directly from VCF:                      :  look straight to consensus
      //  2) Not multiple sample and no overlapping variants :  look straight to consensus
      //  3) Not multiple sample but w/ overlapping variants :  ??check -- individual SampleVariant
      //  4) Multiple samples but flat alignment             :  look at consensus
      //  5) Multiple samples but non-flat alignment         :  look at individual SampleVariant

      // i.e. if we have a non-flat alignment, or infile_overlap, then we need to consider the original SVs. 
      // otherwise, always go to the consensus
      

      // 1) Under simple scenarios, we can work directly from the consensus
      
      if (  svar_id == -1 ||  ( flat()  && ! infile_overlap() )  )
	{
	  return & consensus.calls.genotype( i );
	}
      
      // 2) Get the file this SV belong to 
      
      const int file_id = svtof[ svar_id ];
 
      // 3) Get the position in the file/SV given the individual's consensus slot 
      
      const int j = align->sample_slot( i , file_id );
      
      if ( j == -1 ) return NULL;

      // 4) Return request genotype
      
      return &svar[ svar_id ].calls.genotype( j );
      
    }



  //
  // Creating and populating new Variants
  //

  /// When directly reading a VCF, store the line for (potential) later parsing. The VCFReader
  /// will have created a single SVAR for this, so use that

  void set_vcf_buffer( const Helper::char_tok & tokens ) { svar[0].vcf_direct_buffer = tokens; }

  void set_vcf_buffer( const int i , std::vector<meta_index_t*> * f )
  { 
    consensus.vcf_direct = true; // need to set for both source and target 
    svar[0].vcf_direct = true;   // so that decode_BLOB functions do not try to use PB
    svar[0].vcf_gt_field = i;
    svar[0].vcf_formats = f;
  }

  
  /// Add a new Genotype to the consensus SampleVariant (when reading from VCF, and so straight to consensus)
  void add( const Genotype & g ) { consensus.calls.add(g); } 
  void add( const Genotype & g , const int i ) { consensus.calls.add(g,i); } 
  

  
  


  //
  // Some basic access functions for Sample Variants. 
  // By default these point to the consensus variant  
  //
  
  /// Get reference allele (from consensus SampleVariant)
  std::string reference() const { return consensus.ref; }
  
  /// Get string of alternate allele(s), as comma-delimited string list
  std::string alternate() const { return consensus.alt; }
  
  
  /// Pretty-print versions of the above (for large events)
  /// Get reference allele (from consensus SampleVariant)
  std::string pp_reference() const;
  std::string pp_alternate() const;
  
  /// Misc.

  void valid( bool b ) { is_valid = b; } 
 
  bool valid() { return is_valid; } 

  bool invalid() { return ! is_valid; } 

  
  //
  // Access to genotype data and helper functions
  //

  /// Test for concordance between two genotypes
  bool concordant( int s1, const Genotype * g1, int s2, const Genotype * g2 ) const;
  
  /// Test for concordance between two genotypes
  bool concordant( const SampleVariant * s1, const Genotype * g1, 
		   const SampleVariant * s2, const Genotype * g2 ) const;


  /// Printing access to meta-information 

  std::string print_samples( const std::string & delim = " " ) const;
  std::string print_meta( const std::string & key , const std::string & delim = " ") const;
  std::string print_meta_filter( const std::string & delim = " ") const;
  std::set<std::string> meta_filter( ) const;
  
  //
  // Access to individuals (assumes an attached IndidivualMap)
  //

  /// Number of individuals with data for this variant  
  int size() const;
  
  /// Numner of individuals in svar 'si'
  int size(const int si) const;

  /// Set number of individuals with data for this variant  
  void resize(const int n);

  /// Reserve space for calls
  void reserve(const int n);

  /// Pointer to an individual (0..N-1) in the consensus/indmap, given N
  Individual * ind(const int) const;
  
  /// Pointer to an individual (0..N-1) in the consensus/indmap, given ID  
  Individual * ind(const std::string &) const;
  
  /// All IDs from the consensus indmap
  std::vector< std::string> ind_id() const { return align->ind_id(); }

  /// return N for individual given their ID
  int ind_n( const std::string & ) const;
  
  /// does the indmap contain more than one sample?
  bool multi_sample() const { return align->multi_sample(); }
  
  /// is the indmap alignment flat (i.e. no individual seen >1 time)
  bool flat() const { return align->flat() ; }

  /// show which fileset (1..NS) an individual belongs to, if they belong to a single fileset (else returns 0) 
  int ind_sample( const int i ) const { return align->sample(i); }
  

  //
  // Primary functions to access a (consensus) genotype, via a returned pointer
  //
  
  /// Overload () to provide references to consensus genotypes  
  Genotype & operator()(const int i) { return consensus.calls.genotype(i); }

  /// Overload () to provide const references to consensus genotypes  
  const Genotype & operator()(const int i) const { return consensus.calls.genotype(i); }

  /// All genotypes for this Variant for a given individual (SV->genotype map)
  std::map<int, Genotype *> all_genotype( const int );
  
  /// All genotypes for a given individual (const version)
  std::map<int,const Genotype *> all_genotype(const int) const;

  /// From a specific sample (but which will be consensus under a flat alignment)  
  const Genotype * genotype( const SampleVariant * svar , const int i) const;
  

  //
  // Sample-level access functions
  //


  /// Is any one sample repeated for this Variant (1:many mapping of sample:sample-variant)
  bool infile_overlap() const { return svtof.size() > ftosv.size(); }

  /// Return list of all SV IDs
  std::vector<int> samples() const;
  
  /// Return set of only unique file-IDs
  std::set<int> unique_files() const ;
  
  /// Return svar-slot given a file ID (or -1 if no uniq mapping)
  int unique_svar_slot( int f ) const;

  /// Is this file-ID present in the Variant?
  bool file_present( const int f ) const;

  
  /// Return a pointer to a given sample
  SampleVariant * psample(const int s) 
    {
      if ( s == -1 ) return &consensus;
      return s < 0 || s >= (int)svar.size() ? NULL : &(svar[s]);
    }
  
  const SampleVariant * psample(const int s) const
  {
    if ( s == -1 ) return &consensus;
    return s < 0 || s >= (int)svar.size() ? NULL : &(svar[s]);
  }
  
  const SampleVariant & sample(const int s) const 
  {
    return s == -1 ? consensus : svar[ s ] ;
  }
  
  SampleVariant & sample(const int s) 
    {
      return s == -1 ? consensus : svar[ s ] ;
    }

  /// Point to svar that holds genotypes 
  SampleVariant & sample_genotypes( const int s ) const
    {
      return flat() ? (SampleVariant&)consensus : (SampleVariant&)svar[s] ; 
    }
  
  SampleVariant & sample_genotypes( const SampleVariant & sv ) const
    {
      return flat() ? (SampleVariant&)consensus : (SampleVariant&)sv ; 
    }

  SampleVariant * sample_genotypes( const SampleVariant * sv ) const
    {
      return flat() ?  (SampleVariant*)&consensus : (SampleVariant*)sv ;
    }


/*   SampleVariant & sample_metainformation( const int s ) const */
/*     { */
/*       return multi_sample() ? (SampleVariant&)svar[s] : (SampleVariant&)consensus; */
/*     } */

  
  SampleVariant & sample_metainformation( const SampleVariant & sv ) const 
    { 
      return multi_sample() ? (SampleVariant&)sv : (SampleVariant&)consensus; 
    } 
  
  SampleVariant * sample_metainformation( const SampleVariant * sv ) const 
  { 
    return multi_sample() ? (SampleVariant*)sv : (SampleVariant*)&consensus; 
  } 


  /// As above, but via external file codes
  /// Return the number of SVs in a file

  int fsample_svar_counts(const int f) const 
  {
    std::map<int,std::vector<int> >::const_iterator i = ftosv.find(f);
    return  i == ftosv.end() ? 0 : i->second.size(); 
  }
  
  /// Return the i'th SV from the f'th file
  std::vector<SampleVariant *> fsample( const int file_id ) 
    {
      std::map<int,std::vector<int> >::iterator i = ftosv.find( file_id );
      std::vector<SampleVariant*> s;
      if ( i == ftosv.end() ) return s;
      std::vector<int>::iterator j = i->second.begin();
      while ( j != i->second.end() )
	{
	  s.push_back( psample( *j ) );
	  ++j;
	}
      return s;
    }
  
  std::vector<const SampleVariant *> fsample( const int file_id ) const 
    {
      std::map<int,std::vector<int> >::const_iterator i = ftosv.find( file_id );
      std::vector<const SampleVariant*> s;
      if ( i == ftosv.end() ) return s;
      std::vector<int>::const_iterator j = i->second.begin();
      while ( j != i->second.end() )
	{
	  s.push_back( psample( *j ) );
	  ++j;
	}
      return s;
    }



  //
  // Allele-level summaries/queries
  //

  /// The number of specified alleles, from the REF/ALT
  int n_alleles() const;

  /// Does allele-count equal 2?
  bool biallelic() const;

  /// Does allele-count equal 1?
  bool monomorphic() const;

  /// Is allele-count greater than 2?
  bool multiallelic() const;

  /// Is this a bi-allelic single-base substitution? (i.e. basic SNP)
  bool simple_snp() const;
  
  /// Is this a simple SNP and a transition (A/G or C/T SNPs)
  bool transition() const;

  /// Is this a simple SNP and a transversion (not A/G or C/T SNPs)
  bool transversion() const;

  /// Is this a biallelic insertion (ALT length > REF )
  bool simple_ins() const;

  /// A biallelic deletion? (ALT length < REF )
  bool simple_del() const;
  
  /// Is this an indel?
  bool indel() const;

  /// Is this an MNP?
  bool mnp() const;

  /// Is this a SNP?  (could be multi-allelic)
  bool snp() const;

  /// Return the 'Allele' object for allele 'a'
  const Allele & allele(const int a) const;

  
  //
  // Genotype-level summaries
  //
  
  
  std::vector<int> indiv_mask( const int file_id ) const;

  bool has_nonreference_by_file( const int fide_id ) const;

  bool has_nonreference( const SampleVariant & svar ) const;

  /// Return a map of 
  std::map<std::string,int> allele_counts( const affType & aff = UNKNOWN_PHE ) const;
  
  /// Return the number of individuals with a non-reference (non-missing) genotype
  int n_nonreference() const;


  std::map<std::string,int> genotype_counts( const SampleVariant & , const affType & aff , bool unphased = true ) const;

  std::map<std::string,int> genotype_counts( const int si , const affType & aff , bool unphased = true ) const;

  std::map<std::string,int> genotype_counts( const affType & aff , bool unphased = true ) const;


  /*!
  Return number of
       @param m Reference: returns number of alternate alleles
       @param n Reference: returns total number of non-missing alleles
   */

  void n_alt_allele( int * m = NULL , int * n = NULL , double * aaf = NULL , const affType & aff = UNKNOWN_PHE ) const;

  /*!
    Return number of 
    @param m Reference: returns number of minor alleles 
    @param n Reference: returns total number of non-missing alleles
    @return True if alternate allele is minor allele
  */

  bool n_minor_allele( int * m = NULL , int * n = NULL , double * maf = NULL , const affType & aff = UNKNOWN_PHE ) const;
  
  /*! 
    Queries about genotype data
    @return Number of null genotypes
  */
  
  int n_null() const;

  /*! 
    Queries about genotype data
    @return Number of non-null genotypes
  */

  int n_notnull() const;


  /*!
    Given a mask, return T/F for the frequency filter  
  */
  
  bool frequency_filter(Mask *);
  

  /*!
    Given a mask, return T/F for the null-genotype filter
  */
  
  bool null_filter(Mask *);

  /*!
    Given a mask, return T/F for case/control count filter
  */
  
  bool case_control_filter(Mask *);



  /*!
    @brief Label for a genotype
    @param g Reference to Genotype object
    This function assumes that parse_alleles() has been called
    
    @return String value for name of genotype, e.g. A/C
  */
  

  std::string sample_label( const int, const std::string & delim = "," ) const;
  std::string geno_label( const Genotype & ) const;
  std::string geno_label( const int, const Genotype & ) const;
  std::string phased_geno_label( const Genotype & ) const;
  std::string phased_geno_label( const int, const Genotype & ) const;
  std::string label( const int, const std::string & delim = "," ) const;
  std::string gmeta_label( const int, const std::string & delim = "," ) const;
  std::string alternate_label( const int, const std::string & delim = "," ) const;
  std::string sample_variant_index( const int , const std::string & delim = "," ) const;
  std::string allele1_label( const Genotype & ) const;
  std::string allele2_label( const Genotype & ) const;

  
  /*!
    @brief PED-stype genotype string, for SNPs only
    @param g Reference to Genotype object
    @param Alternate delimiter
    This function assumes that parse_alleles() has been called
    
    @return String value for name of genotype, e.g. A C
  */
  
  std::string print_PED(const Genotype & g, const std::string & delim = " " ) const;
  

  // 
  // Encoding/decoding
  //
    
  /*!
    @brief Write variant and genotypes as a VCF string

    @return Single line in VCF representing variant and calls
  */
  
  std::string VCF() const;


  std::string displaycore() const 
    { 
      std::stringstream ss;
      ss << Helper::chrCode( chr ) << ":" << bp;
      if ( bp2 != 0 && bp2 != bp ) ss << ".." << bp2;
      ss << ":" << vname;			      
      return ss.str();
    }


  friend std::ostream & operator<<( std::ostream & out, const Variant & v);
  
  bool single_sample() const 
  {
    return svar.size() == 0 ;
  }
  
  int n_samples() const 
  {
    return svar.size();
  }
  
  int n_uniq_samples() const 
  {
    return ftosv.size();
  }



  /// Create a new SampleVariant, from a given fileset
  
  SampleVariant & add(const int f)  
    { 
      SampleVariant sample; 
      sample.fileset(f); 
      svar.push_back( sample );
      svtof.push_back(f);
      ftosv[f].push_back( svar.size()-1 );
      return svar.back(); 
    } 

  /// Second function to add a new SampleVariant to a Variant

  SampleVariant & add( SampleVariant & sample )
    {
      svar.push_back( sample );
      svtof.push_back( sample.fileset() );
      ftosv[ sample.fileset() ].push_back( svar.size() -1 );
      return svar.back();      
    }
  
  bool remove( int s ) ;

  
  //
  // Some secondary classes embedded here
  //


  
  
  /// Helper class for E-M estimate of allele frequency and posterior genotype probabilities
  EM em;


 private:

  //
  // Core, invariant information for a variant
  //

  /// Variant ID
  std::string      vname;

  /// Chromosome code: as integer, decoded by Helper::chrCode()
  int              chr;
  
  /// Reference base-pair position (base-1)
  int              bp;

  /// Optional second position (also base-1); if 0, implies bp2 == bp
  int              bp2;
  
  /// Is this a valid Variant (i.e. was it properly defined in the VCF)?
  bool             is_valid;


  /// list of sample variants that contain most of the data
  std::vector<SampleVariant> svar;
  

  //
  // Because we allow overlapping or repeated variants in the same file, 
  // there is a one-to-many mapping of file to sample-variant. These two
  // structures keep track of this for this Variant
  //
  
  /// map SV-IDs (0..NS-1) to file-IDs (any 1+ int)
  std::vector<int> svtof;
  
  /// file-ID --> list of SV-IDs (file-to-SampleVariant)
  std::map<int,std::vector<int> > ftosv; 

  /// Simple flag to indicate if the Sample Variants come from more than one File
  /// QUESTION: should this be SV?
  bool             is_multi_sample;

  /// 'iterator' for looping over sample-variants
  int si;

  
  /// Pointer to individual alignment class
  IndividualMap * align;
  

  //
  // Private member functions
  //


  /// New variants constructed as missing 
  void init();
  
    
 public:

  bool operator<( const Variant & rhs ) const
  {
    if ( chr == rhs.chr )
      {
	if ( bp == rhs.bp )
	  {
	    if ( bp2 == rhs.bp2 ) 
	      {
		return consensus.alt < rhs.consensus.alt;
	      }
	    return bp2 < rhs.bp2;
	  }
	return bp < rhs.bp;
      }
    else
      return chr < rhs.chr;
  }
  


};



#endif
  


