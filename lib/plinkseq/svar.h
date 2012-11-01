#ifndef __PSEQ_SVAR_H__
#define __PSEQ_SVAR_H__


#include "plinkseq/meta.h"
#include "plinkseq/variant.pb.h"
#include "allele.h"
#include "plinkseq/genotype.h"

#include <ostream>

typedef unsigned int uint;

class Genotype;
class Mask;
class BCF;
class VarDBase;
class IndividualMap;


/*!
  
  This is a key class, which represents a variant. This means basic
  information such as chromosomal location, name, alleles, etc, as
  well as extensible meta-information in the MetaInformation member
  class. A variant can also hold genotype calls on 1 or more individuals,
  in the calls members class.
  
  @class SampleVariant
  
  @brief Core class representing individual variants for one file

*/


class SampleVariant {

  friend class Variant;
  friend class VarDBase;
  friend class Eval;
  friend class BCF;

 public:
     

  SampleVariant() { init(); }

    
  //
  // Publicly-visible meta-information classes
  //

  /// SampleVariant meta-information class  
  MetaInformation<VarMeta>        meta;
  

  /// Filter information
  MetaInformation<VarFilterMeta>  meta_filter;


  
  //
  // Public access functions
  //

  /// File index
  
  void index( uint64_t v ) { vindex = v; }
  uint64_t index() const { return vindex; }

  
  /// Set fileset ID
  void fileset(const int f) { fset = f; }
  
  /// Get fileset ID
  const int fileset() const { return fset; }
  
  /// Set reference allele (string, 1+ characters)
  void reference(const std::string & s) { ref = s; }
  
  /// Get reference allele (string, 1+ characters)
  std::string reference() const { return ref; }
  
  /// Set alternate allele(s), comma-delimited list, e.g. A,T
  void alternate(const std::string & s) { alt = s; }
  
  /// Get string of alternate allele(s), as comma-delimited string list
  std::string alternate() const { return alt; }
  
  /// Return number of alleles (one reference + alternates):
  int nalleles() const {
	  return alleles.size();
  }

  /// Return number of alternate alleles:
  int nalt() const {
	  return nalleles() - 1;
  }

  /// Return allele label of k-th allele (0 is reference)
  std::string alternate( const uint k ) const
    {
      if ( k >= alleles.size() ) return ".";
      return alleles[k].name();
    }
  
  /// pretty print versions of the above
  std::string pp_reference() const;
  std::string pp_alternate() const;


  /// If a VARDB attached, get file-name given fileset() value
  std::string file_name() const;


  
  
  //
  // Basic get/set functions
  //
  
  /// Set quality score, expected positive floating-point value
  void quality(double d) { qual = d; }

  /// Get quality score, 
  double quality() const { return qual; }
  
  /// Set meta-information (info and meta class; assume semi-colon delimited)
  void info( const std::string & s , VarDBase * vardb = NULL , int file_id = 0 , Variant * parent = NULL ) ;

  /// Get meta-information string as specified by void info(string)
  std::string info() const { return other_info; }
  
  /// Set filter string; currently no further processing 
  void filter( const std::string & s , VarDBase * vardb = NULL , int file_id = 0);

  /// Get filter string
  std::string filter() const { return filter_info; }

  /// Get set of filter strings
  std::vector< std::string > filters() const { return meta_filter.get_flags(); }

  /// Has as specific filter value?
  bool has_filter( const std::string & s ) const;
  
  /// Has any filter value (i.e. non-PASS in FILTER)?
  bool any_filter() const { return filter_info != PLINKSeq::PASS_FILTER(); }
		   

  //
  // Allow adding genotypes directly, but all reading must be via the Variant interface
  // (to simplify the possible redirection to consensus). In the case of writing, the 
  // calling function will have taken care of what is target for meta-information vs. 
  // genotypes, etc
  //
  
  //
  // If reading from BCF, here is pointer to that class. The BCF class does all the work of 
  // parsing, and is also given then IndividualMap and Mask, so it will only extract what 
  // is needed.  (Note: this is in contrast to BGZF-compressed VCF, where the fine-grained
  // extraction logic is done here in SamplVariant, because the whole text line needs to be 
  // read from disk first).
  //

  void set_pointer_to_bcf( BCF * p , int64_t offset ) { bcf = p; bcf_offset = offset; } 
  
  

  //
  // As above, for BGZF-compressed VCFs (note; when procesing a single VCF from the
  // command line (i.e. that isn't indexed in the VARDB, this is handled differently, from
  // vcfiterate.cpp and the setting of the buffer is done via the parent Variant 
  //

  void set_vcfz_buffer( const Helper::char_tok & buffer ,      // genotypes and meta-information
			int gt_field ,                         // from FORMAT, slot containing GT
			std::vector<meta_index_t*> * formats ) // parsed meta-key vector for FORMAT
  {
    vcf_direct = true;
    vcf_direct_buffer = buffer;
    vcf_gt_field = gt_field;
    vcf_formats = formats;
  }
  


  //
  // Convert BLOB to PB format to SampleVariant structure
  //

  // Variant --> PB --> BLOB
  
  blob encode_var_BLOB() const;
  blob encode_vmeta_BLOB() const;
  blob encode_geno_BLOB() const;
  blob encode_gmeta_BLOB() const;


  // BLOB --> unparsed PB

  void store_BLOBs( blob * , blob * , blob * , blob * );
  
 private:


  /// Unique index from VARDB, used in construction

  uint64_t              vindex;

  
  // Allelic REF/ALT encoding
  
  std::string           ref;
  
  std::string           alt;
  
    
  
  // 
  // Sample-specific meta-information (as specified in VCF)
  // 
  
  /// QUAL field from VCF
  double           qual;
  

  /// FILTER field from VCF --> meta_filter
  std::string      filter_info;


  /// INFO field from VCF --> meta
  std::string      other_info; 

    
  //
  // Other information generated upon loading
  //
  
  
  void              set_allelic_encoding();
 
  
  /// REF base-offset: i.e. REF starts 'offset' bp after the Variant BP1
  int               offset;


  /*! 
    Construct internal alleles[] vector from reference() allele value(s)
    @return Number of alleles 
  */
  
  int parse_alleles();

  //
  // Representation of all alleles observed 
  //


  std::vector<Allele>   alleles;
  
  //
  // Collapse all alternate alleles
  //
  
  void collapse_alternates( const Variant * , int altcode = 0 );
  

  //
  // Helper function to recode overlapping alleles
  //
  
  static bool align_reference_alleles( SampleVariant & s1 , 
				       SampleVariant & s2 , 
				       bool add_alt = false );


  int              fset;  

    
  // Buffers (BLOB/PB -> variant)
  
  VariantBuffer      var_buf;
  VariantMetaBuffer  vmeta_buf;
  GenotypeBuffer     geno_buf;
  GenotypeMetaBuffer gmeta_buf;
  
  
  // Or, direct from a single VCF (and so no BLOB)
  
  bool                         vcf_direct;
  Helper::char_tok             vcf_direct_buffer;
  int                          vcf_gt_field;
  std::vector<meta_index_t*> * vcf_formats;
  
  // Or, from a BCF (and so no BLOB to decode)

  BCF *                 bcf;
  int64_t               bcf_offset;

  /// Genotypes on 1 or more individuals for this variant

  GenotypeSet      calls;

  /// Overload () to access genotype calls for individual i
  
  Genotype & operator()(const int i) { return calls.genotype(i); }

  const Genotype & operator()(const int i) const { return calls.genotype(i); }
     
  /// Initialise a SV as null

  void init()
  {
    vindex = 0 ;
    fset = 0;
    alt = ref = ".";
    filter_info = ".";
    other_info = ".";
    qual = -1;    
    meta.clear();
    calls.clear();     
    alleles.clear();
    offset = 0;
      
    // for VCF --> SV
    vcf_direct = false;
    vcf_direct_buffer.clear();
    vcf_gt_field = 0;
    vcf_formats = NULL;
    
    // for BCF --> SV
    bcf = NULL;
    
  }
  
  
  /// Has at least one non-reference call
  bool has_nonreference( const bool also_poly = false , const std::vector<int> * remap = NULL ) const;
  

  //
  // Input/output and recoding of genotypes
  //


  /// Recall a genotype: i.e. new numeric coding, given a new SV REF/ALT codes  
  void recall( Genotype & g , SampleVariant * );
  
  
  /// Write core variant information to stream (no meta or genotype information)  
  friend std::ostream & operator<<( std::ostream & out, const SampleVariant & v);

  /// ACGT encoding of a genotype (by default, collapse phased/unphased)
  std::string label( const Genotype & g , bool phased = false ) const;
  
  /// Numeric (0/0, etc) encoding of a genotype
  std::string num_label( const Genotype & g ) const;

  /// Similar to the above, but for individual alleles
  std::string allele1_label( const Genotype & g ) const;
  std::string allele2_label( const Genotype & g ) const;


  //
  // Allele/genotype count functions
  //

  std::map<std::string,int> allele_counts( const affType & aff , const Variant * parent) const;
  
  std::map<std::string,int> allele_count(const int ) const;
  
  std::map<std::string,int> genotype_counts( const affType & aff , const Variant * parent , bool phased = false ) const;
  

  

  // PB --> Variant
  
  bool decode_BLOB( Variant * , IndividualMap * , Mask * );
  
  bool decode_BLOB_basic( SampleVariant * target );

  void decode_BLOB_alleles();
 
  bool decode_BLOB_vmeta( Mask * mask ,             // for filters
			  Variant * parent        , // for pop/static meta-info
			  SampleVariant * target ); // for straight-to-consensus
  
  bool decode_BLOB_genotype( IndividualMap * align , 
			     Mask * mask ,
			     Variant * parent , 
			     SampleVariant * source ,
			     SampleVariant * vtarget ,
                             SampleVariant * target );
  


  // ---------------------------------------------------------------------------
  // private SampleVariant members
  
 private:
  
  /// Helper function for reading in int genotype meta-information
  inline int addIntGenMeta( int j ,                        // ?
			    int f ,                        // ?
			    const GenotypeMetaBuffer & v,  // data source
			    IndividualMap * align,         // indiv-alignment  
			    int k,                         // meta-info slot
			    int idx,                       // current counter 		     
			    int l );                       // length arg 
  
  /// Helper function for reading in float genotype meta-information  
  inline int addFloatGenMeta( int j , int f , const GenotypeMetaBuffer & v, 
			      IndividualMap * align, int k, int idx, int l );    
  
  /// Helper function for reading in string genotype meta-information
  inline int addStringGenMeta( int j , int f , const GenotypeMetaBuffer & v, 
			       IndividualMap * align, int k, int idx, int l );  
    
  /// Helper function for reading in bool genotype meta-information
  inline int addBoolGenMeta( int j , int f , const GenotypeMetaBuffer & v, 
			     IndividualMap * align, int k, int idx, int l ); 
  

};




#endif
