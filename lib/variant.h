
#ifndef __VARIANT_H__
#define __VARIANT_H__

#include <vector>
#include <string>
#include <sstream>
#include <iostream>
#include <set>
#include <iterator>

#include "variant.pb.h"

#include "genotype.h"
#include "allele.h"
#include "helper.h"
#include "defs.h"
#include "meta.h"
#include "filemap.h"
#include "em.h"
#include "indmap.h"
#include "spec.h"

class Individual;
class Mask;
class Genotype;
class BCF;

/*!
  
  This is a key class, which represents a variant. This means basic
  information such as chromosomal location, name, alleles, etc, as
  well as extensible meta-information in the MetaInformation member
  class. A variant can also hold genotype calls on 1 or more individuals,
  in the calls members class.
  
  @class SampleVariant
  
  @brief Core class representing individual variants for one file

*/

struct SampleVariant {
  

  //
  // A particular sample/variant combination
  //
  
  uint64_t         vindex;
  
  
  //
  // Allelic, genotypic encoding
  //
  
  std::string           vencoding;
  
  bool                  simple;

  std::string           ref;
  
  std::string           alt;
  
  std::string           vstrand;  
  
  
  // 
  // Sample-specific meta-information (as specified in VCF)
  // 
  
  double           qual;
  
  std::string      filter_info;
   
  bool             is_filtered;
  
  vType            vtype;
  
  std::string      other_info;  // converted -> meta()

    
  //
  // Other information generated upon loading
  //
  

  VariantSpec*      spec;
  
  void              set_allelic_encoding();
 

  // If this is >0, means this REF starts 'offset'
  // bp after the Variant bp

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
  
  void collapse_alternates( int altcode = 0 );
  

  //
  // Helper function to recode overlapping alleles
  //

  static bool align_reference_alleles( SampleVariant & s1 , SampleVariant & s2 , bool add_alt = false );


  int              fset;
  
    
  // Buffer (BLOB/PB -> variant)
  
  VariantBuffer buf;
  
  // Or, direct from a single VCF (and so no BLOB)
  
  bool                       vcf_direct;
  std::vector<std::string>   vcf_direct_buffer;
  bool                       vcf_expand_buffer( Variant * );
  
  // Or, from a BCF (and so no BLOB to decode)

  BCF *                 bcf;
  std::string           bcf_format;
  std::vector<uint8_t>  bcf_genotype_buf;
  


 public:  
     
  SampleVariant() 
    { 
      init(); 
    }

  
  SampleVariant(bool b) 
    { 
      init(); 
    }
  
  SampleVariant(uint64_t idx)
    {		
      init();      
      vindex = idx;
    }
  
  /*!
    Unless otherwise stated, variants are initialised as completely missing.
  */

  void init()
    {
      vindex = 0;
      fset = 0;
      alt = ref = ".";
      filter_info = ".";
      other_info = ".";
      qual = -1;
      is_filtered = false;
      vstrand = ".";
      meta.clear();
      calls.clear();       
      spec = NULL;
      simple = true;
      alleles.clear();

      vcf_direct = false;
      vcf_direct_buffer.clear();

      // for BCF
      offset = 0;
      bcf = NULL;
      bcf_genotype_buf.clear();
    }
 
  
  /// Variant specification 
  
  static specDecoder      decoder;
  
  
  /// SampleVariant meta-information class
  
  MetaInformation<VarMeta>  meta;


  /// Filter information

  MetaInformation<VarFilterMeta> meta_filter;


  /// Genotypes on 1 or more individuals for this variant

  GenotypeSet      calls;

  /// Overload () to access genotype calls for individual i
  
  Genotype & operator()(const int i) { return calls.genotype(i); }
  
  const Genotype & operator()(const int i) const { return calls.genotype(i); }
  
  /// Has at least one non-reference call
  bool has_nonreference( const bool also_poly = false ) const;
  
  /// Set variant index
  void index(uint64_t i) { vindex = i; }  

  /// Get variant index
  const uint64_t index() const { return vindex; }  


  /// Set fileset ID (probably redundant now)
  void fileset(int f) { fset = f; }

  /// Get fileset ID
  const int fileset() const { return fset; }

  /// Set reference allele (string, 1+ characters)
  void reference(std::string s) { ref = s; }

  /// Get reference allele (string, 1+ characters)
  std::string reference() const { return ref; }

  /// Set alternate allele(s), comma-delimited list, e.g. A,T
  void alternate(std::string s) { alt = s; }
  
  /// Get string of alternate allele(s), as comma-delimited string list
  std::string alternate() const { return alt; }
  
  /// Return allele label of k-th allele (0 is reference)
  std::string alternate( const int k ) const
  {
    if ( k < 0 || k >= alleles.size() ) return ".";
    return alleles[k].name();
  }
  

  std::string file_name() const;

  std::string label( const Genotype & g , bool unphased = false ) const;

  std::string num_label( const Genotype & g ) const;

  std::map<std::string,int> allele_counts( const affType & aff , const Variant * parent) const;

  std::map<std::string,int> allele_count(const int ) const;

  std::map<std::string,int> genotype_counts( const affType & aff , const Variant * parent , bool unphased = true ) const;
  
  /// Set strand ("+", "-" or "?")
  void strand(std::string s) { vstrand = s; }

  /// Get strand
  std::string strand() const { return vstrand; }

  /// Set quality score, expected positive floating-point value
  void quality(double d) { qual = d; }

  /// Get quality score, 
  double quality() const { return qual; }
  
  /// Set meta-information: keep as string, but also set in MetaInformation class meta
  /// Uses semi-colon only as field separator

  void info( const std::string & s , VarDBase * vardb = NULL , int file_id = 0 ) ;

  /// Get meta-information string as specified by void info(string)
  std::string info() const { return other_info; }

  
  /// Set filter string; currently no further processing 
  void filter( const std::string & s , VarDBase * vardb = NULL , int file_id = 0);

  /// Get filter string
  std::string filter() const { return filter_info; }

  /// Get set of filter strings
  std::vector< std::string > filters() const 
  {   
    return meta_filter.get_flags();
  }

  /// Has a filter value?
  bool has_filter( const std::string & s ) const;
  
  bool any_filter() const { return filter_info != PLINKSeq::PASS_FILTER(); }

		   

  /// Set boolean values on whether variant is filtered out (if T)
  void filtered(bool b) { is_filtered = b; }

  /// Get boolean value, T is variant is to be filtered out
  bool filtered() const { return is_filtered; }



  // Check which of the following are redundant...

  void type(vType t) { vtype = t; }

  vType type() { return vtype; }


  // Get/set variant specification 

  void specification(VariantSpec *s) { spec = s; }

  VariantSpec * specification() const { return spec; }  
  
  std::string encoding() { return vencoding; }

  void encoding(std::string s) { vencoding = s; }
  

  
  /// Write core variant information to stream (no meta or genotype information)
  
  friend std::ostream & operator<<( std::ostream & out, const SampleVariant & v);
/*     {  */
/*       out << GP->vardb.file_tag( v.fset ) */
/* 	  << ":" << v.ref << "/" << v.alt */
/* 	  << ":" << v.filter_info;       */
/*       return out; */
/*     } */
  
      
  blob encode_BLOB() const;
  void store_BLOB(blob&);

/*   bool decode_BLOB( blob & ,  */
/* 		    Variant * parent = NULL); */

  bool decode_BLOB( Variant * ,
		    IndividualMap * , 
		    Mask * );

  bool decode_BLOB_basic( SampleVariant * target );

  bool decode_BLOB_vmeta( Mask * mask ,             // for filters
			  Variant * parent        , // for pop/static meta-info
			  SampleVariant * target ); // for straight-to-consensus
  
  bool decode_BLOB_genotype( IndividualMap * align , 
			     Mask * mask ,
			     Variant * parent , 
			     SampleVariant * source ,
			     SampleVariant * vtarget ,
                             SampleVariant * target );
  

 private:


  //
  // Helper functions for reading in genotype meta-information
  //
  

 inline int addIntGenMeta( int j , int f , 
			   const VariantBuffer & v, 
			   IndividualMap * align, 
			   int k,   // meta-info slot
			   int idx, // current counter 		     
			   int l )  // length arg 
  {
    
    if ( align )
      {
	j = align->sample_remapping( f , j);
	
	if ( align->flat() )
	  j = align->get_slot( f , j );
      }
    
    if( j == -1 ) 
      return idx + l; 
          
    MetaInformation<GenMeta> & gmeta = calls.genotype( j ).meta ;
    
    if ( l == 1 ) 
      {
	gmeta.set( v.gmeta(k).name() , 
		   v.gmeta(k).int_value( idx++ ) );
      }
    else
      {
	std::vector<int> t(l);
	for ( int i = 0 ; i < l; i++)
	  t[i] = v.gmeta(k).int_value( idx++ );
	gmeta.set( v.gmeta(k).name() , t );
      }
    
    return idx;      
  }
  



  
  inline int addFloatGenMeta( int j , int f , 
			      const VariantBuffer & v, 
			      IndividualMap * align, 
			      int k,   // meta-info slot
			      int idx, // current counter 		     
			      int l )  // length arg
    {
      
      if ( align )
	{
	  j = align->sample_remapping( f , j);
	  
	  if ( align->flat() )
	    j = align->get_slot( f , j );
	}

      if( j == -1 ) return idx + l;
      MetaInformation<GenMeta> & gmeta = calls.genotype( j ).meta ;

      if ( l == 1 ) 
	gmeta.set( v.gmeta(k).name() , 
		  v.gmeta(k).double_value( idx++ ) );
      else
	{
	  std::vector<double> t(l);
	  for ( int i = 0 ; i < l; i++)
	    t[i] = v.gmeta(k).double_value( idx++ );
	  gmeta.set( v.gmeta(k).name() , t );
	}
      
      return idx;
      
    }
  
  inline int addStringGenMeta( int j , int f , 
			       const VariantBuffer & v, 
			       IndividualMap * align, 
			       int k,   // meta-info slot
			       int idx, // current counter 		     
			       int l )  // length arg
    {
      
      if ( align )
	{
	  j = align->sample_remapping( f , j);
	  
	  if ( align->flat() )
	    j = align->get_slot( f , j );
	}

      if( j == -1 ) return idx + l;
      MetaInformation<GenMeta> & gmeta = calls.genotype( j ).meta ;
      
      if ( l == 1 ) 
	gmeta.set( v.gmeta(k).name() , 
		   v.gmeta(k).string_value( idx++ ) );
      else
	{
	  std::vector<std::string> t(l);
	  for ( int i = 0 ; i < l; i++)
	    t[i] = v.gmeta(k).string_value( idx++ );
	  gmeta.set( v.gmeta(k).name() , t );
	}  
      return idx;  
    }
  
  inline int addBoolGenMeta( int j , int f , 
			    const VariantBuffer & v, 
			    IndividualMap * align, 
			    int k,   // meta-info slot
			    int idx, // current counter 		     
			    int l )  // length arg
   {

      if ( align )
	{
	  j = align->sample_remapping( f , j);
	  
	  if ( align->flat() )
	    j = align->get_slot( f , j );
	}

      if( j == -1 ) return idx + l;
      MetaInformation<GenMeta> & gmeta = calls.genotype( j ).meta ;

     if ( l == 1 ) 
       gmeta.set( v.gmeta(k).name() , 
		  v.gmeta(k).bool_value( idx++ ) );
     else
       {
	 std::vector<bool> t(l);
	 for ( int i = 0 ; i < l; i++)
	   t[i] = v.gmeta(k).bool_value( idx++ );
	 gmeta.set( v.gmeta(k).name() , t );
       }
     return idx;
   }
  
 
};




class Variant {
  
  

  // 
  // Core, invariant information for a variant
  //

  std::string      vname;
  
  int              chr;
  
  int              bp;
  
  int              bp2;
  


  //
  // Storage of sample-specific variants
  //

  /// vector of sample variants
  std::vector<SampleVariant> svar;
  
  /// map SV-IDs (0..NS-1) to file-IDs (any 1+ int)
  std::vector<int> svtof;

  /// file-ID --> list of SV-IDs (file-to-SampleVariant)
  std::map<int,std::vector<int> > ftosv; 
    
  /// 'iterator' for looping over samples
  int si;
  
  

  //
  // Pointer to individual alignment class
  //
  
  IndividualMap * align;
  

  //
  // Misc.
  //

  bool             is_valid;

  bool             is_multi_sample;
  
 public:


  Variant() 
    { 
      init(); 
    }


  Variant(bool b) 
    { 
      init(); 
      is_valid = b;  
    }
  
  Variant(std::string n, int c, int b)
    {		
      init();      
      vname = n;
      chr = c;
      bp = bp2 = b;      
    }
  

  /*!
    Unless otherwise stated, variants are initialised as completely missing.
  */

  void init()
    {
      chr = bp = bp2 = 0;
      vname = ".";
      meta.clear();
      align = NULL;
      is_multi_sample = false;
     }
  

  void attach( IndividualMap * a)
  {
    align = a;
  }
  
  EM em;
  

  //
  // A consensus/merged variant
  //
  
  SampleVariant consensus;
  
  bool  make_consensus(IndividualMap*);
  
  
  //
  // Population-level meta-information
  //

  MetaInformation<VarMeta>  meta;
  
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

  /// 
  std::string coordinate() const { return Helper::coordinate( chr, bp, bp2 ); } 
  

  //
  // Some basic access functions for the consensus variant
  //

  /// Get reference allele (string, 1+ characters)
  
  std::string reference() const 
    { 
      return consensus.ref; 
    }
  
  /// Get string of alternate allele(s), as comma-delimited string list
  std::string alternate() const 
    { 
      return consensus.alt; 
    }

  // Misc.

  void valid(bool b) { is_valid = b; } 

  bool valid() { return is_valid; } 

  bool invalid() { return ! is_valid; } 

  bool concordant( int s1, const Genotype * g1, int s2, const Genotype * g2 ) const;

  bool concordant( SampleVariant * s1, const Genotype * g1, SampleVariant * s2, const Genotype * g2 ) const;
  
  /// Printing access to meta-information 

  std::string print_samples( const std::string & delim = " " ) const;
  std::string print_meta( const std::string & key , const std::string & delim = " ") const;
  std::string print_meta_filter( const std::string & delim = " ") const;
  

  //
  // Core individual access function -- look at the indmap, *align
  //

  
  /// Return number of individuals with data for this variant
  
  int size() const;
  

  /// Pointer to the k'th individual (0..N-1) in the consensus/indmap, given N or ID
  
  Individual * ind(const int) const;
  
  Individual * ind(const std::string &) const;
  
  /// vector of all IDs of for the consensus indmap
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
  
  Genotype * genotype(const int);
  
  const Genotype * genotype(const int) const;

  /// Overload () to provide references to consensus genotypes
  
  Genotype & operator()(const int i) { return consensus.calls.genotype(i); }

  const Genotype & operator()(const int i) const { return consensus.calls.genotype(i); }

  /// All genotypes for a given individual
  std::map<int, Genotype *> all_genotype(const int);
  
  /// All genotypes for a given individual (const version)
  std::map<int,const Genotype *> all_genotype(const int) const;


  /// From a specific sample (but which will be consensus under a flat alignment)
  
  const Genotype * genotype( const SampleVariant * svar , const int i) const;
  

  //
  // Sample-level access functions
  //
  
  bool infile_overlap() const
  {
    // if more SV-ids than file ids for this variant, implies that we
    // saw at least one "merged" variant within a single file.
    
    return svtof.size() > ftosv.size();
  }

  /// Return list of all sv-ids (which index SampleVariants)
  
  std::vector<int> samples() const;
  
  /// Return set of only unique file-IDs

  std::set<int> unique_files() const ;
  
  // Return svar-slot given file ID (or -1 if no uniq mapping)
  int unique_svar_slot( int f ) const;;

  // Is this file-ID present in the Variant?
  bool file_present( const int f ) const;

  /// Set the internal iterator to the first sample

  void set_first_sample();
    
  /// Return a reference for the current sample
  
  SampleVariant & sample();

  /// Return int code to that sample

  int sample_n();

  /// Advance to next sample; return T if next sample is valid
  
  bool next_sample();
  
  
  /// Return a pointer to a given sample

  SampleVariant * psample(const int s) const 
    {
      return s < 0 || s >= svar.size() ? NULL : (SampleVariant*)&(svar[s]);
    }

  /// As above, but via external file codes
  // TODO: will need to return multiple sometimes (samefile can appear>1)
  SampleVariant * fsample(const int s) const 
    {
      std::map<int,std::vector<int> >::const_iterator i = ftosv.find(s);
      return  i == ftosv.end() ? NULL : psample( i->second[0] );
    }

  
  /*!
    Allele count for variant, as specified in ALT/REF
    @return Number of unique alleles, including reference 
  */

  int n_alleles() const;
  
  std::map<std::string,int> allele_counts( const affType & aff = UNKNOWN_PHE ) const;
  
  bool biallelic() const;
  
  bool monomorphic() const;
  
  bool multiallelic() const;

  /*!
    Is this variant a simple SNP?
    @return T if only two alleles observed and single base substitution
  */

  bool simple_snp() const;

  bool transition() const;

  bool transversion() const;

  // 
  bool simple_ins() const;

  bool simple_del() const;

  /*! 
    Queries about genotype data
    @return Number of nonreference individuals 
  */
  
  int n_nonreference() const;


  /*!
    Return number of 
    @param m Reference: returns number of minor alleles 
    @param n Reference: returns total number of non-missing alleles
    @return True if alternate allele is minor allele
  */

  bool n_minor_allele( int & m , int & n , const affType & aff = UNKNOWN_PHE ) const;

  bool n_minor_allele() const
    {
      int dummy1, dummy2;
      return n_minor_allele(dummy1, dummy2);
    } 
  
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
  std::string label( const int, const std::string & delim = "," ) const;
  std::string gmeta_label( const int, const std::string & delim = "," ) const;

  
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
  
  std::string VCF();


  std::string displaycore() const 
    { 
      std::stringstream ss;
      ss << Helper::chrCode( chr ) << ":" << bp;
      if ( bp2 != 0 && bp2 != bp ) ss << ".." << bp2;
      ss << ":" << vname;			      
      return ss.str();
    }

  
  bool operator<( const Variant & rhs ) const
  {
    if ( chr == rhs.chr )
      {
	if ( bp == rhs.bp )
	  {
	    return bp2 < rhs.bp2;
	  }
	return bp < rhs.bp;
      }
    else
      return chr < rhs.chr;
  }


  friend std::ostream & operator<<( std::ostream & out, const Variant & v)
    { 
      out << Helper::chrCode( v.chr ) << ":" << v.bp;
      if ( v.bp2 != 0 && v.bp2 != v.bp ) out << ".." << v.bp2;
      out << ":" << v.vname;
      return out;
    }
  
  bool single_sample() const 
  {
    return svar.size() == 0 ;
  }
  
  int n_samples() const 
  {
    // in single-file mode, the data is only put in the consensus, so not even a single
    // svar -- hmm. should straighten this out.
    return svar.size() == 0 ? 1 : svar.size();
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
  
  bool remove( int s ) 
    {

      // this will invalidate 'si'
      if ( s < 0 || s >= svar.size() ) return false;
      
      // remove from svar (and also svtof and ftosv)
      svar.erase( svar.begin() + s );
      
      svtof.erase( svtof.begin() + s );
      std::map<int,std::vector<int> >::iterator i = ftosv.begin();
      while ( i != ftosv.end() ) 
	{
	  std::vector<int>::iterator j = i->second.begin();
	  while ( j != i->second.end() )
	    {	      
	      if ( *j == s ) 
		j = i->second.erase( j );
	      else
		++j;
	    }	  
	  ++i;
	}      
	
      return true;
    }

  SampleVariant & sample(const int s)
    {
      return svar[ s ];
    }
  
  SampleVariant & first_sample()
    {
      return svar[0];
    }

};



#endif
  


