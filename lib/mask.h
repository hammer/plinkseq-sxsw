#ifndef __MASK_H__
#define __MASK_H__

#include <set>

#include "regions.h"
#include "variant.h"
#include "eval.h"
#include "defs.h"

#include "locdb.h"
#include "refdb.h"

class VarDBase;
class LocDBase;
class RefDBase;

struct mask_command_t 
{ 
  mask_command_t(const std::string & n, 
		 const int no = 0 , 
		 const std::string & g = "", 
		 const int go = 0 , 
		 const std::string & a = "", 
		 const std::string & d = "", 
		 const bool h = false ) 
  : name(n) , name_order(no) , group(g), group_order(go) , argtype(a) , desc(d) , hidden(h) 
  { } 
  
  std::string name;
  std::string group;
  int group_order;
  int name_order;
  std::string desc;
  std::string argtype;
  bool hidden;
  
  bool operator<( const mask_command_t & rhs ) const 
  { 
    if ( group == "" || rhs.group == "" ) 
      return name_order < rhs.name_order;       
    if ( group_order < rhs.group_order ) return true;
    if ( group_order > rhs.group_order ) return false;    
    return name_order < rhs.name_order;       
  }
};


typedef bool(*mask_func_t)(Variant &, void *);

/*!
  
  This class defines which variants, individuals and genotypes are
  pulled out for a given iteration over the data. A Mask can
  contain inclusions, requirements or exclusions based on genomic position, 
  variant-database groups, loci, locus-sets, user-defined functions, 
  variant annotation, file, etc

  @class Mask
  
  @brief Specify which variants are passed to the user's function, 
    as well as how they are grouped and ordered.
*/

class Mask {
    

 public:    
  
  Mask( const std::string & , 
	const std::string & expr = "" , 
	const bool finclude = true , 
	bool groups = false );
  
  Mask(VarDBase * vardb = NULL, 
       LocDBase * locdb = NULL, 
       LocDBase * segdb = NULL, 
       RefDBase * refdb = NULL )
    : vardb(vardb) , locdb(locdb), segdb(segdb), refdb(refdb)
    { 	    
      construct();
    }
    
  int include_loc( int x );    
  int require_loc( int x );    
  int exclude_loc( int x );    
    
  int include_var( int x );   
  int require_var( int x );   
  int exclude_var( int x );   
  
  int include_ref( int x );   
  int require_ref( int x );   
  int exclude_ref( int x );   

  int include_loc( std::string n );
  int require_loc( std::string n );
  int exclude_loc( std::string n );
  
  int include_var( std::string n );
  int require_var( std::string n );
  int exclude_var( std::string n );
  
  int include_ref( std::string n );
  int require_ref( std::string n );
  int exclude_ref( std::string n );
  
  void individual_data(const bool b) { inddata = b; if (b) is_simple=false; }
  bool individual_data() const { return inddata; }
  
  bool loc_seg(const int i ) const 
  { 
    std::map<int,bool>::const_iterator f = loc_or_seg.find( i );
    if ( f != loc_or_seg.end() ) return f->second;
    return false;
  }
  
  bool locset_segset(const int i ) const 
  { 
    std::map<int,bool>::const_iterator f = locset_or_segset.find( i );
    if ( f != locset_or_segset.end() ) return f->second;
    return false;
  }
  
  
  void loc_seg(const int i , bool b) 
  { 
    loc_or_seg[i] = b; 
  }
  
  void locset_segset(const int i , bool b) 
  { 
    locset_or_segset[i] = b; 
  }
  
  const std::set<int> & included_loc() const { return in_locset; }
  const std::set<int> & included_var() { return in_varset; }
  const std::set<Region> & included_reg() const { return in_regions; }
  const std::set<int> & included_ref() const { return in_refset; }

  const std::set<int> & required_loc() const { return req_locset; }
  const std::set<int> & required_locset() const { return req_locset_set; }
  const std::set<int> & required_var() const { return req_varset; }
  const std::set<Region> & required_reg() const { return req_regions; }
  const std::set<int> & required_ref() const { return req_refset; }

  const std::set<int> & excluded_loc() const { return ex_locset; }
  const std::set<int> & excluded_locset() const { return ex_locset_set; }
  const std::set<int> & excluded_var() const { return ex_varset; }
  const std::set<Region> & excluded_reg() const { return ex_regions; }
  const std::set<int> & excluded_ref() const { return ex_refset; }
  
  const std::set<int> appended_ref() const { return app_refset; }
  const std::set<int> appended_var() const { return app_varset; }
  const std::set<int> appended_loc() const { return app_locset; }
  const std::set<int> appended_loc_set() const { return app_locset_set; }
  
  
  // Sets within locdb, or segdb
  
  int include_loc_set( int x );
  int include_loc_set( std::string n , std::string p );
  std::set<int> & included_loc_set() { return in_locset_set; }

  int require_loc_set( int x );
  int require_loc_set( std::string n , std::string p );
  std::set<int> & required_loc_set() { return req_locset_set; }

  int exclude_loc_set( int x );
  int exclude_loc_set( std::string n , std::string p );
  std::set<int> & excluded_loc_set() { return ex_locset_set; }


  //
  // From the broad super-sets specified above, we can elect to
  //
  //  --> only extract certain sets (subset)
  //  --> ignore certain sets (skip)
  
  void subset_loc(const int, const std::string&);
  void subset_loc(const int, const std::vector<std::string>&);
  void skip_loc(const int, const std::string&);
  void skip_loc(const int, const std::vector<std::string>&);
  
  void subset_var(const int, const std::string&);
  void subset_var(const int, const std::vector<std::string>&);
  void skip_var(const int, const std::string&);
  void skip_var(const int, const std::vector<std::string>&);
  
  // As above, but by std::string

  void subset_loc(const std::string & , const std::string&);
  void subset_loc(const std::string & , const std::vector<std::string>&);
  void skip_loc(const std::string & , const std::string&);
  void skip_loc(const std::string & , const std::vector<std::string>&);
  
  void subset_var(const std::string & , const std::string&);
  void subset_var(const std::string & , const std::vector<std::string>&);
  void skip_var(const std::string & , const std::string&);
  void skip_var(const std::string & , const std::vector<std::string>&);
  
  void subset_loc_set(const std::string & , const std::string &, const std::string & );
  void skip_loc_set(const std::string & , const std::string &, const std::string & );
  
  // search by 'alternate name' (i.e. typically gene symbol, where 'name' is transcript ID)
  // and so alternate name is not assumed to be unique
  void subset_loc_altname(const std::string & , const std::string&);
  void subset_loc_altname(const std::string & , const std::vector<std::string>&);
    

  //
  // User-defined regions
  //

  void include_reg( const Region & r ) {  is_simple = false; in_regions.insert(r); }
  void include_reg( const std::vector<std::string> & r );
  void require_reg( const Region & r ) {  is_simple = false; req_regions.insert(r); }
  void exclude_reg( const Region & r ) {  is_simple = false; ex_regions.insert(r); }

  //
  // Variant/file masks (i.e. remove variant if not observed/alt-allele in files x,y,z etc
  //

  void include_var_obs_file( const std::vector<std::string> & );
  void require_var_obs_file( const std::vector<std::string> & );
  void exclude_var_obs_file( const std::vector<std::string> & );
    
  void include_var_alt_file( const std::vector<std::string> & );
  void require_var_alt_file( const std::vector<std::string> & );
  void exclude_var_alt_file( const std::vector<std::string> & );
  
  void include_var_alt_group( const std::vector<std::string> & );
  void require_var_alt_group( const std::vector<std::string> & );
  void exclude_var_alt_group( const std::vector<std::string> & );

  void count_obs_file(int i , int j = 0 );
  void count_alt_file(int i , int j = 0 );

  bool eval_obs_file_filter( const Variant & ) const;
  bool eval_alt_file_filter( Variant & ) const;
  bool eval_alt_group_filter( const Variant & ) const;
  bool eval_file_count( Variant & ) const;  


  //
  // Variant 'FILTER' field
  //

  // 1 of ...
  void include_filter( const std::string & s );
  void include_filter( const std::vector<std::string> & s ); 
  
  // all of ...
  void require_filter( const std::string & s );
  void require_filter( const std::vector<std::string> & s ); 

  // none of ...
  void exclude_filter( const std::string & s );
  void exclude_filter( const std::vector<std::string> & s ); 

  // 
  void include_filter_any();
  void require_filter_any();
  void exclude_filter_any();

  bool eval_filters( const SampleVariant & );

  bool polymorphism_filter( const Variant & );

  //
  // User-defined functions, applied after the variant is extracted
  //

  void func( const mask_func_t f ) { is_simple = false; filterFunctions.insert(f); } 
  void include_func( const mask_func_t f ) { is_simple = false; filterFunctions.insert(f); }
  void require_func( const mask_func_t f ) { is_simple = false; req_filterFunctions.insert(f); } 

  bool eval(Variant & v, void * p = NULL);
  
  //
  // Restriction of maximum number of variants that can be returned (safety-net for R)
  //

  void limit(const int i) { is_simple = false; max_var_count = i; } 
  
  //
  // How to handle multiple sites per genomic position
  //
  
  bool fail_on_sample_variant() const ;

  bool test_fail_on_sample_variant(int n , int m ) const ;
  
  void fail_on_sample_variant( int n , int m );


  //
  // If a problem occurred by a routine trying to make a mask, then
  // set this flag
  //

  bool invalid() const { return invalid_mask; } 
  void invalid(const bool b) { invalid_mask = b; } 


  //
  // Inidividuals
  //

  void include_indiv( const std::string & id );
  void include_indiv( const std::vector<std::string> & id );
  
  void exclude_indiv( const std::string & id );
  void exclude_indiv( const std::vector<std::string> & id );
  
  bool use_indiv( const std::string & ) const;

  bool merge_indiv() const { return merge_ind; }
  void merge_indiv(const bool b) { merge_ind = b; }

  bool pheno_screen( Individual * ) const;
  
  void include_phenotype( const std::vector< std::string > & );
  void require_phenotype( const std::vector< std::string > & );
  void exclude_phenotype( const std::vector< std::string > & );
  void require_nonmissing_phenotype( const std::vector< std::string> & );
  

  //
  // Per-individual segments for SEGDB
  //

  void include_indiv_segment( const uint64_t & );
  void require_indiv_segment( const uint64_t & );
  void exclude_indiv_segment( const uint64_t & );

  void include_indiv_segment( const std::string & );
  void require_indiv_segment( const std::string & );
  void exclude_indiv_segment( const std::string & );
  

  //
  // Appends
  //
    
  int append_loc( int x );
  int append_loc( const std::string & n );
  
  int append_loc_set( int x );
  int append_loc_set( const std::string & n , const std::string & m );
  
  int append_var( int x );      
  int append_var( const std::string & n ); 

  int append_ref( int x );	
  int append_ref( const std::string & n ); 

    
  //
  // Groups (allowed a single group, either variant-set or locus set)
  //
  
  void group_var(const int g);

  void group_var(const std::string & g);

  void group_loc(const int g);

  void group_loc(const std::string & g);    
  
  void group_reg(const std::vector<std::string> & g);

  void group_loc_set(const int g);

  void group_loc_set(const std::string & g, const std::string & h);    

  void process_empty_groups(bool b) { is_simple = false; empty_groups = b; }

  bool process_empty_groups() const { return empty_groups; }

  int group_set() const 
    { 
      return group_var() ? 
	group_variant : 
	group_loc_set() ? group_locus_set : group_locus ; 
    }
  
  
  bool any_grouping() const 
  { 
    return named_grouping(); //  || file_grouping() || file_split() ; 
  } 
  
  bool named_grouping() const 
  { 
    return group_var() || group_loc() || group_loc_set() ; 
  } 

  bool group_var() const { return group_variant != 0; }

  bool group_loc() const { return group_locus != 0; }

  bool group_loc_set() const { return group_locus_set != 0; }
  
  bool group_reg() const { return group_region; }

  void ensure_single_include_group() 
    {
      if ( ! named_grouping() ) return;
      if ( group_loc() ) 
	{
	  in_locset.clear();
	  in_locset.insert( group_set() );
	}
      else if ( group_var() ) 
	{
	    in_varset.clear();
	    in_varset.insert( group_set() );
	}
      else if ( group_loc_set() )
	{
	  in_locset_set.clear();
	  in_locset_set.insert( group_set() );
	}
    }
  
  //
  // Simple file-based inclusions, exclusions
  //
  
  int include_file( const std::string & filetag );
  int exclude_file( const std::string & filetag );

  bool use_file( const int ) const;

  
  //
  // Are we looking at a single VCF?
  //

  bool external_vcf_iteration() const { return ext_vcf; }
  void external_vcf_iteration( const std::string & );
  std::string external_vcf_filename() const { return ext_vcffile; } 

  //
  // User-defined meta-variant info.
  //
  
  
  bool attach_meta() const { return will_attach_meta; }
  void attach_meta(const bool b) { will_attach_meta = b; }
  
  bool attach_all_meta() const { return will_attach_all_meta; }
  void attach_all_meta(const bool b) { will_attach_all_meta = b; if (b) attach_meta(true); }

  std::set<std::string> attach_meta_fields() const { return meta_fields; }
  void attach_meta_fields(const std::set<std::string> & s ) { meta_fields = s; attach_meta(true); }
  void attach_meta_fields(const std::string & s ) { meta_fields.insert(s); attach_meta(true); }
  
  //
  // Simple/biallelic variants versus not
  //

  bool require_biallelic() const { return req_biallelic; }
  void require_biallelic(const bool b ) { req_biallelic = b; }

  bool exclude_biallelic() const { return exc_biallelic; }
  void exclude_biallelic(const bool b ) { exc_biallelic = b; }

  bool require_monomorphic() const { return req_monomorphic; }
  void require_monomorphic(const bool b ) { req_monomorphic = b; }

  bool exclude_monomorphic() const { return exc_monomorphic; }
  void exclude_monomorphic(const bool b ) { exc_monomorphic = b; }
  
  // 
  // QUAL filter
  //

  void qual_filter( const std::string & s ) { qual.set(s); use_qual_filter = true; }

  //
  // Variant meta-information masks
  //

  int meta_equals( const std::string & key , int value );
  int meta_not_equals( const std::string & key , int value );

  int meta_equals( const std::string & key , const std::string & value );
  int meta_not_equals( const std::string & key , const std::string & value );
  
  int meta_greater( const std::string & key , double value );
  int meta_greater_equal( const std::string & key , double value );

  int meta_less( const std::string & key , double value );
  int meta_less_equal( const std::string & key , double value );
  
  // as above, for requires

  int req_meta_equals( const std::string & key , int value );
  int req_meta_not_equals( const std::string & key , int value );

  int req_meta_equals( const std::string & key , const std::string & value );
  int req_meta_not_equals( const std::string & key , const std::string & value );
  
  int req_meta_greater( const std::string & key , double value );
  int req_meta_greater_equal( const std::string & key , double value );

  int req_meta_less( const std::string & key , double value );
  int req_meta_less_equal( const std::string & key , double value );

  bool meta_includes() const
    {
      return meta_eq.size() > 0 || 
	meta_ne.size() > 0 || 
	meta_has_text.size() > 0 || 
	meta_has_not_text.size() > 0 || 
	meta_gt.size() > 0 || 
	meta_ge.size() > 0 || 
	meta_lt.size() > 0 || 
	meta_le.size() > 0 ;
    }

  bool meta_requires() const
    {
      return req_meta_eq.size() > 0 || 
	req_meta_ne.size() > 0 || 
	req_meta_has_text.size() > 0 || 
	req_meta_has_not_text.size() > 0 || 
	req_meta_gt.size() > 0 || 
	req_meta_ge.size() > 0 || 
	req_meta_lt.size() > 0 || 
	req_meta_le.size() > 0 ;
    }



  void determine_variant_mask()
    {
      meta_mask = meta_includes() || meta_requires() || eval_expr_set;
    }
  
  bool variant_mask() const { return meta_mask; }

  bool eval( SampleVariant & svar );



  //
  // Genotype masks
  //

  int geno_equals( const std::string & key , int value );
  int geno_not_equals( const std::string & key , int value );

  int geno_equals( const std::string & key , const std::string & value );
  int geno_not_equals( const std::string & key , const std::string & value );
  
  int geno_greater( const std::string & key , double value );
  int geno_greater_equal( const std::string & key , double value );

  int geno_less( const std::string & key , double value );
  int geno_less_equal( const std::string & key , double value );
  
  // as above, for requires

  int req_geno_equals( const std::string & key , int value );
  int req_geno_not_equals( const std::string & key , int value );

  int req_geno_equals( const std::string & key , const std::string & value );
  int req_geno_not_equals( const std::string & key , const std::string & value );
  
  int req_geno_greater( const std::string & key , double value );
  int req_geno_greater_equal( const std::string & key , double value );

  int req_geno_less( const std::string & key , double value );
  int req_geno_less_equal( const std::string & key , double value );

  bool geno_includes() const
    {
      return geno_eq.size() > 0 || 
	geno_ne.size() > 0 || 
	geno_has_text.size() > 0 || 
	geno_has_not_text.size() > 0 || 
	geno_gt.size() > 0 || 
	geno_ge.size() > 0 || 
	geno_lt.size() > 0 || 
	geno_le.size() > 0 ;
    }

  bool geno_requires() const
    {
      return req_geno_eq.size() > 0 || 
	req_geno_ne.size() > 0 || 
	req_geno_has_text.size() > 0 || 
	req_geno_has_not_text.size() > 0 || 
	req_geno_gt.size() > 0 || 
	req_geno_ge.size() > 0 || 
	req_geno_lt.size() > 0 || 
	req_geno_le.size() > 0 ;
    }
  
  
  
  void determine_genotype_mask()
    {
      geno_mask = geno_includes() || geno_requires() ;
    }
  
  bool genotype_mask() const { return geno_mask; }

  bool eval( const Genotype & g ) const;


  //
  // Auto-conversion of Null to Reference
  //

  void assuming_null_is_reference( const bool b )
  { is_simple = false; assume_missing_is_ref = b; } 
  
  bool assuming_null_is_reference() const
  { return assume_missing_is_ref; } 
  

  //
  // Variant merging and allele-downcoding
  //
  
  bool exact_merge() const { return exact_vmerge; } 
  void exact_merge( const bool b ) { exact_vmerge = b; }

  downcode_mode_t downcode() const { return downcode_mode; } 
  void downcode( const downcode_mode_t d ) { downcode_mode = d; } 

  

  //
  // Invoke EM caller
  //
  
  bool EM_caller() const { return use_em; }
  void EM_caller( const bool b ) { is_simple = false; use_em = b; }
  
  bool EM_replace() const { return em_replace; }
  void EM_replace( const bool b ) { em_replace = b; }
  
  double EM_threshold() const { return em_threshold; } 
  void EM_threshold( const double t ) { em_threshold = t; }
  
  
  //
  // Allele-frequency based variant filters
  //

  bool count_filter() const { return mac_filter; }
  bool frequency_filter() const { return maf_filter; }

  bool count_filter(const int i) const 
  {       
    // if mac threshold is -1, means ignore
    if ( mac_lower >= 0 && i < mac_lower ) return false;
    if ( mac_upper >= 0 && i > mac_upper ) return false;
    return true;
  }
  
  bool frequency_filter(const double f) const 
    { 
      // if maf threshold is -1, means ignore
      if ( maf_lower >= 0 && f < maf_lower ) return false;
      if ( maf_upper >= 0 && f > maf_upper ) return false;
      return true;      
    }

  void minor_allele_count(const int c, const int d ) 
  { 
    is_simple = false; 
    mac_filter = true;
    mac_lower = c; 
    mac_upper = d;
  }
  
  void minor_allele_frequency(const double c, const double d ) 
  { 
    is_simple = false; 
    maf_filter = true;
    maf_lower = c; 
    maf_upper = d;
  }
  
  
  void hwe( double l, double u )
  {
    is_simple = false;
    use_hwe_filter = true;
    hwe_lower = l;
    hwe_upper = u;
  }
  
  bool hwe_filter( const double p ) 
  {
    // if not defined, do not exclude (e.g. monomorphic)
    if ( !Helper::realnum( p ) ) return true;
    return p >= hwe_lower && p <= hwe_upper;
  }
  
  bool hwe_filter() 
  {
    return use_hwe_filter;
  }
  
  bool get_minor_allele_count(int & c, int & d ) 
  {
    if ( ! mac_filter ) return false;
    c = mac_lower; d = mac_upper;
    return true;
  }
  
  bool get_minor_allele_frequency(double & c, double & d ) 
  {
    if ( ! maf_filter ) return false;
    c = maf_lower; d = maf_upper;
    return true;
  }
  

  //
  // Null genotype filter
  //

  void null_filter( const int_range & r );
  bool null_filter( ) const;
  bool null_filter( const int ) const;


  //
  // Case/control counts
  //

  void case_control_filter( const std::string & a , const std::string & u );
  bool case_control_filter( ) const;
  bool case_control_filter( const int , const int ) const;

  //
  // Sequence annotations
  //
  
  void include_annotation_nonsyn();
  void require_annotation_nonsyn();
  void exclude_annotation_nonsyn();
  
  void include_annotation( const std::string & );
  void include_annotation( const std::vector<std::string> & );
  
  void require_annotation( const std::string & );
  void require_annotation( const std::vector<std::string> & );
  
  void exclude_annotation( const std::string & );
  void exclude_annotation( const std::vector<std::string> & );
  
  void append_annotation() { is_simple = false; annot = true; }
  
  //
  // Get a dump of all mask options
  //

  static std::string describe_options();     

  //
  // Helper functions
  //
        
  void reset()
    {
      
      is_simple = true;

      inddata = true;

      in_locset.clear();
      in_varset.clear();
      in_refset.clear();
      
   
      // LOCDB or SEGDB? (for normal loci, and locus-sets)
      loc_or_seg.clear();
      locset_or_segset.clear();
      
      req_locset.clear();
      req_varset.clear();
      req_refset.clear();

      ex_locset.clear();
      ex_varset.clear();
      ex_refset.clear();

      subset_varset.clear();
      subset_locset.clear();
      subset_locset_set.clear();

      skip_varset.clear();
      skip_locset.clear();
      skip_locset_set.clear();
      
      // Var/file masks

      obs_file_filter = false;
      obs_file_count = 0;
      alt_file_count = 0;
      obs_file_max = 0;
      alt_file_max = 0;

      alt_file_filter = false;
      alt_group_filter = false;
      
      inc_obs_file.clear();
      req_obs_file.clear();
      exc_obs_file.clear();
      
      inc_alt_file.clear();
      req_alt_file.clear();
      exc_alt_file.clear();
      
      inc_alt_group.clear();
      req_alt_group.clear();
      exc_alt_group.clear();

      // Filters
      
      inc_filter.clear();
      req_filter.clear();
      exc_filter.clear();

      inc_filter_any = false;
      req_filter_any = false;
      exc_filter_any = false;

      // var Meta-information 

      meta_eq.clear();
      meta_ne.clear();
      
      meta_has_text.clear();
      meta_has_not_text.clear();
      
      meta_gt.clear();
      meta_ge.clear();
      
      meta_lt.clear();
      meta_le.clear();


      req_meta_eq.clear();
      req_meta_ne.clear();
      
      req_meta_has_text.clear();
      req_meta_has_not_text.clear();
      
      req_meta_gt.clear();
      req_meta_ge.clear();
      
      req_meta_lt.clear();
      req_meta_le.clear();

      // Genotypes

      geno_eq.clear();
      geno_ne.clear();
      
      geno_has_text.clear();
      geno_has_not_text.clear();
      
      geno_gt.clear();
      geno_ge.clear();
      
      geno_lt.clear();
      geno_le.clear();


      req_geno_eq.clear();
      req_geno_ne.clear();
      
      req_geno_has_text.clear();
      req_geno_has_not_text.clear();
      
      req_geno_gt.clear();
      req_geno_ge.clear();
      
      req_geno_lt.clear();
      req_geno_le.clear();

      // variant/allele merging/splitting

      exact_vmerge = true;
      downcode_mode = DOWNCODE_MODE_ALL_ALT;
      
      // EM caller
      
      use_em = false;
      em_replace = false;
      em_threshold = 0;

      // Files
      
      in_files.clear();
      ex_files.clear();
      
      // External VCF mode?
      
      ext_vcf = false;
      ext_vcffile = "";

      // Individuals
      
      in_indset.clear();
      ex_indset.clear();

      // Phenotypes
      
      in_phe.clear();
      req_phe.clear();
      ex_phe.clear();

      
      // User-specified regions

      in_regions.clear();
      req_regions.clear();
      ex_regions.clear();

      // Grouping
	    
      group_region = false;
      group_locus = group_variant = group_locus_set = 0;

      // Functions

      filterFunctions.clear();
      req_filterFunctions.clear();


      // Appends
      
      app_locset.clear();
      app_varset.clear();
      app_refset.clear();

      // Locus-sets (and parallele SEGDB)
	    
      in_locset_set.clear();
      req_locset_set.clear();
      ex_locset_set.clear();
      app_locset_set.clear();
      
      // Allele frequency filters
      // either absolute counts 
      // or frequencies (that respect missing data)
      
      mac_filter = false;
      maf_filter = false;
      
      has_null_filter = false;
      has_case_control_filter = false;
	
    }
  

  bool simple( const int n_files = 2 ) const
  {
    if ( ! is_simple ) return false;

    // is this next line needed now??
    if ( n_files > 1 ) return false; 

    return true;
  }
  
  
  bool loc() const { return in_locset.size() > 0 ; }
  bool var() const { return in_varset.size() > 0 ; }
  bool reg() const { return in_regions.size() > 0 ; }
  
  bool ref() const { return in_refset.size() > 0; }
  bool loc_set() const { return in_locset_set.size() > 0; }

  bool rloc() const { return req_locset.size() > 0; }
  bool rvar() const { return req_varset.size() > 0; }
  bool rreg() const { return req_regions.size() > 0; }
  bool rref() const { return req_refset.size() > 0; }
  
  bool xloc() const { return ex_locset.size() > 0; }
  bool xvar() const { return ex_varset.size() > 0; }
  bool xreg() const { return ex_regions.size() > 0; }
  bool xref() const { return ex_refset.size() > 0; }
  
  bool requires() const 
  { 
    return req_locset.size() > 0
      || req_varset.size() > 0
      || req_regions.size() > 0 ;
      }
  
  bool excludes() const { return ex_locset.size() > 0 
			    || ex_varset.size() > 0 
			    || ex_regions.size() > 0 ;  }
  
  
  bool loc_exceptions() const { return subset_locset.size() > 0 
				  || skip_locset.size() > 0; } 
  
/*   bool xloc_exceptions() const { return subset_xlocset.size() > 0  */
/* 				   || skip_xlocset.size() > 0; }  */
  
  bool var_exceptions() const { return subset_varset.size() > 0 
				  || skip_varset.size() > 0; } 
  
/*   bool xvar_exceptions() const { return subset_xvarset.size() > 0  */
/* 				   || skip_xvarset.size() > 0; }  */
  
  bool locset_exceptions() const { return subset_locset_set.size() > 0 
      || skip_locset_set.size() > 0 ; } 
  
  bool func() const 
    { 
      return filterFunctions.size() > 0 
	|| req_filterFunctions.size() > 0 
	|| annot;
    } 
  
  //    bool rfunc() const { return RfilterFunctions.size() > 0; } 
  
  bool files() const { return in_files.size() > 0; }
  bool xfiles() const { return ex_files.size() > 0; }
  
  
  bool build_temporary_db() const
    {
      return var_exceptions() 
	|| named_grouping() 
	|| loc_exceptions() 
	|| locset_exceptions()
	|| loc_set() 
	|| loc_set_append() 
	|| requires() 
	|| excludes() 
	|| reg() 
	|| xreg() 
	|| files() 
	|| xfiles();
    }
  
  
  bool var_append() const { return app_varset.size() > 0; }
  bool loc_append() const { return app_locset.size() > 0; }
  bool loc_set_append() const { return app_locset_set.size() > 0; }
  bool ref_append() const { return app_refset.size() > 0; }
  
  bool append() const 
    { return var_append() || loc_append() || ref_append(); }
  
  int loc_size() const { return in_locset.size(); }
  int var_size() const { return in_varset.size(); }
  
  int variant_limit() const { return max_var_count; }
  
  
  bool loc_any() const { return loc() || loc_append() || rloc() || xloc() || loc_set(); }
  bool var_any() const { return var() || var_append() || rvar() || xvar(); }
  bool ref_any() const { return ref_append(); }
  bool reg_any() const { return reg() || rreg() || xreg(); } 
  bool loc_set_any() const { return loc_set(); }
  

  std::set<std::string> subset_loc(int j) const
    {       
      std::map<int, std::set<std::string> >::const_iterator i = subset_locset.find(j);
      if ( i != subset_locset.end() ) return i->second;
      std::set<std::string> t;      
      return t;
    }
  
  std::set<std::string> subset_loc_set(int j) const
    {       
      std::map<int, std::set<std::string> >::const_iterator i = subset_locset_set.find(j);
      if ( i != subset_locset_set.end() ) return i->second;      
      std::set<std::string> t;
      return t;
    }

  bool insert_locset( const int j , const std::string & n ) const;


  std::set<std::string> subset_var(int j) const
    {       
      std::map<int, std::set<std::string> >::const_iterator i = subset_varset.find(j);
      if ( i != subset_varset.end() ) return i->second;
      std::set<std::string> t;
      return t;
    }
  
    ///

    std::set<std::string> skip_loc(int j) const
      { 
	
	std::map<int, std::set<std::string> >::const_iterator i = skip_locset.find(j);
	if ( i != skip_locset.end() ) return i->second;
	std::set<std::string> t;
	return t;
      }
    
    std::set<std::string> skip_loc_set(int j) const
      { 	
	std::map<int, std::set<std::string> >::const_iterator i = skip_locset.find(j);
	if ( i != skip_locset.end() ) return i->second;
	std::set<std::string> t;
	return t;
      }
    
    std::set<std::string> skip_var(int j) const
      { 	
	std::map<int, std::set<std::string> >::const_iterator i = skip_varset.find(j);
	if ( i != skip_varset.end() ) return i->second;
	std::set<std::string> t;
	return t;
      }
    

    ///

/*     std::set<std::string> subset_xloc(int j) const */
/* 	{  */
/* 	    std::set<std::string> t; */
/* 	    std::map<int, std::set<std::string> >::const_iterator i = subset_xlocset.find(j); */
/* 	    if ( i != subset_xlocset.end() ) */
/* 		t = i->second; */
/* 	    return t; */
/* 	} */

/*     std::set<std::string> subset_xvar(int j) const */
/* 	{  */
/* 	    std::set<std::string> t; */
/* 	    std::map<int, std::set<std::string> >::const_iterator i = subset_xvarset.find(j); */
/* 	    if ( i != subset_xvarset.end() ) */
/* 		t = i->second; */
/* 	    return t; */
/* 	} */

///

/*     std::set<std::string> skip_xloc(int j) const */
/* 	{  */
/* 	    std::set<std::string> t; */
/* 	    std::map<int, std::set<std::string> >::const_iterator i = skip_xlocset.find(j); */
/* 	    if ( i != skip_xlocset.end() ) */
/* 		t = i->second; */
/* 	    return t; */
/* 	} */

/*     std::set<std::string> skip_xvar(int j) const */
/* 	{  */
/* 	    std::set<std::string> t; */
/* 	    std::map<int, std::set<std::string> >::const_iterator i = skip_xvarset.find(j); */
/* 	    if ( i != skip_xvarset.end() ) */
/* 		t = i->second; */
/* 	    return t; */
/* 	} */

///

    std::string set2str(const std::set<int> & s) const
      {
	std::string t="";
	std::set<int>::iterator j = s.begin();
	bool first = true;
	while ( j != s.end() )
	  {
	    if ( ! first ) t += ", ";
	    first = false;
	    t += Helper::int2str( *j );
	    ++j;
	  }	
	return t;
      }
    

    std::string var_include_string() const
      { 
	return set2str( in_varset );
      }

    std::string var_require_string() const
      { 
	return set2str( req_varset );
      }

    std::string var_exclude_string() const
      { 
	return set2str( ex_varset );
      }


    std::string loc_include_string() const
      { 
	return set2str( in_locset );
      }

    std::string loc_require_string() const
      { 
	return set2str( req_locset );
      }

    std::string loc_exclude_string() const
      { 
	return set2str( ex_locset );
      }
    
    
    std::string loc_set_include_string() const
      {
	return set2str( in_locset_set );
      }


    std::string var_append_string() const
      { 
	return set2str( app_varset );
      }

    std::string loc_append_string() const
      { 
	return set2str( app_locset );	
      }

    std::string loc_set_append_string() const
      { 
	return set2str( app_locset_set );	
      }	

    std::string ref_append_string() const
      { 
	return set2str( app_refset );	
      }	


    std::string files_include_string() const
      {
	return set2str( in_files );		
      }
    
    std::string files_exclude_string() const
      {
	return set2str( ex_files );		
      }	


    std::string locdb_name() const { return locdb ? locdb->filename() : ""; }
    std::string segdb_name() const { return segdb ? segdb->filename() : ""; }
    std::string refdb_name() const { return refdb ? refdb->filename() : ""; }
    RefDBase * refdb_pointer() const { return refdb; }
    VarDBase * vardb_pointer() const { return vardb; }
    LocDBase * locdb_pointer() const { return locdb; }
    LocDBase * segdb_pointer() const { return segdb; }

 private:
    
    void construct()
      {
	searchDB();
	group_variant = group_locus = group_locus_set = 0;
	group_region = false;
	group_mode = false; 
	empty_groups = false;
	annot = false;
	max_var_count = 0; // no limit for number of variants processed
	invalid_mask = false; // used by R API
	req_nonmissing_phenotype.clear(); 
	merge_ind = true;
	meta_mask = false;
	geno_mask = false;
	assume_missing_is_ref = false;
	mac_filter = false;
	maf_filter = false;
	has_null_filter = false;
	has_case_control_filter = false;
	inc_filter_any = false;
	req_filter_any = false;
	exc_filter_any = false;
	mac_lower = mac_upper = -1;
	maf_lower = maf_upper = -1.0;
	use_hwe_filter = false;
	hwe_lower = 0;
	hwe_upper = 1;
	use_qual_filter = false;
	qual.set(".");
	null_fltr.reset();
	case_fltr.reset();
	control_fltr.reset();
	is_simple = true;
	fail_on_sample_variant_allow = -1; // allow inf.
	fail_on_sample_variant_require = -1; // no req.
	use_em = false;
	em_replace = false;
	em_threshold = 0;
	will_attach_meta = false;
	will_attach_all_meta = false;
	req_biallelic = false;
	exc_biallelic = false;
	req_monomorphic = false;
	exc_monomorphic = false;
	eval_expr_set = false;
	eval_filter_includes = true;
	var_eval_expr_set = false;
	var_eval_filter_includes = true;
	ext_vcf = false;
	ext_vcffile = "";
	obs_file_filter = false;
	obs_file_count = 0;
	obs_file_max = 0;
	alt_file_filter = false;
	alt_file_count = 0;
	alt_file_max = 0;
	alt_group_filter = false;
      } 
    
    void searchDB();
    
    static std::set<mask_command_t> known_commands; // mask commands
    static std::set<std::string> known_commands_str; // simple str duplicate
           
    bool group_mode;  // will the mask be used in a group-iteration context?
    bool group_region;

    VarDBase * vardb;
    LocDBase * locdb;
    LocDBase * segdb;
    RefDBase * refdb;

    // Any data-base specified groups (regions, variant sets, or positions (from ref-db)
    // These are used in constructing the original SQL query
 
    std::set<int> in_locset;
    std::set<int> req_locset;
    std::set<int> ex_locset;  

    std::set<int> in_varset; 
    std::set<int> req_varset; 
    std::set<int> ex_varset; 

    std::set<int> in_refset; 
    std::set<int> req_refset; 
    std::set<int> ex_refset; 

    // Locus-sets

    std::set<int> in_locset_set;
    std::set<int> req_locset_set;
    std::set<int> ex_locset_set;
    std::set<int> app_locset_set;
    
    // Indicator of whether a group should be taken from 
    // LOCDB or SEGDB

    std::map<int,bool> loc_or_seg;
    std::map<int,bool> locset_or_segset;

    // Var/file masks
    
    std::set<int> inc_obs_file;
    std::set<int> req_obs_file;
    std::set<int> exc_obs_file;
      
    std::set<int> inc_alt_file;
    std::set<int> req_alt_file;
    std::set<int> exc_alt_file;

    std::set<int> inc_alt_group;
    std::set<int> req_alt_group;
    std::set<int> exc_alt_group;

    bool obs_file_filter;
    int  obs_file_count;
    int  obs_file_max;
    bool alt_file_filter;
    int  alt_file_count;
    int  alt_file_max;
    bool alt_group_filter;
    

    // Filter-strings

    std::set<std::string> inc_filter;
    std::set<std::string> req_filter;
    std::set<std::string> exc_filter;

    bool inc_filter_any;
    bool req_filter_any;
    bool exc_filter_any;


    // Groups

    int group_variant;
    int group_locus;  // specifies either LOCDB or SEGDB
    int group_locus_set;
    bool empty_groups;

    // From the above sets, are we only interested in a named subset of
    // entries, that have been specified directly to this Mask, i.e. not
    // explicitly represented in the database as of yet?  These can be 
    // both inclusion and exlcusion parameters. The keys here should match
    // the return value from the insertion
    
    std::map<int,std::set<std::string> > subset_locset;
    std::map<int,std::set<std::string> > subset_locset_set;
    std::map<int,std::set<std::string> > subset_varset;

    std::map<int,std::set<std::string> > skip_locset;
    std::map<int,std::set<std::string> > skip_locset_set;
    std::map<int,std::set<std::string> > skip_varset;
    
    // Files

    std::set<int> in_files;
    std::set<int> ex_files;


    // External VCF mode

    bool ext_vcf;
    std::string ext_vcffile;

    // User-defined data, not in database
    
    //
    // Regions:
    //

    std::set<Region> in_regions;
    std::set<Region> req_regions;
    std::set<Region> ex_regions;

    
    // Functions that can be applied to the full variant, after extracting
    // from the database but before sending on the the calling function. All
    // standard attributes and meta-information for the variant will be
    // available.
    
    std::set< mask_func_t > filterFunctions  ; 
    std::set< mask_func_t > req_filterFunctions  ; 

    //
    // Misc. valeus
    //
    
    int max_var_count;
    
    bool invalid_mask;
    
    bool is_simple;

    bool inddata;

    int fail_on_sample_variant_allow;
    int fail_on_sample_variant_require;

    //
    // Individuals
    //

    std::set<std::string> in_indset;
    std::set<std::string> ex_indset;
    bool merge_ind;

    //
    // Phenotypes
    //
    
    std::map< std::string , std::set<std::string> > in_phe;
    std::map< std::string , std::set<std::string> > req_phe;
    std::map< std::string , std::set<std::string> > ex_phe;
    std::set<std::string> req_nonmissing_phenotype;

    //
    // Appends from a given source (i.e. do not filter, but annotate)
    //

    std::set<int> app_locset;
    std::set<int> app_varset; 
    std::set<int> app_refset;
    

    //
    // User-defined meta-info
    //

    bool will_attach_meta;
    bool will_attach_all_meta;
    std::set<std::string> meta_fields;

    //
    // Bialleic SNP status
    //

    bool req_biallelic;
    bool exc_biallelic;
    bool req_monomorphic;
    bool exc_monomorphic;

    //
    // Meta masks
    //

    std::map<std::string,int> meta_eq;
    std::map<std::string,int> meta_ne;

    std::map<std::string,std::set<std::string> > meta_has_text;
    std::map<std::string,std::set<std::string> > meta_has_not_text;

    std::map<std::string,double> meta_gt;
    std::map<std::string,double> meta_ge;

    std::map<std::string,double> meta_lt;
    std::map<std::string,double> meta_le;


    std::map<std::string,int> req_meta_eq;
    std::map<std::string,int> req_meta_ne;

    std::map<std::string,std::set<std::string> > req_meta_has_text;
    std::map<std::string,std::set<std::string> > req_meta_has_not_text;

    std::map<std::string,double> req_meta_gt;
    std::map<std::string,double> req_meta_ge;

    std::map<std::string,double> req_meta_lt;
    std::map<std::string,double> req_meta_le;

    bool meta_mask;


    //
    // Genotype masks
    //

    std::map<std::string,int> geno_eq;
    std::map<std::string,int> geno_ne;

    std::map<std::string,std::string> geno_has_text;
    std::map<std::string,std::string> geno_has_not_text;

    std::map<std::string,double> geno_gt;
    std::map<std::string,double> geno_ge;

    std::map<std::string,double> geno_lt;
    std::map<std::string,double> geno_le;


    std::map<std::string,int> req_geno_eq;
    std::map<std::string,int> req_geno_ne;

    std::map<std::string,std::string> req_geno_has_text;
    std::map<std::string,std::string> req_geno_has_not_text;

    std::map<std::string,double> req_geno_gt;
    std::map<std::string,double> req_geno_ge;

    std::map<std::string,double> req_geno_lt;
    std::map<std::string,double> req_geno_le;

    bool geno_mask;
    
    //
    // VarDB functions
    //

    bool exact_vmerge;
    downcode_mode_t downcode_mode;
    bool assume_missing_is_ref;

    //
    // EM caller
    //

    bool use_em;
    bool em_replace;
    double em_threshold;
    
    //
    // Frequency filters
    //

    bool mac_filter;
    int mac_lower;
    int mac_upper;

    bool maf_filter;
    double maf_lower;    
    double maf_upper;

    //
    // HWE filters
    //

    bool use_hwe_filter;
    double hwe_lower;
    double hwe_upper;

    //
    // QUAL filters
    //

    bool use_qual_filter;
    dbl_range qual;
    

    //
    // Null filters, case/control filrers
    //

    bool has_null_filter;
    int_range null_fltr;

    bool has_case_control_filter;
    int_range case_fltr;
    int_range control_fltr;

    //
    // Annotation filters and appends
    //

    bool annot;
    bool annot_append;
    std::vector<std::string> in_annotations;
    std::vector<std::string> req_annotations;
    std::vector<std::string> ex_annotations;

    bool f_include_annotation( const Variant & );
    bool f_require_annotation( const Variant & );
    bool f_exclude_annotation( const Variant & );


    //
    // Meta-information expression filters
    //
    
    Eval eval_expr;
    bool eval_expr_set;
    void set_filter_expression(const std::string & );
    bool eval_filter_includes;

    Eval var_eval_expr;
    bool var_eval_expr_set;
    void var_set_filter_expression(const std::string & );
    bool var_eval_filter_includes;
    
 public:
    bool filter_expression() const { return eval_expr_set; } 
    bool filter_expression_requires_genotypes() const;
    bool calc_filter_expression( SampleVariant & );
    bool calc_filter_expression( SampleVariant & , SampleVariant & );    

    bool var_filter_expression() const { return var_eval_expr_set; } 
    bool var_filter_expression_requires_genotypes() const;
    bool var_calc_filter_expression( Variant & );

 private:
    
    //
    // Display
    //
    
    std::string onoff(const bool b) const { return b ? "on" : "off" ; } 

    friend std::ostream & operator<<( std::ostream & out , const Mask & m )
      {
	if ( m.simple() ) 
	  {
	    out << "simple mask";
	    return out;
	  }
	
	

	
	if ( m.inc_filter.size() > 0 ) out << "  w/ filter-includes\n";
	if ( m.req_filter.size() > 0 ) out << "  w/ filter-requires\n";
	if ( m.exc_filter.size() > 0 ) out << "  w/ filter-excludes\n";
	
	if ( m.geno_eq.size() > 0 || 
	     m.geno_ne.size() > 0 || 
	     m.geno_has_text.size() > 0 || 
	     m.geno_has_not_text.size() > 0 || 
	     m.geno_gt.size() > 0 || 
	     m.geno_ge.size() > 0 || 
	     m.geno_lt.size() > 0 || 
	     m.geno_le.size() > 0 ) out << "  w/ genotype-filters\n";
	
	
	if ( m.in_indset.size() > 0 ) out << "  w/ individual-includes\n";
	if ( m.ex_indset.size() > 0 ) out << "  w/ individual-excludes\n";
	
	
	// Grouping
	
	if ( m.any_grouping() )
	  {
	    out << "grouping\n";
	    out << "  group-by      : ";
	    
	    if ( m.group_locus ) out << "loc\n";
	    else if ( m.group_variant ) out << "var\n";
	    else if ( m.group_locus_set ) out << "loc_set\n";
	    
/* 	    out << "  file-grouping : " << m.onoff( m.file_grouping() ) << "\n" */
/* 		<< "  file-split    : " << m.onoff( m.file_split() ) << "\n"; */
	  }

      
      
      if ( m.loc() ) out << "  w/ loc-includes\n";
      if ( m.var() ) out << "  w/ var-includes\n";
      if ( m.reg() ) out << "  w/ reg-includes\n";
      if ( m.ref() ) out << "  w/ ref-includes\n";
      if ( m.loc_set() ) out << "  w/ loc_set-includes\n";
  
      if ( m.rloc() ) out << "  w/ loc-requires\n";
      if ( m.rvar() ) out << "  w/ var-requires\n";
      if ( m.rreg() ) out << "  w/ reg-requires\n";
      if ( m.rref() ) out << "  w/ ref-requires\n";

      if ( m.xloc() ) out << "  w/ loc-excludes\n";
      if ( m.xvar() ) out << "  w/ var-excludes\n";
      if ( m.xreg() ) out << "  w/ reg-excludes\n";
      if ( m.xref() ) out << "  w/ ref-excludes\n";
      
      if ( m.loc_exceptions() ) out << "  w/ loc-exceptions\n";
      if ( m.var_exceptions() ) out << "  w/ var-exceptions\n";  
      
      
      if ( m.func() ) out << "  w/ functions\n";

      if ( m.filterFunctions.size() > 0 ) out << "  w/ include-functions\n";
      if ( m.req_filterFunctions.size() > 0 ) out << "  w/ reqiure-functions\n";
      if ( m.annot ) out << "  w/ annotation\n";

      if ( m.files() ) out << "  w/ file-includes\n";
      if ( m.xfiles() ) out << "  w/ file-excludes\n";

      if ( m.var_append() ) out << "  w/ var-appends\n";
      if ( m.loc_append() ) out << "  w/ loc-appends\n";
      if ( m.ref_append() ) out << "  w/ ref-appends\n";
      if ( m.loc_set_append() ) out << "  w/ locset-appends\n";

      if ( m.variant_limit() ) out << "  w/ variant limit " << m.variant_limit() << "\n";

      return out;
      
      }



};


#endif
