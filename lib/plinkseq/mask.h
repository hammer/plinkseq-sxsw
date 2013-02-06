#ifndef __MASK_H__
#define __MASK_H__

#include <set>

#include "regions.h"
#include "plinkseq/variant.h"
#include "plinkseq/eval.h"
#include "plinkseq/defs.h"

#include "locdb.h"
#include "refdb.h"

class VarDBase;
class LocDBase;
class RefDBase;

struct mask_group_t 
{
  mask_group_t( const std::string & g , const std::string & d )
  : name(g) , desc(d) { } 
  std::string name;
  std::string desc;
  bool operator<( const mask_group_t & rhs ) const { return name < rhs.name; } 
};

  
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
  
    bool operator<( const mask_command_t & rhs ) const ;

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
    
  int include_seg( int x );    
  int require_seg( int x );    
  int exclude_seg( int x );    

  int include_var( int x );   
  int require_var( int x );   
  int exclude_var( int x );   
  
  int include_varset( int x );   
  int require_varset( int x );   
  int exclude_varset( int x );   

  int include_ref( int x );   
  int require_ref( int x );   
  int exclude_ref( int x );   

  int include_loc( const std::string & n );
  int require_loc( const std::string & n );
  int exclude_loc( const std::string & n );
  
  int include_seg( const std::string & n );
  int require_seg( const std::string & n );
  int exclude_seg( const std::string & n );

  int include_var( const std::string & n );
  int require_var( const std::string & n );
  int exclude_var( const std::string & n );
  
  int include_varset( const std::string & n );
  int require_varset( const std::string & n );
  int exclude_varset( const std::string & n );

  int include_ref( const std::string & n );
  int require_ref( const std::string & n );
  int exclude_ref( const std::string & n );
  
  void individual_data(const bool b) { inddata = b; }
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
  const std::set<int> & included_varset() { return in_varset_set; }
  const std::set<Region> & included_reg() const { return in_regions; }
  const std::set<Region> & included_ereg() const { return in_eregions; }
  const std::set<std::string> & included_id() const { return in_ids; }
  const std::set<int> & included_ref() const { return in_refset; }


  const std::set<int> & required_loc() const { return req_locset; }
  const std::set<int> & required_locset() const { return req_locset_set; }
  const std::set<int> & required_var() const { return req_varset; }
  const std::set<int> & required_varset() const { return req_varset_set; }
  const std::set<Region> & required_reg() const { return req_regions; }
  const std::set<Region> & required_ereg() const { return req_eregions; }
  const std::set<std::string> & required_id() const { return req_ids; }
  const std::set<int> & required_ref() const { return req_refset; }

  const std::set<int> & excluded_loc() const { return ex_locset; }
  const std::set<int> & excluded_locset() const { return ex_locset_set; }
  const std::set<int> & excluded_var() const { return ex_varset; }
  const std::set<int> & excluded_varset() const { return ex_varset_set; }
  const std::set<Region> & excluded_reg() const { return ex_regions; }
  const std::set<Region> & excluded_ereg() const { return ex_eregions; }
  const std::set<std::string> & excluded_id() const { return ex_ids; }
  const std::set<int> & excluded_ref() const { return ex_refset; }
  
  const std::set<int> appended_ref() const { return app_refset; }
  const std::set<int> appended_var() const { return app_varset; }
  const std::set<int> appended_var_set() const { return app_varset_set; }
  const std::set<int> appended_loc() const { return app_locset; }
  const std::set<int> appended_loc_set() const { return app_locset_set; }
  
  
  //
  // Allele-specific variant-set masks
  //
  
  bool apply_vset_allelemap() { return using_allelemap; } 

  bool apply_vset_allelemap_excludes() { return using_allelemap_excludes; } 
  
  bool attach_vset_allelemap();

  bool attach_vset_allelemap_excludes();
  
  std::vector<std::string> fetch_vset_allelemap( const uint64_t & idx , bool * okay )
    {
      std::map<uint64_t,std::vector<std::string> >::iterator ii = allelemap.find( idx );
      if ( ii == allelemap.end() ) { *okay = false; std::vector<std::string> dummy; return dummy; }
      *okay = true; return ii->second;
    }
  
  std::vector<std::string> fetch_vset_allelemap_excludes( const uint64_t & idx , bool * okay )
    {
      std::map<uint64_t,std::vector<std::string> >::iterator ii = allelemap_excludes.find( idx );
      if ( ii == allelemap_excludes.end() ) { *okay = false; std::vector<std::string> dummy; return dummy; }
      *okay = true; return ii->second;
    }
  
  bool vset_allelemap_any_excludes( const uint64_t & idx )
  {
    return allelemap_excludes.find( idx ) != allelemap_excludes.end() ;
  }

  //
  // Sets within locdb, or segdb
  //

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
  // and so alternate name is not assumed to be unique; return number of transcripts added
  int subset_loc_altname(const std::string & , const std::string&);
  int subset_loc_altname(const std::string & , const std::vector<std::string>&);
    

  //
  // User-defined regions
  //

  void include_reg( const Region & r ) {  in_regions.insert(r); }
  void include_reg( const std::vector<std::string> & r );
  void require_reg( const Region & r ) {  req_regions.insert(r); }
  void exclude_reg( const Region & r ) {  ex_regions.insert(r); }

  //
  // Simular, but exact match (not overlap) regions (eregions)
  //

  void include_ereg( const Region & r ) {  in_eregions.insert(r); }
  void include_ereg( const std::vector<std::string> & r );
  void require_ereg( const Region & r ) {  req_eregions.insert(r); }
  void exclude_ereg( const Region & r ) {  ex_eregions.insert(r); }


  //
  // Variant IDs
  //

  void include_id( const std::vector<std::string> & r );
  void require_id( const std::vector<std::string> & r );
  void exclude_id( const std::vector<std::string> & r );

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

  void func( const mask_func_t f ) { filterFunctions.insert(f); } 
  void include_func( const mask_func_t f ) { filterFunctions.insert(f); }
  void require_func( const mask_func_t f ) { req_filterFunctions.insert(f); } 

  bool eval(Variant & v, void * p = NULL);
  
  //
  // Restriction of maximum number of variants that can be returned (safety-net for R)
  //

  void limit(const int i) { max_var_count = i; } 
  
  //
  // How to handle multiple sites per genomic position
  //
  
  bool fail_on_sample_variant() const ;

  bool test_fail_on_sample_variant(int n , int m ) const ;
  
  void fail_on_sample_variant( int n , int m );

  //
  // Specify what meta-information and genotype/meta-genotype information to load
  // (only load what is needed)
  //
  
  bool load_genotype_data() const { return load_genotypes; } 
  void load_genotype_data(const bool b) { load_genotypes = b; } 
  
  bool load_variant_meta() const { return load_vmeta; }
  void load_variant_meta(const bool b) { load_vmeta = b; }

  bool load_variant_meta( const std::string & m ) const;
  void set_load_variant_meta( const std::string & m );
  
  bool load_genotype_meta() const { return load_gmeta; }
  void load_genotype_meta(const bool b) { load_gmeta = b; }

  bool load_genotype_meta( const std::string & m ) const;
  void set_load_genotype_meta( const std::string & m );

  //
  // Directly attach meta-information from a file
  //

  bool onthefly_append() const { return otf_meta_append; }  
  void onthefly_attach_from_file( const std::string & filename );
  void onthefly_attach_to_variant( Variant & parent );

  

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

/*   void include_indiv_segment( const uint64_t & ); */
/*   void require_indiv_segment( const uint64_t & ); */
/*   void exclude_indiv_segment( const uint64_t & ); */

/*   void include_indiv_segment( const std::string & ); */
/*   void require_indiv_segment( const std::string & ); */
/*   void exclude_indiv_segment( const std::string & ); */
  

  //
  // Appends
  //
    
  int append_loc( int x );
  int append_loc( const std::string & n );
  
  int append_loc_set( int x );
  int append_loc_set( const std::string & n , const std::string & m );
  
  int append_var( int x );      
  int append_var( const std::string & n ); 

  int append_var_set( int x );      
  int append_var_set( const std::string & n ); 

  int append_ref( int x );	
  int append_ref( const std::string & n ); 

    
  //
  // Groups (allowed a single group, either variant-set or locus set)
  //
  
  void group_var(const int g);

  void group_var(const std::string & g);

  void group_var_set(const int g);

  void group_var_set(const std::string & g);

  void group_loc(const int g);

  void group_loc(const std::string & g);    
  
  void group_reg(const std::vector<std::string> & g);

  void group_loc_set(const int g);

  void group_loc_set(const std::string & g, const std::string & h);    

  void process_empty_groups(bool b) { empty_groups = b; }

  bool process_empty_groups() const { return empty_groups; }

  int group_set() const 
    { 
      if ( group_var_set() ) return group_variant_set;
      if ( group_loc_set() ) return group_locus_set;
      if ( group_loc() ) return group_locus;
      return 0;
    }  
  
  bool any_grouping() const 
  { 
    return all_group || named_grouping(); //  || file_grouping() || file_split() ; 
  } 
  
  bool all_grouping() const 
  {
    return all_group;
  }

  bool named_grouping() const 
  { 
    return group_var() || group_loc() || group_var_set() || group_loc_set() ; 
  } 

  bool group_var() const { return group_variant != 0; }

  bool group_var_set() const { return group_variant_set != 0; }

  bool group_loc() const { return group_locus != 0; }

  bool group_loc_set() const { return group_locus_set != 0; }
  
  bool group_reg() const { return group_region; }

  void ensure_single_include_group() 
    {
      
      if ( ! named_grouping() ) return;
      
      if ( group_loc() ) 
	{
	  if ( in_varset.size() || in_eregions.size() || in_regions.size() || in_locset_set.size() || in_varset_set.size() ) 
	    Helper::halt( "you cannot specify other includes in the mask with loc.group" );

	  in_locset.clear();
	  in_locset.insert( group_set() );
	}
      else if ( group_var() ) 
	{
	  if ( in_locset.size() || in_eregions.size() || in_regions.size() || in_locset_set.size() || in_varset_set.size() ) 
	    Helper::halt( "you cannot specify other includes in the mask with var.group" );

	    in_varset.clear();
	    in_varset.insert( group_set() );
	}
      else if ( group_reg() )
	{
	  if ( in_locset.size() || in_eregions.size() || in_varset.size() || in_locset_set.size() || in_varset_set.size() ) 
	    Helper::halt( "you cannot specify other includes in the mask with reg.group" );
	}
      else if ( group_loc_set() )
	{
	  if ( in_locset.size() || in_varset.size() || in_eregions.size() || in_regions.size() || in_varset_set.size() ) 
	    Helper::halt( "you cannot specify other includes in the mask with locset.group" );

	  in_locset_set.clear();
	  in_locset_set.insert( group_set() );
	}
      else if ( group_var_set() )
	{

	  if ( in_locset.size() || in_varset.size() || in_eregions.size() || in_regions.size() || in_locset_set.size() ) 
	    Helper::halt( "you cannot specify other includes in the mask with varset.group" );

	  in_varset_set.clear();
	  in_varset_set.insert( group_set() );
	}

    }
  
  //
  // Simple file-based inclusions, exclusions
  //
  
  int include_file( const std::string & filetag );
  int exclude_file( const std::string & filetag );

  bool use_file( const int ) const;
  
  // Keep track of site-only VCFs too
  
  bool site_only(const int i) const { return site_only_samples.find(i) != site_only_samples.end() ; } 
  void set_site_only( const int i ) { site_only_samples.insert(i); }
  
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

  void only_indels( const bool b ) { req_indels = b; }
  bool only_indels() const { return req_indels; } 

  void skip_indels( const bool b ) { exc_indels = b; } 
  bool skip_indels() const { return exc_indels; } 
  
  void only_mnps( const bool b ) { req_mnps = b; }
  bool only_mnps() const { return req_mnps; } 

  void skip_mnps( const bool b ) { exc_mnps = b; } 
  bool skip_mnps() const { return exc_mnps; } 

  void only_snps( const bool b ) { req_snps = b; }
  bool only_snps() const { return req_snps; } 

  void skip_snps( const bool b ) { exc_snps = b; } 
  bool skip_snps() const { return exc_snps; } 

  void only_novel( const bool b ) { req_novel = b; }
  bool only_novel() const { return req_novel; } 

  void skip_novel( const bool b ) { exc_novel = b; } 
  bool skip_novel() const { return exc_novel; } 
  
  void add_allele( const std::vector<std::string> & a );
  void add_allele_ex( const std::vector<std::string> & a );

  bool allele_filter() const { return do_allele_match; }
  bool allele_ex_filter() const { return do_allele_ex_match; }

  bool allele_match( const std::string & , const std::string & );

  bool test_allele( const std::string & , const std::string & );
  bool test_allele_ex( const std::string & , const std::string & );
  
  // 
  // QUAL filter
  //

  void qual_filter( const std::string & s ) { qual.set(s); use_qual_filter = true; }

  //
  // Variant meta-information masks
  //

  int meta_set( const std::string & key );
  int meta_not_set( const std::string & key );
  
  int meta_equals( const std::string & key , int value );
  int meta_not_equals( const std::string & key , int value );

  int meta_equals( const std::string & key , const std::string & value );
  int meta_not_equals( const std::string & key , const std::string & value );
  
  int meta_greater( const std::string & key , double value );
  int meta_greater_equal( const std::string & key , double value );

  int meta_less( const std::string & key , double value );
  int meta_less_equal( const std::string & key , double value );
  
  // as above, for requires

  int req_meta_set( const std::string & key );
  int req_meta_not_set( const std::string & key );

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
      return meta_is_set.size() > 0 ||
	meta_is_not_set.size() > 0 ||
	meta_eq.size() > 0 || 
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
      return req_meta_is_set.size() > 0 || 
	req_meta_is_not_set.size() > 0 ||
	req_meta_eq.size() > 0 || 
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
  // Borders on locdb
  //

  int border_3prime() const { return loc_border_3p; }
  void border_3prime(const int a) { loc_border_3p = a; }
  
  int border_5prime() const { return loc_border_5p; }
  void border_5prime(const int a ) { loc_border_5p = a; }

  
  //
  // Genotype masks
  //

  int geno_set( const std::string & key );
  int geno_not_set( const std::string & key );

  int geno_equals( const std::string & key , int value );
  int geno_not_equals( const std::string & key , int value );

  int geno_equals( const std::string & key , const std::string & value );
  int geno_not_equals( const std::string & key , const std::string & value );
  
  int geno_greater( const std::string & key , double value );
  int geno_greater_equal( const std::string & key , double value );

  int geno_less( const std::string & key , double value );
  int geno_less_equal( const std::string & key , double value );
  
  // as above, for requires

  int req_geno_set( const std::string & key );
  int req_geno_not_set( const std::string & key );

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
      return geno_is_set.size() > 0 ||
	geno_is_not_set.size() > 0 || 
	geno_eq.size() > 0 || 
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
      return req_geno_is_set.size() > 0 || 
	req_geno_is_not_set.size() > 0 || 
	req_geno_eq.size() > 0 || 
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
      geno_segmask = seg() || rseg() || xseg(); 
      if ( geno_segmask ) prep_segmask();
    }
  
  bool genotype_mask() const { return geno_mask; }

  bool eval( const Genotype & g ) const;

  
  bool genotype_segmask() const { return geno_segmask; } 

  bool eval_segmask( const int i , const Region & );

  bool in_any_segmask( const Region & var , const std::vector<Region> & segs );
  
  bool in_all_segmask( const Region & var , const std::map<int,std::vector<Region> > & segs );

  void prep_segmask();
  

  //
  // Auto-conversion of Null to Reference
  //

  void assuming_null_is_reference( const bool b )
  { assume_missing_is_ref = b; } 
  
  bool assuming_null_is_reference() const
  { return assume_missing_is_ref; } 
  

  //
  // Treatment of REF '0' allele when merging variants -- should '0' indicate only the absence of ALT?
  // i.e. and so might not be inconsistent with another SampleVariant containing a 'contradictory' ALT

  //  e.g. if person really C/G genotype, but represented as follows: 

  //          REF   ALT   P0001
  //            A     C     0/1
  //            A     G     1/1
  
  // OR, more likely, a combination of an indel and SNP
  
  //   true genotype C/AA
  //          REF    ALT   P0001
  //            A    C     0/1
  //            A    AA    0/1
  
  //  This is handled in make_consensus() function in variant.cpp
  //  Make this phase/phase-set aware too

  void soft_reference_calls( const bool b ) { soft_ref = b; }
  bool soft_reference_calls() const { return soft_ref; } 

  
  //
  // X chromosome encoding
  //

  void ploidy( const std::string & chr , ploidy_t t );
  ploidy_t ploidy( const std::string & chr ) const;

  void set_pseudo_autosomal( const Region & region );
  bool pseudo_autosomal( const Variant & var ) const;
  
  void set_genotype_scoring_model( const std::string & );

  //
  // Variant merging and allele-downcoding
  //
  
  merge_mode_t mergemode() const { return merge_mode; } 
  void mergemode( const merge_mode_t m ) { merge_mode = m; }
 
  downcode_mode_t downcode() const { return downcode_mode; } 
  void downcode( const downcode_mode_t d ) { downcode_mode = d; } 

  

  //
  // Invoke EM caller
  //
  
  bool EM_caller() const { return use_em; }
  void EM_caller( const bool b ) { use_em = b; }
  
  double EM_threshold() const { return em_threshold; } 
  void EM_threshold( const double t ) { em_threshold = t; }
  
  
  //
  // Allele-frequency based variant filters
  //

  bool count_filter() const { return mac_filter; }
  bool frequency_filter() const { return maf_filter; }

  bool alt_count_filter() const { return aac_filter; }
  bool alt_frequency_filter() const { return aaf_filter; }

  bool count_filter(const int i) const 
  {       
    // if mac threshold is -1, means ignore
    if ( mac_lower >= 0 && i < mac_lower ) return false;
    if ( mac_upper >= 0 && i > mac_upper ) return false;
    return true;
  }
  
  bool alt_count_filter(const int i) const 
  {       
    // if mac threshold is -1, means ignore
    if ( aac_lower >= 0 && i < aac_lower ) return false;
    if ( aac_upper >= 0 && i > aac_upper ) return false;
    return true;
  }

  bool frequency_filter(const double f) const 
    { 
      // if maf threshold is -1, means ignore
      if ( maf_lower >= 0 && f < maf_lower ) return false;
      if ( maf_upper >= 0 && f > maf_upper ) return false;
      return true;      
    }

  bool alt_frequency_filter(const double f) const 
    { 
      // if maf threshold is -1, means ignore
      if ( aaf_lower >= 0 && f < aaf_lower ) return false;
      if ( aaf_upper >= 0 && f > aaf_upper ) return false;
      return true;      
    }

  void minor_allele_count(const int c, const int d ) 
  { 
    mac_filter = true;
    mac_lower = c; 
    mac_upper = d;
  }
  
  void minor_allele_frequency(const double c, const double d ) 
  { 
    maf_filter = true;
    maf_lower = c; 
    maf_upper = d;
  }
  
  void alt_allele_count(const int c, const int d ) 
  { 
    aac_filter = true;
    aac_lower = c; 
    aac_upper = d;
  }
  
  void alt_allele_frequency(const double c, const double d ) 
  { 
    aaf_filter = true;
    aaf_lower = c; 
    aaf_upper = d;
  }

  
  void hwe( double l, double u )
  {
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
  

  bool get_alt_allele_count(int & c, int & d ) 
  {
    if ( ! aac_filter ) return false;
    c = aac_lower; d = aac_upper;
    return true;
  }
  
  bool get_alt_allele_frequency(double & c, double & d ) 
  {
    if ( ! aaf_filter ) return false;
    c = aaf_lower; d = aaf_upper;
    return true;
  }

  //
  // Null genotype filter
  //

  void null_filter( const int_range & r );
  bool null_filter( ) const;
  bool null_filter( const int ) const;

  void null_prop_filter( const dbl_range & r );
  bool null_prop_filter( ) const;
  bool null_prop_filter( const double ) const;


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
  
  void append_annotation() { annot = true; }
  
  //
  // Get a dump of all mask options
  //

  static std::string describe_options();     
  static std::string list_groups( bool verbose = false );     
  static std::string list_masks( const std::string & g );


  //
  // Helper functions
  //
        
  void reset()
    {
      
      inddata = true;

      in_locset.clear();
      in_segset.clear();
      in_varset.clear();
      in_refset.clear();
      
   
      // LOCDB or SEGDB? (for normal loci, and locus-sets)
      loc_or_seg.clear();
      locset_or_segset.clear();
      
      req_locset.clear();
      req_segset.clear();
      req_varset.clear();
      req_refset.clear();

      ex_locset.clear();
      ex_segset.clear();
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

      meta_is_set.clear();
      meta_is_not_set.clear();

      meta_eq.clear();
      meta_ne.clear();
      
      meta_has_text.clear();
      meta_has_not_text.clear();
      
      meta_gt.clear();
      meta_ge.clear();
      
      meta_lt.clear();
      meta_le.clear();


      req_meta_is_set.clear();
      req_meta_is_not_set.clear();

      req_meta_eq.clear();
      req_meta_ne.clear();
      
      req_meta_has_text.clear();
      req_meta_has_not_text.clear();
      
      req_meta_gt.clear();
      req_meta_ge.clear();
      
      req_meta_lt.clear();
      req_meta_le.clear();

      // Genotypes

      geno_is_set.clear();
      geno_is_not_set.clear();

      geno_eq.clear();
      geno_ne.clear();
      
      geno_has_text.clear();
      geno_has_not_text.clear();
      
      geno_gt.clear();
      geno_ge.clear();
      
      geno_lt.clear();
      geno_le.clear();


      req_geno_is_set.clear();
      req_geno_is_not_set.clear();

      req_geno_eq.clear();
      req_geno_ne.clear();
      
      req_geno_has_text.clear();
      req_geno_has_not_text.clear();
      
      req_geno_gt.clear();
      req_geno_ge.clear();
      
      req_geno_lt.clear();
      req_geno_le.clear();

      // variant/allele merging/splitting
      
      // default to keep as separate sites
      merge_mode = MERGE_MODE_NONE;           
      
      // default to treat as separate alleles
      downcode_mode = DOWNCODE_MODE_EACH_ALT;  

      // EM caller      
      use_em = false;
      em_threshold = 0;

      // Files
      
      in_files.clear();
      ex_files.clear();
      
      site_only_samples.clear();

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

      // User-specified exact-match regions

      in_eregions.clear();
      req_eregions.clear();
      ex_eregions.clear();

      // Variant IDs

      in_ids.clear();
      req_ids.clear();
      ex_ids.clear();

      
      // Grouping
	    
      group_region = false;
      all_group = false;
      group_locus = group_variant = group_variant_set = group_locus_set = 0;

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
      
      // Var Super-sets
      in_varset_set.clear();
      req_varset_set.clear();
      ex_varset_set.clear();
      app_varset_set.clear();


      // Allele frequency filters
      // either absolute counts 
      // or frequencies (that respect missing data)
      
      mac_filter = false;
      maf_filter = false;
      
      aac_filter = false;
      aaf_filter = false;

      has_null_filter = false;
      has_null_prop_filter = false;
      has_case_control_filter = false;
	
    }
  

  bool loc() const { return in_locset.size() > 0 ; }
  bool seg() const { return in_segset.size() > 0 ; } 
  bool var() const { return in_varset.size() > 0 ; }
  bool reg() const { return in_regions.size() > 0 ; }  
  bool ereg() const { return in_eregions.size() > 0 ; }  
  bool id() const  { return in_ids.size() > 0; }
  bool ref() const { return in_refset.size() > 0; }

  bool loc_set() const { return in_locset_set.size() > 0; }
  bool var_set() const { return in_varset_set.size() > 0; }

  bool rloc() const { return req_locset.size() > 0; }
  bool rseg() const { return req_segset.size() > 0; }
  bool rvar() const { return req_varset.size() > 0; }
  bool rreg() const { return req_regions.size() > 0; }
  bool rereg() const { return req_eregions.size() > 0; }
  bool rid() const  { return req_ids.size() > 0; }
  bool rref() const { return req_refset.size() > 0; }
  
  bool xloc() const { return ex_locset.size() > 0; }
  bool xseg() const { return ex_segset.size() > 0; }
  bool xvar() const { return ex_varset.size() > 0; }
  bool xreg() const { return ex_regions.size() > 0; }
  bool xereg() const { return ex_eregions.size() > 0; }
  bool xid() const { return ex_ids.size() > 0 ; } 
  bool xref() const { return ex_refset.size() > 0; }
  
  bool requires() const { 
    return req_locset.size() > 0
      || req_varset.size() > 0
      || req_regions.size() > 0
      || req_eregions.size() > 0
      || req_ids.size() > 0;
  }
  
  bool excludes() const { 
    return ex_locset.size() > 0 
      || ex_varset.size() > 0 
      || ex_regions.size() > 0 
      || ex_eregions.size() > 0 
      || ex_ids.size() > 0;  
  }
  
  
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
	|| rreg()
	|| ereg()
	|| rereg()
	|| id()
	|| xreg() 
	|| xereg() 
	|| files() 
	|| xfiles();
    }
  
  
  bool var_append() const { return app_varset.size() > 0; }
  bool var_set_append() const { return app_varset_set.size() > 0; }
  bool loc_append() const { return app_locset.size() > 0; }
  bool loc_set_append() const { return app_locset_set.size() > 0; }
  bool ref_append() const { return app_refset.size() > 0; }

  void ref_allelic( const bool b ) { ref_allelic_match = b; }
  bool ref_allelic() const { return ref_allelic_match; }

  bool append() const 
  { return var_append() || loc_append() || ref_append() || var_set_append() || loc_set_append() ; }
  
  int loc_size() const { return in_locset.size(); }
  int var_size() const { return in_varset.size(); }
  int ref_size() const { return in_refset.size(); } 

  int variant_limit() const { return max_var_count; }
    
  bool loc_any() const { return loc() || loc_append() || rloc() || xloc() || loc_set(); }
  bool var_any() const { return var() || var_append() || rvar() || xvar() || var_set(); }
  bool ref_any() const { return ref_append(); }
  bool reg_any() const { return reg() || rreg() || xreg(); } 
  bool ereg_any() const { return ereg() || rereg() || xereg(); } 

  // 'force' masks not yet implemented
  bool reg_force() const { return using_force_vlist; } 
  void reg_force(const Region & r ) { using_force_vlist = true; force_vlist.insert(r); } 
  bool forced( int,int,int,int,int,int,Region * ) const ;

  bool loc_set_any() const { return loc_set(); }
  bool var_set_any() const { return var_set(); }
  
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

    // TO ADD
/*     std::string id_include_string() cosnt */
/*       { */
	
/*       } */

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
	group_variant = group_variant_set = group_locus = group_locus_set = 0;
	all_group = false;
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
	geno_segmask = false;
	segs.clear();
	req_segs.clear();
	ex_segs.clear();
	assume_missing_is_ref = false;
	using_allelemap = false;
	using_allelemap_excludes = false;
	ref_allelic_match = false;
	soft_ref = false;
	genotype_model = GENOTYPE_MODEL_ALLELIC;
	hard_call = false;
	hard_call_label = "";
	hard_call_threshold = 0.1;
	fixxy_mode = false;
	mac_filter = false;	
	maf_filter = false;
	aac_filter = false;	
	aaf_filter = false;
	has_null_filter = false;
	has_null_prop_filter = false;
	has_case_control_filter = false;
	inc_filter_any = false;
	req_filter_any = false;
	exc_filter_any = false;
	mac_lower = mac_upper = -1;
	maf_lower = maf_upper = -1.0;
	aac_lower = aac_upper = -1;
	aaf_lower = aaf_upper = -1.0;
	use_hwe_filter = false;
	hwe_lower = 0;
	hwe_upper = 1;
	use_qual_filter = false;
	qual.set(".");
	null_fltr.reset();
	null_prop_fltr.reset();
	case_fltr.reset();
	control_fltr.reset();	
	fail_on_sample_variant_allow = -1; // allow inf.
	fail_on_sample_variant_require = -1; // no req.
	use_em = false;
	merge_mode = MERGE_MODE_EXACT;
	downcode_mode = DOWNCODE_MODE_ALL_ALT;
	em_threshold = 0;
	otf_meta_append = false;
	otf_meta.clear();
	will_attach_meta = false;
	will_attach_all_meta = false;
	req_biallelic = false;
	exc_biallelic = false;
	req_indels = false;
	exc_indels = false;
	req_snps = false;
	exc_snps = false;
	req_mnps = false;
	exc_mnps = false;
	do_allele_match = false;
	do_allele_ex_match = false;
	allele_match_ref.clear();
	allele_match_alt.clear();
	allele_exclude_ref.clear();
	allele_exclude_alt.clear();
	req_novel = false;
	exc_novel = false;
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
	loc_border_3p = 0;
	loc_border_5p = 0;
	using_force_vlist = false;
	force_vlist.clear();

	load_genotypes = true;
	load_vmeta = true;
	load_gmeta = true;
	load_partial_vmeta = false;
	load_partial_gmeta = false;
	load_vmeta_list.clear();
	load_gmeta_list.clear();
      } 
    
    void searchDB();
    
    static std::set<mask_command_t> known_commands; // mask commands
    static std::set<mask_group_t> known_groups; // mask groups

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

    std::set<int> in_segset;
    std::set<int> req_segset;
    std::set<int> ex_segset;  

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
    
    // Variant super-sets

    std::set<int> in_varset_set;
    std::set<int> req_varset_set;
    std::set<int> ex_varset_set;
    std::set<int> app_varset_set;

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
    int group_variant_set;
    int group_locus;  // specifies either LOCDB or SEGDB
    int group_locus_set;
    bool empty_groups;
    bool all_group;

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
    
    std::set<int> site_only_samples;

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

    //
    // Regions: exact matches
    //

    std::set<Region> in_eregions;
    std::set<Region> req_eregions;
    std::set<Region> ex_eregions;


    //
    // Force variants to be output
    //
    
    bool using_force_vlist;
    std::set<Region> force_vlist;
    
    
    //
    // Variant IDs
    //
    
    std::set<std::string> in_ids;
    std::set<std::string> ex_ids;
    std::set<std::string> req_ids;


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
    // On the-fly meta-information
    //
    
    bool otf_meta_append;
    std::map<RefVariant,MetaInformation<RefMeta> > otf_meta;
    
    //
    // REFDB matching
    //
    
    bool ref_allelic_match;
    

    //
    // VARSET allele matching
    //

    std::map<uint64_t,std::vector<std::string> > allelemap;
    std::map<uint64_t,std::vector<std::string> > allelemap_excludes;
    bool using_allelemap;    
    bool using_allelemap_excludes;    
    

    //
    // Bialleic SNP status
    //

    bool req_biallelic;
    bool exc_biallelic;

    bool req_monomorphic;
    bool exc_monomorphic;
    
    bool req_indels;
    bool exc_indels;

    bool req_mnps;
    bool exc_mnps;

    bool req_snps;
    bool exc_snps;

    bool req_novel;
    bool exc_novel;

    bool do_allele_match;
    bool do_allele_ex_match;

    std::vector<std::string> allele_match_ref;
    std::vector<std::string> allele_match_alt;

    std::vector<std::string> allele_exclude_ref;
    std::vector<std::string> allele_exclude_alt;

    //
    // Meta masks
    //

    std::set<std::string> meta_is_set;
    std::set<std::string> meta_is_not_set;

    std::map<std::string,int> meta_eq;
    std::map<std::string,int> meta_ne;

    std::map<std::string,std::set<std::string> > meta_has_text;
    std::map<std::string,std::set<std::string> > meta_has_not_text;

    std::map<std::string,double> meta_gt;
    std::map<std::string,double> meta_ge;

    std::map<std::string,double> meta_lt;
    std::map<std::string,double> meta_le;


    std::set<std::string> req_meta_is_set;
    std::set<std::string> req_meta_is_not_set;

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
    // Locus borders
    //
    
    int loc_border_3p;
    int loc_border_5p;


    //
    // Genotype masks
    //

    std::set<std::string> geno_is_set;
    std::set<std::string> geno_is_not_set;

    std::map<std::string,int> geno_eq;
    std::map<std::string,int> geno_ne;

    std::map<std::string,std::string> geno_has_text;
    std::map<std::string,std::string> geno_has_not_text;

    std::map<std::string,double> geno_gt;
    std::map<std::string,double> geno_ge;

    std::map<std::string,double> geno_lt;
    std::map<std::string,double> geno_le;


    std::set<std::string> req_geno_is_set;
    std::set<std::string> req_geno_is_not_set;

    std::map<std::string,int> req_geno_eq;
    std::map<std::string,int> req_geno_ne;

    std::map<std::string,std::string> req_geno_has_text;
    std::map<std::string,std::string> req_geno_has_not_text;

    std::map<std::string,double> req_geno_gt;
    std::map<std::string,double> req_geno_ge;

    std::map<std::string,double> req_geno_lt;
    std::map<std::string,double> req_geno_le;

    // filter on GT meta fields
    bool geno_mask;
    
    // filter GTs on SEGDB (i.e. set to missing potentially)
    bool geno_segmask;


    //
    // cache for storing all individual segments
    // (populated in prep_segmask)
    //

    std::map<int,std::vector<Region> > segs;
    std::map<int,std::map<int,std::vector<Region> > > req_segs;
    std::map<int,std::vector<Region> > ex_segs;


    //
    // VarDB functions
    //

    merge_mode_t      merge_mode;
    downcode_mode_t   downcode_mode;
    bool              assume_missing_is_ref;
    bool              soft_ref;

    genotype_model_t                 genotype_model;
    bool                             fixxy_mode;
    std::map<std::string,ploidy_t>   fixxy_map_chr;
    std::vector<Region>              fixxy_map_par;

    bool                             hard_call;
    std::string                      hard_call_label;
    double                           hard_call_threshold;
    
 public:

    void fixxy( const bool b ) { fixxy_mode = b; }
    bool fixxy() const { return fixxy_mode; }

    bool make_hard_calls() const { return hard_call; } 
    bool make_hard_calls( const std::string & l , double x ) 
    { 
      hard_call = true; 
      hard_call_label = l;
      hard_call_threshold = x ; 
    }
    
    void revise_hard_call( Genotype & g );

 private:

    //
    // EM caller
    //

    bool use_em;
    double em_threshold;
    

    // Frequency filters

    bool         mac_filter;
    int          mac_lower;
    int          mac_upper;

    bool         maf_filter;
    double       maf_lower;    
    double       maf_upper;

    // as above, for alternate, not minor, allele

    bool         aac_filter;
    int          aac_lower;
    int          aac_upper;

    bool         aaf_filter;
    double       aaf_lower;    
    double       aaf_upper;

    // HWE filters

    bool         use_hwe_filter;
    double       hwe_lower;
    double       hwe_upper;


    // QUAL filters

    bool         use_qual_filter;
    dbl_range    qual;
    

    // Null filters, case/control filrers

    bool         has_null_filter;
    int_range    null_fltr;

    bool         has_null_prop_filter;
    dbl_range    null_prop_fltr;

    bool         has_case_control_filter;
    int_range    case_fltr;
    int_range    control_fltr;


    // Annotation filters and appends

    bool                       annot;
    bool                       annot_append;

    std::vector<std::string>   in_annotations;
    std::vector<std::string>   req_annotations;
    std::vector<std::string>   ex_annotations;
 

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
    

    //
    // Load indicators (do we need to get from disk?)
    //
    
    bool load_genotypes;
    bool load_vmeta;
    bool load_gmeta;
    bool load_partial_vmeta;
    bool load_partial_gmeta;
    std::set<std::string> load_vmeta_list;
    std::set<std::string> load_gmeta_list;

    
 public:
    bool filter_expression() const { return eval_expr_set; } 
    bool filter_expression_requires_genotypes() const;
    bool calc_filter_expression( Variant & , SampleVariant & );
    bool calc_filter_expression( Variant & , SampleVariant & , SampleVariant & );    

    bool var_filter_expression() const { return var_eval_expr_set; } 
    bool var_filter_expression_requires_genotypes() const;
    bool var_calc_filter_expression( Variant & );


};


#endif
