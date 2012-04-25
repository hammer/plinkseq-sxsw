#ifndef __VARDB_H__
#define __VARDB_H__

#include "iterrep.h"
#include "sqlwrap.h"
#include "variant.h"
#include "regions.h"
#include "mask.h"
#include "vgroup.h"

#include <string>
#include <map>
#include <set>
#include <vector>


class Mask;
class IndividualMap;
class BCF;

class VarDBase {


 public:
  

 VarDBase(IndividualMap & indmap) : indmap(indmap)
  {

    using_compression = false;
    
    // only combine SampleVariants into a single Variant if they 
    // specify *exactly* the same underlying variation (i.e. not 
    // just overlapping)
    
    merge_mode = MERGE_MODE_EXACT;

    
    //-------------------------------------------------------------
    // mutli-allelic downcoding 
    
    //      REF=C  ALT=G,T

    //  given the above, we can either merge all alt-alleles

    //  --> REF=C  ALT=G|T  

    // or downcode to K-1 biallelic variants

    //  --> REF=C   ALT=G
    //      REF=C   ALT=T
    
    // when doing the latter, we have a choice about how to handle the other 
    // alternate alleles -- either to set as reference, or to make missing
    
    // 0. send multi-allelic variant to user function
    // I. merge all non-ref alleles to make single biallelic 
    // II. send one biallelic per non-ref allele (other non-ref set to missing)
    // III. send one biallelic per non-ref allele (other non-ref set to reference)
    
    downcode_mode = DOWNCODE_MODE_ALL_ALT;

    // track explicitly whether temporary databases are attached
    tmpdb_attached = locdb_attached = refdb_attached = false;
    
  }
  
  ~VarDBase()
    {
      dettach();
    }
  
  bool newDB(std::string name);
  bool attach(std::string name);
  bool attached() { return sql.is_open(); }
  bool check_version();
  bool dettach();
  bool init();
  bool release();
  bool index();
  bool drop_index();

  std::string summary( Mask * mask = NULL , bool ugly = false );
  int2 make_summary(std::string);
  int2 make_summary(int);

  int variant_count();
  int variant_count(uint64_t);
  int set_count(uint64_t);
  int indiv_count(uint64_t);
    
  
  
  //---------------------------------------
  

  // 
  // Insertions
  // 
  
  void insert_header( uint64_t file_id , const std::string & name , std::string value );
  void insert_metatype( uint64_t file_id , const std::string & name , 
			mType mt , int num, int mgrp, const std::string & description );
  uint64_t insert( const std::string & filename , const std::string & tagname );
  uint64_t insert( uint64_t file_id , const Individual & );
  uint64_t insert_consensus( uint64_t file_id , const Variant & );


  //
  // BCF indexing/queries
  //

  void populate_bcf_map();
  void insert_bcf_index( uint64_t file_id , const Variant & , int64_t );
  void store_bcf_n( uint64_t , const std::string & , int t , int s );
  // primary query: get a Variant from a BCF
  bool fetch_bcf( Variant & , const uint64_t & offset );
  

  //
  // Meta-information
  //

  uint64_t insertMetaField( const std::string & , mType, const std::string & );
  void populateMetaFields();


  //
  // Variant groups
  //
  
  uint64_t add_set( const std::string & name , const std::string & desc = "" , bool donotadd = false );
  bool add_var_to_set( const std::string & , Variant & , bool allelic = false );
  
  uint64_t add_superset( const std::string & name , const std::string & desc = "" , bool donotadd = false );
  bool add_set_to_superset( const std::string & , const std::string & );
  
  void drop_set( const std::string & );
  void drop_superset( const std::string & );

  std::vector<std::string> get_sets();
  std::vector<std::string> get_supersets();
  std::vector<std::string> get_sets( const std::string & superset );
  int get_set_size( const std::string & );

  void add_set_description( const std::string & name , const std::string & desc );
  void add_superset_description( const std::string & name , const std::string & desc );
  std::string get_set_description( const std::string & name );
  std::string get_superset_description( const std::string & name );
  
  //
  // Add/remove meta-information on variants (independent of VCF load)
  //
  
  int loader_indep_meta( const std::string & , int , const std::string & );
  int flush_indep_meta( const std::string & );
  void flush_indep_meta( );
  bool attach_indep_metadata( const uint64_t & svar_id , SampleVariant & t , const std::set<std::string> * keys = NULL );


  //
  // Chromosome codes
  //

  bool chr_code( const int , const std::string & , const ploidy_t ); 
  void chr_name( const int , const std::string & );
  int chr_code( const std::string & , ploidy_t * p = NULL );
  std::string chr_name( const int ) ;
  ploidy_t ploidy( const int c );

  bool chr_known( const std::string & n );


  //
  // Queries
  // 
  
  std::string print_headers( uint64_t file_id );	
  void set_mask_metatypes( const Mask & mask );
  void set_metatypes(bool clear = false );
  void set_file_metatypes(uint64_t file_id, bool clear = false );
  
  void append_metainformation( Variant & v , const std::set<int> & grp );

  // File, indiv and meta-information
    
  std::vector<std::string> fetch_individuals(uint64_t file_id);
  std::vector<std::map<std::string,std::string> > fetch_headers(uint64_t file_id);
  std::vector<std::map<std::string,std::string> > fetch_metatypes(uint64_t file_id);
  std::map<int,std::string> fetch_files( Mask * mask = NULL );
  int n_files( Mask * mask = NULL );
  int fileID(const std::string & );

  bool replace_individual_id( const std::string & old_id , const std::string & new_id );
  
  // Fetch single variant at single (1bp) position 
  
  int populate_individual_alignment( IndividualMap &, Mask & );
  int populate_individual_map(PhenotypeMap &, const std::vector<Individual> & inds );
  
  void addMetaFields(Variant &, sqlite3_stmt *, Mask &);
  
  // Redundant:
  int get_individuals( uint64_t file_id );
  
  MetaInformation<VarMeta> metaVariant( uint64_t file_id , uint64_t var_id );
  std::map<uint64_t, MetaInformation<GenMeta> > metaGenotype( uint64_t file_id , uint64_t var_id );
  
  void addGenotypes( uint64_t file_id , uint64_t var_id , Variant & var , bool );
  
  // primary function make a new variant
  SampleVariant & construct( Variant & , sqlite3_stmt * , IndividualMap * align = NULL ); 
  

  //
  // Fetch a single variant, given it's IDX
  //

  Variant fetch(uint64_t idx);

  Variant fetch(int chr, int bp);
  std::set<Variant> key_fetch( const Region & region );

  std::set<Variant> fetch( const Region & );

  //
  // Fetch a single Region/Variant, given it's ID
  //

  Region get_position_from_id( const std::string & id1 , const std::string & id2 ) ;

  //
  // Deletions
  //
  
  void drop(int g);

  void vacuum();


  //
  // Single variant iteration
  //
  
  
  IterationReport iterate( void (*f)(Variant&, void *) , 
			   void * data ,
			   Mask & mask ); 
  
  IterationReport iterate( void (*f)(Variant&, void *) , 
			   void * data = NULL )
    {
      Mask simple;
      return iterate( f, data, simple );
    }
  
  /*!
    Group iteration
    
  */
  
  IterationReport iterate( void (*f)(VariantGroup &, void *) , 
			   void * data ,
			   Mask & mask );
  
  IterationReport generic_iterate( void (*f)(Variant&, void *) ,
				   void (*g)(VariantGroup&, void *) ,
				   void * data ,	
				   Mask & mask );
  
  IterationReport vcf_iterate( void (*f)(Variant&, void *) ,			       
			       void * data ,	
			       Mask & mask );

  bool vcf_iterate_read_header( Mask & mask );

  
  //
  // Other functions
  // 

  void build_temporary_db( Mask & );

  bool eval_and_call( Mask & mask , 
		      IndividualMap * align ,
		      Variant & var , 
		      void (*f)(Variant&, void *) ,
		      void * data );
  
  
  //
  // Temporary in-memory table
  //
  
  void attachMemoryDB();
  void insertMemoryDB(const std::string &);
  void detachMemoryDB();
  

  //
  // Helpers
  //
  
  void begin() { sql.begin(); }
  void commit() { sql.commit(); }
  
  uint64_t lookup_file_id( const std::string & );
  void insert_file_tag( uint64_t , const std::string & );
  std::string file_tag( uint64_t );
  uint64_t file_tag( const std::string & );
  
  bool compression( const bool b ) { using_compression=b; }

  //
  // Misc.
  //


 private:
  
  SQL sql;

  bool tmpdb_attached;
  bool locdb_attached;
  bool refdb_attached;

  bool using_compression; 

  // misc options controlling how variants are fetched/merged
  // these are populated via the Mask

  merge_mode_t    merge_mode;
  downcode_mode_t downcode_mode;


  // SQL queries

  sqlite3_stmt * stmt_vcount;
  sqlite3_stmt * stmt_totvcount;
  sqlite3_stmt * stmt_indcount;
  sqlite3_stmt * stmt_setcount;

  sqlite3_stmt * stmt_insert_header;
  sqlite3_stmt * stmt_fetch_headers;

  sqlite3_stmt * stmt_insert_metatype;
  sqlite3_stmt * stmt_fetch_metatypes;
  sqlite3_stmt * stmt_fetch_files;
  sqlite3_stmt * stmt_fetch_file_id;
  sqlite3_stmt * stmt_fetch_file_summary;

  sqlite3_stmt * stmt_insert_bcf_n;
  sqlite3_stmt * stmt_fetch_bcf;
  sqlite3_stmt * stmt_fetch_bcfs;
  sqlite3_stmt * stmt_insert_bcf_idx;
  

  sqlite3_stmt * stmt_fetch_sets;

  sqlite3_stmt * stmt_insert_chr_code;
  sqlite3_stmt * stmt_insert_chr_name;
  sqlite3_stmt * stmt_fix_chr_code;
  sqlite3_stmt * stmt_fetch_chr_code;
  sqlite3_stmt * stmt_fetch_chr_all_codes;
  sqlite3_stmt * stmt_fetch_chr_all_codes_2;
  sqlite3_stmt * stmt_fetch_chr_name;

  sqlite3_stmt * stmt_insert_file;
  sqlite3_stmt * stmt_insert_file_summary;
  sqlite3_stmt * stmt_insert_individual;
  sqlite3_stmt * stmt_fetch_individual;
  sqlite3_stmt * stmt_fetch_individuals;
  sqlite3_stmt * stmt_replace_individual_id;

  sqlite3_stmt * stmt_fetch_file_from_tag;
  sqlite3_stmt * stmt_fetch_tag_from_file;
  sqlite3_stmt * stmt_insert_file_tag;

  // Core insertion/access statements:

  sqlite3_stmt * stmt_insert_variant_key; 
  sqlite3_stmt * stmt_insert_variant_data; 

  sqlite3_stmt * stmt_fetch_variant_key;
  sqlite3_stmt * stmt_fetch_variant_key_from_id;
  sqlite3_stmt * stmt_fetch_variant_pos;
  sqlite3_stmt * stmt_fetch_variant_range;
	
  sqlite3_stmt * stmt_fetch_variant_data_all;
  sqlite3_stmt * stmt_fetch_variant_data_vmeta_geno;  
  sqlite3_stmt * stmt_fetch_variant_data_vmeta;
  sqlite3_stmt * stmt_fetch_variant_data_geno;

  
  // ID-lookups

  sqlite3_stmt * stmt_fetch_var_from_position;
  sqlite3_stmt * stmt_fetch_var_from_position2;
  sqlite3_stmt * stmt_fetch_var_from_name;
  
  // Aux. meta-information
  
  sqlite3_stmt * stmt_insert_indep_meta_group;
  sqlite3_stmt * stmt_fetch_indep_meta_group;
  sqlite3_stmt * stmt_dump_indep_meta_group;
  sqlite3_stmt * stmt_insert_indep_meta_type;
  sqlite3_stmt * stmt_fetch_indep_meta_type;
  sqlite3_stmt * stmt_insert_indep_meta_value;
  sqlite3_stmt * stmt_fetch_indep_meta_value;
  
  
  // Basic iteration
  
  sqlite3_stmt * stmt_iterate_variants;


  // Sets & supersets

  sqlite3_stmt * stmt_insert_set;
  sqlite3_stmt * stmt_insert_superset;
  sqlite3_stmt * stmt_insert_set_variant;
  sqlite3_stmt * stmt_attach_set_to_superset;
  sqlite3_stmt * stmt_lookup_set;
  sqlite3_stmt * stmt_lookup_superset;
  sqlite3_stmt * stmt_lookup_set_name;
  sqlite3_stmt * stmt_lookup_superset_name;
  sqlite3_stmt * stmt_lookup_set_names;
  sqlite3_stmt * stmt_lookup_set_desc;
  sqlite3_stmt * stmt_lookup_superset_desc;
  sqlite3_stmt * stmt_dump_all_set_names;
  sqlite3_stmt * stmt_dump_all_superset_names;
  sqlite3_stmt * stmt_fetch_set_variants;
  sqlite3_stmt * stmt_fetch_superset_variants;

  // Temporary table

  sqlite3_stmt * stmt_tmp_insert;

  // For any given variant, keep track of what we've added so far
  
  std::map<std::string,int> added;

  // Look-up tables

  std::vector<uint64_t> indiv;


  IndividualMap & indmap;

  void clear()
    {
      indiv.clear();
      varset_map.clear();
      varsuperset_map.clear();     
    }
  

  //
  // Aux. meta-data
  //

  int process_indep_meta_header( const std::string & , const int );
  
  void populate_indep_metadata_map();
  
  std::map<std::string,int> indep_metamap;
  std::map<int,std::string> reverse_indep_metamap;

  std::map<std::string,int2> indep_group_metamap;

  std::map<int,std::string> file_tag_map;
  std::map<std::string,int> reverse_file_tag_map;

  std::map<std::string,int> chr_code_map;
  std::map<int,std::string> chr_name_map;
  std::map<int,ploidy_t> chr_ploidy_map;
  
  std::map<std::string,int> varset_map;
  std::map<std::string,int> varsuperset_map;		       

  std::map<int,BCF*> bcfmap;
  std::map<int,VCFZ*> vcfzmap;
  
  enum fetch_mode_t { ALL = 0 , NO_GMETA , ONLY_VMETA , ONLY_GENO , ONLY_CORE } fetch_mode;
  
};

#endif
