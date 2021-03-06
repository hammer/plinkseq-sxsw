
#ifndef __GSTORE_H__
#define __GSTORE_H__


#include "plinkseq/variant.h"
#include "plinkseq/genotype.h"
#include "plinkseq/individual.h"
#include "regions.h"
#include "permute.h"
#include "phmap.h"

// Main SQL data-stores

#include "vardb.h"
#include "locdb.h"
#include "refdb.h"
#include "seqdb.h"
#include "protdb.h"

// In memory stores / interfaces

#include "filemap.h"
#include "inddb.h"
#include "vcf.h"

// Helper functions, etc

#include "plinkseq/defs.h"

#include <vector>
#include <map>


class Genotype;
class Variant;
class File;
class Variant;

extern Log plog;
extern Log debug;

class GStore { 
    
 public:
  
  //
  // Core databases
  //
  
  // Locus (reference) database
  
  LocDBase         locdb;

  // Protein database
  ProtDBase        protdb;

  // Segment database (project-specific)

  LocDBase         segdb;
  
  
  // Variant/genotype database
  
  VarDBase         vardb;


  // Reference variant database
  
  RefDBase         refdb;

  
  // Sequence database
  
  SeqDBase         seqdb;

  
  // Individual database
  
  IndDBase         inddb;
  
  
  //
  // Helper classes
  //
  

  // Global Individual Map

  IndividualMap    indmap;
  

  // Phenotype Map

  PhenotypeMap    phmap;
  
  
  // Permutation class
  
  Permute          perm;

  
  // Filemap
  
  FileMap          fIndex;

  // Pointer to single-file BCF iteration mode
  
  BCF *            bcf;


  // Password
  
  void set_pwd( const std::string & p ) { _pwd = p; } 
  bool pwd( const std::string & p) const { return _pwd == "" || _pwd == p; }
  

  // Other misc. params set/get functions go here

  void show_id( const bool b ) { param_var_show_ids = b; } 
  bool show_id() const { return param_var_show_ids; }

  void set_param( const std::string & );


  /////////////////////////////////////////////////
  // Initialise a genotype store

  GStore(bool r = false);
  
  
  /////////////////////////////////////////////////
  // Delete a genotype store
  
  ~GStore() 
    { 
      plog.print_warnings();
      cleanUp();

    };
  

  void cleanUp() 
  {

      fIndex.reset();
      
      if ( bcf ) 
	{
	  delete bcf;
	  bcf = NULL;
	}
      
      

  }
  
  
//////////////////////////////////////////////////////
//                                                  //
//  C/C++ API for GStore class                      //
//                                                  //
//////////////////////////////////////////////////////

  //
  // Primary library project specificiations
  //

  bool set_project( const std::string & filename, bool verbose = false);

  void reset();
  
  bool register_mask( Mask & );

  void buildIndex();
  
  std::string summary( bool );

  void show_version() const;

  std::map<std::string,std::string> version() const;
  
  bool single_file_mode() const { return in_single_file_mode; }
  void single_file_mode( const bool b ) { in_single_file_mode = b; }

  bool single_file_bcf() const { return single_file_bcf_mode; }
  void single_file_bcf( const bool b ) { single_file_bcf_mode = b; }

  bool has_project_file() const { return has_projfile; }
  void has_project_file( const bool b ) { has_projfile = b; }

  
  //
  // Some helper functions (of SEGDB, LOCDB or neither)
  //
  
  LocDBase * resolve_locgroup( const std::string & );
  
  
  // Note most functions below are now redundant and can be removed.
  
  //
  // SNP/main variant database functions
  //
  
  void vardb_new(std::string filename) { vardb.newDB(filename); }
  void vardb_attach(std::string filename) { vardb.attach(filename); }
  void vardb_dettach() { vardb.dettach(); }
  bool vardb_load_vcf( Mask & mask , 
		       const std::set<std::string> & , 
		       const std::set<std::string> & , 
		       const std::string * region_mask = NULL );
  bool vardb_load_vcf( const std::string & , const std::string & , const std::string & , 
		       Mask & mask , 
		       const std::set<std::string> & , 
		       const std::set<std::string> & , 
		       const std::set<Region> * = NULL );
  void vardb_write_vcf(std::string filename);
  void vardb_load_plink(std::string filename);
  void vardb_write_plink(std::string filename);
  void vardb_load_variant_information(std::string filename);
  void vardb_iterate( void (*)(Variant&, void *) , void * data = NULL );
  void vardb_giterate( void (*)(std::vector<Variant>&, void *) , int grp_id , void * data = NULL );

      

  // Return a variant given it's index
  Variant vardb_variant(uint64_t);
  
  //
  // Individuals
  //

  void inddb_count();
  void inddb_display();
  
  //
  // Files
  //

  void filedb_display();


  //
  // Data use / application of functions (tests)
  //
    
  void vardb_allele_freq();

  // Mask out or in varints, genotypes, individuals, etc, for iteration
  
  void vardb_mask();
  void inddb_mask();
  

  //
  // Locus-database/segment related functions
  //
  
  void locdb_load_names( const std::string & file , const std::string & name );
  void locdb_summary();
  void locdb_remove_group( ID_t id );
  void locdb_subregions(bool t);
  void locdb_meta(bool t);
  void locdb_rename( const std::string & group , const std::string & alias, const std::string & new_label );
  void locdb_extract_intersection( const std::string & group1, const std::string & group2 , const std::string & newLabel);
  void locdb_display_regions( const std::string & name );
 
  void locdb_overlap_analysis( const std::string & target, 
			       const std::string & preload,
			       const std::string & alias , 
			       const std::string & listmode );
  
  // Protein database functions
  void protdb_attach(std::string);
  void protdb_summary();
  void protdb_init();
  ProtFeatureSet protdb_lookup(Variant &);

  //
  // Reference database functions
  //
  
  void refdb_new(std::string);
  void refdb_attach(std::string);
  void refdb_summary();
  RefVariant refdb_lookup(Variant &, int);


  //
  // Meta-information
  //

  void summary_metatypes() { Helper::metatype_summary(); }


  //
  // Misc. parameters
  //
  
  void gseq_tracking( const std::string &  h , const std::string & n ) { gseq=true; ghist=h; gjobn=n; }
  
  bool gseq_mode() const { return gseq; }
  std::string gseq_history() const { return ghist; }
  std::string gseq_job() const { return gjobn; }  
  
  void R_mode( const bool b ) { r_mode = b; } 

 private:
  
  bool r_mode;
  
  bool gseq;
  
  std::string ghist;

  std::string gjobn;

  bool in_single_file_mode;

  bool single_file_bcf_mode;
  
  BCF * singe_bcf;

  bool has_projfile;

  std::string _pwd;

  bool param_var_show_ids;

};


#endif
