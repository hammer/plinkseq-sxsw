#ifndef __RINT_H__
#define __RINT_H__

#include "gstore.h"
#include "vgroup.h"

#include <iostream>
#include <vector>
#include <string>

class Rdisplay_options;

//#define R_NO_REMAP 1


extern "C" { 

#include <R.h>
#include <Rversion.h>
#include <Rinternals.h>
#include <R_ext/Rdynload.h>
        
  void R_init_pseqr(DllInfo *info);
  
  // Attach/dettach a database
  
  SEXP Rattach(SEXP d);
  
  SEXP Rdettach();

  //
  // Report an error, print any output
  //

  void R_error( const std::string & );
  void R_warning( const std::string & );
  void R_flush_plog(); 

  //
  // Set a project, databases
  //
  
  SEXP Rset_project(SEXP name);
  SEXP Rsummary();
  
  
  //
  // Variant database
  //

  SEXP Rvardb_new(SEXP filename);
  SEXP Rvardb_attach(SEXP filename);
  SEXP Rvardb_dettach();
    
  
  //
  // Load VCF file directly into R space
  //
  
  SEXP Rdirect_load_vcf( SEXP , SEXP , SEXP );
    
  
  //
  // Make a set of variants, given a set of regions
  //
  
  SEXP Rvardb_make_set(SEXP grp, SEXP name);
  
    
  //
  // Load GTF file
  //
  
  SEXP Rlocdb_attach(SEXP name);
  SEXP Rlocdb_load_gtf( SEXP file, SEXP name);
  SEXP Rlocdb_load_alias( SEXP file , SEXP name );
  
  
  //
  // Convert regions to R list
  //
  
  SEXP Rlocdb_summary();
  SEXP Rfetch_regions(SEXP grp);
  SEXP Rlocdb_collapse_subregions(SEXP x, SEXP y);
  
  SEXP Rfetch_set_names(SEXP x, SEXP y);
  SEXP Rfetch_set_members( SEXP x, SEXP y, SEXP z);


  //
  // Create a Variant(Group) object in R
  //
  
  SEXP Rvariant(Variant & var, Rdisplay_options & );
  SEXP Rsample_variant(SampleVariant & sample, Variant&, Rdisplay_options &);
  SEXP Rvariant_group(VariantGroup & vars, Rdisplay_options &);
  
  
  //
  // Iterate over all variants, applying a function
  //
  
  SEXP Riterate(SEXP fn, SEXP msk, SEXP ret, SEXP rho);
  
  //
  // Get meta-information data-frame 
  //

  SEXP Rgetmeta(SEXP fn, SEXP msk, SEXP rho);
  
  //
  // Helper functions for iteration
  //
  
  void R_iterate_func1( Variant & v ,             void * p );
  
  void R_iterate_func2( VariantGroup & v ,        void * p );
  
  bool R_filter_func( SEXP fn );
  

  //
  // Construct a mask object from a list
    //
  
  Mask R_make_mask(SEXP m);
  
  
  
  // Get list of individuals (per file)
  
  SEXP Rind_list(SEXP f, SEXP p);
  
  // Get headers for a file
  
  SEXP Rhdr_list_int(int f);
  SEXP Rhdr_list(SEXP f);
  
  // Get meta-information for a file
  
  SEXP Rmeta_list_int(int f , VarDBase * vardb );
  SEXP Rmeta_list(SEXP f);
  
  
  // Get list of files
  
  SEXP Rfile_list();

  
  //
  // Reference-database functions
  //
    
  SEXP Rrefdb_load(SEXP x);
  SEXP Rrefdb_attach(SEXP x);
  SEXP Rrefdb_summary();
  SEXP Rrefdb_lookup(SEXP x,SEXP y);
  SEXP Rrefdb_index_lookup(SEXP x);


  //
  // Sequence database functions
  //
  
  SEXP Rseqdb_loadFASTA( SEXP filename );
  SEXP Rseqdb_attach( SEXP filename );
  SEXP Rseqdb_lookup( SEXP pos );
  SEXP Rseqdb_annotate_load( SEXP loc_id );
  SEXP Rseqdb_annotate( SEXP pos, SEXP alleles );
  
  // 
  // Helper functions
  //

  SEXP Rmake_string_vector( std::vector<std::string> & r );

  std::string Rversion();

}


struct Rdisplay_options{ 

  Rdisplay_options() 
    {
      show_consensus = true;
      show_multi_sample = true;
      show_genotypes = true;
      show_genotype_metainformation = true;
    }

  bool show_consensus;
  bool show_multi_sample;
  bool show_genotypes;
  bool show_genotype_metainformation;

};

#endif
