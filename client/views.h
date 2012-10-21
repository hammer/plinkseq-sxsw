#ifndef __VIEWS_H__
#define __VIEWS_H__

#include <set>
#include <map>
#include <string>
#include <vector>
#include <sstream>
#include "plinkseq.h"

class Variant;
class VariantGroup;
class Individual;
class GStore;
class int2;
class ProtDBase;

void f_view( Variant & , void * p );

void g_view( VariantGroup & , void * p );

void i_view( Individual & );


// Specialist views

void g_loc_view( VariantGroup & , void * p );

void g_geneseq( VariantGroup & , void * p );

void f_view_lik( Variant & , void * p );

void f_view_matrix( Variant & , void * p );

void f_view_meta_matrix( Variant & , void * p );

void f_view_var_meta_matrix( Variant & , void * p );

void f_view_gene_matrix( VariantGroup & , void * p );

void f_view_gene_meta_matrix( VariantGroup & , void * p );

void f_view_tped( Variant & , void * p );

void f_extra_qc_metrics( Variant & , void * p);

void f_simple_counts( Variant & , void * p );

struct XQC_varid { 
    XQC_varid( const Variant & v ) 
	: chr(v.chromosome()), bp1(v.position()), alt(v.alternate()) { } 
    int chr;
    int bp1;
    std::string alt;
    bool operator<( const XQC_varid & rhs ) const {
	if ( chr < rhs.chr ) return true;
	if ( chr > rhs.chr ) return false;
	if ( bp1 < rhs.bp1 ) return true;
	if ( bp1 > rhs.bp1 ) return false;
	return  alt < rhs.alt;
    }    
};

struct XQCstats {
    
  XQCstats() 
	{
	  curr_chr = -1;
	  em_stats = false;
	}
  
  int curr_chr;
  bool em_stats; // do extra EM-based stats (entropy, etc)
  
  template<class T> void add( const XQC_varid & pos , const T & t , int k )
	{
	    std::stringstream ss;
	    ss << t;
	    data[pos].push_back( ss.str() );
	    datak[pos].push_back( k );
	}

    std::map<XQC_varid,std::vector<std::string> > data;
    std::map<XQC_varid,std::vector<int> > datak;
    std::map<XQC_varid,std::set<int> > neighbours;
    
    std::map<XQC_varid,std::set<int> >::iterator flush( std::map<XQC_varid,std::set<int> >::iterator n );    
    void flush();
    
};


struct OptSimpleCounts
{
  bool apply_annot;
  bool apply_full_annot;
  bool dichot_pheno;
  bool qt_pheno;
  bool show_filter;
  std::set<std::string> meta;
  bool genotypes;
};

struct OptVView {

  OptVView() : vars( dummy )
  {
    vmeta = false;
    geno = false;
    gmeta = false;
    vexpand = false;
    show_samples = false;
    show_nonmissing_geno = true;
    show_only_minor = false;
    show_only_alt = false;
    mview = false;
    simple = false; 
  }
  
  bool vmeta;
  bool geno;
  bool gmeta;
  bool vexpand;
  bool show_samples;
  bool show_nonmissing_geno;
  bool show_only_minor;
  bool show_only_alt;
  bool mview;
  bool simple;

  // for m-view
  Mask dummy;
  VariantGroup vars;
};


struct OptGView {

  OptGView() 
  {
    vmeta = false;
    vexpand = false;
    geno = false;
    gmeta = false;
    transpose = false;    
    rarelist = false;
    show_phenotype = false;
  }
  
  bool vmeta;
  bool vexpand;
  bool geno;
  bool gmeta;
  bool transpose;
  bool rarelist;
  bool show_phenotype;
};

class Feature;

struct Opt_geneseq {
  Opt_geneseq()   
  { 
    pheno = false; 
    ref = 0; 
    protdb = NULL;
    protdom.clear();
    all_cds = true;
    only_variant_sites = false;
    R_plot = false;
  }

  bool all_cds;  // show *all* CDS codons in output
  bool only_variant_sites;

  bool R_plot;

  bool pheno;
  int ref;
  ProtDBase * protdb;
  std::set<std::string> protdom;
  
};

struct OptGMatrix {

  OptGMatrix(GStore * g ) : g(g)  
  {
    hide_zero_variance = false;
    collapse_01 = false;
  }
  GStore * g;
  bool hide_zero_variance;
  bool collapse_01;
};

struct OptGMetaMatrix {

  OptGMetaMatrix()
  {
    name = "";
    show_mean = false;
    show_sum = false;
    show_flag01 = false;
    show_min = false;
    show_max = false;
    show_na = false;
  }

  bool show_mean;  
  bool show_sum;
  bool show_flag01;
  bool show_min;
  bool show_max;
  bool show_na;
  std::string name;
};

struct OptUniq {
  OptUniq()
  {
    ingroup_req = -1;
    outgroup_allow = 0;
  }
  std::set<Individual*> indiv;
  int ingroup_req;
  int outgroup_allow;
};


#endif
