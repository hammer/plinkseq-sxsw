#ifndef __VIEWS_H__
#define __VIEWS_H__

#include <set>
#include <map>
#include <string>
#include <vector>
#include <sstream>

class Variant;
class VariantGroup;
class Individual;
class GStore;
class int2;

void f_view( Variant & , void * p );

void g_view( VariantGroup & , void * p );

void i_view( Individual & );


// Specialist views

void f_view_lik( Variant & , void * p );

void f_view_matrix( Variant & , void * p );

void f_view_meta_matrix( Variant & , void * p );

void f_view_var_meta_matrix( Variant & , void * p );

void f_view_gene_matrix( VariantGroup & , void * p );

void f_view_tped( Variant & , void * p );

void f_extra_qc_metrics( Variant & , void * p);

void f_simple_counts( Variant & , void * p );

struct XQCstats {
  
  XQCstats() 
  {
    curr_chr = -1;
    em_stats = false;
  }

  int curr_chr;
  bool em_stats; // do extra EM-based stats (entropy, etc)
  
  template<class T> void add( const int2 & pos , const T & t , int k )
  {
    std::stringstream ss;
    ss << t;
    data[pos].push_back( ss.str() );
    datak[pos].push_back( k );
  }

  std::map<int2,std::vector<std::string> > data;
  std::map<int2,std::vector<int> > datak;
  std::map<int2,std::set<int> > neighbours;
  
  std::map<int2,std::set<int> >::iterator flush( std::map<int2,std::set<int> >::iterator n );    
  void flush();

};

struct OptSimpleCounts
{
  bool apply_annot;
  bool dichot_pheno;
  bool show_filter;
  std::set<std::string> meta;
};


struct OptVView {

  OptVView() 
  {
    vmeta = false;
    geno = false;
    gmeta = false;
    vexpand = false;
    show_samples = false;
  }
  
  bool vmeta;
  bool geno;
  bool gmeta;
  bool vexpand;
  bool show_samples;
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



struct OptGMatrix {

  OptGMatrix(GStore * g ) : g(g)  
  {
    hide_zero_variance = false;
  }
  GStore * g;
  bool hide_zero_variance;
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
