
#ifndef __SUMMARIES_H__
#define __SUMMARIES_H__

#include "func.h"

#include <map>
#include <set>
#include <string>

class GStore;

namespace Pseq 
{

  //
  // 1) Accumulators
  //
  
  //
  // Variant summary statistics that are tracked across the entire
  // (masked) dataset
  //
  
  struct VStat {

  VStat( GStore * g = NULL ) : g(g) 
    {
      // init all variables
      
      nvar = 0;
      ti = tv = 0;
      ti_singleton = tv_singleton = 0;
      n_filtered = n_filtered_singleton = 0;
      dbSNP = dbSNP_singleton = 0;
      thousandG = thousandG_singleton = 0;
      nssnp = nssnp_singleton = 0;
     
      var_depth_label = PLINKSeq::META_DP();
      ind_depth_label = PLINKSeq::META_GENO_DP();
      func_str = PLINKSeq::META_ANNOT();
      ind_qual_label = PLINKSeq::META_GENO_GQ();
      
      dbSNP_label = PLINKSeq::DEFAULT_DBSNP();
      thousandG_label = PLINKSeq::DEFAULT_G1K();
      
      hwe_failure[ 0.01 ] = 0;
      hwe_failure[ 0.0001 ] = 0;
      
      altmin = true;
      single_minor_allele_count = 0;
 }

    void reset()
    {
      nvar = 0;
      ti = tv = 0;
      ti_singleton = tv_singleton = 0;
      n_filtered = n_filtered_singleton = 0;
      dbSNP = dbSNP_singleton = 0;
      thousandG = thousandG_singleton = 0;
      nssnp = nssnp_singleton = 0;
      altmin = true;
      single_minor_allele_count = 0;

      // some maps we can completley clear (keys+values)
      // others, only values (i.e. keys are pre-set)

      call_rate.clear();     
      minor_allele_count.clear();
      total_allele_count.clear();
      minor_allele_freq_bin100.clear();
      filter.clear();
      
      // flagged functional SNPs
      {
	std::map<std::string,int>::iterator i = func.begin();
	while ( i != func.end() )
	  {
	    i->second = 0;
	    ++i;
	  }
      }
    
      // HWE
      {
	std::map<double,int>::iterator i = hwe_failure.begin();
	while ( i != hwe_failure.end() )
	  {
	    i->second = 0;
	    ++i;
	  }
      }
      
      // Variant depth
      {
	std::map<int,int>::iterator i = depth.begin();
	while ( i != depth.end() )
	  {
	    i->second = 0;
	    ++i;
	  }
	
      }

      // Mean per-genotype depth
      {
	std::map<double,int>::iterator i = ind_depth.begin();
	while ( i != ind_depth.end() )
	  {
	    i->second = 0;
	    ++i;
	  }
      }
      
      // Geno GQ
      {
	std::map<double,int>::iterator i = ind_qual.begin();
	while ( i != ind_qual.end() )
	  {
	    i->second = 0;
	    ++i;
	  }	
      }

      
      // per variant quality score
      {
	std::map<double,int>::iterator i = qual.begin();
	while ( i != qual.end() )
	  {
	    i->second = 0;
	    ++i;
	  }		
      }
      
    }


    void set_depth_label( const std::string & s ) { var_depth_label = s; } 
    void add_depth( const int i ) { depth[i] = 0; }   
  
    void add_qual( const double d ) { qual[d] = 0; }   
  
    void set_indiv_depth_label( const std::string & s ) { ind_depth_label = s; } 
    void add_indiv_depth( const double d ) { ind_depth[d] = 0; }   

    void set_indiv_geno_qual_label( const std::string & s ) { ind_qual_label = s; } 
    void add_indiv_geno_qual( const double d ) { ind_qual[d] = 0; }   

    void set_func_label( const std::string & s ) { func_str = s; }
    void add_func( const std::string & s ) { func[s] = 0; }

    void add_filter( const std::string & s ) { filter_out.insert(s) ; }
    
    void set_dbSNP_label( const std::string & s ) { dbSNP_label = s; }
    void set_1KG_label( const std::string & s ) { thousandG_label = s; }

    void add_hwe_p( const double d ) { hwe_failure[d] = 0; }

    //
    // Display function
    //
    
    void report();

    //
    // Data members
    //

    GStore * g ;
    

    // This could apply to one individual or to entire sample
    // The sample will be used to assess things such as MAF, of course
    
    // Simple # of variants, inds

    int nvar;    

    //
    // Call-rate
    //

    // key   = number of genotypes called per variant
    // value = number of sites with this call-rate

    std::map<int,int> call_rate;
    
    
    //
    // Frequency and HWE statistics
    //
    
    bool altmin;
    
    int single_minor_allele_count;
    
    std::map<int,int> minor_allele_count;

    std::map<int,int> total_allele_count;

    std::map<int,int> minor_allele_freq_bin100;
    
    std::map<double,int> hwe_failure;
     
    //
    // Depth at called sites
    //
    
    std::string var_depth_label;
    std::map<int,int> depth;
    
    // mean per-individual depth

    std::string ind_depth_label;
    std::map<double,int> ind_depth;

    // mean GQ

    std::string ind_qual_label;
    std::map<double,int> ind_qual;
    
    // per variant quality score

    std::map<double,int> qual;


    // 
    // Ti/Tv (all SNPs, and at singletons)
    //
    
    int ti, tv;
    int ti_singleton, tv_singleton;

    
    //
    // Filtered-out variants
    //
    
    std::map<std::string,int> filter;
    std::set<std::string> filter_out;
    
    int n_filtered;
    int n_filtered_singleton;

    

    //
    // dbSNP percentages
    //

    std::string dbSNP_label;
    std::string thousandG_label;

    int dbSNP;
    int dbSNP_singleton;

    int thousandG;
    int thousandG_singleton;

    
    //
    // Functional class
    //

    std::string func_str; 
    std::map<std::string,int> func; 

    int nssnp;
    int nssnp_singleton;
            
  };


  //
  // Individual report: each individual has VStat tracked separately
  //
  
  struct IStat {
    
  IStat( GStore * g ) 
  : g(g) , vstat(g) 
    {
      nvar = 0; 
      Pseq::Util::set_default( vstat );
    }
    
    GStore * g;

    VStat vstat;
    
    std::map<std::string,VStat> stat;

    int nvar;
    std::map<std::string,int> nalt;
    std::map<std::string,int> nhet;
    std::map<std::string,int> nobs;

    void report();
  };

  
  struct GStat {    
  GStat( GStore * g, int group_id , VStat & vstat ) 
  : g(g) , group_id(group_id) , vstat(vstat) { }
    GStore * g;    
    VStat & vstat;
    int group_id;
  };


  struct AuxVDist {
    bool use_binary_phenotype;
    bool match_on_strata;
    bool within_stratum_counts;
    std::map< int2 , int > counts;
    std::map< std::string, std::map<int2,int> > phe_counts;
  };
  
}


//
// 2) Iteration functions
//

class Variant;
class VariantGroup;
void f_vstat( Variant & , void * p);
void g_gstat( VariantGroup & , void * p );
void f_istat( Variant & , void * p);
void f_vdist( Variant & , void * p);  

#endif
