
#ifndef __SUMMARIES_H__
#define __SUMMARIES_H__

#include "func.h"

#include <map>
#include <set>
#include <string>

class GStore;


namespace Pseq 
{
  
  // forward declare IStat
  struct IStat;
  

  // Variant summary statistics that are tracked across the entire
  // (masked) dataset
  
  struct VStat {
    
    template<typename T, typename U> void zero( std::map<T,U> & x ) 
    { 
      typename std::map<T,U>::iterator ii = x.begin();  
      while ( ii != x.end() )  
	{  
	  ii->second = 0;  
	  ++ii;  
	}  
    }

    template<typename T, typename U> void zero2( std::map<std::string,std::map<T,U> > & x ) 
    { 
      typename std::map<std::string,std::map<T,U> >::iterator ii = x.begin();  
      while ( ii != x.end() )
	{
	  typename std::map<T,U>::iterator k = ii->second.begin();  
	  while ( k != ii->second.end() )
	    {
	      k->second = 0;  
	      ++k;
	    }
	  ++ii;
	}
    }
  

    template<typename T, typename U> std::string display_mean( const T & t, const U & u )
    {
      std::stringstream ss;
      if ( u ) ss << t /(double)u;
      else ss << "NA";
      return ss.str();
    }
    
    
    template<typename T, typename U> void display( const std::string & label , std::map<T,U> & x , double denom = 0 ) 
    { 
      typename std::map<T,U>::iterator ii = x.begin();  
      while ( ii != x.end() )  
	{  
	  plog << label << "|"
	       << ii->first << "\t"
	       << ( denom ? ii->second / denom : ii->second ) << "\n";
	  ++ii;  
	}  
    }
    
    template<typename T, typename U> void display( const std::string & label , std::map<T,U> & x , std::map<std::string,long int> & mn_cnt )
    { 
      typename std::map<T,U>::iterator ii = x.begin();  
      while ( ii != x.end() )  
	{  
	  double denom = mn_cnt[ ii->first ];
	  if ( denom )
	    plog << label << "|"
		 << ii->first << "\t"
		 << ii->second / denom << "\n";
	  else
	    plog << label << "|"
		 << ii->first << "\t"
		 << "NA" << "\n";
	  ++ii;  
	}  
    }


    template<typename T, typename U> void headers2( const std::string & label , std::map<std::string,std::map<T,U> > & x ) 
    {
      typename std::map<std::string, std::map<T,U> >::iterator ii = x.begin();
      while ( ii != x.end() )
	{
	  headers( label + "|" + ii->first , ii->second );
	  ++ii;
	}
    }

    template<typename T, typename U> void headers( const std::string & label , std::map<T,U> & x ) 
    { 
      typename std::map<T,U>::iterator ii = x.begin();  
      while ( ii != x.end() )  
	{  
	  plog << "\t" << label << "|"<< ii->first;
	  ++ii;  
	}  
    }


    template<typename T, typename U> void row_display( const std::string & label , std::map<T,U> & x , double denom = 0 ) 
    { 
      typename std::map<T,U>::iterator ii = x.begin();  
      while ( ii != x.end() )  
	{  
	  plog << "\t" << ( denom ? ii->second / denom : ii->second );
	  ++ii;  
	}  
    }
    
    template<typename T, typename U> void row_display( const std::string & label , std::map<T,U> & x , std::map<std::string,long int> & mn_cnt )
    { 
      typename std::map<T,U>::iterator ii = x.begin();  
      while ( ii != x.end() )  
	{  
	  double denom = mn_cnt[ ii->first ];
	  if ( denom )
	    plog << "\t" << ii->second / denom;
	  else
	    plog << "\tNA";
	  ++ii;  
	}  
    }

    template <class T , class U, class V  > void score( std::map<T,U> & x , V v ) 
    { 
      typename std::map<T,U>::iterator i = x.begin(); 
      while ( i != x.end() ) 
 	{ 
 	  if ( i->first.in( v ) ) i->second++; 
 	  ++i; 
 	} 
    } 


    void score_table( std::map<std::string,long int> & x , Variant & var , const std::string & key ) 
    { 
      // assume meta is a single int, flag or textual field that can sensibly be 
      // enumertaed; handle for FLAGs a '1' 
      std::string val;
      if ( var.meta.has_field( key ) )
	{ val = var.meta.as_string( key ); if ( val == "" ) val = "1"; }
      else if ( var.consensus.meta.has_field( key ) )
	{ val = var.consensus.meta.as_string( key ); if ( val == "" ) val = "1"; } 
      if ( val != "" ) x[val]++;
    } 
    

    
    VStat( GStore * g = NULL ) : g(g) 
    {
      // init all variables
      
      nvar = 0;                              // NVAR
      call_rate = 0;                         // RATE
      
      n_singleton = 0;                       // SING
      n_mono = 0;                            // MONOMORPHIC  

      n_pass = 0;                            // FILTER PASS
      n_filter.clear();                      // FILTER xxx
      zero(n_istat_filter);                  // FILTER xxx
      
      n_pass_singleton = 0;                  // FILTER_SING PASS
      n_filter_singleton.clear();            // FILTER_SING xxx

      ti = tv = 0;                           // TITV      
      ti_singleton = tv_singleton = 0;       // TITV_SING
      
      zero( refgrp );                        // REF_xxx
      zero( locgrp );                        // LOC_xxx
      zero( mean_tag );                      // MEAN_xxx
      
      mean_qual = 0;
      cnt_qual = 0;
      mean_dp = 0;
      cnt_dp = 0;
      mean_mac = 0;
      mean_maf = 0;
      cnt_ma = 0;
      
      // G-tags
      zero( mean_gtag );
      zero( cnt_gtag );
      zero( mean_nrgtag );
      zero( cnt_nrgtag );

      zero2( in_range_int_gtag );
      zero2( in_range_dbl_gtag );
      zero2( in_range_int_nrgtag );
      zero2( in_range_dbl_nrgtag );
      
      zero( hwe_failure );                   // HWE

      var_depth_label = PLINKSeq::META_DP();

      // Misc.
      altmin = true;
      single_minor_allele_count = 0;

      // i-stats mode, we have pointer to IStat object to populate per-genotype things
      istat = NULL;
    }
    
    void set_istat( Pseq::IStat * p ) 
    { 
      istat = p; 
    } 
    
   
    void reset()
    {
      
      nvar = 0;                              // NVAR
      call_rate = 0;                         // RATE
      
      n_singleton = 0;                       // SING
      n_mono = 0;                            // MONO

      n_pass = 0;                            // FILTER PASS
      n_filter.clear();                      // FILTER xxx      
      zero(n_istat_filter);                  // FILTER xxx      
      
      n_pass_singleton = 0;                  // FILTER_SING PASS
      n_filter_singleton.clear();            // FILTER_SING xxx

      ti = tv = 0;                           // TITV      
      ti_singleton = tv_singleton = 0;       // TITV_SING
      
      zero( refgrp );                        // REF_xxx
      zero( locgrp );                        // LOC_xxx

      zero( hwe_failure );                   // HWE
      zero( n_mac );                         // MAC
      zero( n_maf );                         // MAF
      
      // Mean TAGs
      zero( mean_tag );                       // MEAN_xxx
      zero( mean_cnt_tag );                   // MEAN_xxx
      
      // Prop. TAGs
      std::map<std::string, std::map<int_range,long int> >::iterator i1 = cnt_itag.begin();
      while ( i1 != cnt_itag.end() )
	{
	  zero( i1->second );                       // MEAN_xxx
	  ++i1;
	}

      // Prop. TAGs
      std::map<std::string, std::map<dbl_range,long int> >::iterator i2 = cnt_ftag.begin();
      while ( i2 != cnt_ftag.end() )
	{
	  zero( i2->second );                       // MEAN_xxx
	  ++i2;
	}
      
      // Table TAGs
      std::map<std::string, std::map<std::string,long int> >::iterator i3 = table_tag.begin();
      while ( i3 != table_tag.end() )
	{
	  zero( i3->second );
	  ++i3;
	}
    
      mean_qual = 0;
      cnt_qual = 0;
      
      mean_dp = 0;
      cnt_dp = 0;
      
      mean_mac = 0;
      mean_maf = 0;
      cnt_ma = 0;

      // G-tags
      zero( mean_gtag );
      zero( cnt_gtag );
      zero2( in_range_int_gtag );
      zero2( in_range_dbl_gtag );

      zero( mean_nrgtag );
      zero( cnt_nrgtag );
      zero2( in_range_int_nrgtag );
      zero2( in_range_dbl_nrgtag );

      // ?? Misc
      altmin = true;
      single_minor_allele_count = 0;
      
    }


    void add_refgroup( const std::string & g ) { refgrp[ g ] = 0; } 
    void add_locgroup( const std::string & g ) { locgrp[ g ] = 0; }

    void add_qual( const std::string & s ) { qual[ dbl_range(s) ] = 0; }
    void add_depth( const std::string & s ) { depth[ int_range(s) ] = 0; }
    
    void add_count_tag( const std::string & s , const std::string & r ) 
    { 
      mType mt = MetaInformation<VarMeta>::type( s );
      if ( mt == META_INT )
	cnt_itag[ s ][ int_range(r) ] = 0;
      else if ( mt == META_FLOAT ) 
	cnt_ftag[ s ][ dbl_range(r) ] = 0;	
    }
    
    void add_mean_tag( const std::string & s ) 
    { 
      mType mt = MetaInformation<VarMeta>::type( s );
      if ( mt == META_FLAG || mt == META_INT || mt == META_BOOL || mt == META_FLOAT )
	{
	  mean_tag[ s ] = 0;
	  mean_cnt_tag[ s ] = 0;
	}
    }

    
    void add_table_tag( const std::string & s ) 
    { 
      std::map<std::string,long int> tbl;
      table_tag[ s ] = tbl;
    }
    
    
    // g_TAGS

    void add_mean_gtag( const std::string & s ) 
    { 
      mType mt = MetaInformation<GenMeta>::type( s );
      if ( mt == META_FLAG || mt == META_INT || mt == META_BOOL || mt == META_FLOAT )
	{
	  mean_gtag[ s ] = 0;
	  cnt_gtag[ s ] = 0;

	  mean_nrgtag[ s ] = 0;
	  cnt_nrgtag[ s ] = 0;

	}
    }

    void add_count_gtag( const std::string & s , const std::string & r ) 
    { 
      mType mt = MetaInformation<GenMeta>::type( s );

      if ( mt == META_INT )
	{
	  in_range_int_gtag[ s ][ int_range(r) ] = 0;
	  cnt_gtag[ s ]  = 0;
	  in_range_int_nrgtag[ s ][ int_range(r) ] = 0;
	  cnt_nrgtag[ s ] = 0;
	}
      else if ( mt == META_FLOAT )
        {
	  in_range_dbl_gtag[ s ][ dbl_range(r) ] = 0;
	  cnt_gtag[ s ] = 0;
	  in_range_dbl_nrgtag[ s ][ dbl_range(r) ] = 0;
	  cnt_nrgtag[ s ] = 0;
	}
    }
    
    void add_hwe_p( const std::string & s ) { hwe_failure[ dbl_range(s) ] = 0; }
    void add_maf( const std::string & s ) { n_maf[ dbl_range(s) ] = 0; }
    void add_mac( const std::string & s ) { n_mac[ int_range(s) ] = 0; }
          
    void add_filter( const std::string & s ) { n_filter[s] = 0; }
    void add_istat_filter( const std::string & s ) { n_istat_filter[s] = 0; }
    

    //
    // Display function
    //
    
    void report();
    void row_headers();
    void row_report( const int , const int , const bool show_genic = false );

    //
    // Data members
    //

    GStore * g ;
    

    // This could apply to one individual or to entire sample
    // The sample will be used to assess things such as MAF, of course
    
    // Simple # of variants, inds
    
    long int nvar;    
    double call_rate;

    long int n_singleton, n_mono;

    // QUAL, DP

    double mean_qual;
    long int cnt_qual;
    
    double mean_dp;
    long int cnt_dp;


    
    
    // FILTERs

    long int n_pass;
    long int n_pass_singleton;
    std::map<std::string,long int> n_filter;
    std::map<std::string,long int> n_filter_singleton;
    std::map<std::string,long int> n_istat_filter;

    // HWE, MAC and MAF

    std::map<dbl_range,long int> hwe_failure;
    std::map<int_range,long int> n_mac;
    std::map<dbl_range,long int> n_maf;
    double mean_maf;
    long int mean_mac;
    long int cnt_ma;
    
    bool altmin;
    int single_minor_allele_count;    
 
    // Depth at called sites

    std::string var_depth_label;
    std::map<int_range,int> depth;
    
    
    // per variant quality score
    
    std::map<dbl_range,long int> qual;


    // Generic variant TAGs 

    std::map<std::string, double> mean_tag;
    std::map<std::string,long int> mean_cnt_tag;

    std::map<std::string, std::map<int_range,long int> > cnt_itag;
    std::map<std::string, std::map<dbl_range,long int> > cnt_ftag;
    
    std::map<std::string, std::map<std::string,long int> > table_tag;

    // G-tags
    std::map<std::string, double> mean_gtag;
    std::map<std::string,long int> cnt_gtag;
    std::map<std::string, std::map<int_range,long int> > in_range_int_gtag;
    std::map<std::string, std::map<dbl_range,long int> > in_range_dbl_gtag;    

    // the following only used in ISTAT as accumulators -- i.e. 
    //  the sums/mans for non-reference gtags only

    std::map<std::string, double> mean_nrgtag;
    std::map<std::string,long int> cnt_nrgtag;
    std::map<std::string, std::map<int_range,long int> > in_range_int_nrgtag;
    std::map<std::string, std::map<dbl_range,long int> > in_range_dbl_nrgtag;    


    // Ti/Tv (all SNPs, and at singletons)
    
    long int ti, tv;
    long int ti_singleton, tv_singleton;
    
    
    // Misc. group counts
    
    std::map<std::string,long int> refgrp;
    std::map<std::string,long int> locgrp;
    
    // Optional pointer to an IStat (for i-stats mode)
    // this will contain an map if ID->vstat for indiv-specific 
    // accumulators

    IStat * istat;
                
  };


  //
  // Individual report: each individual has VStat tracked separately
  //
  
  template<typename T, typename U> void istat_score( std::map<T,U> & v , std::map<T,U> & i )  
  { 
    typename std::map<T,U>::iterator vv = v.begin();  
    while ( vv != v.end() )  
      {  
	i[ vv->first ] += vv->second;
	++vv;  
      }  
  }

  template<typename T, typename U> void istat_score2( std::map<std::string,std::map<T,U> > & v , std::map<std::string,std::map<T,U> > & u )  
  {
    typename std::map<std::string,std::map<T,U> >::iterator i = v.begin();
    while ( i != v.end() )
      {
	Pseq::istat_score( i->second , u[ i->first ] );
	++i;
      }
  }

  
  struct IStat {

    
    IStat( GStore * g ); 
    
    GStore * g;

    VStat vstat;
    
    std::map<std::string,VStat> stat;
    
    bool alt_not_min;

    int nvar;
    std::map<std::string,int> nalt;
    std::map<std::string,int> nmin;
    std::map<std::string,int> nhet;
    std::map<std::string,int> nobs;
    
    // Mean of a particular tag

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
