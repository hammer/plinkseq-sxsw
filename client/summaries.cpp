
#include "summaries.h"
#include "pseq.h"
#include "func.h"

#include <cmath>

extern GStore g;
extern Pseq::Util::Options options;

void f_vstat( Variant & v , void * p) 
{

  Pseq::VStat * vstat = (Pseq::VStat*)p;
  
  // Number of varaits ( denom in RATE )
  
  vstat->nvar++;
  

  // Genotype call rate 
  
  int call_rate = 0;
  const int n = v.size();
  for (int i = 0 ; i < n ; i++)
    if ( v(i).notnull() ) ++call_rate;
  
  vstat->call_rate += (double) call_rate / (double) n;

  
  // Minor allele frequency
  
  int minor_allele_count = 0;
  int total_allele_count = 0;
  
  // track whether alternate is minor allele (used in IStat)
  
  vstat->altmin = v.n_minor_allele( minor_allele_count, total_allele_count );
  
  // track minor allele count (for i-stats)

  vstat->single_minor_allele_count = minor_allele_count; 
  

  // Singleton count 

  bool singleton = minor_allele_count == 1;  
  if ( singleton ) vstat->n_singleton++;
  
  if ( minor_allele_count == 0 ) vstat->n_mono++;

  // HWE test
  
  double hwe_p = Helper::hwe( v );  
  vstat->score( vstat->hwe_failure , hwe_p );

  

  // Read-depth 
  
  bool has_depth = v.consensus.meta.has_field( vstat->var_depth_label );
  if ( has_depth ) 
    {
      int dp = v.consensus.meta.get1_int( vstat->var_depth_label );
      vstat->mean_dp += dp;
      vstat->cnt_dp++;      
      vstat->score( vstat->depth , dp ); 
    }
  
  
  // Quality score
  
  bool has_qual = v.consensus.quality() >= 0;
  if ( has_qual )
    {
      vstat->mean_qual += v.consensus.quality();
      vstat->cnt_qual++;
      vstat->score( vstat->qual , v.consensus.quality() );
    }
  

  // Ti/Tv
  
  if ( v.transition() ) 
    {
      vstat->ti++;
      if ( singleton ) vstat->ti_singleton++;
    }
  else 
    {
      vstat->tv++;
      if ( singleton ) vstat->tv_singleton++;
    }

  
  // FILTERs
  
  bool pass = true;
  std::vector<std::string> fltrs = v.consensus.filters();  
  for (int ft=0; ft < fltrs.size(); ft++)
    {      
      vstat->n_filter[ fltrs[ft] ]++;    
      // this is used only in i-stats (from --options filter=StrandBias, for example
      if ( vstat->n_istat_filter.find( fltrs[ft] ) != vstat->n_istat_filter.end() )
	vstat->n_istat_filter[ fltrs[ft] ]++;
      if ( fltrs[ft] != "PASS" ) pass = false;
    }

  if ( pass )
    {
      ++( vstat->n_pass );
      if ( singleton ) ++(vstat->n_pass_singleton );
    }
  
  // REF groups
  
  std::map<std::string,long int>::iterator ir = vstat->refgrp.begin();
  while ( ir != vstat->refgrp.end() )
    {
      RefVariant refvar = vstat->g->refdb.lookup( v , ir->first );
      if ( refvar.observed() ) ir->second++;
      ++ir;
    }
  

  // LOC groups
  
  ir = vstat->locgrp.begin();
  while ( ir != vstat->locgrp.end() )
    {
      if ( vstat->g->locdb.contains( ir->first , v.chromosome() , v.position() , v.stop() ) )
	ir->second++;
      ++ir;
    }
  
  
  // MAC and MAF

  if ( total_allele_count )
    {
      double maf = minor_allele_count / (double)total_allele_count;   
      vstat->mean_maf += maf;
      vstat->mean_mac += minor_allele_count;      
      vstat->cnt_ma++;	
      vstat->score( vstat->n_mac , minor_allele_count );
      vstat->score( vstat->n_maf , maf );
    }


  // ** Arbitrary TAGs ** 
  
  // Means
  std::map<std::string,double>::iterator i_mean = vstat->mean_tag.begin();
  while ( i_mean != vstat->mean_tag.end() )
    {
      // FLAGs handled separately (i.e. as not present means '0' effectively)
      mType mt = MetaInformation<VarMeta>::type( i_mean->first );
      
      if ( mt == META_FLAG ) 
	{
	  
	  if ( v.meta.flag( i_mean->first ) || v.consensus.meta.flag( i_mean->first ) ) 
	    i_mean->second++;
	  vstat->mean_cnt_tag[ i_mean->first ]++;
	}
      else
	{
	  double val = 0;      
	  bool obs = false;
	  MetaInformation<VarMeta> * m = v.meta.has_field( i_mean->first ) ? &v.meta : NULL;
	  if ( ! m ) m = v.consensus.meta.has_field( i_mean->first ) ? &v.consensus.meta : NULL;
	  if ( m ) 
	    {
	      obs = true;
	      if ( mt == META_INT ) val = m->get1_int( i_mean->first );
	      else if ( mt == META_FLOAT ) val = m->get1_double( i_mean->first );
	      else if ( mt == META_BOOL ) val = m->get1_bool( i_mean->first );
	      else obs = false;
	      if ( obs ) 
		{
		  i_mean->second += val;
		  vstat->mean_cnt_tag[ i_mean->first ]++;
		}
	    }
	}
      ++i_mean;
    }
	  

  // itags
  std::map<std::string,std::map<int_range,long int> >::iterator i_int = vstat->cnt_itag.begin();
  while ( i_int != vstat->cnt_itag.end() )
    {
      int val;
      bool obs = v.meta.has_field( i_int->first );
      if ( obs ) val = v.meta.get1_int( i_int->first );
      if ( (!obs) && v.consensus.meta.has_field( i_int->first ) ) 
	{ obs = true; val = v.consensus.meta.get1_int( i_int->first ); } 
      if ( obs ) vstat->score( i_int->second , val );
	
      ++i_int;
    }

  
  // ftags
  std::map<std::string,std::map<dbl_range,long int> >::iterator i_dbl = vstat->cnt_ftag.begin();
  while ( i_dbl != vstat->cnt_ftag.end() )
    {
      double val;
      bool obs = v.meta.has_field( i_dbl->first );
      if ( obs ) val = v.meta.get1_double( i_dbl->first );
      if ( (!obs) && v.consensus.meta.has_field( i_dbl->first ) ) 
	{ obs = true; val = v.consensus.meta.get1_double( i_dbl->first ); } 
      if ( obs ) vstat->score( i_dbl->second , val );
      ++i_dbl;
    }

    
  // Table tags
  std::map<std::string,std::map<std::string,long int> >::iterator i_table = vstat->table_tag.begin();
  while ( i_table != vstat->table_tag.end() )
    {
      vstat->score_table( i_table->second , v , i_table->first );
      ++i_table;
    }


  // Generic Genotypic Tags, if any 
  // get cnts here; but ensure do not double cnt in itag or ftag by use of cnted
  
  for (int i=0; i<v.size(); i++)
    {
      
      // In VStat mode, we accumulate across all people in single variant
      // In IStat mode, we point to each person's own VStat immediately
      // (and so then no need to sort out later in IStat)
      
      Pseq::VStat * pvstat = vstat->istat ? &( vstat->istat->stat[v.ind(i)->id()])  : vstat ; 
      
      bool has_nonref = v(i).nonreference();


      // gMEANs
      
      std::map<std::string,double>::iterator i_gmean = pvstat->mean_gtag.begin();
      while ( i_gmean != pvstat->mean_gtag.end() )
	{
	  
	  mType mt = MetaInformation<GenMeta>::type( i_gmean->first );
	  
	  if ( v(i).meta.has_field( i_gmean->first ) )
	    {
	      if ( mt == META_INT )
		{
		  i_gmean->second += v(i).meta.get1_int( i_gmean->first );
		  pvstat->cnt_gtag[ i_gmean->first ]++;
 		  if ( has_nonref ) 
		    {
		      pvstat->mean_nrgtag[ i_gmean->first ] += v(i).meta.get1_int( i_gmean->first );
		      pvstat->cnt_nrgtag[ i_gmean->first ]++;
		    }
		}
	      else if ( mt == META_FLOAT ) 
		{
		  i_gmean->second += v(i).meta.get1_double( i_gmean->first );
		  pvstat->cnt_gtag[ i_gmean->first ]++;
 		  if ( has_nonref ) 
		    {
		      pvstat->mean_nrgtag[ i_gmean->first ] += v(i).meta.get1_double( i_gmean->first );
		      pvstat->cnt_nrgtag[ i_gmean->first ]++;
		    }
		}
	      else if ( mt == META_BOOL ) 
		{
		  i_gmean->second += v(i).meta.get1_bool( i_gmean->first );
		  pvstat->cnt_gtag[ i_gmean->first ]++;
		  if ( has_nonref ) 
		    {
		      pvstat->mean_nrgtag[ i_gmean->first ] += v(i).meta.get1_bool( i_gmean->first );
		      pvstat->cnt_nrgtag[ i_gmean->first ]++;
		    }
		}
	      else if ( mt == META_FLAG ) 
		{
		  i_gmean->second += v(i).meta.flag( i_gmean->first );
 		  if ( has_nonref ) 
		    {
		      pvstat->mean_nrgtag[ i_gmean->first ] += v(i).meta.flag( i_gmean->first );		      
		    }

		}	      
	    }
	  
	  if ( mt == META_FLAG ) 
	    {
	      pvstat->cnt_gtag[ i_gmean->first ]++;
	      if ( has_nonref ) 
		pvstat->cnt_nrgtag[ i_gmean->first ]++;
	    }
	  ++i_gmean;
	}


      
      // gTAGs: int_ranges 
      
      std::map<std::string,std::map<int_range,long int> >::iterator i_int = pvstat->in_range_int_gtag.begin();
      while ( i_int != pvstat->in_range_int_gtag.end() )
	{	  
	  if ( v(i).meta.has_field( i_int->first ) )
	    {
	      int x = v(i).meta.get1_int( i_int->first );
	      pvstat->score( i_int->second , x );  	      
	      if ( has_nonref ) 
		pvstat->score( pvstat->in_range_int_nrgtag[ i_int->first ] , x );
	    }	  
	  ++i_int;	  
	}
      
      std::map<std::string,std::map<dbl_range,long int> >::iterator i_dbl = pvstat->in_range_dbl_gtag.begin();
      while ( i_dbl != pvstat->in_range_dbl_gtag.end() )
	{	  
	  if ( v(i).meta.has_field( i_dbl->first ) )
	    {
	      double x = v(i).meta.get1_double( i_dbl->first );
	      pvstat->score( i_dbl->second , x );  	      
	      if ( has_nonref ) 
		pvstat->score( pvstat->in_range_dbl_nrgtag[ i_dbl->first ] , x );
	    }
	  ++i_dbl; 
	}
    }
  
}




  

void g_gstat( VariantGroup & vars, void * p )
{
  
  // Calculate and display relevant per-gene statistics


  Pseq::GStat * aux = (Pseq::GStat*)p;
    
  plog << vars.name() << "\t";
  
  bool empty = vars.size() == 0;

  // get chromosomal range (will be sorted within)

  const int n = vars.size();

  
  // Get original gene -- this assumes that the "grouping" variable
  // was in fact a gene from the LOCDB.  We'll need some flag to show
  // if this is not the case (and so not attempt this step).
  
  Region reg = aux->g->locdb.get_region( aux->group_id , vars.name() );
  int chr = reg.chromosome();
  int bp1 = reg.start.position();
  int bp2 = reg.stop.position();

  int n_exons = reg.subregion.size();
  int n_bp = 0;
  
  // If subregions exist, use these to get BP (i.e. assumes correspond to exons)

  if ( reg.subregion.size() > 0 ) 
    {
      for ( int x = 0 ; x < n_exons ; x++) 
	n_bp += reg.subregion[x].stop.position() - reg.subregion[x].start.position() + 1;    
    }
  else
    {
      n_bp = bp2 - bp1 + 1;
    }
  

  plog << Helper::chrCode( chr ) << ":" 
       << bp1 << ".." << bp2 << "\t";
  
  plog << (bp2-bp1+1)/1000.0 << "\t";
  
  plog << n_bp << "\t"
       << n_exons << "\t";
  
  // rate of variants per bp
  
  plog << n_bp/(double)n << "\t"; 
  

  // Use the variant counter for each gene; also gives 
  // thresholds etc
  
  aux->vstat.reset();


  // some basic means

  double mean_qual = 0;
  double mean_depth = 0;
  
  for (int v = 0 ; v < n ; v++ ) 
    {
      
      // Accumulate statistics for all variants in this gene
      
      f_vstat( vars(v) , &(aux->vstat) );
      
      // means
      
      mean_qual += vars(v).consensus.quality();
      mean_depth += vars(v).consensus.meta.get1_int( aux->vstat.var_depth_label );
      
    }
  
  mean_qual /= (double)n;
  mean_depth /= (double)n;
  

  //
  // Report 
  //

  plog << n << "\t";    

  // # of singletons
  plog << aux->vstat.n_singleton << "\t";
  
  // mean call rate
  plog << aux->vstat.call_rate / (double)n << "\t";

  // mean QUAL, Depth
  
  plog << aux->vstat.mean_qual << "\t"
       << mean_depth << "\t";

  // Ti/Tv
  if ( aux->vstat.tv > 0 )     
    plog << aux->vstat.ti /(double)(aux->vstat.tv) << "\t";
  else
    plog << "NA" << "\t";
  
  // FILTER PASSing
  
  plog << aux->vstat.n_pass << "\t"
       << (double)aux->vstat.n_pass/(double)aux->vstat.nvar << "\t";
  

  // REFGROUPss..

//   plog << aux->vstat.dbSNP << "\t" 

//        << (double)aux->vstat.dbSNP/(double)aux->vstat.nvar << "\t";
  


  // HWE groups

  aux->vstat.display( "HWE" , aux->vstat.hwe_failure );



  
  
  plog << "\n";

  plog.flush();
}


bool Pseq::VarDB::gene_stats_header(Mask & m)
{
  
  plog << "NAME" << "\t"
       << "POS" << "\t"
       << "KB" << "\t"
       << "BP" << "\t"
       << "N_EXON" << "\t"
       << "VRATE" << "\t"
       << "NVAR" << "\t"
       << "NSING" << "\t"
       << "CALLS" << "\t"
       << "QUAL" << "\t"
       << "DEPTH" << "\t"
       << "TI:TV" << "\t"
       << "QC_FAIL" << "\t"
       << "QC_PCT" << "\t"
       << "DBSNP" << "\t"
       << "DBSNP_PCT" << "\t"
       << "THOU_G" << "\t"
       << "THOU_G_PCT" << "\t"
       << "HWE" << "\t"
       << "\n";
  
  plog.flush();
  
 return true;
}


void Pseq::VStat::report()
{
  
  plog << "NVAR" << "\t"
       << nvar  << "\n";
  
  plog << "RATE" << "\t"
       << display_mean( call_rate , nvar ) << "\n";
  
  plog << "MAC\t" 
       << display_mean( mean_mac , cnt_ma ) << "\n";
  display( "MAC" , n_mac );

  plog << "MAF\t" 
       << display_mean( mean_maf , cnt_ma ) << "\n";

  display( "MAF" , n_maf );

  plog << "SING" << "\t"
       << n_singleton << "\n";
  
  plog << "MONO" << "\t"
       << n_mono << "\n";

  if ( tv )
    plog << "TITV" << "\t"
	 << (double)ti/(double)tv << "\n";
  else
    plog << "TITV\tNA\n";

  if ( tv_singleton ) 
    plog << "TITV_S" << "\t"
	 << (double)ti_singleton/(double)tv_singleton << "\n";
  else
    plog << "TITV_S\tNA\n";

  plog << "DP" << "\t"
       << display_mean( mean_dp , cnt_dp ) << "\n";
  
  display( "DP", depth );
  
  plog << "QUAL\t"
       << display_mean( mean_qual , cnt_qual ) << "\n";

  display( "QUAL", qual );
  
  display( "HWE", hwe_failure );

  plog << "PASS\t"
       << n_pass / ( nvar ? (double)nvar : 1.0 ) << "\n";

  display( "FILTER", n_filter , nvar );

  plog << "PASS_S\t"
       << n_pass_singleton / ( n_singleton ? (double)n_singleton : 1.0 ) << "\n";

  display( "FILTER_S", n_filter_singleton , n_singleton );
  
  // REF groups
  
  display( "REF" , refgrp , nvar );
  display( "LOC" , locgrp , nvar );

  
  // Generic tags: means
  display( "MEAN" , mean_tag , mean_cnt_tag );
  
  // Generic tags: int_ranges
  std::map<std::string,std::map<int_range,long int> >::iterator i_int = cnt_itag.begin();
  while ( i_int != cnt_itag.end() )
    { 
      display( i_int->first , i_int->second );
      ++i_int;
    }

  // Generic tags: dbl_ranges
  std::map<std::string,std::map<dbl_range,long int> >::iterator i_dbl = cnt_ftag.begin();
  while ( i_dbl != cnt_ftag.end() )
    { 
      display( i_dbl->first , i_dbl->second );
      ++i_dbl;
    }
  
  // Generic tags: tables
  std::map<std::string,std::map<std::string,long int> >::iterator i_tab = table_tag.begin();
  while ( i_tab != table_tag.end() )
    { 
      display( i_tab->first , i_tab->second );
      ++i_tab;
    }
  

  // Generic tags: genotypic tag means
  
  display( "G" , mean_gtag , cnt_gtag );
  display( "NRG" , mean_nrgtag , cnt_nrgtag );

  // Generic tags: genotypic table-counts; reuse iterators

  // TODO -- make the below four things work.
  
//   i_int = in_range_int_gtag.begin();
//   while ( i_int != in_range_int_gtag.end() )
//     {
//       display( "G|" + i_int->first , i_int->second , cnt_gtag );
//       ++i_int;
//     }
  
//   i_int = in_range_int_nrgtag.begin();
//   while ( i_int != in_range_int_nrgtag.end() )
//     {
//       display( "NRG|" + i_int->first , i_int->second );
//       ++i_int;
//     }

//   i_dbl = in_range_dbl_gtag.begin();
//   while ( i_dbl != in_range_dbl_gtag.end() )
//     {
//       display( "G|" + i_dbl->first , i_dbl->second );
//       ++i_dbl;
//     }

//   i_dbl = in_range_dbl_nrgtag.begin();
//   while ( i_dbl != in_range_dbl_nrgtag.end() )
//     {
//       display( "NRG|" + i_dbl->first , i_dbl->second );
//       ++i_dbl;
//     }

  
}


Pseq::IStat::IStat( GStore * g ) 
  : g(g) , vstat(g) 
{
  nvar = 0; 
  Pseq::Util::set_default( vstat );
  Pseq::IStat * p = this;
  vstat.set_istat( p );
  //copy all all headers etc over
  if (g)
    {
      for (int i=0; i<g->indmap.size(); i++)
	stat[ g->indmap(i)->id() ] = vstat; 	
    }
}


void f_istat( Variant & v , void * p)
{
  
  Pseq::IStat * istat = (Pseq::IStat*)p;
  
  Pseq::VStat * vstat = &(istat->vstat);
  
  // Keep track of total variants looked at
  
  istat->nvar++;
  
  // Calculate statistics for this one variant
  
  vstat->reset();
  
  f_vstat( v , vstat ) ; 
  
  bool singleton = vstat->single_minor_allele_count == 1 ;
  
  
  //
  // Compile equivalent terms into per-individual report
  //

  const int n = v.size();
  
  for (int i=0; i<n; i++)
    {

      // Individual ID
      const std::string id = v.ind(i)->id();

      // Genotpe
      const Genotype & genotype = v(i);
      
      // Accumulate over variants in this individual's slot
      Pseq::VStat & s = istat->stat[ id ];
            
            
      if ( genotype.notnull() )
	{
	
	  // Genotype call-rate 
	  s.call_rate++;
  
	  // 1) For all variant-level statistics, calculate per individual, only 
	  //    considering the variants for which the individual has non-reference 
	  //    a non-reference genotype
	  
	  if ( genotype.minor_allele( vstat->altmin ) )
	    {
	      
	      // HET and HOM alternate genotype counts
	      
	      istat->nalt[ id ]++;
	      if ( genotype.heterozygote() ) istat->nhet[id]++;
	      
	      
	      // SING

	      if ( vstat->n_singleton ) 
		s.n_singleton++;
	      
	      // QUAL

	      if ( vstat->cnt_qual )
		{
		  s.mean_qual += vstat->mean_qual;
		  s.cnt_qual++;
		}
	      
	      // DP
	      
	      if ( vstat->cnt_dp )
		{
		  s.mean_dp += vstat->mean_dp;
		  s.cnt_dp++;
		}
	      
	      // TITV
	      
	      if ( vstat->tv ) s.tv++;
	      else if ( vstat->ti ) s.ti++;

	      // FILTERs
	      
	      if ( vstat->n_pass ) s.n_pass++;
	      if ( vstat->n_pass_singleton ) s.n_pass_singleton++;	      
	      Pseq::istat_score( vstat->n_istat_filter , s.n_istat_filter );

	      // REF & LOC groups	      
	      
	      Pseq::istat_score( vstat->refgrp , s.refgrp );
	      Pseq::istat_score( vstat->locgrp , s.locgrp );

	      // FREQ (HWE, MAC, MAF, etc)

	      Pseq::istat_score( vstat->n_mac, s.n_mac );
	      Pseq::istat_score( vstat->n_maf, s.n_maf );
	      Pseq::istat_score( vstat->hwe_failure, s.hwe_failure );

	      // Generic VTAGS: means, itags, ftags (but not and tables)
	      
	      Pseq::istat_score( vstat->mean_tag , s.mean_tag );
	      Pseq::istat_score( vstat->mean_cnt_tag , s.mean_cnt_tag );
	      Pseq::istat_score2( vstat->cnt_itag , s.cnt_itag );
	      Pseq::istat_score2( vstat->cnt_ftag , s.cnt_ftag );
	      
	    }

	}
      
      
    } // Next individual

}


void Pseq::IStat::report() 
{

  // N=0. nothing to do
  if ( stat.size() == 0 ) return;
  
  // All keys will be 0 or 1, i.e. seen or not seen
  
  plog << "ID" 
       << "\t" << "NVAR" 
       << "\t" << "RATE" 
       << "\t" << "NALT" 
       << "\t" << "NHET" 
       << "\t" << "SING";

  vstat.headers( "MAC", vstat.n_mac );
  vstat.headers( "MAF", vstat.n_maf );
  vstat.headers( "HWE", vstat.hwe_failure );

  plog << "\tTITV";

  plog << "\tPASS"
       << "\tPASS_SING";

  // For v-stats, we will tabulate all FILTERs
  // To make i-stats more manageable, look at a user-specified list only

  vstat.headers( "FILTER", vstat.n_istat_filter );

  plog << "\tQUAL";
  vstat.headers( "QUAL" , vstat.qual );

  plog << "\tDP";
  vstat.headers( "DP" , vstat.depth );

  vstat.headers( "REF", vstat.refgrp );
  vstat.headers( "LOC" , vstat.locgrp );
  vstat.headers( "MEAN" , vstat.mean_tag );
  
  // Generic tags: int_ranges
  vstat.headers2( "" , vstat.cnt_itag );
  vstat.headers2( "" , vstat.cnt_ftag );

  // G-TAGs: means
  vstat.headers( "G" , vstat.mean_gtag );  
  vstat.headers( "NRG" , vstat.mean_nrgtag );  
  
  // G-TAGs: range counts

  vstat.headers2( "G" , vstat.in_range_int_gtag );
  vstat.headers2( "G" , vstat.in_range_dbl_gtag );

  vstat.headers2( "NRG" , vstat.in_range_int_nrgtag );
  vstat.headers2( "NRG" , vstat.in_range_dbl_nrgtag );

  plog << "\n";


  // One row of output per individual

  std::map<std::string,VStat>::iterator i = stat.begin();
  while ( i != stat.end() )
    {

      plog << i->first;
      
      Pseq::VStat & s = i->second;
      
      int actnvar = nalt[i->first]; 
      
      plog << "\t" << s.call_rate 
	   << "\t" << s.call_rate / (double)nvar 
	   << "\t" << actnvar
	   << "\t" << nhet[i->first]
	   << "\t" << s.n_singleton;

      // FREQs 

      s.row_display( "MAC" , s.n_mac , actnvar );
      s.row_display( "MAF" , s.n_maf , actnvar );
      s.row_display( "HWE" , s.hwe_failure , actnvar );

      // Ti/Tv
      if ( s.tv > 0 ) plog << "\t" << (double)s.ti/(double)s.tv ;
      else plog << "\tNA";

      // FILTERs
      plog << "\t" << s.n_pass 
	   << "\t" << s.n_pass_singleton ;

      s.row_display( "FILTER", s.n_istat_filter , actnvar );
      
      // QUAL scores
      
      plog << "\t" << s.display_mean( s.mean_qual , s.cnt_qual );
      s.row_display( "QUAL", s.qual );

      // Read DP

      plog << "\t" << s.display_mean( s.mean_dp , s.cnt_dp );
      s.row_display( "DP", s.qual );

      // REF groups

      s.row_display( "REF" , s.refgrp , nalt[i->first] );
      s.row_display( "LOC" , s.locgrp , nalt[i->first] );
      
      // Generic tags: Means

      s.row_display( "MEAN" , s.mean_tag , s.mean_cnt_tag );
      
      // Generic tags: int_ranges
      
      std::map<std::string,std::map<int_range,long int> >::iterator i_int = s.cnt_itag.begin();
      while ( i_int != s.cnt_itag.end() )
	{ 
	  s.row_display( i_int->first , i_int->second );
	  ++i_int;
	}
      
      // Generic tags: dbl_ranges

      std::map<std::string,std::map<dbl_range,long int> >::iterator i_dbl = s.cnt_ftag.begin();
      while ( i_dbl != s.cnt_ftag.end() )
	{ 
	  s.row_display( i_dbl->first , i_dbl->second );
	  ++i_dbl;
	}


      //
      // 2) Mean per-individual gtags
      //

      s.row_display( "G" , s.mean_gtag , s.cnt_gtag );
      s.row_display( "NRG" , s.mean_nrgtag , s.cnt_nrgtag );

      // Range-counts
      
      i_int = s.in_range_int_gtag.begin();
      while ( i_int != s.in_range_int_gtag.end() )
	{ 
	  s.row_display( i_int->first , i_int->second );
	  ++i_int;
	}
      
      i_dbl = s.in_range_dbl_gtag.begin();
      while ( i_dbl != s.in_range_dbl_gtag.end() )
	{ 
	  s.row_display( i_dbl->first , i_dbl->second );
	  ++i_dbl;
	}

      i_int = s.in_range_int_nrgtag.begin();
      while ( i_int != s.in_range_int_nrgtag.end() )
	{ 
	  s.row_display( i_int->first , i_int->second );
	  ++i_int;
	}
      
      i_dbl = s.in_range_dbl_nrgtag.begin();
      while ( i_dbl != s.in_range_dbl_nrgtag.end() )
	{ 
	  s.row_display( i_int->first , i_int->second );
	  ++i_dbl;
	}

      plog << "\n";
      
      ++i;
    }
}



void f_vdist( Variant & v , void * p)
{
  Pseq::AuxVDist * d = (Pseq::AuxVDist*)p;

  int c     = 0; // minor allele
  int c_tot = 0; // total counts	  
  bool altmin = v.n_minor_allele( c , c_tot );
  
  if ( ( ! d->within_stratum_counts ) &&  c < 1 || c > 4 ) return;
  
  std::map<std::string,int> ca;
  std::map<std::string,int> cu;
  std::set<std::string> grps;

  if ( ! d->match_on_strata ) 
    grps.insert( "." );

  // ignore copy-number for now
  
  std::set<int> carriers;
  for( int i = 0; i< v.size(); i++)
    {
      if ( altmin ? v(i).nonreference() : v(i).reference() )
	{
	  carriers.insert(i);
	  
	  int ac = v(i).minor_allele_count( altmin );
	  
	  if ( d->use_binary_phenotype )
	    {
	      if ( d->match_on_strata ) 
		{

		  std::string label = v.ind(i)->group_label();
		  
		  if ( v.ind(i)->affected() == CASE ) 
		    {		      
		      if ( label != "." ) 
			{
			  ca[ label ] += ac;
			  grps.insert(label);
			}
		    }
		  else if ( v.ind(i)->affected() == CONTROL )
		    {	      
		      if ( label != "." ) 
			{
			  cu[ label ] += ac;
			  grps.insert(label);
			}
		    }
		  
		}	      
	      else // no stratification
		{
		  // Total
		  if ( v.ind(i)->affected() == CASE ) ca["."] += ac;
		  else if ( v.ind(i)->affected() == CONTROL ) cu["."] += ac;
		}
	    }
	}
    }
  
  // track pairwise sharing (including self-self)
  if ( c == 2 ) 
    {
      std::set<int>::iterator i = carriers.begin();
      if ( carriers.size() == 1 ) 
	d->counts[ int2( *i, *i ) ]++;
      else if ( carriers.size() == 2 ) 
	{
	  int i1 = *i;
	  d->counts[ int2( i1, *(++i) ) ]++;
	}
      else if ( carriers.size() > 2 )
	plog.warn( "internal problem in v-dist counting, 2 != 2..." , Helper::int2str( carriers.size() ) );      
    }
  

  // track phenotypic sharing
  
  if ( d->use_binary_phenotype )
    {

      // real-denom here is ca+cu, in case there is missing phenotype data
      // this also ignores doubleton homozygotes
      // i.e. condition on seeing allele (at least once) in N different people
      // note that might will introduce slight bias potentially, check later
      
      if ( d->within_stratum_counts || ( ! d->match_on_strata ) )
	{
	  
	  std::set<std::string>::iterator i = grps.begin();
	  while ( i != grps.end() )
	    {
	      if ( ca[ *i ] + cu[ *i ] >= 1 )
		{
		  d->phe_counts[ *i ][ int2( ca[ *i ] , ca[*i] + cu[ *i ] ) ]++;
		  if ( *i != "." )
		    d->phe_counts[ "." ][ int2( ca[ *i ] , ca[*i] + cu[ *i ] ) ]++;
		}
	      ++i;
	    }
	}

      // Are we defining a singleton in terms of seen once in whole sample, 
      // or seen once in stratum ?   If the former, we need to revise the counts
      
      else 
	{

	  int tot = 0;

	  std::set<std::string>::iterator i = grps.begin();
	  while ( i != grps.end() )
	    {
	      if ( *i != "." )
		tot += ca[*i] + cu[ *i ] ;		
	      ++i;
	    }
	  
	  if ( tot >= 1 ) 
	    {
	      std::set<std::string>::iterator i = grps.begin();	  
	      while ( i != grps.end() )
		{
		  const int stot = ca[ *i ] + cu[ *i ];		  
		  if ( stot == tot  ) // only consider if, e.g. 4 variants, and all 4 are in this strata
		    {	    
		      d->phe_counts[ *i ][ int2( ca[ *i ] , stot ) ]++;
		      if ( *i != "." )
			d->phe_counts[ "." ][ int2( ca[ *i ] , stot ) ]++;
		    }	  
		  ++i;
		}
	    }
	}
    }

}


bool Pseq::VarDB::vdist_summary( Mask & mask )
{
  
  // report # of pairwise sharing of doubletons for each pair, 
  //  report pairs with excessive sharing 

  // { IGNORE FOR NOW... }
  

  // assume a phenotype has been specified; then break down for variants
  // seen 2,3 and 4 times the breakdown of case/case, case/control, control/control
  // counts, and contrast to expectation;  potentially perform this within strata
  
  Pseq::AuxVDist aux;
  aux.use_binary_phenotype = g.phmap.type() == PHE_DICHOT ;
  aux.match_on_strata = g.phmap.strata_set();
  aux.within_stratum_counts = aux.match_on_strata ? ! options.key( "whole-sample-counts" ) : false; 
  
  g.vardb.iterate( f_vdist , &aux , mask );
  
  
  // Report summary
  
  if ( aux.use_binary_phenotype )
    {

      plog.precision(4);

      std::map<std::string, std::map<int2,int> >::iterator s = aux.phe_counts.begin();

      while ( s != aux.phe_counts.end() )
	{
	  
	  std::string strata = s->first;
		
	  // 1,2,3,4
	  int cnt_0 = 0, cnt_1 = 0;
	  int cnt_00 = 0, cnt_01 = 0, cnt_11 = 0;
	  int cnt_000 = 0, cnt_001 = 0, cnt_011 = 0, cnt_111 = 0;
	  int cnt_0000 = 0, cnt_0001 = 0, cnt_0011 = 0, cnt_0111 = 0, cnt_1111 = 0;
	  std::map<int2,int>::iterator i = s->second.begin();
	  while ( i != s->second.end() )
	    {
	      
	      if (  i->first == int2(0,1) ) cnt_0 += i->second;
	      else if (  i->first == int2(1,1) ) cnt_1 += i->second;
	      
	      else if ( i->first == int2(0,2) ) cnt_00 += i->second;
	      else if (  i->first == int2(1,2) ) cnt_01 += i->second;
	      else if (  i->first == int2(2,2) ) cnt_11 += i->second;
	      
	      else if (  i->first == int2(0,3) ) cnt_000 += i->second;
	      else if (  i->first == int2(1,3) ) cnt_001 += i->second;
	      else if (  i->first == int2(2,3) ) cnt_011 += i->second;
	      else if (  i->first == int2(3,3) ) cnt_111 += i->second;
	      
	      else if (  i->first == int2(0,4) ) cnt_0000 += i->second;
	      else if (  i->first == int2(1,4) ) cnt_0001 += i->second;
	      else if (  i->first == int2(2,4) ) cnt_0011 += i->second;
	      else if (  i->first == int2(3,4) ) cnt_0111 += i->second;
	      else if (  i->first == int2(4,4) ) cnt_1111 += i->second;
	      
	      ++i;

	    }
	  
	  // Get expectation
	  
	  // total counts (that might be strata-specific)
	  
	  int ta = 0, tu = 0;
	  for (int i = 0 ; i < g.indmap.size(); i++)
	    {
	      if ( strata == "." || g.indmap(i)->group_label() == strata )
		{
		  if ( g.indmap(i)->affected() == CASE ) ++ta;
		  else if ( g.indmap(i)->affected() == CONTROL ) ++tu;
		}
	    }
	  
 	double pa = ta / (double)(ta + tu);
 	double pu = tu / (double)(ta + tu);
 	int n = ta + tu; // revise N if missing phenotypes
	
 	double e_0 = pu*(cnt_0+cnt_1);
 	double e_1 = pa*(cnt_0+cnt_1);
	
 	double e_00 = pu*pu*(cnt_00+cnt_01+cnt_11); 
 	double e_01 = 2*pa*pu*(cnt_00+cnt_01+cnt_11); 
 	double e_11 = pa*pa*(cnt_00+cnt_01+cnt_11); 
	
 	double e_000 = pu*pu*pu*(cnt_000+cnt_001+cnt_011+cnt_111);
 	double e_001 = 3*pa*pu*pu*(cnt_000+cnt_001+cnt_011+cnt_111);
 	double e_011 = 3*pa*pa*pu*(cnt_000+cnt_001+cnt_011+cnt_111);
 	double e_111 = pa*pa*pa*(cnt_000+cnt_001+cnt_011+cnt_111);
	
 	double e_0000 = pu*pu*pu*pu*(cnt_0000+cnt_0001+cnt_0011+cnt_0111+cnt_1111);
 	double e_0001 = 4*pu*pu*pu*pa*(cnt_0000+cnt_0001+cnt_0011+cnt_0111+cnt_1111);
 	double e_0011 = 6*pu*pu*pa*pa*(cnt_0000+cnt_0001+cnt_0011+cnt_0111+cnt_1111);
 	double e_0111 = 4*pu*pa*pa*pa*(cnt_0000+cnt_0001+cnt_0011+cnt_0111+cnt_1111);
 	double e_1111 = pa*pa*pa*pa*(cnt_0000+cnt_0001+cnt_0011+cnt_0111+cnt_1111);
	
 	double chi1 = (cnt_0-e_0)*(cnt_0-e_0) / e_0
 	  +  (cnt_1-e_1)*(cnt_1-e_1) / e_1;
	
 	double chi2 = (cnt_00-e_00)*(cnt_00-e_00) / e_00
 	  +  (cnt_01-e_01)*(cnt_01-e_01) / e_01
 	  +  (cnt_11-e_11)*(cnt_11-e_11) / e_11;
	
 	double chi3 = (cnt_000-e_000)*(cnt_000-e_000) / e_000
 	  +  (cnt_001-e_001)*(cnt_001-e_001) / e_001
 	  +  (cnt_011-e_011)*(cnt_011-e_011) / e_011
 	  +  (cnt_111-e_111)*(cnt_111-e_111) / e_111;
	
 	double chi4 = (cnt_0000-e_0000)*(cnt_0000-e_0000) / e_0000
 	  +  (cnt_0001-e_0001)*(cnt_0001-e_0001) / e_0001
 	  +  (cnt_0011-e_0011)*(cnt_0011-e_0011) / e_0011
 	  +  (cnt_0111-e_0111)*(cnt_0111-e_0111) / e_0111
 	  +  (cnt_1111-e_1111)*(cnt_1111-e_1111) / e_1111;
	
	
 	plog << "--------------------------- " 
	     << ( strata == "." ? "All cases & controls" : ( g.phmap.strata() + " = " + strata ) )
	     << " ( " 
	     << ta << " cases, " 
	     << tu << " controls ) ---------------------------\n\n";
	
	if ( cnt_0 + cnt_1 > 0 )
	  {
	    std::string w = e_0 < 5 || e_1 < 5 ? " ** low expected counts, p-value not valid ** " : ""; 
	    plog << "Singletons : chi-sq (1df) = " << chi1 << " p = " << Statistics::chi2_prob( chi1 , 1 ) << w << "\n";
	    
	    plog << "  A / U  \tExp\tObs\tRatio\n"
		 << "  0 / 1  \t" << e_0 << "\t" << cnt_0 << "\t" << cnt_0 / e_0 << "\n"
		 << "  1 / 0  \t" << e_1 << "\t" << cnt_1 << "\t" << cnt_1 / e_1 << "\n\n";
	  }

	if ( cnt_00 + cnt_01 + cnt_11 > 0 ) 
	  {
	    std::string w = e_00 < 5 || e_01 < 5 || e_11 < 5 ? " ** low expected counts, p-value not valid ** " : ""; 
	    plog << "Doubletons : chi-sq (2df) = " << chi2 << " p = " << Statistics::chi2_prob( chi2 , 2 ) << w << "\n";

	    plog << "  A / U  \tExp\tObs\tRatio\n"
		 << "  0 / 2  \t" << e_00 << "\t" << cnt_00 << "\t" << cnt_00 / e_00 << "\n"
		 << "  1 / 1  \t" << e_01 << "\t" << cnt_01 << "\t" << cnt_01 / e_01 << "\n"
		 << "  2 / 0  \t" << e_11 << "\t" << cnt_11 << "\t" << cnt_11 / e_11 << "\n\n";
	  }

	if ( cnt_000 + cnt_001 + cnt_011 + cnt_111 > 0 )
	  {
	    std::string w = e_000 < 5 || e_001 < 5 || e_011 < 5 || e_111 < 5 ? " ** low expected counts, p-value not valid ** " : ""; 
	    plog << "Tripletons : chi-sq (3df) = " << chi3 << " p = " << Statistics::chi2_prob( chi3 , 3 ) << w << "\n";
	    
	    plog << "  A / U  \tExp\tObs\tRatio\n"
		 << "  0 / 3  \t" << e_000 << "\t" << cnt_000 << "\t" << cnt_000/e_000 << "\n"
		 << "  1 / 2  \t" << e_001 << "\t" << cnt_001 << "\t" << cnt_001/e_001 << "\n"
		 << "  2 / 1  \t" << e_011 << "\t" << cnt_011 << "\t" << cnt_011/e_011 << "\n"
		 << "  3 / 0  \t" << e_111 << "\t" << cnt_111 << "\t" << cnt_111/e_111 << "\n\n";
	  }

	if ( cnt_0000 + cnt_0001 + cnt_0011 + cnt_0111 + cnt_1111 > 0 )
	  {
	    std::string w = e_0000 < 5 || e_0001 < 5 || e_0011 < 5 || e_0111 < 5 || e_1111 < 5 ? " ** low expected counts, p-value not valid ** " : ""; 
	    plog << "Quadruples : chi-sq (4df) = " << chi4 << " p = " << Statistics::chi2_prob( chi4 , 4 ) << w << "\n";
	    
	    plog << "  A / U  \tExp\tObs\tRatio\n"
		 << "  0 / 4  \t" << e_0000 << "\t" << cnt_0000 << "\t" << cnt_0000/e_0000 << "\n"
		 << "  1 / 3  \t" << e_0001 << "\t" << cnt_0001 << "\t" << cnt_0001/e_0001 << "\n"
		 << "  2 / 2  \t" << e_0011 << "\t" << cnt_0011 << "\t" << cnt_0011/e_0011 << "\n"
		 << "  3 / 1  \t" << e_0111 << "\t" << cnt_0111 << "\t" << cnt_0111/e_0111 << "\n"
		 << "  4 / 0  \t" << e_1111 << "\t" << cnt_1111 << "\t" << cnt_1111/e_1111 << "\n\n";
	  }

 	++s; // next strata
       }
     }

}

