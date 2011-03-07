
#include "summaries.h"
#include "pseq.h"
#include "func.h"

#include <cmath>

extern GStore g;
extern Pseq::Util::Options options;

void f_vstat( Variant & v , void * p) 
{

  Pseq::VStat * vstat = (Pseq::VStat*)p;

  vstat->nvar++;

  //
  // Genotype call rate 
  //

  int call_rate = 0;
  const int n = v.size();
  for (int i = 0 ; i < n ; i++)
    if ( v(i).notnull() ) ++call_rate;


  //
  // Minor allele frequency
  //

  int minor_allele_count = 0;
  int total_allele_count = 0;

  // track whether alternate is minor allele (used in IStat)
  vstat->altmin = v.n_minor_allele( minor_allele_count, 
				    total_allele_count );
  
  // track minor allele count
  vstat->single_minor_allele_count = minor_allele_count; 
  
  int minor_allele_freq_bin100 = total_allele_count > 0 ? 
    (int)(floor( ( (double)minor_allele_count 
		   / (double)total_allele_count ) 
		 * 100)) : 0 ; 
  
  bool singleton = minor_allele_count == 1;

  vstat->call_rate[ call_rate ]++;
  vstat->minor_allele_count[ minor_allele_count ]++;
  vstat->total_allele_count[ total_allele_count ]++;
  if ( total_allele_count > 0 ) 
    vstat->minor_allele_freq_bin100[ minor_allele_freq_bin100 ]++;
  

  //
  // HWE test
  //

  double hwe_p = Helper::hwe( v );  
  std::map<double,int>::iterator i_hwe = vstat->hwe_failure.begin();
  while ( i_hwe != vstat->hwe_failure.end() )
    {
      if ( hwe_p < i_hwe->first ) 
	i_hwe->second++;
      ++i_hwe;
    }


  //
  // Depth at called-variant sites
  //
  
  int depth = v.consensus.meta.get1_int( vstat->var_depth_label );
  std::map<int,int>::iterator i = vstat->depth.begin();
  while ( i != vstat->depth.end() )
    {
      if ( depth >= i->first ) i->second++;
      ++i;
    }  
  
  
  //
  // Quality score per variant
  //

  std::map<double,int>::iterator iq = vstat->qual.begin();
  while ( iq != vstat->qual.end() )
    {
      if ( v.consensus.quality() >= iq->first ) iq->second++;
      ++iq;
    }  
  

  //
  // Genotypic depth and GQ
  //

  // Note: assumes DP present on all individuals -- change this...

  if ( vstat->ind_depth_label != "" )
    {
      double mean_geno_depth = 0;
      for (int i=0; i<n; i++)
	mean_geno_depth += v(i).meta.get1_int( vstat->ind_depth_label );	  
      mean_geno_depth /= (double)n;

      std::map<double,int>::iterator i = vstat->ind_depth.begin();
      while ( i != vstat->ind_depth.end() )
	{
	  if ( mean_geno_depth < i->first ) i->second++;
	  ++i;
	}  
    }

  if ( vstat->ind_qual_label != "" )
    {
      double mean_geno_qual = 0;
      for (int i=0; i<n; i++)
	mean_geno_qual += v(i).meta.get1_double( vstat->ind_qual_label );	  
      mean_geno_qual /= (double)n;

      std::map<double,int>::iterator i = vstat->ind_qual.begin();
      while ( i != vstat->ind_qual.end() )
	{
	  if ( mean_geno_qual >= i->first ) i->second++;
	  ++i;
	}  
    }



  //
  // Ti/Tv
  //

  bool ti =
    ( v.reference() == "A" && v.alternate() == "G" ) ||
    ( v.reference() == "G" && v.alternate() == "A" ) ||
    ( v.reference() == "C" && v.alternate() == "T" ) ||
    ( v.reference() == "T" && v.alternate() == "C" ) ;

  if ( ti ) vstat->ti++;
  else vstat->tv++;

  if ( singleton )
    {
      if ( ti ) vstat->ti_singleton++;
      else vstat->tv_singleton++;
    }

  
  //
  // Filters
  //

  bool f_out = false;
  
  std::vector<std::string> fltrs = v.consensus.filters();

  for (int ft=0; ft < fltrs.size(); ft++)
    {      
      vstat->filter[ fltrs[ft] ]++;      
      std::set<std::string>::iterator j = vstat->filter_out.begin();
      while ( j != vstat->filter_out.end() )
	{
	  if ( fltrs[ft] == *j ) { f_out = true; break; }
	  ++j;
	}
    }

  if ( f_out )
    {
      ++(vstat->n_filtered);
      if ( singleton ) ++(vstat->n_filtered_singleton);
    }


  //
  // dbSNP percentages, 1000G
  //

  RefVariant refdbSNP = vstat->g->refdb.lookup( v , vstat->dbSNP_label );
  RefVariant ref1KG = vstat->g->refdb.lookup( v , vstat->thousandG_label );

  if ( refdbSNP.observed() )
    {
      ++(vstat->dbSNP);
      if ( singleton ) ++(vstat->dbSNP_singleton);
    }

  if ( ref1KG.observed() )
    {
      ++(vstat->thousandG);
      if ( singleton ) ++(vstat->thousandG_singleton);
    }


  //
  // Functional class
  //

  std::string func_str = v.meta.get1_string( vstat->func_str );
  bool func = vstat->func.find( func_str ) != vstat->func.end();
  
  if ( func ) 
    {
      vstat->func[ func_str ]++;
      vstat->nssnp++;
      if ( singleton ) vstat->nssnp_singleton++;    
    }
    
}




  

void g_gstat( VariantGroup & vars, void * p )
{
  
  Pseq::GStat * aux = (Pseq::GStat*)p;
  

  //
  // Calculate and display relevant per-gene statistics
  //
  
  plog << vars.name() << "\t";
  
  
  bool empty = vars.size() == 0;

  // get chromosomal range (will be sorted within)

  const int n = vars.size();

  
  // Get original gene -- this assumes that the "grouping" variable
  // was in fact a gene from the LOCDB.  We'll need some flag to show
  // if this is not the case, and not do it.

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

  
  // TODO: an exon-count, 
  // std::map<int,int> exon_count;
  
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

  double numer = 0, denom = 0;
  
  std::map<int,int>::iterator i_call_rate = aux->vstat.call_rate.begin();
  while ( i_call_rate != aux->vstat.call_rate.end() )
    {
      numer += i_call_rate->first * i_call_rate->second;
      denom += i_call_rate->second;
      ++i_call_rate;
    }

  double mean_callrate = (double)numer/(double)denom;

  //
  // Report 
  //

  plog << n << "\t";

  // # of singletons
  plog << aux->vstat.minor_allele_count[1] << "\t";

  // call rate
  plog << mean_callrate << "\t";

  // mean QUAL, Depth
  plog << mean_qual << "\t"
	    << mean_depth << "\t";

  // Ti/Tv
  if ( aux->vstat.tv > 0 )     
    plog << aux->vstat.ti /(double)(aux->vstat.tv) << "\t"
	      << aux->vstat.ti << "\t"
	      << aux->vstat.tv << "\t";
  else
    plog << "NA" << "\t"
	      << aux->vstat.ti << "\t"
	      << aux->vstat.tv << "\t";

  // Filtered in/out
  
  plog << aux->vstat.n_filtered << "\t"
	    << (double)aux->vstat.n_filtered/(double)aux->vstat.nvar << "\t";

  // dbSNP / 1KG

  plog << aux->vstat.dbSNP << "\t" 
	    << (double)aux->vstat.dbSNP/(double)aux->vstat.nvar << "\t";
  
  plog << aux->vstat.thousandG << "\t"
	    << (double)aux->vstat.thousandG/(double)aux->vstat.nvar << "\t";


  plog << aux->vstat.nssnp << "\t"
	    << aux->vstat.nssnp_singleton << "\t";


  // HWE failures (just print first)

  std::map<double,int>::iterator i_hwe = aux->vstat.hwe_failure.begin();
  while ( i_hwe != aux->vstat.hwe_failure.end() )
    {
      plog << i_hwe->second << "\t";
      break;
    }
  
  
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
       << "TI" << "\t"
       << "TV" << "\t"
       << "QC_FAIL" << "\t"
       << "QC_PCT" << "\t"
       << "DBSNP" << "\t"
       << "DBSNP_PCT" << "\t"
       << "THOU_G" << "\t"
       << "THOU_G_PCT" << "\t"
       << "FUNC" << "\t"
       << "SING_FUNC" << "\t"
       << "HWE" << "\t"
       << "\n";
  
  plog.flush();
  
 return true;
}


void Pseq::VStat::report()
{
  
  plog << "NVAR" << "\t"
       << "ALL" << "\t"
       << nvar  << "\n";
  
  plog << "NVAR" << "\t"
       << "SING" << "\t"
       << minor_allele_count[1] << "\n";
  
  
  // read-depth
  
  std::map<int,int>::iterator i = depth.begin();
  while ( i != depth.end() )
    {
      plog << "DEPTH" << "\t"
	   << i->first << "\t"
	   << i->second << "\n";
      ++i;
    }
  
  std::map<double,int>::iterator i_hwe = hwe_failure.begin();
  while ( i_hwe != hwe_failure.end() )
    {
      plog << "HWE" << "\t"
	   << i_hwe->first << "\t"
	   << i_hwe->second << "\n";
      ++i_hwe;
    }
  
  
  plog << "TI:TV" << "\t"
       << "ALL" << "\t"
       << (double)ti/(double)tv << "\n";
  
  plog << "TI" << "\t"
       << "ALL" << "\t"
       << ti << "\n";
  
  plog << "TV" << "\t"
       << "ALL" << "\t"
       << tv << "\n";
  
  plog << "TI:TV" << "\t"
       << "SING" << "\t"
       << (double)ti_singleton/(double)tv_singleton << "\n";
  
  plog << "TI" << "\t"
       << "SING" << "\t"
       << ti_singleton << "\n";
  
  plog << "TV" << "\t"
       << "SING" << "\t"
       << tv_singleton << "\n";
  
  
  std::map<std::string,int>::iterator i_filter = filter.begin();
  while ( i_filter != filter.end() )
    {
      plog << "FILTER" << "\t"
	   << i_filter->first << "\t"
	   << i_filter->second << "\n";
      ++i_filter;
    }
  
  
  
  plog << "QC_FAIL" << "\t"
       << "ALL" << "\t"
       << n_filtered << "\n";
  
  plog << "QC_FAIL_PCT" << "\t"
       << "ALL" << "\t"
       << (double)n_filtered/(double)nvar << "\n";
  
  
  plog << "QC_FAIL" << "\t"
       << "SING" << "\t"
       << n_filtered_singleton << "\n";
  
  plog << "QC_FAIL_PCT" << "\t"
       << "SING" << "\t"
       << (double)n_filtered_singleton/(double)minor_allele_count[1] << "\n";
  
  
  plog << dbSNP_label << "\t"
       << "ALL" <<"\t"
       << dbSNP << "\n";
  
 plog << dbSNP_label + "_PCT" << "\t"
      << "ALL" << "\t"
      << (double)dbSNP/(double)nvar << "\n";
 
 plog << thousandG_label << "\t"
      << "ALL" << "\t"
      << thousandG << "\n";
 
 plog << thousandG_label +"_PCT" << "\t"
      << "ALL" << "\t"
      << (double)thousandG/(double)nvar << "\n";
 

 plog << dbSNP_label << "\t"
      << "SING" << "\t"
      << dbSNP_singleton << "\n";
 
 plog << dbSNP_label + "_PCT" << "\t"
      << "SING" << "\t"
      << (double)dbSNP_singleton/(double)minor_allele_count[1] << "\n";
 
 plog << thousandG_label << "\t"
      << "SING" << "\t"
      << thousandG_singleton << "\n";
 
 plog << thousandG_label + "_PCT" << "\t"
      << "SING" << "\t"
      << (double)thousandG_singleton/(double)minor_allele_count[1] << "\n";
 
 std::map<std::string,int>::iterator i_func = func.begin();
 while ( i_func != func.end() )
   {
     plog << "FUNC" << "\t"
	  << i_func->first << "\t"
	  << i_func->second << "\n";
     ++i_func;
   }
 
 plog << "FUNC_FLAG" << "\t"
      << "ALL" << "\t"
      << nssnp << "\n";
 
 plog << "FUNC_FLAG" << "\t"
      << "SING" << "\t"
      << nssnp_singleton << "\n";
 
 
 // quality
 
 std::map<double,int>::iterator i_qual = qual.begin();
 while ( i_qual != qual.end() )
   {
     plog << "QUAL" << "\t"
	  << i_qual->first << "\t"
	  << i_qual->second << "\n";
     ++i_qual;
   }

 
 // mean ind. depth
 
 if ( ind_depth_label != "" )
   {
     std::map<double,int>::iterator i = ind_depth.begin();
     while ( i != ind_depth.end() )
       {
	 plog << "MEAN(IND_DP) < " << "\t"
	      << i->first << "\t"
	      << i->second << "\n";
	 ++i;
       }      
     plog << "\n";
   }
 
 // mean ind. GQ
 
 if ( ind_qual_label != "" )
   {
     std::map<double,int>::iterator i = ind_qual.begin();
     while ( i != ind_qual.end() )
       {
	 plog << "MEAN(GQ) >= " << "\t"
	      << i->first << "\t"
	      << i->second << "\n";
	 ++i;
       }            
   }
 
 
 // call rate
 
 std::map<int,int>::iterator i_call_rate = call_rate.begin();
 while ( i_call_rate != call_rate.end() )
   {
     plog << "CALL_RATE" << "\t"
	  << i_call_rate->first << "\t"
	  << i_call_rate->second << "\n";       
     ++i_call_rate;
   }
 
 

 std::map<int,int>::iterator i_minor_allele = minor_allele_count.begin();
 while ( i_minor_allele != minor_allele_count.end() )
   {
     plog << "MA_COUNT" << "\t"
	  << i_minor_allele->first << "\t"
	  << i_minor_allele->second << "\n";       
     ++i_minor_allele;
   }
 
 std::map<int,int>::iterator i_total_allele = total_allele_count.begin();
  while ( i_total_allele != total_allele_count.end() )
    {
      plog << "TA_COUNT" << "\t"
	   << i_total_allele->first << "\t"
	   << i_total_allele->second << "\n";       
      ++i_total_allele;
    }
  
  std::map<int,int>::iterator i_minor_allele_freq = minor_allele_freq_bin100.begin();
  while ( i_minor_allele_freq != minor_allele_freq_bin100.end() )
    {
      plog << "MA_FREQ" << "\t"
	   << i_minor_allele_freq->first << "\t"
	   << i_minor_allele_freq->second << "\n";
      ++i_minor_allele_freq;
    }
   
  
}



void f_istat( Variant & v , void * p)
{
  
  Pseq::IStat * istat = (Pseq::IStat*)p;

  Pseq::VStat * vstat = &(istat->vstat);

  // Keep track of total variants looked at

  istat->nvar++;
  
  //
  // Calculate statistics for this one variant
  //
  
  vstat->reset();
  
  f_vstat( v , vstat ) ; 
  
  bool singleton = vstat->single_minor_allele_count == 1 ;
  bool nssnp = vstat->nssnp = 1;
  
  //
  // Compile equivalent terms into per-individual report
  //

  const int n = v.size();
  
  for (int i=0; i<n; i++)
    {
      const std::string id = v.ind(i)->id();
      
      const Genotype & genotype = v(i);
      
      Pseq::VStat & s = istat->stat[ id ];
      
      //
      // Obtain information on this single variant (for whole sample)
      //
      
      // Call rate
      s.call_rate[ genotype.notnull() ? 1 : 0 ]++;

      // For all observed genotypes
      if ( genotype.notnull() )
	{

	  istat->nobs[id]++;

	  // How we have the minor allele?
	  if ( genotype.minor_allele( vstat->altmin ) )
	    {
	      	      
	      istat->nalt[ id ]++;
	      
	      if ( genotype.heterozygote() )
		istat->nhet[id]++;
	      
	      // Record that this person had a minor allele
	      s.minor_allele_count[ vstat->single_minor_allele_count ]++;
	      
	      // Record this person had a nsSNP
	      if ( vstat->nssnp == 1 ) s.nssnp++;
	      if ( vstat->nssnp_singleton == 1 ) s.nssnp_singleton++;

	      // Record this person had a minor allele of a filtered-out variant
	      if ( vstat->n_filtered == 1 ) { s.n_filtered++; }
	      if ( vstat->n_filtered_singleton == 1 ) s.n_filtered_singleton++;

	      // dbSNP/1kG
	      if ( vstat->dbSNP == 1 ) s.dbSNP++;
	      if ( vstat->thousandG == 1 ) s.thousandG++;
	      
	      // Count number of Ti/Tv
	      if ( vstat->tv == 1 ) s.tv++;
	      else if ( vstat->ti == 1 ) s.ti++;
	      	      
	      
	    }
	  
	}
      
      
      // mean ind-GQ -- track separately
      // mean ind-DP


    }
}

void Pseq::IStat::report() 
{

  // All keys will be 0 or 1, i.e. seen or not seen
  
  plog << "ID" << "\t"
	    << "NVAR" << "\t"
	    << "NOBS" << "\t"
	    << "RATE" << "\t"	    
	    << "NALT" << "\t"
	    << "NHET" << "\t"
	    << "SING" << "\t"
	    << "QCFAIL" << "\t"
	    << "QCFAIL_SING" << "\t"
	    << "TI" << "\t"
	    << "TV" << "\t"
	    << "TI:TV" << "\t"
	    << "DBSNP" << "\t"
	    << "THOU_G" << "\t"
	    << "FUNC" << "\t"
	    << "FUNC_SING" << "\n";

  std::map<std::string,VStat>::iterator i = stat.begin();
  while ( i != stat.end() )
    {
      plog << i->first << "\t";
      
      Pseq::VStat & s = i->second;
      
      plog << nvar << "\t"
		<< nobs[i->first] << "\t"
		<< nobs[i->first] / (double)nvar << "\t"	       				
		<< nalt[i->first] << "\t"		
		<< nhet[i->first] << "\t"
		<< s.minor_allele_count[1] << "\t" 		
		<< s.n_filtered << "\t"
		<< s.n_filtered_singleton << "\t"
		<< s.ti << "\t"
		<< s.tv << "\t";

      if (s.tv > 0 ) plog << (double)s.ti/(double)s.tv << "\t";
      else plog << "NA\t";
      
      plog << s.dbSNP << "\t"
		<< s.thousandG << "\t"
		<< s.nssnp << "\t"
		<< s.nssnp_singleton << "\n";
      
      
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

