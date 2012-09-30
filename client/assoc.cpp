#include "assoc.h"
#include "plinkseq.h"
#include "util.h"
#include "genic.h"
#include "plinkseq/prob.h"

#include <iostream>
#include <vector>
#include <string>
#include <sstream>

extern GStore g;
extern Pseq::Util::Options args;

bool Pseq::Assoc::variant_assoc_test( Mask & m , 
				      Pseq::Assoc::Aux_vassoc_options & aux , 
				      const Pseq::Util::Options & args )
{


  //
  // Single variant association testsing; assumes a dichotomous or QT
  //
  
  if ( ! ( g.phmap.type() == PHE_DICHOT || g.phmap.type() == PHE_QT ) )
    Helper::halt("v-assoc requires a dichotomous phenotype or quantitive trait");
  

  //
  // Allelic, dominant, recessive tests
  //
  
  const int ntests = 3 ;
  

  //
  // Display options
  //
  
  bool show_ival = args.has( "info" );
  bool show_meta = args.has( "vmeta" );
  
  // Use Yates-chisq (instead of standard, or instead of Fisher's exact if no perms)

  aux.yates_chisq = args.has( "yates" );
  
  // QT testing?
  
  aux.qt = g.phmap.type() == PHE_QT;
  

  //
  // Header row
  //
  
  // T 0 means 'variant level info', e.g. VMETA
  //     that we only list once, and reference
  //     and omnibus test statistic, if any
  //   1 is alternate allele 1
  //   2 is alternate allele 2, etc

  Out & pout = Out::stream( "vassoc" );

  pout.data_reset();
  
  pout.data_group_header( "VAR" );
  pout.data_header( "REF" );
  pout.data_header( "SAMPLES" );
  pout.data_header( "FILTER" );
  pout.data_header( "VMETA" );
  pout.data_header( "CONMETA" );
  pout.data_header( "ALT" );

  pout.data_header( "MAF" );
  pout.data_header( "HWE" );

  if ( aux.qt )
    {
      pout.data_header( "REF" );
      pout.data_header( "HET" );
      pout.data_header( "HOM" );
      
      pout.data_header( "REFN" );
      pout.data_header( "HETN" );
      pout.data_header( "HOMN" );

    }
  else // binary trait
    {
      pout.data_header( "MINA" );
      pout.data_header( "MINU" );
      
      pout.data_header( "OBSA" );
      pout.data_header( "OBSU" );
      
      pout.data_header( "REFA" );
      pout.data_header( "HETA" );
      pout.data_header( "HOMA" );
      
      pout.data_header( "REFU" );
      pout.data_header( "HETU" );
      pout.data_header( "HOMU" );
    }

  pout.data_header( "P" );
  pout.data_header( aux.qt ? "BETA" : "OR" );
  if ( aux.nrep ) pout.data_header( "I" );
  
  pout.data_header( "PDOM" );
  pout.data_header( aux.qt ? "BETADOM" : "ORDOM" );
  if ( aux.nrep ) pout.data_header( "IDOM" );

  pout.data_header( "PREC" );
  pout.data_header( aux.qt ? "BETAREC" : "ORREC" );
  if ( aux.nrep ) pout.data_header( "IREC" );

  pout.data_header_done();

  
  //
  // Set up permutation class, etc
  //
  
  Pseq::Assoc::Aux a;
  a.g     = &g;
  a.rseed = PLINKSeq::DEFAULT_RNG_SEED();
  a.show_info = args.has( "info" );
  g.perm.initiate( aux.nrep , ntests );
  //  if ( args.has("aperm") ) g.perm.adaptive( args.as_int_vector( "aperm" ) );
  a.fix_null_genotypes = args.has( "fix-null" );
  
  //
  // Apply function
  //

  g.vardb.iterate( f_variant_association , &a , m );


  //
  // Post-processing to obtain corrected p-values
  //

  //  pout << Helper::sw( "TEST" , 10 ) 
  // 	    << Helper::sw( "VAR" , 10 )
  // 	    << Helper::sw( "PCORR" , 10 )
  // 	    << "\n";
  
  //   for (int t=0; t < g.perm.n_tests(); t++)
  //     for (int s=0; s < g.perm.n_stats(); s++)
  //       pout << Helper::sw( s , 10 ) 
  // 		<< Helper::sw( t , 10 ) 
  // 		<< Helper::sw( g.perm.max_pvalue(s,t) , 10 ) 
  // 		<< "\n";
  
  return true;
  
}


void f_variant_association( Variant & v , void * p )
{
  
  //
  // Implements single-variant association tests; explicitly testing alternate versus 
  // reference alleles; always presents the genotypic as well as allelic tests; for 
  // binary traits;  uses permutation 
  

  //
  // For now, skip anything other than a basic SNP
  //
  
  if ( ! v.biallelic() ) return;


  //
  // TODO: Copy-number variable markers? 
  //
  
  
  
  //
  // Set up permutations
  //
  
  Pseq::Assoc::Aux * data = (Pseq::Assoc::Aux*)p;
  
  GStore * g = data->g;
  
  const int R = 1 + g->perm.replicates();
  
  g->perm.seed( data->rseed );
  
  g->perm.reset();
  
  if ( data->fix_null_genotypes ) 
    g->perm.fix( VarFunc::missing_genotype_mask( v ) );

  

  //
  // Begin permutations
  //
  
  int obs_a = 0 , obs_u = 0 , obs_tota = 0, obs_totu = 0;

  const int n = v.size();
  
  int tota = 0, totu = 0;
  
  std::set<int> het_carriers;
  std::set<int> hom_carriers;
  std::set<int> missing;

  int obs_refa, obs_heta , obs_homa;
  int obs_refu, obs_hetu , obs_homu;

  double obs_maf;
  double obs_hwe;

  double obs_statistic;
  double obs_statistic_dom;
  double obs_statistic_rec;
      
  double obs_odds;
  double obs_odds_dom;
  double obs_odds_rec;

  double fisher_pv0, fisher_pv1, fisher_pv2;

  for (int p = 0; p < R ; p++ )
    {

      std::vector<double> statistics(3,0);
      
      double & statistic     = statistics[0];
      double & statistic_dom = statistics[1];
      double & statistic_rec = statistics[2];
      
      // Output original data, calculate OR, track alt.allele carriers
      
      if ( p == 0 )
	{

	  // Genotype counts
	  int refa = 0 , refu = 0;
	  int heta = 0 , hetu = 0;
	  int homa = 0 , homu = 0;

	  for (int i=0; i<n; i++)
	    {
	      
	      affType aff = v.ind( g->perm.pos( i ) )->affected();
	      
	      if ( v(i).null() )
		{
		  missing.insert( i );
		  if ( aff == CASE ) ++tota; 
		  else if ( aff == CONTROL ) ++totu;
		}
	      else
		{
		  
		  const int ac = v(i).allele_count( );	      
		  
		  if ( aff == CASE ) 
		    {
		      
		      ++tota;
		      
		      if ( ac == 1 ) 
			{
			  ++heta;
			  het_carriers.insert(i);
			}
		      else if ( ac == 2 ) 
			{
			  ++homa;
			  hom_carriers.insert(i);
			}
		      else if ( ac == 0 ) 
			{
			  ++refa;
			}		      
		    }
		  else if ( aff == CONTROL )
		    {
		      
		      ++totu;
		      
		      if ( ac == 1 ) 
			{
			  ++hetu;
			  het_carriers.insert(i);
			}
		      else if ( ac == 2 ) 
			{
			  ++homu;
			  hom_carriers.insert(i);		      
			}
		      else if ( ac == 0 )
			{
			  ++refu; 
			}
		    }
		}
	    }
      
	  //
	  // Calculate X2 and OR
	  //

	  // Allele counts

	  const int allele_refa = refa*2 + heta;
	  const int allele_refu = refu*2 + hetu;
	  
	  const int allele_alta = homa*2 + heta;
	  const int allele_altu = homu*2 + hetu;
	  
	  if ( R == 1 ) 
	    {
	      if ( ! fisher( Table( allele_refa, allele_refu, allele_alta, allele_altu ) , &fisher_pv0 ) ) fisher_pv0 = 1;
	    }
	  else
	    statistic = Helper::chi2x2( allele_refa ,
					allele_refu , 
					allele_alta , 
					allele_altu );
	  
	  bool zero_count = 
	    allele_refa == 0 || allele_refu == 0 || 
	    allele_alta == 0 || allele_altu == 0;
	  
	  double odds = 0;

	  if ( zero_count ) 
	    obs_odds = ( ( allele_alta + 0.5 ) 
			 * ( allele_refu + 0.5 ) ) 
	      / ( ( ( allele_refa + 0.5 ) 
		    * ( allele_altu + 0.5 ) ) ) ;
	  else
	    obs_odds = (double)( allele_alta * allele_refu ) 
	      / (double)( allele_altu * allele_refa );
	  
	  
	  // Dominant test

	  const int dom_alta = heta + homa;
	  const int dom_altu = hetu + homu;
	  
	  if ( R == 1 ) 
	    {
	      if (! fisher( Table( refa, refu, dom_alta, dom_altu ) , &fisher_pv1 ) ) fisher_pv1 = 1;
	    }
	  else
	    statistic_dom = Helper::chi2x2( refa , 
					    refu , 
					    dom_alta , 
					    dom_altu );
	  
	  zero_count = 
	    refa == 0 || refu == 0 || 
	    dom_alta == 0 || dom_altu == 0;
	
	  double odds_dom = 0;
	  if ( zero_count ) 
	    obs_odds_dom = ( ( dom_alta + 0.5 ) * ( refu + 0.5 ) ) 
	      / ( ( ( refa + 0.5 ) * ( dom_altu + 0.5 ) ) ) ;
	  else
	    obs_odds_dom = (double)( dom_alta * refu ) 
	      / (double)( dom_altu * refa );

	  // Recessive test
	  
	  const int rec_refa = refa + heta;
	  const int rec_refu = refu + hetu;
	  
	  if ( R == 1 ) 
	    {
	      if ( ! fisher( Table( rec_refa, rec_refu, homa, homu ) , &fisher_pv2 ) ) fisher_pv2 = 1;
	    }
	  else
	    statistic_rec = Helper::chi2x2( rec_refa , 
					    rec_refu , 
					    homa , 
					    homu );
	  
	  zero_count = 
	    rec_refa == 0 || rec_refu == 0 || 
	    homa == 0 || homu == 0;

	  double odds_rec = 0;
	  if ( zero_count ) 
	    obs_odds_rec = ( ( homa + 0.5 ) * ( rec_refu + 0.5 ) ) 
	      / ( ( ( rec_refa + 0.5 ) * ( homu + 0.5 ) ) ) ;
	  else
	    obs_odds_rec = (double)( homa * rec_refu ) 
	      / (double)( homu * rec_refa );
	  
	  obs_hwe = Helper::hwe( v );  

	  obs_maf = ( allele_alta + allele_altu ) 
	    / (double) ( allele_alta + allele_altu 
			 + allele_refa + allele_refu );
	  

	  // track for output
	  obs_a = allele_altu < allele_refu ? allele_alta : allele_refa ;
	  obs_u = allele_altu < allele_refu ? allele_altu : allele_refu ;
	  obs_tota = refa + heta + homa;
	  obs_totu = refu + hetu + homu;

	  obs_statistic = statistic;
	  obs_statistic_dom = statistic_dom;
	  obs_statistic_rec = statistic_rec;
	  
	  obs_refa = refa;
	  obs_heta = heta;
	  obs_homa = homa;
	  obs_refu = refu;
	  obs_hetu = hetu;
	  obs_homu = homu;

	}
      else
	{
	  
	  // For permuted datasets, we've already tracked who is an
	  // alternate allele carrier, and who has a missing genotype
	  // -- so only look at these people now

	  int heta = 0 , hetu = 0;
	  int homa = 0 , homu = 0;
	  int misa = 0 , misu = 0;
	  
	  std::set<int>::iterator i = het_carriers.begin();
	  while ( i != het_carriers.end() )
	    {
	      affType aff = v.ind( g->perm.pos( *i ) )->affected();
	      if ( aff == CASE ) ++heta;
	      else if ( aff == CONTROL ) ++hetu;	      
	      ++i;
	    }
	  
	  i = hom_carriers.begin();
	  while ( i != hom_carriers.end() )
	    {
	      affType aff = v.ind( g->perm.pos( *i ) )->affected();
	      if ( aff == CASE ) ++homa;
	      else if ( aff == CONTROL ) ++homu;	      
	      ++i;
	    }

	  i = missing.begin();
	  while ( i != missing.end() )
	    {
	      affType aff = v.ind( g->perm.pos( *i ) )->affected();
	      if ( aff == CASE ) ++misa;
	      else if ( aff == CONTROL ) ++misu;
	      ++i;
	    }
	  
	  
	  int refa = tota - misa - heta - homa;
	  int refu = totu - misu - hetu - homu;
	  
	  // Construct the three tests: allelic, dominant and recessive
	  
	  statistic =     Helper::chi2x2( refa*2+heta , refu*2+hetu , 
					  homa*2+heta , homu*2+hetu );
	  statistic_dom = Helper::chi2x2( refa        , refu        , 
					  homa+heta   , homu+hetu   );
	  statistic_rec = Helper::chi2x2( refa+heta   , refu+hetu   , 
					  homa        , homu        );

	  
	}
         
      //
      // Keep track of statistic, permute
      //

      if ( ! g->perm.score( statistics ) ) break; 
      
    } // Next permutation
  


  //
  // Output empirical p-value (and minimum obtainable, optionally)
  //

  Out & pout = Out::stream( "vassoc" );
  
  pout.data_group( v );

  pout.data( v.reference() );
  pout.data( v.print_samples() );
  pout.data( v.print_meta_filter(";") );
  pout.data( v.meta );
  pout.data( v.consensus.meta );
  
  pout.data( v.alternate() , 1 );
  pout.data( obs_maf , 1 );
  pout.data( obs_hwe , 1 );

  pout.data( obs_a , 1 );
  pout.data( obs_u , 1 );

  pout.data( obs_tota , 1 );
  pout.data( obs_totu , 1 );
  
  pout.data( obs_refa , 1 );
  pout.data( obs_heta , 1 );
  pout.data( obs_homa , 1 );

  pout.data( obs_refu , 1 );
  pout.data( obs_hetu , 1 );
  pout.data( obs_homu , 1 );

  pout.data( R == 1 ? fisher_pv0 : g->perm.pvalue(0) , 1 );
  pout.data( obs_odds , 1 );  
  if ( R != 1 ) pout.data( g->perm.min_pvalue(0) , 1 );
  
  pout.data( R == 1 ? fisher_pv1 : g->perm.pvalue(1) , 1 );
  pout.data( obs_odds_dom , 1 );
  if ( R != 1 ) pout.data( g->perm.min_pvalue(1) , 1 );

  pout.data( R == 1 ? fisher_pv2 : g->perm.pvalue(2) , 1 );
  pout.data( obs_odds_rec , 1 );
  if ( R != 1 ) pout.data( g->perm.min_pvalue(2) , 1 );

  pout.print_data_group();

}



void f_variant_qtassociation( Variant & v , void * p )
{
  

//   //
//   // Implements single-variant association tests; explicitly testing alternate versus 
//   // reference alleles; always presents the genotypic as well as allelic tests; for 
//   // QTs;  uses permutation 
  

//   //
//   // For now, skip anything other than a basic SNP
//   //
  
//   if ( ! v.biallelic() ) return;


//   //
//   // TODO: Copy-number variable markers? 
//   //
  
    
//   //
//   // Set up permutations
//   //
  
//   Pseq::Assoc::Aux * data = (Pseq::Assoc::Aux*)p;
  
//   GStore * g = data->g;
  
//   const int R = 1 + g->perm.replicates();
  
//   g->perm.seed( data->rseed );
  
//   g->perm.reset();
  
//   if ( data->fix_null_genotypes ) 
//     g->perm.fix( VarFunc::missing_genotype_mask( v ) );

  

//   //
//   // Begin permutations
//   //
  
//   int    obs_ref  = 0 , obs_het  = 0 , obs_hom  = 0;
//   double mean_ref = 0 , mean_het = 0 , mean_hom = 0;
//   double var_ref  = 0 , var_het  = 0 , var_hom  = 0;
  
//   const int n = v.size();
  
//   double obs_maf;
//   double obs_hwe;

//   // t-test values
//   // use both empirical and asympotic p-values
//   double obs_statistic;
//   double obs_statistic_dom;
//   double obs_statistic_rec;
  
//   double obs_beta;
//   double obs_beta_dom;
//   double obs_beta_rec;

//   for (int p = 0; p < R ; p++ )
//     {
      
//       std::vector<double> statistics(3,0);
      
//       double & statistic     = statistics[0];
//       double & statistic_dom = statistics[1];
//       double & statistic_rec = statistics[2];
      

//       // Genotype counts
//       int ref = 0; 
//       int het = 0;
//       int hom = 0;
	  
//       // Means
//       double refm = 0.0;
//       double hetm = 0.0;
//       double homm = 0.0;
      
//       for (int i=0; i<n; i++)
// 	    {
	      
// 	      int pos = p == 0 ? i : g->perm.pos( i );
	      
// 	      if ( v.ind( pos )->missing() ) continue;

// 	      if ( v(i).null() )

// 	      double qt = v.ind( g->perm.pos( i ) )->qt();
	      
// 	      if ( v(i).null() )
// 		{
// 		  missing.insert( i );
// 		  if ( aff == CASE ) ++tota; 
// 		  else if ( aff == CONTROL ) ++totu;
// 		}
// 	      else
// 		{
		  
// 		  const int ac = v(i).allele_count( );	      
		  
// 		  if ( aff == CASE ) 
// 		    {
		      
// 		      ++tota;
		      
// 		      if ( ac == 1 ) 
// 			{
// 			  ++heta;
// 			  het_carriers.insert(i);
// 			}
// 		      else if ( ac == 2 ) 
// 			{
// 			  ++homa;
// 			  hom_carriers.insert(i);
// 			}
// 		      else if ( ac == 0 ) 
// 			{
// 			  ++refa;
// 			}		      
// 		    }
// 		  else if ( aff == CONTROL )
// 		    {
		      
// 		      ++totu;
		      
// 		      if ( ac == 1 ) 
// 			{
// 			  ++hetu;
// 			  het_carriers.insert(i);
// 			}
// 		      else if ( ac == 2 ) 
// 			{
// 			  ++homu;
// 			  hom_carriers.insert(i);		      
// 			}
// 		      else if ( ac == 0 )
// 			{
// 			  ++refu; 
// 			}
// 		    }
// 		}
// 	    }
      
// 	  //
// 	  // Calculate X2 and OR
// 	  //

// 	  // Allele counts

// 	  const int allele_refa = refa*2 + heta;
// 	  const int allele_refu = refu*2 + hetu;
	  
// 	  const int allele_alta = homa*2 + heta;
// 	  const int allele_altu = homu*2 + hetu;
	  
// 	  if ( R == 1 ) 
// 	    {
// 	      if ( ! fisher( Table( allele_refa, allele_refu, allele_alta, allele_altu ) , &fisher_pv0 ) ) fisher_pv0 = 1;
// 	    }
// 	  else
// 	    statistic = Helper::chi2x2( allele_refa ,
// 					allele_refu , 
// 					allele_alta , 
// 					allele_altu );
	  
// 	  bool zero_count = 
// 	    allele_refa == 0 || allele_refu == 0 || 
// 	    allele_alta == 0 || allele_altu == 0;
	  
// 	  double odds = 0;

// 	  if ( zero_count ) 
// 	    obs_odds = ( ( allele_alta + 0.5 ) 
// 			 * ( allele_refu + 0.5 ) ) 
// 	      / ( ( ( allele_refa + 0.5 ) 
// 		    * ( allele_altu + 0.5 ) ) ) ;
// 	  else
// 	    obs_odds = (double)( allele_alta * allele_refu ) 
// 	      / (double)( allele_altu * allele_refa );
	  
	  
// 	  // Dominant test

// 	  const int dom_alta = heta + homa;
// 	  const int dom_altu = hetu + homu;
	  
// 	  if ( R == 1 ) 
// 	    {
// 	      if (! fisher( Table( refa, refu, dom_alta, dom_altu ) , &fisher_pv1 ) ) fisher_pv1 = 1;
// 	    }
// 	  else
// 	    statistic_dom = Helper::chi2x2( refa , 
// 					    refu , 
// 					    dom_alta , 
// 					    dom_altu );
	  
// 	  zero_count = 
// 	    refa == 0 || refu == 0 || 
// 	    dom_alta == 0 || dom_altu == 0;
	
// 	  double odds_dom = 0;
// 	  if ( zero_count ) 
// 	    obs_odds_dom = ( ( dom_alta + 0.5 ) * ( refu + 0.5 ) ) 
// 	      / ( ( ( refa + 0.5 ) * ( dom_altu + 0.5 ) ) ) ;
// 	  else
// 	    obs_odds_dom = (double)( dom_alta * refu ) 
// 	      / (double)( dom_altu * refa );

// 	  // Recessive test
	  
// 	  const int rec_refa = refa + heta;
// 	  const int rec_refu = refu + hetu;
	  
// 	  if ( R == 1 ) 
// 	    {
// 	      if ( ! fisher( Table( rec_refa, rec_refu, homa, homu ) , &fisher_pv2 ) ) fisher_pv2 = 1;
// 	    }
// 	  else
// 	    statistic_rec = Helper::chi2x2( rec_refa , 
// 					    rec_refu , 
// 					    homa , 
// 					    homu );
	  
// 	  zero_count = 
// 	    rec_refa == 0 || rec_refu == 0 || 
// 	    homa == 0 || homu == 0;

// 	  double odds_rec = 0;
// 	  if ( zero_count ) 
// 	    obs_odds_rec = ( ( homa + 0.5 ) * ( rec_refu + 0.5 ) ) 
// 	      / ( ( ( rec_refa + 0.5 ) * ( homu + 0.5 ) ) ) ;
// 	  else
// 	    obs_odds_rec = (double)( homa * rec_refu ) 
// 	      / (double)( homu * rec_refa );
	  
// 	  obs_hwe = Helper::hwe( v );  

// 	  obs_maf = ( allele_alta + allele_altu ) 
// 	    / (double) ( allele_alta + allele_altu 
// 			 + allele_refa + allele_refu );
	  

// 	  // track for output
// 	  obs_a = allele_altu < allele_refu ? allele_alta : allele_refa ;
// 	  obs_u = allele_altu < allele_refu ? allele_altu : allele_refu ;
// 	  obs_tota = refa + heta + homa;
// 	  obs_totu = refu + hetu + homu;

// 	  obs_statistic = statistic;
// 	  obs_statistic_dom = statistic_dom;
// 	  obs_statistic_rec = statistic_rec;
	  
// 	  obs_refa = refa;
// 	  obs_heta = heta;
// 	  obs_homa = homa;
// 	  obs_refu = refu;
// 	  obs_hetu = hetu;
// 	  obs_homu = homu;

// 	}
//       else
// 	{
	  
// 	  // For permuted datasets, we've already tracked who is an
// 	  // alternate allele carrier, and who has a missing genotype
// 	  // -- so only look at these people now

// 	  int heta = 0 , hetu = 0;
// 	  int homa = 0 , homu = 0;
// 	  int misa = 0 , misu = 0;
	  
// 	  std::set<int>::iterator i = het_carriers.begin();
// 	  while ( i != het_carriers.end() )
// 	    {
// 	      affType aff = v.ind( g->perm.pos( *i ) )->affected();
// 	      if ( aff == CASE ) ++heta;
// 	      else if ( aff == CONTROL ) ++hetu;	      
// 	      ++i;
// 	    }
	  
// 	  i = hom_carriers.begin();
// 	  while ( i != hom_carriers.end() )
// 	    {
// 	      affType aff = v.ind( g->perm.pos( *i ) )->affected();
// 	      if ( aff == CASE ) ++homa;
// 	      else if ( aff == CONTROL ) ++homu;	      
// 	      ++i;
// 	    }

// 	  i = missing.begin();
// 	  while ( i != missing.end() )
// 	    {
// 	      affType aff = v.ind( g->perm.pos( *i ) )->affected();
// 	      if ( aff == CASE ) ++misa;
// 	      else if ( aff == CONTROL ) ++misu;
// 	      ++i;
// 	    }
	  
	  
// 	  int refa = tota - misa - heta - homa;
// 	  int refu = totu - misu - hetu - homu;
	  
// 	  // Construct the three tests: allelic, dominant and recessive
	  
// 	  statistic =     Helper::chi2x2( refa*2+heta , refu*2+hetu , 
// 					  homa*2+heta , homu*2+hetu );
// 	  statistic_dom = Helper::chi2x2( refa        , refu        , 
// 					  homa+heta   , homu+hetu   );
// 	  statistic_rec = Helper::chi2x2( refa+heta   , refu+hetu   , 
// 					  homa        , homu        );

	  
// 	}
         
//       //
//       // Keep track of statistic, permute
//       //

//       if ( ! g->perm.score( statistics ) ) break; 
      
//     } // Next permutation
  


//   //
//   // Output empirical p-value (and minimum obtainable, optionally)
//   //

//   Out & pout = Out::stream( "vassoc" );
  
//   pout.data_group( v );

//   pout.data( v.reference() );
//   pout.data( v.print_samples() );
//   pout.data( v.print_meta_filter(";") );
//   pout.data( v.meta );
//   pout.data( v.consensus.meta );
  
//   pout.data( v.alternate() , 1 );
//   pout.data( obs_maf , 1 );
//   pout.data( obs_hwe , 1 );

//   pout.data( obs_a , 1 );
//   pout.data( obs_u , 1 );

//   pout.data( obs_tota , 1 );
//   pout.data( obs_totu , 1 );
  
//   pout.data( obs_refa , 1 );
//   pout.data( obs_heta , 1 );
//   pout.data( obs_homa , 1 );

//   pout.data( obs_refu , 1 );
//   pout.data( obs_hetu , 1 );
//   pout.data( obs_homu , 1 );

//   pout.data( R == 1 ? fisher_pv0 : g->perm.pvalue(0) , 1 );
//   pout.data( obs_odds , 1 );  
//   if ( R != 1 ) pout.data( g->perm.min_pvalue(0) , 1 );
  
//   pout.data( R == 1 ? fisher_pv1 : g->perm.pvalue(1) , 1 );
//   pout.data( obs_odds_dom , 1 );
//   if ( R != 1 ) pout.data( g->perm.min_pvalue(1) , 1 );

//   pout.data( R == 1 ? fisher_pv2 : g->perm.pvalue(2) , 1 );
//   pout.data( obs_odds_rec , 1 );
//   if ( R != 1 ) pout.data( g->perm.min_pvalue(2) , 1 );

//   pout.print_data_group();

}







bool Pseq::Assoc::set_assoc_test( Mask & m , const Pseq::Util::Options & args )
{

  //
  // Helper class
  //
  
  Pseq::Assoc::AuxGenic a;

  a.g     = &g;
  a.rseed = PLINKSeq::DEFAULT_RNG_SEED();
  
  int nrep = args.has("perm") ? args.as_int( "perm" ) : -1 ;
    
  //
  // Implement the following gene-based tests
  //

  if ( args.has( "midpoint" ) ) a.show_midbp = true;
  
  //
  // Standard output mode, or dumping a matrix of null-statistics
  // (with originals, and labels, as first row)
  //

  a.dump_stats_matrix = args.has( "dump-null-matrix" );



  //
  // Which tests to apply?
  //

  a.vanilla       =   args.has( "tests" , "sumstat" );
  a.burden        =   args.has( "tests" , "burden" );
  a.uniq          =   args.has( "tests" , "uniq" );
  a.site_burden   =   args.has( "tests" , "site-burden" );
  a.mhit          =   args.has( "tests" , "mhit" );
  a.vt            =   args.has( "tests" , "vt" );
  a.fw            =   args.has( "tests" , "fw" );
  a.calpha        =   args.has( "tests" , "calpha" );
  a.cancor        =   args.has( "tests" , "cancor" );
  a.hoffman_witte =   args.has( "tests" , "stepup" );
  a.kbac          =   args.has( "tests" , "kbac" );
  a.two_hit       =   args.has( "tests" , "two-hit");
  a.skat          =   args.has( "tests" , "skat" );
  a.skato         =   args.has( "tests" , "skato" );


  
  int ntests = a.n_tests() ;
  
  // if no tests specified, default is SKAT

  if ( ntests == 0 ) { a.skat = true;  ntests = 1; } 
  
  // Convenience specification for asymptotic-only tests 

  // If only assymptotic tests, and no permutation requested
  // by default, then add the equivalent of "--perm 0" to the 
  // command line
  
  const int n_asymptotic_tests = a.skat + a.skato;
  
  if ( ntests == n_asymptotic_tests && nrep == -1 ) nrep = 0;

  // labels   vanilla  SUMSTAT
  //          burden   BURDEN
  //          burden2  BURDEN2
  //          uniq     UNIQ
  //          mhit     MHIT
  //          vt       VT
  //          fw       FRQWGT
  //          calpha   CALPHA
  //          cancor   CANCOR
  //          stepup   STEPUP
  //          kbac     KBAC 
  //          two-hit  TWO-HIT
  //          skat     SKAT
  //          skato    SKAT-O


  //
  // Write list of tests, if in matrix-dump mode
  //

  if ( a.dump_stats_matrix ) 
    {
      Out & pmat = Out::stream( "matrix" );

      pmat << nrep ;
      if ( a.burden ) pmat << "\t" << "BURDEN"; 
      if ( a.uniq )   pmat << "\t" << "UNIQ"; 
      if ( a.skat )   pmat << "\t" << "SKAT"; 
      if ( a.skato )  pmat << "\t" << "SKAT-O"; 
      pmat << "\n";
    }

  //
  // Do we have an appropriate pheotype specified?
  //
  
  if ( g.phmap.type() == PHE_DICHOT )
    { 
      // fine for now
    }
  else if ( g.phmap.type() == PHE_QT )
    {
      // for now, only SKAT can handle QTs
      if ( ! ( a.skat || a.skato ) ) 
	Helper::halt( "only SKAT/SKAT-O can handle quantitative traits" );
    }
  else      
    Helper::halt("no dichotomous phenotype specified");
  
  //
  // Get main output stream
  //

  Out & pout = Out::stream( "assoc" );
  Out & pdet = Out::stream( "assoc.det" );


  //
  // Write header for 'assoc' file
  //

  pout.data_reset();
      
  pout.data_group_header( "LOCUS" );
  
  pout.data_header( "POS" );
  
  if ( a.show_midbp ) 
    {
      pout.data_header( "MID" );
      pout.data_header( "BP" );
    }
  
  pout.data_header( "ALIAS" );
  
  pout.data_header( "NVAR" );
  
  pout.data_header( "TEST" );
  
  pout.data_header( "P" );
  pout.data_header( "I" );
  pout.data_header( "DESC" );
  
  pout.data_header_done();

  
  //
  // Set up permutation class
  //
  
  g.perm.initiate( nrep , ntests );
  a.fix_null_genotypes = args.has( "fix-null" );
 
  
  //
  // Covariates (for SKAT)
  //
  
  if ( args.has( "covar" ) )
    {
      Pseq::Assoc::Aux_skat::has_covar = true;
      Pseq::Assoc::Aux_skat::covars = args.as_string_vector( "covar" );
    }


  //
  // Weights?
  //
  
  if ( args.has( "weights" ) )
    {
      a.weights = args.has("weights");
      a.weight_tag = args.as_string( "weights" );
      
      plog >> "Applying variant-weights to BURDEN, UNIQ and SUMSTAT tests, from " >> a.weight_tag >> "\n"; 
      
      Pseq::Assoc::Aux_skat::has_weights = true;
      Pseq::Assoc::Aux_skat::use_freq_weights = false;
      Pseq::Assoc::Aux_skat::weights = args.as_string( "weights" );
    }
  
  else if ( args.has( "skat-weights" ) )
    {
      Pseq::Assoc::Aux_skat::use_freq_weights = true;
      std::vector<double> w = args.as_float_vector( "skat-weights" );
      if ( w.size() != 2 ) Helper::halt( "expecting --skat-weights a b" );
      Pseq::Assoc::Aux_skat::a1 = w[0];
      Pseq::Assoc::Aux_skat::a2 = w[1]; 
    }

  //
  // Initialize TWO-HIT test
  //

  if ( a.two_hit )
    Pseq::Assoc::Aux_two_hit::initialize();

  //
  // Apply tests to dataset
  //
 

  g.vardb.iterate( g_set_association , &a , m );


  //
  // Post-processing to obtain corrected p-values
  //

  if ( false ) // skip for now
    {
      if ( ! a.dump_stats_matrix ) 
	{
	  for (int t=0; t < g.perm.n_tests(); t++)
	    for (int s=0; s < g.perm.n_stats(); s++) 
	      {      
		pout << "_PCORR\t" 
		     << s << "\t"
		     << t << "\t"
		     << g.perm.max_pvalue(s,t) << "\n";      
	      }
	}
    }


  // All done

  return true; 

}



//
// Gene/pathway-based association tests
//

void g_set_association( VariantGroup & vars , void * p )
{


  //
  // Is this a suitable group?
  //
  
  if ( vars.n_individuals() < 2 ) return;  
  if ( vars.size() < 1 ) return;  
  

  //
  // Get auxiliary data
  //

  if ( ! p ) return; 
  Pseq::Assoc::AuxGenic * data = (Pseq::Assoc::AuxGenic*)p;
  GStore * g = data->g;

  //
  // Set up permutations
  //

  const int R = 1 + g->perm.replicates();
  g->perm.seed( data->rseed );
  g->perm.reset();

  if ( data->fix_null_genotypes ) g->perm.fix( VarFunc::missing_genotype_mask( vars ) ) ;

  //
  // Metainfo transport structs
  //


  Pseq::Assoc::Aux_prelim aux_prelim;
  Pseq::Assoc::Aux_burden aux_burden(data->vanilla,data->burden,data->uniq,data->mhit,data->site_burden);
  Pseq::Assoc::Aux_fw_vt  aux_fw_vt(data->fw,data->vt);
  Pseq::Assoc::Aux_calpha aux_calpha;
  Pseq::Assoc::Aux_cancor aux_cancor( 1 , vars.size() , vars.n_individuals() ) ;

  Pseq::Assoc::prelim( vars , &aux_prelim , data );

  
  //  
  // External tests
  //

  Pseq::Assoc::Aux_hoffman_witte aux_hoffman_witte( data->hoffman_witte , vars, &aux_prelim );
  Pseq::Assoc::Aux_kbac aux_kbac;
  Pseq::Assoc::Aux_two_hit aux_two_hit( 1 , vars.size() , vars.n_individuals() );
  Pseq::Assoc::Aux_skat aux_skat;

  //
  // If in dump-stats-matrix mode, this is the second header row -- output vars.nam()
  //

  Out * pmat = NULL;
  Out & pout = Out::stream( "assoc" );

  if ( data->dump_stats_matrix ) 
    {
      pmat = &Out::stream( "matrix" );
      *pmat << vars.name() ;
    }
  
  //
  // Apply tests to original dataset
  //
  
  std::vector<double> test_statistic;
  std::vector<std::string> test_name;
  std::map<std::string,std::string> test_text;



  if ( data->vanilla || data->burden || data->uniq || data->mhit ) 
    {
      
      Pseq::Assoc::stat_burden( vars , 
 				&aux_prelim, 
 				&aux_burden , 
 				&test_text , 
 				true );
      
      if ( data->vanilla ) 
	{ 
	  test_name.push_back("SUMSTAT");
	  test_statistic.push_back( aux_burden.stat_vanilla );
	}
      
      if ( data->burden ) 
	{ 
	  test_name.push_back("BURDEN");
	  test_statistic.push_back( aux_burden.stat_burden );
	  if ( data->dump_stats_matrix ) *pmat << "\t" << aux_burden.stat_burden ;
	}
      
      if ( data->uniq ) 
	{ 
	  test_name.push_back("UNIQ");
	  test_statistic.push_back( aux_burden.stat_uniq );	  
	  if ( data->dump_stats_matrix ) *pmat << "\t" << aux_burden.stat_uniq ;
	}
      
      if ( data->mhit )
	{ 
	  test_name.push_back("MHIT");
	  test_statistic.push_back( aux_burden.stat_mhit );
	}
      
    }
  else
    {

      // at least ensure this is run once, to get the .assoc.det file produced
      // even if burden/uniq/etc is not a specified test; but no need to repeat
      // under the null replicates, etc.

      Pseq::Assoc::stat_burden( vars , 
 				&aux_prelim, 
 				&aux_burden , 
 				&test_text , 
 				true );
    }
  

  if ( data->vt || data->fw )
    {

      aux_fw_vt.fw = data->fw;
      aux_fw_vt.vt = data->vt;

      Pseq::Assoc::stat_fw_vt( vars , 
			       &aux_prelim, 
			       &aux_fw_vt , 
			       &test_text , 
			       true );
      
      if ( data->vt ) 
	{
	  test_name.push_back("VT");
	  test_statistic.push_back( aux_fw_vt.stat_vt );
	}
      
      if ( data->fw )
	{
	  test_name.push_back("FRQWGT");
	  test_statistic.push_back( aux_fw_vt.stat_fw );
	}
      
    }
  

  if ( data->calpha )
    {

      test_name.push_back("CALPHA");

      double statistic = Pseq::Assoc::stat_calpha( vars , 
						   &aux_prelim, 
						   &aux_calpha , 
						   &test_text , 
						   true );
      test_statistic.push_back( statistic ) ;

    }


  if ( data->cancor )
    {
      Pseq::Assoc::stat_cancor( vars , &aux_prelim, &aux_cancor , &test_text , true );
      test_name.push_back( "CANCOR" );
      test_statistic.push_back( aux_cancor.stat );
    }
  

  if ( data->hoffman_witte )
    {
      test_name.push_back( "STEPUP" );
      double statistic = Pseq::Assoc::stat_hoffman_witte( vars , &aux_hoffman_witte , &test_text , true );
      test_statistic.push_back( statistic );
   }
  
  
  if ( data->kbac )
    {
      test_name.push_back( "KBAC" );
      double statistic = Pseq::Assoc::stat_kbac( vars , &aux_prelim , &aux_kbac , &test_text , true );
      test_statistic.push_back( statistic );
    }


  if ( data->two_hit )
    {
      test_name.push_back( "TWO-HIT" );
      double statistic = Pseq::Assoc::stat_two_hit( vars , &aux_prelim , &aux_two_hit , &test_text, true );
      test_statistic.push_back( statistic );
    }
  
  
  if ( data->skat ) 
    {
      test_name.push_back( "SKAT" );
      aux_skat.set_optimal_mode( false );
      double statistic = Pseq::Assoc::stat_skat( vars , &aux_prelim , &aux_skat , &test_text , true ); 
      test_statistic.push_back( statistic );
      if ( data->dump_stats_matrix ) *pmat << "\t" << -log10( statistic );
    }


  if ( data->skato ) 
    {
      test_name.push_back( "SKAT-O" );
      aux_skat.set_optimal_mode( true );
      double statistic = Pseq::Assoc::stat_skat( vars , &aux_prelim , &aux_skat , &test_text , true ); 
      test_statistic.push_back( statistic );
      if ( data->dump_stats_matrix ) *pmat << "\t" << -log10( statistic );
    }



  //
  // Register original test statistics with the permutation class
  //
  
  g->perm.score( test_statistic );

    
  //
  // Begin permutations
  //
  
  
  for (int p = 1; p < R ; p++ )
    {
      std::vector<double> test_statistic;      
      
      if ( data->vanilla || data->burden || data->uniq || data->mhit )
	{	  
	  Pseq::Assoc::stat_burden( vars , &aux_prelim, &aux_burden , NULL , false );	  
	  if ( data->vanilla ) test_statistic.push_back( aux_burden.stat_vanilla );
	  if ( data->burden ) test_statistic.push_back( aux_burden.stat_burden );
	  if ( data->uniq ) test_statistic.push_back( aux_burden.stat_uniq );
	  if ( data->mhit ) test_statistic.push_back( aux_burden.stat_mhit );

	  if ( data->dump_stats_matrix )
	    {
	      if ( data->burden ) *pmat << "\t" << aux_burden.stat_burden;
	      if ( data->uniq ) *pmat << "\t" << aux_burden.stat_uniq;
	    }

	}

      
      if ( data->vt || data->fw )
	{	  
	  Pseq::Assoc::stat_fw_vt( vars , &aux_prelim, &aux_fw_vt , NULL , false ) ;
	  if ( data->vt ) test_statistic.push_back( aux_fw_vt.stat_vt );
	  if ( data->fw ) test_statistic.push_back( aux_fw_vt.stat_fw );
	}
      
      
      if ( data->calpha ) 
	{ 
	  double statistic = Pseq::Assoc::stat_calpha( vars , 
						       &aux_prelim, 
						       &aux_calpha, 
						       NULL , 
						       false );
	  test_statistic.push_back( statistic );
	}
      

      if ( data->cancor )
	{
	  Pseq::Assoc::stat_cancor( vars , &aux_prelim, &aux_cancor , NULL , false );	  
	  test_statistic.push_back( aux_cancor.stat );
	}
      
      
      if ( data->hoffman_witte )
	{
	  double statistic = Pseq::Assoc::stat_hoffman_witte( vars , &aux_hoffman_witte , NULL , false );
	  test_statistic.push_back( statistic );
	}

      
      if ( data->kbac )
	{
	  double statistic = Pseq::Assoc::stat_kbac( vars , &aux_prelim , &aux_kbac , NULL , false );
	  test_statistic.push_back( statistic );
	}
      if ( data->two_hit )
        {
	  double statistic = Pseq::Assoc::stat_two_hit( vars , &aux_prelim , &aux_two_hit , NULL , false );
	  test_statistic.push_back( statistic );
	}

      if ( data->skat )
	{	  
	  double statistic = Pseq::Assoc::stat_skat( vars , &aux_prelim , &aux_skat , NULL , false  );
	  test_statistic.push_back( statistic );
	  if ( data->dump_stats_matrix ) *pmat << "\t" << statistic;
	}
      


      //
      // Store all statistics; permute phenotype labels
      //      

      if ( ! g->perm.score( test_statistic ) ) break;
    
      
    } // Next permutation
  
  


  //
  // Output and return for next gene
  //


  pout.data_group( vars.name() );
  pout.data( vars.coordinate() );
  pout.data( g->locdb.alias( vars.name() , false ) );
  pout.data( vars.size() );
  
  if ( data->show_midbp)
    {
      pout.data( vars.midposition() );
      pout.data( vars.span() );
    }
  
  for (int t=0; t<test_name.size(); t++ )
    {      
      std::string output = test_text[ test_name[t] ];      
      pout.data( test_name[t] , "TEST" , test_name[t] );      
      
      if ( R != 1 ) 
	{
	  pout.data( g->perm.pvalue(t) , "P", test_name[t] );
	  pout.data( g->perm.min_pvalue(t) , "I" , test_name[t] ); 
	}
      else // asymptotic p-values for SKAT, two-hit
	{ 
	  if ( test_name[t] == "SKAT" || test_name[t] == "SKAT-O" ) 
	    {
	      pout.data( aux_skat.returned_pvalue , "P", test_name[t] );
	      pout.data( "." , "I", test_name[t] );
	    }
	  else
	    {
	      if ( test_name[t] == "TWO-HIT" ){
		pout.data( aux_two_hit.returned_pvalue , "P", test_name[t] );
		pout.data( "." , "I", test_name[t] );
	      }
	      else
		{
		  pout.data( "." , "P", test_name[t] );
		  pout.data( "." , "I" , test_name[t] ); 
		}
	    }
	}
      pout.data( output == "" ? "." : output , "DESC" , test_name[t] );
    }
  
  pout.print_data_group();
  
  // if writing null-matrix, add new-line here
  if ( data->dump_stats_matrix ) *pmat << "\n";
  
  return;

}




//
// Per-individual set-enrichment scan
//

struct vidx_t { 
  vidx_t( int chr , int bp ) : chr(chr) , bp(bp) { } 
  int chr;
  int bp;
  bool operator<( const vidx_t & rhs ) const
  {
    if ( chr < rhs.chr ) return true;
    if ( chr > rhs.chr ) return false;
    return bp < rhs.bp;
  }
};

struct aux_var_in_set {
  
  vidx_t   vidx;   // variant index
  int      setid;  // set-name
  
  bool operator<( const aux_var_in_set & rhs ) const
  {
    if ( vidx     <  rhs.vidx ) return true;
    if ( rhs.vidx <  vidx     ) return false;
    return setid < rhs.setid;
   } 
};


struct aux_indiv_enrichment{

  // counts per set/per individual
  std::map<int,double> total_cnts;
  std::map<int,double> total_sqrs;

  std::map<std::string,std::map<int,double> > set_cnts;
  std::map<std::string,std::map<int,double> > set_sqrs;
  

  // per individual, keep track of the rare alleles
  // i.e. for a given set, only count each rare allele once
  //      but keep track of the genes it hits
  //      use variant index to track
  
  // for this individual, have we already seen this variant?
  //  but same variant can be in different genes, that have different set memberships
  //  so need to track what actually gets added?
  
  // for each individual, track which genes have been mapped to which
  // sets i.e. all overlapping genes also get 'tracked' as having
  // mapped to the sets that the index gene gets mapped to.  This
  // means that that gene cannot subsequently be counted for the same
  // set in the same individual. (but it could count towards different
  // sets, i.e. those other than already implicated by the index gene)
  
  // tracked_overlap{indiv}->{gene}->{sets} 
  std::map<int,std::map<uint64_t,std::set<int> > > track_overlap;
  
  // track hits (that also reports overlap), from index
  // track_hits{indiv}->{set}->{vars/genes}
  std::map<int,std::map<int,std::vector<vidx_t> > > track_var_hits;
  std::map<int,std::map<int,std::vector<int> > > track_gene_hits;
  
  // for determining which sets to add to
  std::map<std::string,std::set<int> > * g2s;
  std::vector<std::string> * set_slot;
  std::map<std::string,int> * set_map;

  // locdb group
  int gid;
  std::map<uint64_t,std::string> gene_names;

};



void g_set_enrichment( VariantGroup & vars , void * p )
{
  
  if ( vars.size() == 0 ) return;
  
  int chr = 0, bp1 = 0, bp2 = 0;

  // get location, ensuring maps to a single chromosome interval in this case
  
  if ( ! vars.location( &chr, &bp1, &bp2 ) ) return;
  
  aux_indiv_enrichment * aux = (aux_indiv_enrichment*)p;


  // get ID for this region
  uint64_t index_gene_id = g.locdb.get_region_id( aux->gid , vars.name() );

  // and track the actual name with this ID
  aux->gene_names[ index_gene_id ] = vars.name();

  // for this gene, pull all overlapping genes IDs (which will include 'gene_id' index)
  std::vector<uint64_t> gene_ids = g.locdb.get_region_ids( aux->gid , chr , bp1 , bp2 );

//   std::cout << vars.name() << "\t"
// 	    << index_gene_id << "\t"
// 	    << gene_ids.size() << "\n";

  //
  // Get MAFs for each variant, given sample
  //

  int a = 0;
  const int n = g.indmap.size();
  const int s = vars.size();
  
  std::vector<bool> altmin(s);
  std::vector<double> maf(s);
  for (int j = 0 ; j < s ; j++  )
    {
      int c, c_tot;
      altmin[j] = vars(j).n_minor_allele( &c , &c_tot );
      maf[j] = (double)c/(double)c_tot;
      if ( maf[j] > 0.5 ) maf[j] = 1-maf[j];
      if ( maf[j] == 0 ) maf[j] = 1.0;
    }


  //
  // Primary loop over each individual
  //
  
  for ( int i = 0 ; i < n ; i++ )
    {

      // For this individual, get list of sets that have already been implicated by this gene, 
      // if this gene happened to overlap previous index genes
      
      const std::set<int> already_implicated = aux->track_overlap[ i ][ index_gene_id ]; 
      
      bool observed = false;
      double w = 1; 
      int added = 0;
      
      for (int j = 0 ; j < s ; j++  )
	{	  	  
	  if ( vars(j,i).minor_allele( altmin[j] ) )
	    {	      
	      observed = true;
	      
	      if ( maf[j] < w ) { w = maf[j]; added = j; } 
	      	      
	    }
	}

      //
      // for an individual, we get a (weighted) score which is the 
      // 1/MAF for the lowest frequency variant that have 
      //

      if ( observed ) 
	{

	  // MAF weighted
	  //w = 1.0 / w;
	  w = 1.0;

	  std::set<int> & ss = (*aux->g2s)[ vars.name() ];
	  std::set<int>::iterator ii = ss.begin();
	  while ( ii != ss.end() )
	    {
	      // has this set already been implicated by this gene?
	      
	      if ( already_implicated.find( *ii ) != already_implicated.end() )
		{
// 		  std::cout << "not counting " << index_gene_id << "(" << vars.name() 
// 			    << ") as already implicated for " << *ii << "\n";
		  ++ii;
		  continue; // then just skip
		}

	      //
	      // otherwise, count the 'best' variant towards this set 
	      //
	      
	      aux->set_cnts[ (*aux->set_slot)[ *ii ]][ i ]  += w;
	      aux->set_sqrs[ (*aux->set_slot)[ *ii ]][ i ]  += w * w ;	      	      
	      
	      //
	      // and record that this index, and all overlaps, have
	      // contributed for this person
	      //
	      
	      std::vector<uint64_t>::iterator jj = gene_ids.begin();
	      while ( jj != gene_ids.end() )
		{
		  //		  std::cout << "acknowledging that " << *jj << " now implicates " << *jj << " for " << *ii << "\n"; 
		  aux->track_overlap[i][ *jj ].insert( *ii );
		  ++jj;
		}


	      //
	      // Track for output, the actual 'index variant' hit, and genes for this person/set
	      //
	      
	      aux->track_var_hits[i][*ii].push_back( vidx_t( vars(added).chromosome() , vars(added).position() ) );
	      aux->track_gene_hits[i][*ii].push_back( index_gene_id );	      

	      ++ii;
	    }
	  
	  // and for the background, track exactly the same quantity, but irrespective of set-membership
	  
	  // add as background
	  
	  aux->total_cnts[i] += w ;
	  aux->total_sqrs[i] += w * w ;
	  
	  // track the variants observed here, in relation to the sets they've hit (via this gene)
	  
	}
      
    }
  
}


bool Pseq::Assoc::set_enrich_wrapper( Mask & mask , const Pseq::Util::Options & args )
{
  
  Out & pout = Out::stream( "indiv.enrich" );
  
  //
  // Usage: pseq proj indiv-enrich --mask loc.req={background} 
  //                                      reg.req={background} 
  //                               --loc refseq
  //                               --locset synapse
  //

  // Specifying a background is option.

  // We assume the initial scan has to be based on genes --mask loc.group=refseq
  
  // We need a geneset list to be specified by the argment --locset go
  // and that this relates to a locset that has been loaded into the LOCDB
  // (and that will match in terms of 'refseq')
  
  //
  // Populate a list of individual / gene initially
  //
  
  std::string locset = args.as_string( "locset" );
  std::string loc = args.as_string( "loc" );
  
  
  //
  // Speed-up by not pulling meta-info or subregions here
  //

  bool orig_get_meta = g.locdb.get_meta();
  bool orig_get_subregions = g.locdb.get_subregions();

  g.locdb.get_meta( false );
  g.locdb.get_subregions( false );


  //
  // 'Exclusion mask'
  //

  // the 'background' defines the parts of the experiment in which a
  // variant could feasibly have occurred (i.e. if variants are only
  // called in coding regions, etc; otherwise, the implied region by
  // 'loc' alone would be (for example) all regions from transcription
  // start to stop, not just coding exons, etc.  In this way, we allow
  // for a more flexible joint specification of genes (that link to
  // locus-sets) and the 'background' (that could include
  // experimenter- defined regions of high-coverage, etc) and/or
  // restriction of coding exons (with an optional 50bp flank, for
  // example).  When calculating the effective extent of the gene
  // groups, we extract out anything not in the background.
  
  // As noted above, this is defined by --mask loc.req= and reg.req=


  const std::set<int> & locus_requires     = mask.required_loc();
  const std::set<Region> & region_requires = mask.required_reg();

  std::set<Region> requires = region_requires;
  std::set<int>::iterator il = locus_requires.begin();
  while ( il != locus_requires.end() ) 
    {
      if ( requires.size() == 0 ) 
	requires = g.locdb.get_regions( *il );
      else
	{
	  std::set<Region> rtmp = g.locdb.get_regions( *il );
	  requires = RegionHelper::region_require( requires , rtmp ); 
	}
      ++il;
    }
  
  //plog << "Prior to merging, " << requires.size() << " required intervals\n";
  requires = RegionHelper::region_merge_overlap( requires );  
 
 plog << "\nPost merging, " << requires.size() << " required intervals (reg.req & loc.req masks)\n";

  
  //
  // All genes ( background )
  // 
  
  std::set<Region> all_genes = g.locdb.get_regions( loc );


  //
  // Flatten
  //
  
  plog << "Initial gene/transcript list contains " << all_genes.size() << " entries\n";
  
  all_genes = RegionHelper::region_merge_overlap( all_genes );
  
  plog << "After merging, we now have " << all_genes.size() << " entries\n";


  //
  // only include 'required' regions, if specifed (--mask reg.req=  loc.req= )
  //

  if ( requires.size() > 0 ) 
    all_genes = RegionHelper::region_require( all_genes , requires );
  
  // get total span
  
  double span = RegionHelper::region_span( all_genes );
  
  plog << "After intersecting with masks, " << all_genes.size() << " intervals spanning " << span << " bp\n\n";
  
  
  //
  // Gene-set genes
  //
  
  std::map<std::string,std::set<Region> > targets = g.locdb.fetch_set_regions( loc , locset ); 

  // get gene-name --> gene-sets mapping 
  std::vector<std::string> sets = g.locdb.fetch_set_names( loc , locset );

  plog << sets.size() << " sets found\n";
  
  // track gene -> sets
  
  std::map<std::string,std::set<int> > g2s;
  for (int s=0;s<sets.size();s++)
    {
      std::vector<std::string> members = g.locdb.fetch_set_members( loc , locset , sets[s] );
      plog << "  " << members.size() << " members in " << sets[s] << "\n";
      for (int m = 0 ; m < members.size(); m++ ) 
	g2s[ members[m] ].insert( s );      
    }
  
  plog << "\n";
  

  //
  // get % hitting target by chance  (background rate 'f')
  //

  std::map<std::string,double> f;

  std::map<std::string,std::set<Region> >::iterator ii = targets.begin();
  while ( ii != targets.end() )
    {
      
      std::set<Region> target_genes = RegionHelper::region_merge_overlap( ii->second );
      
      if ( requires.size() > 0 ) 
	target_genes = RegionHelper::region_require( target_genes , requires );
      
      double target_span = RegionHelper::region_span( target_genes );
      
      f[ ii->first ] = ( target_span / span ) ; 
      
      ++ii;
    }


  //
  // return LOCDB to original settings w.r.t. meta-info and subregions
  //
  
  g.locdb.get_meta( orig_get_meta );
  g.locdb.get_subregions( orig_get_meta );
  
  //
  // Now that we've calculated a simple estimate of the prior probability of hitting 
  // a given target (given gene size and the total effective exome being used in the study)
  // now get the empirical counts for the # of variants
  //
  // note: in future, you could imagine other things going into the
  // calculation of 'f', such as base-context  
  //


  //
  // For each individual, calculate a) the total number of variants, b) for each set, the 
  // number of variants in a gene.  Count a gene at a time, and collapse 2+ variants to 1
  // (i.e. to try to help avoiding violations of the independence assumption).
  //
  
  mask.group_loc( loc );
  
  aux_indiv_enrichment aux;
  
  aux.set_slot = &sets;

  // assign a unique numeric ID to each set, for easier tracking, that 
  // can be quickly looked-up
  
  std::map<std::string,int> set_map;
  for (int s=0; s<sets.size(); s++ ) set_map[ sets[s] ] = s; 
  aux.set_map = &set_map;

  // Gene --> set mapping

  aux.g2s = &g2s;

  aux.gid = g.locdb.lookup_group_id( loc );

  if ( aux.gid == 0 ) 
    Helper::halt( "could not find locus-group " + loc );

  plog << "Now iterating over all genes...\n";
  
  // header row
  pout << "ID\tPHE\tSET\tP\tF\tOBS1\tOBS0\tEXP1\tEXP0\tGENES\tVARS\n";
  
  g.vardb.iterate( g_set_enrichment , &aux , mask ) ;
  
  plog << "\n...finished, now calculating per-individual enrichment\n";
  
  

  //
  // Now we should have populated for each set, for each individual, the total 
  // number of genes they have with a variant of 'interest'
  //
  

  //
  // we have populated aux.total_cnts and aux.set_cnts[].  For each individual, for each set, 
  // determine the probability that they have n of m (weighted) counts
  //


  // Use weighted binominal heuristic test, as described here:
  //   http://www.cv.nrao.edu/adass/adassVI/theilerj.html


  // Li & Ma (1983) ratio of Poisson likelihoods method

  // At = 'size of target genes'
  // Ab = 'size of non-target genes'
  // f = At / ( At + Ab ) 

  // S( sqrt( 2 * ( Nt * ln( ( Nt / Et ) ) + nb * ln( ( Nb / Eb ) ) ) ) ) 
  
  // where Et = f ( Nt + Nb ) 
  //       Rb = (1-f) ( Nt + Nb ) 


  // S(s) = 0.5 * ( 1 - erfc( s / sqrt(2) ) )
  // to get 1-tailed p-value

  // Weighted counts:
  //   given sample frequencies (or any other types of arbitary weights), we wish 

  //  Wt  = sum { w_i )       for all in target points
  //  Qt  = sum { w_i ^2 }              "

  // no assumptions about weights summing or averaging to unity. 

  //
  // heuristic for P value
  //
  // At = 'size of target genes'
  // Ab = 'size of non-target genes'
  // f = At / ( At + Ab ) 

  // S = ( sqrt(  ( ( 2 Wt + Wb ) / ( Qt + Qb )  )   *   ( Wt * ln( Wt / Et )  + Wb * ln ( Wb / Eb ) ) ) ) 

  // where Et = f * ( Wt + Wb ) 
  //       Eb = (1-f) * ( Wt + Wb ) 



  const int n = g.indmap.size();
  
  for (int i = 0 ; i < n ; i++ ) 
    {
      
      double w_total = aux.total_cnts[ i ];
      double q_total = aux.total_sqrs[ i ];
      
      for (int s = 0 ; s < sets.size() ; s++ ) 
	{

	  double w_set = aux.set_cnts[ sets[s] ][ i ];
	  double q_set = aux.set_sqrs[ sets[s] ][ i ];
	  
	  double w_bak = w_total - w_set;
	  double q_bak = q_total - q_set;

	  double e_set = f[ sets[s] ] * w_total;
	  double e_bak = w_total - e_set ; 

	  

	  // alpha = f / ( 1 - f )

	  double alpha =  f[ sets[s] ] / (1- f[ sets[s] ] );
	  
	  // can swap in a different approximation from Li & Ma (1983) here:

	  double x = ( w_set - alpha * w_bak ) /  sqrt( alpha * q_set + alpha * q_bak )  ; 	  
	  
// 	  double x = sqrt( ( ( ( 2 * ( w_set + w_bak ) ) ) / 
// 	   		     ( q_set + q_bak ) ) * ( w_set * log( w_set / e_set ) + w_bak * log( w_bak / e_bak ) ) );
	  
	  // one-sided p-value

	  double pv = x < 0 ? 1.0 : 0.5 * ( 1 - erf( x / 1.414214 ) ) ; 

	  //
	  // primary output
	  // 

	  pout  << g.indmap(i)->id() << "\t"
		<< ( g.indmap(i)->affected() == CASE ? "CASE" : g.indmap(i)->affected() == CONTROL ? "CONTROL" : "." ) << "\t"		
		<< sets[s] << "\t"
		<< pv << "\t"
		<< f[ sets[s] ] << "\t"
		<< w_set << "\t" << w_bak << "\t"		
		<< e_set << "\t" << e_bak << "\t";
	  

	  //
	  // gene list
	  //

	  if ( aux.track_gene_hits[i].find(s) != aux.track_gene_hits[i].end() )
	    {
	      std::vector<int> & genes = aux.track_gene_hits[i][s];
	      std::vector<int>::iterator gii = genes.begin();
	      while ( gii != genes.end() )
		{
		  pout << aux.gene_names[ *gii ] << ";";
		  ++gii;
		}
	  
	    }
	  else
	    pout << ".";
	  
	  pout << "\t";
	  
	  
	  //
	  // var list
	  //

	  if ( aux.track_var_hits[i].find(s) != aux.track_var_hits[i].end() )
	    {
	      std::vector<vidx_t> & vars = aux.track_var_hits[i][s];
	      std::vector<vidx_t>::iterator vii = vars.begin();
	      while ( vii != vars.end() )
		{
		  pout << Helper::chrCode( vii->chr ) << ":" << vii->bp << ";";
		  ++vii;
		}
	    }
	  else
	    pout << ".";

	  pout << "\n";


	} // next set
      
    } // next individual
		 

  return true ;


}
