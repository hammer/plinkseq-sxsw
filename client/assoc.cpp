#include "assoc.h"
#include "pseq.h"
#include "util.h"
#include "genic.h"

#include <iostream>
#include <vector>
#include <string>
#include <sstream>

extern GStore g;


bool Pseq::Assoc::variant_assoc_test( Mask & m , 
				      Pseq::Assoc::Aux_vassoc_options & aux , 
				      const Pseq::Util::Options & args )
{
  

  //
  // Single variant association testsing; assumes a dichotomous phenotype; QTs 
  // are redirected to a different function (below)
  //
  
  if ( g.phmap.type() == PHE_QT ) 
    return variant_qtassoc_test( m , aux , args );
  else if ( ! g.phmap.type() == PHE_DICHOT ) 
    Helper::halt("basic association tests assumes a dichotomous phenotype");
  

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
  

  //
  // Header row
  //
  
  // T 0 means 'variant level info', e.g. VMETA
  //     that we only list once, and reference
  //     and omnibus test statistic, if any
  //   1 is alternate allele 1
  //   2 is alternate allele 2, etc

  plog.data_reset();
  
  plog.data_group_header( "VAR" );
  plog.data_header( "REF" );
  plog.data_header( "SAMPLES" );
  plog.data_header( "FILTER" );
  plog.data_header( "VMETA" );
  plog.data_header( "CONMETA" );
  plog.data_header( "ALT" );

  plog.data_header( "MAF" );
  plog.data_header( "HWE" );

  plog.data_header( "MINA" );
  plog.data_header( "MINU" );

  plog.data_header( "OBSA" );
  plog.data_header( "OBSU" );

  plog.data_header( "REFA" );
  plog.data_header( "HETA" );
  plog.data_header( "HOMA" );

  plog.data_header( "REFU" );
  plog.data_header( "HETU" );
  plog.data_header( "HOMU" );

  plog.data_header( "P" );
  plog.data_header( "OR" );
  if ( aux.nrep ) plog.data_header( "I" );

  plog.data_header( "PDOM" );
  plog.data_header( "ORDOM" );
  if ( aux.nrep ) plog.data_header( "IDOM" );

  plog.data_header( "PREC" );
  plog.data_header( "ORREC" );
  if ( aux.nrep ) plog.data_header( "IREC" );

  plog.data_header_done();

  
  //
  // Set up permutation class, etc
  //
  
  Pseq::Assoc::Aux a;
  a.g     = &g;
  a.rseed = time(0);
  a.show_info = args.has( "info" );
  g.perm.initiate( aux.nrep , ntests );
  //  if ( args.has("aperm") ) g.perm.adaptive( args.as_int_vector( "aperm" ) );
  a.fix_null_genotypes = args.has("fix-null");


  //
  // Apply function
  //

  g.vardb.iterate( f_variant_association , &a , m );


  //
  // Post-processing to obtain corrected p-values
  //

  //  plog << Helper::sw( "TEST" , 10 ) 
  // 	    << Helper::sw( "VAR" , 10 )
  // 	    << Helper::sw( "PCORR" , 10 )
  // 	    << "\n";
  
  //   for (int t=0; t < g.perm.n_tests(); t++)
  //     for (int s=0; s < g.perm.n_stats(); s++)
  //       plog << Helper::sw( s , 10 ) 
  // 		<< Helper::sw( t , 10 ) 
  // 		<< Helper::sw( g.perm.max_pvalue(s,t) , 10 ) 
  // 		<< "\n";
  
  return true;
  
}


void f_variant_association( Variant & v , void * p )
{
  
  // Implements single-variant association tests; assumes a dominant
  // model based on non-ref alleles (i.e. rather than minor allele)

  
  //
  // For now, skip anything other than a basic SNP
  //
  
  // Handle multi-allelic markers in separate funcion
    
  if ( ! v.biallelic() ) return;

  //
  // Copy-number variable markers? 
  //

  // TODO ... (?) 
  
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
  
  plog.data_group( v );

  plog.data( v.reference() );
  plog.data( v.print_samples() );
  plog.data( v.print_meta_filter() );
  plog.data( v.meta );
  plog.data( v.consensus.meta );
  
  plog.data( v.alternate() , 1 );
  plog.data( obs_maf , 1 );
  plog.data( obs_hwe , 1 );

  plog.data( obs_a , 1 );
  plog.data( obs_u , 1 );

  plog.data( obs_tota , 1 );
  plog.data( obs_totu , 1 );
  
  plog.data( obs_refa , 1 );
  plog.data( obs_heta , 1 );
  plog.data( obs_homa , 1 );

  plog.data( obs_refu , 1 );
  plog.data( obs_hetu , 1 );
  plog.data( obs_homu , 1 );

  plog.data( R == 1 ? fisher_pv0 : g->perm.pvalue(0) , 1 );
  plog.data( obs_odds , 1 );  
  if ( R != 1 ) plog.data( g->perm.min_pvalue(0) , 1 );
  
  plog.data( R == 1 ? fisher_pv1 : g->perm.pvalue(1) , 1 );
  plog.data( obs_odds_dom , 1 );
  if ( R != 1 ) plog.data( g->perm.min_pvalue(1) , 1 );

  plog.data( R == 1 ? fisher_pv2 : g->perm.pvalue(2) , 1 );
  plog.data( obs_odds_rec , 1 );
  if ( R != 1 ) plog.data( g->perm.min_pvalue(2) , 1 );

  plog.print_data_group();

}




//
// Quantitative traits
//


bool Pseq::Assoc::variant_qtassoc_test( Mask & m , 
					Pseq::Assoc::Aux_vassoc_options & aux , 
					const Pseq::Util::Options & args )
{
  

//   //
//   // Allelic, dominant, recessive tests
//   //
  
//   const int ntests = 3 ;
  

//   //
//   // Display options
//   //
  
//   bool show_meta = args.has("vmeta");


//   //
//   // Header row
//   //
  
//   // T 0 means 'variant level info', e.g. VMETA
//   //     that we only list once, and reference
//   //     and omnibus test statistic, if any

//   //   1 is alternate allele 1

//   //   2 is alternate allele 2, etc


//   plog.data_reset();
  
//   plog.data_group_header( "VAR" );
//   plog.data_header( "REF" );
//   plog.data_header( "SAMPLES" );
//   plog.data_header( "FILTER" );
//   plog.data_header( "VMETA" );
//   plog.data_header( "CONMETA" );
//   plog.data_header( "ALT" );

//   plog.data_header( "MAF" );
//   plog.data_header( "HWE" );

//   plog.data_header( "REFMEAN" );
//   plog.data_header( "REFSD" );
//   plog.data_header( "REFOBS" );
  
//   plog.data_header( "HETMEAN" );
//   plog.data_header( "HETSD" );
//   plog.data_header( "HETOBS" );
  
//   plog.data_header( "HOMMEAN" );
//   plog.data_header( "HOMSD" );
//   plog.data_header( "HOMOBS" );

//   plog.data_header( "P" );
//   plog.data_header( "BETA" );

//   plog.data_header( "PDOM" );
//   plog.data_header( "BETADOM" );

//   plog.data_header( "PREC" );
//   plog.data_header( "BETAREC" );
  
//   plog.data_header_done();

  
//   //
//   // Set up permutation class, etc
//   //
  
//   Pseq::Assoc::Aux a;
//   a.g     = &g;
//   a.rseed = time(0);
//   g.perm.initiate( aux.nrep , ntests );
//   a.fix_null_genotypes = args.has("fix-null");

  
  
//   // Set phenotype

//   Pseq::Assoc::Aux_glm aux;
  
//   const int n = g.indmap.size();
//   aux.y.resize( n );
//   aux.mask.resize( n , false );
  
//   for (int i=0; i < n; i++)
//     {      

//       Individual * person = g.indmap.ind(i);
      
//       if ( person->missing() )
// 	aux.mask[i] = true; //mask out
//       else
// 	aux.y[i] = person->qt();

//     }


//   //
//   // Apply function
//   //

//   g.vardb.iterate( f_variant_qtassociation , &aux , m );


//   //
//   // Post-processing to obtain corrected p-values
//   //

//   //  plog << Helper::sw( "TEST" , 10 ) 
//   // 	    << Helper::sw( "VAR" , 10 )
//   // 	    << Helper::sw( "PCORR" , 10 )
//   // 	    << "\n";
  
//   //   for (int t=0; t < g.perm.n_tests(); t++)
//   //     for (int s=0; s < g.perm.n_stats(); s++)
//   //       plog << Helper::sw( s , 10 ) 
//   // 		<< Helper::sw( t , 10 ) 
//   // 		<< Helper::sw( g.perm.max_pvalue(s,t) , 10 ) 
//   // 		<< "\n";
  
  return true;
  
}



void f_variant_qtassociation( Variant & v , void * p )
{
  
  // Implements single-variant quantitative trait association tests
  // Uses univariante linear regression to calculate a quick beta and P, etc

  //
  // For now, skip anything other than a basic SNP.
  //
  
  if ( ! v.biallelic() ) return;
  
  //
  // Copy-number variable markers? 
  //

  // TODO ... (?) 

  
//   Pseq::Assoc::Aux_glm * data = (Pseq::Assoc::Aux_glm *)p;
  

//   // Extract out non-missing phenotypes
//   // (in future, set up a MASK on the MATRIX)
  
//   // For a given variant, get non-missing n and 
//   // create phenotypes
  
//   const int n = v.size();
//   const int np = 1;  // num of (genetic) parameters

//   // actual number of non-missing individuals

//   int an = n;

//   std::vector<bool> mask( n , false ); // F means include
  
//   for (int i=0; i<n; i++)
//     {
//       if ( data->mask[i] ) { --an; mask[i] = true; }
//       else if ( data->c.masked(i) ) { --an; mask[i] = true; }
//       else if ( data->use_dosage && ! v(i).meta.has_field( data->softtag ) ) { --an; mask[i] = true; }
//       else if ( data->use_postprobs && ( ( ! v(i).meta.has_field( data->softtag ) ) || v(i).meta.get_double( data->softtag ).size() != 3 ) ) { --an; mask[i] = true; }
//       else if ( v(i).null() ) { --an; mask[i] = true; }
//     }

//   Data::Vector<double> y( an );

//   Data::Matrix<double> x( an , 2 );  // 2 = intercept plus basic encoding

//   // Populate

//   int ni = 0;

//   for (int i=0; i<n; i++)
//     {
  
//       if ( ! mask[i] ) 
// 	{
	
// 	  // DV
// 	  y[ni] = data->y[i];
	  
// 	  // Intercept
// 	  x(ni,0) = 1;
	  
// 	  // Genotype
// 	  if ( data->use_dosage ) // dosage of alt-allele(s)
// 	    {
// 	      x(ni,1) = v(i).meta.get1_double( data->softtag );
// 	    }	  
// 	  else if ( data->use_postprobs ) // post-probs (assume biallelic)
// 	    {
// 	      std::vector<double> pp = v(i).meta.get_double( data->softtag ); 
// 	      x(ni,1) = pp[1] + 2 * pp[2]; 		
// 	    }
// 	  else // use Genotype::genotype_model to score()
// 	    {
// 	      //x(ni,1) = v(i).minor_allele_count( true );
// 	      x(ni,1) = v(i).score();
// 	    }

// 	  ++ni;
// 	}
//     }


//   // Perform test

//   GLM glm( GLM::LINEAR );
  
//   glm.set( y , x ); 

//   glm.fit(); // should use fast univariate fit() 
  
//   bool valid = glm.valid();
  
//   // Aux. output
  
//   // calc genotype freqs taking phenotype status into account
  
//   // minor allele frequency
//   double maf = 0;

//   // genotype frequencies (reference-based)
//   double f_ref = 0, f_het = 0 , f_hom = 0;
  
//   int tot_allele = 0, tot_genotype = 0;
  
//   for (int i=0; i<n; i++)
//     {

//       if ( ! mask[i] ) 
// 	{

// 	  int c = v(i).copy_number(); 

// 	  int a = v(i).minor_allele_count( true );

// 	  maf          += a;
// 	  tot_allele   += c;
// 	  tot_genotype += 1;
	  
// 	  if ( a == 0 ) ++f_ref;
// 	  else if ( a == 1 ) ++f_het;
// 	  else ++f_hom;
	  
// 	}
      
//     }
  
//   maf /= tot_allele ? (double)tot_allele : 1;
//   f_ref /= tot_genotype ? (double)tot_genotype : 1 ;
//   f_het /= tot_genotype ? (double)tot_genotype : 1 ;
//   f_hom = 1 - f_het - f_ref;

//   double pval = valid ? glm.test_pval() : -1 ;

//   if ( pval < 0 ) valid = false;
  
//   // Output results

//   if ( valid )
//     {
      
//       double se = glm.test_se();
//       double statistic = glm.test_statistic();
//       double coef = glm.test_coef();

//       // For now, let's ignore long-format potential for output. 
     
//       plog.data_group( v );

//       plog.data( v.reference() );
      
//       plog.data( v.alternate() , 1 );  
            
//       plog.data( an , 1 );

//       plog.data( maf , 1 );
          
//       if ( data->dichot_pheno )
// 	{
// 	  if ( ! data->has_covar ) 
// 	    {
// 	      plog.data( mafa );
// 	      plog.data( mafu );
// 	    }
// 	  plog.data( coef );
// 	}
//       else
// 	{
// 	  plog.data( coef );      
// 	}
      
//       plog.data( se );  
//       plog.data( statistic );
//       plog.data( pval );
      
//       plog.print_data_group();
      
//     }
//   else
//     {
     
//       plog.data_group( v );

//       plog.data( v.reference() );
//       plog.data( v.alternate() , 1 );  

//       plog.data( an , 1 );      

//       plog.data( maf , 1 );    

//       if ( data->dichot_pheno )
// 	{
// 	  if ( data->has_covar ) 
// 	    plog.data( maf );
// 	  else
// 	    {
// 	      plog.data( mafa );
// 	      plog.data( mafu );
// 	    }
// 	  plog.data( "NA" );
// 	}
//       else
// 	{
// 	  plog.data( maf );
// 	  plog.data( "NA" );      
// 	}
      
//       plog.data( "NA" );  
//       plog.data( "NA" );
//       plog.data( "NA" );
      
//       plog.print_data_group();
      
//     }




//   // ---------------------------------------------------------------------------------------





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
  
//   std::vector<double>  


//   int obs_a = 0 , obs_u = 0 , obs_tota = 0, obs_totu = 0;

//   const int n = v.size();
  
//   int tota = 0, totu = 0;
  
//   std::set<int> het_carriers;
//   std::set<int> hom_carriers;
//   std::set<int> missing;

//   int obs_refa, obs_heta , obs_homa;
//   int obs_refu, obs_hetu , obs_homu;

//   double obs_maf;
//   double obs_hwe;

//   double obs_statistic;
//   double obs_statistic_dom;
//   double obs_statistic_rec;
      
//   double obs_odds;
//   double obs_odds_dom;
//   double obs_odds_rec;

//   double fisher_pv0, fisher_pv1, fisher_pv2;

//   for (int p = 0; p < R ; p++ )
//     {

//       std::vector<double> statistics(3,0);
      
//       double & statistic     = statistics[0];
//       double & statistic_dom = statistics[1];
//       double & statistic_rec = statistics[2];
      
//       // Output original data, calculate OR, track alt.allele carriers
      
//       if ( p == 0 )
// 	{

// 	  // Genotype counts
// 	  int refa = 0 , refu = 0;
// 	  int heta = 0 , hetu = 0;
// 	  int homa = 0 , homu = 0;

// 	  for (int i=0; i<n; i++)
// 	    {
	      
// 	      affType aff = v.ind( g->perm.pos( i ) )->affected();
	      
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
  
//   plog.data_group( v );

//   plog.data( v.reference() );
//   plog.data( v.print_samples() );
//   plog.data( v.print_meta_filter() );
//   plog.data( v.meta );
//   plog.data( v.consensus.meta );
  
//   plog.data( v.alternate() , 1 );
//   plog.data( obs_maf , 1 );
//   plog.data( obs_hwe , 1 );

//   plog.data( obs_a , 1 );
//   plog.data( obs_u , 1 );

//   plog.data( obs_tota , 1 );
//   plog.data( obs_totu , 1 );
  
//   plog.data( obs_refa , 1 );
//   plog.data( obs_heta , 1 );
//   plog.data( obs_homa , 1 );

//   plog.data( obs_refu , 1 );
//   plog.data( obs_hetu , 1 );
//   plog.data( obs_homu , 1 );

//   plog.data( R == 1 ? fisher_pv0 : g->perm.pvalue(0) , 1 );
//   plog.data( obs_odds , 1 );  
//   if ( R != 1 ) plog.data( g->perm.min_pvalue(0) , 1 );
  
//   plog.data( R == 1 ? fisher_pv1 : g->perm.pvalue(1) , 1 );
//   plog.data( obs_odds_dom , 1 );
//   if ( R != 1 ) plog.data( g->perm.min_pvalue(1) , 1 );

//   plog.data( R == 1 ? fisher_pv2 : g->perm.pvalue(2) , 1 );
//   plog.data( obs_odds_rec , 1 );
//   if ( R != 1 ) plog.data( g->perm.min_pvalue(2) , 1 );

//   plog.print_data_group();

}







bool Pseq::Assoc::set_assoc_test( Mask & m , const Pseq::Util::Options & args )
{

  //
  // Helper class
  //
  
  Pseq::Assoc::AuxGenic a;

  a.g     = &g;
  a.rseed = time(0);

  int nrep = args.has("perm") ? args.as_int( "perm" ) : -1 ;
    
  //
  // Implement the following gene-based tests
  //

  if ( args.has( "midpoint" ) ) a.show_midbp = true;
  
  plog.data_reset();
  
  plog.data_group_header( "LOCUS" );

  plog.data_header( "POS" );

  if ( a.show_midbp ) 
    {
      plog.data_header( "MID" );
      plog.data_header( "BP" );
    }

  plog.data_header( "ALIAS" );

  plog.data_header( "NVAR" );

  plog.data_header( "TEST" );

  plog.data_header( "P" );
  plog.data_header( "I" );
  plog.data_header( "DESC" );
  
  plog.data_header_done();
  

  //
  // Which tests to apply?
  //
  
  a.vanilla       =   args.has( "tests" , "sumstat" );
  a.burden        = ! args.has( "tests" , "no-burden" );
  a.uniq          =   args.has( "tests" , "uniq" );
  a.site_burden   =   args.has( "tests" , "site-burden" );
  a.mhit          =   args.has( "tests" , "mhit" );
  a.vt            =   args.has( "tests" , "vt" );
  a.fw            =   args.has( "tests" , "fw" );
  a.calpha        =   args.has( "tests" , "calpha" );
  a.cancor        =   args.has( "tests" , "cancor" );
  a.hoffman_witte =   args.has( "tests" , "stepup" );
  a.kbac          =   args.has( "tests" , "kbac" );


  const int ntests = a.n_tests();
  
  if ( ntests == 0 ) Helper::halt( "no assoc tests specified" );

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

  //
  // Set up permutation class
  //
   
  g.perm.initiate( nrep , ntests );
  a.fix_null_genotypes = args.has("fix-null");
 
  
  //
  // Apply tests to dataset
  //
  
  g.vardb.iterate( g_set_association , &a , m );


  //
  // Post-processing to obtain corrected p-values
  //

  if ( false ) // skip for now
    {
      for (int t=0; t < g.perm.n_tests(); t++)
	for (int s=0; s < g.perm.n_stats(); s++) 
	  {      
	    plog << "_PCORR\t" 
		 << s << "\t"
		 << t << "\t"
		 << g.perm.max_pvalue(s,t) << "\n";      
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
  // Get auxiliary daya
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
  Pseq::Assoc::Aux_fw_vt aux_fw_vt(data->fw,data->vt);
  Pseq::Assoc::Aux_calpha aux_calpha;
  Pseq::Assoc::Aux_cancor aux_cancor( 1 , vars.size() , vars.n_individuals() ) ;
  Pseq::Assoc::prelim( vars , &aux_prelim );
  
  // External tests
  Pseq::Assoc::Aux_hoffman_witte aux_hoffman_witte( vars, &aux_prelim );
  Pseq::Assoc::Aux_kbac aux_kbac;


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
	}
      
      if ( data->uniq ) 
	{ 
	  test_name.push_back("UNIQ");
	  test_statistic.push_back( aux_burden.stat_uniq );	  
	}
      
      if ( data->mhit )
	{ 
	  test_name.push_back("MHIT");
	  test_statistic.push_back( aux_burden.stat_mhit );
	}
      
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
      test_name.push_back("CANCOR");
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


      //
      // Store all statistics; permute phenotype labels
      //      

      if ( ! g->perm.score( test_statistic ) ) break;
      
      
    } // Next permutation
  
  
  //
  // Output and return for next gene
  //

  plog.data_group( vars.name() );
  plog.data( vars.coordinate() );
  plog.data( g->locdb.alias( vars.name() , false ) );
    //plog.data( "." );
  plog.data( vars.size() );

  if ( data->show_midbp)
    {
      plog.data( vars.midposition() );
      plog.data( vars.span() );
    }

  for (int t=0; t<test_name.size(); t++ )
    {      
      std::string output = test_text[ test_name[t] ];      
      plog.data( test_name[t] , "TEST" , test_name[t] );      
      plog.data( g->perm.pvalue(t) , "P", test_name[t] );
      plog.data( g->perm.min_pvalue(t) , "I" , test_name[t] ); 
      plog.data( output == "" ? "." : output , "DESC" , test_name[t] );
    }
  
  plog.print_data_group();

  return;

}





//void f_set_

// // Per-individual set-enrichment scan

// bool Pseq::Assoc::set_enrich_wrapper( Mask & mask , const Pseq::Util::Options & args )
// {
  

//   if ( args.has( "make-gene-map" ) )
//     {
      
//       std::string filename = args.as_string( "make-gene-map" );
      
//       std::map<std::string,int> counts;
      
//       // iterate through each set; write number of individuals/genes that have 1+ allele
      
//       g.vardb.iterate( g_set_enrich_make_count_map , &counts , );

//       return true;
//     }
  
// //   if ( args.has( "make-set-count-map" ) ) 
// //     {

// //       return true;
// //     }


//   //
//   // otherwise assume that we are going to score individuals
//   //

//   if ( ! args.has( "site-map" ) ) Pseq::halt( "no --site-map {file} specified" );
//   if ( ! args.has( "set-count" ) ) Pseq::halt( "no --set-count {file} specified" );
  
//   std::string sitefile = args.as_string( "site-map" );
//   std::string setfile = args.as_string( "set-count" );
  
//   Helper::checkFileExists( sitefile );
//   Helper::checkFileExists( setfile );
  
//   InFile sites( sitefile );
//   InFile sets( setfile );

//   std::map<std::string,int> gene_counts;
  
//   while ( ! sites.eof() ) 
//     {
//       std::string gene;
//       int count;
//       sites >> gene >> count;
//       gene_counts[ gene ] = count;
//     }
//   sites.close();
  
// //   while ( ! sets.eof() )
// //     {
// //       std::sintr 
// //     }
// //   sets.close();
  

//   // Step 0 : create site map

//   // Step 1 : create set/count map
  
  
//   // total number of genes/individuals per set with 1+ rare-variant per individual
  
//   std::map<std::string,int> stot;
//   int t

//   // For a given individual we have N variants of interest, e.g. LoF variants
//   // Of those, n fall in genes of interest, m do not (N=n+m)

//   // Reframe as, in a set with G genes, G1 have at least one variant of interest, G0 have no
//   // 
//   //  // |   |  |  |  |  |  |  |  |  |  |  |  | 
//   //     20  
//   //     

//   //  // Step 1 -- for each set, create a table of Set-Variant-Count, where Count is collapsed entry per gene

//   //  // Step 2 -- go back, for each set, for each individual, permute, count up # of genes hit 

//   //     For a set, the overall ratio , calcuate T / 

//   //   Now we will have N variants, n1 in sets, n0 not; 

//   //   PERM-i : pick at random N variants from the entire list;  calculate all SET stats;  
  

//   // Create genic-scoring matrix; the total 

//   // Sets will be specified by a locset.group Mask, we assume 
//   // iterate for each set;   
  
//   return true;
// }


