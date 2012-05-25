#include "assoc.h"
#include "pseq.h"
#include "util.h"
#include "genic.h"

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
  
  //
  // Standard output mode, or dumping a matrix of null-statistics
  // (with originals, and labels, as first row)
  //

  a.dump_stats_matrix = args.has( "dump-null-matrix" );


  //
  // Write header, if in stanard OUTPUT mode
  //

  if ( ! a.dump_stats_matrix )
    {
  
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
    }


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
  a.two_hit       =   args.has( "tests" , "two_hit");
  a.skat          =   args.has( "tests" , "skat" );

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
  //          two-hit  TWO-HIT
  //          skat     SKAT



  //
  // Write list of tests, if in matrix-dump mode
  //

  if ( a.dump_stats_matrix ) 
    {
      plog << nrep ;
      if ( a.burden ) plog << "\t" << "BURDEN"; 
      if ( a.uniq )  plog << "\t" << "UNIQ"; 
      if ( a.skat )  plog << "\t" << "SKAT"; 
      plog << "\n";
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
      if ( ntests > 1 || ! a.skat ) 
	Helper::halt( "only SKAT can handle quantitative traits" );      
    }
  else      
    Helper::halt("no dichotomous phenotype specified");


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
		plog << "_PCORR\t" 
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

  Pseq::Assoc::prelim( vars , &aux_prelim );

  
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

  if ( data->dump_stats_matrix ) 
    {
      plog << vars.name() ;
    }



  //
  // Apply tests to original dataset
  //
  
  std::vector<double> test_statistic;
  std::vector<std::string> test_name;
  std::map<std::string,std::string> test_text;
  

  double prev = .006;
  if ( args.has( "prev" ) )
    prev =  Helper::str2dbl(args.as_string( "prev" ));

  bool mhit = args.has( "mhit" );


  std::map< std::string, int > func_inc;
  std::map< std::string, int > func_exc;

  if( args.has( "func-inc" ) ){
    std::vector<std::string> inc = args.as_string_vector( "func-inc" );
    for( int i = 0; i < inc.size(); i++ )
      func_inc[inc[i]] = i;
  }

  if( args.has( "func-exc" ) ){
    std::vector<std::string> inc = args.as_string_vector( "func-exc" );
    for( int i = 0; i < inc.size(); i++ )
      func_exc[inc[i]] = i;
  }



  




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
	  if ( data->dump_stats_matrix ) plog << "\t" << aux_burden.stat_burden ;
	}
      
      if ( data->uniq ) 
	{ 
	  test_name.push_back("UNIQ");
	  test_statistic.push_back( aux_burden.stat_uniq );	  
	  if ( data->dump_stats_matrix ) plog << "\t" << aux_burden.stat_uniq ;
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


  if ( data->two_hit )
    {
      test_name.push_back( "TWO-HIT" );
      
      double statistic = Pseq::Assoc::stat_two_hit( vars , &aux_prelim , &aux_two_hit , &test_text , true , func_inc, func_exc , prev, mhit );
      test_statistic.push_back( statistic );
    }
  
  
  if ( data->skat ) 
    {
      test_name.push_back( "SKAT" );
      double statistic = Pseq::Assoc::stat_skat( vars , &aux_prelim , &aux_skat , &test_text , true ); 
      test_statistic.push_back( statistic );
      if ( data->dump_stats_matrix ) plog << "\t" << statistic;
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
	      if ( data->burden ) plog << "\t" << aux_burden.stat_burden;
	      if ( data->uniq ) plog << "\t" << aux_burden.stat_uniq;
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
	  double statistic = Pseq::Assoc::stat_two_hit( vars , &aux_prelim , &aux_two_hit , NULL , false , func_inc, func_exc , prev, mhit );
	  test_statistic.push_back( statistic );
	}
      if ( data->skat )
	{	  
	  double statistic = Pseq::Assoc::stat_skat( vars , &aux_prelim , &aux_skat , NULL , false  );
	  test_statistic.push_back( statistic );
	  if ( data->dump_stats_matrix ) plog << "\t" << statistic;
	}



      //
      // Store all statistics; permute phenotype labels
      //      

      if ( ! g->perm.score( test_statistic ) ) break;
    
      
    } // Next permutation
  
  


  //
  // Output and return for next gene
  //

  if ( ! data->dump_stats_matrix ) 
    {

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
	  if ( R != 1 ) 
	    {
	      plog.data( g->perm.pvalue(t) , "P", test_name[t] );
	      plog.data( g->perm.min_pvalue(t) , "I" , test_name[t] ); 
	    }
	  else // asymptotic p-values
	    { 
	      plog.data( "." , "P", test_name[t] );
	      plog.data( "." , "I" , test_name[t] ); 
	    }
	  plog.data( output == "" ? "." : output , "DESC" , test_name[t] );
	}
      
      plog.print_data_group();
    }
  else // in dump-stats-mode
    {
      plog << "\n";
    }

  return;

}




// Per-individual set-enrichment scan

struct aux_set_enrichment{
  std::vector<int> indiv;
  std::vector<int> gene;
  std::map<std::string,int> * genes_slot;
};


void g_set_enrichment( VariantGroup & vars , void * p )
{

  aux_set_enrichment * aux = (aux_set_enrichment*)p;
  
  // only consider scoring for genes that are in the gene-map we will be using
  if ( aux->genes_slot->find( vars.name() ) == aux->genes_slot->end() ) return;
  
  // get numeric ID 
  const int gn = (*aux->genes_slot)[ vars.name() ];
  
  int a = 0;
  const int n = g.indmap.size();
  const int s = vars.size();
  
  std::vector<bool> altmin(s);
  for (int j = 0 ; j < s ; j++  )
    {
      int c, c_tot;      
      altmin[j] = vars(j).n_minor_allele( &c , &c_tot );      
    }

  for ( int i = 0 ; i < n ; i++ )
    {
      for (int j = 0 ; j < s ; j++  )
	{
	  if ( vars(j,i).minor_allele( altmin[j] ) )
	    {
	      aux->gene.push_back( gn );
	      aux->indiv.push_back( i );
	      ++a;
	      break;
	    }
	}
    }

  //  std::cout << " added " << a << " for " << vars.name() << "\n";

}


bool Pseq::Assoc::set_enrich_wrapper( Mask & mask , const Pseq::Util::Options & args )
{
  
  // We assume the initial scan has to be based on genes --mask loc.group=refseq

  // We need a geneset list to be specified by the argment --locset go
  // and that this relates to a locset that has been loaded into the LOCDB
  // (and that will match in terms of 'refseq')
  
  // Populate a list of individual / gene initially
  
  aux_set_enrichment aux;
  
  std::string locset = args.as_string( "locset" );
  std::string loc = args.as_string( "loc" );
  
  //
  // First, extract set information
  //

  plog << "pulling all set information\n";

  std::map<int,std::string> sets;
  std::map<int,std::string> genes;
  std::map<int,std::set<int> > s2g;
  std::map<int,std::set<int> > g2s;
  
  bool okay = g.locdb.populate_set_structures( loc , locset , &genes , &sets , &s2g , &g2s );

  //
  // Create gene-string --> numeric-ID mapping, and attach to aux
  //
  
  std::map<std::string,int> genes_slot;
  std::map<int,int> genes_slot2map;  
  std::map<int,std::string>::iterator i = genes.begin();

  int slot = 0;
  while ( i != genes.end() )
    {
      genes_slot[ i->second ] = slot;
      genes_slot2map[ slot ] = i->first;
      ++slot;
      ++i;
    }
  
  aux.genes_slot = &genes_slot;
  
  
  
  //
  // ensure that we iterate by gene initially, to build the gene table
  // we want to do this for all genes (whether in set or not)
  //

  mask.group_loc( loc );
  

  //
  // iterate through each set; write number of individuals/genes that have 1+ allele
  //
  
  g.vardb.iterate( g_set_enrichment , &aux , mask ) ;



  //
  // Set ID --> Set Slot mapping
  //

  std::map<int,std::string>::iterator si = sets.begin();
  std::map<int,int> s2s;
  slot = 0;
  while ( si != sets.end() )
    {
      s2s[ si->first ] = slot++;
      ++si;
    }


  //
  // Create large score matrix of indiv x set
  //
  
  const int n = g.indmap.size();  
  std::vector<std::vector<int> > scores(n);  
  for (int k=0;k<n;k++) scores[k].resize( sets.size() , 0 );
  

  //
  // Create per individual counts
  //
  
  std::vector<int> indcnt( n , 0 );
  std::map<int,int> gcnt;
  
  //
  // Populate counts
  //
  
  for (int k = 0 ; k < aux.indiv.size() ; k++)
    {
      const int i = aux.indiv[k];
      const int gslot = aux.gene[k];                // col position in SCORES matrix
      const int gmap  = genes_slot2map[ gslot ];    // ID from LOCDB that connects genes to sets
      
      //      std::cout << "Individual/gene " << k << " of " << aux.indiv.size() << "\t" << i << " " << gslot << " " << gmap << "\n";
      
      // track total number of mutant genes per individual
      indcnt[i]++;
      
      // similarly, for each gene
      gcnt[gmap]++;
      
      // track number of mutants per set per individual
      std::set<int> & gsets = g2s[ gmap ];
      std::set<int>::iterator ii = gsets.begin();
      while ( ii != gsets.end() ) 
	{
	  scores[i][ s2s[ *ii ] ]++;
	  ++ii;
	}
      
    }

  
  //
  // Calculate enrichment statistic per individual for each set
  //


  //           | Set | Not-Set
  // ----------|-----|----------
  // Indiv V   |     |
  // ----------|-----|----------
  // Rest  V-1 |     |
  // ----------|-----|----------


  const int numgenes = g2s.size();
  const int numsets = s2g.size();
  
  for ( int i = 0 ; i < n ; i++ ) 
    plog << "ICNT" << "\t"
	 << g.indmap.ind(i)->id() << "\t"
	 << indcnt[i] << "\n";
  
  // For each gene, number of mutations (per person)
  std::map<int,int> genic_hits;
  std::map<int,std::set<int> >::iterator gi = g2s.begin();
  while ( gi != g2s.end() )
    {
      plog << "GCNT" << "\t"
	   << gi->second.size() << "\t"
	   << gcnt.find( gi->first )->second << "\t"
	   << genes[ gi->first ] << "\n";
	   ++gi;
    }
      

  // For each set, number of genes & number of mutations (per person)   
  si = sets.begin();
  while ( si != sets.end() )
    {

      plog << "SCNT" << "\t"
	   << s2g[ si->first ].size() << "\t";
      
      int scnt = 0;
      std::set<int> & gs = s2g.find( si->first )->second;
      std::set<int>::iterator gi = gs.begin();
      while ( gi != gs.end() )
	{
	  scnt += gcnt[ *gi ];
	  ++gi;
	}
      plog << scnt << "\t"
	   << si->second << "\n";
      ++si;
    }



  // 
  // Test each individual for each set
  //

  ///  std::cout << "n = " << n << "\n";
  
  for ( int i = 0 ; i < n ; i++ ) 
    {
      
      int g0s0 = 0 , g0s1 = 0 , g1s0 = 0 , g1s2 = 0;

      int slot = 0;

      //      std::cout << "set N = " << s2g.size() << "\n";

      std::map<int,std::set<int> >::iterator si = s2g.begin();
      while ( si != s2g.end() )
	{
	  
 	  int g1s1 = scores[i][slot++];
 	  int g1s0 = indcnt[i]  - g1s1;
 	  int g0s1 = si->second.size() - g1s1; 
	  int g0s0 = numgenes - g0s1 - g1s0 - g1s1;
	  
	  plog << g.indmap.ind(i)->id() << "\t";

	  plog << g0s0 << " " << g0s1 << " / " << g1s0 << " " << g1s1 << "\t";
	  
	  // basic Fisher's test
 	  double pval = 0;
 	  if ( ! fisher( Table( g0s0 , g0s1 , g1s0 , g1s1 ) , &pval ) ) pval = -1;
	  if ( pval < 0 ) plog << "NA";
	  else plog << pval;

	  plog << "\t" << sets[ si->first ] << "\n";

	  ++si;

 	}
    }
  


  //
  // Permutation procedure is simply to shuffle the phenotype on the list of (i.e. disconnect aux.indiv[] and aux.gene[]
  //
  
  
  
  //
  //  
  // 
  
  return true;

}
