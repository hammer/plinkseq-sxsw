#include "assoc.h"
#include "pseq.h"
#include "func.h"
#include "genic.h"

#include <iostream>
#include <vector>
#include <string>
#include <sstream>

extern GStore g;

extern Pseq::Util::Options args;

void f_glm_association( Variant & v , void * p );


bool Pseq::Assoc::glm_assoc_test( Mask & m , 
				  Pseq::Assoc::Aux_glm & aux )
  
{

  plog.data_reset();

  plog.data_group_header( "VAR" );
  
  // if will be running many tests, we need to be explicit...
  
  plog.data_header( "REF" );
  plog.data_header( "ALT" );
  
  plog.data_header( "N" );

  plog.data_header( "F" );
  

  if ( aux.test_list )
    {
      plog.data_header("SET");  
      plog.data_header("PHE");      
    }

  if ( aux.test_list || aux.show_all_covar ) 
    plog.data_header( "TERM" );  

  if ( aux.test_list )
    {
      plog.data_header( "BETA" );
    }
  else if ( aux.dichot_pheno )
    {
      if ( ! aux.has_covar ) 
	{
	  plog.data_header( "F_A" );
	  plog.data_header( "F_U" );
	}
      plog.data_header( "OR" );
    }
  else
    {
      plog.data_header( "BETA" );      
    }

  plog.data_header( "SE" );
  plog.data_header( "STAT" );
  plog.data_header( "P" );
  
  plog.data_header_done();
  
  //
  // Either just iterate through each SNP in the mask, given a single phenotype and covariate set
  // OR, read a list of 'tests' from a file
  //
  
  if ( aux.test_list )
    {
      return Pseq::Assoc::glm_assoc_testlist( m , aux );
    }

  //
  // Set phenotype, and any covariates
  //

  const int n = g.indmap.size();

  aux.y.resize( n );
  aux.mask.resize( n , false );
  
  for (int i=0; i<g.indmap.size(); i++)
    {      
      Individual * person = g.indmap.ind(i);
      
      if ( person->missing() )
	{ 
	  aux.mask[i] = true; //mask out
	}
      else
	{
	  if ( aux.dichot_pheno ) 
	    aux.y[i] = person->affected() == CASE ? 1 : 0;
	  else
	    aux.y[i] = person->qt();
	}            
    }

  
  // Covariates, that come with their own mask, handled below  
  if ( aux.has_covar )
    aux.c = g.phmap.covariates( aux.covars , g.indmap , &g.vardb );
      
  // Run GLM tests  
  g.vardb.iterate( f_glm_association , &aux , m );
  
  // Done
  return true;
  
}




void f_glm_association( Variant & v , void * p )
{
  
  Pseq::Assoc::Aux_glm * data = (Pseq::Assoc::Aux_glm *)p;
  
  // Extract out non-missing phenotypes
  // (in future, set up a MASK on the MATRIX)
  
  // For a given variant, get non-missing n and 
  // create phenotypes
  
  const int n = v.size();
  const int np = 1;  // num of (genetic) parameters
  const int nc = data->c.dim2(); // num of covariates
  int an = n;

  std::vector<bool> mask( n , false ); // F means include
  
  for (int i=0; i<n; i++)
    {
      if ( data->mask[i] )     
	{ 
	  --an; 
	  mask[i] = true; 
	}
      else if ( data->c.masked(i) ) 
	{ 
	  --an; 
	  mask[i] = true; 
	}
      else if ( data->use_dosage && ! v(i).meta.has_field( data->softtag ) ) 
	{ 
	  --an; 
	  mask[i] = true; 
	}
      else if ( data->use_postprobs ) 
	{
	  if ( ( ! v(i).meta.has_field( data->softtag ) ) 
	       || v(i).meta.get_double( data->softtag ).size() != 3 ) 
	    { 
	      --an; mask[i] = true; 
	    }
	}
      else if ( v(i).null() && ! ( data->use_postprobs || data->use_dosage ) ) 
	{ 
	  --an; 
	  mask[i] = true; 
	}
    }

  if ( an == 0 ) Helper::halt("no valid observations");

  Data::Vector<double> y( an );
  Data::Matrix<double> x( an , 1 + np + nc ) ;

  // Populate
  int ni = 0;
  for (int i=0; i<n; i++)
    {

      if ( ! mask[i] ) 
	{
	  // DV
	  y[ni] = data->y[i];
	  
	  // Intercept
	  x(ni,0) = 1;
	  
	  // Genotype
	  if ( data->use_dosage ) // dosage of alt-allele(s)
	    {	      
	      x(ni,1) = v(i).meta.get1_double( data->softtag );
	    }	  
	  else if ( data->use_postprobs ) // post-probs (assume biallelic)
	    {
	      std::vector<double> pp = v(i).meta.get_double( data->softtag ); 
	      x(ni,1) = pp[1] + 2 * pp[2]; 		
	    }
	  else // use Genotype::genotype_model to score()
	    {
	      //x(ni,1) = v(i).minor_allele_count( true );
	      x(ni,1) = v(i).score();
	    }

	  // Covariates (1,variant,covar)
	  int z = 1; // skip first two slots
	  for (int j=0;j<nc;j++) x(ni,++z) = data->c(i,j);
	  
	  ++ni;
	}
    }

  // Perform test
  
  GLM glm( data->dichot_pheno ? GLM::LOGISTIC : GLM::LINEAR );

//   for (int i=0;i<an;i++)
//     {
//       std::cout << "data = " << y[i] << "\t";
//       for (int j=0;j<x.dim2();j++) 
// 	std::cout << "\t" << x(i,j);
//       std::cout << "\n";
//     }

  glm.set( y , x ); 

  glm.fit();

  bool valid = glm.valid();
  
  // Aux. output
  
  // calc freqs taking phenotype status into account
  double maf = 0;
  double mafa = 0 , mafu = 0;
  int tot = 0, tota = 0 , totu = 0;
  for (int i=0; i<n; i++)
    {

      if ( ! mask[i] ) 
	{

	  int c = v(i).copy_number(); 
	  int a = v(i).minor_allele_count( true );
	  
	  maf += a;
	  tot += c;
	  if ( data->dichot_pheno ) 
	    {
	      if ( v.ind(i)->affected() == CASE ) 
		{
		  tota += c;
		  mafa += a;
		}
	      else 
		{
		  totu += c;
		  mafu += a;
		}
	    }
	}
    }
  
  maf /= tot ? (double)tot : 1;
  mafa /= tota ? (double)tota : 1;
  mafu /= totu ? (double)totu : 1;
  
  double pval = 0;
  if ( valid ) pval = glm.test_pval();
  if ( pval < 0 ) valid = false;
  
  // Output results

  if ( valid )
    {
      
      if ( data->show_all_covar ) 
	{

	  std::vector<bool> mask;
	  Data::Vector<double> beta;
	  Data::Vector<double> se;
	  Data::Vector<double> lowci;
	  Data::Vector<double> uprci;
	  Data::Vector<double> statistic;
	  Data::Vector<double> pvalue;
	
	  glm.display( &beta, &se, &pvalue , &mask, &lowci, &uprci, &statistic );

	  for (int term = 0 ; term < pvalue.size() ; term++ ) 
	    {
	      if ( term == 0 && ! data->show_intercept ) continue;

	      plog.data_group( v );	      
	      plog.data( v.reference() );
	      plog.data( v.alternate() , 1 );              
	      plog.data( an , 1 );

	      if ( term != 1 || data->use_dosage || data->use_postprobs ) 
		plog.data( "." , 1 );
	      else
		plog.data( maf , 1 );
	      
	      if ( data->test_list ) 
		{
		  plog.data( data->test_list );
		  plog.data( g.phmap.phenotype() );
		}

	      if ( term == 0 ) 
		plog.data( "b0" );
	      else if ( term == 1 )
		plog.data( v.name() );
	      else 
		plog.data( data->covars[ term - 2 ] ) ;
	      
	      plog.data( beta[ term ] );      
	      
	      plog.data( se[ term ] );  
	      plog.data( statistic[ term ] );
	      plog.data( pvalue[ term ] );	      
	      plog.print_data_group();	      
	    }
	  
	}
      else
	{

	  double se = glm.test_se();
	  double statistic = glm.test_statistic();
	  double coef = glm.test_coef();
	  
	  // For now, let's ignore long-format potential for output. 
	  //  So everything has 'group' == 1 
	  
	  plog.data_group( v );
	  
	  plog.data( v.reference() );
	  plog.data( v.alternate() , 1 );              
	  plog.data( an , 1 );
	  
	  if ( data->use_dosage || data->use_postprobs ) 
	    plog.data( "." , 1 );
	  else
	    plog.data( maf , 1 );
	  
	  if ( data->show_all_covar )
	    plog.data( v.name() );
	  
	  if ( data->dichot_pheno )
	    {
	      if ( ! data->has_covar ) 
		{
		  plog.data( mafa );
		  plog.data( mafu );
		}
	      plog.data( coef );
	    }
	  else
	    {
	      plog.data( coef );      
	    }
	  
	  plog.data( se );  
	  plog.data( statistic );
	  plog.data( pval );
	  
	  plog.print_data_group();
	}
    }
  else
    {
     
      plog.data_group( v );
      plog.data( v.reference() );
      plog.data( v.alternate() , 1 );  

      plog.data( an , 1 );      

      plog.data( maf , 1 );    

      if ( data->show_all_covar )
	plog.data( v.name() );
      
      if ( data->dichot_pheno )
	{
	  if ( data->has_covar ) 
	    plog.data( maf );
	  else
	    {
	      plog.data( mafa );
	      plog.data( mafu );
	    }
	  plog.data( "NA" );
	}
      else
	{
	  plog.data( maf );
	  plog.data( "NA" );      
	}
      
      plog.data( "NA" );  
      plog.data( "NA" );
      plog.data( "NA" );
      
      plog.print_data_group();
      
    }

}



bool Pseq::Assoc::glm_assoc_testlist( Mask & m , 
				      Pseq::Assoc::Aux_glm & aux )
{

  Helper::checkFileExists( aux.test_list_file );
  
  InFile f( aux.test_list_file );
  
  const int n = g.indmap.size();
  
  std::vector<std::string> base_covars = aux.covars;
  
  while ( ! f.eof() )
    {

      // expecting format : 4 tab-delimited columns

      //  Y  X C1 C2 C2
      //  If -C1 -C2  means residualize Y first 
      
      
      // If X == ".", then iterate over all SNPs in the model with this phenotype/covariate set-up
      
      std::vector<std::string> h = f.tokenizeLine("\t ");
      if ( h.size() == 0 ) continue;
      if ( h[0][0] == '#' ) continue; // skip comments

      if ( h.size() < 2 ) 
	{
	  plog.warn( "row of specified tests has less than two fields" );
	  continue;
	}


      // Phenotype (can use default phenotype if entry is '.')

      std::string phenotype = h[0] != "." ? h[0] : g.phmap.phenotype();
      
      if ( ! g.phmap.set_phenotype ( phenotype ) ) 
	Helper::halt( "problem setting phenotype " + phenotype );
      
      // Variant (will look at all variants in mask, if '.')
      std::string variant = h[1];
      
      // Covariates
      // 1) if any --covar specified on command line always include these in any model
      // 2) if any "-COVAR" listed, means to residualize phenotype first

      aux.covars = base_covars; 

      std::vector<std::string> resid;
      for (int j=2;j<h.size();j++)
	{	  
	  if ( h[j][0] == '-' ) resid.push_back( h[j].substr(1) );
	  else aux.covars.push_back( h[j] );
	}
      
      aux.has_covar = aux.covars.size() > 0 ;
      
      
      //
      // Do we need to calculate a new, residualized phenotype?
      //

      if ( resid.size() > 0 ) 
	Pseq::IndDB::make_residuals( resid );
      
      

      //
      // Create Y and C
      //
      
      aux.dichot_pheno = g.phmap.type() == PHE_DICHOT;
      
      aux.mask.clear();
      aux.mask.resize( n , false );
      aux.y.resize( n );
      
      for (int i=0; i<n; i++)
	{      	  
	  Individual * person = g.indmap.ind(i);	  
	  if ( person->missing() )
	    aux.mask[i] = true; //mask out
	  else
	    {
	      if ( aux.dichot_pheno ) 
		aux.y[i] = person->affected() == CASE ? 1 : 0;
	      else
		aux.y[i] = person->qt();
	    }            
	}
      

      // Covariates?
      
      if ( aux.has_covar )
	aux.c = g.phmap.covariates( aux.covars , g.indmap , &g.vardb );
      else
	aux.c.clear();
      
      // Either run *all* GLM tests, of for a specific variant

      if ( variant == "." ) 
	{
	  g.vardb.iterate( f_glm_association , &aux , m );
	}
      else
	{
	  bool okay = false;
	  Region region( variant , okay );
	  if ( ! okay ) continue;
	  std::set<Variant> vars = g.vardb.fetch( region );
	  std::set<Variant>::iterator v = vars.begin();
	  while ( v != vars.end() )
	    {
	      f_glm_association( (Variant&)*v , &aux );
	      ++v;
	    }
	}

      // set a new unique test-number (to help organise output)
      ++aux.test_list;

    } // next row of input file

  // tidy up
  f.close();

  return true;
}



//
// Helper function to make a new phenotype that is based on residuals, given
// an existing phenotype in the phmap
//

bool Pseq::IndDB::make_residuals( const std::vector<std::string> & covars )
{
  // set up a GLM to regress standard phenotype on X1, X2, X3, etc
  // then create a new phenotype, and set as standard phenotype, 
  // that is the residuals from this model
  
  const int n = g.indmap.size();
  
  int an = 0;  // actual number (after case-wise deletion)
  std::vector<bool> mask( n , false ) ; 
  
  pType ptype = g.phmap.type(); 


  // Get covariates, that come with their own mask, handled below  
  Data::Matrix<double> C = g.phmap.covariates( covars , g.indmap , &g.vardb );
  
  for (int i=0; i<g.indmap.size(); i++)
    {
      if ( g.indmap(i)->missing() ) mask[i] = true; //mask out
      if ( C.masked(i) ) mask[i] = true; 
      if (  ! mask[i] ) ++an;
    }
  
  if ( an == 0 ) Helper::halt("no valid observations");

  const int ncov = C.dim2();
  
  Data::Vector<double> Y(an);
  Data::Matrix<double> X(an, 1 + ncov ); // intercept + covariates
  
  int ni = 0;
  for (int i=0; i<n; i++)
    {      
      
      if ( ! mask[i] ) 
	{
	  // DV
	  Y[ni] = ptype == PHE_DICHOT ? g.indmap(i)->affected() == CASE : g.indmap(i)->qt() ;
	  
	  // Intercept
	  X(ni,0) = 1;
	  
	  // Covariates
	  int z = 0; // i.e. skip first slot (intercept)
	  for (int j=0;j<ncov;j++) 
	      X(ni,++z) = C(i,j);

	  ++ni;
	}
    }

  // Perform test
  
  GLM glm( ptype == PHE_DICHOT ? GLM::LOGISTIC : GLM::LINEAR );

  glm.set( Y , X ); 
  
  glm.fit();
  
  bool valid = glm.valid();
  
  if ( ! valid ) 
    Helper::halt( "unable to fit model to calculate residuals" );
  
  Data::Vector<double> beta;
  if ( ! glm.display(&beta) ) 
    Helper::halt( "unable to fit model to calculate residuals" );

  // by default, glm() returns OR scale; so return...
  if ( ptype == PHE_DICHOT )    
    for (int j=0;j<beta.size();j++) beta[j] = log(beta[j]);
  
  // calculate expected score for each individual
  // and set residual as new phenotype

  ni = 0;
  
  std::vector<double> resid( n , 0 );
  std::vector<double> missing( n , false );

  for (int i=0;i<n;i++)
    {

      Individual * person = g.indmap(i);
      
      if ( ! mask[i] ) 
	{
	  double yr = beta[0];
	  
	  for (int c=0;c<ncov;c++) 
	    yr += X(ni,c+1) * beta[c+1];
	    
	  // Y=logit^-1(Y)
	  if ( ptype == PHE_DICHOT ) 
	    {	      
	      double e = exp( yr );
	      yr = e / ( 1 + e ) ;	      
	    }

	  // set
	  resid[i] = Y[ni] - yr;	  
	  
	  ++ni;

	}
    }

  std::string pheno_name = "Resid(" + g.phmap.phenotype() + "~" ;
  for (int c=0;c<ncov;c++) 
    {
      if ( c ) pheno_name += "+";
      pheno_name += covars[c] ;
    }
  pheno_name += ")";

  if ( ! g.phmap.attach_qt_phenotype( pheno_name  , mask , resid , g.indmap ) ) 
    Helper::halt( "no individuals with valud residual phenotype set" );

}


