#include "assoc.h"
#include "pseq.h"
#include "func.h"
#include "genic.h"

#include <iostream>
#include <vector>
#include <string>
#include <sstream>

extern GStore g;

extern Pseq::Util::Options options;

void f_glm_association( Variant & v , void * p );

bool Pseq::Assoc::glm_assoc_test( Mask & m , 
				  Pseq::Assoc::Aux_glm & aux )
  
{

  plog.data_reset();
  
  plog.data_group_header( "VAR" );
  
  plog.data_header( "REF" );
  plog.data_header( "ALT" );

  plog.data_header( "N" );

  plog.data_header( "MAF" );

  if ( aux.show_all_covar ) 
    plog.data_header( "TEST" );
  
  if ( aux.dichot_pheno )
    {
      if ( aux.has_covar ) 
	plog.data_header( "F" );
      else
	{
	  plog.data_header( "F_A" );
	  plog.data_header( "F_U" );
	}
      plog.data_header( "OR" );
    }
  else
    {
      plog.data_header( "F" );
      plog.data_header( "BETA" );      
    }

  plog.data_header( "SE" );  
  plog.data_header( "STAT" );
  plog.data_header( "P" );
  
  plog.data_header_done();
  
   
  // Set phenotype, and any covariates

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
      
      // Covariates,  comes with it's own mask handled below
      
      if ( aux.has_covar )
	{
	  aux.c = g.phmap.covariates( aux.covars );
	}
      
    }

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
      if ( v(i).null() ) { --an; mask[i] = true; }
      else if ( data->mask[i] ) { --an; mask[i] = true; }
      else if ( data->c.masked(i) ) { --an; mask[i] = true; }
    }
  
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
	  
	  // Hard-coding
	  x(ni,1) = v(i).minor_allele_count( true );

	  // Covariates
	  int z = 1;
	  for (int j=0;j<nc;j++) x(ni,++z) = data->c(i,j);

	  ++ni;
	}
    }

  // Perform test

  GLM glm( data->dichot_pheno ? GLM::LOGISTIC : GLM::LINEAR );

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

  double pval;
  if ( valid ) pval = glm.test_pval();
  if ( pval < 0 ) valid = false;

  // Output results

  if ( valid )
    {
      
      double se = glm.test_se();
      double statistic = glm.test_statistic();
      double coef = glm.test_coef();

      // For now, let's ignore long-format potential for output. 
     
      plog.data_group( v );
      plog.data( v.reference() );
      
      plog.data( v.alternate() , 1 );  
      
      
      plog.data( an , 1 );

      plog.data( maf , 1 );
    
      if ( data->show_all_covar )
	plog.data( "VAR" );
      
      if ( data->dichot_pheno )
	{
	  if ( data->has_covar ) 
	    plog.data( maf );
	  else
	    {
	      plog.data( mafa );
	      plog.data( mafu );
	    }
	  plog.data( coef );
	}
      else
	{
	  plog.data( maf );
	  plog.data( coef );      
	}
      
      plog.data( se );  
      plog.data( statistic );
      plog.data( pval );
      
      plog.print_data_group();
      
    }
  else
    {
     
      plog.data_group( v );
      plog.data( v.reference() );
      plog.data( v.alternate() , 1 );  

      plog.data( an , 1 );      

      plog.data( maf , 1 );    

      if ( data->show_all_covar )
	plog.data( "VAR" );
      
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
