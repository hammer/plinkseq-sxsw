#include "vardb.h"
#include "vcf.h"
#include "gstore.h"

extern GStore * GP;

void f_vcfiterate_inserter( Variant & var , void * data )
{
  //  std::cout << "inserter... " << var << "\n";
}

IterationReport VarDBase::vcf_iterate( void (*f)(Variant&, void *) , void * data , Mask & mask )
{
  
  IterationReport irep( true , mask.any_grouping() , mask.variant_limit() );

  // VCF file name is kept in the Mask, by the 'ex-vcf' attribute
  
  std::string filename = mask.external_vcf_filename();
  Helper::checkFileExists( filename );

  //
  // Use VCFReader, into a temporary :memory: database
  //

  IndividualMap imap;

  VarDBase tmpdb( imap );
  tmpdb.attach( ":memory:" );

  // Load, parse VCF file; store variant and genotype information, and
  // meta-information, in vardb
  
  File vcffile( filename , VCF );

  VCFReader v( &vcffile , "" , &tmpdb ,  NULL );

  // Selectively filter in/out meta-information?
  // or, add a region filter?

  // We might not want to load the entire VCF into memory; & thus allow
  // the includes/excludes and filters to bring regions/tags into view

//   std::set<Region> filter;
//   std::string locinc = mask.loc_include_string();
//   if ( locinc != "" ) 
//     {
//       filter = GP->locdb.get_regions( locinc );
//       v.set_region_mask( &filter );
//     }
  

//   // filters on meta-fields?
//   std::set<std::string> includes, excludes;
//   if ( options.key( "meta" ) ) 
//     {
//       includes = options.get_set( "meta" );
//       if ( includes.size() > 0 ) v.get_meta( includes );
//     }
  
//   if ( options.key( "meta.ex" ) ) 
//     {
//       excludes = options.get_set( "meta.ex" );
//       if ( excludes.size() > 0 ) v.ignore_meta( excludes );
//     }




  downcode_mode = mask.downcode();


  // 
  // Work through VCF
  //

  tmpdb.begin();

  int inserted = 0;

  v.return_variant( true );

  Variant * pv = NULL;
  
  while ( 1 ) 
    { 

      VCFReader::line_t l = v.parseLine( &pv );
      
      if ( l == VCFReader::VCF_EOF ) break;
      
      if ( l == VCFReader::VCF_INVALID ) 
	{
	  continue;
	}

      if ( l == VCFReader::VCF_HEADER ) 
	{
	  int n = imap.populate( tmpdb , GP->phmap , mask );
	}
      
      
      // If a variant line has been processed and meets criteria, pv
      // will be non-NULL which also implies that a Variant has been
      // created and we are reponsible for cleaning up afterwards
      
      if ( l == VCFReader::VCF_VARIANT && pv )
	{
	  
	  // So that the Variant functions know not to look for data
	  // in a BLOB	  
	  pv->consensus.vcf_direct = true;
	  
	  // Apply all mask filters, and decide whether to call function
	  if ( eval_and_call( mask, &imap , *pv , f , data ) ) 
	    {
	      if ( ! irep.accepted_variant() ) break;
	    }
	  else
	    {
	      irep.rejected_variant();
	    }
	  
	  // and now clean up 	  
	  delete pv;
	  
	}
      
    }

  // Wrap up
  
  tmpdb.commit();
  
  // If we had to use any positional filters or grouping, now run the 
  // actual iterate on the temporary database we've created
     
  return irep;
}

