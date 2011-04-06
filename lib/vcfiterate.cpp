#include "vardb.h"
#include "vcf.h"
#include "gstore.h"

extern GStore * GP;

bool VarDBase::vcf_iterate_read_header( Mask & mask )
{

  // Read meta-information, and header row.
  // Populate the vardb, indmap, etc
  
  std::string filename = mask.external_vcf_filename();
  Helper::checkFileExists( filename );
  
  //
  // Use VCFReader, into a temporary :memory: database
  //

//   IndividualMap imap;
//   VarDBase tmpdb( imap );
//   tmpdb.attach( ":memory:" );

  // Load, parse VCF file; store variant and genotype information, and
  // meta-information, in vardb
  
  File vcffile( filename , VCF );

  VCFReader v( &vcffile , "" , &(GP->vardb) ,  NULL );

  // 
  // Work through VCF
  //

  GP->vardb.begin();

  while ( 1 ) 
    { 
      
      VCFReader::line_t l = v.parseLine();
      
      if ( l == VCFReader::VCF_EOF ) break;
      
      if ( l == VCFReader::VCF_INVALID ) 
	{
	  continue;
	}
      
      if ( l == VCFReader::VCF_HEADER )
	{
	  int n = GP->indmap.populate( GP->vardb , GP->phmap , mask );
	  break;
	}           
    }

  // Wrap up  
  GP->vardb.commit();

  
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

  // Load, parse VCF file; store variant and genotype information, and
  // meta-information, in vardb
  
  File vcffile( filename , VCF );

  VCFReader v( &vcffile , "" , &(GP->vardb) ,  NULL );

  // Selectively filter in/out meta-information?
  // or, add a region filter?

  // We might not want to load the entire VCF into memory; & thus allow
  // the includes/excludes and filters to bring regions/tags into view

  // Respect 'reg' and 'loc' from command line.
  // But not loc.subset; loc.req, loc.ex, etc
  
  std::set<Region> filter;
  std::string locinc = mask.loc_include_string();
  if ( locinc != "" ) 
    filter = GP->locdb.get_regions( locinc );

  std::set<Region> reginc = mask.included_reg();
  std::set<Region>::iterator ii = reginc.begin();
  while ( ii != reginc.end() ) 
    {
      filter.insert( *ii );
      ++ii;
    }

  // Add other "reg" from mask? 
  if ( filter.size() > 0 ) 
    v.set_region_mask( &filter );  
  

  //
  // Misc. settings.
  //
  
  downcode_mode = mask.downcode();


  // 
  // Work through VCF
  //

//  tmpdb.begin();
  GP->vardb.begin();

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

//       if ( l == VCFReader::VCF_HEADER ) 
// 	{
// //	  int n = imap.populate( tmpdb , GP->phmap , mask );
// 	  int n = GP->indmap.populate( GP->vardb , GP->phmap , mask );
// 	}
      
      
      // If a variant line has been processed and meets criteria, pv
      // will be non-NULL which also implies that a Variant has been
      // created and we are reponsible for cleaning up afterwards
      
      if ( l == VCFReader::VCF_VARIANT && pv )
	{
	  
	  // So that the Variant functions know not to look for data
	  // in a BLOB	  
	  pv->consensus.vcf_direct = true;
	  
	  // Apply all mask filters, and decide whether to call function
//	  if ( eval_and_call( mask, &imap , *pv , f , data ) ) 
	  if ( eval_and_call( mask, &(GP->indmap) , *pv , f , data ) ) 
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
  
//  tmpdb.commit();
  GP->vardb.commit();
  
  // If we had to use any positional filters or grouping, now run the 
  // actual iterate on the temporary database we've created
     
  return irep;
}

