#include "plinkseq/vardb.h"
#include "plinkseq/vcf.h"
#include "plinkseq/gstore.h"

extern GStore * GP;

bool VarDBase::vcf_iterate_read_header( Mask & mask )
{

  // Read meta-information, and header row.
  // Populate the vardb, indmap, etc

  std::string filename = mask.external_vcf_filename();

  if ( filename != "-" ) Helper::checkFileExists( filename );
  
  //
  // Use VCFReader, into a temporary :memory: database
  //

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
  
  //
  // In FIX/XY mode when in single-VCF mode, we need to populate the Sex codes in the indmap
  //

  for (int i = 0; i< GP->indmap.size() ; i++)
    {
      sType s = GP->inddb.sex( GP->indmap.ind(i)->id() );
      GP->indmap.ind(i)->sex( s );
    }
  
  
}


IterationReport VarDBase::vcf_iterate( void (*f)(Variant&, void *) , void * data , Mask & mask )
{
  

  IterationReport irep( true , mask.any_grouping() , mask.variant_limit() );


  // VCF file name is kept in the Mask, by the 'ex-vcf' attribute
  
  std::string filename = mask.external_vcf_filename();
  
  if ( filename != "-" ) Helper::checkFileExists( filename );


  // Use VCFReader, into a temporary :memory: database

  // Load, parse VCF file; store variant and genotype information, and
  // meta-information, in vardb


  File vcffile( filename , VCF );

  VCFReader v( &vcffile , "" , &(GP->vardb) ,  NULL );

  // hack 
  if ( filename == "-" ) 
  {
      v.observed_header( true );      
      v.set_number_individuals( GP->indmap.size() );
  }

      
  if ( mask.fixxy() )
    {
      v.set_fixxy( &mask , &(GP->locdb), &(GP->inddb) );
    }


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
  
  if ( mask.id() )
    v.add_id_filter( mask.included_id() );
    

  // Misc. settings.
  
  downcode_mode = mask.downcode();


  // Work through VCF
  
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
      

      // If a variant line has been processed and meets criteria, pv
      // will be non-NULL which also implies that a Variant has been
      // created and we are reponsible for cleaning up afterwards
      
      if ( l == VCFReader::VCF_VARIANT && pv )
	{
	  
	  // bad line, or failed a loc mask filter

	  if ( ! pv->valid() ) 
	    {
	      irep.rejected_variant();
	      delete pv;
	      continue;
	    }
	  	  
	  // So that the Variant functions know not to look for data
	  // in a BLOB; also, they they know how to parse it downstream
	  
	  pv->set_vcf_buffer( v.gt_field , &v.formats );
	  
	  
	  // Apply all mask filters, and decide whether to call function
	  
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
  
  GP->vardb.commit();
  
  // If we had to use any positional filters or grouping, now run the 
  // actual iterate on the temporary database we've created
     
  return irep;
}

