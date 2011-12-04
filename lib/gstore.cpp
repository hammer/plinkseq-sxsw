
#include "gstore.h"
#include "genotype.h"
#include "variant.h"
#include "vcf.h"
#include "defs.h"
#include "crandom.h"

#include <iostream>

using namespace Helper;
using namespace std;

GStore * GP;

GStore::GStore(bool r)
  : phmap(&inddb) , 
    vardb(indmap) , 
    perm( indmap, phmap )
{       	
  
  // Global pointer to primary class
  GP = this; 

  FileMap::setTypes();

  std::set_new_handler(NoMem);
  
  // Ensure random number generator is set
  CRandom::srand( time(0) );

  // (and also the default C RNG)
  srand( time(0) );
  
  // set some standard meta-types (i.e VCF defaults)
  PLINKSeq::register_standard_metatypes();
  
  // this may get changed later, but assume we have a full project
  in_single_file_mode = false;

  // has a proj specification file been set (or ".")
  has_projfile = false;

  // by default, not in R mode
  r_mode = false;
  
}

bool GStore::set_project( const std::string & filename, bool verbose)
{
  
  // We allow a 'dummy' project to be specified, where we
  // do not attempt to attach any databases; this can be 
  // useful, e.g. if a client subsequently just wants to 
  // operate on one database, e.g. SEQDB, without having 
  // a full project in existence

  if ( filename == "." ) 
    {
      // so we know not to attempt certain other functions
      has_project_file( false );
      
      return true;
    }
  
  has_project_file(true);
  
  // TODO: Reserve "-" to mean read from STDIN

  // Starting a new project, we need an index file that points to the
  // key files in the project.
  
  if ( ! fileExists( filename ) ) return false;
  
  // Read in core files and folder locations

  fIndex.setCoreFiles( filename );  

  // Ensure that a RESOURCES/ folder has been specified

  // Read in all other files specified at this point
  
  fIndex.readFileIndex( filename );  

  // Set up core databases

  vardb.attach( fIndex.file( VARDB )->name() );
  inddb.attach( fIndex.file( INDDB )->name() );  
  locdb.attach( fIndex.file( LOCDB )->name() );
  refdb.attach( fIndex.file( REFDB )->name() );
  seqdb.attach( fIndex.file( SEQDB )->name() );

  if ( verbose ) summary( false );
  
  return true;

}

bool GStore::register_mask( const Mask & m )
{
  // set up the individual map  
  int n = indmap.populate( vardb , phmap , m );

  // register all known meta-types (Var, VarFilter, Gen)
  vardb.set_mask_metatypes( m );

  return n > 0;
}

std::string GStore::summary( bool ugly  )
{
  
  std::stringstream ss;
  
  ss << fIndex.summary( ugly );

  ss << "\n";

  if ( vardb.attached() ) 
    ss << vardb.summary( NULL , ugly ) << "\n";
  else 
    ss << "\n -- VARDB not attached -- \n";

  ss << "\n";
  
  if ( inddb.attached() ) 
    ss << inddb.summary( ugly ) << "\n";
  else
    ss << "\n -- INDDB not attached -- \n";

  ss << "\n";

  if ( locdb.attached(  ) ) 
    ss << locdb.summary( ugly ) << "\n";
  else
    ss << "\n -- LOCDB not attached -- \n";

  ss << "\n";
  
  if ( refdb.attached( ) ) 
    ss << refdb.summary( ugly ) << "\n";
  else
    ss << "\n -- REFDB not attached -- \n";

  ss << "\n";
  
  if ( seqdb.attached( ) ) 
    ss << seqdb.summary( ugly ) << "\n";
  else
    ss << "\n -- SEQDB not attached -- \n";
  
  ss << "\n";
  
  ss << Helper::metatype_summary( ugly );
  
  ss << "\n";
  
  return ss.str();

}




bool GStore::vardb_load_vcf( Mask & mask , 
			     const std::set<std::string> & includes, 
			     const std::set<std::string> & excludes, 
			     const std::string * region_mask )
{
  
  vardb.drop_index();
  
  // Have any filters been specified (i.e. only load variants
  // from VCF that belong in a certain LOCDB group?
  
  std::set<Region> filter;
  if ( region_mask )
    filter = locdb.get_regions( *region_mask );
  
  std::set<Region> * pfilter = region_mask ? &filter : NULL  ;

  std::set<File*> files = fIndex.get( VCF );
  std::set<File*>::iterator i = files.begin();
  while ( i != files.end() )
    {           
      std::string filename = (*i)->name();      
      // Do not reload VCF files already in VARDB
      if ( vardb.fileID( (*i)->name() ) == 0 )
	if ( ! vardb_load_vcf( filename , (*i)->tag() , (*i)->comment() , mask , includes , excludes , pfilter ) )
	  return false;
      ++i;
    }
  
  vardb.index();

  // Insert summary Ni, Nv into database for this file  
  i = files.begin();
  while ( i != files.end() )
    {      
      int2 niv = vardb.make_summary( (*i)->name() ) ;  
      plog << (*i)->name() << " : inserted " << niv.p2 << " variants\n";
      ++i;
    }

  return true;
}


bool GStore::vardb_load_vcf( const std::string & file, 
			     const std::string & tag , 
			     const std::string & comment, 
			     Mask & mask , 
			     const std::set<std::string> & includes, 
			     const std::set<std::string> & excludes, 
			     const std::set<Region> * pfilter )
{
    
  // Load, parse VCF file; store variant and genotype information, and
  // meta-information, in vardb
  
  File * f = fIndex.file( file );
  if ( f == NULL ) return false;
    
  VCFReader v( f , tag , &vardb , &seqdb );

  // Selectively filter in/out meta-information?
  if ( includes.size() > 0 ) v.get_meta( includes );
  if ( excludes.size() > 0 ) v.ignore_meta( excludes );
  
  // add a region filter?
  if ( pfilter ) v.set_region_mask( pfilter );

  // If SEQDB has been disabled
  if ( ! GP->seqdb.attached() ) v.set_seqdb( NULL );
  
  // Do we want to apply up-front fix for X/Y genotypes?
  if ( mask.fixxy() )
    {
      v.set_fixxy( &mask , &locdb, &inddb );
    }

  vardb.begin();
  
  int inserted = 0;

  plog.counter( "parsing..." );

  while ( v.parseLine() ) 
    { 
      if ( ++inserted % 1000 == 0 ) 
	plog.counter( "parsed " + Helper::int2str( inserted ) + " rows" );
    }
  plog.counter("\n");
    
  // Wrap up
  
  vardb.commit();    
  
  return true;
}



void GStore::show_version() const
{
  std::map<std::string,std::string> v = version();
  std::map<std::string,std::string>::iterator i = v.begin();
  while ( i != v.end() )
    {
      plog << i->first << "\t" << i->second << "\n";
      ++i;
    }
}


std::map<std::string,std::string> GStore::version() const
{
  std::map<std::string,std::string> v;
  v[ "PLINKSeq" ] = PLINKSeq::VERSION_NUMBER();
  v[ "SQLITE3_LIBRARY" ] = SQL::library_version();
  v[ "SQLITE3_HEADER" ] = SQL::header_version();  
  return v;
}


LocDBase * GStore::resolve_locgroup( const std::string & g )
{
  
  // does group have explicit specifier? e.g. LOCDB::CCDS

  if ( g.size() > 7 )
    {
      if ( g.substr( 0 , 7 ) == "LOCDB::" ) 
	return locdb.lookup_group_id( g.substr( 7 ) ) == 0 ? NULL : &locdb; 
      
      if ( g.substr( 0 , 7 ) == "SEGDB::" ) 
	return segdb.lookup_group_id( g.substr( 7 ) ) == 0 ? NULL : &segdb;       
    }


  // Give preference to user-defined SEGDB, if no explicit qualifier

  if ( segdb.lookup_group_id( g ) != 0 ) return &segdb;

  if ( locdb.lookup_group_id( g ) != 0 ) return &locdb;
  
  return NULL;
}


