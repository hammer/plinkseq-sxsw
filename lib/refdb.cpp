
#include "plinkseq/refdb.h"
#include "plinkseq/variant.h"
#include "plinkseq/gstore.h"
#include "plinkseq/filemap.h"
#include "plinkseq/output.h"

#include <algorithm>
#include <cmath>

extern GStore * GP;

using namespace Helper;
using namespace std;

bool RefDBase::attach(std::string name)
{

  if ( name == "-" || name == "." ) { dettach(); return false; } 
  
  // If this file does not exist, create it, otherwise attach
  
  sql.open( name );

  fname = name;

  sql.synchronous(false);


  //
  // DB version, and a place for various other meta-information in future
  //
  
  sql.query(" CREATE TABLE IF NOT EXISTS dbmeta("
            "   varname      VARCHAR(20) NOT NULL , "
            "   varvalue    VARCHAR(20) NOT NULL , "
            " CONSTRAINT uMeta UNIQUE (varname ) ); " );

  
  sql.query(" CREATE TABLE IF NOT EXISTS refvariants("
            "   group_id  INTEGER NOT NULL , "
            "   name      VARCHAR(40) , "
            "   chr       INTEGER NOT NULL , "
            "   bp1       INTEGER NOT NULL , "
            "   bp2       INTEGER NOT NULL , "
	    "   ref       VARCHAR(1) , "
	    "   alt       VARCHAR(1) , "
	    "   value     VARCHAR(10) ); " );

  
  // Meta-information table
  
  sql.query( " CREATE TABLE IF NOT EXISTS metatypes("
	     "   field_id     INTEGER PRIMARY KEY , "
	     "   group_id     INTEGER NOT NULL , "
	     "   name         VARCHAR(8) , "
 	     "   type         VARCHAR(8) , "
 	     "   number       INTEGER , "
 	     "   description  VARCHAR(20) ); ");
  

  
  // Group information table
  
  sql.query( " CREATE TABLE IF NOT EXISTS groups("
             "   group_id     INTEGER PRIMARY KEY , "
	     "   count        INTEGER , "
             "   name         VARCHAR(40) NOT NULL , "
             "   temp         CHAR(1) , "
             "   description  TEXT ); " );



  
  //
  // set up queries
  //

  init();
  
  //
  // Check version number 
  //

  check_version();


  //
  // set meta-type information
  //

  set_metatypes();

  index();

}


void RefDBase::set_metatypes( bool clear )
{
  if ( clear ) 
    {
      MetaInformation<RefMeta>::reset();    
    }
  
  while ( sql.step( stmt_fetch_metatypes ) )
    {
      std::string name = sql.get_text( stmt_fetch_metatypes , 0 );
      mType mt = (mType)sql.get_int( stmt_fetch_metatypes , 1 );
      int num = sql.get_int( stmt_fetch_metatypes , 2 );
      std::string desc = sql.get_text( stmt_fetch_metatypes , 3 );	
      registerMetatype( name, mt , num, META_GROUP_REF, desc );
      registerMetatype( name, mt , num, META_GROUP_VAR, desc );
    }
  sql.reset( stmt_fetch_metatypes );
  
  while ( sql.step( stmt_fetch_groups ) )
    {
      std::string gname = sql.get_text( stmt_fetch_groups , 0 );
      int gid = sql.get_int( stmt_fetch_groups , 1 );

      registerMetatype( gname, META_FLAG , -1, META_GROUP_REF, "" );
      registerMetatype( gname, META_FLAG , -1, META_GROUP_VAR, "" );
      
      grpmap[gid] = gname;
    }
  sql.reset( stmt_fetch_groups ) ;
}




void RefDBase::drop_index()
{
  if ( ! attached() ) return;
  sql.query( "DROP INDEX IF EXISTS ind1;");
  sql.query( "DROP INDEX IF EXISTS ind2;");
  release();
  init();
}


void RefDBase::index()
{
  if ( ! attached() ) return;  
  sql.query( "CREATE INDEX IF NOT EXISTS ind1 ON refvariants(group_id,chr, bp1); " );    
  sql.query( "CREATE INDEX IF NOT EXISTS ind2 ON refvariants(group_id,name); " );    
  release();
  init();  
}

bool RefDBase::init()
{
    
  //
  // Core
  //

  stmt_lookup = 
    sql.prepare( "SELECT * FROM refvariants WHERE chr == :chr AND bp1 == :bp1 AND :bp2 == :bp2 AND group_id == :group_id ; "  );

  stmt_lookup_allelic_match = 
    sql.prepare( "SELECT * FROM refvariants WHERE chr == :chr AND bp1 == :bp1 AND bp2 == :bp2 AND group_id == :group_id AND alt == :alt ; "  );

  stmt_lookup_range = 
    sql.prepare( "SELECT * FROM refvariants "
		 " WHERE chr == :chr AND bp1 >= :rstart AND group_id == :group_id AND bp1 <= :rend ; "  );

  stmt_lookup_range_count = 
    sql.prepare( "SELECT COUNT(*) FROM refvariants "
		 " WHERE group_id == :group_id AND chr == :chr AND :bp2 >= bp1 AND :bp1 <= bp2 ; "  );

  stmt_lookup_dbsnp = 
    sql.prepare( "SELECT * FROM refvariants WHERE group_id == :group_id AND name == :name ; "  );

  stmt_insert =
    sql.prepare( "INSERT OR REPLACE INTO refvariants ( name, group_id, chr, bp1, bp2, ref, alt, value ) "
		 " values( :name , :group_id, :chr, :bp1, :bp2, :ref, :alt , :value ) ; " );

  stmt_dump = 
    sql.prepare( "SELECT * FROM refvariants WHERE group_id == :group_id ORDER BY chr,bp1; " );
  

  //
  // Meta-information
  //


  stmt_meta_dump =
    sql.prepare(" SELECT name,field_id,type FROM metatypes ;" );

  stmt_meta_insert_prep =
    sql.prepare(" SELECT field_id,type FROM metatypes WHERE name == :name ;" );

  stmt_meta_insert_prep2 =
    sql.prepare(" INSERT INTO metatypes (name,type,number,description,group_id) values( :name, :type, :number, :description, :group_id ); " );

  stmt_fetch_metatypes = 
    sql.prepare(" SELECT name , type , number, description , group_id "
		  " FROM metatypes ; " );


  //
  // Groups
  //

  stmt_insert_group_name =
      sql.prepare("INSERT OR REPLACE INTO groups ( name, temp, description ) "
		  " values( :name, :temp, :description ) ; ");
  
  stmt_lookup_group_name =
    sql.prepare("SELECT group_id FROM groups WHERE name == :name ; " );

  stmt_remove_group1 =
    sql.prepare("DELETE FROM refvariants WHERE group_id == :group ; " );
  
  stmt_remove_group2 =
      sql.prepare("DELETE FROM groups WHERE group_id == :group ; " );
  
  stmt_update_group_count = 
      sql.prepare("UPDATE groups "
		  "SET count = :count "
		  "WHERE group_id == :group_id ; ");

  stmt_fetch_group_count = 
      sql.prepare("SELECT count FROM groups WHERE group_id == :group_id; ");
  
  stmt_fetch_groups = 
      sql.prepare(" SELECT name , group_id FROM groups ; " );



}


bool RefDBase::release()
{
    
  // Main queries

  sql.finalise( stmt_dump );
  sql.finalise( stmt_lookup );
  sql.finalise( stmt_lookup_allelic_match );
  sql.finalise( stmt_insert );
  sql.finalise( stmt_lookup_dbsnp );
  sql.finalise( stmt_lookup_range );
  sql.finalise( stmt_lookup_range_count );
  
  // Meta-information

  sql.finalise( stmt_meta_dump );
  sql.finalise( stmt_meta_insert_prep );
  sql.finalise( stmt_meta_insert_prep2 );
  sql.finalise( stmt_fetch_metatypes );
  sql.finalise( stmt_fetch_groups );
  
  // Groups
    
  sql.finalise( stmt_insert_group_name );
  sql.finalise( stmt_lookup_group_name );
  sql.finalise( stmt_remove_group1 );
  sql.finalise( stmt_remove_group2 );
  sql.finalise( stmt_update_group_count );
  sql.finalise( stmt_fetch_group_count );
  
}



bool RefDBase::dettach()
{
  if ( attached() )
    {
      release();
      sql.close();  
    }
  return true;
}


//
// Main insertion functions
//

void RefDBase::insert_header( uint64_t file_id , const std::string & name , const std::string & value )
{
  // ignored for now...
  // populate table metaheaders...
}


void RefDBase::insert_metatype( uint64_t file_id ,
				const std::string & name ,
				mType mt ,
				int num ,
				int mgrp ,
				const std::string & description )
{
  
  // pre-prend group name on meta-type name
  
  std::string name2 = grpmap[ file_id ] + "_" + name;
  
  sql.bind_text( stmt_meta_insert_prep , ":name" , name2 );
  
  if ( sql.step( stmt_meta_insert_prep ) )
    {
      // store ID for this meta-tag if already seen
      mtmap[ name2 ] = sql.get_int64( stmt_meta_insert_prep , 0 );
    }
  else
    {	      
      // otherwise, need to store this in DB
      
      sql.bind_text( stmt_meta_insert_prep2 , ":name" , name2 );
      sql.bind_int( stmt_meta_insert_prep2 , ":type" , mt );
      sql.bind_int( stmt_meta_insert_prep2 , ":number" , num );
      sql.bind_text( stmt_meta_insert_prep2 , ":description" , description );
      sql.bind_int( stmt_meta_insert_prep2 , ":group_id" , (int)file_id );
      sql.step( stmt_meta_insert_prep2 );

      // store the newly-created field ID locally
      mtmap[ name2 ] = sql.last_insert_rowid();

      // clean up
      sql.reset(stmt_meta_insert_prep2);
    }
  // clean up
  sql.reset( stmt_meta_insert_prep );
  
  
  //
  // Populate mtmap[] with any FLAGs that might already exist in the DB
  //

  while ( sql.step( stmt_meta_dump ) )
    {
      std::string fname = sql.get_text( stmt_meta_dump , 0 );
      int fid = sql.get_int( stmt_meta_dump , 1 );
      int type = sql.get_int( stmt_meta_dump , 2 );      
      if ( mtmap.find( fname ) == mtmap.end() )
	mtmap[ fname ] = fid;
    }

  sql.reset( stmt_meta_dump ) ;


}


bool RefDBase::load_VCF( const std::string & file, 
			 const std::string & tag , 
			 const std::string & comment, 
			 const std::set<std::string> & includes, 
			 const std::set<std::string> & excludes, 
			 const std::set<Region> * pfilter )
{
  
  // Load, parse VCF file; only store variant data in REFDB. Ignore any genotype data
  
  File * f = GP->fIndex.add( file , VCF , tag , comment );
  
  // Point VCF-reader to a REFDB
  
  VCFReader v( f , tag , this , &(GP->seqdb) );
  
  // Selectively filter in/out meta-information?
  if ( includes.size() > 0 ) v.get_meta( includes );
  if ( excludes.size() > 0 ) v.ignore_meta( excludes );
  
  // add a region filter?
  if ( pfilter ) v.set_region_mask( pfilter );
  
  drop_index();
  sql.begin();

  int inserted = 0;
  
  while ( v.parseLine() ) 
    { 
      if ( ++inserted % 1000 == 0 ) 
	plog.counter1( "parsed " + Helper::int2str( inserted ) + " rows" );
    }
  plog.counter1("\n");
    
  // Wrap up
  
  sql.commit();    
  index();

  // insert how many variants added
  
  sql.bind_int64( stmt_update_group_count, ":group_id" , v.group_id()  );
  sql.bind_int( stmt_update_group_count, ":count" , v.variants_inserted()  );
  sql.step( stmt_update_group_count );
  sql.reset( stmt_update_group_count );
    
  return true;
}




bool RefDBase::insert( int file_id , const Variant & v )
{
  
  // Any variant meta-information (will be in consensus SampleVariant)
  
  std::stringstream ss;
  ss << v.consensus.meta;

  RefVariant rv( file_id , 
		 v.name() , 
		 v.chromosome() , 
		 v.position() , 
		 v.stop() , 
		 v.reference() , 
		 v.alternate() , 
		 ss.str()  
		 );
  
  
  
  //  attach_metainformation( rv , v );
  
  // Insert reference variant

  refInsertion( rv );
  
}


bool RefDBase::refInsertion(const RefVariant & rv)
{


    // Insert variant:
    
    sql.bind_text( stmt_insert , ":name" ,   rv.name() );
    sql.bind_int( stmt_insert , ":group_id", rv.group() );
    sql.bind_int( stmt_insert , ":chr" ,     rv.chromosome() );
    sql.bind_int( stmt_insert , ":bp1" ,     rv.start() );

    sql.bind_int( stmt_insert , ":bp1" ,     rv.start() );
    sql.bind_int( stmt_insert , ":bp1" ,     rv.start() );

    sql.bind_text( stmt_insert , ":ref" ,   rv.reference() );
    sql.bind_text( stmt_insert , ":alt" ,   rv.alternate() );
    
    sql.bind_int( stmt_insert , ":bp2" ,     rv.stop() );

    sql.bind_text( stmt_insert , ":value" ,  rv.value() );
    sql.step( stmt_insert );
    uint64_t id = sql.last_insert_rowid();
    sql.reset( stmt_insert );
            
}


std::map<std::string,mType> RefDBase::populate_metatypes( std::map<std::string,int> * meta , const int group_id )
{
  
  std::map<std::string,mType> mt;
  
  if ( meta ) 
    {
      
      std:map<std::string,int>::iterator i = meta->begin();
      
      while ( i != meta->end() )
	{
	  
	  mt[ i->first ] = MetaInformation<RefMeta>::type( i->first );
	  
	  sql.bind_text( stmt_meta_insert_prep , ":name" , i->first );
	  if ( sql.step( stmt_meta_insert_prep ) )
	    {
	      mtmap[ i->first ] = sql.get_int64( stmt_meta_insert_prep , 0 );
	    }
	  else
	    {	      
	      meta_index_t midx = MetaInformation<RefMeta>::field( i->first );
	      sql.bind_text( stmt_meta_insert_prep2 , ":name" , i->first );
	      sql.bind_int( stmt_meta_insert_prep2 , ":type" , midx.mt );
	      sql.bind_int( stmt_meta_insert_prep2 , ":number" , midx.len );
	      sql.bind_text( stmt_meta_insert_prep2 , ":description" , midx.description );
	      sql.bind_int( stmt_meta_insert_prep2 , ":group_id" , group_id );
	      sql.step( stmt_meta_insert_prep2 );
	      mtmap[ i->first ] = sql.last_insert_rowid();
	      sql.reset(stmt_meta_insert_prep2);
	    }
	  sql.reset( stmt_meta_insert_prep );
	  ++i;
	}
    }
  

  //
  // Populate mtmap[] with any FLAGs that might already exist in the DB
  //

  while ( sql.step( stmt_meta_dump ) )
    {
      std::string name = sql.get_text( stmt_meta_dump , 0 );
      int fid = sql.get_int( stmt_meta_dump , 1 );
      int type = sql.get_int( stmt_meta_dump , 2 );

      if ( mtmap.find( name ) == mtmap.end() )
	mtmap[ name ] = fid;
    }

  sql.reset( stmt_meta_dump ) ;
  
  return mt;
}




void RefDBase::attach_metainformation( RefVariant & rv , const Variant & v )
{
  rv.meta.append( v.consensus.meta , "" );
}


void RefDBase::attach_metainformation( RefVariant & rv , 
				       std::map<std::string,int> * meta , 
				       std::set<int> * flags , 
				       std::set<int> * skip , 
				       std::map<std::string,mType> & mt , 
				       std::vector<std::string> & buffer ) 
{
  
  std::string ms = "";

  if ( meta ) 
    {
      std::map<std::string,int>::iterator m = meta->begin();  
      
      while ( m != meta->end() )
	{		
	  
	  if ( skip )
	    if ( skip->find( m->second ) != skip->end() )
	      {
		++m;
		continue;
	      }
	  
	  mType t = mt[ m->first ];
	  
	  if ( ms != "" ) ms += ";";

	  ms += m->first;
	  
	  if ( t != META_FLAG ) 
	    {	      
	      if ( buffer[m->second].find(" ") != std::string::npos ) 		
		ms += "=\"" + buffer[ m->second ] + "\"";
	      else
		ms += "=" + buffer[ m->second ];
	    }
	  
	  ++m;
	}
    }

  
  //
  // Flag meta-information in VCF format -- just append
  //
  
  if ( flags )
    {      
      // Assume a semi-colon delimited list of flags we wish to set
      // in VCF format.  i.e. this can include 

      // key=value                (normal tag)
      // key=value1,value2        (vector)
      // key                      (flags)
      
      // we do not care about type for now

      std::set<int>::iterator m = flags->begin();  
      while ( m != flags->end() )
	{	
	  if ( ms != "" ) ms += ";";
	  ms += buffer[ *m ];
	  ++m;
	}
    }

  // set the meta-info
  rv.value( ms );
  
}

uint64_t RefDBase::insert( const std::string & filename , const std::string & grp )
{
  return set_group_id( grp );  
}

uint64_t RefDBase::loadRefVariants(const std::string & filename, 
				   const std::string & grp, 
				   int col_name,
				   int col_chr,
				   int col_bp1,
				   map<std::string,int> * meta ,
				   std::set<int> * flags , 
				   std::set<int> * skip , 
				   int col_bp2 ,
				   bool zero1 , 
				   bool zero2 )
{

  // for now, no allele codes allowed as input from flat files...
  std::string ref = ".";
  std::string alt = ".";


  if ( ! attached() ) return 0;

  if ( col_bp1 == -1 || col_chr == -1 ) 
    {
      Helper::halt("no CHR or BP1 specified");
    }

  //
  // Create a new group for this class of reference
  //
  
  if ( ! Helper::valid_name( grp ) ) 
    Helper::halt( grp + " is not a valid name" );

  uint64_t group_id = set_group_id( grp );
  
  bool names = col_name > 0 ;
  
  int maxcol = col_name;
  if ( col_name > maxcol ) maxcol = col_name;
  if ( col_chr > maxcol ) maxcol = col_chr;
  if ( col_bp1 > maxcol ) maxcol = col_bp1;
  if ( col_bp2 > maxcol ) maxcol = col_bp2;

  if ( meta )
    {
      map<std::string,int>::iterator i = meta->begin();
      while ( i != meta->end() )
	{
	  if ( i->second > maxcol )
	    maxcol = i->second;
	  
	  ++i;
	}
    }
  
  if ( flags )
    {
      std::set<int>::iterator i = flags->begin();
      while ( i != flags->end() )
	{
	  if ( *i > maxcol ) maxcol = *i;	  
	  ++i;
	}
    }

  ++maxcol;
  
  
  //                                     
  // Begin SQL transaction               
  //        
  
  sql.begin();
  
  int inserted = 0;


  //
  // Get and write meta-information types
  //
  
  std::map<std::string,mType> mt = populate_metatypes( meta , group_id );
  
  
  //                                     
  // Process each row of input           
  //                                     

  
  Helper::checkFileExists( filename );
  
  InFile file( filename );
  
  drop_index();
  
  while ( ! file.eof() )
    {
      
      // Assume a tab-delimited file
      
      std::vector<std::string> buffer = file.tokenizeLine("\t");
      
      // Comment line?  
      
      if ( buffer.size() > 0 && buffer[0].substr(0,1) == "#" ) 
	continue;
      
      // Invalid formatting?
      
      if( (int)buffer.size() < maxcol )
	{
	  if ( ! file.eof() ) 
	    plog.warn( "too few rows for row", buffer );
	  continue;
	}
      
      int p1,p2;
      
      // Here, allow variants with unknown positions
      // i.e. as we can lookup based on the name, sometimes
      // Ensure that ALL co-ordinates are 1-based
      
      if ( Helper::str2int( buffer[ col_bp1 ] , p1 ) ) 
	{
	  if ( zero1 ) ++p1;
	}
      else
	p1 = 0;
      
      if ( col_bp2 == -1 )
	p2 = p1;
      else
	{
	  if ( str2int( buffer[ col_bp2 ] , p2 ) ) 
	    {
	      if ( zero2 ) ++p2;
	    }
	  else 
	    p2 = 0;
	}
      

      std::string name = names ? buffer[ col_name ] : "";	
      
      int chromosome = chrCode( buffer[ col_chr ] ) ;
      
      RefVariant rv( (int)group_id, 
		     name ,
		     chromosome ,
		     p1 , p2 ,
		     ref , alt , 
		     "" );

      
      //
      // Meta-information ("value" field of RefVariant)
      //
      
      attach_metainformation( rv , meta , flags , skip , mt , buffer );
       
      
	
      //
      // Add entry to database
      //


      refInsertion(rv);
	
      ++inserted;
      
      if ( inserted % 1000 == 0 ) 
	{
	  plog << "inserted " << inserted << " reference-variants              \r";
	  plog.flush();
	}
    }
  
  plog << "\n";

  //                  
  // Finish transaction                  
  //                                     
  
  sql.commit();
  
  
  //
  // Store DB meta-information in reference table
  //
  
  sql.bind_int64( stmt_update_group_count, ":group_id" , group_id  );
  sql.bind_int( stmt_update_group_count, ":count" , inserted  );
  sql.step( stmt_update_group_count );
  sql.reset( stmt_update_group_count );
  
  plog << filename << " : inserted " << inserted << " reference variants\n";
  
  file.close();
  
  // Re-index
  index();
  
  return group_id;
}
   




uint64_t RefDBase::set_group_id(const std::string grp, 
				const bool temp ,
				const std::string desc )
{
  
  uint64_t group_id = 0;

  sql.bind_text( stmt_lookup_group_name, ":name" , grp );

  if ( sql.step( stmt_lookup_group_name ) )
    {
      group_id = sql.get_int64( stmt_lookup_group_name , 0 ) ;
      sql.reset( stmt_lookup_group_name );
    }
  else
    {
      sql.reset( stmt_lookup_group_name );
      
      sql.bind_text( stmt_insert_group_name , ":name" , grp );
      sql.bind_int( stmt_insert_group_name , ":temp" , (int)temp );
      sql.bind_text( stmt_insert_group_name , ":description" , grp + " (default name)" );
      sql.step( stmt_insert_group_name );
      group_id = sql.last_insert_rowid();
      sql.reset( stmt_insert_group_name );
      
      grpmap[ group_id ] = grp;
      
    }

  return group_id;

}

uint64_t RefDBase::lookup_group_id(const std::string grp)
{
    uint64_t group_id = 0;  
    sql.bind_text( stmt_lookup_group_name, ":name" , grp );
    if ( sql.step( stmt_lookup_group_name ) )
	{
	    group_id = sql.get_int64( stmt_lookup_group_name , 0 ) ;
	}
    sql.reset( stmt_lookup_group_name );   
    return group_id;
}


void RefDBase::flush(const uint64_t id)
{
  // Remove all temporary groups
  // (for now, ignore group-specific request)

  sql.query("DELETE FROM refvariants WHERE group_id IN ( SELECT group_id FROM groups WHERE temp == 1 ); ");
  sql.query("DELETE FROM groups WHERE temp == 1 ;");
}



//
//  Queries
//


bool RefDBase::init_iterate( const std::string & grp )
{
  if ( ! attached() ) return false;
  int g = lookup_group_id( grp );
  if ( g == 0 ) return false;
  sql.bind_int( stmt_dump , ":group_id" , g );
  return true;
}

void RefDBase::construct_inplace( sqlite3_stmt * s , RefVariant * rv )
{
  rv->group(      sql.get_int( s , 0 ) );	    
  rv->name(       sql.get_text( s , 1 ) );
  rv->chromosome( sql.get_int( s , 2 ) ) ;
  rv->start(      sql.get_int( s , 3 ) ) ;
  rv->stop(       sql.get_int( s , 4 ) ) ;
  rv->reference(  sql.get_text( s , 5 ) ) ;
  rv->alternate(  sql.get_text( s , 6 ) ) ;
  rv->value(      sql.get_text( s , 7 ) );  
}


bool RefDBase::iterate( RefVariant * rv )
{  
  if ( sql.step( stmt_dump ) )
    {
      construct_inplace( stmt_dump , rv );      
      rv->observed( true ) ;
      return true;
    }
  else
    {      
      rv->observed( false ) ;
      sql.reset( stmt_dump );
      return false;
    }
}


int RefDBase::count( const Region & region , const std::string & grp_name  )
{

  if ( ! attached() ) return -1;
  int grp = lookup_group_id( grp_name );  
  if ( grp == 0 ) return 0;

  sql.bind_int( stmt_lookup_range_count , ":group_id" , grp );
  sql.bind_int( stmt_lookup_range_count , ":chr" , region.chromosome() );
  sql.bind_int( stmt_lookup_range_count , ":bp1" , region.start.position() );
  sql.bind_int( stmt_lookup_range_count , ":bp2" , region.stop.position() );
  
  int c = 0;
  
  if ( sql.step( stmt_lookup_range_count ) )
    c = sql.get_int( stmt_lookup_range_count , 0 );
  
  sql.reset( stmt_lookup_range_count ) ;

  return c;

}


std::set<RefVariant> RefDBase::lookup( const Variant & v , const std::string & grp_name , bool allelic_match )
{
    if ( ! attached() )
    	return std::set<RefVariant>();

    return lookup( v, lookup_group_id( grp_name ) , allelic_match );
}

std::set<RefVariant> RefDBase::lookup( const Variant & v , int grp_id , bool allelic_match )
{
  
  //
  // Search given chromosome position, and optionally, match on ALT allele
  //

  sqlite3_stmt * s = allelic_match ? stmt_lookup_allelic_match : stmt_lookup ;
  
//   std::cout << "looking for " << v.chromosome() << " " 
// 	    << v.position() << " " 
// 	    << v.stop() << " "
// 	    << grp_id << " " 
// 	    << v.alternate() << " " << allelic_match << "\n";

  std::set<RefVariant> refSet;

  std::set<string> altsToMatch;
  if (allelic_match) {
	  altsToMatch.insert(v.alternate()); // try just in case the consensus does not exist
	  for (int k = 1; k <= v.consensus.nalt(); ++k) // exclude the reference [at k = 0]
		  altsToMatch.insert(v.consensus.alternate(k));
  }
  else
	  altsToMatch.insert(".");

  for (std::set<string>::const_iterator altIt = altsToMatch.begin(); altIt != altsToMatch.end(); ++altIt) {
	  sql.bind_int( s , ":chr" , v.chromosome() );
	  sql.bind_int( s , ":bp1" , v.position() );
	  sql.bind_int( s , ":bp2" , v.stop() );
	  sql.bind_int( s , ":group_id" , grp_id );
	  if ( allelic_match )
		  sql.bind_text( s , ":alt" , *altIt );

	  while ( sql.step( s ) )
	  {
		  RefVariant r = construct( s );
		  if ( r.observed() )
			  refSet.insert( r );
	  }
	  sql.reset( s );
  }
  
  return refSet;
}


std::set<RefVariant> RefDBase::lookup( const Region & region , const int grp_id , const int limit )
{
  
  // Search given chromosome position
  
  sql.bind_int( stmt_lookup_range , ":chr" , region.chromosome() );
  sql.bind_int( stmt_lookup_range , ":rstart" , region.start.position() );
  sql.bind_int( stmt_lookup_range , ":rend" , region.stop.position() );
  sql.bind_int( stmt_lookup_range , ":group_id" , grp_id );

  int cnt = 0;
  std::set<RefVariant> s;
  while ( sql.step( stmt_lookup_range ) )
    {      
      RefVariant r = construct( stmt_lookup_range );      
      if ( r.observed() ) s.insert( r );      
      if ( limit ) { ++cnt; if ( cnt == limit ) break; } 
    }
  sql.reset( stmt_lookup_range );
  return s;  
}

std::set<RefVariant> RefDBase::lookup( const Region & region , const std::string & grp_name , const int limit )
{
  std::set<RefVariant> empty;
  if ( ! attached() ) return empty;
  return lookup( region, lookup_group_id( grp_name ) , limit );
}



bool RefDBase::annotate( Variant & v , const int grp_id , bool allelic_match )
{
	std::set<RefVariant> rSet = lookup( v , grp_id , allelic_match );
	bool annotated = false;

	for (std::set<RefVariant>::const_iterator rIt = rSet.begin(); rIt != rSet.end(); ++rIt) {
		RefVariant r = *rIt;
		if (!r.observed())
			continue;

		annotated = true;

		// Attach RefVariant 'name' and 'value' fields

		// the below is likely not terribly efficient
		// but leave as the entire meta-info class
		// needs a revamping

		const std::string gname = grpmap[ grp_id ];

		MetaInformation<VarMeta>::field( gname , META_FLAG );
		v.meta.set( gname );

		if ( r.name() != "" && r.name() != "." )
			v.meta.set( gname + "_ID" , r.name() );

		// attach meta information stored in 'value', but first
		// appending the prefix

		if ( r.value() != "" && r.value() != "." )
		{
			// value() contains non-prefixed versions (e.g. A=1;U=0)
			// passing this extra param, means the prefixes get added,
			// which is necessary, because we store and append them this
			// way

			r.meta.parse( r.value() , ';' , false ,  &gname  );

			v.meta.append( r.meta );           // now already done internally

		}
	}

	return annotated;
}


bool RefDBase::annotate( Variant & v , const std::string & name , bool allelic_match )
{
  int id = lookup_group_id( name );
  if ( id == 0 ) return false;
  return annotate(v,id , allelic_match );
}

RefVariant RefDBase::construct( sqlite3_stmt * s)
{
    
    // Primary refvariants table has the following structure:

    // 0 group_id
    // 1 name
    // 2 chr
    // 3 bp1
    // 4 bp2
    // 5 value
  
  int grp_id = sql.get_int( s , 0 );	    
  std::string name = sql.get_text( s , 1 ) ;
  int chr = sql.get_int( s , 2 ) ;
  int bp1 = sql.get_int( s , 3 ) ;
  int bp2 = sql.get_int( s , 4 ) ;
  std::string ref = sql.get_text( s , 5 );
  std::string alt = sql.get_text( s , 6 );
  std::string value = sql.get_text( s , 7 );
  
  return RefVariant(grp_id,name,chr,bp1,bp2,ref,alt,value);

}


std::vector<std::string> RefDBase::fetch_groups()
{
  std::vector<std::string> results;
  if ( ! attached() ) return results;
  sqlite3_stmt * s = 
    sql.prepare( "SELECT name FROM groups ORDER BY group_id; ");
  while ( sql.step(s) )
    {
      std::string name = sql.get_text( s , 0 );
      results.push_back(name);
    }
  sql.finalise(s); 
  return results;
}

std::string RefDBase::summary( bool ugly )
{
    std::stringstream ss;
    
    sqlite3_stmt * s = 
      sql.prepare( "SELECT group_id,name,description FROM groups; ");
    

    if ( ! ugly ) ss << "---Reference DB summary---\n\n";
    
    bool empty = true;

    while ( sql.step(s) )
      {
	empty = false;

 	int id = sql.get_int( s , 0 );
 	std::string name = sql.get_text( s , 1 );
 	std::string desc =  sql.get_text( s , 2 );
	
 	sql.bind_int64( stmt_fetch_group_count , ":group_id", id );

 	sql.step( stmt_fetch_group_count );
 	int c = sql.get_int( stmt_fetch_group_count , 0 );
 	sql.reset(stmt_fetch_group_count);

	if ( ugly ) 
	  ss << "REFDB\t"
	     << "NAME=" << name << "\t"
	     << "N=" << c << "\t"
	     << "DESC=" << desc 
	     << "\n";
	else
	  ss << "Group : " << name << " (" << c << " entries) : " << desc << "\n";

     }
     sql.finalise(s);
     
     if ( empty ) ss << "(empty)\n";

     return ss.str();
}




std::map< std::string, std::map<std::string,std::string> > RefDBase::get_metatypes()
{
  
  std::map< std::string, std::map<std::string,std::string> > mgrp;
  
  while ( sql.step( stmt_fetch_metatypes ) )
    {
      
      std::string name = sql.get_text( stmt_fetch_metatypes , 0 );
      
      mType mt = (mType)sql.get_int( stmt_fetch_metatypes , 1 );
      std::string num = sql.get_text( stmt_fetch_metatypes , 2 ); // get as text
      std::string desc = sql.get_text( stmt_fetch_metatypes , 3 );	
      
      std::map<std::string,std::string> inf;
      
      inf["GRP"] = "RefVariant";
      inf["NAME"] = name;
      inf["LEN"] = num;
      inf["DESC"] = desc;
      
      switch ( mt ) 
	{
	case META_FLAG :
	  inf[ "TYPE" ] = "Flag";
	  break;
	case META_UNDEFINED :
	  inf[ "TYPE" ] = "Undefined";
	  break;
	case META_TEXT :
	  inf[ "TYPE" ] = "String";
	  break;
	case META_INT :
	  inf[ "TYPE" ] = "Integer";
	  break;
	case META_FLOAT :
	  inf[ "TYPE" ] = "Float";
	  break;
	case META_BOOL :
	  inf[ "TYPE" ] = "Bool";
	  break;
	case META_CHAR :
	  inf[ "TYPE" ] = "Char";
	}

      mgrp[name] = inf;
      
    }

  sql.reset( stmt_fetch_metatypes );

  return mgrp;

}



void RefDBase::dump( const std::string & grp , bool with_meta , bool with_verbose )
{

  Out & pout = Out::stream( "refvars" );
  
  if ( ! init_iterate( grp ) ) return;
  
  RefVariant rv;
  while ( iterate( &rv ) )
    {
      pout << rv;
      if ( with_verbose ) 
	{
	  std::vector<std::string> tok = Helper::char_split( rv.value() , ';' );
	  for ( int i = 0 ; i < tok.size() ; i ++ ) 
	    pout << "\t" << tok[i] << "\n";
	  pout << "\n";
	}
      else if ( with_meta ) pout << "\t" << rv.value();
      
      pout << "\n";
    }
  return;
}


void RefDBase::check_version()
{
  
  if ( ! sql.table_exists( "dbmeta" ) )
    Helper::halt( "old database format, expecting REFDB v"
                  + Helper::int2str( PLINKSeq::REFDB_VERSION_NUMBER() )
                  + " : to fix, remake this REFDB" );
  
  // expected version # is given by  PLINKSeq::REFDB_VERSION_NUMBER()                                                                      
  int v = 0;
  
  sqlite3_stmt * s = sql.prepare( "SELECT varvalue FROM dbmeta WHERE varname == 'VERSION'; " );
  
  if ( sql.step(s) )
    {
      v = sql.get_int( s , 0 );
      sql.finalise(s);
    }
  else // implies a new database, as version not yet set -- so add one                                                                    
    {
      sql.finalise(s);
      sqlite3_stmt * si = sql.prepare("INSERT OR REPLACE INTO dbmeta(varname, varvalue ) values( :x , :y ) ; " );
      std::string vn = "VERSION";
      v = PLINKSeq::REFDB_VERSION_NUMBER();
      sql.bind_text( si , ":x" , vn );
      sql.bind_int( si , ":y" , v );
      sql.step(si);
      sql.finalise(si);
    }

  if ( v != PLINKSeq::REFDB_VERSION_NUMBER() )
    Helper::halt("REFDB version "
                 + Helper::int2str( v ) + " but expected "
                 + Helper::int2str( PLINKSeq::REFDB_VERSION_NUMBER() )
                 + " : to fix, remake this REFDB" );

  return;

}

