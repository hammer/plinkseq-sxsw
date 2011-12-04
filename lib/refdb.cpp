
#include "refdb.h"
#include "variant.h"
#include "gstore.h"
#include "filemap.h"

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
  
  sql.query(" CREATE TABLE IF NOT EXISTS refvariants("
            "   group_id  INTEGER NOT NULL , "
            "   name      VARCHAR(40) , "
            "   chr       INTEGER NOT NULL , "
            "   bp1       INTEGER NOT NULL , "
            "   bp2       INTEGER NOT NULL , "
	    "   value     VARCHAR(10) ); " );

  
  // Meta-information table
  
  sql.query( " CREATE TABLE IF NOT EXISTS metatypes("
	     "   field_id     INTEGER PRIMARY KEY , "
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

  init();
  
  // set meta-type information

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
  release();
  init();
}


void RefDBase::index()
{
  if ( ! attached() ) return;  
  sql.query( "CREATE INDEX IF NOT EXISTS ind1 ON refvariants(group_id,chr, bp1); " );    
  release();
  init();
  
}

bool RefDBase::init()
{
    
  //
  // Core
  //

  stmt_lookup = 
    sql.prepare( "SELECT * FROM refvariants WHERE chr == :chr AND bp1 == :bp1 AND group_id == :group_id ; "  );

  stmt_lookup_range = 
    sql.prepare( "SELECT * FROM refvariants "
		 " WHERE chr == :chr AND bp1 >= :rstart AND group_id == :group_id AND bp1 <= :rend ; "  );

  stmt_lookup_range_count = 
    sql.prepare( "SELECT COUNT(*) FROM refvariants "
		 " WHERE group_id == :group_id AND chr == :chr AND :bp2 >= bp1 AND :bp1 <= bp2 ; "  );

  stmt_lookup_dbsnp = 
    sql.prepare( "SELECT * FROM refvariants WHERE group_id == :group_id AND name == :name ; "  );

  stmt_insert =
    sql.prepare( "INSERT OR REPLACE INTO refvariants ( name, group_id, chr, bp1, bp2, value ) "
		 " values( :name , :group_id, :chr, :bp1, :bp2, :value ) ; " );

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
    sql.prepare(" INSERT INTO metatypes (name,type,number,description) values( :name, :type, :number, :description ); " );

  stmt_fetch_metatypes = 
      sql.prepare(" SELECT name , type , number, description "
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
				int num, 
				int mgrp, 
				const std::string & description )
{
  
  
  sql.bind_text( stmt_meta_insert_prep , ":name" , name );
  
  if ( sql.step( stmt_meta_insert_prep ) )
    {
      // store ID for this meta-tag if already seen
      mtmap[ name ] = sql.get_int64( stmt_meta_insert_prep , 0 );
    }
  else
    {	      
      // otherwise, need to store this in DB
      
      sql.bind_text( stmt_meta_insert_prep2 , ":name" , name );
      sql.bind_int( stmt_meta_insert_prep2 , ":type" , mt );
      sql.bind_int( stmt_meta_insert_prep2 , ":number" , num );
      sql.bind_text( stmt_meta_insert_prep2 , ":description" , description );
      sql.step( stmt_meta_insert_prep2 );

      // store the newly-created field ID locally
      mtmap[ name ] = sql.last_insert_rowid();

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
	plog.counter( "parsed " + Helper::int2str( inserted ) + " rows" );
    }
  plog.counter("\n");
    
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
    
    int bp2 = rv.stop() == rv.start() ? 0 : rv.stop() ; 
    sql.bind_int( stmt_insert , ":bp2" ,     bp2 );

    sql.bind_text( stmt_insert , ":value" ,  rv.value() );
    sql.step( stmt_insert );
    uint64_t id = sql.last_insert_rowid();
    sql.reset( stmt_insert );
            
}


std::map<std::string,mType> RefDBase::populate_metatypes( std::map<std::string,int> * meta )
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
  
  std::map<std::string,mType> mt = populate_metatypes( meta );
  
  
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
  std::cout << " binding = " << g << "\n";
  sql.bind_int( stmt_dump , ":group_id" , g );
  return true;
}

void RefDBase::construct_inplace( sqlite3_stmt * s , RefVariant * rv )
{
  rv->group(      sql.get_int( s , 0 ) );	    
  rv->name(       sql.get_text( s , 1 ) );
  rv->chromosome( sql.get_int( s , 2 ) ) ;
  rv->start(      sql.get_int( s , 3 ) ) ;
  int bp2 =       sql.get_int( s , 4 );
  if ( bp2 ) rv->stop( bp2 ? bp2 : rv->start() );
  rv->value(      sql.get_text( s , 5 ) );  
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


RefVariant RefDBase::lookup( const Variant & v , const std::string & grp_name )
{
    if ( ! attached() ) return RefVariant();
    return lookup( v, lookup_group_id( grp_name ) );
}

RefVariant RefDBase::lookup( const Variant & v , int grp_id )
{
  // Create variant -- by default, has 'observed'
  // flag set to false
  
  RefVariant r;
  
  // Search given chromosome position
  
  sql.bind_int( stmt_lookup , ":chr" , v.chromosome() );
  sql.bind_int( stmt_lookup , ":bp1" , v.position() );
  sql.bind_int( stmt_lookup , ":group_id" , grp_id );
  
  if ( sql.step( stmt_lookup ) )
    r = construct( stmt_lookup );
  
  sql.reset( stmt_lookup );
  
  return r;

  // abandon lookup by name
//   if ( r.observed() ) 
//     return r;
  
//   // Otherwise, search given dbSNP number
//   if ( v.name() != "0" && v.name() != "" && v.name() != "." )
//     {
//       sql.bind_text( stmt_lookup_dbsnp , ":snp" , v.name() );
//       sql.bind_int( stmt_lookup_dbsnp , ":group_id" , grp_id );
//       if ( sql.step( stmt_lookup_dbsnp ) )
// 	r = construct( stmt_lookup_dbsnp );
//       sql.reset( stmt_lookup_dbsnp );
//       if ( r.observed() ) 
// 	return r;	
//     }
//   // Otherwise, return a null result    
//   return r;  
}


std::set<RefVariant> RefDBase::lookup( const Region & region , const int grp_id )
{
  
  // Search given chromosome position
  
  sql.bind_int( stmt_lookup_range , ":chr" , region.chromosome() );
  sql.bind_int( stmt_lookup_range , ":rstart" , region.start.position() );
  sql.bind_int( stmt_lookup_range , ":rend" , region.stop.position() );
  sql.bind_int( stmt_lookup_range , ":group_id" , grp_id );
  
  std::set<RefVariant> s;
  while ( sql.step( stmt_lookup_range ) )
    {      
      RefVariant r = construct( stmt_lookup_range );      
      if ( r.observed() ) s.insert( r );      
    }
  sql.reset( stmt_lookup_range );
  return s;  
}

std::set<RefVariant> RefDBase::lookup( const Region & region , const std::string & grp_name )
{
  std::set<RefVariant> empty;
  if ( ! attached() ) return empty;
  return lookup( region, lookup_group_id( grp_name ) );
}



bool RefDBase::annotate( Variant & v , const int grp_id )
{
  
  RefVariant r = lookup( v , grp_id );
  
  if ( ! r.observed() ) return false;

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
      r.meta.parse( r.value() );
      v.meta.append( r.meta , gname );
    }
  
  return true;
}


bool RefDBase::annotate( Variant & v , const std::string & name )
{
  int id = lookup_group_id( name );
  if ( id == 0 ) return false;
  return annotate(v,id);  
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
  std::string value = sql.get_text( s , 5 );
  
  return RefVariant(grp_id,name,chr,bp1,bp2,value);

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
    

    if ( ! ugly ) ss << "---Reference DB summary---\n";
    
    while ( sql.step(s) )
      {
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

    return ss.str();
}

void RefDBase::dump( const std::string & grp , bool with_meta )
{

  if ( ! init_iterate( grp ) ) return;
  
  RefVariant rv;
  while ( iterate( &rv ) )
    {
      plog << rv;
      if ( with_meta ) plog << "\t" << rv.value();
      plog << "\n";
    }
  return;
}

