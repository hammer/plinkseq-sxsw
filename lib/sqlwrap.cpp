
#include "sqlwrap.h"
#include "defs.h"

using namespace std;
using namespace Helper;

extern Log plog;
extern Log debug;

bool SQL::open(string n)
{
  name = n;
  
  rc = sqlite3_open( name.c_str() , &db);
  
  if ( rc ) Helper::halt("problem opening database: " + name );

  if ( PLINKSeq::SQLITE_SCRATCH_FOLDER() != "" )
    {
      query( "PRAGMA temp_store_directory = '" 
	     + PLINKSeq::SQLITE_SCRATCH_FOLDER() 
	     + "';");
    }
  return rc == 0;
}

void SQL::synchronous(bool b)
{
  if ( !b ) 
    query( "PRAGMA synchronous=OFF;" );
  else
    query( "PRAGMA synchronous=FULL;" );
}

bool SQL::table_exists( const std::string & table_name )
{
  sqlite3_stmt * s = prepare( "SELECT name FROM sqlite_master WHERE type='table' AND name= :table_name ; " );
  bind_text( s , ":table_name" , table_name );
  if ( step(s) ) 
    {
      finalise(s);
      return true;
    }
  finalise(s);
  return false;
}

void SQL::close()
{
  if ( db ) 
    {
      sqlite3_close(db);	    
      db = NULL;
    }
}

bool SQL::query( const std::string & q )
{  
  char * db_err;
  rc = sqlite3_exec( db , q.c_str() , 0 , 0 , &db_err );
  if ( rc ) plog.warn( "database warning: " + std::string(db_err) );
  return rc == 0;
}

sqlite3_stmt * SQL::prepare( const std::string & q )
{   
  sqlite3_stmt * p;
  int rc = sqlite3_prepare_v2( db , q.c_str() , q.size() , &p , NULL );   
  if ( rc ) plog.warn( "database warning preparing query " + std::string( sqlite3_errmsg(db) ) );
  else qset.insert(p);
  return rc ? NULL : p;
}

sqlite3_stmt * SQL::prepare( const std::string & q, const std::string & key )
{   
  sqlite3_stmt * p;
  int rc = sqlite3_prepare( db , q.c_str() , q.size() , &p , NULL );   
  if ( rc ) Helper::halt( db_err );
  else qset.insert( p );
  qmap.insert( make_pair( key , p ) );
  return rc ? NULL : p;
}

sqlite3_stmt * SQL::fetch_prepared( const std::string & key )
{
  std::map<std::string,sqlite3_stmt*>::iterator i = qmap.find(key);
  if ( i == qmap.end() ) return NULL;
  return i->second;
}

void SQL::begin()
{  
  char * db_err;
  string q = "BEGIN;";
  rc = sqlite3_exec( db , q.c_str() , 0 , 0 , &db_err );
  if ( rc ) Helper::halt( db_err );
}

void SQL::finalise(sqlite3_stmt * stmt)
{
  std::set<sqlite3_stmt*>::iterator i = qset.find( stmt );
  if ( stmt && i != qset.end() ) 
    {
      qset.erase( i ); 
      sqlite3_finalize( stmt );  
    }
  stmt = NULL;
}

bool SQL::step(sqlite3_stmt * stmt)
{
  
  rc = sqlite3_step( stmt );

  if ( rc != SQLITE_ROW && rc != SQLITE_DONE )
    {
      reset(stmt);
      Helper::halt( "database (" + name +") error (" + Helper::int2str( sqlite3_errcode(db) ) +") " + sqlite3_errmsg(db) );
    }
  
  return rc == SQLITE_ROW;
}

void SQL::reset( sqlite3_stmt * stmt )
{
  sqlite3_reset( stmt );
}

bool SQL::loadExtension(string libname)
{  
  Helper::halt( "sqlite load-extension not supported" );
  return false;

  // no longer supported
  //   sqlite3_enable_load_extension(db,1);
  //   rc = sqlite3_load_extension( db , libname.c_str() , 0 , &db_err );
  //   return rc == 0;
}

void SQL::bind_int( sqlite3_stmt * stmt , const std::string & index , int value )
{
  sqlite3_bind_int( stmt , 
		     sqlite3_bind_parameter_index( stmt , index.c_str() ) ,
		     value );
}

void SQL::bind_null( sqlite3_stmt * stmt , const std::string & index  )
{
    sqlite3_bind_null( stmt , 
		       sqlite3_bind_parameter_index( stmt , index.c_str() ) );
}

void SQL::bind_int64( sqlite3_stmt * stmt , const std::string & index , uint64_t value )
{
  sqlite3_bind_int64( stmt , 
		     sqlite3_bind_parameter_index( stmt , index.c_str() ) ,
		      value );
}

void SQL::bind_double( sqlite3_stmt * stmt , const std::string & index , double value )
{
  sqlite3_bind_double( stmt , 
		     sqlite3_bind_parameter_index( stmt , index.c_str() ) ,
		     value );
}

void SQL::bind_text( sqlite3_stmt * stmt , const std::string & index , const std::string & value )
{

  sqlite3_bind_text( stmt , 
		     sqlite3_bind_parameter_index( stmt , index.c_str() ) ,
		     value.c_str() , 
		     value.size() , 
		     0 );
}
  

void SQL::bind_blob( sqlite3_stmt * stmt , const std::string & index , blob & value )
{
    rc = sqlite3_bind_blob( stmt , 
			    sqlite3_bind_parameter_index( stmt , index.c_str() ) ,
			    value.p,
			    value.l,
			    0 );

}


int SQL::get_int( sqlite3_stmt * stmt , int idx )
{
  return sqlite3_column_int( stmt , idx );
}

uint64_t SQL::get_int64( sqlite3_stmt * stmt , int idx )
{
  return sqlite3_column_int64( stmt , idx );
}

double SQL::get_double( sqlite3_stmt * stmt , int idx )
{
  return sqlite3_column_double( stmt , idx );
}

string SQL::get_text(  sqlite3_stmt * stmt , int idx )
{
  const unsigned char * s = sqlite3_column_text( stmt , idx );
  if ( s == NULL )
    return "";
  else
    return (const char*)s;
}

blob SQL::get_blob( sqlite3_stmt * stmt , int idx )
{
    blob b;
    b.p = (const char*)sqlite3_column_blob( stmt , idx );
    b.l = sqlite3_column_bytes( stmt , idx );
    return b;
}

void SQL::commit()
{
  query( "COMMIT;" );
}

std::vector<int> SQL::intTable( const std::string & q, int cols)
{  
  return intTable( prepare(q) , cols );
}

std::vector<int> SQL::intTable(sqlite3_stmt * stmt, int cols)
{
  vector<int> res;
  rc = sqlite3_step( stmt );
  while ( rc == SQLITE_ROW )
    {
      for ( int i = 0 ; i < cols ; i++ )
	res.push_back ( sqlite3_column_int( stmt , i ) );           
      rc = sqlite3_step( stmt );
    }
  sqlite3_finalize(stmt);  
  return res;
}

std::vector<uint64_t> SQL::int64Table( const std::string & q, int cols)
{  
  return int64Table( prepare(q) , cols );
}

std::vector<uint64_t> SQL::int64Table(sqlite3_stmt * stmt, int cols)
{
  std::vector<uint64_t> res;
  rc = sqlite3_step( stmt );
  while ( rc == SQLITE_ROW )
    {
      for ( int i = 0 ; i < cols ; i++ )
	res.push_back ( sqlite3_column_int64( stmt , i ) );           
      rc = sqlite3_step( stmt );
    }
  sqlite3_finalize(stmt);  
  return res;
}

int SQL::lookup_int(sqlite3_stmt * stmt)
{
  int r = -1;
  rc = sqlite3_step( stmt );
  if ( rc == SQLITE_ROW )
    r = sqlite3_column_int( stmt , 0 );
  return r;
}

int SQL::lookup_int( const std::string & q )
{
  sqlite3_stmt * s = prepare(q);
  int r = -1;
  rc = sqlite3_step( s );
  if ( rc == SQLITE_ROW )
    r = sqlite3_column_int( s , 0 );
  finalise(s);
  return r;
}

uint64_t SQL::lookup_int64(sqlite3_stmt * stmt)
{
  uint64_t r = 0;
  rc = sqlite3_step( stmt );
  if ( rc == SQLITE_ROW )
    r = sqlite3_column_int64( stmt , 0 );  
  return r;
}

