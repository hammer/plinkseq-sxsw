#include "netdb.h"
#include "locdb.h"
#include "helper.h"

bool NetDBase::attach( const std::string & name )
{
  
  if ( name == "-" ) { dettach(); return false; } 
  
  //if ( Helper::fileExists( name ) ) { sql.open( name ); return true; }
  
  if ( attached() ) dettach();  

  sql.open( name );   

  // Create database

  sql.synchronous(false);
  
  // Edge table ( double-enter n1-n2 and n2-n1 )
  
  sql.query(" CREATE TABLE IF NOT EXISTS edges("
	    "   node1_id   INTEGER , "
	    "   node2_id   INTEGER , "	      
	    "   score      REAL ) ; " );
  
  sql.query( "CREATE TABLE IF NOT EXISTS nodes("
	     "  node_id    INTEGER PRIMARY KEY , " 
	     "  name       VARCHAR(12) ); " );
  
  // Indices
  index();
  
  stmt_insert_node = sql.prepare( " INSERT OR REPLACE INTO nodes( name ) values( :name ); " ); 
  stmt_fetch_node  = sql.prepare( " SELECT node_id FROM nodes WHERE name == :name ; " );
  stmt_insert_edge = sql.prepare( " INSERT OR REPLACE INTO edges( node1_id , node2_id , score ) values( :n1 , :n2 , :score ); " );
  stmt_fetch_connections = sql.prepare( " SELECT name FROM nodes WHERE node_id IN ( SELECT node2_id FROM edges WHERE node1_id == :n ) ; " );

  return true;

}


void NetDBase::set_locdb( LocDBase * l , const std::string & g )
{ 
  locdb = l; 
  grp = 0;
  if ( locdb && locdb->attached() ) 
    grp = locdb->lookup_group_id ( g ) ;
}


int NetDBase::node_id( const std::string & node_name )
{

  // stored in cache?
  std::map<std::string,int>::iterator i = node_id_map.find( node_name );
  if ( i != node_id_map.end() ) return i->second; 
  
  sql.bind_text( stmt_fetch_node , ":name" , node_name );
  int nid = 0;
  if ( sql.step( stmt_fetch_node ) )
    nid = sql.get_int( stmt_fetch_node , 0 );    
  sql.reset( stmt_fetch_node );

  return nid;
}

int NetDBase::add_node( const std::string & node_name )
{
  // **assumes** does not already exist
  sql.bind_text( stmt_insert_node , ":name" , node_name );
  sql.step( stmt_insert_node );
  sql.reset( stmt_insert_node );
  int nid = sql.last_insert_rowid();
  node_id_map[ node_name ] = nid;
  return nid;
}

bool NetDBase::dettach()
{
  sql.finalise( stmt_insert_node );
  sql.finalise( stmt_fetch_node );
  sql.finalise( stmt_insert_edge );
  sql.finalise( stmt_fetch_connections );
  sql.close();
}

void NetDBase::drop_index()
{
  sql.query( "DROP INDEX IF EXISTS nameIndex; " );
  sql.query( "DROP INDEX IF EXISTS nodeIndex; " );
}

void NetDBase::index()
{
  sql.query( "CREATE INDEX IF NOT EXISTS nameIndex ON nodes( name ); " );
  sql.query( "CREATE INDEX IF NOT EXISTS nodeIndex ON edges( node1_id ); " );
}

void NetDBase::load( const std::string & filename )
{

  Helper::checkFileExists( filename );

  InFile F1( filename );

  // expect format g1 , g2 , score 
  // tab-delimited

  drop_index();

  sql.begin();
  
  int edge_cnt = 0 , node_cnt = 0;
  
  while ( ! F1.eof() )
    {
      std::string n1 , n2;
      double sc;
      F1 >> n1 >> n2 >> sc;
      if ( n1 == "" ) break;
      
      int nid1 = node_id( n1 );
      if ( nid1 == 0 ) { nid1 = add_node( n1 ); ++node_cnt; }

      int nid2 = node_id( n2 );
      if ( nid2 == 0 ) { nid2 = add_node( n2 ); ++node_cnt; }
      
      sql.bind_int( stmt_insert_edge , ":n1" , nid1 );
      sql.bind_int( stmt_insert_edge , ":n2" , nid2 );
      sql.bind_double( stmt_insert_edge , ":score" , sc );
      sql.step( stmt_insert_edge );
      sql.reset( stmt_insert_edge );
      
      sql.bind_int( stmt_insert_edge , ":n1" , nid2 );
      sql.bind_int( stmt_insert_edge , ":n2" , nid1 );
      sql.bind_double( stmt_insert_edge , ":score" , sc );
      sql.step( stmt_insert_edge );
      sql.reset( stmt_insert_edge );

      ++edge_cnt;

      if ( edge_cnt % 1000 == 0 ) 
	plog << edge_cnt << " edges\t" << node_id_map.size() << " nodes \n";
      
    }

  plog << "added " << node_cnt << " " << node_id_map.size() << " unique nodes, " << edge_cnt << " edges\n";

  sql.commit();
  index();

}

  
std::set<Region> NetDBase::connections_regions( const std::string & seed , int depth , double thresold  )
{    
  std::set<Region> r;
  int nid = node_id( seed );
  if ( nid == 0 || ! locdb ) return r;  
  sql.bind_int( stmt_fetch_connections , ":n" , nid );
  while ( sql.step( stmt_fetch_connections ) )
    r.insert( locdb->get_region( grp , sql.get_text( stmt_fetch_connections , 0 ) ) );
  sql.reset( stmt_fetch_connections );
  return r;
}

std::set<std::string> NetDBase::connections( const std::string & seed , int depth , double thresold  )
{  
  std::set<std::string> r;
  int nid = node_id( seed );
  if ( nid == 0 ) return r;  
  sql.bind_int( stmt_fetch_connections , ":n" , nid );
  while ( sql.step( stmt_fetch_connections ) )
    r.insert( sql.get_text( stmt_fetch_connections , 0 ) );
  sql.reset( stmt_fetch_connections );
  return r;
}


std::set<int> NetDBase::connections( const std::string & seed , 
				     const std::map<std::string,int> & genemap, 
				     int depth , double thresold  )
{  
  std::set<int> r;
  int nid = node_id( seed );
  if ( nid == 0 ) return r;  
  sql.bind_int( stmt_fetch_connections , ":n" , nid );
  while ( sql.step( stmt_fetch_connections ) )
    {
      std::string g = sql.get_text( stmt_fetch_connections , 0 );
      std::map<std::string,int>::const_iterator i = genemap.find( g );
      if ( i != genemap.end() ) r.insert( i->second );
    }
  sql.reset( stmt_fetch_connections );
  return r;
}
