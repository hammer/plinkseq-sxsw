#ifndef __PSEQ_NETDB_H__
#define __PSEQ_NETDB_H__

#include <string>
#include <map>
#include <set>
#include <vector>
#include <ostream>

#include "sqlwrap.h"
#include "regions.h"

class LocDBase;

class NetDBase {
  
 public:

  /// create a database
  bool attach( const std::string & name );

  bool attached() { return sql.is_open(); }
  
  /// Set LOCDB and group
  void set_locdb( LocDBase * l , const std::string & g );
  
  /// Get all score for a particular variant that meet 'threshold' score, up to 'depth' neighbours away
  std::set<Region> connections_regions( const std::string & seed , int depth = 1 , double thresold = 0 );
  std::set<std::string> connections( const std::string & seed , int depth = 1 , double thresold = 0 );
  std::set<int> connections( const std::string & seed , const std::map<std::string,int> &  genemap, int depth = 1 , double thresold = 0 );
  
  /// Load gene pairs into a table
  void load( const std::string & filename );
  
  bool dettach();
  
 private:

  SQL sql;
  
  sqlite3_stmt * stmt_insert_node;
  sqlite3_stmt * stmt_fetch_node;
  sqlite3_stmt * stmt_insert_edge;
  sqlite3_stmt * stmt_fetch_connections;
  
  std::map<std::string,int> node_id_map;

  int add_node( const std::string & node_name );
  int node_id( const std::string & node_name );
  
  void index();
  void drop_index();
  
  LocDBase * locdb;
  int grp;

};

#endif
