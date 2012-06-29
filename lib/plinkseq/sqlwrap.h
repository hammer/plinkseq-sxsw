#ifndef __SQLWRAP_H__
#define __SQLWRAP_H__

#include <iostream>
#include <string>
#include <vector>
#include <map>

#include "plinkseq/helper.h"

#include "sqlite3.h"

class SQL {

 public:

    SQL() {
      db = NULL;
      name = "";
    }
    
  bool open(std::string n);
  void synchronous(bool);
  void close();
  bool is_open() const { return db; }
  bool query( const std::string & q);
  bool table_exists( const std::string & );

  sqlite3_stmt * prepare(const std::string & q);
  sqlite3_stmt * prepare(const std::string & q, const std::string & key);
  sqlite3_stmt * fetch_prepared(const std::string & key);

  bool step(sqlite3_stmt * stmt);
  void reset( sqlite3_stmt * stmt );
  void finalise(sqlite3_stmt * stmt);

  bool loadExtension(std::string);

  void begin();
  void commit();

  uint64_t last_insert_rowid()
    { return sqlite3_last_insert_rowid(db); }
  
  void bind_int( sqlite3_stmt * stmt , const std::string & index , int value );
  void bind_int64( sqlite3_stmt * stmt , const std::string & index , uint64_t value );
  void bind_double( sqlite3_stmt * stmt , const std::string & index , double value );
  void bind_text( sqlite3_stmt * stmt , const std::string & index , const std::string & value );
  void bind_blob( sqlite3_stmt * stmt , const std::string & index , blob & );
  void bind_null( sqlite3_stmt * stmt , const std::string & index );

  int get_int( sqlite3_stmt *, int );
  uint64_t get_int64( sqlite3_stmt *, int );
  double get_double( sqlite3_stmt *, int );
  std::string get_text( sqlite3_stmt *, int );
  blob get_blob( sqlite3_stmt *, int );

  int lookup_int(sqlite3_stmt *);
  int lookup_int(const std::string & q);
  uint64_t lookup_int64(sqlite3_stmt *);
  std::vector<int> intTable( const std::string & q, int cols);
  std::vector<int> intTable(sqlite3_stmt * stmt, int cols);

  std::vector<uint64_t> int64Table(const std::string & q, int cols);
  std::vector<uint64_t> int64Table(sqlite3_stmt * stmt, int cols);

  static std::string header_version() 
    {
      return sqlite3_libversion();
    }

  static std::string library_version() 
    { 
      return SQLITE_VERSION;
    }
  
  sqlite3 * pointer() { return db; }

 private:
  
  // Keep track of all prepared statements
  std::set<sqlite3_stmt*> qset;
  
  // Map of prepared statements
  std::map<std::string,sqlite3_stmt*> qmap;
  
  // Database
  sqlite3 * db;
  
  // Return code, error msg
  int rc;
  char * db_err;

  // Name of database
  std::string name;
 
};

#endif
