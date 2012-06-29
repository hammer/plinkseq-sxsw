#ifndef __INDDB_H__
#define __INDDB_H__

#include "sqlwrap.h"
#include "plinkseq/individual.h"

#include <string>
#include <vector>

class IndDBase {
  
 public:
  
  IndDBase()
    {
    }
  
  ~IndDBase()
    {
      dettach();
    }
  
  void wipe(const std::string & name);
  bool new_db( const std::string & name);
  bool attach( const std::string & name);
  bool dettach();
  bool init();
  bool release();
  bool index();
  bool drop_index();

  bool attached() { return sql.is_open(); }
    
  std::string summary( bool );  
  
  // 
  // Insertions
  // 
  
  uint64_t insert_phenotype( const std::string & , const std::string & , const std::string &, const std::string & ); 
  uint64_t insert( const Individual & , bool * newone = NULL );
  void insert( const uint64_t i, const uint64_t p , const int );
  void insert( const uint64_t i, const uint64_t p , const double );
  void insert( const uint64_t i, const uint64_t p , const std::string & );
  
  bool load_phenotypes( const std::string & filename );
  bool load_ped_info( const std::string & filename );
  
    
  //
  // Queries
  // 

  std::vector<Individual> fetch( );

  std::vector<Individual> fetch(const std::vector<std::string> & ids )
    {
      std::vector<Individual> ind;
      for ( size_t i = 0 ; i < ids.size(); i++ )
	{
	  int code = fetch_id( ids[i] );
	  if ( code == 0 ) continue;
	  ind.push_back( fetch( code ) );
	}
      return ind;
    }
  
  // append phenotypic information to an existing person (match
  // based on ID)
  bool fetch( Individual * person );

  // get sex of a particular individual
  sType sex( const std::string & id );

  // create a new person
  Individual fetch( uint64_t ind_id );
  
  // actually does the fetching
  bool fetch( Individual * , uint64_t ind_id );
  
  Individual fetch( const std::string & name ) 
  { 
    return fetch( fetch_id(name) );
  }

  int fetch_id( const std::string & name );
  
  int fetch_pheno_id( const std::string & name);
  
  void set_metatypes(bool clear = false );
  
  std::map<std::string, std::vector<std::string> > fetch_phenotype_info();
    
  void load_meta(std::vector<Individual> &, const std::string & ) ;
  
  int size();

  //
  // Individuals
  //
  
  bool replace_individual_id( const std::string & old_id , const std::string & new_id );

  //
  // Helpers
  //
  
  void begin() { sql.begin(); }
  void commit() { sql.commit(); }
  
  //
  // Query functions 
  //
  
      
 private:
  
  SQL sql;

  sqlite3_stmt * stmt_insert_metatype;
  sqlite3_stmt * stmt_fetch_metatypes;
  sqlite3_stmt * stmt_fetch_sex;
  sqlite3_stmt * stmt_fetch_phenotype_list;
  sqlite3_stmt * stmt_fetch_phenotype_values;
  
  sqlite3_stmt * stmt_insert_individual;
  sqlite3_stmt * stmt_update_individual;
  sqlite3_stmt * stmt_replace_individual_id;

  sqlite3_stmt * stmt_lookup_id;
  sqlite3_stmt * stmt_lookup_pheno_id;

  sqlite3_stmt * stmt_fetch_individual;
  sqlite3_stmt * stmt_fetch_individuals;
  
  sqlite3_stmt * stmt_insert_phenotype;
  sqlite3_stmt * stmt_insert_phenotypes;
  sqlite3_stmt * stmt_insert_metaphenotype;

  
  // Cache
  
  std::map<uint64_t,Individual> cache;
  
  
};

#endif
