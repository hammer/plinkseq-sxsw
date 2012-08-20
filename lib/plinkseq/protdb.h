#ifndef __PSEQ_PROTDB_H__
#define __PSEQ_PROTDB_H__

#include <string>
#include <map>
#include <set>
#include <vector>
#include <ostream>

#include "sqlwrap.h"
#include "meta.h"

class LocDBase;
class Variant;
class Region;

struct Feature { 

  std::string protein_id;
  std::string source_id;
  std::string feature_id;
  std::string feature_name;
  
  // 1-based, AA space
  int start;
  int stop;

  // Various meta-information
  std::string mstr;
  MetaInformation<MiscMeta> meta;  

  bool operator<( const Feature & rhs ) const 
  {
    if ( start < rhs.start ) return true;
    if ( start > rhs.start ) return false;
    
    if ( stop < rhs.stop ) return true;
    if ( stop > rhs.stop ) return false;

    if ( source_id < rhs.source_id ) return true;
    if ( source_id > rhs.source_id ) return false;

    return feature_id < rhs.feature_id;
  }
  
  friend std::ostream & operator<<( std::ostream & out, const Feature & f)
    {
      out << f.protein_id << "\t"
	  << f.start << ".." << f.stop << "\t"
	  << f.source_id << "::" << f.feature_id << "\t"
	  << f.feature_name << "\t"	  
	  << f.mstr ;
      return out;
    }
  

};


struct ProtFeatureSet {
  
  ProtFeatureSet()
  {    
    feat.clear();
  }
  
  void add( const std::string & t , const Feature & f ) 
  {
    feat[ t ].insert( f );
  }
  
  std::set<Feature> get( const std::string & t )
  {
    return feat[ t ];
  }

  void clear() 
  {
    feat.clear();
  }

  const std::set<Feature> * features( const std::string & t ) const 
  {
    std::map<std::string,std::set<Feature> >::const_iterator ii = feat.find( t );
    if ( ii == feat.end() ) return NULL;
    return &(ii->second);
  }
  

  // actual store
  std::map<std::string,std::set<Feature> > feat;

};



class ProtDBase {
  
 public:
  
  ProtDBase() { }
  
  /// create a database
  bool attach( const std::string & name );

  bool attached() { return sql.is_open(); }

  /// Load a protein-annotation table into the DB
  void load( const std::string & filename );
  
  /// Map a PROTDB source/feature to a LOCDB (i.e. to make a genomic region list)
  int map_to_genomic( LocDBase * locdb , const std::string & , const std::string & , const std::string & feature = "." , std::map<std::string,Region> * plmap = NULL );

  /// Does this variant have this particular feature or not?
  bool has_feature( const Variant & v , const std::string & transcript , const std::string & group , const std::string & annot );
  
  /// Return all features for a transcript, whether a variant exists there or not
  std::set<Feature> fetch( const std::string & transcript );
    
  
  
 private:
  
  
  //
  // Main dataastore
  //
  
  SQL sql;
  
  // insert main data row 
  sqlite3_stmt * stmt_insert;

  // search
  sqlite3_stmt * stmt_fetch_given_transcript;
  sqlite3_stmt * stmt_fetch_given_transcript_and_annot;

  sqlite3_stmt * stmt_fetch_given_annot;
  sqlite3_stmt * stmt_fetch_given_annot_feature;

  sqlite3_stmt * stmt_insert_mapping;
  sqlite3_stmt * stmt_fetch_mapping;

  //
  // Functions
  //
  
  void init();
  void release();
  void index();
  void drop_index();

  /// Add to database
  void insert( const ProtFeatureSet & );
  
  /// Lookup (from cache first) based on transcript name (RefSeq)
  ProtFeatureSet lookup( const std::string & );
  
};

#endif
