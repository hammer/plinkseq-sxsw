#ifndef __PSEQ_POLYPHEN2_H__
#define __PSEQ_POLYPHEN2_H__

#include <string>
#include <map>
#include <set>
#include <vector>
#include <ostream>

#include "sqlwrap.h"

class LocDBase;
class Variant;

struct PPH2Position {  

  PPH2Position() 
  {
    reference = alternate = "";
    score = 0.0;
    prediction = 0;
  }

  std::string reference;    // amino-acid
  std::string alternate;
  double score;   // pph2_prob
  int prediction;       // prediction : 0,1,2 = benign, possibily, probably damaging

  bool operator<(const PPH2Position & rhs ) const
  {
    if ( reference < rhs.reference ) return true;
    if ( reference > rhs.reference ) return false;
    return alternate < rhs.alternate;      
  }

  friend std::ostream & operator<< ( std::ostream & out , const PPH2Position & s ) 
  {
    out << s.reference << ">" << s.alternate << "\t" << s.score << "\t" << s.prediction;
    return out;
  }

};

struct PPH2Set {  

  PPH2Set()
  {
    protein_name = transcript_name = "";        
    scores.clear();
  }

  std::string protein_name;         // Q01234
  std::string transcript_name;      // NM_012345

  std::map<int,std::map<std::string,PPH2Position> > scores; // per AA position score   

  int size() const { return scores.size(); }
  int max_position() const 
  {
    std::map<int,std::map<std::string,PPH2Position> >::const_iterator i = scores.begin();
    int mx = 0;
    while ( i != scores.end() ) { if ( i->first > mx ) mx = i->first; ++i; }
    return mx;
  }

  const PPH2Position * position( const int p , const std::string & aa1, const std::string & aa2) const 
  {
    std::map<int,std::map<std::string,PPH2Position> >::const_iterator i = scores.find( p );
    if ( i == scores.end() ) return NULL;
    std::map<std::string,PPH2Position>::const_iterator j = i->second.find( aa1 + aa2 );
    return j == i->second.end() ? NULL : &(j->second) ;
  }
  
  
  void reset() 
  {
    transcript_name = protein_name = "";
    scores.clear();
  }
  
  friend std::ostream & operator<< ( std::ostream & out , const PPH2Set & s ) 
  {
    out << "PPH2: " << s.protein_name << " " << s.transcript_name << "\n";
    std::map<int,std::map<std::string,PPH2Position> >::const_iterator i = s.scores.begin();
    while ( i != s.scores.end() )
      {
	std::map<std::string,PPH2Position>::const_iterator j = i->second.begin();
	while ( j != i->second.end() )
	  {
	    out << " " << i->first << "\t"
		<< j->first << "\t"
		<< j->second << "\n";
	    ++j;
	  }
	++i;
      }
    return out;
  }

  
};


class PPH2DBase {
  
 public:

  /// create a database
  bool attach( const std::string & name );
  
  /// Set LOCDB for this PPH2 table
  void set_locdb( LocDBase * l ) { locdb = l; }

  /// Get score for a particular variant  
  bool score( const Variant & v , double & score , int & prediction );
  
  /// Load a PPH2 table into the DB
  void load( const std::string & filename );

  /// Does the PPH2 DB contain a given gene?
  bool present( const std::string & gene );

 private:
  
   static const std::string transcript_set_name;

  // Map transcript name (RefSeq) to PPH2 scores
  std::map<std::string, PPH2Set> cache;
  
  // Main dataastore

  SQL sql;
  
  sqlite3_stmt * stmt_insert_pph2_scores;
  sqlite3_stmt * stmt_fetch_pph2_scores;
  
  sqlite3_stmt * stmt_insert_pph2_id;
  sqlite3_stmt * stmt_fetch_pph2_id;

  //
  // Functions
  //
  
  void init();
  void release();
  void index();
  void drop_index();

  // Accumulate scores for a gene
  void accumulate( PPH2Set & , const std::vector<std::string> & tok );
  
  // Add to database
  void insert( const PPH2Set & );
  
  // Lookup (from cache first) based on transcript name (RefSeq)
  PPH2Set * lookup( const std::string & );
  
  // Transcript database
  LocDBase * locdb;

};

#endif
