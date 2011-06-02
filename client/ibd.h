#ifndef __IBD_SEGS_H__
#define __IBD_SEGS_H__

#include <string>
#include <vector>
#include <map>
#include <set>

#include "pseq.h"

class Mask;

namespace Pseq {

  namespace IBD {
    
    void test_wrapper( const std::string & segment_list , const std::string & gwas_phenotypes , int, Mask & m );
    
    void load_wrapper( const std::string & segment_list , const std::string & ibddb );
    
    void sharing_wrapper( const std::string & database , Mask & m );


    struct IBDPartner {
      IBDPartner(std::string id, int affected ) 
	: id(id) , affected(affected) { } 
      std::string id;
      int affected;
    };
    
    
    class IBDDBase {
      
    public:
      
      IBDDBase( const std::string & db );
      ~IBDDBase();
      
      // Load all pairwise segments      
      void load( const std::string & db );
      
      // For a proband, get IDs of all other individuals who 
      // share at least one segment at this position
      std::vector<IBDPartner> fetch( const std::string & id , const Region & r);
      
      int2 case_control_count( const std::string & id , const Region & r );
      int2 case_control_count( const std::string & id , const Region & r , 
			       std::map<std::string,int> & imap, 
			       std::vector<int> & pmap,
			       std::vector<int> & permed );

      // Given a pair of individuals, return all segments they share 
      std::set<Region> shared_for_pair( const std::string & id1 , const std::string & id2 );
      
    private:
      
      // Store segments in a SQLite database
      
      SQL sql;
      
      sqlite3_stmt * stmt_insert;
      sqlite3_stmt * stmt_fetch;
      sqlite3_stmt * stmt_fetch_pair;
      
    };
    

    struct Aux {
      Aux() { rseed = 0; g = NULL; }
      GStore * g;
      IBDDBase * ibd;
      long int rseed;
      int minm;
      int maxm;
      std::map<std::string,int> * imap; 
      std::vector<int> * pmap;
    };
    
  }

}
    
void g_STEST_association( VariantGroup & vars , void * p );

#endif
