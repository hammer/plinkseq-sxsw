#ifndef __IBD_SEGS_H__
#define __IBD_SEGS_H__

#include <string>
#include <vector>
#include <map>

#include "pseq.h"

class Mask;

namespace Pseq {

  namespace IBD {
    
    void test_wrapper( const std::string & segment_list , const std::string & gwas_phenotypes , int, Mask & m );
    
    struct IBDPartner {
      IBDPartner(std::string id, int affected ) 
	: id(id) , affected(affected) { } 
      std::string id;
      int affected;
    };
    
    
    class IBDSegmentHandler {
      
    public:
      
      IBDSegmentHandler( const std::string & db );
      ~IBDSegmentHandler();
      
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
      
    private:
      
      // Store segments in a SQLite database
      
      SQL sql;
      
      sqlite3_stmt * stmt_insert;
      sqlite3_stmt * stmt_fetch;
      
    };
    

    struct Aux {
      Aux() { rseed = 0; g = NULL; }
      GStore * g;
      IBDSegmentHandler * ibd;
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
