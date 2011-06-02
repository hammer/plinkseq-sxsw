#ifndef __PSEQ_NETASSOC_H__
#define __PSEQ_NETASSOC_H__

#include "func.h"

#include <map>
#include <string>
#include <set>

class NetDBase;

namespace Pseq {
  namespace Assoc {
    namespace NetDB {
      
      struct Aux_netdet { 
	Aux_netdet() { nvar = 0; } 
	std::map<int,double> ind;
	int nvar;
	std::string name;
      };
      
      struct Aux_net { 
	std::map<int,Aux_netdet> gscore;
      };
      

      struct Aux_connection { 
	int extension;
	int parent;
	int acnt;
	int ucnt;
	bool operator< ( const Aux_connection & rhs ) const { return extension < rhs.extension; }
      };

      struct next_node_t { 
      next_node_t( int e , int p = -1 , int d = 1 ) 
      : extension(e) , parent(p) , depth(d) , retain(true) { } 
	int extension;
	int parent;
	int depth;    
	bool retain;
	bool operator< ( const next_node_t & rhs ) const { return extension < rhs.extension; }
      };
      
      bool driver( const std::map<std::string,Aux_netdet> & gscore , const Pseq::Util::ArgMap & args, Mask & m );
      
      void net_test( const int seed , long int , NetDBase & netdb , 
		     const std::set<std::string> & testset , 
		     const std::map<std::string,int> &, 
		     const std::map<int,Aux_netdet> & gscore );
      
      // read/write gene-scores
      
      void write_gscores( const std::string & filename , const std::map<std::string,Aux_netdet> & gscore );
      std::map<std::string,Aux_netdet> read_scores( const std::string & filename );
      
      // actual genic test

      double net_statistic( const int seed, 
			    NetDBase * netdb , 
			    std::set<int> & connections , 
			    std::set<Aux_connection> * endset , 
			    const std::map<int, Aux_netdet > & gscore , 			    
			    const std::map<std::string,int> & genemap, 
			    int * nvar = NULL , double * pa = NULL , double * pu = NULL );
      
      double stat_adder( const int seed , 
			 std::set<int> & inset , 
			 const std::map<int,Aux_netdet > & s , 
			 int * nvar = NULL , double * pa = NULL , double * pu = NULL );
      
    }
    
  }
}

#endif
