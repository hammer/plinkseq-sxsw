#ifndef __PSEQ_IBS_H__
#define __PSEQ_IBS_H__

#include <string>
#include <vector>

#include "pseq.h"
#include "util.h"

class Mask;


namespace Pseq {
  
  namespace IBS {

      void regargs( Pseq::Util::Options * args );    
      
      bool calculate( Mask & );

    class Aux {
      
    public:
      
      Aux( const int n ) 
	{ 
	  ibs_cnt.resize(n);
	  for (int i=0; i<n; i++)
	    ibs_cnt[i].resize(i,0);
	  
	  obs_cnt.resize(n);
	  for (int i=0; i<n; i++)
	    obs_cnt[i].resize(i,0);
	  
	  // ibs[0..n-1][1..n-1]
	  // ibs[i][j] where i>j
	}

      void ibs(const int i, const int j, const int s)
      {
	if ( i > j ) ibs_cnt[i][j] += s;
	else if ( j > i ) ibs_cnt[j][i] +=s;
      }
      
      void obs(const int i, const int j, const int s)
      {
	if ( i > j ) obs_cnt[i][j] += s;
	else if ( j > i ) obs_cnt[j][i] +=s;
      }
      
      int ibs( const int i , const int j ) 
      {
	if ( i > j ) return ibs_cnt[i][j];
	else if ( j > i ) return ibs_cnt[j][i];
	else return 1;
      }

      int obs( const int i , const int j ) 
      {
	if ( i > j ) return obs_cnt[i][j];
	else if ( j > i ) return obs_cnt[j][i];
	else return 1;
      }
      
    private:
      
      std::vector< std::vector<int> > ibs_cnt;

      std::vector< std::vector<int> > obs_cnt;      
      
    };
    
  }

}


/// Worker functions
    
void f_IBS_calculator( Variant & v , void * p );


#endif
