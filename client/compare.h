#ifndef __PSEQ_COMPARE_H__
#define __PSEQ_COMPARE_H__

#include <map>
#include <set>

class Mask;
class Variant;
class GStore;

void f_group_comp( Variant & v , void * p);

namespace Pseq 
{  
  namespace Assoc
    {
            
      bool group_comparison( Mask & );
      
      struct AuxGroupComp {  
	AuxGroupComp( GStore * p) { g = p; } 
	GStore * g;
	std::map<std::set<int>,int> counts;
      };
      
    }
}

#endif
  
