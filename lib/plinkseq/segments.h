#ifndef __SEGMENTS_H__
#define __SEGMENTS_H__


#include "plinkseq/gstore.h"
#include "plinkseq/helper.h"
#include "plinkseq/regions.h"

#include <map>
#include <set>


class OverlapResult {
  
public:
  
  OverlapResult()
    {
      totalLength = 0;
      exonLength = 0;
      nExons = 0;
      nTargets = 0;    
    }
    
  int totalLength;
  int exonLength; 
  int nExons;
  int nTargets;
  
  std::map<int,std::set<int2> > cover;

  // IDs of all overlapping Regions
  std::set<uint64_t> overlapping_target_region_ids;
};


class OverlapResults {
    
 public:    

  std::map<Region,OverlapResult> result;

  OverlapResults(uint64_t t) 
    { target_id = t; }    
  
  uint64_t target_id ;    
  
  void load_regions(std::set<Region> r);
};


#endif
