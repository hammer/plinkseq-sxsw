#ifndef __ITER_REPORT_H__
#define __ITER_REPORT_H__

#include "helper.h"

#include <string>

struct IterationReport {

  IterationReport( bool valid = true , bool grouping = false , int v_limit = 0 ) 
    : valid(valid) , grouping(grouping) , v_limit( v_limit ) , 
       v_consid(0), v_accept(0), g_count(0)
    { } 

  std::string report() 
  {

    if ( ! valid ) 
      return "Invalid specification: no variants or groups processed";

    if ( grouping ) 
      return "Processed " + Helper::int2str( v_accept ) 
	+ " of " + Helper::int2str( v_consid ) 
	+ " variants considered, in " + Helper::int2str( g_count ) 
	+ " distinct groups";

    return "Processed " + Helper::int2str( v_accept ) 
      + " of " + Helper::int2str( v_consid ) + " variants considered";
  }

  bool reached_limit()
  {
    return v_limit != 0 && v_accept >= v_limit ;
  }

  bool accepted_variant() 
  {
    ++v_consid;
    ++v_accept;
    return v_limit == 0 || v_accept < v_limit ;
  }
  
  void rejected_variant() 
  {
    ++v_consid;
  }
  
  void processed_group() 
  { 
    ++g_count; 
  }

  int processed() 
  {
    return v_accept;
  }

  private:
  
  bool valid;
  bool grouping;
  int v_limit;
  int v_consid;
  int v_accept;
  int g_count;

};

#endif
