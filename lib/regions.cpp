#include "regions.h"

using namespace std;
using namespace Helper;

Region::Region(const RefVariant & rv)
{
  construct(0,rv.chromosome(), rv.start(), rv.stop(),rv.name(),"",0);
}

bool Subregion::overlaps(const Region & b) const 
{  
  return stop >= b.start && start <= b.stop; 
} 



bool Region::within( std::set<Region> & s )
{
    
  // Find lower and upper bound
  
  set<Region>::iterator lwr, uppr;
  
  lwr = s.lower_bound( *this );
  
  if ( lwr == s.end() ) return false;
  if ( overlaps( *lwr ) ) return true;
  if ( lwr == s.begin() ) return false;
  --lwr;
  
  return overlaps( *lwr );

}


Region::Region(const std::string & s, bool & flag)
{

  id = 0;
  group = 0;
  name = "";
  altname = "";

  // Assume in format chr1:100..200
  //  or chr1:100
  //  or chr1
  //  or chr1:100:rs1001 

  // *or* whitespace delimited, e.g. chr1 100 200
    
  flag = false;

  size_t p = s.find(":");
  
  if ( p == string::npos ) 
    {

      // is this a single chromosome?
      // set a large value (appropriate for human
      // chromosomes) to get whole thing...

      if ( Helper::chr_known( s ) ) 
	{
	  int chr = Helper::chrCode( s );
	  if ( chr == 0 ) return; // should not happen, but in case
	  start.chromosome( chr );
	  start.position( 1 );
	  stop.chromosome( chr );
	  stop.position( 300000000 );      
	  flag = true;
	  return;
	}
    }


  // Get chromosome 
  if ( ! Helper::chr_known( s.substr( 0,p ) ) ) return;
    
  int chr = Helper::chrCode( s.substr( 0,p ) );
  if ( chr == 0 ) return;
  
  
  // In form chr1:1234:rs1234?  remove rs component
  
  std::string r = s.substr(p+1);
  size_t w = r.find(":");
  std::string spos = w != string::npos ? r.substr(0,w) : r ;
  
  // Are one or two positions specified?
  
  size_t q = spos.find("..");
  
  if ( q == string::npos )
    {
      start.chromosome( chr );
      start.position( str2int( spos ) );
      stop.chromosome( chr );
      stop.position( str2int( spos ) );
      flag = true;
      return;
	}
  
  start.chromosome( chr );
  start.position( str2int(spos.substr(0,q)));
  
  stop.chromosome( chr );
  stop.position( str2int(spos.substr(q+2)));
  flag = true;
  return;
}

bool Region::overlaps(const Region& b) const
{ 
  return stop >= b.start && start <= b.stop;
}



////////////////////////////////////////////////////////////////

  
int RangeIntersector::intersect( const Region & a , void * p )
{
  
  // are appropriate functions set?

  if ( f_next == NULL ) return 0;

  // track how many intersections
  
  int cnt = 0;

  // First, compare against the existing cache, cleaning up as we go
  // For now, only clean on an all-or-nothing basis, which should typically 
  // be fine for most data-types
  
  bool any_overlap = false;
  
  std::set< Region >::iterator i = current.begin(); 
  while ( i != current.end() )
    {
      if ( a == *i ) 
	{
	  if ( f_report ) f_report( a , *i , p );
	  ++cnt;
	  any_overlap = true;
	}
      ++i;
    }
  
  // Clear cache?
  
  if ( ! any_overlap ) 
    current.clear();
  
  // Now, look to the database to get the next set of regions
  
  while ( ! finished ) 
    {

      Region b;      
      
      if ( ! f_next( b , p ) ) 
	{
	  finished = true;
	  break;
	}
      
      if ( b.after( a ) ) 
	{
	  current.insert(b);
	  break;
	}
      
      if ( b == a ) // overlap? 
	{
	  if ( f_report ) f_report( a , b , p );
	  ++cnt;
	  current.insert( b );	  
	}
            
    }
  
  return cnt;
}


