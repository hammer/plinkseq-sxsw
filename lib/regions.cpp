#include "plinkseq/regions.h"
#include "plinkseq/gstore.h"

extern GStore * GP;

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

  // Parse a region string into a Region class
  
  // chr1
  // chr1:100
  // 1:100
  // chr1:100..200
  // chr1:100-200
  // chr1:150+50  --> chr1:100..200

  // lookup from VARDB if can't interpret as a number (or chromosome)
  // rs12345      
  // rs12345+50
  // rs12345..rs67890

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
  if ( ! Helper::chr_known( s.substr( 0,p ) ) ) 
    {
      // is this a rs-ID (or two rs-IDs) instead?
      if  ( GP->vardb.attached() && ! GP->single_file_mode() ) 
	{
	  // could be in one of four forms:
	  // rs1234
	  // rs1234..rs5678
	  // rs1234-rs5678
	  // rs1234+10000
	  
	  // Are one or two positions specified?
	  
	  size_t q1 = s.find("..");
	  size_t q2 = s.find("-");
	  size_t q3 = s.find("+");
  
	  // Single ID
	  if ( q1 == string::npos && q2 == std::string::npos && q3 == std::string::npos )
	    {
	      bool okay ;
	      Region r = GP->vardb.get_position_from_id( s , "" , &okay );
	      if ( ! okay ) return; // could not find a position	      
	      start.chromosome( r.start.chromosome() );
	      stop.chromosome( r.stop.chromosome() );
	      start.position( r.start.position() );
	      stop.position( r.stop.position() );
	      flag = true;
	      return;
	    }
	  
	  if       ( q1 != std::string::npos ) // '..'
	    { 
	      bool okay ;
	      Region r = GP->vardb.get_position_from_id( s.substr(0,q1) , s.substr( q1+2 ) , &okay );
	      if ( ! okay ) return; // could not find a position	      
	      start.chromosome( r.start.chromosome() );
	      stop.chromosome( r.stop.chromosome() );
	      start.position( r.start.position() );
	      stop.position( r.stop.position() );
	      flag = true;
	      return;
	    }
	  else if  ( q2 != std::string::npos )   // '-'
	    { 
	      bool okay ;
	      Region r = GP->vardb.get_position_from_id( s.substr(0,q2) , s.substr( q2+1 ) , &okay );
	      if ( ! okay ) return; // could not find a position	      
	      start.chromosome( r.start.chromosome() );
	      stop.chromosome( r.stop.chromosome() );
	      start.position( r.start.position() );
	      stop.position( r.stop.position() );
	      flag = true;
	      return;
	    }
	  else if  ( q3 != std::string::npos )  // '+'
	    { 

	      bool okay ;
	      int pmid = 0, w = 0;
	      if ( ! Helper::str2int( s.substr( q3+1 ) , w ) ) return;
	      Region r = GP->vardb.get_position_from_id( s.substr(0,q3) , "" , &okay );
	      if ( ! okay ) return; // could not find a position	      
	      pmid = r.start.position();

	      start.chromosome( r.start.chromosome() );
	      stop.chromosome( r.stop.chromosome() );
	      start.position( pmid - w < 0 ? 0 : pmid - w );  // impose positive bounds
	      stop.position( pmid + w );
	      flag = true;
	      return;
	    }
	}
    }

  // if the above didn't work, and we still do not recognize this region, bail
  int chr = Helper::chrCode( s.substr( 0,p ) );
  if ( chr == 0 ) return;
    
  // In form chr1:1234:rs1234?  remove rs component  
  std::string r = s.substr(p+1);
  size_t w = r.find(":");
  std::string spos = w != string::npos ? r.substr(0,w) : r ;
  
  // Are one or two positions specified?
  
  size_t q1 = spos.find("..");
  size_t q2 = spos.find("-");
  size_t q3 = spos.find("+");

  // Single position, e.g. chr1:12345 
  if ( q1 == string::npos && q2 == std::string::npos && q3 == std::string::npos )
    {
      int pstart = 0;
      if ( ! Helper::str2int( spos , pstart ) ) return;
      start.chromosome( chr );
      start.position( pstart );
      stop.chromosome( chr );
      stop.position( pstart );
      flag = true;
      return;
    }

  if       ( q1 != std::string::npos ) // '..'
    { 
      int pstart = 0, pend = 0;
      if ( ! Helper::str2int( spos.substr(0,q1) , pstart ) ) return; 
      if ( ! Helper::str2int( spos.substr( q1+2 ) , pend ) ) return;
      if ( pend < pstart ) return;
      start.chromosome( chr );
      start.position( pstart );
      stop.chromosome( chr );
      stop.position( pend );
    } 
  else if  ( q2 != std::string::npos )   // '-'
    { 
      int pstart = 0, pend = 0;
      if ( ! Helper::str2int( spos.substr(0,q2) , pstart ) ) return; 
      if ( ! Helper::str2int( spos.substr( q2+1 ) , pend ) ) return;
      if ( pend < pstart ) return;
      start.chromosome( chr );
      start.position( pstart );
      stop.chromosome( chr );
      stop.position( pend );    
    } 
  else if  ( q3 != std::string::npos )  // '+'
    { 
      int pmid = 0, w = 0;
      if ( ! Helper::str2int( spos.substr(0,q3) , pmid ) ) return; 
      if ( ! Helper::str2int( spos.substr( q3+1 ) , w ) ) return;

      start.chromosome( chr );
      start.position( pmid - w < 0 ? 0 : pmid - w );  // impose positive bounds
      stop.chromosome( chr );
      stop.position( pmid + w );          
    }
  
  //  std::cout << "region = [" << *this << "] \n";

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





std::set<Region> RegionHelper::region_merge_overlap( const std::set<Region> & r )
{
  std::set<Region> m;

  if ( r.size() == 0 ) return m;

  Region current = *r.begin();

  std::set<Region>::iterator ii = r.begin();
  while ( ii != r.end() ) 
    {
      if ( ii->chromosome() != current.chromosome() ) 
	{
	  m.insert( current );
	  current = *ii;
	}
      else if ( ii->start.position() > current.stop.position() ) 
	{
	  m.insert( current ) ;
	  current = *ii;
	}
      else
	{
	  if ( ii->stop.position() > current.stop.position() )
	    current.stop.position( ii->stop.position() ) ;
	}
      ++ii;
    }
  m.insert( current );
  return m;
}

double RegionHelper::region_span( const std::set<Region> & r)
{
  // assumes non-overlapping regions
  double p = 0;
  std::set<Region>::iterator ii = r.begin();
  while ( ii != r.end() ) 
    {
      p += ii->length();
      ++ii;
    }
  return p;
}


std::set<Region> RegionHelper::region_remove( const std::set<Region> & init , 
					      const std::set<Region> & ex )
{
  
  // we assume two 'flattened' lists (i.e. both sent via
  // region_merge_overlap() above
  
  // given regions init, return a new list that removes any parts in
  // 'ex' ie. this may involve returning partial regions, or splitting
  // a single region into multiple regions, if it is partially
  // disrupted by one or more segments in 'ex'

  std::set<Region> r;
  
  std::set<Region> cache;
  std::set<Region>::iterator xi = ex.begin();
  int cache_chr = xi->chromosome();
  int cache_bp  = xi->stop.position();
  
  std::set<Region>::iterator ii = init.begin();
  while ( ii != init.end() ) 
    {
      
      // chance to clear cache?
      
      if ( ii->chromosome() != cache_chr ||
	   ii->start.position() > cache_bp ) 
	{
	  cache.clear();
	}

      // add items to the cache until we are past this segment
      
      while ( xi != ex.end() ) 
	{
	  if ( xi->before( *ii ) ) { ++xi; continue; }
	  if ( xi->after( *ii ) ) break;
	  // implies overlap
	  cache.insert( *xi );
	  if ( xi->stop.position() > cache_bp ) 
	    cache_bp = xi->stop.position() ;
	  ++xi;
	}

      // now the cache should contain any potential overlaps

      Region curr = *ii;

      bool finished = false;

      std::set<Region>::iterator ci = cache.begin();
      while ( ci != cache.end() )
	{
	  // 0) if X starts before the region, just chop of the start (or all)
	  if ( ci->start.position() <= curr.start.position() ) 
	    {
	      // if whole rest of interval covered, we are done
	      if ( ci->stop.position() >= curr.stop.position() ) { finished = true; break;  }
	      
	      // else, update definition of 'current' segment (i.e. rightmost split)
	      curr.start.position( ci->stop.position() + 1 ) ; 
	      
	    }

	  // 1) if X starts after start of region, peel off that region as a separate one, and make 
	  
	  if ( ci->start.position() > curr.start.position() && ci->start.position() <= curr.stop.position() ) 
	    {

	      r.insert( Region( cache_chr , curr.start.position() , ci->start.position() - 1 ) );

	      // if whole rest of interval covered, we are done
	      if ( ci->stop.position() >= curr.stop.position( ) ) { finished = true; break; } 
	      
	      // else, update definition of 'current' segment (i.e. rightmost split)
	      curr.start.position( ci->stop.position() + 1 ) ; 
	      
	    }

	  ++ci;
	  
	}

      // is anything left to add (could be entire segment)
      
      if ( ! finished ) r.insert( curr );

      // consider next 

      ++ii;
    }

  return r;
}



std::set<Region> RegionHelper::region_require( const std::set<Region> & init , 
					      const std::set<Region> & req )
{
  
  // we assume two 'flattened' lists (i.e. both sent via
  // region_merge_overlap() above
  
  // as above, except this returns the complement of region_remove()
  //  i.e. only sub-regions that overlap something in req.

  std::set<Region> r;
  
  std::set<Region> cache;
  std::set<Region>::iterator xi = req.begin();
  int cache_chr = xi->chromosome();
  int cache_bp  = xi->stop.position();
  
  std::set<Region>::iterator ii = init.begin();
  while ( ii != init.end() ) 
    {
      
      // chance to clear cache?
      
      if ( ii->chromosome() != cache_chr ||
	   ii->start.position() > cache_bp ) 
	{
	  cache.clear();
	}

      // add items to the cache until we are past this segment
      
      while ( xi != req.end() ) 
	{
	  if ( xi->before( *ii ) ) { ++xi; continue; }
	  if ( xi->after( *ii ) ) break;
	  // implies overlap
	  cache.insert( *xi );
	  if ( xi->stop.position() > cache_bp ) 
	    cache_bp = xi->stop.position() ;
	  ++xi;
	}

      // now the cache should contain any potential overlaps

      Region curr = *ii;
      
      bool finished = false;

      std::set<Region>::iterator ci = cache.begin();
      while ( ci != cache.end() )
	{
	  // 0) if X starts before the region and overlaps
	  if ( ci->start.position() <= curr.start.position() ) 
	    {
	      // if whole rest of interval covered, add whole region
	      if ( ci->stop.position() >= curr.stop.position() ) 
		r.insert( curr );
	      
	      // else, add partial, and change definition of current segment
	      else 
		{
		  r.insert( Region( cache_chr , curr.start.position() , ci->stop.position() ) );
		  curr.start.position( ci->stop.position() + 1 ) ; 
		}	      
	    }

	  // 1) if X starts after start of region, peel off that region as a separate one, and make 
	  
	  if ( ci->start.position() > curr.start.position() && ci->start.position() <= curr.stop.position() ) 
	    {
	      
	      // if whole rest of interval covered, we are done
	      if ( ci->stop.position() >= curr.stop.position( ) )
		{
		  r.insert( Region( cache_chr , ci->start.position() , curr.stop.position() ) );
		  finished = true;
		  break;
		}
	      
	      // else, insert partial segment and update definition of 'current' segment (i.e. rightmost split)
	      else 
		{
		  r.insert( Region( cache_chr , ci->start.position() , ci->stop.position() ) );
		  curr.start.position( ci->stop.position() + 1 ) ; 
		}
	      
	    }

	  if ( finished ) break;

	  ++ci;
	  
	}

      // consider next 

      ++ii;
    }

  return r;
}
