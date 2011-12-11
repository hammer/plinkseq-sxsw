
#include "segments.h"

using namespace std;
using namespace Helper;

void OverlapResults::load_regions(set<Region> r)
{
  
  std::set<Region>::iterator i = r.begin();
  
  while ( i != r.end() )
    {

      OverlapResult olap;

      olap.totalLength = i->stop.position() - i->start.position() + 1;
      olap.nExons = i->subregion.size();
      
      for (int s=0; s<i->subregion.size(); s++)
	olap.exonLength += i->subregion[s].stop.position() - i->subregion[s].start.position() + 1;
      
      result.insert(make_pair( *i , olap ) );
      
      ++i;
    }
}


void func_locdb_process_overlap(Region & r1, Region & r2, int v_int, int v_union, void * d)
{
  
  OverlapResults * data = (OverlapResults*)d;
  
  // Ignore if tested gene overlaps tested gene  
    if ( r1.group == r2.group ) return;
    
  // We now have 1 tested gene (with exons) and 1 other region
  
  Region & tested   = r1.group == data->target_id ? r1 : r2; 
  
  Region & target   = r1.group == data->target_id ? r2 : r1; 
  
  // Have we seen this tested gene before? 
  
  std::map<Region,OverlapResult>::iterator i = data->result.find( tested );
  
  if ( i != data->result.end() )
    {
      
      OverlapResult & r = i->second;
      
      r.overlapping_target_region_ids.insert( target.id );
      
      // Which subregions does this region overlap with? 
      
	for (int s=0; s< tested.subregion.size(); s++)
	  {	  
	    for (int s2 = 0; s2 < target.subregion.size(); s2++)
	      {
		if ( target.subregion[s2].overlaps( tested.subregion[s] ) )
		{
		    
		  // Count another target that at least partially spans this gene
		  
		  r.nTargets++;
		  
		  
		  std::map<int,set<int2> >::iterator i = r.cover.find(s);		  
		  if ( i == r.cover.end() )
		    {
		      set<int2> t;
		      t.insert( int2( target.subregion[s2].start.position(), 
				      target.subregion[s2].stop.position() ) );
		      r.cover.insert(make_pair(s,t));
		    }
		  else
		    {
		      i->second.insert( int2( target.subregion[s2].start.position(), 
					      target.subregion[s2].stop.position() ) ) ;
		    }
		}
	      }
	  }
    }
}


void GStore::locdb_load_names( string file , string name )
{
    uint64_t candidate_id = locdb.load_regions( file , name , -1, -1, -1 , 0 );
}


void GStore::locdb_summary()
{
    set<GroupInfo> ginfo = locdb.group_information();    
    set<GroupInfo>::iterator i = ginfo.begin();
    while ( i++ != ginfo.end() )
    {
      plog << *i << "\t"
	   << locdb.count( i->idx ) << " records; "
	   << locdb.span( i->idx )/1000.0 << "kb span\n";	    
    }
    
    // TODO: output # of subregions
    
    plog << "\n";
}


void GStore::locdb_remove_group( ID_t id )
{
    locdb.flush( id );
}

void GStore::locdb_subregions(bool t)
{
    locdb.get_subregions(t);
}

void GStore::locdb_meta(bool t)
{
    locdb.get_meta(t);
}



void GStore::locdb_rename( string group , string alias, string new_label )
{
    ID_t group1 = locdb.lookup_group_id( group );
    ID_t alias1 = locdb.lookup_group_id( alias);
    if ( group1 == 0 || alias1 == 0 ) return;
    ID_t renamed_id = locdb.rename( group1, alias1, new_label );
}


void GStore::locdb_extract_intersection(string group1, string group2 , string newLabel)
{
    // Check groups exist
    ID_t g1 = locdb.lookup_group_id( group1 );
    ID_t g2 = locdb.lookup_group_id( group2 );
    if ( g1 == 0 || g2 == 0 ) return;
    uint64_t extract1_id = locdb.extract( g1 , g2 , "newLabel" );    
} 
  



void GStore::locdb_display_regions( string name )
{
    ID_t id = locdb.lookup_group_id( name );

    if ( id == 0 ) 
	return;

    set<Region> regions = locdb.get_regions( id );
    
    set<Region>::iterator i = regions.begin();
    
    while ( i != regions.end() )
    {
	plog << *i << "\n";
	++i;
    }
}


void GStore::locdb_overlap_analysis( const std::string & target, 
				     const std::string & preload , 
				     const std::string & alias , 
				     const std::string & listmode )
{
  
  ID_t target_id = locdb.lookup_group_id( target );
  ID_t preload_id = locdb.lookup_group_id( preload );
  ID_t alias_id = locdb.alias_id( alias );

  if ( target_id == 0 || preload_id == 0 ) return;

  plog << "TARGET" << "\t"
       << ( alias != "" ? "ALIAS" : "" ) << "\t"
       << "POS1" << "\t"
       << "POS2" << "\t"
       << "N_SUB" << "\t"
       << "LEN" << "\t"
       << "LEN_SUB" << "\t"
       << "N" << "\t"
       << "OLAP1" << "\t"
       << "OLAP2" << "\n";

  bool osub = locdb.get_subregions();
  bool ometa = locdb.get_meta();
  
  // Turn off subregion and meta-data reporting to speed things up

  locdb.get_subregions( false );   
  locdb.get_meta( false );
  
  locdb.add_overlap_table( target_id , preload_id );

  locdb.get_subregions( osub );
  locdb.get_meta( ometa );


  // Pull the set of target regions into memory
  
  std::set<Region> regions = locdb.get_regions( preload_id );

  // Helper class to store overlap results
  
  OverlapResults res( target_id );

  res.load_regions( regions );

  locdb.get_regions_and_overlap( &func_locdb_process_overlap , &res );

  std::map<Region,OverlapResult>::iterator j = res.result.begin();
  
  while( j != res.result.end() )
    {
      
      // Name of target gene
      
      const Region & target = j->first; 
      
      // Exon-based overlap results
      
      OverlapResult & olap = j->second;

      std::string aliases = locdb.alias( target.name , preload_id , alias_id ); 

      std::set<uint64_t>::iterator ii = olap.overlapping_target_region_ids.begin();
      
      if ( ii == olap.overlapping_target_region_ids.end() )
	{
	  plog << ".\n";
	}
      else
	{
	  while ( ii != olap.overlapping_target_region_ids.end() )
	    {

	      Region oregion = locdb.get_region( *ii );

	      plog << target.id << "\t"
		   << target.name << "\t"
		   << aliases << "\t"
		   << target.coordinate() << "\t"
		   << oregion.id << "\t"
		   << oregion.name << "\t"
		   << oregion.coordinate() << "\n";

	      ++ii;
	    }

	}

      // skip actual olap calcs for now.


      ++j;
      continue;

      // SKIP THE BELOW FOR NOW
      

      
      plog << target.name << "\t"
	   << aliases << "\t"
	   << target.coordinate() << "\t"
	   << target.start.chromosome() << "\t"
 	   << target.start.position() << "\t"
 	   << target.stop.position() << "\t";
      
      plog << j->second.nExons << "\t"
	   << (double)(j->second.totalLength)/1000.0 << "\t"
	   << j->second.exonLength << "\t"
	   << j->second.nTargets << "\t";
      

      // Calculate overlap metrics
      
      int tc = 0;
      
      for ( int s = 0 ; s < target.subregion.size(); s++ )
	{
	  
	  std::map<int,std::set<int2> >::iterator i = olap.cover.find( s ); 
	  
	  int exonStart = target.subregion[s].start.position();
	  int exonStop = target.subregion[s].stop.position();
	  
	  // No overlap?
	  
	  if ( i == olap.cover.end() ) 
	    {
	      tc += exonStop-exonStart+1;
	      continue;
	    }

	  // Overlap
	  
	  std::set<int2> & segments = i->second;
	  
	  std::set<int2>::iterator si = segments.begin();
	  
	  int c = 0;
	  int p = exonStart;
	  
	  while ( si != segments.end() )
	    {
	      if ( si->p1 > p )
		c += si->p1 - p;
	      if ( si->p2 > p )
		p = si->p2;
	      ++si;
	      }
	  tc += c;
	}
      
      plog << olap.exonLength - tc << "\t"
	   << 1 - ( (double)tc / olap.exonLength ) << "\n";
      
        // Next target gene
        ++j;
	
      }
  
  
  // Clear up the overlap table
  
  locdb.clear_overlaps();
    
}

