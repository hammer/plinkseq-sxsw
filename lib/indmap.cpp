
#include "indmap.h"
#include "phmap.h"
#include "vardb.h"

using namespace std;
using namespace Helper;

class Individual;

void IndividualMap::reset()
{
  ialign.clear();
  ids.clear();
  id2pos.clear();
  wsample.clear();
  uniq.clear();
  mult.clear();
  indiv.clear();
  idvec.clear();
  sample_indiv.clear();
  sample_idvec.clear();
  n_uniq = 0;
  is_flat = false;
  unflattering = false;
}


void IndividualMap::insert( const int file_id, 
			    const int supposed_j , 
			    const int slot_j, 
			    const int slot_k, 
			    const std::string & id )
{
  int2 a( file_id , slot_j );
  int_string_pair  b( slot_k , id );
  ialign.insert( std::make_pair( a,b ) );
  wsample[ file_id][supposed_j] = slot_j;  
  con2svar[ int2(file_id,slot_k) ] = slot_j;
}

int IndividualMap::sample_slot( const int i , const int file_id ) const
{
  // if f==0 then return 'i' as this is the consensus sample variant
  if ( file_id == 0 ) return i;
  std::map<int2,int>::const_iterator ii = con2svar.find( int2(file_id , i) );
  return ii == con2svar.end() ? -1 : ii->second;
}

std::vector<int> * IndividualMap::svar2consensus( int f )
{
  if ( wsint.find(f) == wsint.end() ) 
    Helper::halt("internal error: trying to reference non-existent svar, svar2consensus()");
  return &wsint[f];
}


void IndividualMap::construct_wsint_vector()
{

  // After insertion, given we have wsample made as maps, also 
  // make some basic std::vector<int> for each file, to be used
  // when mapping svar->consensus slots in multi-sample, flat alignments

  // probably a really dumb way to do the below, but who cares we only 
  // do it once

  std::map<int,int> ns;  
  std::map<int2,int_string_pair>::const_iterator i = ialign.begin();  
  while ( i != ialign.end() )
    {
      ns[ i->first.p1 ]++;
      ++i;
    }
  
  std::map<int,int>::iterator j = ns.begin();
  while ( j != ns.end() )
    {
      wsint[ j->first ].resize( j->second );
      ++j;
    }
  
  i = ialign.begin();  
  while ( i != ialign.end() )
    {
      wsint[ i->first.p1 ][ i->first.p2 ] = i->second.i ;
      ++i;
    }
}


std::map<int,std::string> IndividualMap::map_slot_to_id() const
{
  std::map<int,std::string> m;
  std::map<int2,int_string_pair>::const_iterator k = ialign.begin();
  while ( k != ialign.end() )
    {
      m[ k->second.i ] = k->second.s;
      ++k;
    }
  return m;
}



//
// Create an individual-map (given a VARDB and a Mask)
//

int IndividualMap::populate( VarDBase & vardb, PhenotypeMap & phmap , const Mask & m )
{

  // First clear the map
  reset();
  
  // check the necessary databases are attached
  if ( ! vardb.attached() ) return 0;
  
  std::map<std::string,int> id_list;
  std::map<int, std::map<int,int> > tracker;
  std::map<int, Individual*> ptracker;

  // Get a list of all individuals from the VARDB, apply Mask, then
  // load into the database

  std::map<int,std::string> files = vardb.fetch_files();
  std::map<int,std::string>::iterator f = files.begin();

  // track files actually appearing in this map, to 
  // assign 'multi-sample' status or not

  std::set<int> obs_files;
  
  //
  // Start a transaction
  //
  
  phmap.begin();
  
  while ( f != files.end() )
    {      
      
      // File-level screen
      if ( ! m.use_file( f->first ) ) 
	{
	  ++f;
	  continue;
	}

      int actual_j = 0;

      std::vector<std::string> ind = vardb.fetch_individuals( f->first );
      
      for (int j = 0 ; j < ind.size(); j++ ) 
	{
	  
	  std::string id = ind[j];

	  // Are we interested in this individual?  
	  
  	  if ( ! m.use_indiv( id ) ) 
 	    continue;

	  
	  // Pull individual into the phenotype-map

	  Individual * person = phmap.new_individual( id );	  
	  
	  
	  // Phenotype-based mask-screen (missingness, value)
	  
	  if ( ! m.pheno_screen( person ) )
	    continue;

	  
	  //
	  // Track the number of unique files observed
	  //
	  
	  obs_files.insert( f->first );


	  // Store mapping as appropriate
	  
	  if ( id_list.find( id ) == id_list.end() ) 
	    {
	      int n = id_list.size();
	      id_list.insert( make_pair( id , n ) );
	      tracker[n].insert( make_pair(f->first, actual_j ) );
	      ptracker[n] = person;
	    }
	  else
	    {
	      tracker[ id_list[id] ].insert( make_pair(f->first, actual_j) );
	      ptracker[ id_list[id] ] = person;
	    }

	  // Store mapping of {file/slot} --> { slot/id-string }
	  
	  insert( f->first , j, actual_j++ , id_list[id] , id );
	  
	  // Simply list of people in map
	  
	  ids.insert(id);

	}

      // REMOVED -- this is a bad idea; and means that we don' see missing 
      // genotypes etc for people, e.g. from g-view

      // If no individuals were entered from this file, remove it from the Mask
      // NOTE: need to see how this interacts with other uses of VARDB, eg. from ACDB
      // in which there are no individual-level data...  or, downstream, if we have 
      // some optimisations which only consider variant level data

//       if ( obs_files.find( f->first ) == obs_files.end() )
// 	const_cast<Mask&>(m).exclude_file( vardb.file_tag( f->first ) );

      
      // Move on to next file.

      ++f;
    }



  //
  // End INDDB transaction
  //
  
  phmap.commit();
  

  //
  // Number of files in VARDB used
  //
  
  is_multi_sample = obs_files.size() > 1;

  //
  // Number of unique individuals to be extracted
  //

  const int n = id_list.size();

  //
  // Handle n==0 cases specially
  //

  //  if ( n == 0 ) { is_multi_sample = true; }

  //
  // (-1,-1) in uniq implies this person seen more than once, and so
  // lives in mult[]
  //
  
  uniq.resize( n , int2(-1,-1) );   
  

  //
  // Track people who appear more than once. Ultimately, these will be 
  // resolved into the consensus/uniq set, but let the user do this 
  // how they wish, downstream.
  //

  mult.resize( n );


  //
  // Record whether or not any merging/resolving needs to be performed
  // (i.e. is the same genotype ever seen more than once?) If not, we 
  // can make some savings when constructing a Variant
  //
  
  is_flat = true;
  

  //
  // Handle a special case: if no genotype-level data are supplied, then 
  // do not treat as a flat alignment.  Need to revisit this, but should 
  // make life easier with ACDB
  //

  //  if ( n == 0 ) is_flat = false;


  //
  // Look for a non-flat alignment
  //
  
  std::map<int,map<int,int> >::iterator i = tracker.begin();
  
  while ( i != tracker.end() )
    {

      // track position in files

      std::map<int,int> & m = i->second;

      if ( m.size() == 1 ) 
	{
	  std::map<int,int>::iterator j = m.begin();
	  uniq[ i->first ] = int2( j->first, j->second );	  
	}
      else
	{
	  
	  is_flat = false;
	  
	  std::map<int,int>::iterator j = m.begin();
	  while ( j != m.end() )
	    {
	      mult[ i->first ].insert( int2( j->first, j->second ) );
	      ++j;
	    }
	}

      ++i;
    }
  
  // Store the number of unique individuals
  
  n_uniq = uniq.size();
  


  //
  // Allow appropriate size for pointers to actual data
  //

  indiv.resize( n_uniq , NULL );
  
  // And similarly for ID strings

  idvec.resize( n_uniq , "" );

  std::map<int,Individual*>::iterator pi = ptracker.begin();
  while ( pi != ptracker.end() )
    {

      indiv[ pi->first ] = pi->second;
      idvec[ pi->first ] = pi->second->id();
      
      // Add to ID -> n map
      
      id2pos[ pi->second->id() ] = pi->first;
      
      
      // And add sample-level trackers too
      
      if ( uniq[ pi->first ] == int2(-1,-1) )
	{
	  int2 k = uniq[ pi->first ];
	  sample_indiv[ k.p1 ][ k.p2 ] = pi->second;
	  sample_idvec[ k.p1 ][ k.p2 ] = pi->second->id();
	}
      else
	{
	  std::set<int2>::iterator j = mult[ pi->first].begin();
	  while ( j != mult[ pi->first].end() )
	    {
	      sample_indiv[ j->p1 ][ j->p2 ] = pi->second;
	      sample_idvec[ j->p1 ][ j->p2 ] = pi->second->id();
	      ++j;
	    }
	}
      ++pi;
    }


  
  // Align with PhenotypeMap
  
  phmap.align( ids );
  
  // Populate some more convenient structures used in 
  // Svar genotype look-ups for flat alignments

  construct_wsint_vector();
  
  // return number of people in map

  return id_list.size();
}





//
// Create an simple flat map from a ID list (no VARDB, Mask or INDDB assumed)
//

int IndividualMap::populate( const std::vector<std::string> & ind )
{

  // NOTE -- this is a quick and dirty hack to allow a simple population
  //         of the map via a specification of IDs -- i.e. for case 
  //         where people do not come from the DB, e.g. simulation, etc.
  
  // As such, not all structures (e.g. those dealing with how individuals
  // are distributed across files, or with phenotypes) will necessarily be
  // fully / validly populated.  Need to check.

  reset();
 
  n_uniq = ind.size();

  indiv.resize( n_uniq , NULL );
  idvec.resize( n_uniq , "" );

  for (int i = 0 ; i < ind.size(); i++ ) 
    {            
      Individual * person = new Individual( ind[i] );
      indiv[ i ] = person;
      idvec[ i ] = person->id();      
      id2pos[ person->id() ] = i;
    }

  is_multi_sample = false;
  is_flat = true;
  return n_uniq;
}
