#include "compare.h"
#include "pseq.h"

extern GStore g;


//
// Group comparison function: wrapper 
//

bool Pseq::Assoc::group_comparison( Mask & m )
{
  
  // Assume an n-level group (PHE_FACTOR) has been specified: tabulate
  // variants by these groups
  
  if ( g.phmap.type() != PHE_FACTOR )
    Helper::halt("requires a factor-level phenotype");
  
  std::map<std::string,int> groups = g.phmap.summarise_phenotype();
  
  std::map<std::string,int>::iterator i = groups.begin();
  while ( i != groups.end() )
    {
      plog << i->first << "\t" 
		<< i->second << "\n";
      ++i;
    }
  
  
  //
  // Aux structure for params, results
  // 
  
  Pseq::Assoc::AuxGroupComp aux(&g);
    
    
  //
  // Get counts
  //
  
  g.vardb.iterate( f_group_comp , &aux , m );
    
    
  //
  // Display results
  //
    
  plog << "\n\n";
  
  for (unsigned int s = 1; s< groups.size(); s++)
    {
      std::map< std::set<int> , int >::iterator j = aux.counts.begin();
      
      while ( j !=  aux.counts.end() )
	{
	  if ( j->first.size() != s ) 
	    {
	      ++j;
	      continue;
	    }
	  
	  plog << j->second << "\t";
	  
	  std::set<int>::iterator k = j->first.begin();
	  while ( k != j->first.end() )
	    {
	      if ( k != j->first.begin() ) plog << ",";
	      plog << Individual::label( *k );
	      ++k;
	    }
	  plog << "\n";
	  ++j;
	}
    }
  
  return true;
    
}



//
// Group comparison function: worker 
//

void f_group_comp( Variant & v , void * p)
{
  
  Pseq::Assoc::AuxGroupComp * data = (Pseq::Assoc::AuxGroupComp*)p;
  
  std::set<int> obs;
  
  const int n = v.size();
  bool altmin = v.n_minor_allele( );
  
  for (int i = 0; i < n; i++)
    if ( v(i).minor_allele( altmin ) )
      obs.insert( v.ind(i)->group() );	  
  
  
  // (shouldn't happen but...) do not save empty group

  if ( obs.size() == 0 ) return;
      

  // Record that we've seen this
  
  data->counts[ obs ]++;
  
}

