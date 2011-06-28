#include "ibs.h"
#include "func.h"

extern GStore g;
extern Pseq::Util::Options options;

bool Pseq::IBS::calculate( Mask & m ) 
{
  const int n = g.indmap.size();
  
  class Pseq::IBS::Aux aux( n ) ;
  
  g.vardb.iterate( f_IBS_calculator , &aux , m );
  
  bool matrix = ! options.key("long-format");
  bool two_counts = (!matrix) && options.key("two-counts");
  
  if ( matrix ) 
    {
      plog << "IBS";
      for (int i=0; i<n; i++) plog << "\t" << g.indmap(i)->id();
      plog << "\n";
      
      for (int i=0; i<n; i++)
	{
	  plog << g.indmap(i)->id() << "\t";
	  for (int j=0; j<n; j++)
	    plog << ( aux.obs(i,j) > 0 ? (double)aux.ibs(i,j) / (double)aux.obs(i,j) : 0 ) << "\t";
	  plog << "\n";
	}
    }
  else
    {
      for (int i=0; i<n; i++)
	for (int j=i+1; j<n; j++)
	  {
	    plog << g.indmap(i)->id() << "\t" 
		 << g.indmap(j)->id() << "\t";

	    if ( two_counts ) 
	      plog << aux.ibs(i,j) << "\t" << aux.obs(i,j) << "\n";
	    else 
	      plog << ( aux.obs(i,j) > 0 ? (double)aux.ibs(i,j) / (double)aux.obs(i,j) : 0 ) << "\n";
	  }
    }

}



/// Worker functions
    
void f_IBS_calculator( Variant & v , void * p )
{

  Pseq::IBS::Aux * aux = (Pseq::IBS::Aux*)p;
  
  if ( ! v.biallelic() ) return;

  bool altmin = v.n_minor_allele();
  
  const int n = v.size();
  
  for (int i=1; i<n; i++)
    {
      if ( v(i).null() ) continue;
      
      const int s = v(i).minor_allele_count( altmin );
      
      for (int j=0; j<i; j++)
	{
	
	  // count only similar non-ref alleles at one of two indiv	  
	  if ( ! v(j).null() )
	    {	      
	      // genotype counting (1=both non-reference)
 	      aux->obs(i,j,1);	      
	      const int t = v(j).minor_allele_count( altmin );
	      aux->ibs(i , j , s > 0 && t > 0 );
	    }
	}
    }
}
