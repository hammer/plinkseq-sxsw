#include "ibs.h"
#include "func.h"

extern GStore g;
extern Pseq::Util::Options args;

void Pseq::IBS::regargs( Pseq::Util::Options * args )
{
  args->reg( "long-format" , Pseq::Util::Options::NONE , "output one row per line (instead of n-by-n matrix)" );
  args->reg( "two-counts"  , Pseq::Util::Options::NONE , "only consider variants seen twice in the sample" );
}


bool Pseq::IBS::calculate( Mask & m ) 
{
  
  Out & pout = Out::stream( "ibs" );
  
  const int n = g.indmap.size();
  
  class Pseq::IBS::Aux aux( n ) ;
  
  g.vardb.iterate( f_IBS_calculator , &aux , m );
  
  bool matrix = ! args.has("long-format");
  bool two_counts = (!matrix) && args.has("two-counts");
  
  if ( matrix ) 
    {
      pout << "IBS";
      for (int i=0; i<n; i++) pout << "\t" << g.indmap(i)->id();
      pout << "\n";
      
      for (int i=0; i<n; i++)
	{
	  pout << g.indmap(i)->id() << "\t";
	  for (int j=0; j<n; j++)
	    pout << ( aux.obs(i,j) > 0 ? (double)aux.ibs(i,j) / (double)aux.obs(i,j) : 0 ) << "\t";
	  pout << "\n";
	}
    }
  else
    {
      for (int i=0; i<n; i++)
	for (int j=i+1; j<n; j++)
	  {
	    pout << g.indmap(i)->id() << "\t" 
		 << g.indmap(j)->id() << "\t";
	    
	    if ( two_counts ) 
	      pout << aux.ibs(i,j) << "\t" << aux.obs(i,j) << "\n";
	    else 
	      pout << ( aux.obs(i,j) > 0 ? (double)aux.ibs(i,j) / (double)aux.obs(i,j) : 0 ) << "\n";
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
