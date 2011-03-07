
#include "clusters.h"
#include "func.h"
#include "pseq.h"
#include <cmath>

extern GStore g;

bool Pseq::VarDB::cluster_scan( Mask & m )
{
  ClusterAux aux;
  IterationReport r = g.vardb.iterate( g_clusterscan , &aux , m );
  return true;
}


void g_clusterscan( VariantGroup & vars , void * p )
{

  // checks
  
  //  1) report any homozygote singletons (i.e. clusters within individuals)
  //     perhaps give an option to turn this off, i.e. for small samples
  
  //  2) any variants within 2bp of each other, where an individual
  //  has two minor alleles

  //  3) a region with high number of variants; also, individuals with
  //  high numbers of variants in a region.  Regions defined in terms
  //  of kb windows.

  
  const int n = vars.n_individuals();
  const int nv = vars.size();

  // Report any individuals where the only two alleles are in one
  // individual
  
  std::vector<bool> altmin( vars.size() , false );
  
  for (int v = 0 ; v < vars.size(); v++ )
    {
      
      int c     = 0; // minor allele
      int c_tot = 0; // total counts
      altmin[v] = vars(v).n_minor_allele( c , c_tot );      
           
      // Is this a unique doubleton? 
      
      if ( c == 2 ) 
	{
	  for (int i=0; i<n; i++ )
	    if ( vars.geno(v,i).minor_allele_count( altmin[v] )  == 2 )
	      { 
		plog << "UNIQ_HOMOZYG" << "\t"
		     << vars.name() << "\t"
		     << vars(v) << "\t"
		     << vars.ind(i)->id() << "\n";		
	      }
	}
    }
  

  // 
  // Number of variants per individual
  //
  
  std::vector<int> vcnt(n,0);
  for ( int i=0; i<n; i++ )
    for ( int v=0; v<nv; v++ )      
      if ( vars(v,i).minor_allele( altmin[v] ) ) ++vcnt[i];
  
  double mean = 0 , var = 0;
  for ( int i=0; i<n; i++ )
    {
      mean += vcnt[i];
      var += vcnt[i] * vcnt[i];
    }  
  mean /= (double)n;
  var /= (double)n;
  double sd = sqrt( var - mean*mean );
  double thresh = mean + 5 * sd;
  
  // any individuals above 5SDs in count?
  for ( int i=0; i<n; i++ )
    if ( vcnt[ i] > thresh )
      plog << "EXCESS_VARS" << "\t"
	   << vars.name() << "\t"
	   << vars.size() << "\t"
	   << vars.coordinate() << "\t"
	   << vars.ind(i)->id() << "\t"
	   << vcnt[i] << "\t"
	   << mean << "\t"
	   << sd << "\t"
	   << thresh << "\n";
  
}
