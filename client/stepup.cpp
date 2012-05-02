
#include "genic.h"

extern GStore g;

Pseq::Assoc::Aux_hoffman_witte::Aux_hoffman_witte( bool run_this , const VariantGroup & vars, Aux_prelim * p )
{ 

  // Only do all this if we'll actually be using this test... (messy, I know)

  if ( ! run_this ) return;


  // Populate helper struct for the STEPUP test

  clear();


  // Contains an indicator of whether the 'alternate' is the minor
  // allele, altmin[];

  aux_prelim = p;

  const int K = vars.size(); 
  const int I = vars.n_individuals();
  
  // Phenotype type

  if ( g.phmap.type() == PHE_DICHOT ) dichot = true;
  else if ( g.phmap.type() == PHE_QT ) dichot = false;
  else Helper::halt( "no dichotomous or quantitative phenotype specified" );	
  
    
  // Find cut-points
  
  const int STEP_UP_CUT = 15;
  
  if( K > STEP_UP_CUT ) 
    {
      int ncuts = (int)(K / STEP_UP_CUT) + 1;
      int cutsize  = (int)(K / ncuts) + 1;
      int cur_size = 0;
      while( cur_size < K )
	{
	  piece_begin.push_back( cur_size );
	  cur_size += cutsize;
	  if(cur_size > K)
	    cur_size = K;
	  piece_end.push_back( cur_size-1 );
	}
    }
  else
    {
      // it's only 1 piece
      piece_begin.push_back( 0 );
      piece_end.push_back( K-1 );
    }
  
  // calculate the trait and genotype means (invariant to permutation)
  // either handle dichot or QTs
  
  // assumes people with missing phenotypes have been removed (which, by
  // default, they always are)
  
  tbar = 0;
  if ( dichot ) 
    for(int i=0; i<I; i++)
      tbar += vars.ind(i)->affected() == CASE;
  else 
    for(int i=0; i<I; i++)
      tbar += vars.ind(i)->qt(); 

  tbar /= (double)I;  


  // Genotypes

  xbar.resize(K,0);
  
  for( int k=0; k<K; k++ )
    {
      
      const Variant & var = vars(k);
      
      unsigned int Iobs = 0;
      
      for(int i=0; i<I; i++)
	{
	  if ( ! var(i).null() )
	    {
	      xbar[k] += var(i).minor_allele_count( aux_prelim->altmin[k] ); 
	      ++Iobs;
	    }
	}
      xbar[k] /= (double)Iobs;
      
    }    

}
  

double Pseq::Assoc::Aux_hoffman_witte::calc_stat( const VariantGroup & vars,      //   genotypes/phenotypes from plink
						  const std::vector<int> & m,     //   'mask'
						  const std::vector<double> & s,  //   sign
						  const std::vector<double> & w ) //   weight
  
{
  
  const int K = vars.size();     
  const int I = vars.n_individuals();

  
  // Test statistic
  double num=0.0, den=0.0;    
  for(int i=0; i<I; i++)
    {
      double temp = 0;
      for(int k=0; k<K; k++)
	if( m[k] == 1 )
	  {
	    if ( ! vars.geno(k,i).null() )
	      {
		const double x = vars.geno(k, i).minor_allele_count( aux_prelim->altmin[k] );
		const double t = dichot ? 
		  ( vars.ind( g.perm.pos(i) )->affected() == CASE ? 1.0 : 0.0 ) :
		  vars.ind( g.perm.pos(i) )->qt() ;		
		temp += s[k] * w[k] * ( x - xbar[k] ) * ( t - tbar ); 
	      }
	  }
      
      num += temp;
      den += temp * temp;
    }

  return (num*num) / den;
} 



double Pseq::Assoc::stat_hoffman_witte( const VariantGroup & vars,  
					Pseq::Assoc::Aux_hoffman_witte * aux,    
					std::map<std::string,std::string> * output , 
					const bool original ) 
  
{
  
  // Primary driver of the association; call R+1 times per gene; 
  // For the first call, original == T
  
  const int K = vars.size();     
  const int I = vars.n_individuals(); 
  
  // Compute the signs
  
  std::vector<double> s;
  for(int k=0; k<K; k++)
    {
      double covk = 0;

      const Variant & var = vars(k);
      
      for(int i=0; i<I; i++)
	{
	  if ( ! var(i).null() )
	    {
	      const double x = var(i).minor_allele_count( aux->aux_prelim->altmin[k] ) ; 
	      const double y = aux->dichot ? ( var.ind( g.perm.pos(i) )->affected() == CASE ? 1.0 : 0.0 ) : var.ind( g.perm.pos(i) )->qt() ;	
	      covk += ( x - aux->xbar[k]) * ( y  - aux->tbar); 
	    }
	}

      s.push_back( covk ? ( covk < 0 ? -1 : 1 ) : 0 );
    }
  

  // allow weights 
  std::vector<double> w(K , 1.0 );
  
  // now construct the masks for each piece
  std::vector<std::vector<int> > masks;
  int P = aux->piece_begin.size();
  masks.resize(P);
  std::vector<double> z(P);

  for(int p=0; p<P; p++)
    {

      masks[p].resize(K,0);
      z[p] = 0;

      bool keep_going = true;
      
      while(keep_going)
	{
	  int best_loc = -1;
      
	  // try all of the locations
	  for( int k=aux->piece_begin[p]; k<=aux->piece_end[p]; k++)
	    {
	      if( masks[p][k] == 0 )
		{
		  // toggle it on
		  masks[p][k] = 1;
		  
		  // compute the test statistic, note if better
		  double newz = aux->calc_stat( vars, masks[p], s, w );
		  
		  if( newz > z[p])
		    {
		      z[p] = newz; // store it if better
		      best_loc = k;
		    }
		  
		  // toggle it back off
		  masks[p][k] = 0;
		} 
	    }
	  
	  // store the new best if there was one...
	  if( best_loc == -1 ) keep_going = false;
	  else masks[p][best_loc] = 1;
	}
    }
  
  // finally step-up on the pieces!
  
  bool keep_going = true;
  double bestz = 0;
  std::vector<int> finalmask(K,0);

  std::vector<int> gene_used(P,0);

  while( keep_going )
    {
      int best_gene = -1;
      double best_gene_z = 0;

      for(int p=0; p<P; p++)
	{

	  if( gene_used[p] == 0 )
	    {
	      // toggle the gene on
	      for(int k=0; k<K; k++)
		finalmask[k] = finalmask[k] || masks[p][k];
	      
	      double newz = aux->calc_stat(vars, finalmask, s, w);

	      // compute the test statistic
	      if( newz > best_gene_z )
		{
		  best_gene_z = newz;
		  best_gene = p;
		}
	      	 
	      // toggle it off
	      for(int k=0; k<K; k++)
		finalmask[k] = finalmask[k] && !masks[p][k];
	    }
	}
      
      if( best_gene_z > bestz )
	{
	  bestz = best_gene_z;
	  gene_used[best_gene] = 1;
	  for(int k=0; k<K; k++)
	    finalmask[k] = finalmask[k] || masks[best_gene][k];
	}
      else
	keep_going = false;      
    }
  

  if ( original ) 
    {
      (*output)["STEPUP"] = Helper::dbl2str( bestz );
    }

  // should be +ve?  perhaps constrain to 0 if neg.
  return bestz < 0 ? -bestz : bestz ; 
  
}
  


