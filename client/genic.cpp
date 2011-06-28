#include "genic.h"

#define MATHLIB_STANDALONE
#include "Rmath.h"

extern GStore g;

void   Pseq::Assoc::prelim( const VariantGroup & vars , Aux_prelim * aux )  
{
  
  // Track observed min/max minor allele counts
  
  aux->minm = -1;
  aux->maxm = -1;
  
  // Frequency-weights
  aux->fweights.resize( vars.size() , 0 );  
  aux->acounts.resize( vars.size() , 0 );
  aux->altmin.resize( vars.size() , true );
  
  for ( int v = 0 ; v < vars.size(); v++ )
    {
      
      int c     = 0; // minor allele
      int c_tot = 0; // total counts	  
      aux->altmin[v] = vars(v).n_minor_allele( &c , &c_tot );      
      if ( ! aux->altmin[v] ) aux->refmin.insert(v);
      if ( aux->minm == -1 || c < aux->minm ) aux->minm = c;
      if ( aux->maxm == -1 || c > aux->maxm ) aux->maxm = c;
      
      aux->mc[ c ]++;
      
      // Frequency weights
      double f = (double)( 1 + c ) / (double)( 2 + c_tot );
      aux->fweights[ v ] = f > 0 && f < 1 ? 1.0 / sqrt( f * (1-f) ) : 0 ;
      aux->acounts[ v ] = c;      

      // Track # of case/control-alleles (for output)
      const int n = vars.n_individuals();
      int c_a = 0 , c_u = 0;
      for (int i = 0; i < n; i++)
	{	  	  
	  if ( ! vars.geno(v,i).null() )
	    {
	      int ac = vars.geno(v,i).minor_allele_count( aux->altmin[v] );
	      if ( ac ) 
		{
		  affType aff = vars(v).ind(i)->affected();
		  if ( aff == CASE ) c_a += ac;	      
		  else if ( aff == CONTROL ) c_u += ac;
		}
	    }
	}      
      std::string t = Helper::int2str(c_a) + "/" + Helper::int2str(c_u);      
      aux->mc_a[t]++;	
    }

  
  // Identify carriers of 1+ rare alleles (typically, but not necessarily non-reference)

  aux->carriers.clear();
  
  const int n = vars.n_individuals();  
  aux->n_a = 0;
  aux->n_t = 0;
  
  for (int i = 0; i < n; i++)
    {
      int j = g.perm.pos( i );
      
      affType aff = vars.ind( j )->affected();
      if ( aff == CASE ) aux->n_a++;
      if ( aff == CASE || aff == CONTROL ) aux->n_t++;
      
      for (int v=0; v<vars.size(); v++ )
	if ( ! vars.geno(v,i).null() && vars.geno(v,i).minor_allele( aux->refmin.find(v) == aux->refmin.end() ) )
	  aux->carriers[i].insert( v );		  	  
    }
}



double Pseq::Assoc::stat_calpha( const VariantGroup & vars , 
				 Aux_prelim * pre , 
				 Aux_calpha * aux , 
				 std::map<std::string,std::string> * output , 
				 bool original )
{
  
  if ( original )
    {
      aux->p_a = pre->n_a / (double)(pre->n_t);
      aux->p_ap_u = aux->p_a * (1-aux->p_a);
    }
  
  std::vector<int> y( vars.size() , 0);
  std::vector<int> n( vars.size() , 0);
  std::map<int, std::set<int> >::iterator i1 = pre->carriers.begin();
  
  while ( i1 != pre->carriers.end() )
    {
      int j = g.perm.pos(  i1->first  );
      affType aff = vars.ind( j )->affected();
	  
      std::set<int>::iterator k = i1->second.begin();
      while ( k != i1->second.end() )
	{
	  if ( aff == CASE ) y[ *k ]++;
	  n[ *k ]++;	      
	  ++k;
	}
      ++i1;
    }
  
  double score = 0;
  
  for (int i = 0 ; i < vars.size(); i++ )
    {          
      double t = y[i] - n[i] * aux->p_a;
      t *= t;
      score += t - n[i] * aux->p_ap_u;
    }
  

  if ( original ) 
    {
      // calculate denominator only once

      for (int m = pre->minm ; m <= pre->maxm ; m++ )
	{
	  double t = 0;
	  for ( int u = 0 ; u <= m ; u++ )
	    {
	      double s = ( u - m * aux->p_a );
	      s *= s;
	      s -= m * aux->p_ap_u;
	      s *= s; 	      
	      s *= dbinom( u , m , aux->p_a , false );
	      t += s;
	    }	      
	  aux->variance += pre->mc[m] * t;
	}
    }


  // Frequency breakdown output

  if ( original ) 
    {
      (*output)["CALPHA"] = "Z=" + Helper::dbl2str( score / sqrt( aux->variance ) )  ;
      std::map<std::string,int>::iterator i = pre->mc_a.begin();
      while ( i != pre->mc_a.end() )
	{
	  (*output)["CALPHA"] += ";" + i->first  
	    + "=" + Helper::int2str( i->second ) ;
	  ++i;
	}	  
    }
  
  
  // C-alpha statistic
  return score / sqrt( aux->variance ) ;      
    
}



void Pseq::Assoc::stat_burden( const VariantGroup & vars , 
			       Aux_prelim * pre ,
			       Aux_burden * aux , 
			       std::map<std::string,std::string> * output , 
			       bool original )
{

  const int n = vars.n_individuals();

  int cnt_a = 0, cnt_u = 0;    // basic counts
  int unq_a = 0;               // uniq case counts 
  double chi2 = 0;             // 2-by-2 X^2 statistics sum
  std::vector<int> ghits(n,0); // multi-hit instances
  int multi_a = 0 , multi_u = 0; // multiple-hit test counts

  for (int v=0; v<vars.size(); v++ ) 
    {
      
      int alta = 0 , altu = 0; 
      int refa = 0 , refu = 0;
      
      for (int i = 0; i < n; i++)
	{

	  // missing genotype?
	  if( vars.geno(v,i).null() ) continue;
	  
	  // get (permuted) phenotype
	  int j = g.perm.pos( i );	  
	  affType aff = vars.ind( j )->affected();
	  
	  // allele-count
	  if ( vars.geno(v,i).minor_allele( pre->refmin.find(v) 
					    == pre->refmin.end() ) )
	    {	      
	      if ( aff == CASE ) alta++; 
	      else if ( aff == CONTROL ) altu++;
	      if ( aux->mhit ) ghits[i]++;
	    }
	  else
	    {
	      if ( aff == CASE ) refa++; 
	      else if ( aff == CONTROL ) refu++;	
	    }	      
	}


      // 2x2 tables
      chi2 += Helper::chi2x2( alta , altu , refa , refu );
      
      // counting sites, not alleles (collapse 1+ --> 1)
      if ( aux->site_burden ) 
	{
	  alta = alta ? 1 : 0; 
	  altu = altu ? 1 : 0;
	}

      // case-burden
      cnt_a += alta;
      cnt_u += altu;
      
      // case-unique burden
      if ( altu == 0 ) unq_a += alta;
      


    } // next variant

  
  if ( aux->mhit ) 
    {
      for (int i=0; i < vars.n_individuals(); i++)
	{
	  if ( ghits[i] > 1 ) 
	    {	      
	      int j = g.perm.pos( i );
	      affType aff = vars.ind( j )->affected();	      
	      if ( aff == CASE ) ++multi_a;
	      else ++multi_u;	      
	    }
	}
    }

  // accumulate statistics
  aux->stat_vanilla = chi2;  
  aux->stat_burden = cnt_a - cnt_u;  
  aux->stat_uniq = unq_a; 
  aux->stat_mhit = multi_a - multi_u;

  // any output?
  if ( original ) 
    {      

      (*output)["BURDEN"] = Helper::int2str( cnt_a ) + "/" + Helper::int2str( cnt_u );
      (*output)["UNIQ"] = Helper::int2str( unq_a ) + "/0";
	
      if ( aux->mhit ) 
 	output->insert(make_pair( "MHIT" , 
 				 Helper::int2str( multi_a ) 
 				 + "/" + Helper::int2str( multi_u ) ) );
    }
  
  
  return;
}

  


void Pseq::Assoc::stat_fw_vt( const VariantGroup & vars , 
			      Aux_prelim * pre , 
			      Aux_fw_vt * aux , 
			      std::map<std::string,std::string> * output , 
			      bool original )
{
  
  const int n = vars.n_individuals();
  
  // calculate phenotypic mean (only for original)
  
  if ( original )  
    {
      aux->pmean = 0; 
      int n1 = 0;
      for ( int i = 0 ; i < n ; i++ )
	{
	  if ( vars.ind( i )->affected() == CASE ) 
	    {
	      ++aux->pmean;
	      ++n1;
	    }
	  else if ( vars.ind( i )->affected() == CONTROL ) 
	    {
	      ++n1;
	    }
	}
      aux->pmean /= (double)n1;
    }



  // combined function for FW and/or VT test
  
  // key = sample count
  // value = component of statistic for vt-test

  std::map<int,double> sumx;
  std::map<int,double> sum1;
  std::map<int,double> sum2;
  
  aux->stat_vt = aux->stat_fw = 0;

  int bestk = 0;
  
  std::map<int, std::set<int> >::iterator i = pre->carriers.begin();
  while ( i != pre->carriers.end() )
    {
      int j = g.perm.pos(  i->first  );	  
      affType aff = vars.ind( j )->affected();
      double ph = aff == CASE ? 1 : 0;
      
      std::set<int>::iterator k = i->second.begin();
      while ( k != i->second.end() )
	{
	  int na = vars(*k,i->first).minor_allele_count( pre->refmin.find(*k) 
							 == pre->refmin.end() );
	  
	  // FW-test
	  if ( aux->fw ) 
	    aux->stat_fw += ph * na * pre->fweights[ *k ];
	  
	  // VT-test
	  if ( aux->vt )
	    {
	      sumx[ pre->acounts[*k] ] += ph * na;
	      sum1[ pre->acounts[*k] ] += na;
	      sum2[ pre->acounts[*k] ] += na * na;	      	      
	    }

	  ++k; // next variant
	}
      ++i; // next carrier
    }
       

  // Find VT-threshold
  
  if ( aux->vt )
    {
      double tsumx = 0, tsum1 = 0, tsum2 = 0;
      
      std::map<int,double>::iterator z = sumx.begin();
      
      // Consider each possible threshold
      
      while ( z != sumx.end() )
	{
	  
	  int th = z->first;
	  
	  tsumx += z->second;
	  tsum1 += sum1[th];
	  tsum2 += sum2[th];
	  
	  double tscore = ( tsumx - tsum1 * aux->pmean ) / sqrt( tsum2 );
	  
	  if ( tscore > aux->stat_vt )
	    {
	      aux->stat_vt = tscore;
	      bestk = th; // and track threshold
	    }
	  ++z;
	}
      
      if ( original ) (*output)["VT"] = Helper::int2str( bestk );
    }
  
  return;

}


void Pseq::Assoc::stat_cancor( const VariantGroup & vars , 
			       Aux_prelim * pre ,
			       Aux_cancor * aux , 
			       std::map<std::string,std::string> * output , 
			       bool original )
{

  const int n = vars.n_individuals();

  Data::Matrix<double> PP( n , 1 );
  
  if ( original ) 
    {
      // populate P and G matrices

      for (int i=0; i < n; i++)
	aux->P(i,0) = (double)( vars.ind(i)->affected() == CASE );

      PP = aux->P;

      for (int v=0; v<vars.size(); v++)
	for (int i=0; i < n ; i++)
	  aux->G(i,v) = vars(v,i).null() ? 
	    0 : 
	    vars(v,i).minor_allele_count( pre->refmin.find(v) == pre->refmin.end() ); 
    }    
  else
    {
      // P and G already exist; permute P      
      for (int i=0; i<n; i++)
	PP(i,0) = aux->P( g.perm.pos(i) , 0 );      
    }
 
  double pvalue = 0;
  
  std::vector<double> cc = Statistics::canonical_correlation( PP , aux->G , &pvalue );
  
  if ( original ) (*output)["CANCOR"] = Helper::print( cc , false , false , "," );

  // either use 1st canonical correlation, or 1-p value for Bartletts (swap to X2 for Bartletts)
  aux->stat = 1 - pvalue;

}

