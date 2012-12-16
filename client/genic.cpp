#include "genic.h"
#include "util.h"
#include "plinkseq/defs.h"

#include <iostream>
#include <algorithm>
#include "util.h"

const double EPS = 1e-8;

extern GStore g;
extern Pseq::Util::Options args;

void   Pseq::Assoc::prelim( const VariantGroup & vars , Aux_prelim * aux , AuxGenic * data )  
{
  
  // track AuxGenic here too
  aux->paux = data;

  // Track observed min/max minor allele counts
  
  aux->minm = -1;
  aux->maxm = -1;
  
  // Frequency-weights
  aux->fweights.resize( vars.size() , 0 );  
  aux->acounts.resize( vars.size() , 0 );
  aux->altmin.resize( vars.size() , true );
  aux->maf.resize( vars.size() , 0 );

  // Generic weights
  aux->use_wgt = data->weights;
  if ( data->weights )
    {

      aux->wgt.resize( vars.size() , 0 );

      meta_index_t midx = MetaInformation<VarMeta>::field( data->weight_tag );
      
      if ( midx.mt != META_FLOAT && midx.len != 1 ) 
	Helper::halt( "expecting --weights tag to be a scalar Float" );
      
      for ( int v = 0 ; v < vars.size(); v++ )
	{
	  if ( vars(v).meta.has_field( data->weight_tag ) )
	    aux->wgt[v] = vars(v).meta.get1_double( data->weight_tag );
	}
    }
  
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
      aux->maf[ v ] = c / (double)c_tot; 
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
	{
	  if ( ( ! vars.geno(v,i).null() ) && vars.geno(v,i).minor_allele( aux->refmin.find(v) == aux->refmin.end() ) )
	    aux->carriers[i].insert( v );		  	  	      
	}
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
	  const int ac = vars.geno( *k , i1->first ).minor_allele_count( pre->refmin.find(*k) == pre->refmin.end() ) ;
	  if ( aff == CASE ) y[ *k ] += ac ;
	  n[ *k ] += ac ;	      
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
	  
	  if ( pre->mc[m] == 0 ) continue;
	  for ( int u = 0 ; u <= m ; u++ )
	    {	      
	      double s = ( u - m * aux->p_a );
	      s *= s;
	      s -= m * aux->p_ap_u;
	      s *= s; 	      
	      s *= Statistics::dbinom( u , m , aux->p_a );	      
	      t += s;
	    }	      
	  
	  aux->variance += pre->mc[m] * t;
	  
	}
    }


  // Frequency breakdown output

  if ( original ) 
    {
      (*output)["CALPHA"] = "";
      std::map<std::string,int>::iterator i = pre->mc_a.begin();
      while ( i != pre->mc_a.end() )
	{
	  if ( i != pre->mc_a.begin() ) (*output)["CALPHA"] += ";";
	  (*output)["CALPHA"] += i->first + "(" + Helper::int2str( i->second ) + ")";
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


  // for original data, represent the breakdown of variants in a nice manner
  Out * pdet = original ? &Out::stream( "assoc.det" ) : NULL ;
  
  // for original data, optionally output actual carriers
  Out * pcar = original &&  pre->paux->dump_carriers ? &Out::stream( "assoc.carriers" ) : NULL ;

  const int n = vars.n_individuals();

  // Given we are evaluating by permutation, just use a Z-test to get
  // 1-sided case vs. control comparison of (weighted) allele counts
  
  // Z = ( mA - mU ) / sqrt( SA^2/NA + SU^2/NU ) 

  std::vector<double> x(n,0);

  double cnt_a = 0, cnt_u = 0;      // basic counts
  double all_a = 0, all_u = 0;      // all non-missing allele-counts

  double unq_a = 0;                 // uniq case counts 
  double unq_u = 0;                 // likewise for controls

  double chi2 = 0;                  // 2-by-2 X^2 statistics sum

  // 'mhit' test is unweighted
  std::vector<int> ghits(n,0);      // multi-hit instances
  int multi_a = 0 , multi_u = 0;    // multiple-hit test counts

  for (int v=0; v < vars.size(); v++ ) 
    {
      
      double alta = 0 , altu = 0; 
      double refa = 0 , refu = 0;
      
      double w = pre->use_wgt ? pre->wgt[v] : 1;
      
      bool adir = pre->refmin.find(v) == pre->refmin.end();
      
      bool obsu = false; // this variant seen in any controls?
      bool obsa = false; // or cases

      for (int i = 0; i < n; i++)
	{
	  
	  // missing genotype?
	  if ( vars.geno(v,i).null() ) continue;
	  
	  // get (permuted) phenotype
	  int j = g.perm.pos( i );	  
	  affType aff = vars.ind( j )->affected();
	  
	  // allele-count (domainant model)
	  if ( vars.geno(v,i).minor_allele( adir ) )
	    {	      
	      
	      if      ( aff == CASE )    { alta += w; obsa = true; } 
	      else if ( aff == CONTROL ) { altu += w; obsu = true; }

	      // for Z-test
	      x[i] += w;
	      
	      if      ( aux->mhit ) ghits[i]++;
	      
	    }
	  else
	    {
	      if      ( aff == CASE )    refa += w; 
	      else if ( aff == CONTROL ) refu += w;	
	      
	    }	      
	}


      // 2x2 tables
      chi2 += Helper::chi2x2( alta , altu , refa , refu );
      
      // counting sites, not alleles (collapse 1+ --> 1)
      // TODO -- need to fix to make this weighted

      if ( aux->site_burden ) 
	{
	  Helper::halt("not implemented");
	  alta = alta ? w : 0; 
	  altu = altu ? w : 0;
	}

      
      // case-burden
      cnt_a += alta;
      cnt_u += altu;
      
      all_a += alta + refa;
      all_u += altu + refu;
      
      // case/control-unique burden 
      if ( obsa != obsu )
	{
	  unq_a += alta;
	  unq_u += altu;
	}
      

      // optional output, for real data

      if ( pdet ) 
	*pdet << vars.name() << "\t"
	      << vars(v) << "\t"
	      << vars(v).reference() << "/" << vars(v).alternate() << "\t"
	      << "W="<< w << "\t"
	      << alta << ":" << altu << "\n";
      
      
    } // next variant
  
  
      // carriers of ALT alleles; cycle through again
  if ( pcar )
    {
      
      for (int v=0; v < vars.size(); v++ ) 
	{
	  double w = pre->use_wgt ? pre->wgt[v] : 1;      
	  bool adir = pre->refmin.find(v) == pre->refmin.end();
	  
	  for (int i = 0; i < n; i++)
	    {
	      if ( vars.geno(v,i).null() ) continue;		  
	      affType aff = vars.ind( i )->affected();	  
	      
	      int ac = vars.geno(v,i).minor_allele( adir );
	      if ( vars.geno(v,i).minor_allele( adir ) )
		{	      
		  *pcar << vars.name() << "\t"
			<< vars(v) << "\t"
			<< vars.ind(i)->id() << "\t"
			<< w << "\t";
		  
		  if      ( aff == CASE )    *pcar << "CASE\t";
		  else if ( aff == CONTROL ) *pcar << "CONTROL\t";
		  else *pcar << ".\t";
		  
		  *pcar << ac << "\n";
		}
	    }
	}
    }



  // For Z-test, re-loop to get dominator

  double mA = 0 , mU = 0;
  int nA = 0 , nU = 0;
  
  for (int i=0;i<n;i++)
    {
      int j = g.perm.pos( i );
      affType aff = vars.ind( j )->affected();
      if      ( aff == CASE )    { nA++; mA += x[i]; }
      else if ( aff == CONTROL ) { nU++; mU += x[i]; }
    }
  
  if ( nA ) mA /= (double) nA;
  if ( nU ) mU /= (double) nU;

  double xxA = 0 , xxU = 0;

  for (int i=0;i<n;i++)
    {
      int j = g.perm.pos( i );
      affType aff = vars.ind( j )->affected();
      if      ( aff == CASE )    xxA += ( x[i] - mA ) * ( x[i] - mA );
      else if ( aff == CONTROL ) xxU += ( x[i] - mU ) * ( x[i] - mU );
    }
  
  if ( nA ) xxA /= (double)(nA-1);
  if ( nU ) xxU /= (double)(nU-1);
   
  //  std::cout << "DET " << vars.name() << " " << nA << " " << nU <<" " << mA << " " << mU << " " << xxA << " " << xxU << " " << ( mA - mU ) / sqrt( xxA/(double)nA + xxU/(double)nU ) << "\n";

  // simple multi-hit test
  
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

  //
  // accumulate statistics
  //

  aux->stat_vanilla = chi2;  


  // primary BURDEN test (Z-test, 1-sided)

  if ( nA == 0 || nU == 0 || mA < mU || ( xxA + xxU == 0 ) )
    aux->stat_burden = 0;
  else
    aux->stat_burden = ( mA - mU ) / sqrt( xxA/(double)nA + xxU/(double)nU ) ; 
  
  // track just case/control counts (dom. geno model)
  aux->alta = mA; aux->altu = mU; 
  aux->na = nA; aux->nu = nU;
  
  // use -log10(p) from Fisher's exact test; but 1-sided (so require higher alt-allele freq. in cases
  // if ( all_a == 0 || all_u == 0 || cnt_a / (double)all_a <= cnt_u / (double)all_u ) 
  //   aux->stat_burden = 0; 
  // else 
  //   aux->stat_burden = Helper::chi2x2( cnt_a , cnt_u , all_a - cnt_a , all_u - cnt_u );
  
  //  std::cout << "stats = " << aux->stat_burden << "\t" << cnt_a << " " << cnt_u << " " << all_a - cnt_a << " " << all_u - cnt_u << "\n";

  // uniq-burden test

  if ( all_a == 0 || all_u == 0 || unq_a / (double)all_a <= unq_u / (double)all_u ) 
    aux->stat_uniq = 0; 
  else 
    aux->stat_uniq = Helper::chi2x2( unq_a , unq_u , all_a - unq_a , all_u - unq_u );    


  // MHIT test (not really used now)

  aux->stat_mhit = multi_a - multi_u;
  

  // any output?

  if ( original ) 
    {      

      (*output)["BURDEN"] = Helper::dbl2str( cnt_a ) + "/" + Helper::dbl2str( cnt_u );
      
      (*output)["UNIQ"] = Helper::dbl2str( unq_a ) + "/0";
	
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


