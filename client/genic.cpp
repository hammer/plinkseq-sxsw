#include "genic.h"
#include "util.h"

#include <iostream>

extern GStore g;
extern Pseq::Util::Options args;

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
      //(*output)["CALPHA"] = "Z=" + Helper::dbl2str( score / sqrt( aux->variance ) )  ;
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

double
Pseq::Assoc::stat_two_hit(const VariantGroup & vars, Aux_prelim * pre, Aux_two_hit * aux, std::map<std::string, std::string> * output, bool original, std::map<std::string, int> var_class, double prev)
{

  const int n = vars.n_individuals();
  /*  double prev = .006;
  if ( args.has( "prev" ) )
    prev =  Helper::str2dbl(args.as_string( "prev" ));
  */

  // see stat_cancor above - you may be able to just copy that code for P and G
  // populate P and G matrices

  int ncases = 0;
  int ncontrols = 0;

  for (int i = 0; i < n; i++)
    {
      if(vars.ind(i)->affected() == CASE)
        ncases++;
      else
        ncontrols++;
    }

  std::vector<double> hets;
  std::vector<std::string> ann;
  int mat[3][3];
  int ahom = 0; int uhom = 0;
  int achet = 0; int uchet = 0;
  int ahet = 0; int uhet = 0;

  int ahomi = 0; int uhomi = 0;
  int acheti = 0; int ucheti = 0;
  int aheti = 0; int uheti = 0;
  

  for (int i = 0; i < n; i++)
    {
      // permuted individual index i->j
      int j = g.perm.pos(i);

      std::string id = vars.ind( j )->id();
      ahomi = uhomi = acheti = ucheti = aheti = uheti = 0;
      
      hets.clear();
      ann.clear();
      for( int j = 0; j < 3; j++ )
        for( int k = 0; k < 3; k++ )
          mat[j][k] = 0;
      
      for (int v = 0; v < vars.size(); v++){
        double d = vars(v, i).null() ? 0 : vars(v, i).minor_allele_count(pre->refmin.find(v) == pre->refmin.end());
	std::vector<int> pl;
	std::vector<int> ad;
        double ab  = 0;

        const Genotype & g = vars.geno(v,i);
        
        if ( vars.geno(v,i).meta.has_field( "AD" ) )
          {
            ad = vars.geno(v,i).meta.get_int( "AD" );   
            ab = (ad[0] * 1.0) / (ad[0] + ad[1]);
          }

        // count non-ref homozygotes
	std::string annot = "";

        if( vars(v).consensus.meta.has_field( "transcript" ) && vars(v).consensus.meta.has_field( "func" )){
	  //std::cout << "c\n";
	  std::vector<std::string> func = vars(v).consensus.meta.get_string( "func" );
	  std::vector<std::string> transcript = vars(v).consensus.meta.get_string( "transcript");
	  std::vector<std::string> func_split;
	  std::vector<std::string> trans_split;

	  std::string::size_type i1 = 0;
	  std::string::size_type i2 = transcript[0].find(",");
          while (i2 != std::string::npos) {
            trans_split.push_back(transcript[0].substr(i1, i2-i1));
            i1 = ++i2;
            i2 = transcript[0].find(",", i2);
            if (i2 == std::string::npos)
              trans_split.push_back(transcript[0].substr(i1, transcript[0].length( )));
          }
          
          i1 = 0;
          i2 = func[0].find(",");
          while (i2 != std::string::npos) {
            func_split.push_back(func[0].substr(i1, i2-i1));
            i1 = ++i2;
            i2 = func[0].find(",", i2);
            if (i2 == std::string::npos)
              func_split.push_back(func[0].substr(i1, func[0].length( )));
          }

	  std::string annot1 = "";
          
          if( trans_split.size() == 0 )
            annot = func[0];
          else{
            for( int ti = 0; ti < trans_split.size(); ti++)
              if( trans_split[ti].compare(vars.name()) == 0 ){
                if( annot.compare("") == 0 ){
                  annot = func_split[ti];
                }
                else{
		  std::string test = func_split[ti].substr(0,7);
                  if( test.compare("esplice") == 0  && annot.compare("nonsense") != 0 )
                    annot = func_split[ti];             
                }
              }
          }
        }

	//        bool pass = false;
	//        if( annot.compare("esplice3") == 0 || annot.compare("esplice5") == 0 || annot.compare("nonsense") == 0 ) // || annot.compare("missense") == 0 )
	//          pass = true;


	bool pass = true;
	if( var_class.size() > 0 ){
	  if( var_class.count(annot) > 0 )
	    pass = true;
	  else
	    pass = false;
	}

        if( d == 2 && ab < .1 && pass )
          {         
            if( vars.ind( j )->affected() == CASE){
              if(original)
		std::cout << "HOM: " << vars(v).chromosome() << ":" << vars(v).position() << " " << id << " case " << annot << "\n"; 
              ahomi++;
            }
            else{
              if(original)
		std::cout << "HOM: " << vars(v).chromosome() << ":" << vars(v).position() << " " << id << " control " << annot << "\n";
              uhomi++;
            }
          }

        // found het, store for later
        if( d == 1 && ab > .3 && ab < .7 && pass){
          hets.push_back(v);
          ann.push_back(annot);
          
          if(vars.ind(j)->affected() == CASE){
            aheti++;
          }
          else{
            uheti++;
          }
        }
      }

      if( hets.size() > 0 ){
        for(int z = 0; z < hets.size()-1; z++){
          for(int k = z+1; k < hets.size(); k++){

            // fill in genotype matrix to test two hits
            for (int l = 0; l < n; l++){
              int var1 = (int) vars(hets[z], l).null() ? 0 : vars(hets[z], l).minor_allele_count(pre->refmin.find(hets[z]) == pre->refmin.end());
              int var2 = (int) vars(hets[k], l).null() ? 0 : vars(hets[k], l).minor_allele_count(pre->refmin.find(hets[k]) == pre->refmin.end());         
              mat[var1][var2]++;
            }
            
            // test two hit
            int rest = mat[0][1] + mat[0][2] + mat[1][0] + mat[1][2] + mat[2][0] + mat[2][1] + mat[2][2];
            
            if( ((mat[0][1] + mat[0][2]) > 0 && (mat[1][0] + mat[2][0]) > 0 && (mat[1][2] + mat[2][1] + mat[2][2]) == 0) || (mat[1][1] == 1 && rest == 0)){
              if(vars.ind(j)->affected() == CASE){
                if(original)
		  std::cout << "CHET: " << vars(hets[z]).chromosome() << ":" << vars(hets[z]).position() << "-" << vars(hets[k]).chromosome() << ":" << vars(hets[k]).position() << " " << id << " case " << ann[z] << "-" << ann[k] << "\n";
                acheti++;
              }
              else{
                if(original)
		  std::cout << "CHET: " << vars(hets[z]).chromosome() << ":" << vars(hets[z]).position() << "-" << vars(hets[k]).chromosome() << ":" << vars(hets[k]).position() << " " << id << " control  " << ann[z] << "-" << ann[k] << "\n";
                ucheti++;
              }
            }   
          }
        }
      }
      
      // count individual instances
      /*ahom += ahomi;
    uhom += uhomi;
    achet += acheti;
    uchet += ucheti;
    ahet += aheti;
    uhet += uheti;
      */

      // count total number of cases and controls
      if( ahomi >0 )
	ahom++;
      if( uhomi >0 )
	uhom++;
      if( acheti >0 )
	achet++;
      if( ucheti >0 )
	uchet++;
      if( aheti >0 )
	ahet++;
      if( uheti >0 )
	uhet++; 
    }
  
  int case_allele = (2 * ncases) + ahet;
  
  double total_hets = ahet + uhet;
  double a = 2 * n;
  double b = -1 * 2 * n;
  double c = total_hets;
  double descrim = (b*b) - (4*a*c);

  double f = 1;  
  if( descrim >= 0 )
    f = ((-1*b) - sqrt(descrim)) / (2 * a);
  else
    f = ((-1*b) - sqrt(-1 * descrim)) / (2 * a);

  double arec = ahom + achet;
  double urec = uhom + uchet;
  
  double pAA = f * f;
  double pAA_1 = urec / ncontrols;
  if ( pAA < pAA_1) {
    pAA = pAA_1;
  }
  if ( pAA == 0 ) {
    pAA = 0.0001;
  }
  
  double l0 = (arec * log(pAA)) + (ncases - arec) * log( 1.0 - pAA );
  l0 += (urec * log(pAA)) + (ncontrols - urec) * log( 1.0 - pAA );
  double l1 = l0;
  
  double denom1 = 0;
  double this_l1 = 0;
  double x = 0;
  double pAA_nodis = 0;
  double grrmax = 0;


  for (double gr = 1.01; gr <= 50; gr += 0.01) {
    
    denom1 = gr * pAA + ( 1 - pAA );
    
    this_l1 = arec * log(gr*pAA/denom1) + (ncases-arec)*log(1.0 - (gr*pAA/denom1));
    x = prev/denom1;
    pAA_nodis = (1.0 - gr*x)*pAA / ( 1.0 - prev );
    this_l1 += urec*log(pAA_nodis) + (ncontrols-urec)*log(1.0 - (pAA_nodis));
    
    if (this_l1 > l1) { 
      l1 = this_l1; 
      grrmax = gr; 
    }
  }

  double diff = (l1-l0) / log(10);
  double diff2 = l1-l0;
  
  double chisqVal = 2 * diff2;

  double pvalue = Statistics::chi2_prob(chisqVal, 1);

  // set the statistic on the Aux_two_hit struct.
  // gseq will read this value each permutation (this also includes adaptive permutation)

  if (original)
    {
      (*output)["TWO-HIT"] = "P=" + Helper::flt2str( pvalue ) + ";AF=" + Helper::flt2str(f) + ";CASES=" + Helper::dbl2str(ahom) + "," + Helper::dbl2str(achet) + "," + Helper::dbl2str(ahet) + "," + Helper::dbl2str(ncases) + ";CONTROLS=" + Helper::dbl2str(uhom) + "," + Helper::dbl2str(uchet) + "," + Helper::dbl2str(uhet) + "," + Helper::dbl2str(ncontrols);
    }
  
  return chisqVal;

   
}

