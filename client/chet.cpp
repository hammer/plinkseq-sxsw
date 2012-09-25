#include "genic.h"
#include <iostream>
#include <functional>
#include <cmath>
#include <cstdlib>
#include <cstddef>
#include <iomanip>
#include <list>
#include <algorithm>
#include "util.h"


extern GStore g;
extern Pseq::Util::Options args;
using namespace std;

double Pseq::Assoc::Aux_two_hit::prev = PLINKSeq::DEFAULT_PREV();
std::map< std::string, int > Pseq::Assoc::Aux_two_hit::func_inc;
std::map< std::string, int > Pseq::Assoc::Aux_two_hit::func_exc;

int Pseq::Assoc::Aux_two_hit::ncases = 0;
int Pseq::Assoc::Aux_two_hit::ncontrols = 0;
int Pseq::Assoc::Aux_two_hit::nmissing = 0;

bool Pseq::Assoc::Aux_two_hit::singles = false;
bool Pseq::Assoc::Aux_two_hit::mhit = false;

void Pseq::Assoc::Aux_two_hit::initialize(){

  prev = PLINKSeq::DEFAULT_PREV();
  singles = false;
  mhit = false;

  if ( args.has( "prev" ) )
    prev =  Helper::str2dbl(args.as_string( "prev" ));

  if ( args.has( "singles" ) )
    singles = true;

  if ( args.has( "mhit" ) )
    mhit = true;


  if( args.has( "func-inc" ) ){
    std::vector<std::string> inc = args.as_string_vector( "func-inc" );
    for( int i = 0; i < inc.size(); i++ )
      func_inc[inc[i]] = i;
  }
  
  if( args.has( "func-exc" ) ){
    std::vector<std::string> inc = args.as_string_vector( "func-exc" );
    for( int i = 0; i < inc.size(); i++ )
      func_exc[inc[i]] = i;
  }

  
}


double Pseq::Assoc::stat_two_hit( const VariantGroup & vars, 
				  Aux_prelim * pre, 
				  Aux_two_hit * aux, 
				  std::map<std::string, 
				  std::string> * output, 
				  bool original )
{
  Out & pout = Out::stream( "twohit.vars" );

  if( aux->ncases == 0 ) {

    const int n = vars.n_individuals();

    for (int i = 0; i < n; i++)
      {
	if ( vars.ind(i)->affected() == CASE )
	  aux->ncases++;
	if ( vars.ind(i)->affected() == CONTROL )
	  aux->ncontrols++;
	if ( vars.ind(i)->missing() )
	  aux->nmissing++;
      }
  }
  
  std::vector<double> hets;
  std::vector<std::string> ann;
  std::map<std::string, int> valid_chet;

  int mat[3][3];
  int ahom = 0; int uhom = 0;
  int achet = 0; int uchet = 0;
  int ahet = 0; int uhet = 0;

  int ahomi = 0; int uhomi = 0;
  int acheti = 0; int ucheti = 0;
  int aheti = 0; int uheti = 0;

  const int n = vars.n_individuals();
  
  
  for (int i = 0; i < n; i++)
    {
      
      // permuted individual index i->j
      int j = g.perm.pos(i);
      
      ahomi = uhomi = acheti = ucheti = aheti = uheti = 0;
      
      hets.clear();
      ann.clear();
      
      for( int j = 0; j < 3; j++ )
        for( int k = 0; k < 3; k++ )
          mat[j][k] = 0;

      
      for (int v = 0; v < vars.size(); v++)
	{
	  
	  const Genotype & gt = vars.geno(v,i);
	  double d = gt.null() ? 0 : gt.minor_allele_count(pre->refmin.find(v) == pre->refmin.end());
	  
	  std::vector<int> pl;
	  std::vector<int> ad;
	  double ab  = -1;
	  
	  if ( vars.geno(v,i).meta.has_field( PLINKSeq::DEFAULT_AD() ) )  
	    {
	      ad = vars.geno(v,i).meta.get_int( PLINKSeq::DEFAULT_AD() );   
	      ab = (ad[0] * 1.0) / (ad[0] + ad[1]);
	    }
	  else
	    if( d != 0 )
	      plog.warn( "NO ALLELIC DEPTH INFORMATION FOR VARIANT LEVEL QC");

	  // count non-ref homozygotes
	  std::vector< std::string > annot;
	  
	  if( vars(v).meta.has_field( PLINKSeq::DEFAULT_TRANS() ) && vars(v).meta.has_field( PLINKSeq::DEFAULT_FUNC() )){

	    std::vector<std::string> func = vars(v).meta.get_string( PLINKSeq::DEFAULT_FUNC() );
	    std::vector<std::string> transcript = vars(v).meta.get_string( PLINKSeq::DEFAULT_TRANS() );

	    //  fast string tokenizer:
	    int n = 0;
	    Helper::char_tok tok_trans( transcript[0] , &n , ',' );
            
	    n = 0;
	    Helper::char_tok tok_func( func[0] , &n , ',' );
	    if( tok_trans.size() == 0 )
              annot.push_back(func[0]);
	    else{
	      for( int ti = 0; ti < tok_trans.size(); ti++)
                if( vars.name().compare(tok_trans[ti]) == 0 )
                  annot.push_back(tok_func[ti]);      
	    }
	  }

	  std::string annot1 = "";
	  int cnt = 0;
	  
	  // only include defined annotations unless none given then include all
	  bool pass = false;	
       	  if (aux->func_inc.size() == 0)
	    pass = true;
	  

	  for( int ti = 0; ti < annot.size(); ti++){
	    if( cnt > 0 )
	      annot1 += "/";
	    annot1 += annot[ti];
	    cnt++;
	    if( aux->func_inc.count(annot[ti]) > 0 || aux->func_inc.size() == 0 )
	      pass = true;
	    
	    if( aux->func_exc.count(annot[ti]) > 0 ){
	      pass = false;
	      break;
	    } 
	  }

	  if( annot1.length() == 0 )
	    annot1 = ".";
	  if( d == 2 && (ab < PLINKSeq::DEFAULT_AB_HOMMAX() || ab == -1) && pass )
	    {         
	      if( vars.ind( j )->affected() == CASE){
		if(original)
		  pout << "HOM\t" \
		       << vars.name() << "\t" \
		       << g.locdb.alias( vars.name() , false ) << "\t" \
		       << vars(v) << "\t" \
		       << vars.ind( j )->id() \
		       << "\tcase\t" \
		       << annot1 << "\n"; 
		ahomi++;
	      } 
	      if ( vars.ind(j)->affected() == CONTROL ){
		if(original)
		  pout << "HOM\t" \
		       << vars.name() << "\t" \
                       << g.locdb.alias( vars.name() , false ) << "\t" \
		       << vars(v) << "\t" \
		       << vars.ind( j )->id() \
		       << "\tcontrol\t" \
		       << annot1 << "\n";
		uhomi++;
	      }
	    }
	  // found het, store for later
	  if( d == 1 && ((ab > PLINKSeq::DEFAULT_AB_HETMIN() && ab < PLINKSeq::DEFAULT_AB_HETMAX()) || ab == -1) && pass){
	    hets.push_back(v);	  	  
	    ann.push_back(annot1);
	    if ( vars.ind(j)->affected() == CASE )
	      aheti++;
	    if ( vars.ind(j)->affected() == CONTROL )
	      uheti++;
	  }
	}

      if ( hets.size() > 0 ) {
        for (int z = 0; z < hets.size()-1; z++){
          for (int k = z+1; k < hets.size(); k++){
	    
	  std:string pair = Helper::int2str(vars(hets[z]).position()) + "," + Helper::int2str(vars(hets[k]).position()); 
	    
	    // clear matrix
	    
	    //perform population check only when pair hasn't been seen before
	    //	    bool chet = false;
	    if (  valid_chet.count(pair) == 0 ){
	      //clear matrix
	      for (int u = 0; u < 3; u++ )
		for (int v = 0; v < 3; v++ )
		  mat[u][v] = 0;
	      
	      
	      // fill in genotype matrix to test two hits
	      for (int l = 0; l < n; l++){
		int var1 = (int) vars(hets[z], l).null() ? 0 : vars(hets[z], l).minor_allele_count(pre->refmin.find(hets[z]) == pre->refmin.end());
		int var2 = (int) vars(hets[k], l).null() ? 0 : vars(hets[k], l).minor_allele_count(pre->refmin.find(hets[k]) == pre->refmin.end());         
		mat[var1][var2]++;
	      }
	      
	      // test two hit
	      int rest = mat[0][1] + mat[0][2] + mat[1][0] + mat[1][2] + mat[2][0] + mat[2][1] + mat[2][2];            
	      if( ((mat[0][1] + mat[0][2]) > 0 && (mat[1][0] + mat[2][0]) > 0 && ( mat[1][2] + mat[2][1] + mat[2][2] == 0 ) ) || (mat[1][1] == 1 && rest == 0 && aux->singles) || aux->mhit )
		//	chet = true;
		valid_chet[pair] = 1;
	      else
		valid_chet[pair] = 0;
	      }
	    
	    if ( valid_chet[pair] == 1 ){ // chet || ( valid_chet.count(pair) > 0 && valid_chet[pair] == 1 ) ){

	      //sort two annotations for consistancy
	      std::vector<std::string> vtmp;
	      if( ann[z].length() > 0 )
		vtmp.push_back(ann[z]);
	      else
		vtmp.push_back(".");
	      if( ann[k].length() > 0 )
		vtmp.push_back(ann[k]);
	      else
		vtmp.push_back(".");

	      sort(vtmp.begin(), vtmp.end());
	      

              if ( vars.ind(j)->affected() == CASE )
		{
		  if ( original )
		    pout << "CHET\t"		\
			 << vars.name() << "\t" \
			 << g.locdb.alias( vars.name() , false ) << "\t" \
			 << vars(hets[z]) << "," << vars(hets[k]) << "\t" \
			 << vars.ind( j )->id() \
			 << "\tcase\t" \
			 << vtmp[0] << "," << vtmp[1] << "\n";
		  acheti++;
		}
              if ( vars.ind(j)->affected() == CONTROL )
		{
		  if( original )
		    pout << "CHET\t" \
		         << vars.name() << "\t" \
                         << g.locdb.alias( vars.name() , false ) << "\t" \
			 << vars(hets[z]) << "," << vars(hets[k]) << "\t" \
			 << vars.ind( j )->id() \
			 << "\tcontrol\t" \
			 << vtmp[0] << "," << vtmp[1] << "\n";
		  ucheti++;
		}
            }
          }
        }
      }

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
  


  // skip chrX and Y 
  //  haploid ??? chrX ?? 

  double total_hets = ahet + uhet;
  double a = 2 * aux->ncontrols;
  double arec = ahom + achet;
  double urec = uhom + uchet;
  double f = (( 2 * urec ) + uhet) / a;

  //double f = (uhet + (2*urec))/(2*aux->ncontrols);

  double pAA = f * f;
  double pAA_1 = urec / aux->ncontrols;
  if ( pAA < pAA_1) {
    pAA = pAA_1;
  }
  if ( pAA == 0 ) {
    pAA = 0.0001;
  }
  
  double l0 = (arec * log(pAA)) + (aux->ncases - arec) * log( 1.0 - pAA );
  l0 += (urec * log(pAA)) + (aux->ncontrols - urec) * log( 1.0 - pAA );
  double l1 = l0;
  
  double denom1 = 0;
  double this_l1 = 0;
  double x = 0;
  double pAA_nodis = 0;
  double grrmax = 0;


  for (double gr = 1.01; gr <= 50; gr += 0.01) {
    
    denom1 = gr * pAA + ( 1 - pAA );
    
    this_l1 = arec * log(gr*pAA/denom1) + (aux->ncases-arec)*log(1.0 - (gr*pAA/denom1));
    x = aux->prev/denom1;
    pAA_nodis = (1.0 - gr*x)*pAA / ( 1.0 - aux->prev );
    this_l1 += urec*log(pAA_nodis) + (aux->ncontrols-urec)*log(1.0 - (pAA_nodis));
    
    if (this_l1 > l1) { 
      l1 = this_l1; 
      grrmax = gr; 
    }
  }

  double diff = (l1-l0) / log(10);
  double diff2 = l1-l0;
  
  double chisqVal = 2 * diff2;

  double pvalue = Statistics::chi2_prob(chisqVal, 1);
  aux->returned_pvalue = pvalue < 0 ? 1.00 : pvalue ;
  

  // set the statistic on the Aux_two_hit struct.
  // gseq will read this value each permutation (this also includes adaptive permutation)

  if (original)
    {
      (*output)["TWO-HIT"] = "P=" + Helper::dbl2str( pvalue ) 
	+ ";AF=" + Helper::dbl2str(f) 
	+ ";CASES=" + Helper::int2str(ahom) + "," + Helper::int2str(achet) + "," + Helper::int2str(ahet) + "," + Helper::int2str(aux->ncases) 
	+ ";CONTROLS=" + Helper::int2str(uhom) + "," + Helper::int2str(uchet) + "," + Helper::int2str(uhet) + "," + Helper::int2str(aux->ncontrols);
    }
  return chisqVal;   
}

