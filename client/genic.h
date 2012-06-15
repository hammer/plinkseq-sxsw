#ifndef __PSEQ_GENIC_H__
#define __PSEQ_GENIC_H__

#include <map>
#include "pseq.h"
#include <cmath>

class VariantGroup;
class SKAT_param;

namespace Pseq 
{
  
  namespace Assoc
  {
    
    struct Aux_prelim;
    struct Aux_burden;
    struct Aux_calpha;
    struct Aux_fw_vt;
    struct Aux_cancor;
    struct Aux_hoffman_witte;
    struct Aux_kbac;
    struct Aux_two_hit;
    struct Aux_skat;

    void  prelim( const VariantGroup & vars , Aux_prelim * aux );
    
    void stat_burden( const VariantGroup & vars , 
		      Aux_prelim * pre , 
		      Aux_burden * aux , 
		      std::map<std::string,std::string> * output , 
		      bool original );
    
    double stat_calpha( const VariantGroup & , 
			Aux_prelim * , 
			Aux_calpha * , 
			std::map<std::string,std::string> * ,
			bool );
        
    void stat_fw_vt( const VariantGroup & , 
		     Aux_prelim * , 
		     Aux_fw_vt *  , 
		     std::map<std::string,std::string> *  , 
		     bool original );
    
    void stat_cancor( const VariantGroup & , 
		      Aux_prelim * , 
		      Aux_cancor * , 
		      std::map<std::string,std::string> *  , 
		      bool original );
    
    double stat_hoffman_witte( const VariantGroup & , 
			       Aux_hoffman_witte * , 
			       std::map<std::string,std::string> *  , 
			       bool original );
    
    
    double stat_kbac( const VariantGroup & , 
		      Aux_prelim * , 
		      Aux_kbac * , 
		      std::map<std::string,std::string> *  , 
		      bool original );

    double stat_two_hit( const VariantGroup & ,
                         Aux_prelim * ,
                         Aux_two_hit * ,
                         std::map<std::string,std::string> *  ,
                         bool original, std::map< std::string, int >, std::map< std::string, int >, double, bool );
    
    
    double stat_skat( const VariantGroup & , 
		      Aux_prelim * , 
		      Aux_skat * , 
		      std::map<std::string,std::string> *  , 
		      bool original );
    
    struct AuxGenic 
    {
      
      AuxGenic() 
      { 
	g=NULL;
	rseed=0; 
	fix_null_genotypes = true;
	dump_stats_matrix = false;
	show_info = false; 
	show_midbp = false;
	vanilla = true;
	burden = true;
	site_burden = false;
	uniq = true;
	mhit = false;
	vt = false; 
	fw = false;
	calpha = false;
	cancor = false;
	hoffman_witte = false;
	kbac = false;
	two_hit = false;
	skat = false;
      }
      

      int n_tests() const 
      {
	return vanilla 
	  + burden 
	  + uniq 
	  + mhit 
	  + vt 
	  + fw 
	  + calpha 
	  + cancor 
	  + hoffman_witte
	  + kbac
	  + two_hit
	  + skat 
	  + skato;
      }
      
      GStore * g;
      long int rseed;
      bool fix_null_genotypes;
      bool show_info;
      bool show_midbp;
      bool dump_stats_matrix;
      
      // genic tests      
      bool mhit;
      bool vanilla;
      bool burden;
      bool site_burden; // modifier, not a test
      bool uniq;
      bool fw;
      bool vt;
      bool calpha;
      bool cancor;
      bool hoffman_witte;
      bool kbac;
      bool two_hit;
      bool skat;
      bool skato; // optimal SKAT
    };
 
    
    struct Aux_prelim {
      
      int minm;
      int maxm;

      bool dichot;
      int n_a;
      int n_u;
      int n_t;
	
      std::set<int> refmin;            // track alleles at which reference is minor 
      std::vector<bool> altmin;        // same info but as vector

      std::map<int,int> mc;            // key=#minor alleles,value=#variants   
      std::map<std::string,int> mc_a;  // count in affecteds
      
      std::vector<double> wgt;         // generic weights

      std::vector<double> maf;         // MAF (of minor allele always)
      std::vector<double> fweights;    // frequency weights
      std::vector<int> acounts;        // allele counts (used in VT)
      
      std::map<int,std::set<int> > carriers; 
      
    };
    

    struct Aux_burden {
      
      Aux_burden( bool vanilla , 
		  bool burden , 		  
		  bool uniq , 
		  bool mhit , 
		  bool site_burden  )
	: vanilla(vanilla) , burden(burden) , 
	  uniq(uniq) , mhit(mhit) , site_burden(false) 
      { 
      } 
      
      bool vanilla;
      bool burden;
      bool site_burden;
      bool uniq;
      bool mhit;

      double stat_vanilla;
      double stat_burden;
      double stat_uniq;
      double stat_mhit;

    };


    struct Aux_fw_vt { 
      Aux_fw_vt(bool fw , bool vt) : fw(fw), vt(vt) { } 
      double pmean; // phenotypic mean
      bool fw;
      bool vt;
      double stat_fw;
      double stat_vt;
    };
    

    struct Aux_calpha {
      double p_a;
      double p_ap_u;      
      double variance;
    };


    struct Aux_cancor { 
      Aux_cancor(int np, int ng, int ni) : P(ni,np) , G(ni,ng) { } 
      double stat;
      Data::Matrix<double> P; 
      Data::Matrix<double> G;
    };
    

    struct Aux_hoffman_witte 
    {      
      
      Aux_hoffman_witte( bool , const VariantGroup & vars, Aux_prelim * p );      
      
      std::vector<int> piece_begin; // cutpoint begin points
      std::vector<int> piece_end;   // cutpoint end points
      double tbar;                  // mean of the trait
      std::vector<double> xbar;     // mean of the genotypes      
      Aux_prelim * aux_prelim;
      bool dichot;
      
      double calc_stat( const VariantGroup & , 
			const std::vector<int> & , 
			const std::vector<double> & , 
			const std::vector<double> & );
      
      void clear() 
      {
	piece_begin.clear();
	piece_end.clear();
	tbar = 0.0;
	xbar.clear();
      }      
      
    };

    
    struct Aux_kbac 
    {
      static double lnfact_table[]; 
      double gw_ln_factorial( const double ) const;
      double gw_lnchoose( const double , const double ) const;
      double gw_hypergeometric_pmf(const unsigned int , const unsigned int , const unsigned int , const unsigned int ) const;
      double gw_hypergeometric_cmf(const unsigned int , const unsigned int , const unsigned int , const unsigned int ) const;      
    };


    //
    // SKAT test
    //
    
    struct Aux_skat 
    {    
      
      Aux_skat() 
      { 
	optimal = false;
      }
      
      // when first run, precalculate (only ever once) the 
      // adjust phenotypic expected values for each individual.

      static bool precalculated;

      // note -- having both 'y' and 'Y' is redundant
      //     Y only exists because don't have a vector * matrix operation...

      static Data::Vector<double> y; // phenotype --> modified (y_i - u_i)
      static Data::Matrix<double> Y; // phenotype --> modified (y_i - u_i) (KLUDGE for matrix mult)
      static Data::Vector<double> u; // phenotype --> modified (u_i)
      
      // original values (for permutation)
      static Data::Vector<double> y_orig; // phenotype --> modified (y_i - u_i)
      static Data::Matrix<double> Y_orig; // phenotype --> modified (y_i - u_i) (KLUDGE for matrix mult)
      static Data::Vector<double> u_orig; // phenotype variance --> modified (u_i)

      static Data::Matrix<double> X; // covariates
      static std::vector<bool> mask; // included? (i.e. non-missing pheno/covar?)
      static int n_actual; 
      static bool logistic_model;   

      static std::vector<double> rho; // grid of correlation values
      static int rho_est;
      static int nr;

      // main functions
      static void fit_null();
      
      void populate_G( const VariantGroup & , Aux_prelim * );

      void populate_K();

      double calculate_Q( Data::Matrix<double> * );

      double calculate_optimal_Q( Data::Matrix<double> * );
      Data::Vector<double> sub_optimal_get_Q( const Data::Matrix<double> & );
      double sub_optimal_get_P( const Data::Vector<double> & , const Data::Matrix<double> & );
      Data::Vector<double> get_lambda( const Data::Matrix<double> & );
      void get_optimal_param( const Data::Matrix<double> & Z1 ,  
 			      double * muQ, double * varQ, double * kerQ,  
 			      double * varRemain, double * df,  
 			      Data::Vector<double> * tau, Data::Vector<double> * lambda ); 

      Data::Vector<double> sub_optimal_P_foreach_Q( const Data::Vector<double> & Q , 
						    const std::vector<std::vector<double> > & lambda );
						   
      
      double sub_optimal_P_Davies( const Data::Vector<double> & pminq , const SKAT_param & param );
      double sub_optimal_P_Liu( const Data::Vector<double> & pminq , const SKAT_param & param );
      
      double calculate_pvalue( double , Data::Matrix<double> & );

      double calculate_optimal_pvalue( double , Data::Matrix<double> & );
      
      double get_liu_pval( double , Data::Matrix<double> & );
      
      void get_liu_param( bool mod , // use modified form?
			  double c1 , double c2 , double c3 , double c4 , 
			  double * muX , double * sigmaX , 
			  double * muQ , double * sigmaQ , 
			  double * l , double * d ) ;

      
      // SKAT-O specific functions
      
      void set_optimal_mode( const bool b = true )
      {
	optimal = b;
      }
      
      bool optimal_mode() const { return optimal; } 
      
      void set_optimal_rcorr( const std::vector<double> & r )
      {
	rho = r;
	// for numerical reasons, following original SKAT R implementation
	for (int i=0;i<rho.size();i++) if ( rho[i] > 0.999 ) rho[i] = 0.999; 
      }

      void set_optimal_rcorr()
      {
	// default 0.00 , 1.00 grid
	rho.clear();
	for (int i=0;i<=10;i++) rho.push_back( (double)i/(double)10.0 );
	rho[10] = 0.999;
      }
            
      
      // run as SKAT-O
      bool optimal; 

      //
      static bool has_covar;
      static std::vector<std::string> covars;

      static bool has_weights;
      static std::string weights;

      // Beta-MAF weight parameters      
      static bool use_freq_weights;
      static double a1;
      static double a2;

      Data::Vector<double> w; // weights
      Data::Matrix<double> G; // genotype-data
      Data::Matrix<double> K; // kernel
      
      // use slot as return value for asymptotic p-value from each test
      double returned_pvalue;
  };

    

    
    struct Aux_two_hit 
    {
    Aux_two_hit(int np, int ng, int ni) : P(ni,np) , G(ni,ng), LD(ng, ng) { }
      double stat;
      Data::Matrix<double> P; // matrix of phenotypes (1 x num individuals)
      Data::Matrix<double> G; // genotypes of each individual
      Data::Matrix<double> LD; // LD matrix (num genotypes x num genotypes)
    };


  }  
  
}


#endif
