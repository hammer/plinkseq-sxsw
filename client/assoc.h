#ifndef __PSEQ_ASSOC_H__
#define __PSEQ_ASSOC_H__

#include "matrix.h"
#include "func.h"

#include <map>
#include <set>

class Mask;
class VariantGroup;
class Variant;
class GStore;

namespace Pseq 
{

  namespace Util {
    class Options;
  }

  namespace Assoc
    {
      
      struct Aux_vassoc_options;
      struct Aux_glm;

      bool variant_assoc_test( Mask & , Aux_vassoc_options &  , const Pseq::Util::Options & options );

      bool glm_assoc_test( Mask & , Aux_glm &  );
      
      bool set_assoc_test( Mask & , const Pseq::Util::ArgMap & , const Pseq::Util::Options & options );

      bool net_assoc_test( Mask & , const Pseq::Util::ArgMap & );
      
      struct Aux 
      {
	Aux() { g=0; rseed=0; show_info = false; fix_null_genotypes = true; } 
	GStore * g;
	long int rseed;
	bool show_info;
	bool fix_null_genotypes;
      };
      
      struct Aux_vassoc_options {
	Aux_vassoc_options()
	{
	  show_istat = false;
	  show_meta = false;
	  separate_chr_bp = false;
	  nrep = 0;
	  yates_chisq = false;
	}
	
	bool show_istat;
	bool show_meta;
	bool separate_chr_bp;
	bool yates_chisq;
	int nrep;
      };
      

      struct Aux_glm
      {
	
	Aux_glm()
	{
	  show_meta = false;
	  nrep = 0;
	  dichot_pheno = true;
	  show_all_covar = false;
	  has_covar = false;
	  use_dosage = false;
	  use_postprobs = false;
	  softtag = "";
	}
	
	bool show_meta;
	bool dichot_pheno;
	bool show_all_covar;
	bool has_covar;
	bool use_dosage;
	bool use_postprobs;
	std::string softtag;
	int nrep;

	Data::Vector<double> y; 
	Data::Matrix<double> c; 
	
	std::vector<bool> mask; // ultimately, add mask into Data::Vector	
	std::vector<std::string> covars;
      };
      


    }
}

// Iteration functions outside of namespace 

void f_variant_association( Variant & v , void * p );

void g_set_association( VariantGroup & vars , void * p );

void g_net_assoc_collector( VariantGroup & vars , void * p );

#endif
