#ifndef __PSB_H__
#define __PSB_H__

#include "pseq.h"

namespace ExomeBrowser {
  
  //
  // Query types
  //
  
  enum QType { Q_ERROR    = 0 ,
	       Q_VARIANT  = 1 , 
	       Q_INDIV    = 2 , 
	       Q_GENE     = 3 , 
	       Q_REGION   = 4 ,
	       Q_GENELIST = 5 ,
	       Q_METALIST = 6 ,
	       Q_PHELIST  = 7 , 
	       Q_LOCSETLIST = 8 , 
	       Q_PROJSUMMARY = 9 };
  

  //
  // Core iteration functions
  //
  
  void f_display(Variant & vars, void *p);
  
  void f_display_indiv(Variant & var, void *p);
  
  void g_display_indiv(VariantGroup & vars, void *p);

  
  //
  // Helper functions
  //
  
  // Indicate which exon a variant is in
  
  int exon_overlap( const Region & reg , int pos );

  // Write a link to dbSNP, is name starts with rs#####

  std::string rs_link(const std::string & label );

  // Pretty print long strings

  std::string pp( const std::string & str , const int len = 15 );
  
  struct Aux {
    Aux() 
    {
      add_annot = false;
      g = NULL;
      show_phenotype = false;
      phenotype_name = "";
      multi_transcripts = false;
      reg_list = reg_list_url = "";
    }

    GStore * g;


    // Phenotype information

    std::string phenotype_name;
    bool show_phenotype;
    
    
    // Gene info

    Region region;
    bool multi_transcripts;
    std::string loc_set;
    std::string genename;
    std::string reg_list;
    std::string reg_list_url;
    
    // Optional variant meta-fields
    
    std::vector<std::string> mf;
    std::string mf_print() 
    { 
      std::string s; 
      for (int i=0;i<mf.size(); i++) 
	{ 
	  if (i>0) s += ","; 
	  s += mf[i]; 
	} 
      return s; 
    }

    
    // Mask

    std::vector<std::string> msk;
    std::string msk_print() 
    { 
      std::string s; 
      for (int i=0;i<msk.size(); i++) 
	s += ( i ? "," : "" ) + msk[i]; 
      return s; 
    }
    
    bool add_annot;

    // Regions, and other genes to add

    std::vector<Region> regions;
    std::vector<std::string> other_genes;
    bool extended_search;


    // Form information

    std::map<std::string,std::string> form;

    bool has_form_value(const std::string & s) const 
    { return form.find(s) != form.end(); }

    std::string print_form_value(const std::string & s) 
    { return has_form_value(s) ? s+"="+form[s] : ""; }

    void add_form_value(const std::string & k, const std::string & v) 
    { form[k] = v; }
    

    // Individual mode

    bool indiv_mode;
    std::string indiv_id;

    
    // Structure to build table

    std::map<int,std::string> table_row;
    std::string headers;
    int vcnt;
  };



  //
  // Helper list functions
  //
  
  void make_gene_list(Aux * a);
  
  void make_phe_list(Aux * a);
  
  void make_mf_list(Aux * a);

  void make_locset_list(Aux * a);

  void make_proj_summary(Aux * a);
  
};

#endif
