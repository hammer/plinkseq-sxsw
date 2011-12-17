#ifndef __PSB_H__
#define __PSB_H__

#include "pseq.h"

namespace Helper {
  std::string url_encode(std::string);
}

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
  
  struct BrowserURL {

    BrowserURL(
        std::string project="",
        std::string q="",
        std::string gene="",
        std::string masks="",
        std::string meta="",
        std::string pheno="",
        std::string regs=""
    )
    {
      fields["proj"] = project;
      fields["q"] = q;
      fields["gene"] = gene;
      fields["masks"] = masks;
      fields["meta"] = meta;
      fields["pheno"] = pheno;
      fields["regs"] = regs;

    }

    // todo: why doesn't this work?
//    BrowserURL(Aux * a)
//    {
//      fields["project"] = a->print_form_value("proj");
//      fields["meta"] = a->mf_print();
//      fields["masks"] = a->msk_print();
//      fields["pheno"] = a->phenotype_name;
//      fields["project"] = a->reg_list_url;
//    }

    std::map<std::string, std::string> fields;

    BrowserURL * addField(std::string key, std::string val) {
      fields[key] = val;
      return this;
    }

    BrowserURL * removeField(std::string key) {
      fields[key] = "";
      return this;
    }

    std::string printURL()
    {
      std::string s = "/pbrowse.cgi?";
      for (std::map<std::string, std::string>::iterator i = fields.begin(); i != fields.end(); i++)
        {
          if (i->second != "")
            s += i->first + "=" + Helper::url_encode(i->second) + "&";
        }
      return s;
    }

    std::string printLink(std::string text)
    {
      return "<a href=\"" + printURL() + "\">" + text + "</a>";
    }
  
  };
  

  // ID for transcript --> gene symbol (LOCDB alias)

  std::string symbol = "symbol";

  
  struct Aux {
    Aux() 
    {
      add_annot = false;
      g = NULL;
      show_phenotype = false;
      phenotype_name = "";
      multi_transcripts = false;
      region_search = false;
      reg_list = reg_list_url = "";
      url = NULL;
    }

    ~Aux() {
      delete url;
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

    
    // Masks

    std::vector<std::string> msk;
    std::string msk_print() 
    { 
      std::string s; 
      for (int i=0;i<msk.size(); i++) 
	s += ( i ? "," : "" ) + msk[i]; 
      return s; 
    }
    
    // Include filters

    std::string inc_fltr;
    std::string vinc_fltr;



    bool add_annot;

    // Regions, and other genes to add
    
    std::vector<Region> regions;
    std::vector<std::string> genes;
    bool single_transcript;
    bool region_search;

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

    BrowserURL * url;
    BrowserURL * getURL()
    {
      if ( url == NULL)
        {
          url = new BrowserURL(
              has_form_value("proj") ? form["proj"] : "",
              "",
              "",
              msk_print(),
              mf_print(),
              phenotype_name,
              reg_list_url
              );
        }
      return url;
    }
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

namespace Helper {
	
	std::string html_encode(std::string & data) {

	  //	  return data;

	  std::string buffer;
	  buffer.reserve(data.size());
	  for(size_t pos = 0; pos != data.size(); ++pos) {
	    switch(data[pos]) {
	    case '&':  buffer.append("&amp;");       break;
	    case '\"': buffer.append("&quot;");      break;
	    case '\'': buffer.append("&apos;");      break;
	    case '<':  buffer.append("&lt;");        break;
	    case '>':  buffer.append("&gt;");        break;
	    default:   buffer.append(1, data[pos]); break;
	    }
	  }
	  return buffer;
	}
	
	std::string url_encode(std::string data) {
	  
	  //return data;
	  
	  std::string buffer="";
	  buffer.reserve(data.size());
	  for(size_t pos = 0; pos != data.size(); ++pos)
	    {
	      if ( (48 <= data[pos] && data[pos] <= 57) ||//0-9
		   (65 <= data[pos] && data[pos] <= 90) ||//abc...xyz
		   (97 <= data[pos] && data[pos] <= 122) || //ABC...XYZ
		   (data[pos]=='~' || data[pos]=='!' || data[pos]=='*' || 
		    data[pos]=='(' || data[pos]==')' || data[pos]=='\'')
		   )
	        {
		  buffer.append( &data[pos], 1);
	        }
	        else
		  {
	            buffer.append("%");
		    char s[2];
		    sprintf(s, "%x", data[pos]);
	            buffer.append( s, 2 );//converts char 255 to string "ff"
		  }
	    }
	  return buffer;
	}
	
};

#endif
