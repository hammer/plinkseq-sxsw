#ifndef __PSB_H__
#define __PSB_H__

#include "pseq.h"

namespace Helper {
  std::string url_encode(const std::string & );
}

namespace ExomeBrowser {

  //
  // Query types
  //
  
  enum QType { Q_ERROR    = 0 ,
	       Q_VARIANT  ,
	       Q_INDIV    ,
	       Q_GENE     ,
	       Q_REGION   ,
	       Q_INDGRID  ,
	       Q_PHENOGRID,
	       Q_GENELIST ,
	       Q_METALIST ,
	       Q_PHELIST  ,
	       Q_VARSETLIST , 
	       Q_INDPHENO , 
	       Q_LOCSETLIST ,
	       Q_PROJSUMMARY ,
               Q_GRAPHICAL_VIEW };
  
  //
  // 
  //


  //
  // Core iteration functions
  //
  
  void f_display(Variant & vars, void *p);
  
  void f_display_indiv(Variant & var, void *p);
  
  void g_display_indiv(VariantGroup & vars, void *p);

  struct Aux;
  
  void show_graphical_view( GStore & g , 
			    const std::string & loc_set,
			    const QType & q, 
			    Aux & a, 
			    const std::string & pheno,
			    const std::string & pwd , 
			    const std::string & project_path , 
			    Mask & );

  // correspnds to "pgrid"
  void show_indiv( Aux & a );

  void show_varsets( Aux & a );

  // corresponds to indgrid;
  void index_grid( Aux & a );

  //
  // Helper functions
  //
  
  // Write HTML head

  void write_html_header( const std::string & head_preamble ="" );

  // Write top panel for browser page

  void write_start_page( const GStore & g , 
			 const std::string & loc_set,
			 const QType & q, 
			 Aux & a, 
			 const std::string & pheno,
			 const std::string & pwd , 
			 const std::string & project_path );


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
        std::string regs="",
	std::string varset="",
	std::string ref_append="",
	std::string loc_append="",
	std::string indiv_list="")
    {
      fields["proj"] = project;
      fields["q"] = q;
      fields["gene"] = gene;
      fields["masks"] = masks;
      fields["meta"] = meta;
      fields["pheno"] = pheno;
      fields["regs"] = regs;
      fields["ref_append"] = ref_append;
      fields["loc_append"] = loc_append;
      fields["varset"] = varset;
      fields["indiv_list"] = indiv_list;
    }


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
      std::string s = "pbrowse.cgi?";
      for (std::map<std::string, std::string>::iterator i = fields.begin(); i != fields.end(); i++)
        {
          if (i->second != "")
            s += i->first + "=" + Helper::url_encode(i->second) + "&";
        }
      return s;
    }

    std::string printLink(const std::string & text = "" )
    {
      return text == "" ? printURL() : "<a href=\"" + printURL() + "\">" + text + "</a>";
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
      single_transcript = false;
      region_search = false;
      reg_list = reg_list_url = "";
      indiv_list_url = "";
      url = NULL;
      indiv_genogrid = false;
    }
    
    ~Aux() {
      if ( url ) 
	{
	  delete url;
	  url = NULL;
	}
    }

    GStore * g;


    // Phenotype information

    std::string phenotype_name;
    bool show_phenotype;
    
    bool indiv_genogrid;
    
    // Gene info

    Region region;
    bool single_transcript;
    std::string loc_set;
    std::string genename;
    std::string reg_list;
    std::string reg_list_url;
    
    // Individual includes
    std::set<std::string> indiv_list;
    std::vector<std::string> indiv_list_vec;
    std::string indiv_list_url;
    
    // Optional variant meta-fields
    
    std::vector<std::string> mf;
    std::map<std::string,bool> mfpp; // skip pretty-print
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


    //
    // Appends
    //

    std::vector<std::string> ref_append;
    std::vector<std::string> loc_append;
    
    std::string ref_append_url;
    std::string loc_append_url;

    // Variant sets
    
    std::string varset_url;
    std::vector<std::string> varset;
    
    bool add_annot;

    // Regions, and other genes to add
    
    std::vector<Region> regions;
    std::vector<std::string> genes;
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
      if ( url == NULL )
        {
          url = new BrowserURL(
			       has_form_value("proj") ? form["proj"] : "", // project
			       "", // q
			       "", // gene (?redundant?)
			       msk_print(), // mask 
			       mf_print(),  // meta-inf
			       phenotype_name, // pheno
			       reg_list_url,   // main region list
			       varset_url,     // variant sets
			       ref_append_url, // ref-appends
			       loc_append_url, // loc-appends
			       indiv_list_url  // indiv mask
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
  
  std::string html_encode( const std::string & data) 
    {
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
  
  
  std::string url_encode( const std::string & data) 
    {
      
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
	      char s[3]; // SMP. fix, 2 + 1 for \0
	      sprintf(s, "%x", data[pos]);
	      buffer.append( s, 2 );//converts char 255 to string "ff"
	    }
	}
      
      return buffer;
    }
  
};

#endif
