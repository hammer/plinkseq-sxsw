#ifndef __PSEQ_DOSE_H__
#define __PSEQ_DOSE_H__

#include <string>
#include <vector>

#include "plinkseq/helper.h"

class ifstream;
class SeqDBase;
class IndDBase;
class VarDBase;
class Variant;

class DoseReader {
  
 public:
  
  static double Rsq( const Variant & v , double * maf = NULL );

  DoseReader( VarDBase * v )
    {
      vardb = v;     
      seqdb = NULL;
      filetag = "1";
      tagname = "";
      separate_map = false;
      separate_fam = false;
      dose_filename = "";
      map_filename = "";
      fam_filename = "";
      has_meta = false;
      map_pos_only = false;
      map_allele_only = false;
      make_hard_call = true;
      hard_call_prob_threshold = 0.9;
      hard_call_dosage_threshold = 0.1;
      skip_header = false;
      spaced = false;
      vardelim = false;
      has_rsid = true;
      maf = 0.0;
      mac = 0;
      maf_threshold = false;
      mac_threshold = false;
    }  

  void set_mac( const int m )
  {
    mac_threshold = true;
    mac = m;
  }

  void set_maf( const double m )
  {
    maf_threshold = true;
    maf = m;
  }

  void set_skip_header( const bool b ) 
  {
    skip_header = b;
  }

  void set_spaced() 
  {
    spaced = true;
  }

  void set_variable_delimiter()
  {
    vardelim = true;
  }

  void set_id( const bool b )
  {
    has_rsid = b;
  }

  void set_mapfile( const std::string & n )
  {
    map_filename = n;
    separate_map = true;    
  }
  
  void set_mapfile_only_pos() 
  {
    map_pos_only = true;
  }
  
  void set_mapfile_only_alleles()
  {
    map_allele_only = true;
  }

  void set_famfile( const std::string & n )
  {
    fam_filename = n;
    separate_fam = true;
  }
    
  void use_seqdb( SeqDBase * s ) 
  {
    seqdb = s;
  }
  
  void set_filetag( const std::string & t ) 
  {
    filetag = t;
  }


  void set_metatag( const std::string & t , const std::string & ttype ) 
  {
    tagname = t;
    storage_format = ttype;
  }
  
  void set_meta( const std::string & t )
  {
    Helper::checkFileExists( t );
    meta_filename = t;
    has_meta = true;
  }

  void set_input_format( const std::string & t)
  {
    format_prob2 = t == "prob2" ;
    format_prob3 = t == "prob3" ;
    format_dose2 = t == "dose2" ;
    format_dose1 = t == "dose1" ;
  } 

  bool read_dose( const std::string & n );
  
  void set_hard_call_threshold( const double d ) 
  {
    // if being set by user, they know which storage mode is being
    // used, thus set both to user-specified value
    hard_call_prob_threshold = d;  
    hard_call_dosage_threshold = d;
  }
  
  
 private:
  
  int read_fam( const int );
  int read_bim();
  void read_meta( const int );

  std::string fam_filename;
  std::string map_filename;
  std::string dose_filename;

  std::string meta_filename;
  std::vector<std::string> metas;
  
  // formats
  
  std::string input_format;
  std::string storage_format;
  std::string tagname;  
  std::string filetag;
  
  // extras
  
  bool separate_map;
  bool separate_fam;  
  bool map_pos_only;
  bool map_allele_only;

  bool format_prob2;
  bool format_prob3;
  bool format_dose1;
  bool format_dose2;

  bool skip_header;
  bool spaced;
  bool vardelim;
  bool has_rsid;

  bool has_meta;
  bool make_hard_call;
  double hard_call_prob_threshold;
  double hard_call_dosage_threshold;
  
  bool maf_threshold;
  bool mac_threshold;
  double maf;
  int mac;


  // DB pointers
  
  VarDBase * vardb;
  SeqDBase * seqdb;
  
};


#endif
