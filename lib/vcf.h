#ifndef __VCF_H__
#define __VCF_H__

#include "vardb.h"

class RefDBase;
class SeqDBase;

// Supported versions = 3.3, 4.0

enum VCF_version { VCF_UNKNOWN , 
		   VCF_3_3 , 
		   VCF_4_0 , 
		   VCF_4_1  };


class VCFReader {
  
 public:

  
  enum line_t { VCF_EOF = 0 ,
		VCF_INVALID = 1 ,
		VCF_META = 2 , 
		VCF_HEADER = 3 , 
		VCF_VARIANT = 4 } ;
		
    
 VCFReader(File * f,
	   const std::string & tag , 
	   VarDBase * v , 
	   SeqDBase * s = NULL )
   : vardb(v)
   {       
     set(f);       
     file_id = vardb->insert( f->name() , tag );
     version = VCF_UNKNOWN;
     pfilter = NULL;
     explicit_meta = false;
     current_format = "";
     refdb = NULL;
     vcnt = 0;
     icnt = 0;
     obs_header = false;
     set_seqdb(s);
     return_var = false;
   }

 VCFReader(File * f,
	   const std::string & tag , 
	   RefDBase * p , 
	   SeqDBase * s = NULL )
   : refdb(p)
   {       
      set(f);       
      file_id = refdb->insert( f->name() , tag );
      version = VCF_UNKNOWN;
      pfilter = NULL;
      explicit_meta = false;
      current_format = "";
      vardb = NULL;
      vcnt = 0;
      icnt = 0;
      obs_header = false;
      set_seqdb(s);
      return_var = false;
   }
   
   bool set(File * f)
   {    
     fileinfo = f;
     file.open( f->name().c_str() );
     return true;
   }
   
   void return_variant( const bool b ) { return_var = b; } 

   uint64_t group_id() const { return file_id; } 

   void set_seqdb( SeqDBase * );
   
   void to_refdb( RefDBase * p )
   {
     vardb = NULL;
     refdb = p;
   }
	       
 void set_region_mask( const std::set<Region> * );
 
 line_t parseLine( Variant ** pvar = NULL );
 
 void summary() const;

 bool explicit_meta;
 
 void get_meta( const std::set<std::string> & );
 
 void ignore_meta( const std::set<std::string> & );
 
 int variants_inserted() const { return vcnt; } 

 void insert_meta( const std::string & s ) { getMetaInformation(s); } 
 
 std::map<std::string,std::string> last_meta();

 std::vector<std::string>          last_header();
 
 private:
  
  uint64_t file_id;

  bool return_var;
 
  void getMetaInformation(const std::string &);
  
  bool getHeader(const std::string &);
  
  Variant getVariant(const std::string &);
  
  VarDBase * vardb;
  
  // if loading only site-data into a REFDB instead of VARDB

  RefDBase * refdb;

  // if set, then use SEQDB to check REF alleles
  
  SeqDBase * seqdb;
  
  // For the FORMAT field

  bool set_format( const std::string & );
  std::string current_format;
  std::vector<meta_index_t*> formats;

  // Meta-fields to ignore/read

  std::set<std::string> meta_want;
  std::set<std::string> meta_ignore;

  std::set<std::string> meta_read_var;
  std::set<std::string> meta_read_geno;
  std::set<std::string> meta_read_filter;

  // Filter-group filter

  std::set<Region> * pfilter;
  int largest_region;
  bool contains( int chr , int bp1, int bp2 );

  // File information

  File * fileinfo; 
  InFile file; 
  

  // Misc.

  int icnt;
  int vcnt;
  VCF_version version;
  bool obs_header;

  bool processVCF(std::vector<std::string>::iterator , 
		  int & a);
  
  bool processVCF(std::vector<std::string>::iterator , 
		  double & a);
  
  bool processVCF(std::vector<std::string>::iterator , 
		  std::string & a);
  
  
};


#endif
