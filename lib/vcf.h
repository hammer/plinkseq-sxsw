#ifndef __VCF_H__
#define __VCF_H__

#include "vardb.h"

class RefDBase;
class SeqDBase;
class LocDBase;
class IndDBase;

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
     fixxy_mode = false;
     locdb = NULL;
     inddb = NULL;
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
      fixxy_mode = false;
      locdb = NULL;
      inddb = NULL;
      obs_header = false;
      set_seqdb(s);
      return_var = false;
   }
   
   bool set(File * f)
   {    
       if ( f->name() == "-" ) { from_stdin = true; return true; }  
       from_stdin = false;
       fileinfo = f;
       file.open( f->name().c_str() );
       return true;
   }

   ~VCFReader();
   
   void return_variant( const bool b ) { return_var = b; } 

   uint64_t group_id() const { return file_id; } 

   void set_seqdb( SeqDBase * );
   
   void to_refdb( RefDBase * p )
   {
     vardb = NULL;
     refdb = p;
   }
	       
 void set_region_mask( const std::set<Region> * );

 void set_fixxy( Mask * mask , LocDBase * p , IndDBase * pi );
 
 
 line_t parseLine( Variant ** pvar = NULL );
 
 void summary() const;

 bool explicit_meta;
 
 void get_meta( const std::set<std::string> & );
 
 void ignore_meta( const std::set<std::string> & );
 
 int variants_inserted() const { return vcnt; } 
 
 void insert_meta( const std::string & s ) { getMetaInformation(s); } 
 
 std::map<std::string,std::string> last_meta();

 std::vector<std::string>          last_header();


 // public, so that SampleVariant can see these

  std::vector<meta_index_t*> formats;
  int gt_field;
  
  // say we've already got the header (fix for STDIN mode, when we use 
  // two separate VCFReaders in vcfiterate.cpp
  
  void observed_header( const bool b ) { obs_header = b; } 
  void set_number_individuals( const int i ) { icnt = i; } 

 private:

  uint64_t file_id;

  bool return_var;
 
  void getMetaInformation(const std::string &);
  
  bool getHeader(const std::string &);
  
  Variant getVariant(const std::string &);
  
  VarDBase * vardb;
  
  // if loading only site-data into a REFDB instead of VARDB

  RefDBase * refdb;

  // if applying an X/Y mask/recoding
  
  bool       fixxy_mode;
  Mask     * mask;
  LocDBase * locdb;
  IndDBase * inddb;
  std::vector<sType> sex;

  // if set, then use SEQDB to check REF alleles
  
  SeqDBase * seqdb;
  
  // For the FORMAT field

  bool set_format( const std::string & );
  std::string current_format;

  
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
  bool   from_stdin;

  // Misc.

  int icnt;
  int vcnt;
  VCF_version version;
  bool obs_header;

  inline bool processVCF( const char * , int * a ) ;
  
  inline bool processVCF( const char * , double * a);
  
  inline bool processVCF( const char * , std::string * a);
  
  
};


#endif
