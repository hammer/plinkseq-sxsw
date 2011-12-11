#ifndef __PSEQ_VCFZ_H__
#define __PSEQ_VCFZ_H__

#include "bgzf.h"
#include "helper.h"
#include "meta.h"

#include <vector>
#include <string>
#include <map>

class Variant;
class SampleVariant;
class VarDBase;
class Mask;

class VCFZ {
  
 public:
  
  // rmode 0=read 1=write
  VCFZ( const std::string & f , VarDBase * v = NULL , int m = 0 )
    {      
      file = NULL;   // BGZF handle
      file_id = 0;   // refers to VADDB 
      vardb = v;
      filename = f;
      readmode = m;
    }

  ~VCFZ()
      {
	  if ( file ) 
	  {
	      bgzf_close( file );
	      file = NULL;
	  }
      }

  void set_vardb( VarDBase * v ) { vardb = v; } 
  bool open();
  void close();
  void reading() { readmode = 1; }
  void writing() { readmode = 0; }
  

  // Main client functions

  void read_header( Mask & mask );
  bool index_record( );
  bool read_record( Variant & , SampleVariant & , SampleVariant & , SampleVariant & );
  bool read_record( Variant & , SampleVariant & , SampleVariant & , SampleVariant & , int64_t offset );

  bool write_header();
  bool write_record( const Variant & );
  
  // convert to VCF functions
  std::string vcf_header();

  std::string name() const { return filename; } 

 private:
  
  BGZF *      file;    
  int         readmode;  
  VarDBase *  vardb;  
  int         file_id;
  std::string filename;

  // parsing the FORMAT tag
  bool   set_format( const std::string & f );
  static std::string                current_format;
  static int                        gt_field;
  static std::vector<meta_index_t*> formats;


  // I/O functions
  
  // read a single line into a char_tok object
  bool read_line( std::vector<char> * line );
    
  // writing
  // ( not yet )
  
};

#endif


