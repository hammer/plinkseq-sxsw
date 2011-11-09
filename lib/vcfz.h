#ifndef __PSEQ_VCFZ_H__
#define __PSEQ_VCFZ_H__

#include "bgzf.h"

#include <vector>
#include <string>
#include <map>

class Variant;
class SampleVariant;
class VarDBase;

struct BCF_header {
  std::vector<std::string> seq_names;
  std::map<std::string,int> seq_map; // for reverse lookup
  std::vector<std::string> sample_names;
  std::vector<std::string> meta_text;
  void clear() { seq_names.clear(); sample_names.clear(); meta_text.clear(); } 
};


class VCFZ {
  
 public:
  
    enum bcf_type_t {
      BCF_undef = 0 , 
    BCF_genotype ,
    BCF_uint8 , 
    BCF_uint16 , 
    BCF_int32 , 
    BCF_uint32 , 
    BCF_uint64 ,
    BCF_float , 
    BCF_double , 
    BCF_flag , 
    BCF_char ,
    BCF_string };   
  
  struct bcf_meta_t {
    
    bcf_meta_t( bcf_type_t type , int len ) : type(type) , len(len) { } 

    bcf_meta_t( ) 
    {
      type = BCF_undef;
      len = 0;
    }

    bcf_type_t type;

    int        len;   // +ve number = length; 0 = variable-number; -1 = #alt-alleles, -2 #alleles; -3 = #genotypes
  };


  std::map<std::string,bcf_meta_t> bcftype;

  
  // rmode 0=read 1=write
  BCF( const std::string & f , int m = 1 )
    {
      endian = determine_endian();      
      file = NULL;   // BGZF handle
      file_id = 0;   // refers to VADDB 
      filename = f;  
      readmode = m;
      n = 0;
    }

  void set_types();
  bool open();
  void close();
  void reading() { readmode = 1; }
  void writing() { readmode = 0; }
  int sample_n() { return n; }

  // Main client functions
  void set_nind( const int i ) { n = i; }
  void read_header( VarDBase * vardb = NULL );
  bool index_record( );
  bool read_record( Variant & , SampleVariant & , SampleVariant & );
  bool read_record( Variant & , SampleVariant & , SampleVariant & , int64_t offset );
  bool write_header();
  bool write_record( const Variant & );

  template<class T> 
    void set_size( const int nelem , std::vector<T> & tmp , const int nallele , const int ngen )
    {
      if ( nelem == 0 ) write( (int32_t)tmp.size() );
      else if ( nelem > 0 && tmp.size() != nelem ) tmp.resize( nelem ); 
      else if ( nelem == -1 && tmp.size() != nallele - 1 ) tmp.resize( nallele - 1 );
      else if ( nelem == -2 && tmp.size() != nallele ) tmp.resize( nallele );
      else if ( nelem == -3 && tmp.size() != ngen ) tmp.resize( ngen );
    }
  
  // convert to VCF functions
  std::string vcf_header();
  
  // VCF --> BCF wrapper
  bool create_header();
  
 private:
  
  BGZF * file;
    
  int readmode;
  
  VarDBase * vardb;
  
  int file_id;
  
  //
  // I/O functions
  //

  // read a single line


  // writing
  // ( not yet )
  
  // random access

  inline void seek(int64_t offset);
  inline int64_t tell();

};

#endif
