#ifndef __PSEQ_BCF_H__
#define __PSEQ_BCF_H__

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

class BCF {
  
 public:
  
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
 
  // convert to VCF functions
  std::string vcf_header();

  // VCF --> BCF wrapper
  bool vcf2bcf( const std::string & vcfname , const std::string & bcfname );
  
 private:
  
  BGZF * file;
  
  BCF_header hdr;

  int n;

  std::string filename; 
  
  int readmode;
  
  VarDBase * vardb;
  
  int file_id;
  
  //
  // I/O functions
  //


  // reading

  inline bool read( char & c );
  inline bool read( std::vector<char> & buf, int l);
  inline bool read( uint8_t & i ); 
  inline bool read( uint8_t * i , int32_t l); 
  inline bool read( uint16_t & i );
  inline bool read( int32_t & i );  
  inline bool read( float & f );  
  inline bool read( double & i );  
  inline bool read( std::string & s, int l);  
  inline bool read( std::string & s);
  inline bool read( std::vector<std::string> & s );

  // writing

  inline void write( char c );
  inline void write( const std::vector<char> & buf, int l);
  inline void write( uint8_t );
  inline void write( uint16_t );
  inline void write( uint32_t );
  inline void write( uint64_t );
  inline void write( int32_t );
  inline void write( float );
  inline void write( double );
  inline void write( const std::string & s );
  void write( const std::vector<std::string> & s );

  // random access

  inline void seek(int64_t offset);
  inline int64_t tell();

  // handling endianness
  
  enum endian_t {     
    MACHINE_LITTLE_ENDIAN = 0 ,
    MACHINE_BIG_ENDIAN = 1 };
  
  endian_t endian;
  
  endian_t determine_endian() 
  {
    int i = 1;
    char *p = (char *)&i;
    return p[0] == 1 ? MACHINE_LITTLE_ENDIAN : MACHINE_BIG_ENDIAN ;
  }

  inline uint16_t swap_endian( uint16_t v )
    {
      return (uint16_t)(((v & 0x00FF00FFU) << 8) | ((v & 0xFF00FF00U) >> 8));
    }
  
  inline int32_t swap_endian( int32_t i ) 
  { 
    unsigned char b1, b2, b3, b4;    
    b1 = i & 255;
    b2 = ( i >> 8 ) & 255;
    b3 = ( i >> 16 ) & 255;
    b4 = ( i >> 24 ) & 255;    
    return ((int)b1 << 24) + ((int)b2 << 16) + ((int)b3 << 8) + b4;
  }

  uint64_t swap(double d)
  {
    uint64_t a;
    unsigned char *dst = (unsigned char *)&a;
    unsigned char *src = (unsigned char *)&d;

    dst[0] = src[7];
    dst[1] = src[6];
    dst[2] = src[5];
    dst[3] = src[4];
    dst[4] = src[3];
    dst[5] = src[2];
    dst[6] = src[1];
    dst[7] = src[0];

    return a;
  }

  // unswap using char pointers
  double swap_double(uint64_t a) 
  {

    double d;
    unsigned char *src = (unsigned char *)&a;
    unsigned char *dst = (unsigned char *)&d;

    dst[0] = src[7];
    dst[1] = src[6];
    dst[2] = src[5];
    dst[3] = src[4];
    dst[4] = src[3];
    dst[5] = src[2];
    dst[6] = src[1];
    dst[7] = src[0];

    return d;
  }
  
  uint32_t swap(float d)
  {
    uint32_t a;
    unsigned char *dst = (unsigned char *)&a;
    unsigned char *src = (unsigned char *)&d;
    dst[0] = src[3];
    dst[1] = src[2];
    dst[2] = src[1];
    dst[3] = src[0];
    return a;
  }
  
  // unswap using char pointers
  float swap_float(uint32_t a) 
  {
    float d;
    unsigned char *src = (unsigned char *)&a;
    unsigned char *dst = (unsigned char *)&d;
    dst[0] = src[3];
    dst[1] = src[2];
    dst[2] = src[1];
    dst[3] = src[0];
    return d;
  }

};

#endif
