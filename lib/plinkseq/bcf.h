#ifndef __PSEQ_BCF_H__
#define __PSEQ_BCF_H__

#include "bgzf.h"

#include <vector>
#include <string>
#include <map>

class Variant;
class SampleVariant;
class VarDBase;
class IndividualMap;
class Mask;

struct BCF_header {

  std::vector<std::string> seq_names;
  std::map<std::string,int> seq_map; // for reverse lookup
  std::vector<std::string> sample_names;
  std::vector<std::string> meta_text;
  
  // two string dictionaries : one for contigs, one and INFO/FILTER/FORMAT

  std::map<std::string,int> string2num;
  std::vector<std::string> num2string;

  int dictionary() const 
  {
    return num2string.size();
  }

  int dictionary( const std::string & s ) 
  {
    std::map<std::string,int>::iterator ii = string2num.find( s );
    if ( ii != string2num.end() ) return ii->second;
    const int p = num2string.size();
    string2num[ s ] = p;
    num2string.push_back( s );
    return p;
  }

  bool dictionary( const int p , std::string * str ) const
  {
    if ( p < 0 || p >= num2string.size() ) return false;
    *str = num2string[p];
    return true;
  }

  bool dictionary( const std::string & str , int * c ) const
  {
    std::map<std::string,int>::const_iterator ii = string2num.find( str );
    if ( ii != string2num.end() ) 
      {
	*c = ii->second;
	return true;
      }
    return false;
  }

  // contig dictionary
  std::map<std::string,int> contig_string2num;
  std::vector<std::string> contig_num2string;

  int contig_dictionary() const 
  {
    return contig_num2string.size();
  }

  int contig_dictionary( const std::string & s ) 
  {
    std::map<std::string,int>::iterator ii = contig_string2num.find( s );
    if ( ii != contig_string2num.end() ) return ii->second;
    const int p = contig_num2string.size();
    contig_string2num[ s ] = p;
    contig_num2string.push_back( s );
    return p;
  }

  bool contig_dictionary( const int p , std::string * str ) const
  {
    if ( p < 0 || p >= contig_num2string.size() ) return false;
    *str = contig_num2string[p];
    return true;
  }

  bool contig_dictionary( const std::string & str , int * c ) const
  {
    std::map<std::string,int>::const_iterator ii = contig_string2num.find( str );
    if ( ii != contig_string2num.end() ) 
      {
	*c = ii->second;
	return true;
      }
    return false;
  }

  void clear() 
  { 
    seq_names.clear(); 
    sample_names.clear(); 
    meta_text.clear(); 
    
    num2string.clear();
    string2num.clear();

    contig_num2string.clear();
    contig_string2num.clear();
  } 

};


class BCF {
  
 public:
 
  enum bcf2_typed_value { 
    bcf2_flag    = 0 , 
    bcf2_int8_t  = 1 ,  
    bcf2_int16_t = 2 ,  
    bcf2_int32_t = 3 , 
    bcf2_float   = 5 , 
    bcf2_char    = 7 , 
    bcf2_void    = 99 
  }; 
  
  void dictionary_set( const int k , const std::string & v )
  {
    if ( k == hdr.num2string.size() ) 
      {
	hdr.num2string.push_back( v );
	hdr.string2num[v] = k;
	return;
      }
    else if ( k > hdr.num2string.size() )
      {
	hdr.num2string.resize( k + 1 );
      }    
    hdr.num2string[k] = v;
    hdr.string2num[v] = k;
  }
  
  void contig_dictionary_set( const int k , const std::string & v )
  {
    if ( k == hdr.contig_num2string.size() ) 
      {
	hdr.contig_num2string.push_back( v );
	hdr.contig_string2num[v] = k;
	return;
      }
    else if ( k > hdr.contig_num2string.size() )
      {
	hdr.contig_num2string.resize( k + 1 );
      }    
    hdr.contig_num2string[k] = v;
    hdr.contig_string2num[v] = k;
  }
    
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

  void read_header( VarDBase * vardb = NULL );

  bool index_record( );

  bool read_record( Variant & , SampleVariant & , SampleVariant & , int64_t offset );
  bool read_record( Variant & , SampleVariant & , SampleVariant & );
  bool read_genotypes( int64_t offset , IndividualMap * , SampleVariant * , SampleVariant * , Mask * );
  
  bool write_header();
  bool write_record( const Variant & );
  

  // Read BCF2 typed value

  bool      get_next_typed_value( bcf2_typed_value * type , uint32_t * length );
  uint32_t  write_typing_byte( bcf2_typed_value type , uint32_t length );
  uint32_t  size_typed_value( bcf2_typed_value type , uint32_t length );
  uint32_t  size_typing_byte( bcf2_typed_value type , uint32_t length );
  
 private:
  
  BGZF * file;

  bool have_read_header;

  BCF_header hdr;

  int n;
  
  int last_record_n_fmt;
  int last_record_n_sample;

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

  inline bool read( int8_t & i ); 
  inline bool read( uint8_t & i ); 
  inline bool read( uint8_t * i , int32_t l); 

  inline bool read( int16_t & i );
  inline bool read( uint16_t & i );

  inline bool read( int32_t & i );  
  inline bool read( uint32_t & i );  

  inline bool read( std::vector<int8_t> & i );  
  inline bool read( std::vector<int16_t> & i );  
  inline bool read( std::vector<int32_t> & i );  
  inline bool read( std::vector<float> & i );  


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

  inline void write( int8_t );
  inline void write( int16_t );
  inline void write( int32_t );

  inline void write( float );
  inline void write( double );
  inline void write( const std::string & s );
  void write( const std::vector<std::string> & s );
  
  //
  // BCF2 typed values
  //

  inline bool read_typed_string( std::string & s );
  inline bool read_typed_int( int & i , bool & missing );
  inline bool read_FILTER( std::vector<int> & s );

  // typed-writers return number of bytes used

  inline uint32_t  write_typed_int( int i , bool missing = false );
  inline uint32_t  write_typed_float( float f , bool missing = false );
  inline uint32_t  write_typed_string( const std::string & s , bool missing = false );

  inline uint32_t  write_typed_vec( const std::vector<int8_t>       & , const std::vector<bool> * missing = NULL);
  inline uint32_t  write_typed_vec( const std::vector<int16_t>      & , const std::vector<bool> * missing = NULL);
  inline uint32_t  write_typed_vec( const std::vector<int32_t>      & , const std::vector<bool> * missing = NULL);
  inline uint32_t  write_typed_vec( const std::vector<float>        & , const std::vector<bool> * missing = NULL);

  // note-- current confusion over string as vector of char, vs. how to define a vector of strings
  inline uint32_t write_typed_vec( const std::vector<std::string>  & , const std::vector<bool> * missing = NULL);


  inline void write( const std::vector<int8_t> & );
  inline void write( const std::vector<int16_t> & );
  inline void write( const std::vector<int32_t> & );
  inline void write( const std::vector<double> & );

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
  
  inline uint32_t swap_endian( uint32_t i ) 
  { 
    unsigned char b1, b2, b3, b4;    
    b1 = i & 255;
    b2 = ( i >> 8 ) & 255;
    b3 = ( i >> 16 ) & 255;
    b4 = ( i >> 24 ) & 255;    
    return ((int)b1 << 24) + ((int)b2 << 16) + ((int)b3 << 8) + b4;
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
