#ifndef __PSEQ_BCF_H__
#define __PSEQ_BCF_H__

#include "zfstream.h"
#include "filemap.h"
#include <string>
#include <vector>

class BCF {
  
 public:
 
  BCF( const std::string & f ) { } 
  
 private:
  
  BGZF dat;  
  
};



class BCFReader {

 public:

  BCFReader( const std::string & filename ) 
    {
      endian = determine_endian();
      open( filename );
    }

  bool open( const std::string & filename );
  
  bool parse();
    
  
 private:
  
  InFile BCF;
  
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


  inline float swap_endian( float f )
  {
    union
    {
      float f;
      unsigned char b[4];
    } dat1, dat2;
    
    dat1.f = f;
    dat2.b[0] = dat1.b[3];
    dat2.b[1] = dat1.b[2];
    dat2.b[2] = dat1.b[1];
    dat2.b[3] = dat1.b[0];
    return dat2.f;
  }
  
  // read functions
  bool read( std::string & , int32_t );  
  bool read( std::string & );  
  bool read( std::vector<std::string> & s );
  bool read( int32_t & );
  bool read( uint16_t & );
  bool read( uint8_t & );  
  bool read( double & );

};

#endif
