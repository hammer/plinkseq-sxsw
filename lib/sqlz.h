
#include <cstdlib>
#include <zlib.h>
#include <cassert>

static void compressFunc( sqlite3_context *context, int argc, sqlite3_value **argv)
{
  int nIn, nOut;
  long int nOut2;
  const unsigned char *inBuf;
  unsigned char *outBuf;
  assert( argc==1 );
  nIn = sqlite3_value_bytes(argv[0]);
  inBuf = (const unsigned char*)sqlite3_value_blob(argv[0]);
  nOut = 13 + nIn + (nIn+999)/1000;
  outBuf = (unsigned char *)malloc( nOut+4 );
  outBuf[0] = nIn>>24 & 0xff;
  outBuf[1] = nIn>>16 & 0xff;
  outBuf[2] = nIn>>8 & 0xff;
  outBuf[3] = nIn & 0xff;
  nOut2 = (long int)nOut;
  compress(&outBuf[4], (uLongf*)&nOut2, inBuf, nIn);
  sqlite3_result_blob(context, outBuf, nOut2+4, free);
}

static void uncompressFunc( sqlite3_context *context, int argc, sqlite3_value **argv)
{
  unsigned int nIn, nOut, rc;
  const unsigned char *inBuf;
  unsigned char *outBuf;
  long int nOut2;
  
  assert( argc==1 );
  nIn = sqlite3_value_bytes(argv[0]);
  if( nIn<=4 ){
    return;
  }
  inBuf = (const unsigned char *)sqlite3_value_blob(argv[0]);
  nOut = (inBuf[0]<<24) + (inBuf[1]<<16) + (inBuf[2]<<8) + inBuf[3];
  outBuf = (unsigned char *)malloc( nOut );
  nOut2 = (long int)nOut;
  rc = uncompress(outBuf, (uLongf*)&nOut2, &inBuf[4], nIn);
  if( rc!=Z_OK ){
    free(outBuf);
  }else{
    sqlite3_result_blob(context, outBuf, nOut2, free);
  }
}

