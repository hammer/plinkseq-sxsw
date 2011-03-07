#include "bcf.h"
#include "helper.h"
#include "variant.h"

bool BCFReader::open( const std::string & filename ) 
{

  Helper::checkFileExists( filename );
  
  BCF.open( filename.c_str() , std::ios::in | std::ios::binary );

  // BCF format

  // DP uint16
  // GL float
  
  // Read BCF header
  
  char ch[5];
  BCF.read(ch,4);
  ch[4] = '\0';
  
  if ( std::string(ch) != "BCF\4" ) 
    Helper::halt( "problem with format of BCF file (1) " ); 
  
  // Sequence/chromosome names
  
  std::vector<std::string> seq_names;
  if ( ! read( seq_names ) ) Helper::halt( "problem with format of BCF file" );
  
  std::vector<std::string> sample_names;
  if ( ! read( sample_names ) ) Helper::halt( "problem with format of BCF file" );

  int n = sample_names.size() ;
  
  std::string meta_text;
  if ( ! read( meta_text ) ) Helper::halt( "problem with format of BCF file" );
  
  //
  // Process BCF headers
  //


  //
  // Read each record from BCF
  // 

  while ( ! BCF.eof() )
    {
      
      int32_t seq_id; // 'chromosome' code
      int32_t bp;     // position (BP1)
      double qual;     // quality score (actually reads a float and converts)
      std::vector<std::string> mstr; // variant meta-information
      // ID REF ALT FILTER INFO FORMAT (null padded)

      read( seq_id ) ;
      read( bp );
      read( qual );
      read( mstr );
      
//       std::cout << "det = " 
// 		<< seq_names[ seq_id ] << "\t" 
// 		<< bp << "\t"
// 		<< qual << "\t"
// 		<< mstr.size() << "\n";
            
      //
      // Read genotype; # alleles based on ALT tag, comma-sep values
      //

      std::vector<std::string> alt = Helper::char_split( mstr[2] , ',' );      
      int nalt = alt.size() + 1;
      int ngen = (int) (nalt * (nalt+1) * 0.5);
      

      //
      // Populate Variant
      //

      Variant var;

      var.chromosome( Helper::chrCode( seq_names[ seq_id ] ) );
      var.position( bp );
      var.stop( bp + mstr[1].size() - 1 );
      var.name( mstr[0] == "" ? "." : mstr[0] );

      var.consensus.quality( qual );
      var.consensus.reference( mstr[1] );
      var.consensus.alternate( mstr[2] );
      var.consensus.filter( mstr[3] );
      var.consensus.info( mstr[4] );
      
      var.consensus.calls.size( n );

      // 
      // Genotype fields to expect based on 
      //
      
      // use set_format() function in vcf.cpp

      std::vector<std::string> format = Helper::char_split( mstr[5] , ':' );

      for ( int t = 0 ; t < format.size(); t++ )
	{
	  
	  // BCF recognised tags are:
	  //     DP  uint16_t
	  //     GL  float   (*ngen)
	  //     GT  uint8_t
	  //     GQ  uint8_t 
	  //     HQ  uint8_t (*2)
	  //     PL  uint8_t (*ngen)
	  //     SP  uint8_t  // strand-bias p-value 
	  
	  if ( format[t] == "DP" ) 
	    {

	      MetaInformation<GenMeta>::field( "DP", META_INT , 1, "Read depth" );	      
	      for ( int per = 0 ; per < n ; per++ )
		{
		  uint16_t dp; 
		  read(dp);
		  var(per).meta.set( "DP" , (int)dp );		  
		}
	      
	    }
	  else if ( format[t] == "GT" )
	    {
	      for ( int per = 0 ; per < n ; per++ )
		{
		  uint8_t gt;
		  read( gt );
		  
		  // genotype encoded phase << 6 | allele1 << 3 | allele2 

		  Genotype g( &var );
		  
		  g.null(   ( gt >> 7 ) & 1 );
		  g.phased( ( gt >> 6 ) & 1 );
		  
		  // TMP -- ignore >2 alleles, but casting to bools
		  
		  g.pat( ( gt >> 3 ) & 7 );
		  g.mat( gt & 7 );
		  
		  var(per) = g;

// 		  if ( missing ) std::cout << "gt = " << "./." << "\n";
// 		  else 
// 		    std::cout << "gt = " << allele1 << ( phase ? "|" : "/" ) << allele2 << "\n";

		}

	      
	    }
	  else if ( format[t] == "GQ" )
	    {

	      MetaInformation<GenMeta>::field( "GQ" , META_INT , 1 , "Genotype Quality score (phred-scaled)" );

	      for ( int per = 0 ; per < n ; per++ )
		{
		  uint8_t x;
		  read( x );
		  var(per).meta.set( "GQ" , (int)x );
		}

	    }
	  else if ( format[t] == "HQ" )
	    {
	      MetaInformation<GenMeta>::field( "HQ" , META_INT , 2 , "Haplotype Quality score (phred-scaled)" );

	      for ( int per = 0 ; per < n ; per++ )
		{
		  uint8_t x,y;
		  read( x );
		  read( y );
		  std::vector<int> hq(2);
		  hq[0]=x; hq[1]=y;
		  var(per).meta.set( "HQ" , hq );
		}

	    }
	  else if ( format[t] == "PL" )
	    {
	      
	      MetaInformation<GenMeta>::field( "PL" , META_INT , -1 , "Phred-scaled genotype likelihood" );

	      for ( int per = 0 ; per < n ; per++ )
		{
		  std::vector<int> pl(ngen); 
		  uint8_t pl0;
		  for (int g=0; g<ngen; g++) 
		    {
		      read( pl0 );		  
		      pl[g] = pl0;
		    }
		  var(per).meta.set( "PL" , pl );
		}
	      
	    }
	  else if ( format[t] == "SP" ) 
	    {

	      MetaInformation<GenMeta>::field( "SP" , META_INT , 1 , "Strand Bias p-value (bcftools)" );
	      
	      for ( int per = 0 ; per < n ; per++ )
		{
		  uint8_t x;
		  read( x );
		  var(per).meta.set( "SP" , (int)x );
		}
	    }
	  else
	    plog.warn( "unrecognise tag in BCF", format[t] );
	  
	} // next FORMAT-specified tag


      //
      // Check variant
      //

      std::cout << var << "\t" << var.consensus.meta << "\n";

//       for (int i=0;i<n; i++)
// 	std::cout << sample_names[i] << "\t" << var(i) << "\t" << var(i).meta << "\n";

    }
  
  BCF.close();
}

bool BCFReader::read( std::string & s )
{  
  int32_t l;
  read(l);
  char buf[l+1];
  BCF.read(buf , l );      
  buf[l] = '\0';  
  s = buf;
  return true; 
}

bool BCFReader::read( std::string & s , int32_t l )
{
  char buf[l+1];
  BCF.read(buf,l);      
  buf[l] = '\0';
  std::cout << "read buf [" << buf << "]\n";
  s = buf;
    return true; 
}

bool BCFReader::read( std::vector<std::string> & s  )
{
  // assume 1) int32_t with size of str
  //        2) str is NULL-padded, and so we split tino a vector

  s.clear();
  int32_t l;
  read(l);
  char buf[l+1];
  BCF.read(buf,l);      
  buf[l] = '\0';
  char * p = &(buf[0]);
  for (int j=0;j<l;j++)
    {
      if ( buf[j] == '\0' || j == l-1 ) 
	{
	  s.push_back( p );	  
	  p = &(buf[j+1]);
	}
    }
  return true; 
}

bool BCFReader::read( int32_t & i)
{
  BCF.read((char*)&i,sizeof(uint32_t));
  if ( endian == MACHINE_BIG_ENDIAN ) i = swap_endian(i);
  return true;
}

bool BCFReader::read( uint16_t & i)
{
  BCF.read((char*)&i,sizeof(int16_t));
  if ( endian == MACHINE_BIG_ENDIAN ) i = swap_endian(i);
  return true;
}

bool BCFReader::read ( uint8_t & i)
{
  BCF.read((char*)&i,sizeof(uint8_t));
  return true;
}

bool BCFReader::read(double & d)
{
  float f;
  BCF.read((char*)&f,sizeof(float));
  if ( endian == MACHINE_BIG_ENDIAN ) f = swap_endian(f);
  d = f;
  return true;
}






