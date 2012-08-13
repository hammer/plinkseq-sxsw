#include "plinkseq/bcf.h"
#include "plinkseq/variant.h"
#include "plinkseq/vardb.h"
#include "plinkseq/gstore.h"
#include "plinkseq/meta.h"

#include <iostream>

extern GStore * GP;

bool BCF::get_next_typed_value( bcf2_typed_value * type , uint32_t * length )
 {

   // read the next byte, and determine type and length based on BCF2 spec.
   // return false is problem reading or invalid codes given
   
   // Bit 0-3    Type         Missing value      Description
   
   //  1         int8_t       0x80               signed 8-bit integer
   //  2         int16_t      0x8000             signed 16-bit integer
   //  3         int32_t      0x80000000         signed 32-bit integer 
   //  5         float        0x7F800001         IEEE 32-bit floating point number
   //  7         char         '\0'               character 
   
   char c;  
   if ( ! read( c ) ) return false;  

   uint8_t t = c & 0x0F;  // type
   uint8_t l = ( c & 0xF0 ) >> 4; // length

   //
   // FLAG types are defined by type 0 apparently?
   // 

   // if ( l == 0 ) 
   //   {
   //     *length = 0;
   //     *type = bcf2_flag;
   //     return true;
   //   }

   switch ( t ) { 
   case 0 :
     *type = bcf2_flag; 
     *length = 1;
     return true;
     break;
   case 1 :  
     *type = bcf2_int8_t;
     break;
   case 2 :
     *type = bcf2_int16_t;
     break;
   case 3 :
     *type = bcf2_int32_t;
     break;
   case 5 :
     *type = bcf2_float;
     break;
   case 7 :
     *type = bcf2_char;
     break;    
   default:       
     return false;
   }
   
   
   if ( l < 15 ) 
     {
       *length = l;
       return true;
     }
   
   // Otherwise, next byte determines length

   // Determine length: 0 = missing; 1 = scalar; 2-14 = vector length;
   // 15 = next typed integer gives length
   
   if ( l == 15 ) 
     {

       // expecting a typed integer to give array size

       char lc;       
       if ( ! read( lc ) ) return false;
       
       uint8_t lt =   lc & 0x0F;
       uint8_t ll = ( lc & 0xF0 ) >> 4 ;

       // doesn't make sense if length is not a scalar
       if ( ll != 1 ) return false; 

       if ( lt == bcf2_int8_t ) 
	 {
	   int8_t x;
	   if ( ! read( x ) ) return false;
	   l = (int)x;
	 }
       else if ( lt == bcf2_int16_t )
	 {
	   int16_t x;
	   if ( ! read( x ) ) return false;
	   l = (int)x;
	 }
       else if ( lt == bcf2_int32_t )
	 {
	   int32_t x;
	   if ( ! read( x ) ) return false;
	   l = (int)x;
	 }
       else 
	 return false;
     }
   
   *length = l;
   
   return true;
   
 }

uint32_t BCF::write_typing_byte( bcf2_typed_value type , uint32_t length )
{
  uint8_t t = type;
  if ( length < 15 ) 
    {
      uint8_t l = length;
      write( (char) ( ( l << 4 ) | t ) );
      return 1;
    }
  else 
    {
      write( (char) ( 0xF << 4 ) | t );
      if ( length < 128 ) 
	{
	  write( (char)( ( 0x1 << 4 ) | 0x01 ) );
	  write( (int8_t) length );
	  return 3;
	}
      else if ( length < 32767 ) 
	{
	  write( (char)( ( 0x1 << 4 ) | 0x02 ) );
	  write( (int16_t) length );
	  return 4;
	}
      else
	{
	  write( (char)( ( 0x1 << 4 ) | 0x03 ) );
	  write( (int32_t) length );
	  return 6;
	}
    }
  return 0;
}

uint32_t BCF::size_typed_value( bcf2_typed_value type , uint32_t length )
{  

  int l = 1;
  if ( length >= 15 )
    {
      if ( length < 127 ) l = 3;  
      else if ( length < 32767 ) l = 4;
      else l = 6;
    }

  if ( type == bcf2_int8_t || type == bcf2_char || type == bcf2_flag ) 
    return l + length ;
  
  if ( type == bcf2_int16_t ) 
    return l + 2 * length;

  if ( type == bcf2_int32_t || type == bcf2_float )
    return l + 4 * length;
  
  // void value has zero size other than typing bit
  if ( type == bcf2_void ) return l;

  return 0;
}


uint32_t BCF::size_typing_byte( bcf2_typed_value type , uint32_t length )
{  
  
  // similar to the above, although this is just the size of the typing byte

  int l = 1;
  if ( length >= 15 )
    {
      if ( length < 127 ) l = 3;  
      else if ( length < 32767 ) l = 4;
      else l = 6;
    }
  
  return l;
}

bool BCF::open()
{
  
  if ( file ) close();
  
  if ( readmode )
    file = bgzf_open( filename.c_str(), "r" );
  else
    file = bgzf_open( filename.c_str(), "w" );

  // does this file have correct EOF marker?
  
  if ( file == NULL || ( readmode && bgzf_check_EOF( file ) != 1 ) )
    Helper::halt( "problem opening BCF -- is this a BGZF BCF file?" );
  
  return file == NULL;
}


void BCF::close() 
{

  if ( file ) bgzf_close( file );

  file = NULL;

}


void BCF::read_header( VarDBase * v )
{
  
  //
  // has a VARDB been specified?
  //

  if ( v ) vardb = v;
  
  hdr.clear();

  
  //
  // BCF2 magic number
  //
  
  std::vector<char> ch(4);

  if ( ! read( ch , 4 ) ) 
    Helper::halt( "problem reading from BCF file" );
  
  if ( ch[0] != 'B' || ch[1] != 'C' || ch[2] != 'F' || ch[3] != '\2' )
    Helper::halt( "no magic number BCF/2 in header: is this a BCF2 file?" );

  
  //
  // Expecting standard VCF headers 
  //

  // length of all header text, including null padding

  uint32_t len;

  if ( bgzf_read( file, &len, sizeof(uint32_t) ) == 0 )
    Helper::halt( "problem reading header length from BCF" );
  
  std::string headers;
  
  if ( ! read( headers , len ) ) 
    Helper::halt( "problem reading headers from BCF" );
  
  // parse headers into lines ( F not '.' for empty slot )

  hdr.meta_text = Helper::char_split( headers , '\n' , false ); 
  
  // Assign to VarDB
  
  if ( vardb )
    {
  
      // add chrCodes and sample IDs to VARDB
      // and meta-information. 
      // we'll never want to look at this header
      // again when accessing genotype/variant data 

      // Use VCFReader class functions to process meta-information
      
      // BCF will be in file-index

      File * f = GP->fIndex.file( filename );

      if ( f == NULL ) 
	Helper::halt( "internal error in BCF class, parsing header"  );
      
      // NULL means no SEQDB attached for now, i.e no REF checking
      VCFReader vcf( f , f->tag() , vardb , NULL );
      
      // Get file-ID (would have been created from VARDB)
      file_id = vcf.group_id();
      
      // Individuals; should always be on last line of 

      if ( hdr.meta_text[ hdr.meta_text.size() - 1 ].substr( 0 , 6 ) == "#CHROM" ) 
	{
	  int n = 0;
	  Helper::char_tok tok( hdr.meta_text[ hdr.meta_text.size() - 1 ] , &n , '\t' );
	  for (int i=9;i<tok.size();i++)
	    hdr.sample_names.push_back( tok(i) );	      	    
	}
      else
	Helper::halt( "did not find a valid header row #CHROM ... in BCF: " + headers );
      
      //
      // store number of individuals, and put in VARDB
      //

      n = hdr.sample_names.size();
      
      for (int i=0; i<n; i++)
	{
	  Individual ind( hdr.sample_names[i] );
	  vardb->insert( file_id , ind );
	}
      
      
      //
      // build meta/string dictionary
      //
      
      
      // Register BCF file in VARDB; also the # of individuals in the
      // BCF (i.e. so we do not need to parse the header when
      // accessing variant data from the BCF in future

      // 2 == BCF mode

      vardb->store_bcf_n( file_id , filename , 2 , n );
      
      plog << "added " << hdr.sample_names.size() 
	   << " individuals from BCF " << filename << "\n";
      

      //
      // Meta-fields
      //
      
      // add 'PASS' implicitly in the first offset
      hdr.dictionary( "PASS" );      
      

      // first pass for 'contig' codes
      // then 'dictionary'       
      // then eval as normal meta-information
      
      for (int i=0; i<hdr.meta_text.size(); i++)
	{
	  
	  //	  std::cout << "[" << hdr.meta_text[i] << "]\n";
	  
	  // assume '##contig=<ID='
	  // assume '##contig=<ID=chr20>'
	  // OR     '##contig=<ID=chr20,Description="dndndn">'

	  if ( hdr.meta_text[i].substr( 0 , 9 ) == "##contig=" ) 
	    {
	      size_t pos = hdr.meta_text[i].find( "," ) ;
	      if ( pos == std::string::npos ) 
		pos = hdr.meta_text[i].find( ">" ) ;
	      if ( pos != std::string::npos ) 
		{
		  std::string chr = hdr.meta_text[i].substr( 13 , pos - 13 );
		  plog << "adding to contig dictionary [" << chr << "]\n";
		  hdr.contig_dictionary( chr );							   
		}
	    }
	  
	  // ##dictionary=

	  if ( hdr.meta_text[i].substr( 0 , 13 ) == "##dictionary=" ) 
	    {
	      std::vector<std::string> dic = Helper::parse( hdr.meta_text[i].substr(14) , "," );
	      for (int i = 0 ; i < dic.size() ; i++) 
		{
		  plog << "adding dictionary entry [" << dic[i] << "]\n";
		  hdr.dictionary( dic[i] );
		}
	    }
	  
	  
	  

	  // put in VARDB and parse/store meta-information header as usual
	  
	  std::string name = vcf.insert_meta( hdr.meta_text[i] );

	  if ( name != "" ) 
	    {
	      plog << "adding " << name << " to dictionary\n";
	      hdr.dictionary( name );
	    }

	}
      
      // add all to VARDB, to store

      int ds = hdr.dictionary();
      for (int i = 0 ; i < ds ; i++ ) 
	{
	  std::string s ; 
	  hdr.dictionary( i , &s );
	  vardb->insert_bcf_dictionary( file_id , 1 , i , s );
	}

      // and contigs
      ds = hdr.contig_dictionary();
      for (int i = 0 ; i < ds ; i++ ) 
	{
	  std::string s ; 
	  hdr.contig_dictionary( i , &s );
	  vardb->insert_bcf_dictionary( file_id , 2 , i , s );
	}

      //
      // show dictionary 
      //

      if ( false )
	{
	  int ds = hdr.dictionary();
	  for (int i = 0 ; i < ds ; i++ ) 
	    {
	      std::string s ; 
	      hdr.dictionary( i , &s );
	      std::cout << "dictionary = " << i << " = [" << s << "]\n";
	    }
	  
	  // and contigs
	  ds = hdr.contig_dictionary();
	  for (int i = 0 ; i < ds ; i++ ) 
	    {
	      std::string s ; 
	      hdr.contig_dictionary( i , &s );
	      std::cout << "contig dictionary = " << i << " = [" << s << "]\n";
	    }
	}
      
    }
  
}




bool BCF::index_record( )
{
  
  // requires an attached VARDB
  if ( ! vardb ) return false;
  
  // track offset for this variant
  int64_t offset = tell();
  
  // we do not need to process information here
  
  uint32_t l_shared, l_indiv; 

  // these two give the entire record length
  if ( ! read( l_shared ) ) return false;  // EOF
  if ( ! read( l_indiv ) ) return false;   // should not fail
  

  // read basic chromosome code information   
  int32_t chr;
  int32_t pos;
  int32_t rlen; // length of reference sequence

  float   qual; // ignored
  uint32_t n_allele_info;
  uint32_t n_fmt_sample;

  std::string ID; // variant ID, is a BCF 'typed str'
  
  if ( ! read( chr ) ) Helper::halt( "problem indexing record (chr)" );
  if ( ! read( pos ) ) Helper::halt( "problem indexing record (pos)" );
  if ( ! read( rlen ) ) Helper::halt( "problem indexing record (rlen)" );

  
  if ( ! read( qual ) ) Helper::halt( "problem indexing record (QUAL)" );
  if ( ! read( n_allele_info ) ) Helper::halt( "problem indexing record (n_allele_info)" );
  if ( ! read( n_fmt_sample ) ) Helper::halt( "problem indexing record (n_fmt_sample)" );

  if ( ! read_typed_string( ID ) ) Helper::halt( "problem indexing record (ID) " );
  
  // can ignore the rest of the record now
  
  // get chromosome code
  std::string chrcode = "";
  if ( ! hdr.contig_dictionary( chr , &chrcode ) ) return false;

  
  // note, pos is 0-based.  The pos references an offset into a
  // necessarily-specified reference/contig. Currently we are not
  // enforcing this, so just +1

  ++pos;
  

  //
  // Create a stub Variant
  //

  Variant var;
  var.chromosome( chrcode );
  var.position( pos );
  var.stop( pos + rlen - 1 );
  var.name( ID == "" ? "." : ID );

  //
  // Add as offset to VARDB
  //

  vardb->insert_bcf_index( file_id , var , offset );
  
  //
  // Skip to one past the end of record
  //
  
  // size of rest of record, minus what we've already read 
  int sz = l_shared + l_indiv - 6*4 - size_typed_value( bcf2_char , var.name().size() );

  std::vector<char> dummy( sz  );
  if ( ! read( dummy , sz   ) ) return false;
  
  return true;
  
}

bool BCF::read_record( Variant & var , SampleVariant & svar , SampleVariant & svar_g )
{


  //
  // Reading record (meta-information and/or genotype, depending on
  // mask settings) from BCF record. The variant-level does goes into
  // 'svar'.  The genotype level data goes into svar_g.  These two
  // SampleVariants may or may not point to the same structure.
  //

  //
  // At this point, we should already be pointed at the correct position, and 
  // we will have read the variant stub from the VARDB index. 
  //


  // these two give the entire record length  
  uint32_t l_shared, l_indiv; 

  if ( ! read( l_shared ) ) return false;  
  if ( ! read( l_indiv ) ) return false;  
  
  // basic chromosome code information: check it matches VARDB
  int32_t chr;
  int32_t pos;
  int32_t rlen; // length of reference sequence
  float   qual; 
  uint32_t n_allele_info;
  uint32_t n_fmt_sample;
  std::string ID; 
  
  if ( ! read( chr ) ) Helper::halt( "problem indexing record (chr)" );
  if ( ! read( pos ) ) Helper::halt( "problem indexing record (pos)" );
  if ( ! read( rlen ) ) Helper::halt( "problem indexing record (rlen)" );
  
  if ( ! read( qual ) ) Helper::halt( "problem indexing record (QUAL)" );
  if ( ! read( n_allele_info ) ) Helper::halt( "problem indexing record (n_allele_info)" );
  if ( ! read( n_fmt_sample ) ) Helper::halt( "problem indexing record (n_fmt_sample)" );
  if ( ! read_typed_string( ID ) ) Helper::halt( "problem indexing record (ID) " );
  


  //  std::string chrcode = hdr.dictionary( chr ) ;
  
  // note, pos is 0-based.  The pos references an offset into a
  // necessarily-specified reference/contig. Currently we are not
  // enforcing this, so just +1.

  ++pos;
  
  //
  // Check that Variant spec (chr/bp) from VARDB matches what is in
  // the BCF. Although the chr-coding will be different, so let's just
  // check base-position
  //
  
  if ( pos != var.position() )
    Helper::halt( "mismatching physical position between VARDB and BCF: " +
		  Helper::int2str( pos ) + " vs " + var.coordinate() );

  //
  // QUAL : we use -9 as the missing code for quality scores
  //
  
  
  // TODO ; check for IEEE NaN code as well as BCF2 missing code

  svar.quality( qual == (float)0x7F800001 ? -9.0 : qual );

  //
  // Expand out allele encoding, etc
  //

  // Number of alleles (including REF)
  uint16_t n_allele = n_allele_info >> 16;

  // Number of key-value INFO pairs
  int n_info   = n_allele_info & 0x0000FFFF;;

  // Number of FORMAT fields
  int n_fmt     = n_fmt_sample >> 24;

  // Number of samples
  int n_sample  = n_fmt_sample & 0x00FFFFFF; 
  
  // save, so that read_genotypes() knows about size
  last_record_n_fmt = n_fmt;
  last_record_n_sample = n_sample;

  //
  // Read ALLELE codes
  //

  std::string ref_allele = "";
  std::string alt_allele = "";
  
  if ( ! read_typed_string( ref_allele) ) return false;
  // should not happen, but just in case, set explicitly missing
  // (should fail downstream anyway)
  if ( ref_allele == "" ) ref_allele == ".";
  
  if ( n_allele == 2 ) 
    {
      if ( ! read_typed_string( alt_allele ) ) return false;
    }
  else
    {
      std::vector<std::string> alt_alleles;
      for (int a=1;a<n_allele;a++) 
	{
	  std::string aa;
	  if ( ! read_typed_string( aa ) ) return false;
	  if ( aa == "" ) aa == ".";
	  alt_alleles.push_back( aa );
	}      
      alt_allele = Helper::stringize( alt_alleles , "," );
    }
  
  svar.reference( ref_allele );
  svar.alternate( alt_allele );
  
  
  //
  // Read FILTER : typed vec, where items are defined in the dictionary
  //
  

  std::vector<int> filtervec;
  if ( ! read_FILTER( filtervec ) ) return false;
  
  std::vector<std::string> fltr;  
  if ( filtervec.size() == 0 ) 
    {
      svar.filter( "." );
    }
  else
    {
      
      for (int f=0;f<filtervec.size();f++)
	{
	  
	  // MISSING code (i.e. no offset into dictionary)
	  if ( filtervec[f] == -1 ) 
	    { 
	      fltr.push_back(".");
	  }
	  else
	  {
	    std::string value;      
	    if ( ! hdr.dictionary( filtervec[f] , &value ) )
	      plog.warn( "FILTER value not in dictionary" , Helper::int2str( filtervec[f] ) );
	    else
	      {
		fltr.push_back( value );      
	      }
	  }
	  
	// TODO: silly to merge then parse again later -- so add filter()
	// that accepts a pre-tokenized vector
	  
	  svar.filter( Helper::stringize( fltr , ";" ) );
	}
    }
  

  //
  // INFO fields
  //
 

  // std::map<bcf2_typed_value,std::string> types;
  // types[ bcf2_int8_t ]  = "int8_t";
  // types[ bcf2_int16_t ] = "int16_t";
  // types[ bcf2_int32_t ] = "int32_t";
  // types[ bcf2_float ]   = "float";
  // types[ bcf2_char ]    = "char";
  // types[ bcf2_void ]    = "void";

  std::string info_str;

  for (int i=0;i<n_info;i++)
    {
      
      // 1) Get key into dictionary (a typed INT)
      int d = 0;
      bool m = false; // missing value

      if ( ! read_typed_int( d , m ) ) return false;

      if ( ! m ) 
	{

	  std::string mkey;
	  
	  if ( ! hdr.dictionary( d , &mkey ) ) 
	    Helper::halt( "could not find dictionary entry for BCF that should exist in header" );	  
	  
	  if ( info_str != "" ) info_str += ";";	  
	  info_str += mkey;
	  
	  // 2) Get value (a typed value).  This should already have
	  //    been registered in MetaInformation<VarMeta> by virtue
	  //    of being in the VARDB header.  Nonetheless, for now,
	  //    make as a a simple string and re-parse, just to ensure
	  //    existing checks are all in place. Yes, I know this
	  //    misses half the point of having a nice binary-encoded
	  //    format in the first place...
	  

	  bcf2_typed_value type;
	  uint32_t length;  

	  if ( ! get_next_typed_value( &type , &length ) ) return false;
	  
   	  
   //  1         int8_t       0x80               signed 8-bit integer
   //  2         int16_t      0x8000             signed 16-bit integer
   //  3         int32_t      0x80000000         signed 32-bit integer 
   //  5         float        0x7F800001         IEEE 32-bit floating point number
   //  7         char         '\0'               character 

   // ?9         string                          string (vector of chars)
	  
	  // Flag? 
	  
	  meta_index_t midx = MetaInformation<VarMeta>::field( mkey );

	  // can check length specification here, etc

	  
	  if ( midx.mt == ::META_FLAG ) 
	    {
	      // handle Flag variables differently
	      if ( type != bcf2_int8_t )
		Helper::halt( "expecting BCF implementation to encode Flags as int8_t 0x11 0x01" );
	      int8_t x;
	      if ( ! read(x) ) return false; // read this value, ignore	      
	    }	  
	  else if ( length == 0 ) 
	    {
	      // a missing, non-flag value
	      info_str += "=.";
	    }
	  else 
	    {
	      
	      // for non-missing, non-flag values
	      
	      if ( type == bcf2_int8_t ) 
		{
		  std::vector<int8_t> x( length);
		  if ( ! read(x) ) return false;
		  for (int i=0;i<length;i++)
		    info_str += (i>0?",":"=") + Helper::int2str(x[i]);		     
		}
	      else if ( type == bcf2_int16_t ) 
		{
		  std::vector<int16_t> x( length);
		  if ( ! read(x) ) return false;
		  for (int i=0;i<length;i++)
		    info_str += (i>0?",":"=") + Helper::int2str(x[i]);		     		  
		}
	      else if ( type == bcf2_int32_t )
		{
		  std::vector<int32_t> x( length);
		  if ( ! read(x) ) return false;
		  for (int i=0;i<length;i++)
		    info_str += (i>0?",":"=") + Helper::int2str(x[i]);		     		  
		}
	      else if ( type == bcf2_float )
		{
		  std::vector<float> x( length);
		  if ( ! read(x) ) return false;
		  for (int i=0;i<length;i++)
		    info_str += (i>0?",":"=") + Helper::dbl2str(x[i]);		     		  
		}
	      else if ( type == bcf2_char )
		{
		  std::vector<char> x( length + 1 ); // not null delimited
		  if ( ! read( x , length ) ) return false;
		  x[length] = '\0'; 

		  // now have to parse on commas, to check not a string vector
		  int n = 0;
		  Helper::char_tok tok( &x[0] , &n , ',' );
		  for (int i=0;i<n;i++)
		    info_str += ( i>0 ? "," : "=" ) + std::string( tok(i) ) ; 		  

		}
	      else 
		{
		  return false;
		}
	    }
	}
    }
  


  // TODOL : as above -- make a version of into() that accepts a pre-parsed 
  // list...
  
  svar.info( info_str );
 

  // Let the SampleVariant know where to come back to look for genotype data
  
  svar_g.set_pointer_to_bcf( this , tell() );
  
  // All done for now; leave file pointer at start of genotype records (if any) 
  
  return true;
    
 
}



bool BCF::read_record( Variant & var , 
		       SampleVariant & svar , 
		       SampleVariant & svar_g, 
		       int64_t offset )
{
  seek( offset );
  return read_record( var , svar , svar_g );
}



bool BCF::write_header()
{
  
  // Need to populate BCF dictionaries here too, with four types of entry:
  //  VarMeta, VarFilterMeta, GenMeta and contig/chr-codes.

  // add 'PASS' as first item.
  hdr.dictionary( "PASS" ); 
  
  // Magic keyword
  if ( ! file ) return false;  
  write( 'B' );
  write( 'C' );
  write( 'F' );
  write( '\2' );
  
  // VCF header lines as plain text

  std::stringstream ss;
  
  ss << "##fileformat=" << PLINKSeq::CURRENT_VCF_VERSION() << "\n"
    "##source=pseq\n"    
    "##encoding=BCF2\n"
     << MetaInformation<VarMeta>::headers( ) 
     << MetaInformation<GenMeta>::headers( META_GROUP_GEN ) 
     << MetaInformation<VarFilterMeta>::headers( META_GROUP_FILTER );
  
  std::map<meta_name_t,meta_index_t> m_var = MetaInformation<VarMeta>::dump_types();
  std::map<meta_name_t,meta_index_t> m_gen = MetaInformation<GenMeta>::dump_types();
  std::map<meta_name_t,meta_index_t> m_fltr = MetaInformation<VarFilterMeta>::dump_types();
  
  std::map<meta_name_t,meta_index_t>::iterator ii = m_var.begin();
  while ( ii != m_var.end() ) 
    {
      if ( MetaMeta::display( ii->first ) )
	hdr.dictionary( ii->first );
      ++ii;
    }
  
  ii = m_gen.begin();
  while ( ii != m_gen.end() ) 
    {
      if ( MetaMeta::display( ii->first ) )
	hdr.dictionary( ii->first );
      ++ii;
    }
  
  ii = m_fltr.begin();
  while ( ii != m_fltr.end() ) 
    {
      if ( MetaMeta::display( ii->first ) )
	hdr.dictionary( ii->first );
      ++ii;
    }
  

  // also ensure the GT is present
  if ( m_gen.find( "GT" ) == m_gen.end() )
    {
      hdr.dictionary( "GT" );
      ss << "##FORMAT=<ID=GT,Number=1,Type=String,Description=\"Genotype\">\n";
    }

  
  // Expecting contig dictionary also
  // Get a list of all chromosome strings
  
  std::set<std::string> chrs = GP->vardb.fetch_all_chr_names();
  std::set<std::string>::iterator jj = chrs.begin();
  while ( jj != chrs.end() ) 
    {
      ss << "##contig=<ID=" << *jj << ">\n";
      hdr.contig_dictionary( *jj );
      ++jj;
    }
  
  const int n_sample = GP->indmap.size();

  ss << "#CHROM\tPOS\tID\tREF\tALT\tQUAL\tFILTER\tINFO\tFORMAT";  
  for ( int i=0; i<n_sample; i++) ss << "\t" << GP->indmap(i)->id();
  
  // header text should be null-terminated
  // AND l_text is meant to include any NULL padding
  ss << "\0";

  write( ss.str() );

  // add a padding NULL

  return true;

}


bool BCF::write_record( const Variant & var )
{  
  
  std::cout << "writing record " << var << " at " << tell() << "\n";

  // Write the entire record to a buffer, as we need to calculate length 
  
  // l_shared   uint32_t
  // l_indiv    uint32_t
 
  // CHROM      int32_t
  // POS        int32_t
  // rlen       int32_t
  // QUAL       float
  // n_allele_info uint32_t
  // n_fmt_sample uint32_t

  // 


  // 
  // Calculate l_shared
  //

  uint32_t l_shared = 6 * 4 ;  
  // core, mandatory fields (6, listed above) -- does not include l_shared or l_indiv
  
  uint32_t l_indiv  = 0;


  //
  // Calculate length of record, whilst checking types, etc
  //


  //
  // Prep INFO and FORMAT fields up-front
  //

  const int n_sample = var.size();
  const int n_allele = var.n_alleles();

 	
 
  // Variant ID
  l_shared += size_typed_value( bcf2_char , var.name().size() );
  
  // Number of alleles
  const std::vector<Allele> & alleles0 = var.consensus.alleles; 
  for (int a=0;a<n_allele;a++)
    l_shared += size_typed_value( bcf2_char , alleles0[a].name().size() );


  // FILTERs

  std::set<std::string> filters = var.meta_filter();
  int lrg_fltr = 0;
  std::set<int> dfilters;
  std::set<std::string>::iterator ii = filters.begin();
  while ( ii != filters.end() )
    {
      int d;
      if ( ! hdr.dictionary( *ii , & d ) ) 
	Helper::halt( "FILTER value not found in dictionary: " + *ii );
      dfilters.insert( d );
      if ( d > lrg_fltr ) lrg_fltr = d;
      ++ii;
    }
  
  if ( lrg_fltr < 128 )
    l_shared += size_typed_value( bcf2_int8_t , dfilters.size() );
  else if ( lrg_fltr < 32767 ) 
    l_shared += size_typed_value( bcf2_int16_t , dfilters.size() );
  else 
    l_shared += size_typed_value( bcf2_int32_t , dfilters.size() );


 //
  // Get INFO elements
  //

  std::vector<std::string> keys1 = var.consensus.meta.keys();
  
  std::map<std::string,bcf2_typed_value> infotype;
  
  for (int m=0; m < keys1.size(); m++)
    {
      
      //
      // dictionary key
      //

      int d;
      if ( ! hdr.dictionary( keys1[m], & d ) ) 
	Helper::halt( "INFO value not found in dictionary: " + keys1[m] );
      
      if ( d < 128 )
	l_shared += size_typed_value( bcf2_int8_t , 1 );
      else if ( d < 32767 ) 
	l_shared += size_typed_value( bcf2_int16_t , 1 );
      else 
	l_shared += size_typed_value( bcf2_int32_t , 1 );
      

      //
      // actual value
      //

      mType mt = MetaInformation<VarMeta>::type( keys1[m] );
      
      int sz = 0;
      
      if ( mt == META_INT ) 
	{

	  infotype[ keys1[m] ] = bcf2_int8_t;

	  // get max range
	  std::vector<int> x = var.consensus.meta.get_int( keys1[m] );
	  sz = x.size();
	  for (int ii=0; ii<sz; ii++)
	    {
	      if ( x[ii] < -127   || x[ii] > 127   ) infotype[ keys1[m] ] = bcf2_int16_t;
	      if ( x[ii] < -32767 || x[ii] > 32767 ) infotype[ keys1[m] ] = bcf2_int32_t;
	    }
	  
	  // figure out length
	  l_shared += size_typed_value( infotype[ keys1[m] ] , sz );
	  
	}
      else if ( mt == META_FLAG ) 
	{
	  infotype[ keys1[m] ] = bcf2_flag;
	  l_shared += size_typed_value( bcf2_int8_t , 1 );	  
	}
      else if ( mt == META_BOOL ) 
	{
	  infotype[ keys1[m] ] = bcf2_int8_t;
	  l_shared += size_typed_value( bcf2_int8_t , var.consensus.meta.get_size( keys1[m] ) );
	}
      else if ( mt == META_FLOAT )
	{
	  infotype[ keys1[m] ] = bcf2_float;
	  l_shared += size_typed_value( bcf2_float , var.consensus.meta.get_size( keys1[m] ) );
	}
      else if ( mt == META_TEXT )
	{
	  infotype[ keys1[m] ] = bcf2_char;
	  // ugly -- have to get size of string representing a delimited version of this:
	  int sz = Helper::stringize( var.consensus.meta.get_string( keys1[m] ) , "," ).size();
	  l_shared += size_typed_value( bcf2_char , sz );	  
	}


    }



  //
  // Some basic counts for stuff coming up
  //

  const int n_info   = infotype.size();

   //
  // Get FORMAT elements
  //

  // Scans through all individuals, in case new tags have been added.
  // This is quite inefficient, but not clearly a better way to do
  // this, if we want to be sure that any modifications to the data
  // post-loading are preserved...  
  
  // Actually, in this instance, we also need to determine the most 
  // efficient integer type to use to store the data, so perhaps not 
  // so bad: do this at this step also.  Also, we can calculate the 
  // expected length of the record, l_indiv.

  std::set<std::string>                  allkeys;
  std::map<std::string,bcf2_typed_value> inttype;
  std::map<std::string,int>              maxlength;


  for (int i = 0 ; i < n_sample; i++)
    {

      std::vector<std::string> keys = var(i).meta.keys();

      for (unsigned int j=0; j<keys.size(); j++) 
	{
	  if ( MetaMeta::display( keys[j] ) ) 
	    {

	      // if an integer?
	      if ( inttype.find( keys[j] ) == inttype.end() )
		{
		  if ( MetaInformation<GenMeta>::type( keys[j] ) == META_INT ) 
		    inttype[ keys[j] ] = bcf2_int8_t ;
		}
	      
	      if ( inttype.find( keys[j] ) != inttype.end() ) 
		{
		  std::vector<int> x = var(i).meta.get_int( keys[j] );
		  for (int ii=0;ii<x.size();ii++)
		    {
		      if ( x[ii] < -127   || x[ii] > 127 )   inttype[ keys[j] ] = bcf2_int16_t;
		      if ( x[ii] < -32767 || x[ii] > 32767 ) inttype[ keys[j] ] = bcf2_int32_t;		      
		    }
		}	      
	      
	      // get length; if a variable-length vector, get maximum 
	      // length by checking each individual
	      meta_index_t midx = MetaInformation<GenMeta>::field( keys[j] );

	      if ( midx.len == -1 ) 
		{
		  if ( midx.mt == META_INT ) 
		    {
		      std::vector<int> x = var(i).meta.get_int( keys[j] );
		      if ( x.size() > maxlength[ keys[j] ] )
			maxlength[ keys[j] ] = x.size();
		    }
		  else if ( midx.mt == META_FLOAT ) 
		    {
		      std::vector<double> x = var(i).meta.get_double( keys[j] );
		      if ( x.size() > maxlength[ keys[j] ] )
			maxlength[ keys[j] ] = x.size();
		    }
		  else if ( midx.mt == META_TEXT )
		    {
		      std::vector<std::string> x = var(i).meta.get_string( keys[j] );

		      // NOTE: I believe this is the correct implementation for a string-vector, although 
		      // seems would be better to have an explicit typing to avoid parsing comma-delimiters

		      int sz = Helper::stringize( x , "," ).size();
		      if ( x.size() > maxlength[ keys[j] ] )
			maxlength[ keys[j] ] = sz;
		    }
		  else if ( midx.mt == META_BOOL )
		    {
		      std::vector<bool> x = var(i).meta.get_bool( keys[j] );
		      if ( x.size() > maxlength[ keys[j] ] )
			maxlength[ keys[j] ] = x.size();
		    }
		  else if ( midx.mt == META_FLAG )  
		    {
		      // technically not allowed as FORMAT fields in VCF4.1?
		      // TODO: until checking/resolving the above, treat as scalar		      
		      maxlength[ keys[j] ] = 1;
		    }
		  else 
		    Helper::halt( "all FORMAT tags must be properly defined to make a BCF" );
		}
	      
	      // store key
	      allkeys.insert( keys[j] );
	    }
	}
    }
  
  

  // the above will take a while, but assuming write-once, read-many
  // times will hopefully be worth it

  // final list of FMTs

  std::vector<std::string> fmts;
  std::vector<int>         l_fmts;

  fmts.push_back( "GT" );         // always include a genotype field (?)
  l_fmts.push_back( 1 );

  std::set<std::string>::iterator i = allkeys.begin();
  while ( i != allkeys.end() )
    {
      fmts.push_back( *i );
      if ( maxlength.find( *i ) != maxlength.end() )
	{
	  // if variable length, we should have tracked above,
	  // to get the max per-individual; below, all vectors 
	  // will be padded to give simmilar length
	  l_fmts.push_back( maxlength[ *i ] );
	}
      else // should be a fixed-length tag
	{
	  meta_index_t midx = MetaInformation<GenMeta>::field( *i );
	  l_fmts.push_back( midx.len );
	}

      ++i;
    }
  
  
  // calculate l_indiv

  for (int f=0;f<fmts.size();f++)
    {
      // size of typing INT

      int df;
      if ( ! hdr.dictionary( fmts[f] , &df ) )
	Helper::halt( "undefined FORMAT field " + fmts[f] );      

      if      ( df < 127 )   l_indiv += size_typed_value( bcf2_int8_t , 1 ) ;
      else if ( df < 32767 ) l_indiv += size_typed_value( bcf2_int16_t , 1 ) ;
      else                   l_indiv += size_typed_value( bcf2_int32_t , 1 ) ;

      // 2) type of FORMAT field
      // 3) the actual FORMAT-field data on n_sample individuals

      bool genotype = fmts[f] == "GT";

      meta_index_t midx;
      if ( genotype )
	{
	  midx.mt = META_TEXT;
	  midx.len = 1;
	}
      else 
	midx = MetaInformation<GenMeta>::field( fmts[f] );
      
      int length = l_fmts[f];
      
      // NOTE: this currently hard-codes a biallelic allele w/ < 127 alleles (bcf2_int8_t)
      
      if ( genotype )
	{	  
	  l_indiv += 1;  // typing byte, should ~always be single byte
	  l_indiv += n_sample * 2; 
	}
      else 
	{
	  bcf2_typed_value t = inttype[ fmts[f] ];
	  l_indiv += size_typing_byte( t , length );

	  int s = 1;
	  if ( t == bcf2_int16_t ) s = 2;
	  else if ( t == bcf2_int32_t || t == bcf2_float ) s = 4;
	  if ( t == bcf2_int16_t ) 
	    l_indiv += n_sample * 2 * length;
	  else if ( t == bcf2_int32_t || t == bcf2_float ) 
	    l_indiv += n_sample * 4 * length;
	  else 
	    l_indiv += n_sample * length;	  
	}
    }
 

  //
  // Reset/save counters
  // 


  int e_shared = l_shared;
  int e_indiv = l_indiv;

  l_shared = 0;
  l_indiv = 0;



  //
  // Write out
  // 

  write( e_shared );
  write( e_indiv );

  //
  // Core info
  //

  const int n_fmt    = fmts.size();

  const uint32_t n_allele_info = n_allele << 16 | n_info ;
  const uint32_t n_fmt_sample  = n_fmt << 24    | n_sample ;
    
  int32_t contig;
  if ( ! hdr.contig_dictionary( Helper::chrCode( var.chromosome() ) , &contig ) )
    Helper::halt( "could not find chromosome/contig " + Helper::chrCode( var.chromosome() ) );
  

  write( (int32_t)contig ); 
  l_shared += 4;

  // Note: BCF uses 0-based position; as internally, we always use
  // 1-based need to adjust here

  write( (int32_t)var.position() - 1 );
  l_shared += 4;

  write( (int32_t)var.length() );
  l_shared += 4;
  
  float q = var.consensus.quality();
  if ( q < 0 ) q = (float)0x7F800001;
  write( q );
  l_shared += 4;
  
  
  write( n_allele_info );
  l_shared += 4;

  write( n_fmt_sample );
  l_shared += 4;


  //
  // Variant ID
  //

  l_shared += write_typed_string( var.name() );
  
  
  //
  // REF / ALT allele list
  //

  const std::vector<Allele> & alleles = var.consensus.alleles; 

  for (int a=0;a<n_allele;a++)
    l_shared += write_typed_string( alleles[a].name() );
  

  
  //
  // Write FILTERs
  //

  
  if ( lrg_fltr < 128 )        // int8_t
    {
      std::vector<int8_t> f;
      std::set<int>::iterator ii = dfilters.begin();
      while ( ii != dfilters.end() )
	{
	  f.push_back( *ii );
	  ++ii;
	}
      l_shared += write_typed_vec( f );
    }
  else if ( lrg_fltr < 32767 ) // int16_t 
    {
      std::vector<int16_t> f;
      std::set<int>::iterator ii = dfilters.begin();
      while ( ii != dfilters.end() )
	{
	  f.push_back( *ii );
	  ++ii;
	}
      l_shared += write_typed_vec( f );
    }
  else                         // in32_t
    {
      std::vector<int32_t> f;
      std::set<int>::iterator ii = dfilters.begin();
      while ( ii != dfilters.end() )
	{
	  f.push_back( *ii );
	  ++ii;
	}
      l_shared += write_typed_vec( f );
    }
  


  //
  // INFO fields
  // 
  
  std::map<std::string,bcf2_typed_value>::iterator jj = infotype.begin();
  while ( jj != infotype.end() ) 
    {


      // write as dictionary code
      int d;

      if ( ! hdr.dictionary( jj->first , &d ) ) 
	Helper::halt( "could not find INFO tag in BCF dictionary" );


      //
      // show dictionary 
      //

      if ( false ) 
	{
	  int ds = hdr.dictionary();
	  for (int i = 0 ; i < ds ; i++ ) 
	    {
	      std::string s ; 
	      hdr.dictionary( i , &s );
	      std::cout << "dict " << i << " " << s << "\n";
	    }
	  
	  // and contigs
	  ds = hdr.contig_dictionary();
	  for (int i = 0 ; i < ds ; i++ ) 
	    {
	      std::string s ; 
	      hdr.contig_dictionary( i , &s );
	      std::cout << "contig dict " << i << " " << s << "\n";
	    }
	}
	  
      l_shared += write_typed_int( d );


      bcf2_typed_value t = jj->second;
      
      if ( t == bcf2_int8_t ) 
	{
	  std::vector<int8_t> x;
	  std::vector<int> y = var.consensus.meta.get_int( jj->first );
	  for (int i=0;i<y.size();i++) x.push_back( y[i] );
	  l_shared += write_typed_vec( x ); 
	}
      else if ( t == bcf2_int16_t )
	{
	  std::vector<int16_t> x;
	  std::vector<int> y = var.consensus.meta.get_int( jj->first );
	  for (int i=0;i<y.size();i++) x.push_back( y[i] );
	  l_shared += write_typed_vec( x ); 
	}
      else if ( t == bcf2_int32_t )
	{
	  std::vector<int32_t> x;
	  std::vector<int> y = var.consensus.meta.get_int( jj->first );
	  for (int i=0;i<y.size();i++) x.push_back( y[i] );
	  l_shared += write_typed_vec( x ); 
	}
      else if ( t == bcf2_float )
	{
	  std::vector<float> x;
	  std::vector<double> y = var.consensus.meta.get_double( jj->first );
	  for (int i=0;i<y.size();i++) x.push_back( y[i] );
	  l_shared += write_typed_vec( x ); 
	}
      else if ( t == bcf2_char )
	{	  
	  std::vector<std::string> y = var.consensus.meta.get_string( jj->first );
	  
	  // i.e. current a vector of strings is just a single string (a vector of chars)
	  //      using comma-delimitation, that will have to be parsed upon reading

	  std::string x = Helper::stringize( y , "," );

	  l_shared += write_typed_string( x ); 	  

	}
      else if ( t == bcf2_flag )
	{
	  // Use 0x11 0x01 as righthandside value to represent a FLAG is set
	  write( (char) 0x11 );
	  write( (char) 0x01 );
	  l_shared += 2;
	}
      
      ++jj;
    }

   

  //
  // Finally, write out individual level genotype data
  //
  
  
  for (int f=0;f<fmts.size();f++)
    {
    
      // 1) name of FORMAT field
      
      int df;
      if ( ! hdr.dictionary( fmts[f] , &df ) )
	Helper::halt( "undefined FORMAT field " + fmts[f] );      

      l_indiv += write_typed_int( df );

      
      // 2) type of FORMAT field
      // 3) finally, print the actual FORMAT-field data on n_sample individuals

      bool genotype = fmts[f] == "GT";

      if ( (!genotype) && ! MetaInformation<GenMeta>::exists( fmts[f] ) ) 
	Helper::halt( "undefined FORMAT tag [ " + fmts[f] + " ] -- cannot write to BCF" );

      meta_index_t midx;
      if ( genotype )
	{
	  midx.mt = META_TEXT;
	  midx.len = 1;
	}
      else 
	midx = MetaInformation<GenMeta>::field( fmts[f] );
      
      int length = maxlength.find( fmts[f] ) != maxlength.end() ? maxlength[ fmts[f] ] : midx.len ;
      
      if ( genotype )
	{

	  // TODO: always write as int8_t (for now) and assume diploidy...

	  l_indiv += write_typing_byte( bcf2_int8_t , 2 );

	  for (int i=0;i<n_sample;i++)
	    {
	      
	      const Genotype & g = var(i);	      

	      bool haploid = g.haploid();
	      bool phased = g.phased();
	      
	      if ( g.null() )
		{
		  write( (int8_t) ( 0x00<<1 ) );
		  write( (int8_t) ( haploid ? 0x80 : 0x00<<1 | phased ) );
		}
	      else
		{

		  int allele1 = g.acode1();
		  int allele2 = g.acode2();
	      
		  // Note -- this is hard coded for now, meaning <127 alleles
		  write( (int8_t) ( (allele1+1)<<1 ) );
		  write( (int8_t) ( haploid ? 0x80 : (allele2+1)<<1 | phased ) );
		}
	      
	      l_indiv += 2;
	    }
	}
      else if ( midx.mt == META_INT ) 
	{
	  bcf2_typed_value t = inttype[ fmts[f] ];
	  
	  l_indiv += write_typing_byte( t , length );

	  for (int i=0;i<n_sample;i++)
	    {
	      std::vector<int> d = var(i).meta.get_int( fmts[f] );
	      if ( t == bcf2_int8_t )
		{
		  for (int k=0;k<length;k++)
		    {
		      if ( k < d.size() ) write( (int8_t) d[k] );
		      else write( (int8_t) 0x80 ); // pad w/ MISSING
		    }
		  l_indiv += length;		  
		}
	      else if ( t == bcf2_int16_t )
		{
		  for (int k=0;k<length;k++)
		    {
		      if ( k < d.size() ) write( (int16_t) d[k] );
		      else write( (int16_t) 0x8000 ); 
		    }
		  l_indiv += length * 2 ;
		}
	      else if ( t == bcf2_int32_t )
		{
		  for (int k=0;k<length;k++)
		    {
		      if ( k < d.size() ) write( (int32_t) d[k] );
		      else write( (int32_t) 0x80000000 ); 
		    }
		  l_indiv += length * 4 ;
		}

	    }
	}
      else if ( midx.mt == META_FLOAT )
	{
	  l_indiv += write_typing_byte( bcf2_float , length );

	  for (int i=0;i<n_sample;i++)
	    {
	      std::vector<double> d = var(i).meta.get_double( fmts[f] );	      
	      for (int k=0;k<length;k++)
		{
		  if ( k < d.size() ) write( (float) d[k] );
		  else write( (float) 0x7F800001 ); 
		}
	      l_indiv += length * 4 ;
	    }
	}
      else if ( midx.mt == META_TEXT ) 
	{
	  l_indiv += write_typing_byte( bcf2_char , length );	  
	  
	  for (int i=0;i<n_sample;i++)
	    {
	      std::vector<std::string> d = var(i).meta.get_string( fmts[f] );	      
	      std::string s = Helper::stringize( d , "," );
	      for (int k=0;k<length;k++)
		{
		  if ( k < d.size() ) write( (char) s[k] );
		  else write( (char) '\0' );
		}
	      l_indiv += length;
	    }
	}
      else if ( midx.mt == META_BOOL ) 
	{
	  l_indiv += write_typing_byte( bcf2_int8_t , length );	  

	  for (int i=0;i<n_sample;i++)
	    {
	      std::vector<bool> d = var(i).meta.get_bool( fmts[f] );	      
	      for (int k=0;k<length;k++)
		{
		  if ( k < d.size() ) write( (int8_t) d[k] );
		  else write( (int8_t) 0x80 );
		}
	      l_indiv += length;
	    }
	}
      
    }


  
  
  // 
  // Check that we wrote out correct amount
  // 

  if ( l_shared != e_shared ) 
    Helper::halt( "internal error in formating BCF : l_shared != e_shared, " 
		  + Helper::int2str(l_shared) + " " + Helper::int2str( e_shared ) );
  
  if ( l_indiv != e_indiv ) 
    Helper::halt( "internal error in formating BCF : l_indiv != e_indiv, " 
   		  + Helper::int2str(l_indiv) + " " + Helper::int2str( e_indiv ) );
  
  
  // Have written record

  return true;
  
}


//  Input / output functions 

inline bool BCF::read( char & c )
{
  return bgzf_read(file,&c,sizeof(char)) > 0;
}

inline bool BCF::read( std::vector<char> & buf, int l)
{
  return bgzf_read(file,&buf[0],l) > 0;
}

inline bool BCF::read( int8_t & i)
{
  return bgzf_read(file,&i,sizeof(int8_t)) > 0;
}

inline bool BCF::read( uint8_t & i)
{
  return bgzf_read(file,&i,sizeof(uint8_t)) > 0;
}

inline bool BCF::read( uint8_t * i, int32_t len )
{
  return bgzf_read( file, i, len ) > 0;
}

inline bool BCF::read( uint16_t & i)
{  
  int ret = bgzf_read(file,&i,sizeof(uint16_t));
  if ( endian == MACHINE_BIG_ENDIAN ) i = swap_endian(i);
  return ret > 0;
}

inline bool BCF::read( int16_t & i)
{  
  int ret = bgzf_read(file,&i,sizeof(int16_t));
  if ( endian == MACHINE_BIG_ENDIAN ) i = swap_endian(i);
  return ret > 0;
}

inline bool BCF::read( int32_t & i )
{  
  int ret = bgzf_read(file,&i,sizeof(int32_t));
  if ( endian == MACHINE_BIG_ENDIAN ) i = swap_endian(i);
  return ret > 0 ;
}

inline bool BCF::read( uint32_t & i )
{  
  int ret = bgzf_read(file,&i,sizeof(uint32_t));
  if ( endian == MACHINE_BIG_ENDIAN ) i = swap_endian(i);
  return ret > 0 ;
}

inline bool BCF::read( std::vector<int8_t> & i )
{
  for (int x=0;x<i.size();x++)
    {
      int8_t t;
      if ( ! read( t ) ) return false;
      i[x] = t;
    }
  return true;
}

inline bool BCF::read( std::vector<int16_t> & i )
{
  for (int x=0;x<i.size();x++)
    {
      int16_t t;
      if ( ! read( t ) ) return false;
      i[x] = t;
    }
  return true;
}

inline bool BCF::read( std::vector<int32_t> & i )
{
  for (int x=0;x<i.size();x++)
    {
      int32_t t;
      if ( ! read( t ) ) return false;
      i[x] = t;
    }
  return true;
}

inline bool BCF::read( std::vector<float> & i )
{
  for (int x=0;x<i.size();x++)
    {
      float t;
      if ( ! read( t ) ) return false;
      i[x] = t;
    }
  return true;
}


inline bool BCF::read( float & i)
{
  int ret;
  if ( endian == MACHINE_BIG_ENDIAN )
    {
      uint32_t a;
      ret = bgzf_read(file,&a,sizeof(uint32_t));  
      i = swap_float(a);
    }
  else
    ret = bgzf_read(file,&i,sizeof(float));  
 
  return ret > 0;
}

inline bool BCF::read( double & i)
{
  int ret;
  if ( endian == MACHINE_BIG_ENDIAN )
    {
      uint64_t a;
      ret = bgzf_read(file,&a,sizeof(uint64_t));  
      i = swap_double(a);
    }
  else
    ret = bgzf_read(file,&i,sizeof(double));  
  return ret > 0;
}

inline bool BCF::read( std::string & s, int l )
{
  char buf[l+1];
  int ret = bgzf_read(file,buf,l);
  buf[l] = '\0';
  s = buf;
  return ret > 0;
}

inline bool BCF::read( std::string & s )
{
  uint32_t l;
  read(l);
  char buf[l+1];
  int ret = bgzf_read(file,buf,l);
  buf[l] = '\0';
  s = buf;
  return ret > 0;
}

inline bool BCF::read( std::vector<std::string> & s )
{
  // assume 1) int32_t with size of str
  //        2) str is NULL-padded, and so we split tino a vector
  //        3) convert missing code to '.' (NO)
  
  s.clear();
  uint32_t l;

  if ( ! read(l) ) return false;
  
  char buf[l+1];
  int ret = bgzf_read(file,buf,l);
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
  return ret > 0;
}


// Output functions

inline void BCF::write( char i )
{
  bgzf_write(file,&i,sizeof(char));
}

inline void BCF::write(const std::vector<char> & buf, int l)
{      
  write((int32_t)l);
  bgzf_write(file,&buf[0],l);      
}

inline void BCF::write( uint8_t i)
{
  bgzf_write(file,&i,sizeof(uint8_t));      
}

inline void BCF::write( int8_t i)
{
  bgzf_write(file,&i,sizeof(int8_t));      
}

inline void BCF::write( uint16_t i)
{
  uint16_t v = endian == MACHINE_BIG_ENDIAN ? swap_endian(i) : i ;
  bgzf_write(file,&v,sizeof(uint16_t));      
}

inline void BCF::write( int16_t i)
{
  int16_t v = endian == MACHINE_BIG_ENDIAN ? swap_endian(i) : i ;
  bgzf_write(file,&v,sizeof(int16_t));      
}

inline void BCF::write( int32_t i )
{
  int32_t v = endian == MACHINE_BIG_ENDIAN ? swap_endian(i) : i ;
  bgzf_write(file,&v,sizeof(int32_t));      
}

inline void BCF::write( uint32_t i )
{
  uint32_t v = endian == MACHINE_BIG_ENDIAN ? swap_endian(i) : i ;
  bgzf_write(file,&v,sizeof(uint32_t));      
}

inline void BCF::write( float i )
{
  if ( endian == MACHINE_BIG_ENDIAN )
    {
      uint32_t a = swap(i);
      bgzf_write(file,&a,sizeof(uint32_t));      
    }
  else
    bgzf_write(file,&i,sizeof(float));      
}

inline void BCF::write( double i )
{
  if ( endian == MACHINE_BIG_ENDIAN )
    {
      uint64_t a = swap(i);
      bgzf_write(file,&a,sizeof(uint64_t));      
    }
  else
    bgzf_write(file,&i,sizeof(double));      
}

inline void BCF::write( const std::string & s )
{
  write((uint32_t)s.size());
  bgzf_write(file,s.c_str(),s.size()); // will not include terminal \0 
}

void BCF::write( const std::vector<std::string> & s )
{
  // write NULL-padded string of total length l
  std::string t = "";
  for (int i=0; i<s.size(); i++)
    {
      if ( i > 0 ) t += '\0';
      t += s[i];
    }
  write( (uint32_t)t.size() );
  bgzf_write(file,t.c_str(),t.size());      
}


// Random access function ( seek / tell )

inline void BCF::seek( int64_t offset)
{
  if ( bgzf_seek(file,offset,SEEK_SET) ) 
    Helper::halt( "internal error -- problem seeking to position in BCF/compressed VCF" );
}

inline int64_t BCF::tell()
{
  return bgzf_tell(file);
}  


inline bool BCF::read_typed_string( std::string & s )
{

  // a BCF2 'typed str' read function
  // read typing byte

  bcf2_typed_value type;
  uint32_t length;  
  if ( ! get_next_typed_value( &type , &length ) ) return false;

  // null string
  if ( length == 0 ) { s = ""; return true; }

  // for non-null values, expecting a char type
  if ( type != bcf2_char ) return false;
  
  // read it
  if ( ! read( s , length ) ) return false;  

  return true;
}


inline bool BCF::read_typed_int( int & i , bool & missing )
{
  // a BCF2 'typed int' read function  
  // read typing byte
  bcf2_typed_value type;
  uint32_t length;  

  if ( ! get_next_typed_value( &type , &length ) ) return false;

  // null value
  if ( length == 0 || type == bcf2_void ) 
    { 
      i = 0;
      missing = true;
    }
  
  missing = false;
  if ( type == bcf2_int8_t )
    {
      int8_t x;
      if ( ! read( x ) ) return false;    
      i = x;
    }
  else if ( type == bcf2_int16_t )
    {
      int16_t x;
      if ( ! read( x ) ) return false;    
      i = x;
    }
  else if ( type == bcf2_int32_t )
    {
      int32_t x;
      if ( ! read( x ) ) return false;    
      i = x;
    }
  else 
    return false;

  return true;
}




inline bool BCF::read_FILTER( std::vector<int> & s )
{
  
  s.clear();
  
  // expecting a vector of offsets (i.e. for the FILTER) 
  
  bcf2_typed_value type;
  uint32_t length;  

  if ( ! get_next_typed_value( &type , &length ) ) return false;
  
  // null string means 'missing' (-1)
  if ( length == 0 ) { s.push_back( -1 ); return true; } 
  
  if ( type == bcf2_int8_t )
    {
      std::vector<int8_t> x( length );
      if ( ! read( x ) ) return false;
      for (int i=0;i<length;i++) s.push_back( x[i] );
    }
  else if ( type == bcf2_int16_t )
    {
      std::vector<int16_t> x( length );
      if ( ! read( x ) ) return false;
      for (int i=0;i<length;i++) s.push_back( x[i] );
    }
  else if ( type == bcf2_int32_t )
    {
      s.resize( length );
      if ( ! read( s ) ) return false;      
    }
  else
    {
      // encountered something other than an int-vector, flag failure
      return false;
    }


  return true;
}




bool BCF::read_genotypes( int64_t offset , 
			  IndividualMap * align , 
			  SampleVariant * target , 
			  SampleVariant * source , 
			  Mask * mask )
{


  // move to genotype position

  seek( offset );
  
  //
  // Make room for genotypes
  //
  
  int n_actual = align ? align->size() : last_record_n_sample ;
  int n_buffer = last_record_n_sample;

  target->calls.size( n_actual );

  // Create mapping of slot-to-target
  
  


  //  TODO: doesn't this already exist in IndividualMap ?  Surely we
  //  do not want to recreate this every variant?...
  
  std::vector<int> s2t( n_buffer );
  
  for ( int i=0; i < n_buffer ; i++ )
    {
      int slot = i;	  
      if ( align )
	{
	  // CHECK: don't think both are needed here. (?)
	  slot = align->sample_remapping( source->fileset()  , i ) ;
	  if ( align->flat() ) slot = align->get_slot( source->fileset()  , slot );
	}	  
      
      s2t[ i ] = slot;
    }      
  
  
  //
  // Read FORMAT and GENOTYPE fields
  //

  for (int f=0;f< last_record_n_fmt;f++)
    {
      
      int d = 0;      // dictionary offset  

      bool m = false; // missing value
      
      if ( ! read_typed_int( d , m ) ) 
	return false;
      
      if ( m ) return false;
      
      std::string mkey;
      if ( ! hdr.dictionary( d , &mkey ) ) 
	Helper::halt( "could not find BCF dictionary entry " + Helper::int2str(d) ); 
      
      bool genotype = mkey == "GT";
      
      //
      // next -- type of following data
      //
      
      bcf2_typed_value type;

      uint32_t length;  

      if ( ! get_next_typed_value( &type , &length ) ) return false;
      
      if ( genotype ) // sanity checks
	{
	  
	  if ( length < 1 || length > 2 ) 
	    Helper::halt( "cannot handle ploidy other than 2 or 1 currently" );

	  // TODO -- presumably this could actually be higher, if needed, 
	  //         but for now let's assume it isn't but fail noisily.

	  if ( type != bcf2_int8_t ) 
	    Helper::halt( "nonstandard GT encoding, other than int8_t" );
	}


      //
      // Read individual level data 
      //
      
      for (int i=0;i< last_record_n_sample;i++)
	{
	  
	  bool include_individual = s2t[ i ] != -1 ;
	  
	  if ( genotype ) // read the actual GT field 
	    {
	      
	      int8_t allele1 , allele2; 
	      bool phased;
	      
	      // allele encoding = '(allele+1)<<1|phased’ 
	      // where allele is set to -1 if the allele in GT is a dot ‘.’ 
	      // vector is padded with missing values if the GT has fewer ploidy.
	      
	      if ( ! read( allele1 ) ) return false;
	      allele1 >>= 1;
	      
	      if ( length == 2 )
		{
		  if ( ! read( allele2 ) ) return false;
		  phased = allele2 & 0x01;		  		  
		  allele2 >>= 1;
		}
	      
	      // add Genotype (and meta-information) to SampleVariant?

	      if ( include_individual ) 
		{
		  
		  Genotype g;
		  
		  if      ( allele1 == 0 ) g.null( true );
		  else if ( length == 1 )  g.genotype( allele1-1 );
		  else if ( phased )       g.phased_genotype( allele1-1 , allele2-1 );
		  else                     g.genotype( allele1-1 , allele2-1 );
		  
		  target->calls.add( g , s2t[ i ] );
		  
		}
	      
	    }
	  else // ... read genotype meta-information
	    {

	      if ( type == bcf2_int8_t ) 
		{
		  std::vector<int8_t> x( length );
		  if ( ! read( x ) ) return false;

		  if ( include_individual && length != 0 )
		    {
		      if ( length == 1 && x[0] != (int8_t)0x80 ) 
			{
			  target->calls.genotype( s2t[i] ).meta.set( mkey , x[0] );
			}
		      else 
			{
			  // TODO: clean-up unnecessary copying into temps
			  std::vector<int> t;
			  for (int ii=0;ii<length;ii++)		    
			    if ( x[ii] != (int8_t)0x80 ) t.push_back( x[ii] );
			  if ( t.size() ) 
			    target->calls.genotype( s2t[i] ).meta.set( mkey , t );
			}
		    }
		  
		}
	      else if ( type == bcf2_int16_t )
		{
		  std::vector<int16_t> x( length );
		  if ( ! read( x ) ) return false;
		  
		  if ( include_individual ) 
		    {
		      // TODO: clean-up unnecessary copying into temps
		      std::vector<int> t;
		      for (int ii=0;ii<length;ii++)		    
			if ( x[ii] != (int16_t)0x8000 ) t.push_back( x[ii] );
		      target->calls.genotype( s2t[i] ).meta.set( mkey , t );
		    }

		}
	      else if ( type == bcf2_int32_t )
		{
		  std::vector<int32_t> x( length );
		  if ( ! read( x ) ) return false;

		  if ( include_individual ) 
		    {
		      // TODO: clean-up unnecessary copying into temps
		      std::vector<int> t;
		      for (int ii=0;ii<length;ii++)		    
			if ( x[ii] != (int32_t)0x80000000 ) t.push_back( x[ii] );
		      target->calls.genotype( s2t[i] ).meta.set( mkey , t );
		    }

		}
	      else if ( type == bcf2_float ) 
		{
		  std::vector<float> x( length );
		  if ( ! read( x ) ) return false;
		  if ( include_individual ) 
		    {
		      // TODO: clean-up unnecessary copying into temps
		      if ( length == 1 && x[0] != (float)0x7F800001 )
			target->calls.genotype( s2t[i] ).meta.set( mkey , x[0] );
		      else
			{
			  std::vector<double> t;
			  for (int ii=0;ii<length;ii++)		    
			    if ( x[ii] != (float)0x7F800001 ) t.push_back( x[ii] );
			  target->calls.genotype( s2t[i] ).meta.set( mkey , t );
			}
		    }
		}
	      else if ( type == bcf2_char ) 
		{
		  
		  // TODO: check this handles string vectors appropriately 
		  
		  std::vector<std::string> x( length );
		  if ( ! read( x ) ) return false;
		 
		  if ( include_individual )
		    {
		      target->calls.genotype( s2t[i] ).meta.set( mkey , x );		      
		    }
		}
	      else
		return -1;
	    }
	  
	}
      
    }
  
  return true;
}




inline uint32_t BCF::write_typed_int( int i , bool missing )
{

  if ( missing )                      { write( (char)0x01 ); return 1; } // missing 8-bit int

  // reserve lowest -128 for 'missing'
  if ( i < 128 && i > -128 )          { write( (char)0x11); write( (int8_t)i ); return 2;  }  
  else if ( i < 32767 && i > -32767 ) { write( (char)0x12); write( (int16_t)i ); return 3; }
  else                                { write( (char)0x13); write( (int32_t)i ); return 5; }
  return 0;
}


inline uint32_t BCF::write_typed_float( float f , bool missing )
{
  if ( missing ) { write( (int8_t)0x05 ); return 1; }
  write( (int8_t)0x15 );  // scalar float
  write( (float) f );
  return 5;
}

inline uint32_t BCF::write_typed_string( const std::string & s , bool missing )
{

  if ( missing ) { write( (char)0x07 ); return 1; }

  // string == vector of 'char'

  uint32_t rv = 0;
  if ( s.size() < 15 )
    {
      write( (char)( 0x07 | ( (int8_t)( s.size() << 4 ) ) ) );
      ++rv;
    }
  else
    {
      write( (char)( 0x07 | (int8_t)( 0xF << 4 ) ) );
      ++rv;
      rv += write_typed_int( s.size() );
    }
  
  // no trailing \0 in BCF2
  for (int i=0;i<s.size();i++) 
    write( (char)s[i] );  
  rv += s.size();

  return rv;
}


inline uint32_t BCF::write_typed_vec( const std::vector<int8_t>    & x , const std::vector<bool> * missing )
{

  
  if ( x.size() == 0 ) // null
    {
      write( (char)0x01 ); // void/missing
      return 1;
    }
  
  uint32_t rv = 0;

  if ( x.size() < 15 )
    {
      write( (char)( 0x01 | ( (int8_t)x.size() << 4 ) ) );
      ++rv;
    }
  else
    {
      write( (char)( 0x01 | ( 0xF << 4 ) ) );
      ++rv;
      rv += write_typed_int( x.size() );
    }
  
  if ( missing ) 
    {
      for (int i=0;i<x.size();i++) 
	if ( (*missing)[i] ) write( (int8_t) 0x80 );
	else write( (int8_t)x[i]);
    }
  else
    for (int i=0;i<x.size();i++) 
      {
	write( (int8_t)x[i] );  
      }

  rv += x.size();

  return rv;
}

inline uint32_t BCF::write_typed_vec( const std::vector<int16_t>  &  x , const std::vector<bool> * missing )
{
  if ( x.size() == 0 ) // null
    {
      write( (char)0x02 ); // void/missing
      return 1;
    }
  
  uint32_t rv = 0;

  if ( x.size() < 15 )
    {
      write( (char)( 0x02 | ( (int8_t)x.size() << 4 ) ) );
      ++rv;
    }
  else
    {
      write( (char)( 0x02 | ( 0xF << 4 ) ) );
      ++rv;
      rv += write_typed_int( x.size() );      
    }

  if ( missing ) 
    {
      for (int i=0;i<x.size();i++) 
	if ( (*missing)[i] ) write( (int16_t) 0x8000 );
	else write( (int16_t)x[i]);
    }
  else
    for (int i=0;i<x.size();i++) 
      write( (int16_t)x[i] );  
  
  rv += x.size() * 2 ;
  
  return rv;
}

inline uint32_t BCF::write_typed_vec( const std::vector<int32_t>  & x  , const std::vector<bool> * missing )
{

  if ( x.size() == 0 ) // null
    {
      write( (char)0x03 ); // void/missing
      return 1; 
    }

  uint32_t rv = 0;

  if ( x.size() < 15 )
    {
      write( (char)( 0x03 | ( (int8_t)x.size() << 4 ) ) );
      ++rv;
    }
  else
    {
      write( (char)( 0x03 | ( 0xF << 4 ) ) );
      ++rv;
      rv += write_typed_int( x.size() );
    }

  if ( missing ) 
    {
      for (int i=0;i<x.size();i++) 
	if ( (*missing)[i] ) write( (int32_t) 0x80000000 );
	else write( (int32_t)x[i]);
    }
  else
    for (int i=0;i<x.size();i++) 
      write( (int32_t)x[i] );  

  rv += x.size() * 4;

  return rv;
}

inline uint32_t BCF::write_typed_vec( const std::vector<float>    & x  , const std::vector<bool> * missing )
{
  if ( x.size() == 0 ) // null
    {
      write( (char)0x05 ); // void/missing
      return 1;
    }

  uint32_t rv = 0;
  if ( x.size() < 15 )
    {
      write( (char)( 0x05 | ( (int8_t)x.size() << 4 ) ) );
      ++rv;
    }
  else
    {
      write( (char)( 0x05 | ( 0xF << 4 ) ) );
      ++rv;
      rv += write_typed_int( x.size() );
    }

  if ( missing ) 
    {
      for (int i=0;i<x.size();i++) 
	if ( (*missing)[i] ) write( (float) 0x7F800001 );
	else write( (float)x[i]);
    }
  else
    for (int i=0;i<x.size();i++) 
      write( (float)x[i] );  
  
  rv += x.size() * 4 ; 
  
  return rv;
}


// note-- current confusion over string as vector of char, vs. how to define a vector of strings?
inline uint32_t BCF::write_typed_vec( const std::vector<std::string> & x  , const std::vector<bool> * missing )
{
  // make a single comma-delimited string??
  // TODO: handle individual 'missing' value elements, although this is really a mute point in the context 
  // of a string?  (i.e. zero length char vector == missing).  Alternative is \0 

  std::string s = Helper::stringize( x , "," );
  return write_typed_string( s ); // ignores 'missing' status here...
  
}


