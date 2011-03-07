#include "bcf.h"
#include "variant.h"
#include "vardb.h"
#include "gstore.h"

#include <iostream>

extern GStore * GP;

bool BCF::open()
{
  if ( file ) close();
  if ( readmode )
    file = bgzf_open( filename.c_str(), "r" );
  else
    file = bgzf_open( filename.c_str(), "w" );

  // does this file have correct EOF marker?
  // ignore for now...
  int i = bgzf_check_EOF( file );
  
  return file == NULL;
}

void BCF::close() 
{
  if ( file ) bgzf_close( file );
  file = NULL;
}


void BCF::read_header( VarDBase * v )
{
  // has a VARDB been specified?
  if ( v ) vardb = v;

  hdr.clear();
  
  std::vector<char> ch(4);
  if ( ! read( ch , 4 ) ) Helper::halt( "problem with format of BCF file (1) " );
  if ( ch[0] != 'B' || ch[1] != 'C' || ch[2] != 'F' || ch[3] != '\4' )
    Helper::halt( "problem with format of BCF file (1) " );   
  
  if ( ! read( hdr.seq_names ) ) Helper::halt( "problem with format of BCF header(2) " );   
  if ( ! read( hdr.sample_names ) ) Helper::halt( "problem with format of BCF header(3) " );   
  n = hdr.sample_names.size(); 

  std::string mtext;
  if ( ! read( mtext ) ) Helper::halt( "problem with format of BCF header(4) " );   
  // F not '.' for empty slot  
  hdr.meta_text = Helper::char_split( mtext , '\n' , false ); 
  
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
      if ( f == NULL ) Helper::halt( "internal error in BCF class, parsing header"  );
      
      // NULL means no SEQDB attached for now, i.e no REF checking
      VCFReader v( f , f->tag() , vardb , NULL );
      
      // Get file-ID (would have been created from VARDB)
      file_id = v.group_id();
      
      // Chromosome codes
      // ??
      
      // Individuals
      for (int i=0; i<n; i++)
	{
	  Individual ind( hdr.sample_names[i] );
	  vardb->insert( file_id , ind );
	}


      // Register BCF file in VARDB; also the # of individuals in the
      // BCF (i.e. so we do not need to parse the header when
      // accessing variant data from the BCF in future

      vardb->store_bcf_n( file_id , filename , n );

      plog << "added " << hdr.sample_names.size() << " individuals from BCF " << filename << "\n";

      // Meta-fields
      v.insert_meta( "##format=BCF4.0" );
      for (int i=0; i<hdr.meta_text.size(); i++)
	{	  	  
	  v.insert_meta( hdr.meta_text[i] );	  
	}
    }

}


std::string BCF::vcf_header() 
{
  std::string vcf = "##fileformat=VCFv4.0\n";
  vcf += "##source=plinkseq(BCF->VCF)\n";
  for (int h=0; h<hdr.meta_text.size(); h++)
    vcf += hdr.meta_text[h] + "\n";
  vcf += "#CHROM\tPOS\tID\tREF\tALT\tQUAL\tFILTER\tINFO";
  if ( hdr.sample_names.size() > 0 ) 
    {
      vcf += "\tFORMAT";
      for (int i = 0; i < n ; i++)
	vcf += "\t" + hdr.sample_names[i];
    }
  vcf += "\n";
  return vcf;
}

bool BCF::index_record( )
{
  // requires an attached VARDB
  if ( ! vardb ) return false;
  
  // track offset for this variant
  int64_t offset = tell();

  // we do not need to process information here

  int32_t seq_id; 
  int32_t bp;     
  float qual;     
  std::vector<std::string> mstr;   

  if ( ! read( seq_id ) ) return false;
  if ( ! read( bp ) ) return false;
  if ( ! read( qual ) ) return false;
  if ( ! read( mstr ) ) return false;
  
  std::vector<std::string> alt = Helper::char_split( mstr[2] , ',' );
  int nalt = alt.size() + 1;
  int ngen = (int) (nalt * (nalt+1) * 0.5);
  
  Variant var;
  var.chromosome( Helper::chrCode( hdr.seq_names[ seq_id ] ) );
  var.position( bp );
  var.stop( bp + mstr[1].size() - 1 );
  
  // TODO: change this to calculate a single # of bytes and read/skip those
  std::vector<std::string> format = Helper::char_split( mstr[5] , ':' );
  
  for ( int t = 0 ; t < format.size(); t++ )
    {
      
      // BCF recognised tags are:
      //                                       Length in bytes
      //                                       ----------------------
      //     DP     uint16_t                   sizeof(uint16_t) * n
      //     GL     float   (*ngen)            sizeof(float) * g * n 
      //     GT     uint8_t                    sizeof(uint8_t) * n
      //     GQ     uint8_t                    sizeof(uint8_t) * n
      //     HQ     uint8_t (*2)               sizeof(uint8_t) * 2 * n
      //     PL     uint8_t (*ngen)            sizeof(uint8_t) * 2 * n
      //     SP     uint8_t                    sizeof(uint8_t) * n
      //     misc   int32_t + char*            sizeof(int32_t) + value

      if ( format[t] == "DP" ) 
	{	  
	  for ( int per = 0 ; per < n ; per++ )
	    {
	      uint16_t dp; 
	      if ( ! read(dp) ) return false;
	    }
	  
	}
      else if ( format[t] == "GL" )
	{
	  for ( int per = 0 ; per < n ; per++ )
	    {
	      float gl;
	      for (int g=0; g<ngen; g++) 
		if ( ! read( gl ) ) return false;
	    }	  
	}
      else if ( format[t] == "GT" )
	{
	  for ( int per = 0 ; per < n ; per++ )
	    {
	      uint8_t gt;
	      if ( ! read( gt ) ) return false;
	    }
	}
      else if ( format[t] == "GQ" )
	{
	  for ( int per = 0 ; per < n ; per++ )
	    {
	      uint8_t x;
	      if ( ! read( x ) ) return false;
	    }
	  
	}
      else if ( format[t] == "HQ" )
	{
	  for ( int per = 0 ; per < n ; per++ )
	    {
	      uint8_t x,y;
	      if ( ! read( x ) ) return false;
	      if ( ! read( y ) ) return false;
	    }
	  
	}
      else if ( format[t] == "PL" )
	{	  
	  for ( int per = 0 ; per < n ; per++ )
	    {
	      uint8_t pl0;
	      for (int g=0; g<ngen; g++) 
		{
		  if ( ! read( pl0 ) ) return false;
		}
	    }
	}
      else if ( format[t] == "SP" ) 
	{
	  for ( int per = 0 ; per < n ; per++ )
	    {
	      uint8_t x;
	      if (! read( x ) ) return false;
	    }
	}
      else
	plog.warn( "unrecognise tag in BCF", format[t] );
      
    } // next FORMAT-specified tag
  
  //  std::cout << var << "\t" << offset << "\n";
  
  // Add offset to DB
  vardb->insert_bcf_index( file_id , var , offset );

  return true;
  
}

bool BCF::read_record( Variant & var , SampleVariant & svar , SampleVariant & svar_g )
{
  
  int32_t seq_id;                 // 'chromosome' code  
  int32_t bp;                     // position (BP1)  
  float qual;                     // quality score
  std::vector<std::string> mstr;  // ID REF ALT FILTER INFO FORMAT (null padded)
  
  if ( ! read( seq_id ) ) return false;
  if ( ! read( bp ) ) return false;
  if ( ! read( qual ) ) return false;
  if ( ! read( mstr ) ) return false;  
  if ( mstr.size() < 6 ) return false;

  // check that Variant spec (chr/bp) from VARDB matches what is in
  // the BCF (although chr. coding will be different -- can check that
  // bp matches though)
  
  if ( bp != var.position() )
    plog.warn( "mismatching physical position between VARDB and BCF" , 
	       Helper::int2str(bp) + " vs " + Helper::int2str( var.position() ) );
  
  // Read genotype; # alleles based on ALT tag, comma-sep values
  
  std::vector<std::string> alt = Helper::char_split( mstr[2] , ',' );
  int nalt = alt.size() + 1;
  int ngen = (int) (nalt * (nalt+1) * 0.5);
  

  //
  // Populate Sample Variant -- variant-level information
  //

  // chr, bp1, bp2 and name should already be populated (from VARDB)
  // so we can skip these Variant-level attributes here

  // Add SampleVariant core attributes

  svar.quality( qual );
  svar.reference( mstr[1] );
  svar.alternate( mstr[2] );
  svar.filter( mstr[3] );
  svar.info( mstr[4] );


  //
  // Add in genotype information
  //
  
  // BCF::n (sample size) will have been populated already (from VARDB)
  //                      *but*, we might not want to read it all in... (Masks...)
  

  svar_g.calls.size( n );
  
  // Genotype fields to expect based on 
  // use set_format() function in vcf.cpp
  
  std::vector<std::string> format = Helper::char_split( mstr[5] , ':' );


  // BCF recognised tags are:
  //                                       Length in bytes
  //                                       ----------------------
  //     DP     uint16_t                   sizeof(uint16_t) * n
  //     GL     float   (*ngen)            sizeof(float) * g * n 
  //     GT     uint8_t                    sizeof(uint8_t) * n
  //     GQ     uint8_t                    sizeof(uint8_t) * n
  //     HQ     uint8_t (*2)               sizeof(uint8_t) * 2 * n
  //     PL     uint8_t (*ngen)            sizeof(uint8_t) * g * n
  //     SP     uint8_t                    sizeof(uint8_t) * n
  //     misc   int32_t + char*            sizeof(int32_t) + value
  
  
  // Indicate that this SampleVariant contains a BCF-derived buffer
  // Fill buffer as direct raw byte sequence from BCF, for later parsing

  // set 'target', so that VMETA is also appropriately not handled as a BLOB
  svar.bcf = this;

  // but genotype specific info stays w/ genotypes
  svar_g.bcf_format = mstr[5];
  svar_g.bcf_genotype_buf.clear();

  int32_t buf_sz = 0;
  
  for ( int t = 0 ; t < format.size(); t++ )
    {

      int p = buf_sz;
      
      if ( format[t] == "DP" ) 
	{
	  buf_sz += n * sizeof(uint16_t);
	  svar_g.bcf_genotype_buf.resize( buf_sz );
	  read( &(svar_g.bcf_genotype_buf)[ p ] , n * sizeof(uint16_t) );
	}
      else if ( format[t] == "GL" ) 
	{
	  buf_sz += n * ngen * sizeof(float);
	  svar_g.bcf_genotype_buf.resize( buf_sz );
	  read( &(svar_g.bcf_genotype_buf)[ p ] , n * ngen * sizeof(float) );	  
	}
      else if ( format[t] == "PL" ) 
	{
	  buf_sz += n * ngen * sizeof(uint8_t);
	  svar_g.bcf_genotype_buf.resize( buf_sz );
	  read( &(svar_g.bcf_genotype_buf)[ p ] , n * ngen * sizeof(uint8_t) );	  
	}
      else if ( format[t] == "GT" ) 
	{
	  buf_sz += n * sizeof(uint8_t);
	  svar_g.bcf_genotype_buf.resize( buf_sz );
	  read( &(svar_g.bcf_genotype_buf)[ p ] , n * sizeof(uint8_t) );	  
	}
      else if ( format[t] == "GQ" )
	{
	  buf_sz += n * sizeof(uint8_t);
	  svar_g.bcf_genotype_buf.resize( buf_sz );
	  read( &(svar_g.bcf_genotype_buf)[ p ] , n * sizeof(uint8_t) );	  
	}
      else if ( format[t] == "HQ" ) 
	{
	  buf_sz += 2 * n * sizeof(uint8_t);
	  svar_g.bcf_genotype_buf.resize( buf_sz );
	  read( &(svar_g.bcf_genotype_buf)[ p ] , 2 * n * sizeof(uint8_t) );	  
	}
      else if ( format[t] == "SP" ) 
	{
	  buf_sz += n * sizeof(uint8_t);
	  svar_g.bcf_genotype_buf.resize( buf_sz );
	  read( &(svar_g.bcf_genotype_buf)[ p ] , n * sizeof(uint8_t) );	  
	}
      else 
	{
	  int32_t len;
	  read( len );
	  svar_g.bcf_genotype_buf.resize( buf_sz + sizeof(int32_t) + len );
	  uint8_t * d = (uint8_t*)len;
	  svar_g.bcf_genotype_buf[p] = d[0];
	  svar_g.bcf_genotype_buf[p+1] = d[1];
	  svar_g.bcf_genotype_buf[p+2] = d[2];
	  svar_g.bcf_genotype_buf[p+3] = d[3];
	  read( &(svar_g.bcf_genotype_buf)[ p+4 ] , len );
	}
    }
  
  // now we've finished building the BCF genotype buffer; this may get
  // expanded later (i.e. given a Mask, not everybody will be
  // extracted) necessarily.

  return true;
}


  
// 	  MetaInformation<GenMeta>::field( "DP", META_INT , 1, "Read depth" );	      
// 	  for ( int per = 0 ; per < n ; per++ )
// 	    {
// 	      uint16_t dp; 
// 	      read(dp);
// 	      svar_g.calls.genotype(per).meta.set( "DP" , (int)dp );		  
// 	    }
	  
// 	}
//       else if ( format[t] == "GT" )
// 	{
// 	  for ( int per = 0 ; per < n ; per++ )
// 	    {
// 	      uint8_t gt;
// 	      read( gt );
	      
// 	      // genotype encoded phase << 6 | allele1 << 3 | allele2 
	      
// 	      Genotype g( &var );
	      
// 	      g.null(   ( gt >> 7 ) & 1 );
// 	      g.phased( ( gt >> 6 ) & 1 );
	      
// 	      // TMP -- for now, ignore >2 alleles, but casting to bools
	      
// 	      g.pat( ( gt >> 3 ) & 7 );
// 	      g.mat( gt & 7 );
	      
// 	      svar_g.calls.genotype(per) = g;
	      
// 	    }
// 	}
//       else if ( format[t] == "GQ" )
// 	{
	  
// 	  MetaInformation<GenMeta>::field( "GQ" , META_INT , 1 , "Genotype Quality score (phred-scaled)" );
	  
// 	  for ( int per = 0 ; per < n ; per++ )
// 	    {
// 	      uint8_t x;
// 	      read( x );
// 	      svar_g.calls.genotype(per).meta.set( "GQ" , (int)x );
// 	    }
	  
// 	}
//       else if ( format[t] == "HQ" )
// 	{
// 	  MetaInformation<GenMeta>::field( "HQ" , META_INT , 2 , "Haplotype Quality score (phred-scaled)" );
	  
// 	  for ( int per = 0 ; per < n ; per++ )
// 	    {
// 	      uint8_t x,y;
// 	      read( x );
// 	      read( y );
// 	      std::vector<int> hq(2);
// 	      hq[0]=x; hq[1]=y;
// 	      svar_g.calls.genotype(per).meta.set( "HQ" , hq );
// 	    }
	  
// 	}
//       else if ( format[t] == "PL" )
// 	{
	  
// 	  MetaInformation<GenMeta>::field( "PL" , META_INT , -1 , "Phred-scaled genotype likelihood" );
	  
// 	  for ( int per = 0 ; per < n ; per++ )
// 	    {
// 	      std::vector<int> pl(ngen); 
// 	      uint8_t pl0;
// 	      for (int g=0; g<ngen; g++) 
// 		{
// 		  read( pl0 );		  
// 		  pl[g] = pl0;
// 		}
// 	      svar_g.calls.genotype(per).meta.set( "PL" , pl );
// 	    }
	  
// 	}
//       else if ( format[t] == "SP" ) 
// 	{
	  
// 	  MetaInformation<GenMeta>::field( "SP" , META_INT , 1 , "Strand Bias p-value (bcftools)" );
	  
// 	  for ( int per = 0 ; per < n ; per++ )
// 	    {
// 	      uint8_t x;
// 	      read( x );
// 	      svar_g.calls.genotype(per).meta.set( "SP" , (int)x );
// 	    }
// 	}
//       else
// 	plog.warn( "unrecognised tag in BCF", format[t] );
      
      
//     } // next FORMAT-specified tag
  
// //   for ( int per = 0 ; per < n ; per++ )
// //     std::cout << "in " << per << "\t" << svar_g.calls.genotype(per) << "\t" << svar_g.calls.genotype(per).meta << "\n";
  
//   return true;
  
// }



bool BCF::read_record( Variant & var , SampleVariant & svar , SampleVariant & svar_g, int64_t offset )
{
  seek( offset );
  return read_record( var , svar , svar_g );
}


bool BCF::vcf2bcf( const std::string & vcfname , const std::string & bcfname )
{

  // Get VCF
  if ( ! Helper::fileExists( vcfname ) )
    {
      plog.warn( "could not find VCF" , vcfname );
      return false;
    }
  InFile vcf( vcfname );
  
  // Open a new BCF for writing
  filename = bcfname;
  writing();
  open();
  hdr.clear();

  
  // For now, assume we'll encounter only 'normal' chromosome/sequence codes
  for (int c=1;c<=22;c++)
    {
      hdr.seq_names.push_back( "chr" + Helper::int2str( c ) );
      hdr.seq_names.push_back( Helper::int2str( c ) );
    }
  hdr.seq_names.push_back( "chrX" );
  hdr.seq_names.push_back( "X" );
  hdr.seq_names.push_back( "chrY" );
  hdr.seq_names.push_back( "Y" );
  hdr.seq_names.push_back( "chrM" );
  hdr.seq_names.push_back( "M" );
  
  // lookup code given string
  hdr.seq_map.clear();
  for (int i=0; i<hdr.seq_names.size(); i++)
    hdr.seq_map[ hdr.seq_names[i] ] = i;
  
  std::string mtext = "";

  // Start parsing VCF, first the header

  while ( ! vcf.eof() ) 
    {
      std::string l = vcf.readLine();
      if ( l == "" ) continue;
      if ( l.size() < 2 ) continue;
      // meta text
      if ( l.substr(0,2) == "##" ) mtext += l + "\n";
      // header
      else if ( l.substr(0,1) == "#" ) 
	{
	  
	  // silly fix, but get rid of trailing tab in GATK VCF header
	  const int len = l.size();
	  if ( l.substr( len - 1) == "\t" ) 
	    l = l.substr( 0 , len-1 );

	  std::vector<std::string> t = Helper::char_split( l , '\t' );
	  // skip forst 9 fields (0-8)
	  for (int i=9; i<t.size(); i++)
	    hdr.sample_names.push_back(t[i]);
	  n = hdr.sample_names.size();
	  break;
	}
    }
  
  // F means no empty fields (i.e. last row otherwise)
  hdr.meta_text = Helper::char_split( mtext , '\n' , false );
  
  // Write header info to BCF
  write_header();

  plog << "inserted header into BCF, " << n << " individuals\n";
  int inserted = 0;

  // process each variant row until end of file
  while ( ! vcf.eof() ) 
    {

      // read line from VCF
      std::string s = vcf.readLine();
      if ( s == "" ) continue;

      // Parse to make a Variant

      std::vector<std::string> tok
	= Helper::char_split(  s[ s.size()-1] == '\n' ?
			       s.substr( 0 , s.size() - 1 ) :
			       s , '\t' );

      Variant var( true );
      
      // valid VCF row?
      if ( tok.size() < 8 ) continue;

      // store chr code as per BCF header
      var.chromosome( Helper::chrCode( hdr.seq_map[ tok[0] ] ) );

      int pos;
      Helper::str2int( tok[1] , pos );
      var.position( pos );
      var.stop( tok[3].size() == 1 ? 0 : pos + tok[3].size() - 1 );

      // store as NULL if missing
      var.name( tok[2] == "." ? "" : tok[2] );
      
      var.consensus.reference( tok[3] );
      var.consensus.alternate( tok[4] );
      
      double qual;
      Helper::str2dbl( tok[5] , qual );
      
      var.consensus.quality(qual);
      
      var.consensus.filter( tok[6] == "." ? "" : tok[6] );
      var.consensus.info( tok[7] == "." ? "" : tok[7] );
      
      // no genotype data
      if ( tok.size() > 8 )
	{
      
	  // store format string, to be parsed later
	  var.consensus.bcf_format = tok[8];

	  // Parse the format specifier, and attach the the variant, so that
	  // genotypes can be called.
	  
	  VariantSpec * ps = SampleVariant::decoder.decode( var.consensus.bcf_format + " " + tok[3] + " " + tok[4] );
	  var.consensus.specification( ps );
	  
	  // Call genotypes, add to variant
	  int gcnt = 0;
	  for (int i=9;i<tok.size(); i++)
	    {
	      Genotype g = ps->callGenotype( tok[i] , &var ); 
	      var.consensus.calls.add( g );
	      ++gcnt;
	    }
  
	  
	  // Did we see the correct number of genotypes?
      
	  if ( gcnt != n )
	    {
	      plog.warn( "incorrect number of genotypes: " 
			 + Helper::int2str( gcnt) + " observed, " 
			 + Helper::int2str( n ) + " expected" ) ;
	      var.valid( false );
	    }
	  
	}
      
      // At this point, should be okay to write to BCF
      
      if ( var.valid() ) 
	write_record( var );
	
      if ( ++inserted % 1000 == 0 ) 
	{
	  plog.counter( "inserted " + Helper::int2str( inserted ) + " variants from VCF" ); 
	}
      
    } // next line in VCF
  
  plog << "inserted " << inserted << " variants from VCF\n";

  // close VCF, BCF  
  vcf.close();
  close();
}


bool BCF::write_header()
{
  if ( ! file ) return false;
  write( 'B' );
  write( 'C' );
  write( 'F' );
  write( '\4' );  
  write( hdr.seq_names );
  write( hdr.sample_names );
  std::string mtext = "";
  for ( int m=0;m<hdr.meta_text.size(); m++)
    {
      mtext += hdr.meta_text[m] + "\n";
    }
  write( mtext );
  return true;
}

bool BCF::write_record( const Variant & var )
{  
  
  // We assume this call only comes from the VCF to BCF function
  // i.e. and so the bcf_format variable is set, etc
  
  if ( ! file ) return false;
  
  // variant info

  write( hdr.seq_map[ Helper::chrCode( var.chromosome() ) ] );
  write( var.position() );
  write( (float)var.consensus.quality() );

  // ID REF ALT FILTER INFO FORMAT (null padded)  
  std::vector<std::string> m;
  m.push_back( var.name() == "." ? "" : var.name() );
  m.push_back( var.consensus.reference() );
  m.push_back( var.consensus.alternate() );
  m.push_back( var.consensus.filter() );
  m.push_back( var.consensus.info() );
  m.push_back( var.consensus.bcf_format );
  write( m );
 
  // Read genotype; # alleles based on ALT tag, comma-sep values
  // BCF::n should already be set 

  std::vector<std::string> alt = Helper::char_split( var.alternate() , ',' );
  int nalt = alt.size() + 1;
  int ngen = (int) (nalt * (nalt+1) * 0.5);
  
  // Add in genotype information, from consensus
  
  std::vector<std::string> format = Helper::char_split( var.consensus.bcf_format , ':' );
  

  // BCF recognised tags are:
  //                                       Length in bytes
  //                                       ----------------------
  //     DP     uint16_t                   sizeof(uint16_t) * n
  //     GL     float   (*ngen)            sizeof(float) * g * n 
  //     GT     uint8_t                    sizeof(uint8_t) * n
  //     GQ     uint8_t                    sizeof(uint8_t) * n
  //     HQ     uint8_t (*2)               sizeof(uint8_t) * 2 * n
  //     PL     uint8_t (*ngen)            sizeof(uint8_t) * g * n
  //     SP     uint8_t                    sizeof(uint8_t) * n
  //     misc   int32_t + char*            sizeof(int32_t) + value
  
  
  for ( int t = 0 ; t < format.size(); t++ )
    {

      if ( format[t] == "DP" ) // uint16_t
	{
	  for (int i=0; i<n; i++) write( (uint16_t) var(i).meta.get1_int( "DP" ) );
	}
      else if ( format[t] == "GL" ) // float * G
	{
	  for (int i=0; i<n; i++) 
	    {
	      std::vector<double> gl = var(i).meta.get_double( "GL" );
	      for (int j=0;j<gl.size();j++)
		write( (float)gl[j] ); 
	    }
	}
      else if ( format[t] == "PL" ) // uint8_t * G
	{
	  for (int i=0; i<n; i++) 
	    {
	      std::vector<int> pl = var(i).meta.get_int( "PL" );
	      for (int j=0;j<pl.size();j++)
		write( (uint8_t)pl[j] ); 
	    }
	}
      else if ( format[t] == "GT" ) // uint8_t
	{
	  for (int i=0; i<n; i++) write( var(i).bcf() );
	}
      else if ( format[t] == "GQ" ) // uint8_t
	{
	  for (int i=0; i<n; i++) write( (uint8_t) var(i).meta.get1_int( "GQ" ) );
	}
      else if ( format[t] == "HQ" ) // 2 * uint8_t
	{
	  for (int i=0; i<n; i++) 
	    {
	      std::vector<int> hq = var(i).meta.get_int( "HQ" );
	      for (int j=0;j<hq.size();j++)
		write( (uint8_t)hq[j] ); 
	    }
	}
      else if ( format[t] == "SP" ) // uint8_t
	{
	  for (int i=0; i<n; i++) write( (uint8_t) var(i).meta.get1_int( "SP" ) );
	}
      else 
	{	  
	  // Variable
	  std::vector<std::string> m(n);
	  for (int i=0;i<n;i++)
	    {
	      // get values as comma-delimited list (per individual)
	      m[i] = var(i).meta.as_string( format[t] , "," );
	    }
	  write( m );
	}
    }

  // end of writing this record 

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

inline bool BCF::read( int32_t & i )
{  
  int ret = bgzf_read(file,&i,sizeof(int32_t));
  if ( endian == MACHINE_BIG_ENDIAN ) i = swap_endian(i);
  return ret > 0 ;
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
  int32_t l;
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
  int32_t l;

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

inline void BCF::write( uint16_t i)
{
  uint16_t v = endian == MACHINE_BIG_ENDIAN ? swap_endian(i) : i ;
  bgzf_write(file,&v,sizeof(uint16_t));      
}

inline void BCF::write( int32_t i )
{
  int32_t v = endian == MACHINE_BIG_ENDIAN ? swap_endian(i) : i ;
  bgzf_write(file,&v,sizeof(int32_t));      
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
  write((int32_t)s.size());
  bgzf_write(file,s.c_str(),s.size());      
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
  //  std::cout << t.size() << " t=["<<t<<"]\n";
  write( (int32_t)t.size() );
  bgzf_write(file,t.c_str(),t.size());      
}


// Random access function ( seek / tell )

inline void BCF::seek( int64_t offset)
{
  bgzf_seek(file,offset,SEEK_SET);
}

inline int64_t BCF::tell()
{
  return bgzf_tell(file);
}  
