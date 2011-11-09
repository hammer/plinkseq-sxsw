#include "bcf.h"
#include "variant.h"
#include "vardb.h"
#include "gstore.h"
#include "meta.h"

#include <iostream>

extern GStore * GP;


void BCF::set_types()
{
  
  // Note -- currently no support for FORMAT FLAG items in this implementation of BCF

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

  //   enum bcf_meta_t {
  //     BCF_undef = 0 , 
  //     BCF_uint8 , 
  //     BCF_uint16 , 
  //     BCF_int32 , 
  //     BCF_uint32 , 
  //     BCF_uint64 ,
  //     BCF_float , 
  //     BCF_double , 
  //     BCF_flag , 
  //     BCF_char ,
  //     BCF_string };   
  
  
  bcftype.clear();
  
  // Encode all known meta-fields in this format:
  
  std::vector<std::string> keys = MetaInformation<GenMeta>::field_names();
  for (int k = 0 ; k < keys.size(); k++)
    {
      
      if ( ! MetaMeta::display( keys[k] ) ) continue;

      // BCF   Means                  Internal
      // ------------------------------------------
      // 2     2                      2
      // 1     1                      1
      // ?     ?                      0 (flag)
      // 0     variable-length        -1 (v-length)
      // -1    # alt alleles          -1
      // -2    # alleles              -1
      // -3    # genotypes            -1

      meta_index_t midx = MetaInformation<GenMeta>::field( keys[k] );
      int len = midx.len;
      if ( len ==  0 )      len == 1;  // flag is always 1 bit (0/1)
      else if ( len == -1 ) len == 0;  // 0 is BCF variable length code
      
      if      ( midx.mt == META_INT   ) bcftype[ keys[k] ] = bcf_meta_t( BCF_int32  , len );	
      else if ( midx.mt == META_FLOAT ) bcftype[ keys[k] ] = bcf_meta_t( BCF_double , len );
      else if ( midx.mt == META_TEXT  ) bcftype[ keys[k] ] = bcf_meta_t( BCF_string , len );

      //      else if ( midx.mt == META_FLAG  ) bcftype[ keys[k] ] = bcf_meta_t( BCF_flag   , 1   );
      
    }

  
  // Over-ride with some special, default hard-coded types
  
  bcftype[ "DP" ] = bcf_meta_t( BCF_uint16    , 1 );
  bcftype[ "GL" ] = bcf_meta_t( BCF_float     , -3 );  // -3 is # genotypes
  bcftype[ "GT" ] = bcf_meta_t( BCF_genotype  , 1 );
  bcftype[ "GQ" ] = bcf_meta_t( BCF_uint8     , 1 );
  bcftype[ "HQ" ] = bcf_meta_t( BCF_uint8     , 2 );
  bcftype[ "PL" ] = bcf_meta_t( BCF_uint8     , -3 );
  bcftype[ "SP" ] = bcf_meta_t( BCF_uint8     , 1 );
  bcftype[ "AD" ] = bcf_meta_t( BCF_uint16    , -2 ); // per allele
  bcftype[ "EC" ] = bcf_meta_t( BCF_float     , -1 ); // per alt-allele
  
  // Only the FORMAT fields are encoded here.

}


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
      VCFReader vcf( f , f->tag() , vardb , NULL );
      
      // Get file-ID (would have been created from VARDB)
      file_id = vcf.group_id();
      
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

      // 2 == BCF mode
      vardb->store_bcf_n( file_id , filename , 2 , n );

      plog << "added " << hdr.sample_names.size() << " individuals from BCF " << filename << "\n";
      
      // Meta-fields
      vcf.insert_meta( "##format=BCF4.0" );
      
      for (int i=0; i<hdr.meta_text.size(); i++)
	{
	  vcf.insert_meta( hdr.meta_text[i] );
	}
    }
  
  // refresh bcftypes()
  set_types();
  
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
  int nallele = alt.size() + 1;
  int ngen = (int) (nallele * (nallele+1) * 0.5);
  
  Variant var;
  var.chromosome( Helper::chrCode( hdr.seq_names[ seq_id ] ) );
  var.position( bp );
  var.stop( bp + mstr[1].size() - 1 );
  var.name( mstr[0] == "" ? "." : mstr[0] );
  

  // ** TODO ** change this to calculate a single # of bytes and read/skip those

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
      
      // Look-up BCF type
      
      std::map<std::string,bcf_meta_t>::iterator ii = bcftype.find( format[t] );

      if ( ii == bcftype.end() ) 
	Helper::halt( "could not deal with meta-type " + format[t] );

      int nelem = ii->second.len;
      bool vlen = nelem == 0;
      if ( nelem == -1 ) nelem = nallele - 1;
      else if ( nelem == -2 ) nelem = nallele;
      else if ( nelem == -3 ) nelem = ngen;
      
      // Scan through only
      if ( ii->second.type == BCF_uint16 ) 
	{	  
	  uint16_t dp; 
	  for ( int per = 0 ; per < n ; per++ )
	    {
	      if ( vlen ) read(nelem);
	      for (int g=0; g < nelem; g++) 
		{
		  if ( ! read(dp) ) return false;		  
		}
	    }
	}
//       else if ( ii->second.type == BCF_flag )
// 	{
// 	  // use uint8_t; always 1 element per genotype exactly
// 	  uint8_t x; 
// 	  for ( int per = 0 ; per < n ; per++ )
// 	    if ( ! read(x) ) return false;	  
// 	}      
      else if ( ii->second.type == BCF_float )
	{
	  float gl;

	  for ( int per = 0 ; per < n ; per++ )
	    {
	      if ( vlen ) read(nelem);
	      for (int g=0; g<nelem; g++) 
		if ( ! read( gl ) ) return false;
	    }
	}	 	
      else if ( ii->second.type == BCF_uint8 )
	{
	  uint8_t x; 
	  for ( int per = 0 ; per < n ; per++ )
	    {
	      if ( vlen ) read(nelem);
	      for (int g=0; g < nelem; g++) 
		{
		  if ( ! read(x) ) return false;
		}
	    }
	}
      else if ( ii->second.type == BCF_int32 )
	{
	  int32_t x; 
	  for ( int per = 0 ; per < n ; per++ )
	    {
	      if ( vlen ) read(nelem);
	      for (int g=0; g < nelem; g++) 
		if ( ! read(x) ) return false;
	    }
	}
      else if ( ii->second.type == BCF_double ) 
	{
	  double x; 
	  for ( int per = 0 ; per < n ; per++ )
	    {
	      if ( vlen ) read(nelem);
	      for (int g=0; g < nelem; g++) 
		if ( ! read(x) ) return false;
	    }
	}
      else if ( ii->second.type == BCF_string )
	{
	  std::string x;
	  for ( int per = 0 ; per < n ; per++ )
	    {
	      if ( vlen ) read(nelem);
	      for (int g=0; g < nelem; g++) 
		if ( ! read(x) ) return false;	  
	    }
	}
      else if ( ii->second.type == BCF_genotype )
	{
	  uint8_t x; 
	  for ( int per = 0 ; per < n ; per++ )
	    {
	      if ( vlen ) read(nelem);
	      for (int g=0; g < nelem; g++) 
		{
		  if ( ! read(x) ) return false;
		}
	    }
	}      

    } // next FORMAT-specified tag
  
  // Add offset to DB
  vardb->insert_bcf_index( file_id , var , offset );

  return true;
  
}

bool BCF::read_record( Variant & var , SampleVariant & svar , SampleVariant & svar_g )
{

  // Here we are reading genotpe data for a particular Variant 'var'
  // Variant meta-information goes into 'svar'
  // Genotypes go into 'svar_g' 

  // 'svar' and 'svar_g' may or may not be the same, depending on the structure of
  // the project.


  int32_t                    seq_id;     // 'chromosome' code  
  int32_t                    bp;         // position (BP1)  
  float                      qual;       // quality score
  std::vector<std::string>   mstr;       // ID REF ALT FILTER INFO FORMAT (null padded)
  

  if ( ! read( seq_id )   ) return false;
  if ( ! read( bp )       ) return false;
  if ( ! read( qual )     ) return false;
  if ( ! read( mstr )     ) return false;  
  if (   mstr.size() < 6  ) return false;
  

  // check that Variant spec (chr/bp) from VARDB matches what is in
  // the BCF. Although the chr-coding will be different, so let's just 
  // check base-position
  
  if ( bp != var.position() )
    plog.warn( "mismatching physical position between VARDB and BCF" , 
	       Helper::int2str(bp) + " vs " + var.coordinate() );
  

  // Read genotype; # alleles based on ALT tag, comma-sep values
  
  std::vector<std::string> alt = Helper::char_split( mstr[2] , ',' );
  int nallele = alt.size() + 1;
  int ngen = (int) (nallele * (nallele+1) * 0.5);
  

  // Populate Sample Variant -- variant-level information

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
  

  // Genotype fields to expect based on 
  // use set_format() function in vcf.cpp
  
  std::vector<std::string> format = Helper::char_split( mstr[5] , ':' );


  // Indicate that this SampleVariant contains a BCF-derived buffer
  // Fill buffer as direct raw byte sequence from BCF, for later parsing

  // set 'target', so that VMETA is also appropriately not handled as a BLOB

  // BCF::n (sample size) will have been populated already (from VARDB)
  //                      *but*, we might not want to read it all in... (Masks...)
  
  // but genotype specific info stays w/ genotypes


  svar.set_pointer_to_bcf( this );
    
  svar_g.set_for_bcf_genotypes( n , mstr[5] );



  // Start reading

  int32_t buf_sz = 0;
  
  for ( int t = 0 ; t < format.size(); t++ )
    {
      
      int p = buf_sz;
      
      std::map<std::string,bcf_meta_t>::iterator ii = bcftype.find( format[t] );
      
      if ( ii == bcftype.end() ) 
	Helper::halt( "could not deal with meta-type " + format[t] );
      
      int nelem = ii->second.len;
      bool vlen = nelem == 0;
      if ( nelem == -1 ) nelem = nallele - 1;
      else if ( nelem == -2 ) nelem = nallele;
      else if ( nelem == -3 ) nelem = ngen;
      
      if ( ii->second.type == BCF_uint16 )
	{
	  buf_sz += n * nelem * sizeof(uint16_t);
	  svar_g.bcf_genotype_buf_resize( buf_sz );
	  read(  svar_g.bcf_pointer(p) , n * nelem * sizeof(uint16_t) );
	}
      else if ( ii->second.type == BCF_float )
	{
	  buf_sz += n * nelem * sizeof(float);
	  svar_g.bcf_genotype_buf_resize( buf_sz );
	  read(  svar_g.bcf_pointer(p) , n * nelem * sizeof(float) );
	}
      else if ( ii->second.type == BCF_uint8 )
	{
	  buf_sz += n * nelem * sizeof(uint8_t);
	  svar_g.bcf_genotype_buf_resize( buf_sz );
	  read(  svar_g.bcf_pointer(p) , n * nelem * sizeof(uint8_t) );	  
	}
      else if ( ii->second.type == BCF_double )
	{
	  buf_sz += n * nelem * sizeof(double);
	  svar_g.bcf_genotype_buf_resize( buf_sz );
	  read(  svar_g.bcf_pointer(p) , n * nelem * sizeof(double) );	  
	}
      else if ( ii->second.type == BCF_int32 )
	{
	  buf_sz += n * nelem * sizeof(int32_t);
	  svar_g.bcf_genotype_buf_resize( buf_sz );
	  read( svar_g.bcf_pointer(p) , n * nelem * sizeof(int32_t) );	  
	}
      else if ( ii->second.type == BCF_string ) 
	{
	  for (int i=0;i<n;i++)
	    {
	      int32_t len;
	      read( len );
	      buf_sz += nelem * sizeof(int32_t) + len;
	      
	      svar_g.bcf_genotype_buf_resize( buf_sz );
	      
	      // swap length in
	      uint8_t * d = (uint8_t*)len;

	      *svar_g.bcf_pointer(p)   = d[0];
	      *svar_g.bcf_pointer(p+1) = d[1];
	      *svar_g.bcf_pointer(p+2) = d[2];
	      *svar_g.bcf_pointer(p+3) = d[3];
	      
	      // read text
	      read( svar_g.bcf_pointer(p+4) , len );
	    }
	  
	}
      else if ( ii->second.type == BCF_genotype )
	{
	  
	  // NOTE: currently, assume genotype is uint8_t
	  
	  buf_sz += n * nelem * sizeof(uint8_t);
	  svar_g.bcf_genotype_buf_resize( buf_sz );
	  read( svar_g.bcf_pointer(p) , n * nelem * sizeof(uint8_t) );	  
	}

    }
  
  // now we've finished building the BCF genotype buffer; this may (or
  // may not) get expanded later, as needed given the Mask

  return true;

}



bool BCF::read_record( Variant & var , SampleVariant & svar , SampleVariant & svar_g, int64_t offset )
{
  seek( offset );
  return read_record( var , svar , svar_g );
}


bool BCF::create_header()
{
  
  writing();
  open();
  hdr.clear();

  // Set # of individuals

  n = GP->indmap.size();

  //
  // For now, assume we'll encounter only 'normal' chromosome/sequence
  // codes
  //

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

  
  // Sample IDs
  
  const int n = GP->indmap.size();
  hdr.sample_names.clear(); 
  for ( int i=0; i<n; i++) 
    hdr.sample_names.push_back( GP->indmap(i)->id() );


  // Meta-fields from INFO, FORMAT and FILTER

  std::string mtext = "##source\n" 
    + MetaInformation<VarMeta>::headers( ) 
    + MetaInformation<GenMeta>::headers( META_GROUP_GEN ) 
    + MetaInformation<VarFilterMeta>::headers( META_GROUP_FILTER );

  // F means no empty fields (i.e. last row otherwise)
  // but that should not be a problem in any case
  hdr.meta_text = Helper::char_split( mtext , '\n' , false );  
  
  // Write header info to BCF
  write_header();
  
  set_types();
  
  plog << "inserted header into BCF, " << n << " individuals\n";
  
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
  std::stringstream ss;
  ss << var.consensus.meta;
  m.push_back( ss.str() );
    
  // Get genotype format field; genotype always goes first (this will be *slow*) 
  // look for speed-ups downstream

  std::string format_string = "GT";
  std::set<std::string> allkeys;
  for (int i = 0 ; i < var.size(); i++)
    {
      std::vector<std::string> keys = var(i).meta.keys();
      for (unsigned int j=0; j<keys.size(); j++) allkeys.insert( keys[j] );
    }
  
  std::set<std::string>::iterator i = allkeys.begin();
  while ( i != allkeys.end() )
    {
      format_string += ":" + *i;
      ++i;
    }
  
  m.push_back( format_string );
  write( m );

  // Read genotype; # alleles based on ALT tag, comma-sep values
  // BCF::n should already be set 

  std::vector<std::string> alt = Helper::char_split( var.alternate() , ',' );
  int nallele = alt.size() + 1;
  int ngen = (int) (nallele * (nallele+1) * 0.5);
  
  // Add in genotype information, from consensus
  
  std::vector<std::string> format = Helper::char_split( format_string , ':' );
  
  
  // Default type is  bcf_meta_t( BCF_string , 1 ) 
  // is which implies storage as a single string
  
  for ( int t = 0 ; t < format.size(); t++ )
    {

      const std::string & tag = format[t];
      
      std::map<std::string,bcf_meta_t>::iterator ii = bcftype.find( format[t] );
      
      if ( ii == bcftype.end() ) 
	Helper::halt( "could not deal with meta-type " + format[t] );

      int nelem = ii->second.len;
      if ( nelem == -1 ) nelem = nallele - 1;
      else if ( nelem == -2 ) nelem = nallele;
      else if ( nelem == -3 ) nelem = ngen;

      
      if ( ii->second.type == BCF_uint16 )
	{
	  if ( nelem == 1 ) 
	    for (int i=0; i<n; i++) 
	      write( (uint16_t) var(i).meta.get1_int( tag ) );
	  else
	    {
	      for (int i=0; i<n; i++) 
		{
		  std::vector<int> tmp = var(i).meta.get_int( tag );
		  set_size<int>( nelem , tmp , nallele , ngen );
		  for (int j=0; j<tmp.size(); j++)		    
		    write( (uint16_t) tmp[j] );
		}
	    }
	}

      else if ( ii->second.type == BCF_int32 )
	{
	  if ( nelem == 1 ) 
	    for (int i=0; i<n; i++) 
	      write( (int32_t) var(i).meta.get1_int( tag ) );
	  else
	    {
	      for (int i=0; i<n; i++) 
		{
		  std::vector<int> tmp = var(i).meta.get_int( tag );
		  set_size<int>( nelem , tmp , nallele , ngen );
		  for (int j=0; j<tmp.size(); j++)
		    write( (int32_t) tmp[j] );
		}
	    }
	}

      else if ( ii->second.type == BCF_float )
	{
	  if ( nelem == 1 ) 
	    for (int i=0; i<n; i++) 
	      write( (float) var(i).meta.get1_double( tag ) );
	  else
	    {
	      for (int i=0; i<n; i++) 
		{
		  std::vector<double> tmp = var(i).meta.get_double( tag );
		  set_size<double>( nelem , tmp , nallele , ngen );
		  for (int j=0; j<tmp.size(); j++)
		    write( (float) tmp[j] );
		}
	    }
	}
      
      else if ( ii->second.type == BCF_double )
	{
	  if ( nelem == 1 ) 
	    for (int i=0; i<n; i++) 
	      write( (double) var(i).meta.get1_double( tag ) );
	  else
	    {
	      for (int i=0; i<n; i++) 
		{
		  std::vector<double> tmp = var(i).meta.get_double( tag );
		  set_size<double>( nelem , tmp , nallele , ngen );
		  for (int j=0; j<tmp.size(); j++)
		    write( (double) tmp[j] );
		}
	    }
	}

      else if ( ii->second.type == BCF_uint8 )
	{
	  if ( nelem == 1 ) 
	    for (int i=0; i<n; i++) 
	      write( (uint8_t) var(i).meta.get1_int( tag ) );
	  else
	    {
	      for (int i=0; i<n; i++) 
		{
		  std::vector<int> tmp = var(i).meta.get_int( tag );
		  set_size<int>( nelem , tmp , nallele , ngen );
		  for (int j=0; j<tmp.size(); j++)
		    {
		      // worry about clipping here (i.e. cap at 255 for PL, etc)?
		      if ( tmp[j] < 0 ) tmp[j] == 0;
		      else if ( tmp[j] > 255 ) tmp[j] == 255;
		      write( (uint8_t) tmp[j] );
		    }
		}
	    }
	}
      
      else if ( ii->second.type == BCF_genotype )
	{
	  for (int i=0; i<n; i++) 
	    write( var(i).bcf() );
	}

      else if ( ii->second.type == BCF_string ) 
	{	  
	  // Variable
	  std::vector<std::string> m(n);
	  for (int i=0;i<n;i++)
	    {
	      // store values as comma-delimited list (per individual)
	      m[i] = var(i).meta.as_string( tag , "," );
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
