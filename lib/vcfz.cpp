
#include "vcfz.h"
#include "variant.h"
#include "vardb.h"
#include "gstore.h"
#include "meta.h"

#include <iostream>

extern GStore * GP;

std::string VCFZ::current_format = "";
int VCFZ::gt_field = 0;
std::vector<meta_index_t*> VCFZ::formats;

bool VCFZ::open()
{
  
  if ( file ) close();
  
  if ( readmode ) file = bgzf_open( filename.c_str(), "r" );

  else 
    { 
      Helper::halt( "writing BGZF-VCFs not yet supported" );
      file = bgzf_open( filename.c_str(), "w" );
    }

  return file == NULL;
}


void VCFZ::close() 
{
  if ( file ) bgzf_close( file );
  file = NULL;
}



bool VCFZ::read_line( std::vector<char> * line )
{        
    line->clear();
    while ( 1 ) 
      {
	const int c = bgzf_getc( file );
	if ( c < 0 ) return false;	
	if ( c == '\n' ) 
	  {
	    line->push_back( '\0' );
	    return true;
	  }
	line->push_back( c );
      }
}




void VCFZ::read_header()
{
  if ( ! vardb ) Helper::halt( "no VARDB attached to VCFZ class" );
  
  
  // the header should already have been added when doing the index-vcf
  //  (but implement here: )

//   hdr.clear();
//   hdr.meta_text = Helper::char_split( mtext , '\n' , false ); 
  

  // Assign to VarDB
  
  if ( vardb )
    {

      // and meta-information. 
      // we'll never want to look at this header
      // again when accessing genotype/variant data      
      
      // Use VCFReader class functions to process meta-information
      
      // VCFZ will be in file-index
      
      File * f = GP->fIndex.file( filename );
      if ( f == NULL ) Helper::halt( "internal error in VCFZ class, parsing header"  );
      
      // NULL means no SEQDB attached for now, i.e no REF checking
      VCFReader vcf( f , f->tag() , vardb , NULL );
      
      // Get file-ID (would have been created from VARDB)
      file_id = vcf.group_id();
      
      int nind = 0;

      //       // Individuals
//       for (int i=0; i<n; i++)
// 	{
// 	  Individual ind( hdr.sample_names[i] );
// 	  vardb->insert( file_id , ind );
// 	}


//       // Register BCF file in VARDB; also the # of individuals in the
//       // BCF (i.e. so we do not need to parse the header when
//       // accessing variant data from the BCF in future

//       // Meta-fields
//       vcf.insert_meta( "##format=VCF4.1" );
      
//       for (int i=0; i<hdr.meta_text.size(); i++)
// 	{
// 	  vcf.insert_meta( hdr.meta_text[i] );
// 	}
//     }
  

      // 1 means BGZF-compressed VCF (not BCF)

      //   (only needed to read a BCF)
      
      vardb->store_bcf_n( file_id , filename , 1 , nind );
      
      
      
    }
}


bool VCFZ::index_record( )
{

  // requires an attached VARDB
  if ( ! vardb ) return false;

  // track offset for this variant
  int64_t offset = bgzf_tell( file );

  // error or end of file?
  if ( offset < 0 ) return false;
  
  // read up to the EOL (no need to process the line at all)
  std::vector<char> line;
  if ( ! read_line( &line ) ) return false;

  // skip headers/comments, etc
  if ( line[0] == '#' ) return true;

  Variant var; 
 
  // Need to read chr and bp and put into the Variant
  // rather than parse the whole line, just use a simple search for the first two tabs
  
  int tab1 = 0;
  while ( tab1 < line.size() ) { if ( line[tab1] == '\t' ) break; else ++tab1; }
  int tab2 = tab1 + 1;
  while ( tab2 < line.size() ) { if ( line[tab2] == '\t' ) break; else ++tab2; }
  
  // 012345679
  // xx\xxxx\x
  // tab1 = 2
  // tab2 = 7
  
  if ( tab1 == 0 || tab2 - tab1 == 1 ) Helper::halt( "problem with VCF chr/bp fields" );
  
  std::string c1( line.begin() , line.begin() + tab1 ) ;
  std::string c2( line.begin() + tab1 + 1 , line.begin() + tab2 ) ;

  int bp;
  if ( ! Helper::str2int( c2 , bp ) ) Helper::halt( "could not parse POS field" );

  var.chromosome( Helper::chrCode( c1 ) );
  var.position( bp );

  // Add offset to DB (can use 'BCF' function, does same thing; file_id
  // will identify this as a VCF downstram
  
  vardb->insert_bcf_index( file_id , var , offset );
  
  return true;
  
}



bool VCFZ::read_record( Variant & var , SampleVariant & svar , SampleVariant & svar_g, int64_t offset )
{
  bgzf_seek( file, offset, SEEK_SET);
  return read_record( var , svar , svar_g );
}



bool VCFZ::read_record( Variant & var , SampleVariant & svar , SampleVariant & svar_g )
{

  // Here we are reading genotpe data for a particular Variant 'var'
  // In practice, the wrapper function above will get called, which first moves to 
  // the correct place in the VCF

  // We will receive a Variant spec with chr/bp from the VARDB. We should check that it matches
  // what we see here in the file
  
  // As when parsing a BCF: Variant meta-information goes into 'svar'
  //                        Genotypes go into 'svar_g' 
  
  // 'svar' and 'svar_g' may or may not be the same, depending on the structure of
  // the project.

  //
  // Extract Helper::char_tok from BGZF'ed VCF
  //

  std::vector<char> line;
  read_line( &line );

  int ntok;
  Helper::char_tok tok( &(line[0]) , line.size() , &ntok , '\t' );
  if ( ntok < 8 ) Helper::halt( "invalid VCF entry, less than 8 fields" );
  
  // Sanity check: does BP from VARDB match BP from VCF?
  


  int bp;
  if ( ! Helper::str2int( tok(1) , bp ) ) Helper::halt( "trouble processing POS field in VCF" );
  if ( bp != var.position() ) Helper::halt( "index out of sync with VCF" );
  
  // Populate basic VCF fields: 

  // 1) variants 
  
  // Attempt to populate as much information as possible.
  
  var.name( tok(2) );
  
  // Sample-level information
  
  svar.reference( tok(3) );
  svar.alternate( tok(4) );
  
  double qual;
  if ( Helper::str2dbl( tok(5) , qual ) ) 
    svar.quality( qual ); 
  else 
    svar.quality(-1); 
  
  
  // Filter and meta-information; here we supply the info to be able
  // to put the meta-header info into the VARDB if it hasn't already
  // been declared...
  
  svar.filter( tok(6) , vardb , file_id );
  svar.info( tok(7) , vardb , file_id );
  
  // Let's not bother checking and rechecking the validity of the record; this will
  // be better done with a separate function (or just advocating use of a VCF validator 
  // tool that checks conformance to the spec) 

  var.valid( true );
  
  // If no genotypes, we can leave now

  if ( ntok <= 8 ) return true;
  

  //
  // Set genotype format, using VCFReader to parse the FORMAT field (has cache)
  //
  
  
  set_format( tok(8) );
  

  // Add genotypes to genotype-holding sample-variant; they will be
  // parsed at a later date, if needed; number of genotypes should be
  // checked in svar.cpp when expanding
  
  svar_g.set_vcfz_buffer( tok , gt_field , &formats );

    
  // Return this variant (and any called genotypes)

  return true;
}





bool VCFZ::write_header()
{
  
  std::string mtext = "##source=pseq\n" 
    + MetaInformation<VarMeta>::headers( ) 
    + MetaInformation<GenMeta>::headers( META_GROUP_GEN ) 
    + MetaInformation<VarFilterMeta>::headers( META_GROUP_FILTER );
  
  //  write( mtext );
  return true;
}
 

bool VCFZ::write_record( const Variant & var )
{    
  if ( ! file ) return false;  
  return false;
}


bool VCFZ::set_format( const std::string & f )
{

  // this function basically duplicates the function in VCFReader
  
  if ( f == current_format ) return false;  
  
  // check that all formats are unique
  std::set<std::string> fset; 
  
  std::vector<std::string> tok = Helper::char_split(f,':');
  current_format = f;
  formats.resize( tok.size() , NULL );
  
  gt_field = -1;
  for ( int i = 0 ; i < tok.size() ; i++ )
    {
      
      fset.insert( tok[i] );
      
      // genotype is not treated as meta-information
      // or if no value set in value  
      
      if ( tok[i] == "GT" ) 
      {
        gt_field = i; 
	continue;
      }

      if ( tok[i] == "." ) continue;
      
      // this will return a NULL if the field is not recognized
      // (i.e. does not exist, or has been skipped over and so not
      // registered by this class) NULLs imply the the subsequent
      // set() function will ignore these, as expected
      
      formats[i] = MetaInformation<GenMeta>::index( tok[i] );
      
      if ( ! formats[i] ) 
	{
	  // add to internal VARDB header 
	  MetaInformation<GenMeta>::field( tok[i] , META_TEXT , 1 , "undeclared genotype tag" );
	  vardb->insert_metatype( file_id , tok[i] , META_TEXT , 1 , META_GROUP_GEN , "undeclared genotype tag" );
	  plog.warn( "undefined FORMAT field:", tok[i] );
	}

    }
  
  if ( fset.size() != tok.size() ) 
    Helper::halt( "problem in VCF FORMAT field: repeated tags: " + f ); 
  
  if ( gt_field == -1 ) 
    Helper::halt( "no GT field specified in VCF genotype" );

  return true; // indicates a change was made
}
