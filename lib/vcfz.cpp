
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
  file = bgzf_open( filename.c_str(), readmode ? "r" : "w" );    
  if ( ! file ) Helper::halt( "could not open " + filename ) ;
  return file == NULL;
}


void VCFZ::close() 
{
  if ( file ) { bgzf_close( file ); }
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




void VCFZ::read_header( Mask & mask )
{
  
  if ( ! vardb ) Helper::halt( "no VARDB attached to VCFZ class" );  
  
  // When creating the index for a VCFZ, we first need to read the
  // header and store the INFO in the VARDB, same as if we had used
  // 'load-vcf'.  This function is only called when performing an
  // 'index-vcf' call

  // This is largely the same code as in vcfiterate.cpp, where we use VCFReader to process
  // the VCF header. Because the BGZF can be read like a normal gzipped VCF
  
  // Load, parse VCF file; store variant and genotype information, and
  // meta-information, in vardb
  
  File vcffile( filename , VCF );
  
  VCFReader v( &vcffile , "" , vardb , NULL );  
  
  // Work through VCF

  vardb->begin();

  while ( 1 ) 
    { 
      
      VCFReader::line_t l = v.parseLine();
      
      if ( l == VCFReader::VCF_EOF ) break;
      
      if ( l == VCFReader::VCF_INVALID ) 
	{
	  continue;
	}
      
      if ( l == VCFReader::VCF_HEADER )
	{
	  break;
	}           
    }

  // Wrap up  
  vardb->commit();
  

//   //
//   // In FIX/XY mode when in single-VCF mode, we need to populate the Sex codes in the indmap
//   //
  
//   for (int i = 0; i< GP->indmap.size() ; i++)
//     {
//       sType s = GP->inddb.sex( GP->indmap.ind(i)->id() );
//       GP->indmap.ind(i)->sex( s );
//     }
    

  //
  // Finally, also store in 'bcfs' VARDB table, as type1 (BGZF VCF)
  //

  // Store file-ID that was assigned in VARBD

  file_id = v.group_id();

  
  // 1 == VCF (BGZF'ed)
  // 0 == no need to set N individuals (BCF only)
  
  vardb->store_bcf_n( file_id , filename , 1 , 0 );

  // All done reading header, so can close now
  // File should be closed when VCFReader goes out of scope just below


          
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
 
  // Need to read chr and bp and ID and put into the Variant
  // rather than parse the whole line, just use a simple search for the first three tabs
  
  int tab1 = 0;
  while ( tab1 < line.size() ) { if ( line[tab1] == '\t' ) break; else ++tab1; }
  int tab2 = tab1 + 1;
  while ( tab2 < line.size() ) { if ( line[tab2] == '\t' ) break; else ++tab2; }

  if ( tab1 == 0 || tab2 - tab1 == 1 ) Helper::halt( "problem with VCF chr/bp fields" );
  
  int tab3 = tab2 + 1;
  while ( tab3 < line.size() ) { if ( line[tab3] == '\t' ) break; else ++tab3; }

  std::string c1( line.begin() , line.begin() + tab1 ) ;
  std::string c2( line.begin() + tab1 + 1 , line.begin() + tab2 ) ;
  std::string c3( line.begin() + tab2 + 1 , line.begin() + tab3 ) ;

  int bp;
  if ( ! Helper::str2int( c2 , bp ) ) Helper::halt( "could not parse POS field" );

  var.chromosome( Helper::chrCode( c1 ) );
  var.name( c3 );
  var.position( bp );

  // Add offset to DB (can use 'BCF' function, does same thing; file_id
  // will identify this as a VCF downstram
  
  vardb->insert_bcf_index( file_id , var , offset );
  
  return true;
  
}



bool VCFZ::read_record( Variant & var , SampleVariant & source , SampleVariant & svar , SampleVariant & svar_g, int64_t offset )
{
  bgzf_seek( file, offset, SEEK_SET);
  return read_record( var , source , svar , svar_g );
}



bool VCFZ::read_record( Variant & var , SampleVariant & source , SampleVariant & svar , SampleVariant & svar_g )
{
    
  // Here we are reading genotpe data for a particular Variant 'var'
  // In practice, the wrapper function above will get called, which first moves to 
  // the correct place in the VCF

  // We will receive a Variant spec with chr/bp from the VARDB. We should check that it matches
  // what we see here in the file
  

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
  svar.info( tok(7) , vardb , file_id , &var );
  
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
  
  // NOTE: here we set svar, not svar_g, which is the 'target' for genotypes.
  // Under a single sample, then svar == svar_g == consensus, so no problem.
  // Under any multiple-sample case, we want the VCF-buffers to stay with each Svar at 
  // this point, however, similar to 
  
  source.set_vcfz_buffer( tok , gt_field , &formats );
  
  // Return this variant (and any called genotypes)

  return true;
}





bool VCFZ::write_header()
{
  
  std::string mtext = "##fileformat=" + PLINKSeq::CURRENT_VCF_VERSION() + "\n"
    "##source=pseq\n"    
    + MetaInformation<VarMeta>::headers( ) 
    + MetaInformation<GenMeta>::headers( META_GROUP_GEN ) 
    + MetaInformation<VarFilterMeta>::headers( META_GROUP_FILTER );
  bgzf_write( file , &(mtext[0]) , mtext.size() );

  // individuals
  const int n = GP->indmap.size();
  std::stringstream ss;
  ss << "#CHROM\tPOS\tID\tREF\tALT\tQUAL\tFILTER\tINFO\tFORMAT";  
  for ( int i=0; i<n; i++) ss << "\t" << GP->indmap(i)->id();
  ss << "\n";
  std::string indhdr = ss.str();
  bgzf_write( file , &(indhdr[0]) , indhdr.size() ); 
  
  return true;
}
 

bool VCFZ::write_record( const Variant & var )
{    
  if ( ! file ) return false;  
  std::string vcfline = var.VCF();
  bgzf_write( file , &(vcfline[0]) , vcfline.size() );
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
