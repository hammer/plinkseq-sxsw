
#include "vcf.h"
#include "gstore.h"
#include "genotype.h"
#include "defs.h"

#include <set>
#include <string>

using namespace std;
using namespace Helper;

extern GStore * GP;


void VCFReader::set_seqdb( SeqDBase * s )
{
  seqdb = GP->seqdb.attached() ? s : NULL;
}


void VCFReader::get_meta( const std::set< std::string> & s )
{
  explicit_meta = true;
  std::set<std::string>::const_iterator i = s.begin();
  while ( i != s.end() )
    {
      meta_want.insert( *i );
      ++i;
    }
}


void VCFReader::ignore_meta( const std::set< std::string> & s )
{
  std::set<std::string>::const_iterator i = s.begin();
  while ( i != s.end() )
    {
      meta_ignore.insert( *i );
      ++i;
    }
}


//////////////////////////////////////////////////////////////////////
//  VCF definitions

// Only first 8 fields fixed

// #CHROM
// POS
// ID
// REF
// ALT
// QUAL
// FILTER    
// INFO

// If genotype data, then also:
// FORMAT
// GENOTYPES...


VCFReader::line_t VCFReader::parseLine( Variant ** pvar )
{


  if ( return_var ) *pvar = NULL;

  if ( from_stdin && cin.eof() ) return VCF_EOF;
  if ( (!from_stdin) && file.eof() ) return VCF_EOF;

  //
  // Read next line from file, (or STDIN)
  // 

  std::string s;
  if ( from_stdin ) std::getline( std::cin , s );
  else s = file.readLine();

     
  //
  // If not a valid line, do nothing
  //

  if ( s == "" ) return VCF_EOF;
  
  if ( s.size() < 3 ) 
    {
      plog.warn( "invalid line with fewer than 3 characters in VCF: " + s  );
      return VCF_INVALID;
    }
  
  //
  // Either meta-information, header row, or variant row
  //
  
  if ( s[0] == '#' )
    {
      if ( s[1] == '#' )
	{
	  getMetaInformation(s);
	  return VCF_META;
	}
      else
	{
	  getHeader(s);
	  summary();
	  return VCF_HEADER;
	}
    }

  //
  // Get variant from VCF row, and insert into the database
  //  
  
  if ( return_var ) 
    {
      *pvar = new Variant( getVariant(s) );      
    }
  else
    {
      Variant v = getVariant(s);            
      
      if ( v.valid() )
	{
	  ++vcnt;
	  if ( refdb ) refdb->insert( file_id , v );
	  else vardb->insert_consensus( file_id , v );      
	}
//       else
// 	plog.warn( "read invalid line from VCF" , s );
    }

  return VCF_VARIANT;

}



void VCFReader::summary() const
{
  // return_var implies we are directly processing the VCF, so no output here
  if ( ! GP->single_file_mode() )
    plog << "loading : " << fileinfo->name() 
	 << " ( " << icnt << " individuals )\n";
}


void VCFReader::getMetaInformation(const std::string & s)
{

  // assume key=value pairs

  std::vector<std::string> tokens = Helper::char_split( s.substr(2) , '=' );

  std::vector<std::string> tok(2);
  
  std::vector<std::string>::iterator tok_iter = tokens.begin();
  tok[0] = *tok_iter;
  ++tok_iter;

  bool first = true;
  while ( tok_iter != tokens.end() )
    {
      if ( ! first ) tok[1] += "=";
      tok[1] += *tok_iter;
      ++tok_iter;
      first = false;
    }
  
  
  //
  // If a single keyword, just insert
  //
  
  if ( tok.size() == 1 )
    {
      if ( refdb ) refdb->insert_header( file_id , tok[0] , "" );
      else vardb->insert_header(  file_id , tok[0] , "" ); 
      return;
    }
  

  //
  // Check for VCF version as first line of meta-information
  //
      
  if ( tok[0] == "fileformat" || tok[0] == "format" )
    {
      // match on VCFv3.3, vcf3.3, VCF3.3, etc

      // and also BCFv4.0
      if ( tok[1].size() >= 3 && tok[1].substr( tok[1].size() - 3 ) == "3.3" ) 
	version = VCF_3_3;
      else if ( tok[1].size() >= 3 && tok[1].substr( tok[1].size() - 3  ) == "4.0" )
	version = VCF_4_0;
      else if ( tok[1].size() >= 3 && tok[1].substr( tok[1].size() - 3  ) == "4.1" ) 
	version = VCF_4_1;
      else 
	Helper::halt( "could not recognize VCF version " + tok[1] );
    }
  else if ( version == VCF_UNKNOWN ) 
    {
      Helper::halt( "Version number not specified for VCF" );
    }
  
  //
  // Check for keyword meta-field specification
  //
      
  if ( tok[0] == "INFO" || 
       tok[0] == "FORMAT" )
    {
      
      //
      //  VCF v3.3
      // 
      
      //  name,num,type,description
      
      // VCF 4.0
      
      //  ##INFO=<ID=ID,Number=number,Type=type,Description=description>
      // . for Number implies variable length (== -1 from v3.3)
      // 
      
      // The Number entry is an Integer that describes the
      // number of values that can be included with the INFO
      // field. For example, if the INFO field contains a single
      // number, then this value should be 1. However, if the
      // INFO field describes a pair of numbers, then this value
      // should be 2 and so on.
      
      // If the number of possible values varies, is
      // unknown, or is unbounded, then this value should
      // be '.'.
      
      // Possible Types are: Integer, Float, Character, String
      // and Flag.
      
      // The 'Flag' type indicates that the INFO field
      // does not contain a Value entry, and hence the
      // Number should be 0 in this case.
	
      // Genotype meta-information has similar format:
	  
      // ##FORMAT=<ID=ID,Number=number,Type=type,Description=description> 
	  
      std::string name = "";
      int num = -9;
      std::string type = "";
      mType mt = META_UNDEFINED;
      std::string desc;
      
      if ( version == VCF_3_3 )
	{
	  // check not in 4.0 format (with <ID=id, etc)
	  if ( tok[1].substr(0,1) == "<" ) 
	    plog.warn("likely malformed " + tok[0] + " in VCF3.3" , tok[1] );
	  
	  std::vector<std::string> tok2 = Helper::quoted_parse( tok[1] );
	  if ( tok2.size() == 4 )
	    {
	      name = tok2[0];
	      if ( ! str2int( tok2[1] , num ) ) num = -1;
	      type = tok2[2];
	      desc = Helper::unquote( tok2[3] );
	    }
	  else 
	    plog.warn("malformed VCF3.3 " + tok[0] + " line", tok[1] );
	}
      else if ( version == VCF_4_0 || version == VCF_4_1 )
	{

	  if ( tok[1].substr(0,1) != "<" ) 
	    plog.warn("likely malformed " + tok[0] + " in VCF4.0" , tok[1] );

          // remove any <tags> around KEY=<VALUE>
	  
          tok[1] = Helper::remove_tags( tok[1] );
	  std::map<std::string,std::string> l = Helper::quoted_comma_keypair_split( tok[1] );
	  std::map<std::string,std::string>::iterator i = l.begin();
	  bool obs_name = false , obs_type = false;
	  while ( i != l.end() )
	    {
	      
	      // ##FORMAT=<ID=ID,Number=number,Type=type,Description=description> 
	      if ( i->first == "ID" || i->first == "id" ) { obs_name = true; name = i->second; } 
	      else if ( i->first == "Number" || i->first == "number" || i->first == "num" || i->first == "len" ) 
		{
                  // implicitly makes Number=. be parsed as variable length, 
                  // which I think is spec.
		  if ( ! str2int( i->second , num ) ) num = -1;
		}
	      else if ( i->first == "Type" || i->first == "type" ) { obs_type = true; type = i->second; } 
	      else if ( i->first == "Description" || i->first == "desc" || i->first == "description" ) 
		desc = Helper::unquote( i->second );
	      ++i;
	    }
	  if ( ! ( obs_name && obs_type ) )
	    plog.warn("malformed VCF4.0 " + tok[0] + " line", tok[1] );
	}
   

      //
      // Ignore or read?
      //
      
      if ( explicit_meta && meta_want.find( name ) == meta_want.end() ) return;
      if ( meta_ignore.find( name ) != meta_ignore.end() ) return;

      if ( tok[0] == "INFO" ) 
	meta_read_var.insert( name );
      else 
	meta_read_geno.insert( name );


      //
      // Process and check
      //
      
      if ( Helper::is_int( type ) )  mt = META_INT;
      else if ( Helper::is_float( type ) ) mt = META_FLOAT;
      else if ( Helper::is_text( type ) ) mt = META_TEXT;  // text == String || Char
      else if ( Helper::is_flag( type ) ) mt = META_FLAG;
      
      
       // Does this contain valid information?

       if ( name == "" || mt == META_UNDEFINED || num < -1 ) 
	 return;

       // Insert in DB as header line
       
       if ( refdb ) refdb->insert_header( file_id , tok[0] , tok[1] );
       else vardb->insert_header( file_id , tok[0] , tok[1] ); 
       

       if ( name != "GT" ) 
	 {
	   // Insert in DB as meta-field
	   int mgrp =  tok[0] == "INFO" ? META_GROUP_VAR : META_GROUP_GEN;
	   
	   if ( refdb ) 
	     refdb->insert_metatype( file_id , name , mt, num, mgrp, desc );
	   else 
	     vardb->insert_metatype( file_id , name , mt, num, mgrp, desc );
	   
	   // Register with internal meta-field store (for either variant, or genotype)
	   if ( tok[0] == "INFO" )
	     {
	       if ( refdb ) 
		 {
		   MetaInformation<RefMeta>::field( name , mt , num , desc );
		   MetaInformation<VarMeta>::field( name , mt , num , desc );
		 }
	       else
		 MetaInformation<VarMeta>::field( name , mt , num , desc );
	     }
	   else
	     MetaInformation<GenMeta>::field( name , mt , num , desc );			    
	 }
       else 
	 {	    
	   // just store as header line, as we assume no meta-field has been specified
	   if ( refdb ) 
	     refdb->insert_header( file_id , tok[0] , tok[1] ); 
	   else  
	     vardb->insert_header( file_id , tok[0] , tok[1] ); 
	 }	    

       return;
     }


   if ( tok[0] == "FILTER" ) 
     {

       // Insert in VarDB
       if ( refdb ) 
	 refdb->insert_header( file_id , tok[0] , tok[1] ); 
       else 
	 vardb->insert_header( file_id , tok[0] , tok[1] ); 

       // Parse as meta-information 
       // VCF3.3:  For FILTERs, a simple KEY,"DESCRIPTION" format

       std::string name = "";
       std::string desc = "";

       if ( version == VCF_3_3 )
	 {
	   vector<string> tok2 = Helper::quoted_parse( tok[1] );	  
	   if ( tok2.size() == 2 )
	     {
	       name = tok2[0];
	       desc = tok2[1];
	     }
	   else
	     plog.warn("malformed FILTER line in VCF", tok[1] );
	 }
       else if ( version == VCF_4_0 || version == VCF_4_1 )
	 {
	   // VCF4.0:  ##FILTER=<ID=ID,Description=description>

	   tok[1] = Helper::remove_tags( tok[1] );

	   std::map<std::string,std::string> l = Helper::quoted_comma_keypair_split( tok[1] );
	   std::map<std::string,std::string>::iterator i = l.begin();
	   while ( i != l.end() )
	     {
	       if ( i->first == "ID" || i->first == "id" ) name = i->second;
	       else if ( i->first == "Description" 
			 || i->first == "description" 
			 || i->first == "desc" ) desc = i->second;
	       ++i;
	     }
	 }

       if ( name == "" ) 
	 {
	   plog.warn("malformed FILTER line in VCF" , tok[1] );
	   return;
	 }

      //
      // Ignore or read?
      //
      
      if ( explicit_meta && meta_want.find( name ) == meta_want.end() ) return;
      if ( meta_ignore.find( name ) != meta_ignore.end() ) return;
      meta_read_filter.insert( name );

      // 
      // Okay, insert
      //

      if ( refdb ) 
	refdb->insert_metatype( file_id , name , META_FLAG, 1, META_GROUP_FILTER , desc );
      else  
	vardb->insert_metatype( file_id , name , META_FLAG, 1, META_GROUP_FILTER , desc );
      
      MetaInformation<VarFilterMeta>::field( name , META_FLAG , 1 , desc );

      return;
     }


   // Otherwise, if not INFO, FORMAT or FILTER, treat as whole value
   // (not a meta-field specification)

   if ( refdb ) 
     refdb->insert_header( file_id , tok[0] , tok[1] ); 
   else 
     vardb->insert_header( file_id , tok[0] , tok[1] ); 

   return;

 }




bool VCFReader::getHeader( const std::string & s_ )
{
  
  obs_header = true;
  
  // silly fix, but to get rid of trailing tab issue in GATK VCF
  std::string s = s_;
  const int len = s.size();
  if ( s.substr( len - 1) == "\t" ) 
    s = s.substr( 0 , len-1 );
  
    
  // Insert this into the header database in full
  
  if ( refdb ) 
    refdb->insert_header( file_id , s , "" ); 
  else  
    vardb->insert_header( file_id , s , "" ); 
    
  
  // Skip first # character; space delimited
  
  std::vector<std::string> tok = Helper::char_split( s.substr(1) , '\t' );
  std::vector<std::string>::iterator tok_iter = tok.begin();  
  
  std::vector<std::string> head(8);
  head[0] = "CHROM";
  head[1] = "POS";
  head[2] = "ID";
  head[3] = "REF";
  head[4] = "ALT";
  head[5] = "QUAL";
  head[6] = "FILTER";
  head[7] = "INFO";

  //  head += "CHROM","POS","ID","REF","ALT","QUAL","FILTER","INFO";
  
  for (int i=0; i<8; i++)
    {
      
      if ( tok_iter == tok.end() ) 
	Helper::halt("misformed #header row in VCF (fewer than 8 fields)");
      
      if ( Helper::lexical_cast<std::string>(*tok_iter) != head[i] ) 
	plog.warn( "potential problem with header row: " + *tok_iter + ", expecting " + head[i] + "\n" );
      
      ++tok_iter;
    }
  

  //
  // Genotype do not get read into a REFDB
  //

  if ( refdb ) { icnt = 0; return true;  }
  


  //
  // A genotype FORMAT field and 0+ individuals specified? 
  //
  
  if ( tok_iter != tok.end() && *tok_iter == "FORMAT" ) { ++tok_iter; } 
  else return true;
  

  //
  // Individual IDs
  //

  icnt = 0;
  sex.clear();

  while ( tok_iter != tok.end() )
    {
      
      Individual i( *tok_iter );
      
      // Add to the genotype-database
      
      vardb->insert( file_id , i );

      // Working in fix X/Y mode?       

      if ( fixxy_mode )
	{
	  sType s = inddb->sex( *tok_iter ) ;
	  if ( s == UNKNOWN_SEX ) plog.warn("unknown sex for sample", *tok_iter );
	  sex.push_back( inddb->sex( *tok_iter ) );	  
	  
	}  


      // Keep track of how many individuals to expect per line      

      ++icnt;
      
      ++tok_iter;
    }
  
  return true;
}


Variant VCFReader::getVariant(const std::string & s)
{

  if ( ! obs_header ) 
    Helper::halt("missing header in VCF? error at line:\n" + s );

 
  std::string  chr;  
  int          pos;     
  int          stop = 0;
  std::string  id;   
  std::string  ref;  
  std::string  alt;  
  double       qual;
  std::string  filter;
  std::string  info;
  std::string  format;

  // Tab-delimited (white-space) ; but do not send trailing \n 

  int toksize;

  Helper::char_tok tok( s[ s.size()-1] == '\n' ? s.substr( 0 , s.size() - 1 ) : s , &toksize , '\t' );		

  // When reading in from a VCF file, it will, by definition, only
  // correspond to a single sample; therefore, we want to load up
  // into, and write from, the consensus SampleVariant here
  
  Variant var( false );
  
  if ( toksize < 8 ) return var;

 
  // Get primary identifiers: CHR,POS and ID, etc.  If anything fails
  // here, return a variant in invalid state

  
  // requires a valid chr and position to be accepted as a variant
  if ( ! processVCF( tok(0) , &chr   ) ) return var; 


  if ( ! processVCF( tok(1) , &pos   ) ) 
    {
      bool okay;
      std::string s = tok(0);
      s += ":";
      s += tok(1);
      Region test( s , okay );
      if ( ! okay ) return var; 
      pos = test.start.position();
      stop = test.stop.position();      
    }

  // these return false is non-valid value
  processVCF( tok(2) , &id     );
  processVCF( tok(3) , &ref    );
  processVCF( tok(4) , &alt    );
  processVCF( tok(5) , &qual   );
  processVCF( tok(6) , &filter );
  processVCF( tok(7) , &info   );
  
  //
  // If applying a mask, do we want to keep this variant?
  //

  if ( pfilter )
    {

      if ( ! contains( Helper::chrCode( chr ) , pos , pos ) ) 
      {
	  var.valid( false );
	  return var; 
      }
    }
  
  if ( using_idfilter )
    {
      if ( idfilter.find( id ) == idfilter.end() )
	{
	  var.valid( false );
	  return var;
	}
    }


  // Key variant fields

  var.name( id );

  var.chromosome( Helper::chrCode(chr) );
  var.position( pos );

  
  // SampleVariant (consensus) fields

  var.consensus.reference(ref);   
  var.consensus.alternate(alt); 
  var.consensus.quality(qual); 


  // Here we supply the info to be able to 
  // put the meta-header info into the VARDB
  // if it hasn't already been declared...
  
  var.consensus.filter(filter , vardb , file_id );
  var.consensus.info(info , vardb , file_id);

  //
  // End field
  // 
  
  if ( stop == 0 )
    {
      
      // is this a symbollic REF allele (in which case look for an 'END' INFO field
      // i.e. for SVs.
      
      if ( ref[0] == '<' || ref[0] == '.' ) 
	{
	  if ( var.consensus.meta.has_field( PLINKSeq::VCF_END_FIELD() ) )
	    stop = var.consensus.meta.get1_int( PLINKSeq::VCF_END_FIELD() ) ;	      
	}
      else
	stop = pos + ref.size() - 1 ;  
    }
  
  var.stop( stop );

  
    
  //
  // If SEQDB attached, check that the reference is correct
  //

  if ( seqdb ) 
    {

      std::string sref = seqdb->lookup( Helper::chrCode( chr ) , pos , pos + ref.size() - 1 );

      if ( sref != "" )
	{
	  Helper::str2upper( sref );
	  if ( sref != ref )
	    {
	      // ref is 'N' ? 
	      bool mismatch = false;
	      for (int i=0; i<ref.size(); i++)
		if ( sref[i] != ref[i] && sref[i] != 'N' ) { mismatch = true; break; } 
	      if ( mismatch ) 
		plog.warn( "mismatch in reference sequence between VCF and SEQDB", 
			   var.coordinate() + " " + ref + " vs " + sref );	  
	    }
	}
    }
  var.valid(true);
  

  //
  // In VCF -> REFDB mode, we only care about the variant (not genotype) information
  //

  if ( refdb ) { icnt = 0; return var;  }



  //
  // Do we have any genotypes specified here; if so, we need a format
  // that tells us what to expect
  //

  if ( toksize <= 8 ) 
    {
        //
      // If reading from a single VCF, then create a single SVAR and attach the unparsed 
      // VCF line to that. If needed, it will later be expanded into the consensus.
      //
      
      if ( return_var ) 
	{
	  // Add a sibdummy, single SampleVariant, just to keep everything happy downstream      
	  var.add(1);
	  var.set_vcf_buffer( tok );	
	}
      
      return var;
      
    }

  processVCF( tok(8) , &format );
  

  //
  // Set genotype format 
  //

  set_format( format );

 
  // If we are to add genotypes, we need to have some individuals
  // specified
  
  if ( file_id < 0 ) 
    Helper::halt("No individuals specified but genotypes are present");
  


  //
  // Did we see the correct number of genotypes?
  //
  
  if ( toksize - 9  != icnt )
    {
      plog.warn( "incorrect number of genotypes: " 
		 + Helper::int2str( toksize-9 ) + " observed, " 
		 + Helper::int2str( icnt ) + " expected" ) ;
      
      // if less than we expect/need, do not read line at all
      if ( toksize - 9 < icnt )
	{
	    var.valid( false );
	    return var ;
	}
      // otherwise, okay to read up to end point and ignore rest
    }



  //
  // If reading from a single VCF, then create a single SVAR and attach the unparsed 
  // VCF line to that. If needed, it will later be expanded into the consensus.
  //
  
  if ( return_var ) 
    {
	// Add a sibdummy, single SampleVariant, just to keep everything happy downstream      
	var.add(1);
	var.set_vcf_buffer( tok );
	return var;
    }
  
  
  //
  // Call genotypes, add to variant 
  //
  
  // # of alleles for the svar
  const int na = alt == "." ? 1 : 1 + Helper::char_split( alt , ',' ).size();
  
  int gcnt = 0;
  int idx = 8;

  while ( gcnt != icnt )
  {

      // get next genotype token
      
      Genotype g( tok(++idx) , gt_field , formats , na );
      
      
      // Manually alter genotype to meet X/Y specification?
      
      if ( fixxy_mode ) 
	{
	  
	  // Is this a flagged chromomse?
	  ploidy_t p = mask->ploidy( chr );
	  
	  // AA --> A  | AB --> B  | ./. --> .
	      
	  if ( p == PLOIDY_HAPLOID )
	    {
	      g.make_haploid();
	    }
	  else if (  p == PLOIDY_X 
		     && sex[gcnt] != FEMALE 		     
		     && ! mask->pseudo_autosomal( (const Variant&)var ) )
	    {
	      g.make_haploid();
	    }
	  else if ( p == PLOIDY_Y )
	    {
	      if ( sex[gcnt] != MALE ) 
		g.null( true );
	      else if ( ! mask->pseudo_autosomal( (const Variant&)var ) )
		g.make_haploid();
	    }	  
	}

      var.add( g );

      ++gcnt;
    }
  


  // Return this variant (and any called genotypes)

  return var;
  
}


void VCFReader::set_region_mask( const std::set<Region> * myfilter )
{
  pfilter = (std::set<Region> *)myfilter;
  largest_region = 0;
  std::set<Region>::iterator i = pfilter->begin();
  while ( i != pfilter->end() )
    {
      if ( i->length() > largest_region ) 
	largest_region = i->length();
       ++i;
    }
}

void VCFReader::add_id_filter( const std::set<std::string> & s)
{
  using_idfilter = true;
  std::set<std::string>::iterator ii = s.begin();
  while ( ii != s.end() )
    {
      idfilter.insert( *ii );
      ++ii;
    }
}


bool VCFReader::contains( int chr , int bp1, int bp2 )
{
  Region left( chr, bp1 , bp1 );
  Region right( chr, bp2 , bp2 + largest_region );
    
  std::set<Region>::iterator lwr = pfilter->upper_bound( left );
  std::set<Region>::iterator upr = pfilter->lower_bound( right );
    
  int l = bp1 - largest_region;

  Region left2( chr, l, l );
  
  std::set<Region>::iterator lwr2 = pfilter->upper_bound( left2 );
  
  int cnt = 0 ;
  if ( lwr != lwr2 )
    {	
      while ( lwr2 != lwr )
	{
	  if ( lwr2->stop.position() >= bp1 ) return true;
	  ++lwr2;
	}
    }
  
  while ( lwr != upr ) 
    {
      return true;
      ++lwr;
    }
  
  return false;
}



inline bool VCFReader::processVCF( const char * c , int * a)
{ 
    errno = 0;
    char * endptr;
    *a = strtol( c , &endptr , 10 );
    if ( *endptr == '\0' ) return true;
    *a = 0; 
    return false;
    
}

inline bool VCFReader::processVCF( const char * c , double * a) 
{ 
    char * endptr;
    errno = 0;
    *a = strtod( c , &endptr );
    if ( *endptr == '\0' ) return true;
    *a = 0; 
    return false;    
}

inline bool VCFReader::processVCF( const char * c , std::string * a)
{ 
    *a = c;
    return true;
}

bool VCFReader::set_format( const std::string & f )
{

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



void VCFReader::set_fixxy( Mask * m , LocDBase * p , IndDBase * pi )
{
  mask = m;
  locdb = p;
  inddb = pi;
  fixxy_mode = true;
}

VCFReader::~VCFReader()
{
  if ( (!from_stdin) && file.is_open() ) file.close();
}
