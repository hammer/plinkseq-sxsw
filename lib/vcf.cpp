
#include "vcf.h"
#include "gstore.h"

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

//
// VCFv3.3
//

// Keywords
//   INFO
//   FILTER
//   FORMAT

// Example:
//  ##fileformat=VCFv3.3
//  ##fileDate=20090805
//  ##source=myImputationProgramV3.1
//  ##reference=1000GenomesPilot-NCBI36
//  ##phasing=partial
//  ##INFO=NS,1,Integer,"Number of Samples With Data"
//  ##INFO=DP,1,Integer,"Total Depth"
//  ##INFO=AF,-1,Float,"Allele Frequency"
//  ##INFO=AA,1,String,"Ancestral Allele"
//  ##INFO=DB,0,Flag,"dbSNP membership, build 129"
//  ##INFO=H2,0,Flag,"HapMap2 membership"
//  ##FILTER=q10,"Quality below 10"
//  ##FILTER=s50,"Less than 50% of samples have data"
//  ##FORMAT=GT,1,String,"Genotype"
//  ##FORMAT=GQ,1,Integer,"Genotype Quality"
//  ##FORMAT=DP,1,Integer,"Read Depth"
//  ##FORMAT=HQ,2,Integer,"Haplotype Quality"

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

   
// CHROM  1..22, X, Y, etc

// POS base-1 position 

// ID  rsID, or '.' if none (unlisted)

// REF   Reference base: A,C,G,T,N

// ALT: a comma separated list of alternate non-reference alleles
// called on at least one of the samples. Options are A,C,G,T,Dn
// (for delete n bases starting with the base at POS), I<seq> where
// <seq> is a list of ACGT bases to be inserted just after the base
// at POS, '.' (period character) if there are no alternate alleles.

// QUAL: a phred-scaled quality score for the assertion made in
// ALT. i.e. give -10log_10 prob(call in ALT is wrong). If ALT is
// '.' (no variant) then this is -10log_10 p(variant), and if ALT is
// not '.' this is -10log_10 p(no variant). High QUAL scores
// indicate high confidence calls. Although traditionally people use
// integer phred scores, we agreed to permit floating point scores
// to enable higher resolution for low confidence calls if desired.

// FILTER filter: 0 if this position is not filtered, i.e. a call is
// made at this position. Otherwise a semicolon-separated list of
// codes for filters that fail. e.g. q10;s50 might indicate that
// at this site the quality is below 10 and the number of samples
// with data is below 50% of the total number of samples

// INFO: additional information, encoded as a comma-separated series
// of 2-character keys with optional values in the format:
// <key>=<data>[,data]*. Fields could be e.g.:

// AA ancestral allele, encoded as REF and ALT
// AC allele count in genotypes, for each ALT allele, in the same order as listed
// AN total number of alleles in called genotypes
// AF allele frequency for each ALT allele in the same order as
//   listed: use this when estimated from primary data, not called
//   genotypes    
// DP depth, e.g. D=154    
// MQ RMS mapping quality, e.g. MQ=52
// BQ RMS base quality at this position
// SB strand bias at this position
// DB dbSNP membership  
// H2 membership in hapmap2

// FORMAT: If genotype information is present, then the same types
// of data must be present for all samples. First a FORMAT field is
// given specifying the data types and order. This is followed by
// one field per sample, with the colon-separated data in this field
// corresponding to the types specified in the format. The first
// subfield must always be the genotype.

// GT genotype: 0/0, 0|1, etc  For haploid calls 0.  Missing = .

// GQ genotype quality, encoded as a phred quality
// -10log_10p(genotype call is wrong), max quality 99

// DP read depth at this position for this sample
// HQ haplotype qualities, two phred qualities comma separated

// FT sample genotype filter indicating if this genotype was 
// called (similar in concept to the FILTER record for the
// entire CHROM/POS). Again, use 0 for unfiltered, or a semi-colon
// separated list of codes for filters that fail.

// Additional types could encode probabilities for each genotype
// etc.  If any of the fields is missing, it is replaced by an empty
// string. For example if the format is GT:GQ:DP:HQ then
// A|A::23:23,34 indicates that GQ is missing. Trailing fields can
// be dropped, e.g. A/A:34:23 if HQ is missing.


//---------------------------------------------------------------------------------------------//


// VCF 4.0  (VCF_4_0)

//  1. The REF column can be a String (list of bases).

// 2. The ALT column no longer allows for I/D codes, but does permit angle-bracketed IDs.

// 3. InDels must include the base before the event (and this must
// be reflected by the value in the POS field).

// 4. There is a standard missing "." used.

// In FILTER, "PASS" means is used to signify that a record passes filters instead .
//   --- in 3.4, this was "0", but got automatically changed to "PASS" internally

// 6. The meta information fields in the header use key/value pairs.

//  ##INFO=<ID=NS,Number=1,Type=Integer,Description="Number of Samples With Data">
//  ##INFO=<ID=DP,Number=1,Type=Integer,Description="Total Depth">
//  ##INFO=<ID=AF,Number=-1,Type=Float,Description="Allele Frequency">
//  ##INFO=<ID=AA,Number=1,Type=String,Description="Ancestral Allele">
//  ##INFO=<ID=DB,Number=0,Type=Flag,Description="dbSNP membership, build 129">
//  ##INFO=<ID=H2,Number=0,Type=Flag,Description="HapMap2 membership">
//  ##FILTER=<ID=q10,Description="Quality below 10">
//  ##FILTER=<ID=s50,Description="Less than 50% of samples have data">
//  ##FORMAT=<ID=GT,Number=1,Type=String,Description="Genotype">
//  ##FORMAT=<ID=GQ,Number=1,Type=Integer,Description="Genotype Quality">
//  ##FORMAT=<ID=DP,Number=1,Type=Integer,Description="Read Depth">
//  ##FORMAT=<ID=HQ,Number=2,Type=Integer,Description="Haplotype Quality">

// 7. END and CIGAR have been included in the list of reserved key names for the INFO field. 
//     --- okay, no parsing changes needed

// 8.GL is reserved as a genotype field tag to represent genotype likelihoods. 
//      How to represent this information for multi-allelic sites is still under discussion.
//    --- okay, no parsing changes needed


//---------------------------------------------------------------------------------------------//



VCFReader::line_t VCFReader::parseLine( Variant ** pvar )
{


  if ( return_var ) *pvar = NULL;

  if ( file.eof() ) return VCF_EOF;
  
  // Is this line meta-information, a header, or a variant?
  
  std::string s = file.readLine();

  if ( s == "" ) return VCF_EOF;
  

  //
  // If not a valid line, do nothing
  //
  
  if ( s.size() < 3 ) 
    {
      plog.warn( "invalid line with fewer than 3 characters in VCF: " + s  );
      return VCF_INVALID;
    }
  
  //
  // Either meta-information, header row, or variant row
  //
  
  if ( s.substr(0,2) == "##" ) 
    {
      getMetaInformation(s);
      return VCF_META;
    }
  
  if ( s.substr(0,1) == "#" )
    {
      getHeader(s);
      summary();
      return VCF_HEADER;
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
      else
	plog.warn( "read invalid line from VCF" , s );
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
      else if ( version == VCF_4_0 )
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
       else if ( version == VCF_4_0 )
	 {
	   // VCF4.0:  ##FILTER=<ID=ID,Description=description>

	   tok[1] = Helper::remove_tags( tok[1] );

	   std::map<std::string,std::string> l = Helper::quoted_comma_keypair_split( tok[1] );
	   std::map<std::string,std::string>::iterator i = l.begin();
	   while ( i != l.end() )
	     {
	       if ( i->first == "ID" || i->first == "id" ) name = i->second;
	       else if ( i->first == "Description" || i->first == "description" || i->first == "desc" ) desc = i->second;
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

  while ( tok_iter != tok.end() )
    {
      
      Individual i( *tok_iter );
      
      // Add to the genotype-database
      
      vardb->insert( file_id , i );
      
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

 
  // 1..22, X, Y, etc
  
  std::string chr;  
  

  // base-1 position 

  int pos;     

  
  // rsID, or '.' if none (unlisted)

  std::string id;   


  // Reference base: A,C,G,T,N
  
  std::string ref;  


  // ALT: a comma separated list of alternate non-reference alleles
  // called on at least one of the samples. Options are A,C,G,T,Dn
  // (for delete n bases starting with the base at POS), I<seq> where
  // <seq> is a list of ACGT bases to be inserted just after the base
  // at POS, '.' (period character) if there are no alternate alleles.

  std::string alt;  


  // QUAL: a phred-scaled quality score for the assertion made in
  // ALT. i.e. give -10log_10 prob(call in ALT is wrong). If ALT is
  // '.' (no variant) then this is -10log_10 p(variant), and if ALT is
  // not '.' this is -10log_10 p(no variant). High QUAL scores
  // indicate high confidence calls. Although traditionally people use
  // integer phred scores, we agreed to permit floating point scores
  // to enable higher resolution for low confidence calls if desired.

  double qual;


  // FILTER filter: 0 if this position is not filtered, i.e. a call is
  // made at this position. Otherwise a semicolon-separated list of
  // codes for filters that fail. e.g. ¡Èq10;s50¡É might indicate that
  // at this site the quality is below 10 and the number of samples
  // with data is below 50% of the total number of samples

  std::string filter;

  
  // INFO: additional information, encoded as a comma-separated series
  // of 2-character keys with optional values in the format:
  // <key>=<data>[,data]*. Fields could be e.g.:

  // AA ancestral allele, encoded as REF and ALT
  // AC allele count in genotypes, for each ALT allele, in the same order as listed
  // AN total number of alleles in called genotypes
  // AF allele frequency for each ALT allele in the same order as
  //   listed: use this when estimated from primary data, not called
  //   genotypes    
  // DP depth, e.g. D=154    
  // MQ RMS mapping quality, e.g. MQ=52
  // BQ RMS base quality at this position
  // SB strand bias at this position
  // DB dbSNP membership  
  // H2 membership in hapmap2
  
  std::string info;

  
  // FORMAT: If genotype information is present, then the same types
  // of data must be present for all samples. First a FORMAT field is
  // given specifying the data types and order. This is followed by
  // one field per sample, with the colon-separated data in this field
  // corresponding to the types specified in the format. The first
  // subfield must always be the genotype.


  // GT genotype, encoded as alleles separated by one of /, | and \
  // g. The allele values are 0 for the reference allele (what is in
  // the reference sequence), 1 for the first allele listed in ALT, 2
  // for the second allele list in ALT etc. For haploid calls, e.g. on
  // Y, male X, mitochondrion give only one allele value e.g. 0, for
  // diploid calls examples could be 0/1 or 1|0 etc. The meanings of
  // the separators are:

  //        / : genotype unphased
  //        | : genotype phased
  //        \ : genotype phased but switch probability is high 

  // GQ genotype quality, encoded as a phred quality
  // -10log_10p(genotype call is wrong), max quality 99

  // DP read depth at this position for this sample
  // HQ haplotype qualities, two phred qualities comma separated
  
  // FT sample genotype filter indicating if this genotype was 
  // ¡Ècalled¡É (similar in concept to the FILTER record for the
  // entire CHROM/POS). Again, use 0 for unfiltered, or a semi-colon
  // separated list of codes for filters that fail.

  // Additional types could encode probabilities for each genotype
  // etc.  If any of the fields is missing, it is replaced by an empty
  // string. For example if the format is GT:GQ:DP:HQ then
  // A|A::23:23,34 indicates that GQ is missing. Trailing fields can
  // be dropped, e.g. A/A:34:23 if HQ is missing.
  
  std::string          format;
  

  // Tab-delimited (white-space) ; but do not send trailing \n 
    
  std::vector<std::string> tok = Helper::char_split(  s[ s.size()-1] == '\n' ? 
                                                      s.substr( 0 , s.size() - 1 ) : 
                                                      s , '\t' );  
  
   
  // When reading in from a VCF file, it will, by definition, only
  // correspond to a single sample; therefore, we want to load 
  // up into, and write from, the consensus SampleVariant here

  Variant var( false );
 
  if ( tok.size() < 8 ) return var;

 
  // Get primary identifiers: CHR,POS and ID, etc.  If anything fails
  // here, return a variant in invalid state

  std::vector<std::string>::iterator tok_iter = tok.begin();
  
  if ( ! processVCF( tok_iter , chr ) ) return var; 
  if ( ! processVCF( ++tok_iter , pos ) ) return var; 
  if ( ! processVCF( ++tok_iter , id ) ) return var; 
  if ( ! processVCF( ++tok_iter , ref ) ) return var; 
  if ( ! processVCF( ++tok_iter , alt ) ) return var; 
  if ( ! processVCF( ++tok_iter , qual ) ) return var; 
  if ( ! processVCF( ++tok_iter , filter ) ) return var; 
  if ( ! processVCF( ++tok_iter , info ) ) return var; 
 
  //
  // If applying a mask, do we want to keep this variant?
  //

  if ( pfilter )
    {

      if ( ! contains( Helper::chrCode(chr) , pos , pos ) ) 
	{
	  var.valid( false );
	  return var; 
	}
    }

  // Key variant fields
  var.name( id );
  var.chromosome( Helper::chrCode(chr) );
  var.position( pos );
  var.stop( pos + ref.size() - 1 );

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
  // If SEQDB attached, check that the reference is correct
  //

  if ( seqdb ) 
    {
      std::string sref = seqdb->lookup( chrCode( chr ) , pos , pos + ref.size() - 1 );
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

  if ( tok.size() <= 8 ) return var;
  
  if ( ! processVCF( ++tok_iter , format ) ) return var; 
  




  //
  // Set genotype format 
  //

  set_format( format );


 
  // If we are to add genotypes, we need to have some individuals 
  // specified
  
  if ( file_id < 0 ) 
    Helper::halt("No individuals specified but genotypes are present");
  

  // Parse the format specifier, and attach the the variant, so that
  // genotypes can be called.
  
  VariantSpec * ps = SampleVariant::decoder.decode( format + " " + ref + " " + alt );
  
  var.consensus.specification( ps );


  //
  // Did we see the correct number of genotypes?
  //
  
  if ( tok.size() - 9  != icnt )
    {
      plog.warn( "incorrect number of genotypes: " 
		 + Helper::int2str( tok.size()-9 ) + " observed, " 
		 + Helper::int2str( icnt ) + " expected" ) ;
      var.valid( false );
      return var ;
    }


  //
  // If reading directly from a VCF, do not bother parsing genotypes yet
  // i.e. unless we have to
  //


   if ( return_var ) 
     {
       var.consensus.vcf_direct_buffer = tok;       
       return var;
     }

  
  //
  // Call genotypes, add to variant 
  //
  
  int gcnt = 0;
  while ( ++tok_iter != tok.end() ) 
    {
      Genotype g = ps->callGenotype( *tok_iter , &var ); 
      var.consensus.calls.add( g );
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




bool VCFReader::processVCF( std::vector<std::string>::iterator i , int & a)
{ 
  if ( *i == "." ) { a = -1; return true; } 
  try { a = Helper::lexical_cast<int>( *i ); } 
  catch ( std::exception& e) { return false; } 
  return true;
}


bool VCFReader::processVCF(std::vector<std::string>::iterator i ,
			   double & a) 
{ 
  if ( *i == "." ) { a = -1.0; return true; } 
  try { a = Helper::lexical_cast<double>( *i ); } 
  catch ( std::exception& e) { return false; } 
  return true;
}

bool VCFReader::processVCF(std::vector<std::string>::iterator i , 
			   std::string & a)
{ 
  a = *i;
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
  
  int gi = -1;
  for ( int i = 0 ; i < tok.size() ; i++ )
    {
      
      fset.insert( tok[i] );

      // genotype is not treated as meta-information
      // or if no value set in value  
      
      if ( tok[i] == "GT" ) 
      {
        gi = i; 
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

  // update variant Spec
  VariantSpec::set_format( gi , &formats );

  return true; // indicates a change was made
}

