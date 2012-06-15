#include "meta.h"
#include <iostream>

using namespace std;
using namespace Helper;

std::set<std::string> MetaMeta::pop_static;
std::set<std::string> MetaMeta::show_mask;
std::set<std::string> MetaMeta::hide_mask;
std::set<std::string> MetaMeta::internal_mask;

bool MetaMeta::masked_show = false;
bool MetaMeta::masked_hide = false;
bool MetaMeta::masked_internal = false;
bool MetaMeta::force_consensus_mode = false;


void MetaMeta::load( const std::string & filename)
{

  Helper::checkFileExists( filename );
  
  std::ifstream IN1( filename.c_str() , std::ios::in );
  
  while ( ! IN1.eof() )
    {
      
      std::string a,b;

      IN1 >> a >> b;
      
      // ATTRIB  FLAG
      if ( b == "STATIC" ) pop_static.insert(a);
      
      // Some other defaults
      if ( b == "ANNOT" ) PLINKSeq::META_ANNOT()      = a;
      if ( b == "ANNOT_FLAG" ) PLINKSeq::META_ANNOT_FLAG() = a;

      if ( b == "GENE" )  PLINKSeq::META_GENE()       = a;
      if ( b == "GL" )    PLINKSeq::META_GENO_LIK()   = a;
      if ( b == "PL" )    PLINKSeq::META_GENO_PHRED() = a;

      if ( b == "DP" )    PLINKSeq::META_DP() = a;
      if ( b == "GENO_DP" )    PLINKSeq::META_GENO_DP() = a;
      if ( b == "GENO_GQ" )    PLINKSeq::META_GENO_GQ() = a;

      // Phenotype label
      if ( b == "PHENOTYPE" || b == "PHE" || b == "PHENO" )
	PLINKSeq::DEFAULT_PHENOTYPE() = a;
      
      // default locdb group
      if ( b == "LOCGROUP" ) 
	PLINKSeq::DEFAULT_LOC_GROUP() = a;

      if ( b == "G1K" ) PLINKSeq::DEFAULT_G1K() = a;
      if ( b == "DBSNP" ) PLINKSeq::DEFAULT_DBSNP() = a;
      
      // HIDE or SHOW
      if ( b == "HIDE" )
	{
	  masked_hide = true;
	  hide_mask.insert(a);
	}
      else if ( b == "SHOW" )
	{
	  masked_show = true;
	  show_mask.insert(a);
	}
    }
  
  IN1.close();
  
}


void MetaMeta::clear()
{
  pop_static.clear();
}


bool MetaMeta::static_variant( const std::string & attrib )
{
    return pop_static.find( attrib ) != pop_static.end() ;
}
  

void MetaMeta::hide( const std::string & attrib )
{
  masked_hide = true;
  hide_mask.insert( attrib );
  // ensure that internal fields are kept  
}


void MetaMeta::show( const std::string & attrib )
{
  masked_show = true;  
  show_mask.insert( attrib );
}

void MetaMeta::is_internal( const std::string & attrib )
{
  masked_internal = true;
  internal_mask.insert( attrib );
}


bool MetaMeta::display( const std::string & attrib )
{  
  if ( masked_show     &&  show_mask.find( attrib ) == show_mask.end() ) return false;
  if ( masked_hide     &&  hide_mask.find( attrib ) != hide_mask.end() ) return false;    
  if ( masked_internal &&  internal_mask.find( attrib ) != internal_mask.end() ) return false;
  return true;
}


void registerMetatype( const std::string & name, mType mt, int num , int grp , const std::string & desc )
{
  
  // Helper function -- if we have read back metainformation from a
  // database, the source context might not make it obvious what group
  // it belongs to, e.g. variant vs. genotype Therefore, use the
  // explicit coding in this case
  
  switch ( grp ) { 
  case 1 :     
    MetaInformation<VarMeta>::field( name, mt , num, desc );
    break;
  case 2 :
    MetaInformation<GenMeta>::field( name, mt , num, desc );
    break;
  case 3 :
    MetaInformation<LocMeta>::field( name, mt , num, desc );
    break;
  case 4 :
    MetaInformation<RefMeta>::field( name, mt , num, desc );
    break;
  case 5 :
    MetaInformation<FileMeta>::field( name, mt , num, desc );
    break;
  case 6 :
    MetaInformation<IndivMeta>::field( name, mt , num, desc );
    break;
  case 7 :
    MetaInformation<AlleleMeta>::field( name, mt , num, desc );
    break;
  case 8 :
    MetaInformation<VarFilterMeta>::field( name, mt , num, desc );
    break;   
  default :
    MetaInformation<MiscMeta>::field( name, mt , num, desc );
  }
  
}


std::string Helper::metatype_summary( const bool ugly )
{
  std::stringstream ss;

  if ( ugly ) 
      ss << MetaInformation<VarMeta>::list_fields("META_VARIANT")
	 << MetaInformation<VarFilterMeta>::list_fields("META_FILTER")
	 << MetaInformation<GenMeta>::list_fields("META_GENOTYPE")
	 << MetaInformation<LocMeta>::list_fields("META_LOCUS")
	 << MetaInformation<RefMeta>::list_fields("META_REFVAR")
	 << MetaInformation<FileMeta>::list_fields("META_FILE")
	 << MetaInformation<IndivMeta>::list_fields("META_INDIV")
	 << MetaInformation<AlleleMeta>::list_fields("META_ALLELE");
  else
    {
      ss << "---Meta-information summary---\n\n";
      ss << MetaInformation<VarMeta>::pretty_list_fields("Variants")
	 << MetaInformation<VarFilterMeta>::pretty_list_fields("Filters")
	 << MetaInformation<GenMeta>::pretty_list_fields("Genotypes")
	 << MetaInformation<LocMeta>::pretty_list_fields("Locus")
	 << MetaInformation<RefMeta>::pretty_list_fields("Reference variants")
	 << MetaInformation<FileMeta>::pretty_list_fields("Files")
	 << MetaInformation<IndivMeta>::pretty_list_fields("Individuals")
	 << MetaInformation<AlleleMeta>::pretty_list_fields("Alleles");
    }

  return ss.str();

  // MetaInformation<MiscMeta>::list_fields("MISC");

}
