
#include "mask.h"

#include "locdb.h"
#include "refdb.h"
#include "annot.h"
#include "vardb.h"
#include "gstore.h"

using namespace std;
using namespace Helper;

extern GStore * GP;

bool mask_command_t::operator<( const mask_command_t & rhs ) const
{ 
    if ( group_order < rhs.group_order ) return true;
    if ( group_order > rhs.group_order ) return false;    
    if ( name_order < rhs.name_order ) return true;
    if ( name_order > rhs.name_order ) return false;
    return name < rhs.name;
}

void Mask::searchDB()
{
    // If we do not have explicitly-specified databases, try to attach
    // the generic ones
  
    if ( !vardb ) vardb = GP ? &GP->vardb : NULL ;
    if ( !locdb ) locdb = GP ? &GP->locdb : NULL ;
    if ( !refdb ) refdb = GP ? &GP->refdb : NULL ;
    if ( !segdb ) segdb = GP ? &GP->segdb : NULL ;
}

void mask_add( std::set<mask_command_t> & s , 
	       const int gn , 
	       const int nn , 
	       const std::string & g , 
	       const std::string & n , 
	       const std::string & a , 
	       const std::string & d , 
	       bool h = false ) 
{
  s.insert( mask_command_t(n,nn,g,gn,a,d,h) );
  s.insert( mask_command_t(n) );
}

void mask_group( std::set<mask_group_t> & s , 
		 const std::string & g , 
		 const std::string & d )
{
  s.insert( mask_group_t(g,d) );
}


std::set<mask_command_t> populate_known_commands()
{
  std::set<mask_command_t> s;

  // str             (single string)
  // str-list

  // int
  // int-list        
  // int-range       (single int-range)
  // int-range-list  (comma-sep int-range list)

  // float
  // float-range
  // float-range-list
  
  // flag

  int g = 0;
  int c = 0; 
  std::string gl = "";

  // locus-groups
  ++g; c=0; gl="locus-groups";
  mask_add( s , g , c++ , gl , "loc", "str-list" , "included loci" ); 
  mask_add( s , g , c++ , gl , "loc.inc", "str-list" , "included loci" , true ); //hidden 
  mask_add( s , g , c++ , gl , "loc.group" , "str" , "group variants by locus group" ); 
  mask_add( s , g , c++ , gl , "loc.subset", "str-list" , "subsetted loci" ); 
  mask_add( s , g , c++ , gl , "loc.skip", "str-list" , "skipped loci" ); 
  mask_add( s , g , c++ , gl , "loc.ex", "str-list" , "excluded loci" ); 
  mask_add( s , g , c++ , gl , "loc.req", "str-list" , "required loci");
  mask_add( s , g , c++ , gl , "loc.border" , "float-list" , "3' and 5' borders, bases" );
  mask_add( s , g , c++ , gl , "gene" , "str-list" , "included genes" );
  mask_add( s , g , c++ , gl , "loc.append" , "str-list" , "append meta-information from LOCDB groups" );


  // regions
  ++g; c=0; gl="regions";
  mask_add( s , g , c++ , gl , "reg" , "str-list" , "include region(s)" ); 
  mask_add( s , g , c++ , gl , "reg.inc" , "str-list" , "include regions(s)" , true ); // hidden 
  mask_add( s , g , c++ , gl , "reg.ex" , "str-list" , "require regions(s)" ); 
  mask_add( s , g , c++ , gl , "reg.req" , "str-list" , "exclude regions(s)" ); 
  mask_add( s , g , c++ , gl , "reg.group" , "str-list" , "group variants by list of regions (not implemented)" ); 
  //  mask_add( s , g , c++ , gl , "reg.force" , "str-list" , "force output on variants, whether present or not" ); 

  // variant-IDs
  ++g; c=0; gl="ids";
  mask_add( s , g , c++ , gl , "id" , "str-list" , "include variants based on ID" ); 
  mask_add( s , g , c++ , gl , "id.ex" , "str-list" , "require variants based on ID" ); 
  mask_add( s , g , c++ , gl , "id.req" , "str-list" , "exclude variants based on ID" ); 

  // locus-sets
  ++g; c=0; gl="locus-set-groups";
  mask_add( s , g , c++ , gl , "locset" , "str-list" , "included locus-sets" ); 
  mask_add( s , g , c++ , gl , "locset.inc" , "str-list" , "included locus-sets" , true ); 
  mask_add( s , g , c++ , gl , "locset.group" , "str" , "group variants by locus-set group" ); 
  mask_add( s , g , c++ , gl , "locset.subset", "str-list" , "subsetted locus-sets" ); 
  mask_add( s , g , c++ , gl , "locset.skip" , "str-list" , "skipped locus-sets" ); 
  mask_add( s , g , c++ , gl , "locset.ex" , "str-list" , "excluded locus-sets" ); 
  mask_add( s , g , c++ , gl , "locset.req" , "str-list" , "required locus-sets" );
  mask_add( s , g , c++ , gl , "locset.append" , "str-list" , "append meta-information from LOCDB locus-sets" );


  // variant-groups 
  ++g; c=0; gl="variant-groups";
    
  mask_add( s , g , c++ , gl , "var" , "str-list" , "variant include list(s)" );
  mask_add( s , g , c++ , gl , "var.ex", "str-list" , "excluded variants" ); 
  mask_add( s , g , c++ , gl , "var.req", "str-list" , "required variants" );


//   mask_add( s , g , c++ , gl , "var.inc" , "str-list" , "variant include list(s)" , true ); // hidden

//   mask_add( s , g , c++ , gl , "var.group" , "str" , "group variants by VARDB sets" ); 
   mask_add( s , g , c++ , gl , "var.append" , "str-list" , "append meta-information from VARDB" );

//   mask_add( s , g , c++ , gl , "var.subset", "str-list" , "subsetted variants" , true ); 
//   mask_add( s , g , c++ , gl , "var.skip", "str-list" , "skipped variants" , true ); 

   mask_add( s , g , c++ , gl , "varset", "str-list" , "included variant-supersets" ); 
   mask_add( s , g , c++ , gl , "varset.ex", "str-list" , "excluded variant-supersets");
   mask_add( s , g , c++ , gl , "varset.req", "str-list" , "required variant-supersets");
   mask_add( s , g , c++ , gl , "varset.group" , "str" , "group variants by VARDB super-set" ); 
   mask_add( s , g , c++ , gl , "varset.append" , "str-list" , "append meta-information from VARDB" );

  // Reference databases
  ++g; c = 0 ; gl="ref-variants";
  mask_add( s , g , c++ , gl , "ref" , "str-list" , "include variants in reference variant groups" ); 
  mask_add( s , g , c++ , gl , "ref.ex" , "str-list" , "exclude variants in reference variant groups" );  
  mask_add( s , g , c++ , gl , "ref.req" , "str-list" , "include only variants in all listed reference variant groups" ); 
  mask_add( s , g , c++ , gl , "ref.append" , "str-list" , "append meta-information from REFDB groups" );
  

  // People/sample masks
  ++g; c = 0 ; gl="samples";
  mask_add( s , g , c++ , gl , "file" , "str-list" , "files included" );
  mask_add( s , g , c++ , gl , "file.ex", "str-list" , "files excluded" );
  mask_add( s , g , c++ , gl , "indiv" , "str-list" , "include only these individuals" ); 
  mask_add( s , g , c++ , gl , "indiv.ex" , "str-list" , "exclude these individuals" ); 

  // Variant masks, but on whether present in file, etc  
  ++g; c = 0 ; gl = "files" ;
  mask_add( s , g , c++ , gl , "obs.file" , "str-list" , "include variants observed in file(s)" );
  mask_add( s , g , c++ , gl , "obs.file.req" , "str-list" , "require variants observed in file(s)" );
  mask_add( s , g , c++ , gl , "obs.file.ex", "str-list" , "exclude variants observed in file(s)" );
  
  mask_add( s , g , c++ , gl , "obs.nfile" , "int-range" , "include if observed in [n-m] files" );
  mask_add( s , g , c++ , gl , "alt.nfile" , "int-range" , "include if alt-allele present in [n-m] files" );
  mask_add( s , g , c++ , gl , "fail.nfile" , "int-range" , "include if variant FILTERed in [n-m] files" );

  // as above, but whether alt-allele in ... etc
  mask_add( s , g , c++ , gl , "alt.file" , "str-list" , "include variants with alt-alleles in file(s)" );
  mask_add( s , g , c++ , gl , "alt.file.req" , "str-list" , "require variants have alt-alleles in file(s)" );
  mask_add( s , g , c++ , gl , "alt.file.ex" , "str-list" , "exclude variants with alt-alleles in file(s)" );

  // as above, but with arbitrary groups of people
//   mask_add( s , g , c++ , gl , "alt.group" , "str-list" , "not implemented" , true ); //hidden
//   mask_add( s , g , c++ , gl , "alt.group.req" , "str-list" , "not implemented" , true ); //hidden
//   mask_add( s , g , c++ , gl , "alt.group.ex" , "str-list" , "not implemented" , true ); //hidden
  

  // FILTERs and QUALs
  ++g; c = 0 ; gl="filters";
  mask_add( s , g , c++ , gl , "qual", "float-range" , "include variants with QUAL in [n-m]" );
  mask_add( s , g , c++ , gl , "filter" , "str-list" , "include variants meeting at least one FILTER criterion" ); 
  mask_add( s , g , c++ , gl , "filter.ex" , "str-list" , "exclude variants meeting at least one FILTER criterion" ); 
  mask_add( s , g , c++ , gl , "filter.req" , "str-list" , "include variants meeting all FILTER criteria" );
  mask_add( s , g , c++ , gl , "any.filter" , "flag" , "include variants with any non-PASS FILTER" ); 
  mask_add( s , g , c++ , gl , "any.filter.ex" , "flag" , "exclude variants with any non-PASS FILTER" );
  mask_add( s , g , c++ , gl , "any.filter.req" , "flag" , "include only variants with all non-PASS FILTERs" );
  mask_add( s , g , c++ , gl , "fail.single.file" , "flag" , "only include variants PASSing in all files it is observed in" ); 


  // variant meta-information
  ++g; c = 0 ; gl="vmeta";
  mask_add( s , g , c++ , gl , "include" , "expr" , "filter sample-variants based on logical expression" );
  mask_add( s , g , c++ , gl , "eval" , "expr" , "evaluate sample-variants expression" );
  mask_add( s , g , c++ , gl , "v-include" , "expr" , "filter variants based on logical expression" );
  mask_add( s , g , c++ , gl , "meta" , "str-list" , "include variants passing any meta-field criterion" );
  mask_add( s , g , c++ , gl , "meta.req"  , "str-list" , "require variants passing all meta-field criteria" ); 
  mask_add( s , g , c++ , gl , "meta.attach" , "str-list" , "attach meta-fields uploaded to VARDB" );


  // Presence/frequency masks
  
  // Frequency 
  ++g; c = 0 ; gl="frequency";
  mask_add( s , g , c++ , gl , "mac" , "int-range" , "include variants with minor allele counts between [n-m]" ); 
  mask_add( s , g , c++ , gl , "maf" , "float-range" , "include variants with minor allele frequency in [n-m]" );
  mask_add( s , g , c++ , gl , "biallelic" , "flag" , "include only biallelic sites" ); 
  mask_add( s , g , c++ , gl , "biallelic.ex" , "flag" , "exclude biallelic sites" ); 
  mask_add( s , g , c++ , gl , "monomorphic" , "flag" , "include only monomorphic sites" ); 
  mask_add( s , g , c++ , gl , "monomorphic.ex" , "flag" , "exclude monomorphic sites" );
  mask_add( s , g , c++ , gl , "hwe" , "float-range" , "include variants with HWE p-value in [n-m]" ); 
  mask_add( s , g , c++ , gl , "indel" , "flag" , "only consider indels" );
  mask_add( s , g , c++ , gl , "indel.ex" , "flag" , "exclude indels" );
  mask_add( s , g , c++ , gl , "novel" , "flag" , "only novel SNPs (ID '.')" );
  mask_add( s , g , c++ , gl , "novel.ex" , "flag" , "exclude novel SNPs (ID '.')" );

  // Genotype 
  ++g; c = 0 ; gl="genotype";
  mask_add( s , g , c++ , gl , "geno" , "str-list" , "retain genotypes passing any meta-field criterion" ); 
  mask_add( s , g , c++ , gl , "geno.req" , "str-list" , "retain genotypes passing all meta-field criteria" ); 
  mask_add( s , g , c++ , gl , "null" , "int-range" , "include variants with number of null genotypes in [n-m]" ); 
  mask_add( s , g , c++ , gl , "null.prop" , "float-range" , "include variants with proportion of null genotypes in [n-m]" ); 
  mask_add( s , g , c++ , gl , "assume-ref" , "flag" , "assume null/missing genotypes are reference" );
  mask_add( s , g , c++ , gl , "soft-ref" , "flag" , "0/REF allele initially implies only absence of ALT" );
  mask_add( s , g , c++ , gl , "hard-ref" , "flag" , "0/REF allele always implies presence of REF allele" );
  mask_add( s , g , c++ , gl , "fix-xy" , "str" , "fix X/Y genotypes when loading a VCF" );
  mask_add( s , g , c++ , gl , "genotype-model" , "str" , "set genotype scoring model" );
  mask_add( s , g , c++ , gl , "seg" , "str-list" , "include segment mask(s)" );
  mask_add( s , g , c++ , gl , "seg.req" , "str-list" , "require segment mask(s)" );
  mask_add( s , g , c++ , gl , "seg.ex" , "str-list" , "exclude segment mask(s)" );

  
  // Phenotype
  ++g; c = 0 ; gl="phenotype";
  mask_add( s , g , c++ , gl , "phe" , "str-list" , "include individuals meeting at least one phenotype criterion" ); 
  mask_add( s , g , c++ , gl , "phe.ex" , "str-list" , "exclude individuals meeting at least one phenotype criterion" ); 
  mask_add( s , g , c++ , gl , "phe.req" , "str-list" , "require individuals meet all phenotype criteria" ); 
  mask_add( s , g , c++ , gl , "phe.obs" , "str" , "require phenotype observed" , true ); // hidden  
  mask_add( s , g , c++ , gl , "phe.allow.missing" , "flag" , "do not exclude individuals with missing --pheno" );


  // Case/control
  ++g; c = 0 ; gl="case-control";
  mask_add( s , g , c++ , gl , "case.control" , "int-range-list" , "include variants with alt-alleles [n-m] cases and [i-j] controls" ); 
  mask_add( s , g , c++ , gl , "case.uniq" , "int-range" , "include variants with alt-alleles in [n-m] cases, no controls" ); 
  mask_add( s , g , c++ , gl , "control.uniq" , "int-range" , "include variants with alt-alleles in [n-m] controls, no cases" ); 

 
  // Annotation
  ++g; c = 0 ; gl="annotation";  
  mask_add( s , g , c++ , gl , "annot" , "str" , "include variants with coding annotation");
  mask_add( s , g , c++ , gl , "annot.ex" , "str" , "exclude variants with coding annotations");
  mask_add( s , g , c++ , gl , "annot.req" , "str" , "require variant have coding annotations");
  mask_add( s , g , c++ , gl , "annot.append" , "str" , "append annotations");


  // Skip data
  ++g; c = 0 ; gl="skip";
  mask_add( s , g , c++ , gl , "no-vmeta" , "flag"       , "skip all variant meta-data" ); 
  mask_add( s , g , c++ , gl , "no-geno"  , "flag"       , "skip all genotype data" );  
  mask_add( s , g , c++ , gl , "no-gmeta" , "flag"       , "skip all genotype meta-data" ); 
  mask_add( s , g , c++ , gl , "load-vmeta" , "str-list" , "load only these variant meta-fields" );
  mask_add( s , g , c++ , gl , "load-gmeta" , "str-list" , "load only these genotype meta-fields" );

      
  // Misc.
  ++g; c = 0 ; gl="misc-masks";
  mask_add( s , g , c++ , gl , "em" , "float" , "apply GL/PL-based EM; keep genotypes with prob > n" );
  mask_add( s , g , c++ , gl , "all.group" , "flag" , "all filtered variants as a single group" );
  mask_add( s , g , c++ , gl , "empty.group" , "flag" , "in group-iteration, include groups with no variants" );
  mask_add( s , g , c++ , gl , "limit" , "int" , "limit iteration to first n results" );

  mask_add( s , g , c++ , gl , "downcode" , "flag" , "do not downcode multi-allelic variants as k-1 biallelic variants" );
  mask_add( s , g , c++ , gl , "collapse" , "flag" , "represent multi-allelic variants as one REF/non-REF variant" );
  
  mask_add( s , g , c++ , gl , "no-merge" ,    "flag" , "do not merge any sites at all" );
  mask_add( s , g , c++ , gl , "exact-merge" , "flag" , "only merge variants with same REF/ALT alleles" );
  mask_add( s , g , c++ , gl , "any-merge" ,   "flag" , "combine overlapping variants" );
  
  mask_add( s , g , c++ , gl , "ex-vcf" , "flag" , "name of external VCF" , true ); // hidden 

  return s;
}


std::set<mask_group_t> populate_known_groups()
{
  std::set<mask_group_t> s;
  mask_group( s , "locus-groups" , "Interval-based masks involving LOCDB" );
  mask_group( s , "regions"      , "Interval-based masks specified on the command line" );
  mask_group( s , "locus-set-groups" , "Masks based on sets of loci from a LOCDB" );
  mask_group( s , "ref-variants" , "Masks based on REFDB variants" );
  mask_group( s , "samples"   , "Include/exclude individuals/files" );
  mask_group( s , "files"     , "Include/exclude variants based on presence in one or more files" );  
  mask_group( s , "filters"   , "Masks based on the FILTER and QUAL fields" );
  mask_group( s , "vmeta"     , "Masks based on a variant's meta-information (INFO field)" );
  mask_group( s , "frequency" , "Masks based on variant allele frequency" );
  mask_group( s , "genotype"  , "Per-genotype masks and behaviours" );  
  mask_group( s , "phenotype" , "Individual masks based on phenotypes from INDDB" );
  mask_group( s , "case-control" , "Individual masks based on a disease phenotype from INDDB" );
  mask_group( s , "annotation"  , "Masks based on LOCDB transcript annotation (under revision)" );
  mask_group( s , "skip"        , "Options to skip reading certain things (improves speed)" );
  mask_group( s , "misc-masks"  , "Various other masks" );
  return s;
}


std::set<mask_command_t> Mask::known_commands = populate_known_commands();
std::set<mask_group_t> Mask::known_groups = populate_known_groups();


std::string Mask::describe_options()
{
  std::stringstream ss;
  std::set<mask_command_t>::iterator i = known_commands.begin();
  while ( i != known_commands.end() )
    {
      if ( ! i->hidden )
	ss << i->group << "\t" 
	   << i->name << "\t"
	   << i->argtype << "\t"
	   << i->desc << "\n";
      ++i;
    }
  return ss.str();
}


std::string Mask::list_groups( bool verbose )
{
  std::stringstream ss;

  std::set<mask_group_t>::iterator i = known_groups.begin();
  while ( i != known_groups.end() )
    {

      ss << "\t" << i->name;
      if ( i->name.size() < 8 ) ss << "\t";
      if ( i->name.size() < 16 ) ss << "\t";
      ss << "\t" << i->desc << "\n";
    
      if ( verbose )
	{
	  ss << "\t---------------------------------------------------------\n";
	  ss << list_masks( i->name ) << "\n";
	}

      ++i;
    }
  return ss.str();
}

std::string Mask::list_masks( const std::string & g )
{
  std::stringstream ss;
  std::set<mask_command_t>::iterator i = known_commands.begin();
  while ( i != known_commands.end() )
    {
      if ( ! i->hidden )
	{
	  if ( i->group == g ) 
	    {
	      ss << "\t" << i->name;
	      int len = i->name.size();
	      if ( i->argtype != "flag" ) { len += i->argtype.size() + 5; ss << " { " << i->argtype << " }"; } 
	      if ( len < 8 ) ss << "\t";
	      if ( len < 16 ) ss << "\t";
	      if ( len < 24 ) ss << "\t";
	      ss << "\t" << i->desc << "\n";
	    }
	}
      ++i;
    }
  return ss.str();
}


Mask::Mask( const std::string & d , const std::string & expr , const bool filter_mode , const bool group_mode ) 
    : vardb(NULL) , locdb(NULL) , refdb(NULL) , segdb(NULL) , group_mode( group_mode ) 
{
  

 
  //
  // Swap in any file-lists with @includes
  //

  std::string d2 = Helper::filelist2commalist(d);


  //
  // Attach necessary databases, set sensible defaults, etc
  //
  
  construct();
  

  //
  // Is any expression set? (for SampleVariants)
  //

  if ( expr != "" )
    {
      set_filter_expression( expr );
      eval_filter_includes = filter_mode; // T : include=X  
                                          // F : exclude=X
    }
  
  //
  // Parse command line option as meta-field
  //
  
  MetaInformation<MiscMeta> m;


  // Pass command line ( key=val key key=val1,val2 )

  m.parse( d2, ' ', true ); // automatically add meta-flags 
  
 
  std::vector<std::string> keys = m.keys();
  
  for( int i=0; i<keys.size(); i++)
  {	
      
      if ( keys[i] != "" && known_commands.find( mask_command_t( keys[i] ) ) == known_commands.end() )
      {
	
	// always bail here, i.e. for reg=chr1 syntax rather than --mask chr1
	Helper::halt("Mask option '" + keys[i] + "' not recognised.");
	
	// If we can interpret this as a region, do so

// 	bool okay = false;
// 	Region r( keys[i] , okay );	  
// 	if ( !okay ) Helper::halt("Mask option '" + keys[i] + "' not recognised.");
// 	else keys[i] = "reg=" + keys[i];

      }
  }
  
  if ( m.has_field( "v-include" ) )
    { 
      std::string k = m.get1_string( "v-include" );            
      var_set_filter_expression( k );
    }

  if ( m.has_field( "var" ) ) 
    {
      std::vector<std::string> k = m.get_string( "var" );
      for (int i=0; i<k.size(); i++) include_var(k[i]);
    }
  
  if ( m.has_field( "var.ex" ) ) 
    {
      std::vector<std::string> k = m.get_string( "var.ex" );
      for (int i=0; i<k.size(); i++) exclude_var(k[i]);
    }

  if ( m.has_field( "var.req" ) ) 
    {
      std::vector<std::string> k = m.get_string( "var.req" );
      for (int i=0; i<k.size(); i++) require_var(k[i]);
    }
  

  // variant super-sets

  if ( m.has_field( "varset" ) ) 
    {
      std::vector<std::string> k = m.get_string( "varset" );
      for (int i=0; i<k.size(); i++) include_varset(k[i]);
    }
  
  if ( m.has_field( "varset.ex" ) ) 
    {
      std::vector<std::string> k = m.get_string( "varset.ex" );
      for (int i=0; i<k.size(); i++) exclude_varset(k[i]);
    }

  if ( m.has_field( "varset.req" ) ) 
    {
      std::vector<std::string> k = m.get_string( "varset.req" );
      for (int i=0; i<k.size(); i++) require_varset(k[i]);
    }
  


//   if ( m.has_field( "var.subset" ) ) 
//     {
//       std::vector<std::string> k = m.get_string( "var.subset" );
//       for (int i=1; i<k.size(); i++) subset_var(k[0],k[i]);
//     }
  
//   if ( m.has_field( "var.skip" ) ) 
//     {
//       std::vector<std::string> k = m.get_string( "var.skip" );
//       for (int i=1; i<k.size(); i++) skip_var(k[0],k[i]);
//     }
  
  
  if ( m.has_field( "loc" ) ) 
    {
      std::vector<std::string> k = m.get_string( "loc" );
      for (int i=0; i<k.size(); i++) include_loc(k[i]);
    }
  
  if ( m.has_field( "loc.ex" ) ) 
  {
      std::vector<std::string> k = m.get_string( "loc.ex" );
      for (int i=0; i<k.size(); i++) exclude_loc(k[i]);
    }
  
  if ( m.has_field( "loc.req" ) ) 
    {
      std::vector<std::string> k = m.get_string( "loc.req" );
      for (int i=0; i<k.size(); i++) require_loc(k[i]);
    }


  if ( m.has_field( "loc.border" ) )
  {
      
      // TODO...
      Helper::halt( "loc.border not yet implemented" );
//       std::vector<std::string> k = m.get_string( "loc.border" );
//       if ( k.size() != 1 && k.size() != 2 ) Helper::halt( "expecting one or two values, e.g. loc.border=50,50" );
//       int a, b;
//       if ( k.size() == 2 ) 
//       {
// 	  if ( ! ( Helper::str2int( k[0], a ) || Helper::str2int( k[1], b ) ) )
// 	      Helper::halt( "expecting two values, e.g. loc.border=50,50" );
// 	  loc_border_3prime( a );
// 	  loc_border_5prime( b );
//       }
//       else if ( k.size() == 1 ) 
//       {
// 	  if ( ! ( Helper::str2int( k[0], a ) ) )
// 	      Helper::halt( "expecting one or two values, e.g. loc.border=50,50" );
// 	  loc_border_3prime( a );
// 	  loc_border_5prime( a );
//       }

  }


  if ( m.has_field( "seg" ) ) 
    {
      std::vector<std::string> k = m.get_string( "seg" );
      for (int i=0; i<k.size(); i++) include_seg(k[i]);
    }
  
  if ( m.has_field( "seg.req" ) ) 
    {
      std::vector<std::string> k = m.get_string( "seg.req" );
      for (int i=0; i<k.size(); i++) require_seg(k[i]);
    }

  if ( m.has_field( "seg.ex" ) ) 
    {
      std::vector<std::string> k = m.get_string( "seg.ex" );
      for (int i=0; i<k.size(); i++) exclude_seg(k[i]);
    }

  if ( m.has_field( "loc.subset" ) ) 
    {
      std::vector<std::string> k = m.get_string( "loc.subset" );
      for (int i=1; i<k.size(); i++) subset_loc(k[0],k[i]);      
    }

  if ( m.has_field( "loc.skip" ) ) 
    {
      std::vector<std::string> k = m.get_string( "loc.skip" );
      for (int i=1; i<k.size(); i++) skip_loc(k[0],k[i]);
    }


  // locus-sets

  if ( m.has_field( "locset" ) ) 
    {
      std::vector<std::string> k = m.get_string( "locset" );
      for (int i=1; i<k.size(); i++) include_loc_set(k[0],k[i]);
    }

  if ( m.has_field( "locset.subset" ) ) 
    {
      std::vector<std::string> k = m.get_string( "locset.subset" );
      for (int i=2; i<k.size(); i++) subset_loc_set(k[0],k[1],k[i]);
    }
  
  if ( m.has_field( "locset.skip" ) ) 
    {
      std::vector<std::string> k = m.get_string( "locset.subset" );
      for (int i=2; i<k.size(); i++) skip_loc_set(k[0],k[1],k[i]);
    }
  
  if ( m.has_field( "locset.req" ) ) 
    {
      std::vector<std::string> k = m.get_string( "locset.req" );
      for (int i=1; i<k.size(); i++) require_loc_set(k[0],k[i]);
    }

  if ( m.has_field( "locset.ex" ) ) 
    {
      std::vector<std::string> k = m.get_string( "locset.ex" );
      for (int i=1; i<k.size(); i++) exclude_loc_set(k[0],k[i]);
    }


  // reference-variants

  if ( m.has_field( "ref" ) ) 
    {
      std::vector<std::string> k = m.get_string( "ref" );
      for (int i=0; i<k.size(); i++) include_ref(k[i]);
    }

  if ( m.has_field( "ref.ex" ) ) 
    {
      std::vector<std::string> k = m.get_string( "ref.ex" );
      for (int i=0; i<k.size(); i++) exclude_ref(k[i]);
    }

  if ( m.has_field( "ref.req" ) ) 
    {
      std::vector<std::string> k = m.get_string( "ref.req" );
      for (int i=0; i<k.size(); i++) require_ref(k[i]);
    }
  
  if ( m.has_field( "ex-vcf" ) )
    {
      external_vcf_iteration( m.get1_string( "ex-vcf" ) );
    }


  if ( m.has_field( "file" ) ) 
    {
      std::vector<std::string> k = m.get_string( "file" );
      for (int i=0; i<k.size(); i++) 
	include_file( k[i] );
    }

  if ( m.has_field( "file.ex" ) ) 
    {
      std::vector<std::string> k = m.get_string( "file.ex" );
      for (int i=0; i<k.size(); i++) 
	exclude_file( k[i] );
    }
  
  if ( m.has_field("fail.single.file") )
    {
      fail_on_sample_variant( 0 , 1 ); // (n,m) allow n missing, require m observed
    }

  if ( m.has_field("fail.nfile") ) 
    {
      std::vector<std::string> k = m.get_string( "fail.nfile" );
      if ( k.size() == 2 )
	{
	  int i1, i2;
	  if ( ! ( Helper::str2int( k[0] , i1 ) || Helper::str2int( k[1] , i2 ) ) )
	    plog.warn("incorrect arg for nfile");
	  else
	    fail_on_sample_variant( i1 , i2 ) ; 
	}
      else if ( k.size() == 1 ) // cannot fail more than X, need at least one left
	{
	  int i1;
	  if ( ! Helper::str2int( k[0] , i1 ) )
	    plog.warn("incorrect arg for nfile=n");
	  else
	    fail_on_sample_variant( i1 , 1 ) ; // requires at least 1 file 
	}
      else
	plog.warn("incorrect argument for nfile=n{,m}");       
    }

  if ( m.has_field( "reg" ) ) 
    {
      std::vector<std::string> k = m.get_string( "reg" );
      include_reg( k );
    }

  if ( m.has_field( "reg.ex" ) ) 
    {
      std::vector<std::string> k = m.get_string( "reg.ex" );
      for (int i=0; i<k.size(); i++) 
	{
	  bool okay = false;
	  Region r( k[i] , okay );
	  if ( okay ) exclude_reg( r );
	}
    }

  if ( m.has_field( "reg.req" ) ) 
    {
      std::vector<std::string> k = m.get_string( "reg.req" );
      for (int i=0; i<k.size(); i++) 
	{
	  bool okay = false;
	  Region r( k[i] , okay );
	  if ( okay ) require_reg( r );
	}
    }

  
  // Variant IDs

  if ( m.has_field( "id" ) ) 
    {
      std::vector<std::string> k = m.get_string( "id" );
      include_id( k );
    }
  
  if ( m.has_field( "id.ex" ) ) 
    {
      std::vector<std::string> k = m.get_string( "id.ex" );
      exclude_id( k );
    }
  
  if ( m.has_field( "id.req" ) ) 
    {
      std::vector<std::string> k = m.get_string( "id.req" );
      require_id( k );
    }


  // Variant file masks

  if ( m.has_field( "obs.file" ) )
    include_var_obs_file( m.get_string( "obs.file" ) );
    
  if ( m.has_field( "obs.file.req" ) )
    require_var_obs_file( m.get_string( "obs.file.req" ) );
  
  if ( m.has_field( "obs.file.ex" ) )
    exclude_var_obs_file( m.get_string( "obs.file.ex" ) );


  if ( m.has_field( "alt.file" ) )
    include_var_alt_file( m.get_string( "alt.file" ) );
    
  if ( m.has_field( "alt.file.req" ) )
    require_var_alt_file( m.get_string( "alt.file.req" ) );
  
  if ( m.has_field( "alt.file.ex" ) )
    exclude_var_alt_file( m.get_string( "alt.file.ex" ) );


  if ( m.has_field( "alt.group" ) )
    include_var_alt_group( m.get_string( "alt.group" ) );
    
  if ( m.has_field( "alt.group.req" ) )
    require_var_alt_group( m.get_string( "alt.group.req" ) );
  
  if ( m.has_field( "alt.group.ex" ) )
    exclude_var_alt_group( m.get_string( "alt.group.ex" ) );

  // File counts

  if ( m.has_field( "obs.nfile" ) )
    {
      std::vector<std::string> k = m.get_string( "obs.nfile" );
      if ( k.size() != 1 && k.size() != 2 ) Helper::halt("bad args for obs.nfile");
      int i1, i2 = 0;      
      if ( ! Helper::str2int( k[0] , i1 ) ) 
	plog.warn("invalid arg for obs.nfile=n,{m}");
      if ( k.size() ==2 && ! Helper::str2int( k[1] , i2 ) ) 
	plog.warn("invalid arg for obs.nfile=n,{m}");
      count_obs_file( i1,i2 );
    }

  if ( m.has_field( "alt.nfile" ) )
    {
      std::vector<std::string> k = m.get_string( "alt.nfile" );
      if ( k.size() != 1 && k.size() != 2 ) Helper::halt("bad args for alt.nfile");
      int i1, i2 = 0;      
      if ( ! Helper::str2int( k[0] , i1 ) ) 
	plog.warn("invalid arg for alt.nfile=n,{m}");
      if ( k.size() ==2 && ! Helper::str2int( k[1] , i2 ) ) 
	plog.warn("invalid arg for alt.nfile=n,{m}");
      count_alt_file( i1,i2 );
    }


  // Filters

  if ( m.has_field( "filter" ) )
    {
      include_filter( m.get_string( "filter" ) );
    }

  if ( m.has_field( "filter.ex" ) )
    {
      exclude_filter( m.get_string( "filter.ex" ) );
    }

  if ( m.has_field( "filter.req" ) )
    {
      require_filter( m.get_string( "filter.req" ) );
    }


  //
  // Attach user-defined meta-information
  //

  if ( m.has_field("meta.attach") )
    {
      std::vector<std::string> s = m.get_string( "meta.attach" );
      if ( s.size() == 0 || s[0] == "ALL" ) 
	attach_all_meta(true);
      else 
	for (int i=0; i<s.size(); i++)
	  attach_meta_fields(s[i]);
    }

  
  //
  // Include/exclude biallelic SNPs
  // 

  if ( m.has_field( "biallelic" ) )
    {
      require_biallelic(true);
    }

  if ( m.has_field( "biallelic.ex" ) )
    {
      exclude_biallelic(true);
    }

  if ( m.has_field( "monomorphic" ) )
    {
      require_monomorphic(true);
    }

  if ( m.has_field( "monomorphic.ex" ) )
    {
      exclude_monomorphic(true);
    }

  //
  // SNPs vs indels
  //
  
  if ( m.has_field( "indel" ) )
    {
      only_indels( true );
    }

  if ( m.has_field( "indel.ex" ) )
    {
      skip_indels( true );
    }

  //
  // Novel SNPs
  //
  
  if ( m.has_field( "novel" ) )
    {
      only_novel( true );
    }
    
  if ( m.has_field( "novel.ex" ) )
    {
      skip_novel( true );
    }


  //
  // Variant meta-information
  //

  if ( m.has_field( "meta" ) || m.has_field( "meta.req" ) )
    {
      
      // process includes first, then reqs
      bool req = false;
      while ( 1 ) 
	{
	  
	  if ( ( !req ) && !m.has_field( "meta" ) ) 
	    {
	      req = true;
	      continue;
	    }
	  
	  if ( req && ! m.has_field( "meta.req" ) )
	    {
	      break;
	    }
	  
	  vector<string> s = req ? m.get_string( "meta.req" ) : m.get_string( "meta" );
	  
	  for (int j=0; j<s.size(); j++)
	    {
	      
	      std::vector<std::string> t = parse(s[j],":");
	      
	      // should be multiple of 3
	      if ( t.size() % 3 == 0 ) 
		for (int i=0; i<t.size(); i+=3)
		  {
		    
		    if ( t[i+1] == "eq" )
		      {
			int x;
			if ( str2int(t[i+2],x) )
			  req ? req_meta_equals( t[i] , x ) : meta_equals( t[i] , x );
			else // try as string, is parsed below
			  t[i+1] = "is";
		      }
		    
		    if ( t[i+1] == "ne" )
		      {
			int x;
			if ( str2int(t[i+2],x) )
			  req ? req_meta_not_equals( t[i] , x ) : meta_not_equals( t[i] , x );
			else 
			  t[i+1] == "not";
		      }
		    
		    if ( t[i+1] == "gt" )
		      {
			double x;
			if ( str2dbl(t[i+2],x) )
			  req ? req_meta_greater( t[i] , x ) : meta_greater( t[i] , x );
		      }
		    
		    if ( t[i+1] == "ge" )
		      {
			double x;
			if ( str2dbl(t[i+2],x) )
			  {
			    req ? req_meta_greater_equal( t[i] , x ) : meta_greater_equal( t[i] , x );
			  }
			
		      }
		    
		    if ( t[i+1] == "lt" )
		      {
			double x;
			if ( str2dbl(t[i+2],x) )
			  req ? req_meta_less( t[i] , x ) : meta_less( t[i] , x );
		      }
		    
		    if ( t[i+1] == "le" )
		      {
			double x;
			if ( str2dbl(t[i+2],x) )
			  req ? req_meta_less_equal( t[i] , x ) : meta_less_equal( t[i] , x );
		      }
		    
		    if ( t[i+1] == "is" )
		      {
			string x = t[i+2] == "*" ? "" : t[i+2] ;
			req ? req_meta_equals( t[i] , x ) : meta_equals( t[i] , x );
		      }
		    
		    if ( t[i+1] == "not" )
		      {
			string x = t[i+2] == "*" ? "" : t[i+2] ;
			req ? req_meta_not_equals( t[i] , x ) : meta_not_equals( t[i] , x );
		      }
		    
		  }
	    }

	  if ( ! req ) 
	    req = true;
	  else 
	    break;
	}
    }


  //
  // Genotype meta-information
  //

  if ( m.has_field( "geno" ) || m.has_field( "geno.req" ) )
    {
      
      // process includes first, then reqs
      bool req = false;
      while ( 1 ) 
	{
	  
	  if ( ( !req ) && !m.has_field( "geno" ) ) 
	    {
	      req = true;
	      continue;
	    }
	  
	  if ( req && ! m.has_field( "geno.req" ) )
	    {
	      break;
	    }
	  
	  std::vector<std::string> s = req ? m.get_string( "geno.req" ) : m.get_string( "geno" );
	  
	  for (int j=0; j<s.size(); j++)
	    {
	      std::vector<std::string> t = parse(s[j],":");
	      
	      // should be multiple of 3
	      if ( t.size() % 3 == 0 ) 
		for (int i=0; i<t.size(); i+=3)
		  {
		    
		    if ( t[i+1] == "eq" )
		      {
			int x;
			if ( str2int(t[i+2],x) )
			  req ? req_geno_equals( t[i] , x ) : geno_equals( t[i] , x );
		      }
		    
		    if ( t[i+1] == "ne" )
		      {
			int x;
			if ( str2int(t[i+2],x) )
			  req ? req_geno_not_equals( t[i] , x ) : geno_not_equals( t[i] , x );
		      }
		    
		    if ( t[i+1] == "gt" )
		      {
			double x;
			if ( str2dbl(t[i+2],x) )
			  req ? req_geno_greater( t[i] , x ) : geno_greater( t[i] , x );
		      }
		    
		    if ( t[i+1] == "ge" )
		      {
			double x;
			if ( str2dbl(t[i+2],x) )
			  {
			    req ? req_geno_greater_equal( t[i] , x ) : geno_greater_equal( t[i] , x );
			  }
			
		      }
		    
		    if ( t[i+1] == "lt" )
		      {
			double x;
			if ( str2dbl(t[i+2],x) )
			  req ? req_geno_less( t[i] , x ) : geno_less( t[i] , x );
		      }
		    
		    if ( t[i+1] == "le" )
		      {
			double x;
			if ( str2dbl(t[i+2],x) )
			  req ? req_geno_less_equal( t[i] , x ) : geno_less_equal( t[i] , x );
		      }
		    
		    if ( t[i+1] == "is" )
		      {
			string x = t[i+2] == "*" ? "" : t[i+2] ;
			req ? req_geno_equals( t[i] , x ) : geno_equals( t[i] , x );
		      }
		    
		    if ( t[i+1] == "not" )
		      {
			string x = t[i+2] == "*" ? "" : t[i+2] ;
			req ? req_geno_not_equals( t[i] , x ) : geno_not_equals( t[i] , x );
		      }
		    
		  }
	    }

	  if ( ! req ) 
	    req = true;
	  else 
	    break;
	}
    }

  if ( m.has_field( "em" ) )
    {
      std::string s = m.get1_string("em");
      double t;
      if ( Helper::str2dbl(s,t) )
	{
	  EM_caller( true );
	  EM_threshold( t );
	}
    }
  
  if ( m.has_field( "maf" ) )
    {
      dbl_range r( m.get1_string( "maf" ) , -1 ); // 0.05  means include 0.05 or lower
      minor_allele_frequency( r.has_lower() ? r.lower() : 0 , r.has_upper() ? r.upper() : 1 );
    }
  
  if ( m.has_field( "mac" ) )
    {
      // int_range("2", 0) means 2 -->  *:2
      // int_range("2", 1) means 2 -->  2:*
      // int_range("2",-1) means 2 -->  *:2
      
      int_range r( m.get1_string( "mac" ) , -1 ); // 2 means 2 or *lower*
      minor_allele_count( r.lower() , r.upper() );
    }
 
  if ( m.has_field( "hwe" ) )
    {      
      dbl_range r( m.get1_string( "hwe" ) , 1 ); // x --> include x:*       
      hwe( r.has_lower() ? r.lower() : 0  , r.has_upper() ? r.upper() : 1 );
    }
  

  if ( m.has_field( "assume-ref" ) )
    {
      assuming_null_is_reference( true );
    }


  if ( m.has_field( "soft-ref" ) )
  {
      soft_reference_calls( true );
  }
  else if ( m.has_field( "hard-ref" ) )
  {
      soft_reference_calls( false );
  }


  if ( m.has_field( "fix-xy" ) )
    {
      if ( ! locdb ) 
	{
	  plog.warn("no LOCDB attached for fix-xy");	  
	}
      else
	{
	  std::string g = m.get1_string( "fix-xy" ); 
	  
	  // expected format:

	  //  chrX=X
	  //  chrX:1..3000=PAR
	  //  chrX:150000000..250000000=PAR
	  //  

	  std::vector<std::string> v = locdb->fetch_special( g );

	  if ( v.size() == 0 ) 
	    {
	      plog.warn("no special region/chromosome codes set for fix-xy");	      
	    }
	  else
	    {
	      for ( int i = 0; i < v.size() ; i++ )
		{
		  std::vector<std::string> k = Helper::char_split( v[i] , '=' );
		  if ( k.size() != 2 ) continue;		  

		  fixxy( true );

		  // chr-code or region?
		  if ( k[1] == "PAR" ) 
		    {
		      bool okay = true;
		      Region region( k[0] , okay );
		      if ( okay ) set_pseudo_autosomal( region ); 
		      else plog.warn( "trouble setting PAR" , region.coordinate() );			
		    }		  
		  else if ( k[0].find( ":" ) == std::string::npos )
		    {
		      
		      ploidy_t p = PLOIDY_UNKNOWN;
		      
		      if      ( k[1] == "X" ) ploidy( k[0] , PLOIDY_X );
		      else if ( k[1] == "Y" ) ploidy( k[0] , PLOIDY_Y );
		      else if ( k[1] == "HAPLOID" )  ploidy( k[0] , PLOIDY_HAPLOID );
		      else if ( k[1] == "AUTOSOMAL" )  ploidy( k[0] , PLOIDY_AUTOSOMAL );
		      else 
			{ 
			  plog.warn( "unrecognized ploidy code" , k[1] ); continue; 
			} 
		      		      
		    }
		}
	    }
	}
    }


  if ( m.has_field( "genotype-model" ) )
    {
      std::string model = m.get1_string( "genotype-model" );

      if ( model == "ALLELIC" ) Genotype::model = GENOTYPE_MODEL_ALLELIC;
      else if ( model == "ALLELIC2") Genotype::model = GENOTYPE_MODEL_ALLELIC2;
      else if ( model == "ALLELIC3") Genotype::model = GENOTYPE_MODEL_ALLELIC3;
      else if ( model == "DOM") Genotype::model = GENOTYPE_MODEL_DOM;
      else if ( model == "REC") Genotype::model = GENOTYPE_MODEL_REC;
      else if ( model == "REC2") Genotype::model = GENOTYPE_MODEL_REC2;
      else if ( model == "CN") Genotype::model = GENOTYPE_MODEL_CN;
      else if ( model == "NULL") Genotype::model = GENOTYPE_MODEL_NULL;
      else if ( model == "DOSAGE" ) Genotype::model = GENOTYPE_MODEL_DOSAGE;
      else if ( model == "PROB_REF" ) Genotype::model = GENOTYPE_MODEL_PROB_REF;
      else if ( model == "PROB_HET" ) Genotype::model = GENOTYPE_MODEL_PROB_REF;
      else if ( model == "PROB_HOM" ) Genotype::model = GENOTYPE_MODEL_PROB_REF;
      else Helper::halt( "did not recognize genotype model" );
    }
  

  if ( m.has_field( "downcode" ) )
    downcode( DOWNCODE_MODE_EACH_ALT );
  else if ( m.has_field( "collapse" ) )
    downcode( DOWNCODE_MODE_ALL_ALT );
  else     
    downcode( DOWNCODE_MODE_NONE ); // default -- no downcoding
  
  if ( m.has_field( "no-merge" ) )          // no merging, even if exact match
    mergemode( MERGE_MODE_NONE );
  else if ( m.has_field( "exact-merge" ) )  // only if size of ALTs is similar (is default)
    mergemode( MERGE_MODE_EXACT );
  else if ( m.has_field( "any-merge" ) )
    mergemode( MERGE_MODE_ANY_OVERLAP );    // attempt to merge any calls
  else
    mergemode( MERGE_MODE_EXACT );  


  if ( m.has_field( "null" ) )
    {
      null_filter( m.get1_string( "null" ) );
    }
  
  if ( m.has_field( "null.prop" ) )
    {
      null_prop_filter( m.get1_string( "null.prop" ) );
    }


  if ( m.has_field( "qual" ) )
    {
      qual_filter( m.get1_string( "qual" ) );
    }

  if ( m.has_field( "case.control" ) )
    {
      std::vector<std::string> k = m.get_string( "case.control" );
      if ( k.size() == 2 )
	case_control_filter( k[0] , k[1] );
    }
  
  if ( m.has_field( "case.uniq" ) )
     {
       case_control_filter( "1:" , "0:0" );
     }

   if ( m.has_field( "control.uniq" ) )
     {
       case_control_filter(  "0:0" , "1:" );
     }
   
    
  if ( m.has_field( "var.append" ) ) 
    {
      std::vector<std::string> k = m.get_string( "var.append" );
      for (int i=0; i<k.size(); i++) append_var(k[i]);
    }

  if ( m.has_field( "varset.append" ) ) 
    {
      std::vector<std::string> k = m.get_string( "varset.append" );
      for (int i=0; i<k.size(); i++) append_var_set(k[i]);
    }

  if ( m.has_field( "loc.append" ) ) 
    {
      vector<string> k = m.get_string( "loc.append" );
      for (int i=0; i<k.size(); i++) append_loc(k[i]);
    }

  if ( m.has_field( "locset.append" ) ) 
    {
      vector<string> k = m.get_string( "locset.append" );
      for (int i=1; i<k.size(); i++) append_loc_set(k[0],k[i]);
    }

  if ( m.has_field( "ref.append" ) ) 
    {
      vector<string> k = m.get_string( "ref.append" );
      for (int i=0; i<k.size(); i++) append_ref(k[i]);
    }
  



  if ( m.has_field( "annot" ) ) 
    {
      vector<string> k = m.get_string( "annot" );
      if ( k.size() == 1 && k[0] == "nsSNP" )
	include_annotation_nonsyn();
      else
	for (int i=0; i<k.size(); i++) include_annotation(k[i]);
    }


  if ( m.has_field( "annot.ex" ) ) 
    {
      vector<string> k = m.get_string( "annot.ex" );
      if ( k.size() == 1 && k[0] == "nsSNP" )
	exclude_annotation_nonsyn();
      else
	for (int i=0; i<k.size(); i++) exclude_annotation(k[i]);
    }


  if ( m.has_field( "annot.req" ) ) 
    {
      vector<string> k = m.get_string( "annot.req" );
      if ( k.size() == 1 && k[0] == "nsSNP" )
	{
	  require_annotation_nonsyn();
	}
      else
	for (int i=0; i<k.size(); i++) require_annotation(k[i]);
    }

  if ( m.has_field( "annot.append" ) ) 
    {
      append_annotation();
    }

  if ( m.has_field( "empty.group" ) ) 
    {
      process_empty_groups( true );
    }

  if ( m.has_field( "all.group" ) )
    {
      all_group = true;
    }

//   if ( m.has_field( "var.group" ) ) 
//     {
//       std::vector<std::string> k = m.get_string( "var.group" );
//       if ( k.size() > 0 ) group_var( k[0] );
//     }

  if ( m.has_field( "varset.group" ) ) 
    {
      std::vector<std::string> k = m.get_string( "varset.group" );
      if ( k.size() > 0 ) group_var_set( k[0] );
    }

  if ( m.has_field( "loc.group" ) ) 
    {
      std::vector<std::string> k = m.get_string( "loc.group" );
      if ( k.size() > 0 ) group_loc( k[0] );
      
    }


  // convenience subset function: assumes refseq
  // search on altnames first (i.e. gene symbol to transcript ID)
  // two forms: include and group

  if ( m.has_field( "gene" ) )
    {
      if ( group_mode ) group_loc( PLINKSeq::DEFAULT_LOC_GROUP() );
      std::vector<std::string> k = m.get_string( "gene" );
      if ( k.size() > 0 ) 
	for ( int i = 0 ; i < k.size() ; i++ )
	  subset_loc_altname( PLINKSeq::DEFAULT_LOC_GROUP() , k[i] );
    }

  if ( m.has_field( "reg.group" ) )
    {
      std::vector<std::string> k = m.get_string( "reg.group" );
      if ( k.size() > 0 ) group_reg( k );
    }
  
  if ( m.has_field( "locset.group" ) ) 
    {
      std::vector<std::string> k = m.get_string( "locset.group" );
      if ( k.size() == 2 ) group_loc_set( k[0] , k[1] );
      else plog.warn("invalid arg for locset.group");
    }
 
  if ( m.has_field( "filter" ) ) 
    {
      std::vector<std::string> k = m.get_string( "filter" );
      for (int i=0; i<k.size(); i++) include_filter( k[i] );
    }

  if ( m.has_field( "filter.ex" ) ) 
    {
      std::vector<std::string> k = m.get_string( "filter.ex" );
      for (int i=0; i<k.size(); i++) exclude_filter( k[i] );
    }


  if ( m.has_field( "filter.req" ) ) 
    {
      std::vector<std::string> k = m.get_string( "filter.req" );
      for (int i=0; i<k.size(); i++) require_filter( k[i] );
    }

  if ( m.has_field( "any.filter" ) ) 
    {
      include_filter_any();
    }

  if ( m.has_field( "any.filter.req" ) ) 
    {
      require_filter_any();
    }

  if ( m.has_field( "any.filter.ex" ) ) 
    {
      exclude_filter_any();
    }

  if ( m.has_field( "limit" ) ) 
    {
      int n;
      string k = m.get1_string( "limit" );
      if ( str2int(k,n) ) limit(n);
    }


  if ( m.has_field( "indiv" ) ) 
    {
      include_indiv( m.get_string( "indiv" ) );
    }

  if ( m.has_field( "indiv.ex" ) ) 
    {
      exclude_indiv( m.get_string( "indiv.ex" ) );
    }

  // format phe=batch:2,scz=2:1,

  if ( m.has_field( "phe" ) ) 
    {
      include_phenotype( m.get_string( "phe" ) );
    }

  if ( m.has_field( "phe.ex" ) ) 
    {
      exclude_phenotype( m.get_string( "phe.ex" ) );
    }

  if ( m.has_field( "phe.req" ) ) 
    {
      require_phenotype( m.get_string( "phe.req" ) );
    }

  if ( m.has_field( "phe.obs" ) )
    {      
      // allow over-ride in the mask by phe.allow.missing (i.e. as phe.obs is automatically added if a --phenotype etc is specified)
      if ( ! m.has_field("phe.allow.missing") ) 
	require_nonmissing_phenotype( m.get_string( "phe.obs" ) );
    }

  
  // Options to skip certain types of input from the VARDB

  if ( m.has_field( "no-geno" ) )
    {
      load_genotype_data( false );      
      load_genotype_meta( false );
    }
  
  if ( m.has_field( "no-vmeta" ) )
    {
      load_variant_meta( false );
    }

  if ( m.has_field( "no-gmeta" ) )
    {
      load_genotype_meta( false );
    }

  if ( m.has_field( "load-gmeta" ) )
    {
      std::vector<std::string> t = m.get_string( "load-gmeta" );
      for (int i=0;i<t.size(); i++) set_load_genotype_meta( t[i] );      
    }

  if ( m.has_field( "load-vmeta" ) )
    {
      std::vector<std::string> t = m.get_string( "load-vmeta" );
      for (int i=0;i<t.size(); i++) set_load_variant_meta( t[i] );      
    }


  // TODO: Based on the command and the mask/options/arguments, determine what needs to be loaded and automatically set this.
  //       Will have to be careful/conservative in doing that...

  
}

bool Mask::eval(Variant & v, void * p)
{
  
  //
  // No filters set? 
  //

  if ( (!annot) 
       && filterFunctions.size() == 0  
       && req_filterFunctions.size() == 0 ) return true;
  

  //
  // Required filters (i.e require all must be true )
  //

  set<mask_func_t>::iterator f = req_filterFunctions.begin();
  while ( f != req_filterFunctions.end() )
    {
      if ( ! (*f)(v,p) ) return false;
      ++f;
    }


  //
  // Inclusion filters (i.e. include if at least one to be true )
  //

  bool include = filterFunctions.size() == 0 ? true : false;

  f = filterFunctions.begin();
  while ( f != filterFunctions.end() )
    {
      if ( (*f)(v,p) ) include = true;
      ++f;
    }

  if ( ! include ) return false;
  
  //
  // Any annotation filers?
  //

  if ( annot )
      {

	// This appends as meta-information the annotations to the variant
	Annotate::annotate(v);

	if ( ! f_include_annotation( v ) )
	  return false;
	
	if ( ! f_require_annotation( v ) )
	  return false;

	if ( ! f_exclude_annotation( v ) )
	  return false;
	
      }
  
  //
  // Passed all filters 
  //

  return true;

}


int Mask::include_loc( int x )
{  
  if ( locdb ) 
    {
      in_locset.insert(x); 
    }
  else return 0;
  return x;
}

int Mask::include_seg( int x ) 
{
  if ( segdb ) { in_segset.insert(x); return x; } 
  return 0;
}

void Mask::include_reg( const std::vector<std::string> & k )
{
  for (int i=0; i<k.size(); i++) 
    {
      bool okay = false;
      Region r( k[i] , okay );
      if ( okay ) include_reg( r );
    }
}

// variant IDs

void Mask::include_id( const std::vector<std::string> & r )
{
  for (int i=0; i<r.size(); i++) in_ids.insert(r[i]);
}


void Mask::require_id( const std::vector<std::string> & r )
{
  for (int i=0; i<r.size(); i++) req_ids.insert(r[i]);
}


void Mask::exclude_id( const std::vector<std::string> & r )
{
  for (int i=0; i<r.size(); i++) ex_ids.insert(r[i]);
}


// variant groups

int Mask::include_var( int x )      
{  
    if ( vardb ) 
      {
	in_varset.insert(x); 
      }
    else return 0;
    return x;
}


int Mask::include_loc( const std::string & n )
{
  if ( ! locdb ) return 0;
  int id = locdb->lookup_group_id( n );
  if ( id > 0 ) return include_loc(id);
  else return 0;
}

int Mask::include_seg( const std::string & n )
{
  if ( ! segdb ) return 0;
  int id = segdb->lookup_group_id( n );
  if ( id > 0 ) return include_seg(id);
  else return 0;
}

int Mask::include_var( const std::string & n )
{
  if ( ! vardb ) return 0;
  const bool DO_NOT_ADD = true;
  int id = vardb->add_set( n , "" , DO_NOT_ADD );
  if ( id > 0 ) return include_var(id);
  else return 0;
}


int Mask::include_ref( int x )
{
  if ( refdb ) 
    {
      in_refset.insert(x);
      append_ref(x);
    }
  else return 0;
  return x;
}

int Mask::include_ref( const std::string & n )
{
  if ( ! refdb ) return 0;
  int id = refdb->lookup_group_id( n );
  if ( id > 0 ) return include_ref( id );
  else return 0;
}



//////////////////
// Requires

int Mask::require_loc( int x )
{
  if ( locdb ) 
    {
      req_locset.insert(x); 
    }
  else return 0;
  return x;
}

int Mask::require_seg( int x )
{
  if ( segdb ) 
    {
      req_segset.insert(x); 
    }
  else return 0;
  return x;
}

int Mask::require_loc( const std::string & n )
{
  if ( ! locdb ) return 0;
  int id = locdb->lookup_group_id( n );
  if ( id > 0 ) return require_loc(id);
  else return 0;  
}

int Mask::require_seg( const std::string & n )
{
  if ( ! segdb ) return 0;
  int id = segdb->lookup_group_id( n );
  if ( id > 0 ) return require_seg(id);
  else return 0;  
}

int Mask::require_var( int x )
{
    if ( vardb ) 
      {
	req_varset.insert(x); 
      }
    else return 0;
    return x;
}

int Mask::require_var( const std::string & n )
{
  if ( ! vardb ) return 0;
  const bool DO_NOT_ADD = true;
  int id = vardb->add_set( n , "" , DO_NOT_ADD );
  if ( id > 0 ) return require_var(id);
  else return 0;  
}

int Mask::require_ref( int x )
{
  if ( refdb ) 
    {
      req_refset.insert(x);      
      append_ref(x);
    }
  else return 0;
  return x;
}

int Mask::require_ref( const std::string & n )
{
  if ( ! refdb ) return 0;
  int id = refdb->lookup_group_id( n );
  if ( id > 0 ) return require_ref( id );
  else return 0;
}


//////////////////
// Excludes

int Mask::exclude_loc( int x )
{  
    if ( locdb ) 
      {
	ex_locset.insert(x); 
      }
    else return 0;
    return x;
}

int Mask::exclude_seg( int x )
{  
    if ( segdb ) 
      {
	ex_segset.insert(x); 
      }
    else return 0;
    return x;
}

int Mask::exclude_var( int x )      
{  
  if ( vardb ) 
    {
      ex_varset.insert(x); 
    }
  else return 0;
  return x;
}


int Mask::exclude_loc( const std::string & n )
{
  if ( ! locdb ) return 0;
  int id = locdb->lookup_group_id( n );
  if ( id > 0 ) return exclude_loc(id);
  else return 0;
}

int Mask::exclude_seg( const std::string & n )
{
  if ( ! segdb ) return 0;
  int id = segdb->lookup_group_id( n );
  if ( id > 0 ) return exclude_seg(id);
  else return 0;
}

int Mask::exclude_var( const std::string & n )
{
  if ( ! vardb ) return 0;
  const bool DO_NOT_ADD = true;
  int id = vardb->add_set( n , "" , DO_NOT_ADD );
  if ( id > 0 ) return exclude_var(id);
  else return 0;
}


int Mask::exclude_ref( int x )
{
  if ( refdb ) 
    {
      ex_refset.insert(x);
      append_ref(x);
    }
  else return 0;
  return x;
}

int Mask::exclude_ref( const std::string & n )
{
  if ( ! refdb ) return 0;
  int id = refdb->lookup_group_id( n );
  if ( id > 0 ) return exclude_ref( id );
  else return 0;
}


//
// Subsetting/skipping for locus-inclusions
//

void Mask::subset_loc(const int s, const std::string & d)
{
  include_loc(s);
  if ( in_locset.find(s) == in_locset.end() ) return;
  std::map<int, std::set<std::string> >::iterator i = subset_locset.find( s );
  if ( i == subset_locset.end() ) 
    {
      std::set<string> t;
      t.insert( d );
      subset_locset.insert( make_pair( s , t ) );
    }
  else
    i->second.insert(d);
}

void Mask::skip_loc(const int s, const std::string & d)
{
  include_loc(s);
  if ( in_locset.find(s) == in_locset.end() ) return;
  map<int, set<string> >::iterator i = skip_locset.find( s );
  if ( i == skip_locset.end() ) 
    {
      set<string> t;
      t.insert( d );
      skip_locset.insert( make_pair( s , t ) );
    }
  else
    i->second.insert(d);
}

//
// Subsetting/skipping for locusset-inclusions
//

void Mask::subset_loc_set( const std::string & s1 , 
			   const std::string & s2 , 
			   const std::string & s3 )
{
  int gid = include_loc_set( s1, s2 );
  if ( gid ) subset_locset_set[ gid ].insert( s3 );
}

void Mask::skip_loc_set( const std::string & s1 , 
			 const std::string & s2 , 
			 const std::string & s3 )
{
  Helper::halt(" mask locset.skip not implemented yet ");
}
			   

bool Mask::insert_locset( const int j , const std::string & n ) const
{
  
  // if no subsets specified, this is fine.
  if ( subset_locset_set.size() == 0 ) return true;
  
  std::map<int, std::set<std::string> >::const_iterator i = subset_locset_set.find(j);

  if ( i == subset_locset_set.end() ) 
    Helper::halt( "internal prob in Mask::insert_locset()" );
  
  // If no subset specified, fine to accept all
  if ( i->second.size() == 0 ) return true;

  // Otherwise, T/F for whether this set is in the subset-list
  return i->second.find( n ) != i->second.end();
}



//
// Subsetting/skipping for variant-inclusions
//

void Mask::subset_var(const int s, const std::string & d)
{

  Helper::halt( "var.subset not currently supported" );

  include_var(s);
  if ( in_varset.find(s) == in_varset.end() ) return;
  map<int, set<string> >::iterator i = subset_varset.find( s );
  if ( i == subset_varset.end() ) 
    {
      set<string> t;
      t.insert( d );
      subset_varset.insert( make_pair( s , t ) );
    }
  else
    i->second.insert(d);
}


void Mask::skip_var(const int s, const string & d)
{
  Helper::halt( "var.skip not currently supported" );

  include_var(s);
  if ( in_varset.find(s) == in_varset.end() ) return;
  map<int, set<string> >::iterator i = skip_varset.find( s );
  if ( i == skip_varset.end() ) 
    {
      set<string> t;
      t.insert( d );
      skip_varset.insert( make_pair( s , t ) );
    }
  else
    i->second.insert(d);
}


void Mask::subset_loc(const int s, const std::vector<std::string> & d)
{
  for (int i=0;i<d.size(); i++) subset_loc(s,d[i]);
}

void Mask::subset_var(const int s, const std::vector<std::string> & d)
{
  for (int i=0;i<d.size(); i++) subset_var(s,d[i]);
}


void Mask::skip_loc(const int s, const std::vector<std::string> & d)
{
  for (int i=0;i<d.size(); i++) skip_loc(s,d[i]);
}

void Mask::skip_var(const int s, const std::vector<std::string> & d)
{
  for (int i=0;i<d.size(); i++) skip_var(s,d[i]);
}



//
// Specify individual genes based on altname
//

void Mask::subset_loc_altname(const std::string & group , const std::string & altname )
{
  std::vector<std::string> s = locdb->fetch_name_given_altname( group , altname );
  if ( s.size() > 0 ) 
    subset_loc( group , s );
}

void Mask::subset_loc_altname(const std::string & grp, const std::vector<std::string>& n)
{
  for (int i=0; i<n.size(); i++)
    subset_loc_altname( grp , n[i] );
}



//
// Allow string-idenitification of sets
//

void Mask::subset_loc(const string & n, const string& d)
{
    if ( ! locdb ) return;
    int id = locdb->lookup_group_id( n );
    if ( id > 0 ) return subset_loc(id,d);    
    return;
}

void Mask::subset_loc(const string & n, const vector<string>& d)
{
  if ( ! locdb ) return;
  int id = locdb->lookup_group_id( n );
  if ( id > 0 ) return subset_loc(id,d);
}

void Mask::skip_loc(const string & n, const string& d)
{
    if ( ! locdb ) return;
    int id = locdb->lookup_group_id( n );
    if ( id > 0 ) return skip_loc(id,d);
    else return;
}

void Mask::skip_loc(const string & n, const vector<string>& d)
{
    if ( ! locdb ) return;
    int id = locdb->lookup_group_id( n );
    if ( id > 0 ) return skip_loc(id,d);
    else return;
}


void Mask::subset_var(const string & n, const string& d)
{
  if ( ! vardb ) return;
  const bool DO_NOT_ADD = true;
  int id = vardb->add_set( n , "" , DO_NOT_ADD);
  if ( id > 0 ) return subset_var(id,d);
  else return;
}

void Mask::subset_var(const string & n, const vector<string>& d)
{
  if ( ! vardb ) return;
  const bool DO_NOT_ADD = true;
  int id = vardb->add_set( n , "" , DO_NOT_ADD );
  if ( id > 0 ) return subset_var(id,d);
  else return;
}

void Mask::skip_var(const string & n, const string& d)
{
  if ( ! vardb ) return;
  const bool DO_NOT_ADD = true;
  int id = vardb->add_set( n , "" , DO_NOT_ADD );
  if ( id > 0 ) return skip_var(id,d);
  else return;
}

void Mask::skip_var(const string & n, const vector<string>& d)
{
  if ( ! vardb ) return;
  const bool DO_NOT_ADD = true;
  int id = vardb->add_set( n , "" , DO_NOT_ADD );
  if ( id > 0 ) return skip_var(id,d);
  else return;
}




//
// File inclusion/exclusion
//


int Mask::include_file( const string & filetag )
{
  int id = vardb->file_tag( filetag );
  if ( id != 0 ) 
    {
      in_files.insert(id);
    }
  return in_files.size();
}


// int Mask::exclude_file( const int id )
// {
//   return 
// }


int Mask::exclude_file( const string & filetag )
{
  if ( ! vardb ) return 0;
  int id = vardb->file_tag( filetag );
  if ( id != 0 ) 
    {
      ex_files.insert(id); 
    }
  return ex_files.size();
}


//
// Appends from locus database
//

int Mask::append_loc( int x )
{  
  if ( locdb ) 
    {
      app_locset.insert(x); 
      //include_loc(x);
    }
  else return 0;
  return x;
}

int Mask::append_loc( const string & n )
{
  if ( ! locdb ) return 0;
  int id = locdb->lookup_group_id( n );
  if ( id > 0 ) return append_loc(id);
  else return 0;
}


//
// Appends from variant database
//

int Mask::append_var( int x )      
{  
    if ( vardb ) 
      {
	app_varset.insert(x); 
	// include_var(x);
      }
    else return 0;
    return x;
}

int Mask::append_var( const std::string & n )
{
  if ( ! vardb ) return 0;
  const bool DO_NOT_ADD = true;
  int id = vardb->add_set( n , "" , DO_NOT_ADD );
  if ( id > 0 ) return append_var(id);
  else return 0;
}

int Mask::append_var_set( int x )
{
  Helper::halt("not implemented yet");
  return 0;
}

int Mask::append_var_set( const std::string & n )
{
  Helper::halt("not implemented yet");
  return 0;
}


//
// Appends from reference database
//

int Mask::append_ref( int x )	
{  
  if ( refdb ) 
    {
      app_refset.insert(x); 
    }
  else return 0;
  return x;
}    

int Mask::append_ref( const string & n )
{
  if ( ! refdb ) return 0;
  int id = refdb->lookup_group_id( n );
  if ( id > 0 ) return append_ref(id);
  else return 0;
}


//
// Variant grouping 
//

void Mask::group_var(const int g) 
{ 
  group_variant = g; 
  group_variant_set = 0;
  group_locus = 0; 
  group_locus_set = 0;
  include_var(g);
}

void Mask::group_loc(const int g) 
{ 
  group_locus = g;
  group_variant_set = 0;
  group_locus_set = 0;
  group_variant = 0;
  include_loc(g);
}

void Mask::group_var(const string & g) 
{ 
  if ( ! vardb ) return;
  include_var(g); 
  const bool DO_NOT_ADD = true;
  int id = vardb->add_set( g , "" , DO_NOT_ADD );
  if ( id > 0 ) group_var(id);
  else Helper::halt("could not find var super-set " + g );
}

void Mask::group_var_set( const int g )
{
  group_variant = 0;
  group_variant_set = g;
  group_locus = 0;
  group_locus_set = 0;
}

void Mask::group_var_set( const std::string & g )
{
  if ( ! vardb )return;
  include_varset(g);
  const bool DO_NOT_ADD = true;
  int pid = vardb->add_superset( g , "" , DO_NOT_ADD );
  if ( pid > 0 ) return group_var_set(pid);  
  return;
}


// TODO TODO TODO 
// void Mask::group_loc_set(const string & n, const string & p)
// {  
//   if ( ! locdb ) return;
//   include_loc_set(n,p);
//   int pid = locdb->lookup_set_id( n,p );
//   if ( pid > 0 ) return group_loc_set(pid);  
//   return;
// }


void Mask::group_loc(const string & g) 
{ 
  if ( ! locdb ) return;
  include_loc(g);
  int id = locdb->lookup_group_id( g );
  if ( id > 0 ) group_loc(id);    
}

void Mask::group_reg( const std::vector<std::string> & g )
{
  // allow window specification
  // i.e.  chr1:10..20,chr1:20..30 is fine
  // but also 
  //    w:10:10
  // meaning a window of 10 bases, sliding over 10 bases \
  //  question: how to we know the limits? 
  //  answer: 
  
  Helper::halt("not implemented reg.group yet");
}

void Mask::include_annotation_nonsyn()
{
  annot = true;
  in_annotations.push_back("_MIS");
  in_annotations.push_back("_NON");
}

void Mask::require_annotation_nonsyn()
{
  annot = true;
  req_annotations.push_back("_MIS");
  req_annotations.push_back("_NON");
}

void Mask::exclude_annotation_nonsyn()
{
  annot = true;
  ex_annotations.push_back("_MIS");
  ex_annotations.push_back("_NON");
}


void Mask::include_annotation(const string & s)
{
  annot = true;
  in_annotations.push_back(s);
}

void Mask::require_annotation(const string & s)
{
  annot = true;
  req_annotations.push_back(s);
}

void Mask::exclude_annotation(const string & s)
{
  annot = true;
  ex_annotations.push_back(s);
}

void Mask::include_annotation( const vector<string> & s)
{
  for (int i=0; i<s.size(); i++) include_annotation( s[i] );
}

void Mask::require_annotation( const vector<string> & s)
{
  for (int i=0; i<s.size(); i++) require_annotation( s[i] );
}

void Mask::exclude_annotation( const vector<string> & s)
{
  for (int i=0; i<s.size(); i++) exclude_annotation( s[i] );    
}

bool Mask::f_include_annotation(const Variant & v)
{
  if ( in_annotations.size() == 0 ) return true;
  for (int i=0; i<in_annotations.size(); i++)
    if ( v.meta.has_field( in_annotations[i] ) ) return true;
  return false;
}

bool Mask::f_require_annotation(const Variant & v)
{
  if ( req_annotations.size() == 0 ) return true;
  for (int i=0; i<req_annotations.size(); i++)
    if ( !v.meta.has_field( req_annotations[i] ) ) return false;
  return true;
}

bool Mask::f_exclude_annotation(const Variant & v)
{
  for (int i=0; i<ex_annotations.size(); i++)
    if ( v.meta.has_field( ex_annotations[i] ) ) return false;
  return true;  
}

//
// Functions for variant super-sets
//

int Mask::include_varset( const std::string & g )
{
  // simply add each implied subset individually to the include
  std::vector<std::string> sets = GP->vardb.get_sets( g );
  for (int s=0;s<sets.size();s++) include_var( sets[s] );   
  return 0;
}

int Mask::include_varset( int x )
{
  Helper::halt("Not implemented yet.");
  return 0;
}


int Mask::require_varset( const std::string & g )
{
  // simply add each implied subset individually to the include
  std::vector<std::string> sets = GP->vardb.get_sets( g );
  for (int s=0;s<sets.size();s++) require_var( sets[s] );
  return 0;  
}

int Mask::exclude_varset( const std::string & g )
{
  // simply add each implied subset individually to the include
  std::vector<std::string> sets = GP->vardb.get_sets( g );
  for (int s=0;s<sets.size();s++) exclude_var( sets[s] );   
  return 0;  
}


//
// Functions for locus-sets
//

int Mask::include_loc_set( int x )
{
  if ( locdb ) 
    {
      in_locset_set.insert( x ); 
    }
  else return 0;
  return x;
}

int Mask::include_loc_set( std::string n , std::string p )
{
  int pid = locdb ? locdb->lookup_set_id( n, p ) : 0 ;
  return  pid > 0 ? include_loc_set(pid) : 0 ;    
}

int Mask::append_loc_set( int x )
{
  if ( locdb ) 
    {
      app_locset_set.insert(x); 
    }
  else return 0;
  return x;
}

int Mask::append_loc_set( const string & n , const string & p )
{
  if ( ! locdb ) return 0;
  int pid = locdb->lookup_set_id( n,p );
  if ( pid > 0 ) return append_loc(pid);
  else return 0;

}

void Mask::group_loc_set(const int g)
{
  group_locus_set = g;
  group_locus = 0;
  group_variant = 0;
  group_variant_set = 0;
  include_loc_set(g);
}

void Mask::group_loc_set(const string & n, const string & p)
{  
  if ( ! locdb ) return;
  include_loc_set(n,p);
  int pid = locdb->lookup_set_id( n,p );
  if ( pid > 0 ) return group_loc_set(pid);  
  return;
}


int Mask::require_loc_set( int x )
{
  if ( locdb ) 
    {
      req_locset_set.insert(x); 
    }
  else return 0;
  return x;
}

int Mask::require_loc_set( string n , string p )
{
  if ( ! locdb ) return 0;
  int pid = locdb->lookup_set_id( n, p );
  if ( pid > 0 ) return require_loc_set(pid);
  else return 0;

}


int Mask::exclude_loc_set( int x )
{
  if ( locdb ) 
    {
      ex_locset_set.insert(x); 
    }
  else return 0;
  return x;
}

int Mask::exclude_loc_set( string n , string p )
{
  if ( ! locdb ) return 0;
  int pid = locdb->lookup_set_id( n, p );
  if ( pid > 0 ) return exclude_loc_set(pid);
  else return 0;

}


// filters

void Mask::include_filter( const string & s )
{
  inc_filter.insert(s);
}

void Mask::require_filter( const string & s )
{
  req_filter.insert(s);
}

void Mask::exclude_filter( const string & s )
{
  exc_filter.insert(s);
}

void Mask::include_filter( const vector<string> & s )
{
  for (int i=0; i<s.size(); i++) 
    include_filter( s[i] );
}

void Mask::require_filter( const vector<string> & s )
{
  for (int i=0; i<s.size(); i++) 
    require_filter( s[i] );
}

void Mask::exclude_filter( const vector<string> & s )
{
  for (int i=0; i<s.size(); i++) 
    exclude_filter( s[i] );
}

void Mask::include_filter_any()
{
  inc_filter_any = true;
}

void Mask::require_filter_any()
{
  req_filter_any = true;
}

void Mask::exclude_filter_any()
{
  exc_filter_any = true;
}


bool Mask::eval_filters( const SampleVariant & v )
{
  
  // QUAL filter?

  if ( use_qual_filter && ! qual.in( v.quality() ) ) return false;

  // Needs at least 1 inclusion filter, if 1+ includes
  // Needs all required filters, if 1+ requires
  // Needs to not have any excludes

  bool include = inc_filter.size() == 0 ? true : false;
  
  // Filter of *any* elements?
  
  if ( exc_filter_any && v.any_filter() ) return false;

  if ( ( inc_filter_any || req_filter_any ) && ! v.any_filter() ) return false; 

  // Excludes
    
  set<string>::iterator i = exc_filter.begin();
  while ( i != exc_filter.end() )
    {
      std::vector< std::string > keys = v.meta_filter.keys();
      if ( v.meta_filter.has_field( *i ) ) 
	{
	  return false;
	}
      ++i;
    }

  // Requires
  
  i = req_filter.begin();
  while ( i != req_filter.end() )
    {
      if ( ! v.meta_filter.has_field(*i) ) 
	return false;
      ++i;
    }
  
  if ( req_filter.size() > 0 ) 
    {
      // If here and true, then we must have ssen all the necessary 
      // requires -- therefore, accept, we do not need to further insist 
      // on any includes
      return true;
    }
  

  // Includes 

  i = inc_filter.begin();
  while ( i != inc_filter.end() )
    {
      if ( v.meta_filter.has_field( *i ) ) 
	{
	  include = true;
	  break;
	}
      ++i;
    }
  
  return include;
   
}


bool Mask::polymorphism_filter( const Variant & v )
{
  // also include type of polymorphism here (indels vs snps, initially)
  if ( req_indels || exc_indels )
    {
      bool is_indel = v.indel();
      if ( req_indels && ! is_indel ) return false;
      if ( exc_indels &&   is_indel ) return false;
    }

  if ( req_novel && v.name() != "." ) return false;
  if ( exc_novel && v.name() == "." ) return false;

  if ( require_biallelic() && ! v.biallelic() ) return false;
  if ( exclude_biallelic() &&   v.biallelic() ) return false;
  if ( require_monomorphic() && ! v.monomorphic() ) return false;
  if ( exclude_monomorphic() &&   v.monomorphic() ) return false;
  
  return true;
}


//
// Individuals
//

void Mask::include_indiv( const std::string & id )
{
  in_indset.insert(id);
}

void Mask::include_indiv( const std::vector<std::string> & id )
{
  for (int i=0; i<id.size(); i++ ) include_indiv( id[i] );
}
  
void Mask::exclude_indiv( const std::string & id )
{
  ex_indset.insert(id);
}

void Mask::exclude_indiv( const std::vector<std::string> & id )
{
  for (int i=0; i<id.size(); i++ ) exclude_indiv( id[i] );
}
  
bool Mask::use_indiv( const std::string & id ) const
{
  if ( in_indset.size() > 0 && in_indset.find(id) == in_indset.end() ) return false;
  if ( ex_indset.find(id) != ex_indset.end() ) return false;
  return true;    
}

bool Mask::use_file( const int f ) const
{
  if ( in_files.size() > 0 && in_files.find(f) == in_files.end() ) return false;
  if ( ex_files.find(f) != ex_files.end() ) return false;
  return true;    
}

// expect in format phe=scz:2,dis:1 

void Mask::include_phenotype( const std::vector< std::string > & val )
{
  for (int i=0;i<val.size();i++)
    {
      std::vector< std::string > k = Helper::parse( val[i] , ":" );
      for (int j=1;j<k.size(); j++)
	in_phe[ k[0] ].insert( k[j] );
    }
}

void Mask::require_phenotype( const std::vector< std::string > & val )
{
  for (int i=0;i<val.size();i++)
    {
      std::vector< std::string > k = Helper::parse( val[i] , ":" );
      for (int j=1;j<k.size(); j++)
	req_phe[ k[0] ].insert( k[j] );
    }
}

void Mask::exclude_phenotype( const std::vector< std::string > & val )
{
  for (int i=0;i<val.size();i++)
    {
      std::vector< std::string > k = Helper::parse( val[i] , ":" );
      for (int j=1;j<k.size(); j++)
	ex_phe[ k[0] ].insert( k[j] ) ;
    }
}

void Mask::require_nonmissing_phenotype( const std::vector< std::string> & k)
{
  for (int i=0; i<k.size(); i++)
    req_nonmissing_phenotype.insert( k[i] );
}

bool Mask::pheno_screen( Individual * person ) const
{

  if ( req_nonmissing_phenotype.size() > 0 ) 
    {
      std::set<std::string>::iterator i = req_nonmissing_phenotype.begin();
      while ( i != req_nonmissing_phenotype.end() )
	{
	  if ( ! person->meta.has_field( *i ) ) return false;
	  ++i;
	}
    }

  std::map< std::string, std::set<std::string> >::const_iterator i = ex_phe.begin();
  while ( i != ex_phe.end() )
    {
      std::set<std::string>::const_iterator j = i->second.begin();
      while ( j != i->second.end() )
	{
	  if ( person->meta.as_string( i->first ) == *j ) return false;
	  ++j;
	}
      ++i;
    }
  
  i = req_phe.begin();
  while ( i != req_phe.end() )
    {
      // for each distinct phenotype, need at least value seen.
      bool observed = false;
      std::set<std::string>::const_iterator j = i->second.begin();
      while ( j != i->second.end() )
	{
	  if ( person->meta.as_string( i->first ) == *j ) observed = true;
	  ++j;
	}
      if ( ! observed ) return false;
      ++i;
    }

  if ( in_phe.size() == 0 ) return true;

  bool included = false;
  i = in_phe.begin();
  while ( i != in_phe.end() )
    {
      std::set<std::string>::const_iterator j = i->second.begin();
      while ( j != i->second.end() )
	{	  
	  if ( person->meta.as_string( i->first ) == *j ) included = true;      
	  ++j;
	}
      ++i;
    }
  return included;

}


//
// Variant meta-information filters
//

int Mask::meta_equals( const std::string & key , int value )
{
  meta_eq[key] = value;
}

int Mask::meta_not_equals( const std::string & key , int value )
{
  meta_ne[key] = value;
}

int Mask::meta_equals( const std::string & key , const std::string & value )
{
  meta_has_text[key].insert( value );
}

int Mask::meta_not_equals( const std::string & key , const std::string & value )
{
  meta_has_not_text[key].insert( value );
}

int Mask::meta_greater( const std::string & key , double value )
{
  meta_gt[key] = value;
}

int Mask::meta_greater_equal( const std::string & key , double value )
{
  meta_ge[key] = value;
}

int Mask::meta_less( const std::string & key , double value )
{
  meta_lt[key] = value;
}

int Mask::meta_less_equal( const std::string & key , double value )
{
  meta_le[key] = value;
}



int Mask::req_meta_equals( const std::string & key , int value )
{
  req_meta_eq[key] = value;
}

int Mask::req_meta_not_equals( const std::string & key , int value )
{
  req_meta_ne[key] = value;
}

int Mask::req_meta_equals( const std::string & key , const std::string & value )
{
  req_meta_has_text[key].insert( value );
}

int Mask::req_meta_not_equals( const std::string & key , const std::string & value )
{
  req_meta_has_not_text[key].insert( value );
}

int Mask::req_meta_greater( const std::string & key , double value )
{
  req_meta_gt[key] = value;
}

int Mask::req_meta_greater_equal( const std::string & key , double value )
{
  req_meta_ge[key] = value;
}

int Mask::req_meta_less( const std::string & key , double value )
{
  req_meta_lt[key] = value;
}

int Mask::req_meta_less_equal( const std::string & key , double value )
{
  req_meta_le[key] = value;
}




bool Mask::eval( SampleVariant & svar )
{


  MetaInformation<VarMeta> & m = svar.meta;


  // Other meta=AB:ge:0.75 stype filters:

  // Set a genotype to missing if it fails a filter, or does not have
  // a given meta-field If a text match is "", this means we just need
  // to have the key


  //
  // Evaluate requirements first, then includes
  //
  
  if ( meta_requires() )
    {
      
      std::map<string,int>::const_iterator i = req_meta_eq.begin();
      while ( i != req_meta_eq.end() )
	{
	  if ( m.has_field( i->first ) )
	    {
	      
	      mType mt = MetaInformation<VarMeta>::type( i->first );
	      if ( mt == META_INT ) 
		{
		  if ( m.get1_int( i->first ) != i->second ) return false;
		}
	      else 
		return false;
	    }
	  else 
	    return false;
	  ++i;
	}
      
      i = req_meta_ne.begin();
      while ( i != req_meta_ne.end() )
	{
	  if ( m.has_field( i->first ) )
	    {
	      mType mt = MetaInformation<VarMeta>::type( i->first );
	      if ( mt == META_INT ) 
		{
		  if ( m.get1_int( i->first ) == i->second ) return false;
		}
	      else 
		return false;	  
	    }
	  else 
	    return false;
	  ++i;
	}
      
      
      //
      // Note -- for inequalities, we need to allow for both int and floats..
      //
      
      // Greater than
      
      std::map<string,double>::const_iterator j = req_meta_gt.begin();
      while ( j != req_meta_gt.end() )
	{
	  if ( m.has_field( j->first ) )
	    {
	      mType mt = MetaInformation<VarMeta>::type( j->first );
	      if ( mt == META_INT ) 
		{
		  if ( m.get1_int( j->first ) <= j->second ) return false;
		}
	      else if ( mt == META_FLOAT )
 		{
		  if ( m.get1_double( j->first ) <= j->second ) return false;
		}
	    }
	  else 
	    return false;
	  ++j;
	}
      
      
      // Greater than or equal to 
      
      j = req_meta_ge.begin();
      while ( j != req_meta_ge.end() )
	{
	  if ( m.has_field( j->first ) )
	    {
	      mType mt = MetaInformation<VarMeta>::type( j->first );
	      if ( mt == META_INT ) 
		{
		  if ( m.get1_int( j->first ) < j->second ) return false;
		}
	      else if ( mt == META_FLOAT )
		{
		  if ( m.get1_double( j->first ) < j->second ) return false;
		}
	      
	    }
	  else 
	    return false;
	  ++j;
	}
      
      // Less than
      
      j = req_meta_lt.begin();
      while ( j != req_meta_lt.end() )
	{
	  if ( m.has_field( j->first ) )
	    {	  
	      mType mt = MetaInformation<VarMeta>::type( j->first );
	      if ( mt == META_INT ) 
		{
		  if ( m.get1_int( j->first ) >= j->second ) return false;
		}
	      else if ( mt == META_FLOAT )
		{
		  if ( m.get1_double( j->first ) >= j->second ) return false;
		}
	    }
	  else 
	    return false;
	  ++j;
	}
      
      // Less than or equal to
      
      j = req_meta_le.begin();
      while ( j != req_meta_le.end() )
	{
	  if ( m.has_field( j->first ) )
	    {
	      mType mt = MetaInformation<VarMeta>::type( j->first );
	      if ( mt == META_INT ) 
		{
		  if ( m.get1_int( j->first ) > j->second ) return false;
		}
	      else if ( mt == META_FLOAT )
		{
		  if ( m.get1_double( j->first ) > j->second ) return false;
		}
	    }
	  else 
	    return false;
	  ++j;
	}
      
      // Text matching (equals/set)

      std::map<string,set<string> >::const_iterator k = req_meta_has_text.begin();
      while ( k != req_meta_has_text.end() )
	{
	  
	  bool okay = false;
	  
	  if ( m.has_field( k->first ) )
	    {	      

	      if ( MetaInformation<VarMeta>::type( k->first ) == META_FLAG ) 
		okay = true;
	      else
		{
		  std::vector<std::string> t = m.get_string( k->first );	      
		  for (int z=0; z<t.size(); z++)
		    if ( k->second.find( t[z] ) != k->second.end() ) 
		      okay=true;	      
		}
	    }
	  else 
	    return false;
	  
	  if ( ! okay ) return false;
	  
	  ++k;
	}
      
      // Text matching, not equals/set

      k = req_meta_has_not_text.begin();
      while ( k != req_meta_has_not_text.end() )
	{
	  if ( m.has_field( k->first ) )
	    {
	      if ( MetaInformation<VarMeta>::type( k->first ) == META_FLAG ) 
		return false;

	      if ( k->second.size() > 0 || 
		   k->second.find( m.get1_string( k->first ) ) != k->second.end() ) 
		return false;
	    }
	  else

	  ++k;
	}
    
    }
  
  
  //
  // This genotype has passed all requirements. Now evalaute any
  // includes-- fail if not 1+ meet
  //

  while ( meta_includes() ) 
    {

      // break when we meet our first condition == T

      bool okay = false;

      std::map<string,int>::const_iterator i = meta_eq.begin();
      while ( i != meta_eq.end() )
	{
	  if ( m.has_field( i->first ) )
	    {	      
	      mType mt = MetaInformation<VarMeta>::type( i->first );
	      if ( mt == META_INT ) 
		if ( m.get1_int( i->first ) == i->second ) { okay=true; break; }
	    }

	  ++i;
	}
      if ( okay ) break;
      
      
      i = meta_ne.begin();
      while ( i != meta_ne.end() )
	{
	  if ( m.has_field( i->first ) )
	    {
	      mType mt = MetaInformation<VarMeta>::type( i->first );
	      if ( mt == META_INT ) 
		if ( m.get1_int( i->first ) != i->second ) { okay=true; break; }
	    }
	  ++i;
	}
      if ( okay ) break;
      
      //
      // Note -- for inequalities, we need to allow for both int and floats..
      //
      
      // Greater than
      
      std::map<string,double>::const_iterator j = meta_gt.begin();
      while ( j != meta_gt.end() )
	{
	  if ( m.has_field( j->first ) )
	    {
	      mType mt = MetaInformation<VarMeta>::type( j->first );
	      if ( mt == META_INT ) 
		{
		  if ( m.get1_int( j->first ) > j->second ) { okay=true; break; }
		}
	      else if ( mt == META_FLOAT )
		{
		  if ( m.get1_double( j->first ) > j->second )  { okay=true; break; }
		}
	    }
	  ++j;
	}
      if ( okay ) break;
      
      
      // Greater than or equal to 
      
      j = meta_ge.begin();
      while ( j != meta_ge.end() )
	{
	  if ( m.has_field( j->first ) )
	    {
	      mType mt = MetaInformation<VarMeta>::type( j->first );
	      if ( mt == META_INT ) 
		{
		  if ( m.get1_int( j->first ) >= j->second ) { okay=true; break; }
		}
	      else if ( mt == META_FLOAT )
		{
		  if ( m.get1_double( j->first ) >= j->second ) { okay=true; break; }
		}
	      
	    }
	  ++j;
	}
      if ( okay ) break;


      // Less than
      
      j = meta_lt.begin();
      while ( j != meta_lt.end() )
	{
	  if ( m.has_field( j->first ) )
	    {	  
	      mType mt = MetaInformation<VarMeta>::type( j->first );
	      if ( mt == META_INT ) 
		{
		  if ( m.get1_int( j->first ) < j->second ) { okay=true; break; }
		}
	      else if ( mt == META_FLOAT )
		{
		  if ( m.get1_double( j->first ) < j->second ) { okay=true; break; }
		}
	    }
	  ++j;
	}
      if ( okay ) break;


      // Less than or equal to
      
      j = meta_le.begin();
      while ( j != meta_le.end() )
	{
	  if ( m.has_field( j->first ) )
	    {
	      mType mt = MetaInformation<VarMeta>::type( j->first );
	      if ( mt == META_INT ) 
		{
		  if ( m.get1_int( j->first ) <= j->second ) { okay=true; break; }
		}
	      else if ( mt == META_FLOAT )
		{
		  if ( m.get1_double( j->first ) <= j->second ) { okay=true; break; }
		}
	    }
	  ++j;
	}
      if ( okay ) break;


      std::map<string,set<string> >::const_iterator k = meta_has_text.begin();
      while ( k != meta_has_text.end() )
	{
	  if ( m.has_field( k->first ) )
	    {
	      if ( MetaInformation<VarMeta>::type( k->first ) == META_FLAG ) 
		{ okay=true; break; }

	      if (k->second.size() > 0 && 
		  k->second.find( m.get1_string( k->first )) != k->second.end() ) 
		{ okay = true; break; }
	    }
	  ++k;
	}
      if ( okay ) break;


      k = meta_has_not_text.begin();
      while ( k != meta_has_not_text.end() )
	{
	  if ( m.has_field( k->first ) )
	    {
	      if (k->second.size() > 0 || 
		  k->second.find( m.get1_string( k->first ) ) == k->second.end()  ) 
		{ okay=true; break; }
	    }
	  else
	    {
	      if ( MetaInformation<VarMeta>::type( k->first ) == META_FLAG ) 
		{ okay=true; break; }
	    }
	  ++k;
	}
      if ( okay ) break;

      return false;
    }
  
  
  // Accept this genotype

  return true;

}



//
// Genotype filters
//

int Mask::geno_equals( const std::string & key , int value )
{
  geno_eq[key] = value;
}

int Mask::geno_not_equals( const std::string & key , int value )
{
  geno_ne[key] = value;
}

int Mask::geno_equals( const std::string & key , const std::string & value )
{
  geno_has_text[key] = value;
}

int Mask::geno_not_equals( const std::string & key , const std::string & value )
{
  geno_has_not_text[key] = value;
}

int Mask::geno_greater( const std::string & key , double value )
{
  geno_gt[key] = value;
}

int Mask::geno_greater_equal( const std::string & key , double value )
{
  geno_ge[key] = value;
}

int Mask::geno_less( const std::string & key , double value )
{
  geno_lt[key] = value;
}

int Mask::geno_less_equal( const std::string & key , double value )
{
  geno_le[key] = value;
}



int Mask::req_geno_equals( const std::string & key , int value )
{
  req_geno_eq[key] = value;
}

int Mask::req_geno_not_equals( const std::string & key , int value )
{
  req_geno_ne[key] = value;
}

int Mask::req_geno_equals( const std::string & key , const std::string & value )
{
  req_geno_has_text[key] = value;
}

int Mask::req_geno_not_equals( const std::string & key , const std::string & value )
{
  req_geno_has_not_text[key] = value;
}

int Mask::req_geno_greater( const std::string & key , double value )
{
  req_geno_gt[key] = value;
}

int Mask::req_geno_greater_equal( const std::string & key , double value )
{
  req_geno_ge[key] = value;
}

int Mask::req_geno_less( const std::string & key , double value )
{
  req_geno_lt[key] = value;
}

int Mask::req_geno_less_equal( const std::string & key , double value )
{
  req_geno_le[key] = value;
}



bool Mask::eval( const Genotype & g ) const
{

  // Set a genotype to missing if it fails a filter, or does not have
  // a given meta-field If a text match is "", this means we just need
  // to have the key


  //
  // Evaluate requirements first, then includes
  //
  
  if ( geno_requires() )
    {

      std::map<string,int>::const_iterator i = req_geno_eq.begin();
      while ( i != req_geno_eq.end() )
	{
	  if ( g.meta.has_field( i->first ) )
	    {
	      
	      mType mt = MetaInformation<GenMeta>::type( i->first );
	      if ( mt == META_INT ) 
		{
		  if ( g.meta.get1_int( i->first ) != i->second ) return false;
		}
	      else 
		return false;
	    }
	  else 
	    return false;
	  ++i;
	}
      
      i = req_geno_ne.begin();
      while ( i != req_geno_ne.end() )
	{
	  if ( g.meta.has_field( i->first ) )
	    {
	      mType mt = MetaInformation<GenMeta>::type( i->first );
	      if ( mt == META_INT ) 
		{
		  if ( g.meta.get1_int( i->first ) == i->second ) return false;
		}
	      else 
		return false;	  
	    }
	  else 
	    return false;
	  ++i;
	}
      
      
      //
      // Note -- for inequalities, we need to allow for both int and floats..
      //
      
      // Greater than
      
      std::map<string,double>::const_iterator j = req_geno_gt.begin();
      while ( j != req_geno_gt.end() )
	{
	  if ( g.meta.has_field( j->first ) )
	    {
	      mType mt = MetaInformation<GenMeta>::type( j->first );
	      if ( mt == META_INT ) 
		{
		  if ( g.meta.get1_int( j->first ) <= j->second ) return false;
		}
	      else if ( mt == META_FLOAT )
		{
		  if ( g.meta.get1_double( j->first ) <= j->second ) return false;
		}
	    }
	  else 
	    return false;
	  ++j;
	}
      
      
      // Greater than or equal to 
      
      j = req_geno_ge.begin();
      while ( j != req_geno_ge.end() )
	{
	  if ( g.meta.has_field( j->first ) )
	    {
	      mType mt = MetaInformation<GenMeta>::type( j->first );
	      if ( mt == META_INT ) 
		{
		  if ( g.meta.get1_int( j->first ) < j->second ) return false;
		}
	      else if ( mt == META_FLOAT )
		{
		  if ( g.meta.get1_double( j->first ) < j->second ) return false;
		}
	      
	    }
	  else 
	    return false;
	  ++j;
	}
      
      // Less than
      
      j = req_geno_lt.begin();
      while ( j != req_geno_lt.end() )
	{
	  if ( g.meta.has_field( j->first ) )
	    {	  
	      mType mt = MetaInformation<GenMeta>::type( j->first );
	      if ( mt == META_INT ) 
		{
		  if ( g.meta.get1_int( j->first ) >= j->second ) return false;
		}
	      else if ( mt == META_FLOAT )
		{
		  if ( g.meta.get1_double( j->first ) >= j->second ) return false;
		}
	    }
	  else 
	    return false;
	  ++j;
	}
      
      // Less than or equal to
      
      j = req_geno_le.begin();
      while ( j != req_geno_le.end() )
	{
	  if ( g.meta.has_field( j->first ) )
	    {
	      mType mt = MetaInformation<GenMeta>::type( j->first );
	      if ( mt == META_INT ) 
		{
		  if ( g.meta.get1_int( j->first ) > j->second ) return false;
		}
	      else if ( mt == META_FLOAT )
		{
		  if ( g.meta.get1_double( j->first ) > j->second ) return false;
		}
	    }
	  else 
	    return false;
	  ++j;
	}
      

      std::map<string,string>::const_iterator k = req_geno_has_text.begin();
      while ( k != req_geno_has_text.end() )
	{
	  if ( g.meta.has_field( k->first ) )
	    {
	      if (k->second != "" && 
		  g.meta.get1_string( k->first ) != k->second ) 
		return false;
	    }
	  else 
	    return false;
	  ++k;
	}
      
      k = req_geno_has_not_text.begin();
      while ( k != req_geno_has_not_text.end() )
	{
	  if ( g.meta.has_field( k->first ) )
	    {
	      if (k->second == "" || 
		  g.meta.get1_string( k->first ) == k->second ) 
		return false;
	    }
	  ++k;
	}
    
    }
      
  
  //
  // This genotype has passed all requirements. Now evalaute any includes-- fail if not 1+ meet
  //

  while ( geno_includes() ) 
    {
      // break when we meet our first condition == T
      
      bool okay = false;

      std::map<string,int>::const_iterator i = geno_eq.begin();
      while ( i != geno_eq.end() )
	{
	  if ( g.meta.has_field( i->first ) )
	    {	      
	      mType mt = MetaInformation<GenMeta>::type( i->first );
	      if ( mt == META_INT ) 
		if ( g.meta.get1_int( i->first ) == i->second ) { okay = true; break; }
	    }

	  ++i;
	}
      if ( okay ) break;
      
      i = geno_ne.begin();
      while ( i != geno_ne.end() )
	{
	  if ( g.meta.has_field( i->first ) )
	    {
	      mType mt = MetaInformation<GenMeta>::type( i->first );
	      if ( mt == META_INT ) 
		if ( g.meta.get1_int( i->first ) != i->second ) { okay = true; break; }
	    }
	  ++i;
	}
      if ( okay ) break;

      
      //
      // Note -- for inequalities, we need to allow for both int and floats..
      //
      
      // Greater than
      
      std::map<string,double>::const_iterator j = geno_gt.begin();
      while ( j != geno_gt.end() )
	{
	  if ( g.meta.has_field( j->first ) )
	    {
	      mType mt = MetaInformation<GenMeta>::type( j->first );
	      if ( mt == META_INT ) 
		{
		  if ( g.meta.get1_int( j->first ) > j->second ) { okay = true; break; }
		}
	      else if ( mt == META_FLOAT )
		{
		  if ( g.meta.get1_double( j->first ) > j->second ) { okay = true; break; }
		}
	    }
	  ++j;
	}
      if ( okay ) break;

      
      // Greater than or equal to 
      
      j = geno_ge.begin();
      while ( j != geno_ge.end() )
	{
	  if ( g.meta.has_field( j->first ) )
	    {
	      mType mt = MetaInformation<GenMeta>::type( j->first );
	      if ( mt == META_INT ) 
		{
		  if ( g.meta.get1_int( j->first ) >= j->second ) { okay = true; break; }
		}
	      else if ( mt == META_FLOAT )
		{
		  if ( g.meta.get1_double( j->first ) >= j->second ) { okay = true; break; }
		}
	      
	    }
	  ++j;
	}
      
      if ( okay ) break;


      // Less than
      
      j = geno_lt.begin();
      while ( j != geno_lt.end() )
	{
	  if ( g.meta.has_field( j->first ) )
	    {	  
	      mType mt = MetaInformation<GenMeta>::type( j->first );
	      if ( mt == META_INT ) 
		{
		  if ( g.meta.get1_int( j->first ) < j->second ) { okay = true; break; }
		}
	      else if ( mt == META_FLOAT )
		{
		  if ( g.meta.get1_double( j->first ) < j->second ) { okay = true; break; }
		}
	    }
	  ++j;
	}
      
      if ( okay ) break;


      // Less than or equal to
      
      j = geno_le.begin();
      while ( j != geno_le.end() )
	{
	  if ( g.meta.has_field( j->first ) )
	    {
	      mType mt = MetaInformation<GenMeta>::type( j->first );
	      if ( mt == META_INT ) 
		{
		  if ( g.meta.get1_int( j->first ) <= j->second ) { okay = true; break; }
		}
	      else if ( mt == META_FLOAT )
		{
		  if ( g.meta.get1_double( j->first ) <= j->second ) { okay = true; break; }
		}
	    }
	  ++j;
	}
      
      if ( okay ) break;
      
      
      std::map<string,string>::const_iterator k = geno_has_text.begin();
      while ( k != geno_has_text.end() )
	{
	  if ( g.meta.has_field( k->first ) )
	    {
	      if (k->second != "" && 
		  g.meta.get1_string( k->first ) == k->second ) 
		{ okay = true; break; }
	    }
	  ++k;
	}
      
      if ( okay ) break;


      k = geno_has_not_text.begin();
      while ( k != geno_has_not_text.end() )
	{
	  if ( g.meta.has_field( k->first ) )
	    {
	      if (k->second == "" || 
		  g.meta.get1_string( k->first ) != k->second ) 
		{ okay = true; break; }
	    }
	  ++k;
	}
      
      if ( okay ) break;

      return false;
    }
  
  
  // Accept this genotype

  return true;

}


// Variant/genotype filters

void Mask::null_filter( const int_range & r )
{
  has_null_filter = true;
  null_fltr = r;
}

bool Mask::null_filter( ) const 
{
  return has_null_filter;
}

bool Mask::null_filter( const int j ) const
{
  return null_fltr.in( j );
}


void Mask::null_prop_filter( const dbl_range & r )
{
  has_null_prop_filter = true;
  null_prop_fltr = r;
}

bool Mask::null_prop_filter( ) const 
{
  return has_null_prop_filter;
}

bool Mask::null_prop_filter( const double j ) const
{
  return null_prop_fltr.in( j );
}


void Mask::case_control_filter( const std::string & a , const std::string & u )
{
  has_case_control_filter = true;
  case_fltr.set(a);
  control_fltr.set(u);
}

bool Mask::case_control_filter( ) const
{
  return has_case_control_filter;
}

bool Mask::case_control_filter( const int a, const int u) const
{
  return case_fltr.in(a) && control_fltr.in(u);
}


void Mask::set_filter_expression( const std::string & e )
{
  if ( ! eval_expr.parse(e) ) Helper::halt("could not parse expression: " + e + "\n" + eval_expr.errmsg() );
  eval_expr_set = true;
}

bool Mask::filter_expression_requires_genotypes() const
{
  return eval_expr.requires_genotypes();
}

bool Mask::calc_filter_expression( Variant & var , SampleVariant & svar )
{
  eval_expr.bind( svar );
  eval_expr.evaluate();
  bool passed = false;
  bool valid = eval_expr.value( passed );
  if ( ! valid ) 
    {
      plog.warn( "could not evaluate filter: " , var.coordinate() + ", " + eval_expr.errmsg() );
      return false; //invalid expressions always evaluate as false
    }
  if ( ! eval_filter_includes ) passed = ! passed;
  return passed;
}

bool Mask::calc_filter_expression( Variant & var , SampleVariant & svar , SampleVariant & gvar )
{
  eval_expr.bind( svar , gvar );
  eval_expr.evaluate();
  bool passed = false;
  bool valid = eval_expr.value( passed );
  if ( ! valid ) 
    {
      plog.warn( "could not evaluate filter: " , var.coordinate() + ", " + eval_expr.errmsg() );
      return false; //invalid expressions always evaluate as false
    }
  if ( ! eval_filter_includes ) passed = ! passed;
  return passed;
}



void Mask::var_set_filter_expression(const std::string & e )
{
  if ( ! var_eval_expr.parse(e) ) Helper::halt("could not set filter expression: " + e + "\n" + var_eval_expr.errmsg() );
  var_eval_expr_set = true;
}

bool Mask::var_filter_expression_requires_genotypes() const
{
  // not used
  return var_eval_expr.requires_genotypes();
}

bool Mask::var_calc_filter_expression( Variant & var )
{
  var_eval_expr.bind( var );
  var_eval_expr.evaluate();
  bool passed = false;
  bool valid = var_eval_expr.value(passed);
  if ( ! valid ) Helper::halt( "could not evaluate filter: " + var_eval_expr.errmsg() );
  if ( ! var_eval_filter_includes ) passed = ! passed;
  return passed;
}


void Mask::external_vcf_iteration( const std::string & f )
{
  ext_vcffile = f;
  ext_vcf = true;
}


void Mask::include_var_obs_file( const std::vector<std::string> & f )
{
  obs_file_filter = true;
  for (int i=0;i<f.size(); i++)
    {
      int id = vardb->file_tag( f[i] );
      if ( id != 0 ) inc_obs_file.insert(id);      
    }  
}

void Mask::require_var_obs_file( const std::vector<std::string> & f )
{
  obs_file_filter = true;
  for (int i=0;i<f.size(); i++)
    {
      int id = vardb->file_tag( f[i] );
      if ( id != 0 ) req_obs_file.insert(id);      
    }  
}

void Mask::exclude_var_obs_file( const std::vector<std::string> & f )
{
  obs_file_filter = true;
  for (int i=0;i<f.size(); i++)
    {
      int id = vardb->file_tag( f[i] );
      if ( id != 0 ) exc_obs_file.insert(id);      
    }  
}

void Mask::include_var_alt_file( const std::vector<std::string> & f )
{
  alt_file_filter = true;
  for (int i=0;i<f.size(); i++)
    {
      int id = vardb->file_tag( f[i] );
      if ( id != 0 ) 
	inc_alt_file.insert(id); 
    }  
}

void Mask::require_var_alt_file( const std::vector<std::string> & f )
{ 
  alt_file_filter = true;
  for (int i=0;i<f.size(); i++)
    {
      int id = vardb->file_tag( f[i] );
      if ( id != 0 ) req_alt_file.insert(id);      
    }  
}

void Mask::exclude_var_alt_file( const std::vector<std::string> & f )
{
  alt_file_filter = true;
  for (int i=0;i<f.size(); i++)
    {
      int id = vardb->file_tag( f[i] );
      if ( id != 0 ) exc_alt_file.insert(id);
    }
}


void Mask::include_var_alt_group( const std::vector<std::string> & f )
{
  alt_group_filter = true;
//   for (int i=0;i<f.size(); i++)
//     {      
//       if ( id != 0 ) inc_alt_group.insert(id);      
//     }  
}

void Mask::require_var_alt_group( const std::vector<std::string> & f )
{ 
  alt_group_filter = true;
//   for (int i=0;i<f.size(); i++)
//     {
//       if ( id != 0 ) req_alt_group.insert(id);      
//     }  
}

void Mask::exclude_var_alt_group( const std::vector<std::string> & f )
{
  alt_group_filter = true;
//   for (int i=0;i<f.size(); i++)
//     {
//       if ( id != 0 ) exc_alt_group.insert(id);
//     }
}

void Mask::count_obs_file(int i, int j)
{
  obs_file_count = i;
  obs_file_max = j;
}

void Mask::count_alt_file(int i, int j )
{
  alt_file_count = i;
  alt_file_max = j;
}


bool Mask::eval_obs_file_filter( const Variant & v ) const
{

  if ( ! obs_file_filter ) return true;  
  
  // Needs at least 1 inclusion filter, if 1+ includes
  // Needs all required filters, if 1+ requires
  // Needs to not have any excludes
  
  // Excludes
  
  std::set<int>::iterator i = exc_obs_file.begin();
  while ( i != exc_obs_file.end() )
    {
      if ( v.file_present( *i ) ) return false;
      ++i;
    }

  // Requires
  
  i = req_obs_file.begin();
  while ( i != req_obs_file.end() )
    {
      if ( ! v.file_present( *i ) ) return false;
      ++i;
    }
  
  if ( req_obs_file.size() > 0 ) return true;
  
  
  // Includes 
  
  if ( inc_obs_file.size() == 0 ) return true;
  
  i = inc_obs_file.begin();
  while ( i != inc_obs_file.end() )
    {
      if ( v.file_present( *i ) ) return true;
      ++i;
    } 
  
  return false;

}


bool Mask::eval_alt_file_filter( Variant & v ) const      
{

  // do we see at least one alternate-allele in present files
  // i.e. post-filtering;  doesn't assume that a variant/file is 
  // present only if at least one non-ref
  
  if ( ! alt_file_filter ) return true;  
  
  // Needs at least 1 inclusion filter, if 1+ includes
  // Needs all required filters, if 1+ requires
  // Needs to not have any excludes
  

  // Excludes

  std::set<int>::iterator i = exc_alt_file.begin();
  while ( i != exc_alt_file.end() )
    {      
      if ( v.has_nonreference_by_file( *i ) ) { return false; } 
      ++i;
    }
  
  // Requires
  
  i = req_alt_file.begin();
  while ( i != req_alt_file.end() )
    {
      
      // sample must be present, implicitly
      if ( ! v.fsample_svar_counts( *i ) ) return false;
      if ( ! v.has_nonreference_by_file( *i ) ) return false;
      ++i;
    }
  
  if ( req_alt_file.size() > 0 ) return true;
  

  // Includes 
  
  if ( inc_alt_file.size() == 0 ) return true;
  
  i = inc_alt_file.begin();
  while ( i != inc_alt_file.end() )
    {      
      if ( v.has_nonreference_by_file( *i ) ) return true;      
      ++i;
    } 
  
  return false;

}

bool Mask::eval_alt_group_filter( const Variant & v ) const
{
  // do we see at least one alternate-allele in present files
  // i.e. post-filtering;  doesn't assume that a variant/file is 
  // present only if at least one non-ref

  return true;
}

bool Mask::eval_file_count( Variant & v ) const
{

  // Do we have to see at least N samples with the variant observed?
  if ( obs_file_count || obs_file_max ) 
    {
      if ( v.n_uniq_samples() < obs_file_count ) return false;
      if ( obs_file_max && v.n_uniq_samples() > obs_file_max ) return false;
    }

  // Do we have to see at least N samples with at least one non-reference variant?
  if ( alt_file_count || alt_file_max ) 
    {
      
      const int n = v.n_uniq_samples();
      int c=0;
      for (int s = 0 ; s < n ; s++ )
	{	  	  
	  if ( v.has_nonreference_by_file( s ) ) ++c;
	  if ( c == alt_file_count ) return true;
	  if ( alt_file_max && c >  alt_file_max   ) return false;	  
	}
      return false;
    }
  return true;
}


bool Mask::fail_on_sample_variant() const 
{
  return fail_on_sample_variant_allow == 0;
}

bool Mask::test_fail_on_sample_variant(int n , int m ) const 
{
  if ( n >= fail_on_sample_variant_allow && fail_on_sample_variant_allow != -1 ) return false;
  if ( m < fail_on_sample_variant_require ) return false;
  return true;
}

void Mask::fail_on_sample_variant( int n , int m )
{
  fail_on_sample_variant_allow = n;
  fail_on_sample_variant_require = m;
}

void Mask::set_load_variant_meta( const std::string & m )
{
  load_vmeta = true;
  load_partial_vmeta = true;
  load_vmeta_list.insert( m );
}

bool Mask::load_variant_meta( const std::string & m ) const
{
  return load_partial_vmeta ? load_vmeta_list.find( m ) != load_vmeta_list.end() : load_vmeta ; 
}

void Mask::set_load_genotype_meta( const std::string & m )
{
  load_gmeta = true;
  load_partial_gmeta = true;
  load_gmeta_list.insert( m );
}

bool Mask::load_genotype_meta( const std::string & m ) const
{
  return load_partial_gmeta ? load_gmeta_list.find( m ) != load_gmeta_list.end() : load_gmeta ; 
}


void Mask::ploidy( const std::string & chr , ploidy_t t )
{
  fixxy_map_chr[ chr ] = t;
}

void Mask::set_pseudo_autosomal( const Region & region ) 
{
  fixxy_map_par.push_back( region );
}

ploidy_t Mask::ploidy( const std::string & chr ) const
{
  // default is autosomal.
  std::map<std::string,ploidy_t>::const_iterator i = fixxy_map_chr.find( chr );
  return i == fixxy_map_chr.end() ? PLOIDY_AUTOSOMAL : i->second;
}

bool Mask::pseudo_autosomal( const Variant & var ) const
{
  // for a more focused look, default is UNKNOWN (i.e. in practice, we will
  // only ask this of the X chromosome, to determine PAR regions) 
  std::vector<Region>::const_iterator i = fixxy_map_par.begin();
  while ( i != fixxy_map_par.end() )
    {      
      if ( i->contains( var ) ) return true;     
      ++i;
    }
  return false;
}


void Mask::prep_segmask()
{

  for (int i = 0 ; i < GP->indmap.size() ; i++ )
    {
      
      Individual * person = GP->indmap(i);
      
      std::set<int>::iterator grp = in_segset.begin();
      while ( grp != in_segset.end() )
	{
	  std::set<Region> s = segdb->get_indiv_regions( *grp , person->id() );
	  std::set<Region>::iterator ss = s.begin();
	  while ( ss != s.end() )
	    {
	      segs[ i ].push_back( *ss );
	      ++ss;
	    }
	  ++grp;
	}

      grp = req_segset.begin();
      while ( grp != req_segset.end() )
	{
	  std::set<Region> s = segdb->get_indiv_regions( *grp , person->id() );
	  std::set<Region>::iterator ss = s.begin();
	  while ( ss != s.end() )
	    {
	      req_segs[ i ][ *grp ].push_back( *ss );
	      ++ss;
	    }
	  ++grp;
	}

      grp = ex_segset.begin();
      while ( grp != ex_segset.end() )
	{
	  std::set<Region> s = segdb->get_indiv_regions( *grp , person->id() );
	  std::set<Region>::iterator ss = s.begin();
	  while ( ss != s.end() )
	    {
	      ex_segs[ i ].push_back( *ss );
	      ++ss;
	    }
	  ++grp;
	}
    }
}


bool Mask::in_any_segmask( const Region & var , const std::vector<Region> & segs )  
{
  for ( int s = 0 ; s < segs.size() ; s++ ) 
    if ( var.overlaps( segs[s] ) ) return true;
  return false;
}


bool Mask::in_all_segmask( const Region & var , const std::map<int,std::vector<Region> > & segs )  
{
  if ( segs.size() == 0 ) return true;  
  std::map<int,std::vector<Region> >::const_iterator g = segs.begin();
  while ( g != segs.end() )
    {      
      if ( ! in_any_segmask( var , g->second ) ) return false;
      ++g;
    }
  return true;
}


bool Mask::eval_segmask( const int i , const Region & region )
{

  
  // excludes 
  if ( ex_segs.size() && ex_segs[i].size() ) 
    {
      if ( in_any_segmask( region , ex_segs[i] ) ) return false;
    }

  // requires
  if ( req_segs.size() ) 
    {
      if ( req_segs[i].size() )
	return in_all_segmask( region , req_segs[i] );      
      else 
	return false;
    }

  // includes 
  if ( in_segset.size() == 0 ) return true;  
  if ( in_any_segmask( region , segs[i] ) ) return true;

  return false;
}


bool Mask::forced( int achr , int abp1, int abp2, 
		   int bchr , int bbp1 , int bbp2 , 
		   Region * next ) const 
{
  
  // 'previous' variant is 'a'
  // next in varbd iteration is 'b'

  // just check whether or not we have another variant in the force list that 
  // should come before this. 

  // if a == b , implies that we are at the end of the list w.r.t. vardb variants. 
  
  // of course, could be the case that we have 1 or more force variants anyway
  
  std::set<Region>::iterator upr = force_vlist.upper_bound( Region(achr,abp1,abp2) );
  if ( upr == force_vlist.end() ) return false;
  if ( ! ( *upr < Region(bchr,bbp1,bbp2) ) ) return false;

  // otherwise, implies that this next guy should be inserted
  *next = *upr;
  return true;
}
