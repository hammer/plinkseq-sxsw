
#include "spec.h"
#include "helper.h"
#include "genotype.h"
#include "helper.h"

using namespace std;
using namespace Helper;

extern Log plog;

std::vector<meta_index_t*> * VariantSpec::formats = NULL;
int VariantSpec::gt_field;

void specDecoder::displayCache()
{
  std::map<std::string,CachedSpec>::iterator i = cache.begin();
  while ( i != cache.end() )
    {
      plog << "[" << i->first << "]\t" << i->second.n << "\n";
      ++i;
    }
}


int VariantSpec::allele_count( const int i, const int a )
{
  if ( i < 0 || i >= genotypes.size() ) return 0;
  return genotypes[i]->alleleCount(a);
}


VariantSpec * specDecoder::decode(std::string s )
{


  std::map<std::string,CachedSpec>::iterator i = cache.find(s);
  
  if ( i != cache.end() ) 
    {	  
      ++(i->second.n);
      return i->second.p;
    }
  
  
  // First check if we need to prune the cache 

  if ( cache.size() > 100 && false )
    {
      map<string,CachedSpec>::iterator i = cache.begin();
      while ( i != cache.end() )
	{
	  if ( i->second.n == 1 ) 
	    {
	      delete i->second.p;
	      cache.erase(i);
	    }
	  ++i;
	}      
    }
  
  // Create a new specification, store in cache and return a
  // pointer

  VariantSpec * p = decode_from_string(s);
  if ( p == NULL ) return NULL;

  CachedSpec c;
  c.p = p;
  c.n = 1;
  cache.insert(make_pair( s , c ) );
  return p;
  
}

VariantSpec * specDecoder::decode_from_string(const std::string & encoding )
{
  std::vector<std::string> s = Helper::char_split( encoding , ' ' );  
  if ( s.size() != 3 ) Helper::halt("internal error in VariantSpec");
  return decode_from_string( s[0] ,s[1] ,s[2] );
}


VariantSpec * specDecoder::decode_from_string(const std::string & format, 
					      const std::string & ref,
					      const std::string & alternate)
{


  //
  // Generate listing of possible genotypes
  //

  //
  // That is, given a list of alleles, specify haploid,
  // dipolid(unphased) and diploid (phased) genotype codes.
  //

  std::vector<GenotypeSpec*> genotypes;
  
  //
  // Map of possible GT strings back to actual genotypes
  //

  std::map<std::string,int> genotypeMap;
  std::map<std::string,int> genotypeMap_ACGT;
  
  // Restore defaults
  
  reset();

  
  // Basic structure of genotype/meta-information: use colon-delimited
  // list

  std::vector<std::string> slots = Helper::char_split( format , ':' );
  

  //
  // Reference allele
  //

  alleles.insert( AlleleCode( ref , 0 ) ) ;  
  
  std::vector<std::string> alts = Helper::char_split( alternate , ',' );
  
  for (int i=0; i<alts.size(); i++)
    alleles.insert( AlleleCode( alts[i] , i+1 ) );
  
  phased = true;

  haploid = false;
    

  /////////////////////////////////////////////////////
  // Construct genotypes given alleles
  
  //
  // Haploid data
  //

  std::set<AlleleCode>::iterator a1 = alleles.begin();	  
  while ( a1 != alleles.end() )
    {
      genotypes.push_back( new HaploidSpec( *a1 ) );
      int g = genotypes.size() - 1;
      genotypeMap.insert(make_pair( genotypes[g]->num_print() , g ) );
      genotypeMap_ACGT.insert(make_pair( genotypes[g]->print() , g ) );
      ++a1;
    }	    
  

  //
  // Diploid data
  //

  a1 = alleles.begin();
  while ( a1 != alleles.end() )
    {
      set<AlleleCode>::iterator a2 = a1;
      while ( a2 != alleles.end() )
	{		  
	  genotypes.push_back( new DiploidSpec( *a1 , *a2 , false , false ) );
	  int g = genotypes.size() - 1;
	  
	  genotypeMap.insert(make_pair( genotypes[g]->num_print() , g ) );
	  genotypeMap_ACGT.insert(make_pair( genotypes[g]->print() , g ) );
	  
	  genotypeMap.insert(make_pair( genotypes[g]->num_print_reverse_order() , g ));
	  genotypeMap_ACGT.insert(make_pair( genotypes[g]->print_reverse_order() , g ));
	  
	  ++a2;
	}
      ++a1;
    }


  //
  // Phased diploid data
  //
  
  a1 = alleles.begin();
  while ( a1 != alleles.end() )
    {      
      std::set<AlleleCode>::iterator a2 = alleles.begin();
      while ( a2 != alleles.end() )
	{		  
	  genotypes.push_back( new DiploidSpec( *a1 , *a2 , true , false ) );
	  int g = genotypes.size() - 1;
	  
	  genotypeMap.insert(make_pair( genotypes[g]->num_print() , g ) );
	  genotypeMap_ACGT.insert(make_pair( genotypes[g]->print() , g ) );	  

	  ++a2;
	}
      ++a1;
    }
  

  //
  // Phased diploid data, but w/ high switch prob
  //
  
  a1 = alleles.begin();
  while ( a1 != alleles.end() )
    {      
      std::set<AlleleCode>::iterator a2 = alleles.begin();
      while ( a2 != alleles.end() )
	{		  
	  genotypes.push_back( new DiploidSpec( *a1 , *a2 , true , true ) );
	  int g = genotypes.size() - 1;
	  genotypeMap.insert(make_pair(	genotypes[g]->num_print() , g  )); 
	  genotypeMap_ACGT.insert(make_pair(	genotypes[g]->print() , g  )); 
	  ++a2;
	}
      ++a1;
    }
  

  
  VariantSpec * p = new VariantSpec(genotypes, 
				    genotypeMap, 
				    genotypeMap_ACGT,
				    slots);
  
  return p;
  
}

void VariantSpec::display()
{
  if ( this == NULL ) 
    {
      plog << " [ NULL ] \n";
      return;
    }

  plog << "---------------------\n";
  
  plog << "Meta-information\n";
  
  for (int i=0; i<slots.size(); i++)
    {
      plog << slots[i] << "\n";
    }
    
  plog << "----------\n";
  
  std::map<std::string,int>::iterator i = genotypeMap.begin();

  while ( i != genotypeMap.end() )
    {
      plog << i->first << "\t"
		<< i->second << "\t"
		<< genotypes[i->second]->print() << "\t"
		<< genotypes[i->second]->phased() << "\t"
		<< genotypes[i->second]->haploid() << "\t"
		<< genotypes[i->second]->pswitch() << "\t"
		<< genotypes[i->second]->copyCount() << "\t"
		<< genotypes[i->second]->more() << "\n";
      
      ++i;
    }      
  
  plog << "---------------------\n";

}


Genotype VariantSpec::callGenotype( const std::string & s, const Variant * parent, bool acgt )
{      
  
  Genotype g;
  
  //  g.variant( parent );
  
  std::vector<std::string> tok = Helper::char_split( s , ':' );
  
  int gi = -1 ;

  // If we have a format set, use that 
  if ( formats ) 
    {

      // Get genotype code
      if ( gt_field != -1 )
	{
	  if ( acgt )
	    {
	      std::map<std::string,int>::iterator k = genotypeMap_ACGT.find( tok[ gt_field ] );
	      gi = k != genotypeMap_ACGT.end() ? k->second : -1 ;
	    }
	  else
	    {
	      std::map<std::string,int>::iterator k = genotypeMap.find( tok[ gt_field ] );
	      gi = k != genotypeMap.end() ? k->second : -1 ;
	    }
	}
      
      // Set genotpe meta-fields
      g.meta.set( tok , formats );
      
    } 
  else 
    {

      //
      // Extract out genotype; create string of meta-information 
      // for parsing
      //
      
      for (int i = 0; i < tok.size(); i++)
	{
	  
	  if ( tok[i] == "" ) continue;
	  
	  if ( i >= slots.size() ) break;	   
	  
	  if ( slots[i] == "GT" )
	    {
	      if ( acgt )
		{
		  std::map<std::string,int>::iterator k = genotypeMap_ACGT.find( tok[i] );
		  gi = k != genotypeMap_ACGT.end() ? k->second : -1 ;
		}
	      else
		{
		  std::map<std::string,int>::iterator k = genotypeMap.find( tok[i] );
		  gi = k != genotypeMap.end() ? k->second : -1 ;
		}
	    }
	  else
	    {
	      if ( tok[i] != PLINKSeq::VCF_MISSING_CHAR() )
		g.meta.parse_set( slots[i] , tok[i] );
	    }
	  
	  // Get next token/slot
	}
    }


  //
  // Process genotype code, gi
  // Code of -1 means failed call
  //
  
  g.missing(false);
  
  if ( gi == -1 ) 
    {
      g.failed(true);
    }
  else
    {
      GenotypeSpec * gs = genotypes[gi];
    
      // Set genotype settins
      g.failed(false);
      g.pat( gs->pat() );
      g.mat( gs->mat() );
      g.haploid( gs->haploid() );
      g.phased( gs->phased() );
      g.pswitch( gs->pswitch() );      
      g.more( gs->more() );
      if ( g.more() ) 
	g.code(gi);
    }
  
  return g;
}
