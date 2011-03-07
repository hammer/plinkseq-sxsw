
#include "genotype.h"
#include "spec.h"
#include "variant.h"

#include <bitset>
#include <iostream>
using namespace std;



void GenotypeSet::summarise_meta( std::map<meta_typed_key_t,std::pair<int,int> > & mk ,  
				  std::map<meta_typed_key_t,std::string> & mk_name , 
				  std::set<meta_typed_key_t> & mk_varlen )  const
{

  // mk;          mk key-pair (type/key) --> pair< # times seen , constant length >
  // mk_name;     mk key-pair (type/key) --> string name
  // mk_varlen;   is variable-length?

  mk.clear();
  mk_name.clear();
  mk_varlen.clear();
  
  const unsigned int n = calls.size();
  
  for ( unsigned int i = 0 ; i < n; i++)
    {
      
      //
      // Get/test length of item 
      //
      
      const MetaInformation<GenMeta> & meta = calls[i].meta;

      std::vector<meta_typed_key_t> keys = meta.typed_keys();
      
      for (unsigned int k=0; k < keys.size(); k++) 
	{

	  // Have we seen this key before? If so, what length 
	  
	  std::pair<int,int> & val = mk[ keys[k] ];
	  
	  // A var-length field ?
	  
	  if ( val.first != 0 )
	    {	      	      
	      if ( meta.size( keys[k] ) != val.second ) 
		mk_varlen.insert( keys[k] );
	    }
	  else
	    val.second = meta.size( keys[k] ) ;

	  // Note that we've seen this.
	  
	  val.first++;
	  
	}
    }
  
  // Add string labels
  std::map<meta_typed_key_t,std::pair<int,int> >::iterator i = mk.begin();
  while ( i != mk.end() )
    {
      mk_name[ i->first ] = MetaInformation<GenMeta>::field( i->first.first , i->first.second );
      ++i;
    }

}


// Allele count for allele 'altcode' for a particular allele 
// non-simple genotype

int Genotype::allele_count( const int altcode ) const
{
  VariantSpec * s = parent->consensus.specification();
  if ( !s ) return 0;
  return s->allele_count( gcode , altcode ) ;
}

std::vector<int> Genotype::allele_list() const
{
  if ( more() ) plog.warn("TODO: allele_list() not handling multi-allelic markers yet");
  std::vector<int> ac(2);
  ac[0] = pat();
  ac[1] = mat();
  return ac;
}


// ANy alternate allele count for a simple genotype

int Genotype::allele_count( ) const
{ 
  if ( null() ) 
    {
      return -1;           // missing
    }
  else if ( more() )       // multi-allelic
    {
      VariantSpec * s = parent->consensus.specification();
      if ( !s ) return 0;
      return s->copy_count( gcode ) - s->allele_count( gcode , 0 );
    }
  else // basic biallelic variant
    {
      return (int)pat() + (int)mat();
    }
}


// Alternate allele count for a simple genotype

int Genotype::allele_count( const std::string & acode ) const
{ 
  if ( null() ) return 0; // missing
  else if ( more() )       // multi-allelic
    {
      VariantSpec * s = parent->consensus.specification();
      if ( !s ) return 0;
      return s->allele_count( gcode , acode );
    }
  else // basic biallelic variant
     {
      if ( parent->consensus.ref == acode ) return copy_number() - allele_count();
      else return allele_count();  // hmm, check this is okay
    }
}


uint8_t Genotype::bcf() const
{  
  uint8_t gt;
  gt = null() << 7 | phased() << 6 | pat() << 3 | mat() ;
  return gt;
}


void Genotype::bcf( uint8_t gt )
{
  null( ( gt >> 7 ) & 1 );
  phased( ( gt >> 6 ) & 1 );  
  int pg = ( gt >> 3 ) & 7;
  int mg =  gt & 7;
  
  pat( pg );
  mat( mg );
  
  // more than a SNP?
  
  more( pg > 1 || mg > 1 );
  
  if ( more() )
    {
      // how to set gcode from pg and mg?
      std::cerr << "note-- not handling multi-allelic markers in BCF conversion yet...\n";
    }
  
}
