
#include "genotype.h"
#include "variant.h"

#include <bitset>
#include <iostream>
using namespace std;

std::map<std::string,Genotype> Genotype::gcache;

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


std::vector<int> Genotype::allele_list( const int na ) const
{
  std::vector<int> ac(na);
  if ( is_null ) return ac;
  ac[ allele1 ]++;
  if ( ! is_haploid ) ac[ allele2 ]++;
  return ac;
}



// Alternate allele count for a simple genotype

int Genotype::allele_count( const std::string & acode , const Variant * parent ) const
{ 
  int c = -1;

  int na = parent->n_alleles();
  for ( int a = 0 ; a < na; a++ )
    {
      if ( parent->allele(a).name() == acode ) 
	{
	  c = a; 
	  break;
	}
    }

  if ( c == -1 ) return 0;

  return allele_count( c );

}

uint8_t Genotype::bcf() const
{  
  uint8_t gt;
  gt = is_null << 7 | known_phase << 6 | allele1 << 3 | allele2 ;
  return gt;
}


void Genotype::bcf( uint8_t gt )
{
  is_null     = ( gt >> 7 ) & 1 ;
  known_phase = ( gt >> 6 ) & 1 ;  
  allele1     = ( gt >> 3 ) & 7 ;
  allele2     =   gt        & 7 ;
}

uint32_t Genotype::pack() const
{

  // Use two uint8_t to store the two genotype alleles (PB only has uint32 class)
  // a set of four boolean flags will also be added to VARDB, to indicate
  //  null, etc
  
  uint32_t gt = more() << 19 
    | is_null << 18
    | is_haploid << 17
    | known_phase << 16
    | allele1 << 8 
    | allele2 ; 
  
  return gt;
  
}

bool Genotype::unpack( uint32_t gt )
{
  is_null     = ( gt >> 18 ) & 1 ; 
  is_haploid  = ( gt >> 17 ) & 1 ;
  known_phase = ( gt >> 16 ) & 1 ;
  allele1     = ( gt >> 8  ) & 255 ;
  allele2     =   gt         & 255 ;
  
  // returns T is okay
  // returns F if genotype is encoded by _GT in meta-information

  return ( gt >> 19 ) & 1;
}


bool Genotype::equivalent( const Genotype & a , const Genotype & b )
{
  // an exact match?
  if ( a == b ) return true;

  // ignore phase when comparing heterozygotes here
  if ( ! ( a.is_haploid || a.is_null || b.is_haploid || b.is_null ) ) 
    return a.allele1 == b.allele2 && a.allele2 == b.allele1 ;   
}



//
// GenotypeSet functions: reading commands, which can be redirected
//

int GenotypeSet::size() const 
{ 

//   std::cout  << "X = ";
//   if ( incon == NULL ) std::cout << "incon NULL ";
//   if ( svar ) std::cout << "but svar not";
//   std::cout << "\n";

  return svar ? incon->size() : calls.size(); 
}

// Return Genotype for individual i
Genotype & GenotypeSet::genotype(int i) 
{ 
  return calls[i]; 

//    std::cout << " svar = " << ( svar ? "Y" : "NULL" ) << " "
//   	    << i << " of " << calls.size() << " "
//   	    << ( svar ? svar->calls.size() : -1 ) << "\n";
  
  // return svar ? svar->calls.genotype( (*incon)[i] ) : calls[i]; 
}

/// Return const Genotype for individual i

const Genotype & GenotypeSet::genotype(int i) const {  return calls[i]; } 
//   std::cout << " svarX = " << svar << " "
// 	    << i << " if " << calls.size() << " "
//   	    << ( svar ? svar->calls.size() : -1 ) << "\n";
//   return svar ? svar->calls.genotype( (*incon)[i] ) : calls[i]; 
// }
        

// void GenotypeSet::set_consensus_slotmap( SampleVariant * ps , std::vector<int> * pm )
// { 
//   svar = ps;
//   incon = pm; 
// } 


void Genotype::set_from_string( const std::string & gtok )
{
  
  const Genotype * cached = search_genotype_cache( gtok );
  
  if ( cached ) 
    {
      *this = *cached;
      return;
    }

  
  size_t idx_slash = gtok.find( "/" );
  
  if ( idx_slash != std::string::npos ) // unphased genotype
    {
      int a1,a2;
      bool okay = true;
      if ( ! Helper::str2int( gtok.substr(0,idx_slash) , a1 ) ) okay = false; 
      if ( ! Helper::str2int( gtok.substr(idx_slash+1) , a2 ) ) okay = false; 
      if ( okay ) genotype(a1,a2);
      else set_null(); 
    }
  else 
    {
      size_t idx_pipe  = gtok.find( "|" );
      if ( idx_pipe != std::string::npos ) // phased genotype
	{
	  int a1, a2;
	  bool okay = true;
	  if ( ! Helper::str2int( gtok.substr(0,idx_pipe) ) , a1 ) okay = false;
	  if ( ! Helper::str2int( gtok.substr(idx_pipe+1) ) , a2 ) okay = false;
	  if ( okay ) phased_genotype(a1,a2);		    
	  else set_null();
	}
      else // assume haploid 
	{	      
	  int a;	      
	  if ( gtok == "." || ! Helper::str2int( gtok , a ) ) set_null();
	  else genotype(a);
	}
    }

  // and add to cache
  gcache[ gtok ] = *this;
  
}

Genotype::Genotype( const std::string & str )
{
  set_from_string(str);
}

Genotype::Genotype( const std::string & str , 
		    const int gt_field , 
		    const std::vector<meta_index_t*> & formats )
{ 
  
  std::vector<std::string> tok = Helper::char_split( str , ':' );
  
  if ( gt_field >= tok.size() ) 
    {
      set_null();
    }
  else 
    {      
      set_from_string( tok[ gt_field ] );
    }
    
  // Set genotpe meta-fields (GT field should be NULL and thus skipped)
  
  meta.set( tok , &formats );
  
}


void Genotype::clear_genotype_cache()
{
  gcache.clear();
}
  
const Genotype * Genotype::search_genotype_cache( const std::string & tok )
{
  std::map<std::string,Genotype>::iterator i = gcache.find( tok );
  return i != gcache.end() ? &(i->second) : NULL ;
}

  

