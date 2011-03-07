#include <memory>

#include "variant.h"
#include "genotype.h"
#include "options.h"
#include "helper.h"
#include "mask.h"
#include "bcf.h"

#include "gstore.h"

extern GStore * GP;

using namespace std;
using namespace Helper;

specDecoder SampleVariant::decoder;

std::ostream & operator<<( std::ostream & out, const SampleVariant & v)
{ 
  //  out << GP->vardb.file_tag( v.fset )
  out    << v.ref << "/" << v.alt;
  //      << ":" << v.filter_info;      
  return out;
}

std::string SampleVariant::file_name() const 
{
  return GP->vardb.attached() ? GP->vardb.file_tag( fset ) : ".";
}

bool Variant::make_consensus( IndividualMap * a )
{

  //
  // We are entering this under one of three possible states:
  //
  // 0) not a 'multi-sample' alignment, which is most simple case.
  //    Here the consensus SampleVariant will have exactly what we
  //    need, so we can leave immediately.
  //
  // 1) a 'flat' alignment, in which genotypes have already been
  //    deposited in the consensus slot. This implies >1 file, but no
  //    overlap of samples. MetaInformation from each file will still
  //    be available in each SampleVariant however. Up to this point,
  //    we will not have checked to see whether the multiple
  //    SampleVariants have consistent allele coding; if not, it will
  //    need to be reconciled.
  //
  // 2) a non-'flat' alignment, meaning that at least one genotype
  //    overlaps; here, we need to reconcile everything into a
  //    consensus.  The consensus will not already be populated with
  //    genotypes.
  //


  //
  // Record the alignment used
  //
  
  align = a;

  
  //
  // Also allow for possibility that we have >1 obs of a variant in a
  // single file, and so we will want to force a "non-flat"
  // alignment. If we make it to the end of this function, then
  // infile_overlap() could have returned T, so unset it. Otherwise,
  // no need to worry
  //


  if ( infile_overlap() ) 
    align->force_unflat( true );


  //
  // If this is a flat alignment (no dupe genotypes) then nothing else
  // needs to be done, as we would have previously stored all genotype
  // information in the consensus SampleVariant already.
  //  
  
  if ( ( ! align->multi_sample() ) && svar.size() < 2 ) 
    {
      
      int n_alleles = consensus.parse_alleles();
      
      // for biallelic markers, we can leave now
      if ( n_alleles == 2 ) return true;
      
      consensus.set_allelic_encoding();
      
      return true;
    }
  


  //
  // Align basic allelic informaiton
  //
  
  SampleVariant & first = svar[0];
  
  consensus.ref = first.ref;
  consensus.alt = "";
  
  bool need_to_resolve = false;
  std::set<std::string> obs_alleles;

  bool expanded_ref = false;

  for (int i = 0 ; i < svar.size(); i++ )
    {

      SampleVariant & sv = svar[i];
      
      //
      // Check this is the same reference allele
      //
      
      if ( sv.offset != consensus.offset || sv.ref != consensus.ref ) 
	{

	  // does a simple alteration merge things? e.g  GCCC versus GCC
	  // if not, flag here as problematic...

	  if ( ! SampleVariant::align_reference_alleles( sv , consensus ) )
	    {
	      
	      // Hmm... not sure this is right way to deal with things?
	      // should probably just halt?

	      Helper::halt( "incompatible REF sequences " + coordinate() );	      

	      //need_to_resolve = true;
	      //obs_alleles.insert( sv.ref );
	    }
	  else
	    expanded_ref = true;
	}
      
    }

  
  if ( expanded_ref ) 
    {
      // v. inefficient to redo all this, but hopefully these resolving functions
      // will only be needed very occassionally

      // here we add a flag to also ensure that the possible *alternate* alleles
      
      for (int i = 0 ; i < svar.size(); i++ )
	SampleVariant::align_reference_alleles( svar[i] , consensus , true );
    }
  

  // now do alternate alleles
  
  for (int i = 0 ; i < svar.size(); i++ )
    {
      
      SampleVariant & sv = svar[i];
      
      sv.set_allelic_encoding();

      // 
      // Is the alternate alleles specification exactly the same?
      //
      
      if ( sv.alt != consensus.alt ) 
	{
	  if ( consensus.alt != "" ) 
	    need_to_resolve = true;
	  
	  int n_alleles = sv.parse_alleles();
	  
	  // skip sv.ref here (handled above)
	  for (int a = 1 ; a < n_alleles; a++ )
	    {
	      if ( obs_alleles.find( sv.alleles[a].name() ) == obs_alleles.end() )
		{
		  obs_alleles.insert( sv.alleles[a].name() );
		  if ( consensus.alt != "" ) consensus.alt += ",";
		  consensus.alt += sv.alleles[a].name();
		}
	    }
	}
      
    }
  

  //
  // Parse this allele string for the Variant (i.e. create Allele
  // objects in alleles[] )  
  //
  
  int n_con_alleles = consensus.parse_alleles();   // ultimately, can get rid of this


  //
  // Did we encounter SampleVariant with different allelic-specifications? If so, 
  // we will need to re-call the consensus genotypes
  // 

  if ( need_to_resolve )
    {
      for (int a=1;a<n_con_alleles; a++)
	obs_alleles.insert( consensus.alleles[a].name() );
    }
  

  // i.e. some redundancy now we are using the full VariantSpec, below:
    
  consensus.set_allelic_encoding();


  //
  // Under a flat alignment, we can leave now, as genotype data will
  // already be in consensus
  //
  
  // but what if we need a re-coding (i.e. need_to_resolve=T) ??

  if ( align->flat() ) 
    {
      return true;
    }


  //
  // Genotypes (which might need to be recoded)
  //
  
  const int n = align->size();
  const int ns = n_samples();

  consensus.calls.size( n );


  for (int i=0; i<n; i++)
    {
      
      int2 j = align->unique_mapping(i);      
      
      // Does this individual feature in multiple samples? 
      
      bool multiple_samples = j.p1 != -1 ? ftosv[ j.p1 ].size() > 1 : true;

      if ( ! multiple_samples )
	{ 
	  
	  // not sure this would ever fail, but just incase (perhaps
	  // also need to check ftosv/j.p1 above for range?)
	  
	  if ( ftosv[ j.p1 ].size() == 1 ) 
	    {
	      
	      SampleVariant * p = psample( ftosv[ j.p1 ][0] );
	      
	      if ( p )  
		{		  
		  
		  Genotype g = (*p)(j.p2);
		  
		  // do we need to assign a different underlying encoding
		  // scheme for this variant? if so, recall as is from
		  // textual input
		  
		  if ( need_to_resolve )
		    {
		      std::string label = p->label( g );		  
		      g = consensus.spec->callGenotype( label , this , true );  //T=ACGT mode
		    }	      	      
		  
		  consensus(i) = g;
		  consensus(i).meta = (*p)(j.p2).meta;
		}
	    }
	} 	
      else
	{
	  
	  // This individual potentially features in >1
	  // samples/file. (Although, for this particular variant, we
	  // may only see a single SampleVariant, however.) In
	  // addition, if need_to_resolve=T, then we cannot just copy
	  // the Genotype object: we need to re-call using the
	  // consensus spec.
	  
	  // If more than one is placed in the consensus slot, wipe
	  // any sample-specific meta-information

	  bool slot_filled = false;

	  
	  // Get the list of SampleVariant IDs that this individual features in for this Variant
	  
	  std::vector<int> svids;
	  std::vector<int> svslot;
	  
	  if ( j.p1 != -1 )  // single file, but multi obs within that file
	    {
	      std::vector<int> & xx = ftosv[ j.p1 ];
	      for (int z=0; z<xx.size(); z++)
		{
		  svids.push_back( xx[z] );
		  svslot.push_back( j.p2 );
		}
	    }
	  else
	    {
	      std::set<int2> k = align->multiple_mapping(i);    
	      std::set<int2>::iterator ki = k.begin();
	      while ( ki != k.end() )
		{
		  std::vector<int> & xx = ftosv[ ki->p1 ];
		  for (int z=0; z<xx.size(); z++)
		    {
		      svids.push_back( xx[z] );
		      svslot.push_back( ki->p2 );
		    }
		  ++ki;
		}
	    }
	  
	  
	  //
	  // Now look at the genotypes in each of these samples:
	  //

	  for (int k = 0 ; k < svids.size(); k++ )
	    {
	      
	      SampleVariant * p = psample( svids[k] );
	      
	      if ( p ) 
		{
		  
		  Genotype g = (*p)( svslot[k] );
		  
		  // Do not automatically advance genotype meta-information in  this scenario however
		  // (it will stay with the individual SampleVariants)
		 
		  if ( need_to_resolve )
		    {
		      std::string label = p->label( g );
		      g = consensus.spec->callGenotype( label , this , true ); // T=ACGT mode
		    }
		  
		  // Ignore null genotypes
		  
		  if ( g.null() ) continue; 
		  
		  // If consensus is null, insert here

		  if ( consensus(i).null() ) 
		    {
		      consensus(i) = g;
		      if ( slot_filled ) 
			consensus(i).meta.clear();
		      slot_filled = true;
		      continue;
		    }
		  

		  // otherwise, if a discordant call, set consensus to missing
		  
		  if ( g != consensus(i) )
		    {
		      consensus(i).null( true );
		      consensus(i).meta.clear();	
		      break;
		    }		  
		  
		  
		  // if more than one sample-variant fits the consensus slot, do not try 
		  // to copy meta-information over to the consensus slot

		  if ( slot_filled ) 
		    consensus(i).meta.clear();

		  slot_filled = true;

		}
	      	      
	    }    
	  
	} // end of 'need to resolve multiple instances' part

    } // next individual
      

  // In case we forced a non-flat encoding, revert back to normal now

  align->force_unflat( false );

  return true;
}




std::string Variant::print_PED(const Genotype & g, const std::string & delim ) const
{

  stringstream s; 
  
  if ( g.more() ) // do not display multi-allelic variants for now
    s << "0" << delim << "0";
  else
    {
      if( g.null() )
	s << "0" << delim << "0";      
      else
	{
	  
	  if ( g.pat() ) 
	    s << consensus.alt; else s << consensus.ref;

	  s << delim;
	  
	  if ( g.haploid() )
	    {
	      if ( g.pat() ) 
		s << consensus.alt; else s << consensus.ref;		  
	    }
	  else	    
	    {
	      if ( g.mat() ) 
		s << consensus.alt; else s << consensus.ref;
	    }
	}
    }
  
  return s.str();
}





std::string Variant::label( const int i  , const std::string & delim ) const
{   
  
  // get basic textual representation of genotype from SampleVariant function
  std::string s = consensus.label( consensus(i) );
  
  // under a flat alignment, we are done now
  if ( flat() && ! infile_overlap() ) return s;

  // otherwise, list these other genotype calls from individual SVs
  
  std::map<int, const Genotype *> gm = all_genotype(i);
  std::map<int, const Genotype *>::iterator j = gm.begin();
  if ( gm.size() > 1 )
    {
      s += " {";
      j = gm.begin();	  
      while ( j != gm.end() )
	{
	  SampleVariant * svar = psample( j->first );
	  if ( svar ) 
	    s +=  ( j != gm.begin() ? delim : "" ) + svar->label( *(j->second) );	    
	  ++j;
	}
      s += "}";
    }
  return s;
}


std::string Variant::gmeta_label( const int i , const std::string & delim ) const
{

  std::stringstream ss;
  ss <<  consensus(i).meta; 

  if ( flat() ) return ss.str();

  std::map<int, const Genotype *> gm = all_genotype(i);	  
  std::map<int, const Genotype *>::iterator j = gm.begin();
  if ( gm.size() > 1 )
    {
      ss << "{";
      j = gm.begin();	  
      while ( j != gm.end() )
	{
	  SampleVariant * svar = psample( j->first );
	  if ( svar ) 
	    {
	      
	      ss << ( j != gm.begin() ? delim : "" ) << j->second->meta ;
	      
	    }
	  ++j;
	}
      ss <<  "}";
    }
  return ss.str();
}


std::string Variant::sample_label( const int i, const std::string & delim ) const
{
  std::string s;
  std::map<int, const Genotype *> gm = all_genotype( i );	  
  std::map<int, const Genotype *>::iterator j = gm.begin();	
  while ( j != gm.end() )
    {
      s += ( j == gm.begin() ? "" : delim ) + GP->vardb.file_tag( svtof[ j->first ] ) ;
      ++j;
    }  
  return s == "" ? "." : s;
}


int SampleVariant::parse_alleles()
{
  
  //
  // We can improve performance here by use of a cache
  //
  
  alleles.clear();
  
  //
  // Reference allele
  //
  
  alleles.push_back( Allele( ref ) );
  
  //
  // One or more alternate alleles
  //

  std::vector<std::string> tok2 = Helper::char_split(alt , ',');
  
  for (int i=0; i<tok2.size(); i++)
    {
      alleles.push_back( Allele( tok2[i] ) );
    }

  // 
  // If we only have two alleles, this is fine. If more than two, we need to build the 
  //
  
  return alleles.size();
}



int Variant::n_alleles() const
{
  return consensus.alleles.size();    
}


int Variant::n_nonreference() const
{
  int n = 0;
  for (int i=0; i< consensus.calls.size(); i++) 
    if ( genotype(i)->nonreference() ) ++n;
  return n;
}

int Variant::n_null() const
{
  int n = 0;
  for (int i=0; i< consensus.calls.size(); i++) 
    if ( genotype(i)->null() ) ++n;
  return n;
}

int Variant::n_notnull() const
{
  int n = 0;
  for (int i=0; i< consensus.calls.size(); i++) 
    if ( genotype(i)->notnull() ) ++n;
  return n;
}

bool Variant::n_minor_allele( int & m , int & n , const affType & aff ) const
{
  
  // returns counts of any non-reference allele
  // optionally, if aff is not UNKNOWN_PHE , split this by case/control
  
  m = 0;
  n = 0;

  if ( aff == UNKNOWN_PHE )
    {
      for (int i=0; i< size() ; i++) 
	{
	  const Genotype * g = genotype(i);
	  if ( g->notnull() )
	    {
	      int a = g->allele_count();
	      if ( a >= 0 ) 
		{
		  m += a;
		  n += g->copy_number();
		}
	    }      
	}
    }
  else  // alternate version, in which we look up phenotype
    {
      for (int i=0; i< size() ; i++) 
	{
	  if ( ind(i)->affected() == aff ) 
	    {
	      const Genotype * g = genotype(i);
	      if ( g->notnull() )
		{
		  int a = g->allele_count();
		  if ( a >= 0 ) 
		    {
		      m += a;
		      n += g->copy_number();
		    }
		}      
	    }
	}      
    }


  // If reference allele is minor allele, swap the minor allele count,
  // and let the caller know about this

  if ( (double)m / (double)n > 0.5 ) 
    {
      m = n - m;
      return false;
    }
  
  // Alternate allele is minor allele
  return true;  
}


std::map<std::string,int> SampleVariant::allele_counts( const affType & aff , const Variant * parent ) const 
{
  
   std::map<std::string,int> c;
  
//   for (int i=0; i< calls.size() ; i++) 
//     {
//       const int j = GP->indmap.get_slot( fset , i );
//       Individual * person = GP->indmap(j);

//       bool count = true;
//       if ( aff == CASE && person->affected() != CASE ) count = false;
//       if ( aff == CONTROL && person->affected() != CONTROL ) count = false;
      
//       const Genotype * g = &calls.genotype(i);
      
//       if ( g->null() )
// 	{
// 	  c["."]++;
// 	}
//       else
// 	{
// 	  std::map<std::string,int> a = allele_count(i);
// 	  std::map<std::string,int>::iterator k = a.begin();
// 	  while ( k != a.end() )
// 	    {
// 	      c[ k->first ] += count ? k->second : 0;
// 	      ++k;
// 	    }
// 	}
//     }
  return c;
}


std::map<std::string,int> Variant::allele_counts( const affType & aff ) const 
{
  return consensus.allele_counts( aff, this );
}


std::map<std::string,int> SampleVariant::genotype_counts( const affType & aff , const Variant * parent , bool unphased ) const 
{

  std::map<std::string,int> c;

  const int n = GP->indmap.size( fset );  
  for (int i=0; i < n ; i++)
    {
      bool count = true;

      if ( aff != UNKNOWN_PHE ) 
	{
	  Individual * person = parent->ind( GP->indmap.get_slot( fset , i ) );      	 
	  if ( aff == CASE && person->affected() != CASE ) count = false;
	  if ( aff == CONTROL && person->affected() != CONTROL ) count = false;
	}

      // so we also track unobserved counts (e.g. get 0 for genotype not seen in controls)
      // using the function below will correctly point to consensus if needed

      if ( parent->flat() ) 
	c[ parent->consensus.label( *(parent->genotype( this , i ) ) , unphased ) ] += count ? 1 : 0 ;
      else
	c[ label( *(parent->genotype( this , i ) ) , unphased ) ] += count ? 1 : 0 ;

    }

  return c;
}


bool Variant::biallelic() const
{
  // this function of sensitive to down-coding status: 
  return consensus.alleles.size() == 2;
}

bool Variant::multiallelic() const
{
  return consensus.alleles.size() > 2;
}

bool Variant::monomorphic() const
{
  if ( consensus.alleles.size() == 1 ) return true;
  
  // TODO -- see if this is now needed below?
  //         (and if so, needed above, so perhaps precompute?)
  //          when making the variant in first place? 

  int n=0,m=0;
  n_minor_allele(n,m);
  return n==0 || n==m;
}

bool Variant::simple_snp() const
{
  // not sensitive to down-coding
  return consensus.alleles.size() == 2 && 
    consensus.alleles[1].type() == ALLELE_SUBSTITUTION;
}

bool Variant::transition() const
{
  if ( ! simple_snp() ) return false;
  return (    ( consensus.alt == "A" && consensus.ref == "G" )  
	   || ( consensus.alt == "G" && consensus.ref == "A" ) 
	   || ( consensus.alt == "C" && consensus.ref == "T" )  
	   || ( consensus.alt == "T" && consensus.ref == "C" ) );

}

bool Variant::transversion() const
{
  if ( ! simple_snp() ) return false;
  return ! transition();
}

bool Variant::simple_ins() const
{
  return consensus.alleles.size() == 2 && 
    consensus.alleles[1].type() == ALLELE_INSERTION;
}

bool Variant::simple_del() const
{
  return consensus.alleles.size() == 2 && 
    consensus.alleles[1].type() == ALLELE_DELETION;
}


string Variant::VCF()
{
  // Construct a string that is a VCF format entry
  
  std::ostringstream s;
  
  // VCF is tab-delimited 

  s << Helper::chrCode( chr ) << "\t"
    << bp << "\t"
    << vname << "\t"
    << consensus.ref << "\t"
    << consensus.alt << "\t"
    << consensus.qual << "\t"
    << consensus.filter_info << "\t"
    << consensus.meta << "\t";
  
  // Format field for genotype info:
  // Just take from the first genotype for now (although, 
  // we should improve this subsequently?)

  
  //
  // Format for genotypes (scan all genotypes to get fullest possible)
  //
  
  s << "GT";
  
  std::set<std::string> allkeys;
  for (int i = 0 ; i < size(); i++)
    {
      std::vector<std::string> keys = consensus.calls.genotype(i).meta.keys();
      for (unsigned int j=0; j<keys.size(); j++) allkeys.insert( keys[j] );
    }
  
  std::set<std::string>::iterator i = allkeys.begin();
  while ( i != allkeys.end() )
    {
      s << ":" << *i;
      ++i;
    }
  
  bool withmeta = allkeys.size() > 0;
  
  //
  // Output genotype calls for all individuals
  //
  
  for (int i = 0 ; i < size(); i++)
    {
      const Genotype * g = genotype(i);
      s << "\t" << consensus.num_label( *g );
      if ( withmeta ) s << ":" << g->meta.printValues( allkeys , ":");
    }  
  
  s << "\n";

  return s.str();
  
}



blob SampleVariant::encode_BLOB() const
{

  // Using Protocol Buffers, generate a BLOB object that represents
  // this Variant Specification in variant.proto -> variant.pb.h and
  // variant.pb.cc

  VariantBuffer pbVar;
  
  // Set core variant fields (index, name, chr and position (bp1,bp2)
  // handled separately)
  
  pbVar.set_alt( alt );
  pbVar.set_ref( ref );
  pbVar.set_quality( qual );
  pbVar.set_strand( vstrand );


  // TODO: Just add the filter (string) as a single element (for now)
  //  actually, we should change filter_info to be a parsed list

  pbVar.add_filter( filter_info );
  
  
  //
  // Set variant meta-data
  //
  
  std::vector<std::string> keys = meta.keys();
  
  for (unsigned int k=0; k<keys.size(); k++) 
    {
      
      VarMetaBuffer * m = pbVar.add_vmeta();
      
      meta_index_t midx = MetaInformation<VarMeta>::field( keys[k] );
      
      m->set_name( keys[k] );
      
      // Number of elements in meta-value?
      
      int num = meta.size( keys[k] );
      
      switch ( midx.mt ) {
      case META_INT :
	{
	  m->set_type( VarMetaBuffer::INT );
	  const std::vector<int> * v = meta.ptr_int( midx.key );
	  for (int j=0; j<num; j++)
	    m->add_int_value( (*v)[j] );
	  break;
	}
      case META_FLOAT :
	{
	  m->set_type( VarMetaBuffer::FLOAT );
	  const std::vector<double> * v = meta.ptr_double( midx.key );
	  for (int j=0; j<num; j++)
	    m->add_double_value( (*v)[j] );
	  break;
	}
      case META_TEXT :
	{
	  m->set_type( VarMetaBuffer::TEXT );
	  const std::vector<std::string> * v = meta.ptr_string( midx.key );
	  for (int j=0; j<num; j++)
	    m->add_string_value( (*v)[j] );
	  break;
	}
      case META_BOOL :
	{
	  m->set_type( VarMetaBuffer::BOOL );
	  const std::vector<bool> * v = meta.ptr_bool( midx.key );
	  for (int j=0; j<num; j++)
	    m->add_bool_value( (*v)[j]  );
	  break;
	}
      case META_FLAG :
	{
	  m->set_type( VarMetaBuffer::BOOL );
	  // TODO -- shoudn't something else be set here??
	  m->add_bool_value( true );
	  break;
	}
      default :
	{
	  m->set_type( VarMetaBuffer::TEXT );
	  const std::vector<std::string> * v = meta.ptr_string( midx.key );
	  for (int j=0; j<num; j++)
	    m->add_string_value( (*v)[j] );
	  break;
	}
      }
      
    }
  
 
  //
  // Set all genotypes 
  //
  
  
  unsigned int n = calls.size();
  
  for ( unsigned int i = 0 ; i < n; i++)
    {

      //
      // Add genotype to BLOB for this individual
      // 
 
     const Genotype & g = calls.genotype(i);
      
      pbVar.add_geno1( g.pack() );	
      
      if ( g.more() )
	pbVar.add_geno2(g.code());	

    }


  //
  // Genotype meta-data
  //


  // Also, compile list of all keys for this variant (i.e. allow new ones to
  // have been added since original VCF specification)
  //
  // Also, determine the proportion of missing items, and also the whether or 
  // not the length of the information is variable
  //
  
  std::map<meta_typed_key_t,std::pair<int,int> >  mk;     // key --> # inds seen in / length
  std::map<meta_typed_key_t,std::string>  mk_name;        // key  --> name
  std::set<meta_typed_key_t> mk_varlen;                   // is variable-length?

  // populate the above

  calls.summarise_meta( mk , mk_name , mk_varlen ); 
  
  
  //
  // Write genotype meta-information
  //
  
  std::map<meta_typed_key_t,std::pair<int,int> >::iterator k = mk.begin();

  while ( k != mk.end() )
    {

      GenotypeMetaBuffer * m = pbVar.add_gmeta();

      const std::string mk_label = mk_name[ k->first ];
      const mType       mk_type  = k->first.first;
      const meta_key_t  mk_key   = k->first.second;
      
      m->set_name( mk_label );	
      
      //      meta_index_t midx = MetaInformation<GenMeta>::field( k->first );

           
      // Do all non-missing entries have similar length? If so, do not
      // bother storing a length field
      
      bool constant_length = mk_varlen.find( k->first ) == mk_varlen.end(); 
      
      if ( constant_length ) 
	{
	  m->set_fixed_len( k->second.second );
	}
      
      // Do all individuals have a non-missing entry?
      
      // If inds: 1..10
      //   1 2 3 4 5 6 7 8 9 10   -->  all_nonmissing, no indexes
      //   . . . . . 6 . 8 . .    -->  majority missing, use indiv_index
      //   1 2 3 4 . 6 7 8 9 10   -->  majority not missing, use missing_index
      
      bool all_nonmissing = k->second.first == (int)n;
      
      bool use_missing_index = false;
      
      if ( all_nonmissing ) 
	{
	  m->set_fixed_indiv( true );
	}
      else
	{
	  use_missing_index = (double)( k->second.first ) / ( double)(n) > 0.5;
	}


      //
      // Save value for each individual
      //
      
      for (unsigned int i = 0; i < n ; i++)
	{
	  
	  if ( calls.genotype(i).meta.has_field( k->first ) )
	    {
	      
	      if ( ! ( all_nonmissing || use_missing_index ) )
		{
		  m->add_indiv_index( i );
		}

	      const MetaInformation<GenMeta> & meta = calls.genotype(i).meta;
	      
	      switch ( mk_type )
		{
		case META_INT :
		  {
		    m->set_type( GenotypeMetaBuffer::INT );
		    if ( ! constant_length )
		      m->add_len( meta.size( k->first ) );
		    const std::vector<int> * d = meta.ptr_int( mk_key );
		    for (unsigned int j = 0 ; j < d->size(); j++) 
		      m->add_int_value( (*d)[j] );
		    break;
		  }
		case META_FLOAT :
		  {
		    m->set_type( GenotypeMetaBuffer::FLOAT );
		    if ( ! constant_length )
		      m->add_len( meta.size( k->first) );
		    const std::vector<double> * d = meta.ptr_double( mk_key );
		    for (unsigned int j = 0 ; j < d->size(); j++) 
		      m->add_double_value( (*d)[j] );
		    break;
		  }
		case META_TEXT :
		  {
		    m->set_type( GenotypeMetaBuffer::TEXT );
		    if ( ! constant_length )
		      m->add_len( meta.size( k->first ) );
		    const std::vector<string> * d  = meta.ptr_string( mk_key );
		    for (unsigned int j = 0 ; j < d->size(); j++) 
		      m->add_string_value( (*d)[j] );
		    break;
		  }
		case META_BOOL :
		  {
		    m->set_type( GenotypeMetaBuffer::BOOL );
		    if ( ! constant_length )
		      m->add_len( meta.size( k->first ) );
		    const std::vector<bool> * d = meta.ptr_bool( mk_key );
		    for (unsigned int j = 0 ; j < d->size(); j++) 
		      m->add_bool_value( (*d)[j] );
		    break;
		  }
		default : // add as text if unsure
		  {
		    m->set_type( GenotypeMetaBuffer::TEXT );
		    if ( ! constant_length )			    
		      m->add_len( meta.size( k->first ) );
		    const std::vector<std::string> * d = meta.ptr_string( mk_key );
		    for (unsigned int j = 0 ; j < d->size(); j++) 
		      m->add_string_value( (*d)[j] );
		  }
		  
		}
	    }
	  else if ( use_missing_index )
	    {
	      m->add_missing_index( i );
	    }
	  
	}
	
      // Next genotype meta-field
      ++k;
    }

    

    
  // Return a BLOB 
  
  string s;
  
  pbVar.SerializeToString(&s);
  
  // Debug code:

  //  plog >> pbVar.DebugString() << "\n\n";
  
  return blob(s);
  
}



void SampleVariant::store_BLOB(blob & b)
{
  buf.ParseFromString( b.getString() );
}

// bool SampleVariant::decode_BLOB(blob & b , Variant * parent )
// {
//   store_BLOB(b);
//   return decode_BLOB( parent );
// }

bool SampleVariant::decode_BLOB( Variant * parent , 
				 IndividualMap * align , 
				 Mask * mask )
{
  
  SampleVariant * target = ( ! align->multi_sample() ) ? &(parent->consensus) : this ;
  
  SampleVariant * genotype_target = align->flat() ? &(parent->consensus) : this ;  
  
  // Some of the core fields that aren't already populated
  // extracted either to consensus or self
  
  decode_BLOB_basic( target );
  
  // Extract variant meta-information, either to consensus or self,
  // optionally copying to Variant population-level static meta-fields
  
  if ( ! decode_BLOB_vmeta( mask , parent , target ) ) return false;
    
  // Extract genotypes, either into consensus or self
  // this also returns a variant-level bool, if a include expression w/ a gfunc is in the mask

  if ( ! decode_BLOB_genotype( align , mask , parent , this , target , genotype_target ) ) return false;
  
  return true;
}


bool SampleVariant::decode_BLOB_basic( SampleVariant * svar )
{      

  // If this SampleVariant was created from a BCF, we already will have 
  // populated these fields, nothing to do here

  if ( svar->bcf || svar->vcf_direct ) return true; 

  //
  // The index, chr, name and bp1,2 will be already set
  // for the main Variant. Here we assign SampleVariant
  // specific information.  If we have a non-null pointer
  // to the parent Variant, this imples a single-sample
  // Variant, and so use the Variant's consensus SampleVariant
  // instead of that pointed to here.
  //
  
  //
  // Assign allele codes and other core features
  //
  
  svar->alt = buf.alt();
  svar->ref = buf.ref();
  svar->qual = buf.quality();
  svar->vstrand = buf.strand();

  string my_filter = "";
  unsigned int num = buf.filter().size();
  for (int i=0;i<num; i++)    
    my_filter += buf.filter(i);  
  

  //
  // To help with some old VCF files, translate 0 and . to mean PASS
  //

  if ( my_filter == "0" || my_filter == "." ) 
    my_filter = PLINKSeq::PASS_FILTER();


  //
  // Set and parse into meta_filter
  //

  svar->filter( my_filter );

  return true;
}



bool SampleVariant::decode_BLOB_vmeta( Mask * mask, Variant * parent , SampleVariant * sample )
{    

  // SampleVariant meta-information
  
  // Two functions here:
  //  1) to extract from a PB 
  //  2) to apply some variants (i.e. if things fail here saves time not to have to extract all genotype information)
  

  // For BCF-derived SVs, or those read direct from a VCF,skip the first step
  // i.e. genotypes already extracted

  if ( ! ( sample->bcf || sample->vcf_direct ) ) 
    {

      unsigned int num = buf.vmeta().size();
  
      for (unsigned int k=0; k<num; k++)
	{	
	  
	  // Any fields that are flagged in MetaMeta::pop_static() are 
	  // also copied to the consensus variant:
	  
	  bool incon = parent && MetaMeta::static_variant( buf.vmeta(k).name() ) ; 
	  
	  switch ( buf.vmeta(k).type() ) {
	    
	  case VarMetaBuffer::INT :
	    {
	      vector<int> d( buf.vmeta(k).int_value().size() );
	      for ( unsigned int j = 0 ; j < d.size(); j++ ) 
		d[j] = buf.vmeta(k).int_value(j);
	      sample->meta.set( buf.vmeta(k).name() , d );
	      if ( incon ) parent->meta.set( buf.vmeta(k).name() , d );
	      break;
	    }
	  case VarMetaBuffer::FLOAT :
	    {
	      vector<double> d( buf.vmeta(k).double_value().size() );
	      for ( unsigned int j = 0 ; j < d.size(); j++ ) 
		d[j] = buf.vmeta(k).double_value(j);
	      sample->meta.set( buf.vmeta(k).name() , d );
	      if ( incon ) parent->meta.set( buf.vmeta(k).name() , d );
	      break;
	    }
	  case VarMetaBuffer::BOOL :
	    {
	      // Stored in PB as BOOL, but this could either be a flag or a bool
	      
	      meta_index_t midx = MetaInformation<VarMeta>::field( buf.vmeta(k).name() );
	      if ( midx.mt == META_BOOL ) 
		{
		  vector<bool> d( buf.vmeta(k).bool_value().size() );
		  for ( unsigned int j = 0 ; j < d.size(); j++ ) 
		    d[j] = buf.vmeta(k).bool_value(j);
		  sample->meta.set( buf.vmeta(k).name() , d );
		  if ( incon ) parent->meta.set( buf.vmeta(k).name() , d );
		}
	      else if ( midx.mt == META_FLAG ) 
		{
		  sample->meta.set( buf.vmeta(k).name() );
		  if ( incon ) parent->meta.set( buf.vmeta(k).name() );	      
		}
	      break;
	    }
	  default :
	    {
	      
	      //
	      // NOTE: legacy -- we now save FLAGs as BOOL in PB, so will only need the above
	      //
	      
	      meta_index_t midx = MetaInformation<VarMeta>::field( buf.vmeta(k).name() );
	      if ( midx.mt == META_FLAG ) 
		{	      
		  sample->meta.set( buf.vmeta(k).name() );
		  if ( incon ) parent->meta.set( buf.vmeta(k).name() );	      
		}
	      else
		{
		  vector<string> d( buf.vmeta(k).string_value().size() );
		  for ( unsigned int j = 0 ; j < d.size(); j++ ) 
		    d[j] = buf.vmeta(k).string_value(j);
		  sample->meta.set( buf.vmeta(k).name() , d );
		  if ( incon ) 
		    parent->meta.set( buf.vmeta(k).name() , d );
		}
	      break;
	    }
	    
	  }
	}

    } // end of extracting/assigned from PB/BLOB
  

  //
  // Any variant meta-information level filters?
  //

  if ( mask && mask->variant_mask() ) 
    {
      if ( ! mask->eval( *sample ) ) return false;
    }

  //
  // Any include="expr" filters (that do not require genotype data)? 
  // Ones that require genotype data are deferred until decode_BLOB_genotype() 
  // 

  if ( mask && mask->filter_expression() && ! mask->filter_expression_requires_genotypes() )
    {
      if ( ! mask->calc_filter_expression( *sample ) ) return false; 
    }

  return true;
  
}


bool SampleVariant::decode_BLOB_genotype( IndividualMap * align , 
					  Mask * mask , 
					  Variant * parent , 
					  SampleVariant * source ,
					  SampleVariant * vtarget , 
					  SampleVariant * target ) // genotype_target
{
  
  
  //
  // Sample A : 1 2 3
  // Sample B : 4 5 6
  // Sample C : 3 4 

  // Consensus
  // Uniq      1(A1)  2(A2)                       5(B2)  6(B3)
  // Mult                    3(A3,C1)  4(B1,C2)
  //

  
  //
  // Decode genotype information, for 0+ individuals If an optional
  // alignment is specified, use that to only extract that subset
  //
  
  
  if ( ! ( target->bcf || target->vcf_direct ) ) 
    {
      
      // Number of individuals in PBuffer
      unsigned int n_buffer = buf.geno1().size();
      
      // Number of individuals we actually want
      unsigned int n_variant = align ? align->size() : n_buffer ;
      

      //
      // Allocate space as needed
      //

      target->calls.size( n_variant );

  
      //
      // Add basic genotype information
      //
  
      int j = 0; // geno2 counter
      
      for ( unsigned i=0; i<n_buffer; i++)
	{
	  
	  int slot = i;
	  if ( align )
	    {
	      
	      // get within-sample slot
	      slot = align->sample_remapping( source->fileset() , i ) ;
	      
	      // under a flat alignment, we need to reset to the
	      // consensus slot: check -- or should this be (fset,i) ?
	      
	      if ( align->flat() ) 
		slot = align->get_slot( source->fileset() , slot );
	      
	    }
	  
	  
	  //
	  // Basic genotype
	  //
	  
	  Genotype g( parent );
	  
	  g.unpack( buf.geno1(i) );
      
	  if ( g.more() ) 
	    g.code( buf.geno2( j++ ) );
	  
	  if ( slot != -1 )        	      
	    target->calls.add( g, slot );
	}
      

      //
      // Append genotype meta-information
      //
      
      unsigned int m = buf.gmeta().size();


      for ( unsigned int k = 0 ; k < m ; k++ )
	{
	  
	  // Does this have a set length, how are missing individuals
	  // handlded?
      
	  // Mode: 
	  //  0) constant_length or no? 
	  
	  //  1) all_nonmissing -- no index, just list all values
	  //  2) indiv_index    -- specify IDs for people w/ data
	  //  3) missing_index  -- specify IDs for people w/out data
	  
	  bool constant_length = buf.gmeta(k).has_fixed_len();
	  int  length = 0;
	  if ( constant_length ) 
	    length = buf.gmeta(k).fixed_len();
	  
	  // Should only get one of these:
	  bool missing_index = buf.gmeta(k).missing_index().size() > 0;
	  bool indiv_index = buf.gmeta(k).indiv_index().size() > 0;
	  bool all_nonmissing = ! ( missing_index || indiv_index );
	  
	  // Set buffer sizes

	  int nlen = indiv_index ?
	    (int)buf.gmeta(k).indiv_index().size() :
	    n_buffer;

	  int idx = 0;

	  if ( all_nonmissing ) 
	    {

	      switch (  buf.gmeta(k).type() ) {
		
	      case GenotypeMetaBuffer::INT :
		{
	      
		  for (int j=0;j<nlen; j++)
		    {
		      
		      idx = target->addIntGenMeta( j , source->fileset() , 
						   buf, align, k, idx, 
						   constant_length ? length : buf.gmeta(k).len(j) );		  
		    }
		  break;
		}
	      case GenotypeMetaBuffer::FLOAT :
		{
		  for (int j=0;j<nlen; j++)
		    {
		      idx = target->addFloatGenMeta( j,  source->fileset() ,
						     buf, align, k, idx, 
						     constant_length ? length : buf.gmeta(k).len(j) );
		    }
		  break;
		}
	      case GenotypeMetaBuffer::BOOL :
		{
		  for (int j=0;j<nlen; j++)
		    {
		      idx = target->addBoolGenMeta( j, source->fileset() ,
						    buf, align, k, idx, 
						    constant_length ? length : buf.gmeta(k).len(j) );
		    }		   
		  break;
		}
	      default :
		{
		  for ( int j=0;j<nlen; j++)
		    {
		      idx = target->addStringGenMeta( j, source->fileset() ,
						      buf, align, k, idx, 
						      constant_length ? length : buf.gmeta(k).len(j) );
		    }
		  
		}
	      }
	    }
	  
	  //
	  // Or using an individual-index?
	  //
	  
	  else if ( indiv_index )
	    {
	      switch (  buf.gmeta(k).type() ) {
		
	      case GenotypeMetaBuffer::INT :
		{
		  for ( int j=0;j<nlen; j++)
		    {
		      idx = target->addIntGenMeta( buf.gmeta(k).indiv_index(j) , source->fileset() , 
						   buf, align, k, idx, 
						   constant_length ? length : buf.gmeta(k).len(j) );
		    }
		  break;
		}
	      case GenotypeMetaBuffer::FLOAT :
		{
		  for ( int j=0;j<nlen; j++)
		    {
		      idx = target->addFloatGenMeta( buf.gmeta(k).indiv_index(j) , source->fileset() , 
						     buf, align, k, idx, 
						     constant_length ? length : buf.gmeta(k).len(j) );
		    }
		  break;
		}
	      case GenotypeMetaBuffer::BOOL :
		{
		  for ( int j=0;j<nlen; j++)
		    {
		      idx = target->addBoolGenMeta( buf.gmeta(k).indiv_index(j) ,  source->fileset() ,
						    buf, align, k, idx, 
						    constant_length ? length : buf.gmeta(k).len(j) );
		    }
		  break;
		}
	      default :
		{		
		  for ( int j=0;j<nlen; j++)
		    {
		      idx = target->addStringGenMeta( buf.gmeta(k).indiv_index(j),  source->fileset() ,
						      buf, align, k, idx, 
						      constant_length ? length : buf.gmeta(k).len(j) );
		    }
		}
		
	      }
	      
	    }
	  
	  
	  else // assume missing index
	    {
	      
	      // Missing index 
	      int skip = buf.gmeta(k).missing_index(0);
	      int scnt = 0;
	      int cnt = 0;
	      
	      switch (  buf.gmeta(k).type() ) {
		
	      case GenotypeMetaBuffer::INT :
		{		    
		  for (int j=0;j<nlen; j++)
		    {
		      if ( j == skip )
			{
			  ++scnt;
			  skip = scnt < buf.gmeta(k).missing_index().size() ? 
			    buf.gmeta(k).missing_index(scnt) : -1 ;
			}
		      else
			{
			  idx = target->addIntGenMeta( j, source->fileset() , 
						       buf, align, k, idx, 
						       constant_length ? length : buf.gmeta(k).len(cnt++) );
			}
		    }
		  break;
		}
	      case GenotypeMetaBuffer::FLOAT :
		{
		  for (int j=0;j<nlen; j++)
		    {
		      if ( j == skip )
			{
			  ++scnt;
			  skip = scnt < buf.gmeta(k).missing_index().size() ? 
			    buf.gmeta(k).missing_index(scnt) : -1 ;
			}
		      else
			{
			  idx = target->addFloatGenMeta( j,  source->fileset() ,
							 buf, align, k, idx, 
							 constant_length ? length : buf.gmeta(k).len(cnt++) );
			}
		    }
		  break;
		}
	      case GenotypeMetaBuffer::BOOL :
		{
		  
		  for (int j=0;j<nlen; j++)
		    {
		      if ( j == skip )
			{
			  ++scnt;
			  skip = scnt < buf.gmeta(k).missing_index().size() ? 
			    buf.gmeta(k).missing_index(scnt) : -1 ;
			}
		      else
			{
			  idx = target->addBoolGenMeta( j,  source->fileset() ,
							buf, align, k, idx, 
							constant_length ? length : buf.gmeta(k).len(cnt++) );
			}
		    }
		  break;
		}
	      default :
		{
		  for ( int j=0;j<nlen; j++)
		    {
		      if ( j == skip )
			{
			  ++scnt;
			  skip = scnt < buf.gmeta(k).missing_index().size() ? 
			    buf.gmeta(k).missing_index(scnt) : -1 ;
			}
		      else
			{
			  idx = target->addStringGenMeta( j ,  source->fileset() ,
							  buf, align, k, idx, 
							  constant_length ? length : buf.gmeta(k).len(cnt++) );
			  
			}
		    }
		  
		}
		
	      }
	      
	    }
	  
	  
	} // Next individual/genotype
  

    } // end of extracting all genotypic (meta) information from PB/BLOB
  else if ( target->bcf ) 
    {     
      
      // but, perhaps we need to now extract from the BCF buffer
      // The format of the genotype is a string in bcf_format
      
//       std::string           bcf_format;
//       int                   bcf_n;
//       std::vector<uint8_t>  bcf_genotype_buf;
      

      // Number of individuals in BCF buffer
      unsigned int n_buffer = target->bcf->sample_n();
      
      // Number of individuals we actually want
      unsigned int n_variant = align ? align->size() : n_buffer ;
      

      //
      // Allocate space as needed
      //

      target->calls.size( n_variant );
      
      
      // Add basic genotype information
      
      std::vector<std::string> format = Helper::char_split( target->bcf_format , ':' );

      // really slow ways to get allele count duplicating previous effort, but 
      // keep for now

      int nalt = Helper::char_split( vtarget->alternate() , ',' ).size() + 1;
      int ngen = (int) (nalt * (nalt+1) * 0.5);

      // create mapping of source-to-target slots
      std::vector<int> s2t( n_buffer );
      
      for ( int i=0; i < n_buffer ; i++ )
	{
	  int slot = i;	  
	  if ( align )
	    {	
	      // get within-sample slot
	      slot = align->sample_remapping( source->fileset() , i ) ;	      
	      // under a flat alignment, we need to reset to the
	      // consensus slot: check -- or should this be (fset,i) ?	      
	      if ( align->flat() ) 
		slot = align->get_slot( source->fileset() , slot );		      
	    }	  
	  s2t[ i ] = slot;	  
	}
      
      // posiiton in genotype buffer
      int p = 0;
      
      // consider each format-slot for all individuals      
      for (int t = 0 ; t < format.size(); t++)
	{	  
	  
	  if ( format[t] == "GT" ) // uint8_t
	    {
	      for ( int i=0; i < n_buffer ; i++ )
		{
		  if ( s2t[i] != -1 )
		    {
		      // std::cout << "mapping genotype " << i << " " << s2t[i] << " " << "\n";
		      target->calls.genotype( s2t[i] ).bcf( target->bcf_genotype_buf[ p ] );
		      // std::cout << "geno = " << target->calls.genotype( s2t[i] ) << "\n";
		    }
		  ++p;
		}
	    }
	  else if ( format[t] == "PL" ) // G * uint8_t
	    {
	      for ( int i=0; i < n_buffer ; i++ )
		{		  
		  if ( s2t[i] != -1 ) 
		    {
		      std::vector<int> pl(ngen);		      
		      for (int j=0;j<ngen;j++)
			pl[j] = target->bcf_genotype_buf[p+j];			  
		      target->calls.genotype( s2t[i] ).meta.set( "PL" , pl );
		    }
		  p += ngen; 
		}
	    }	  
	}
    

      ////
      // Need to figure out what is done, what isn't here
      ///


//       //
//       // Append genotype meta-information
//       //
      
//       unsigned int m = buf.gmeta().size();


//       for ( unsigned int k = 0 ; k < m ; k++ )
// 	{
	  
// 	  // Does this have a set length, how are missing individuals
// 	  // handlded?
      
// 	  // Mode: 
// 	  //  0) constant_length or no? 
	  
// 	  //  1) all_nonmissing -- no index, just list all values
// 	  //  2) indiv_index    -- specify IDs for people w/ data
// 	  //  3) missing_index  -- specify IDs for people w/out data
	  
// 	  bool constant_length = buf.gmeta(k).has_fixed_len();
// 	  int  length = 0;
// 	  if ( constant_length ) 
// 	    length = buf.gmeta(k).fixed_len();
	  
// 	  // Should only get one of these:
// 	  bool missing_index = buf.gmeta(k).missing_index().size() > 0;
// 	  bool indiv_index = buf.gmeta(k).indiv_index().size() > 0;
// 	  bool all_nonmissing = ! ( missing_index || indiv_index );
	  
// 	  // Set buffer sizes

// 	  int nlen = indiv_index ?
// 	    (int)buf.gmeta(k).indiv_index().size() :
// 	    n_buffer;

// 	  int idx = 0;

// 	  if ( all_nonmissing ) 
// 	    {

// 	      switch (  buf.gmeta(k).type() ) {
		
// 	      case GenotypeMetaBuffer::INT :
// 		{
	      
// 		  for (int j=0;j<nlen; j++)
// 		    {
		      
// 		      idx = target->addIntGenMeta( j , source->fileset() , 
// 						   buf, align, k, idx, 
// 						   constant_length ? length : buf.gmeta(k).len(j) );		  
// 		    }
// 		  break;
// 		}
// 	      case GenotypeMetaBuffer::FLOAT :
// 		{
// 		  for (int j=0;j<nlen; j++)
// 		    {
// 		      idx = target->addFloatGenMeta( j,  source->fileset() ,
// 						     buf, align, k, idx, 
// 						     constant_length ? length : buf.gmeta(k).len(j) );
// 		    }
// 		  break;
// 		}
// 	      case GenotypeMetaBuffer::BOOL :
// 		{
// 		  for (int j=0;j<nlen; j++)
// 		    {
// 		      idx = target->addBoolGenMeta( j, source->fileset() ,
// 						    buf, align, k, idx, 
// 						    constant_length ? length : buf.gmeta(k).len(j) );
// 		    }		   
// 		  break;
// 		}
// 	      default :
// 		{
// 		  for ( int j=0;j<nlen; j++)
// 		    {
// 		      idx = target->addStringGenMeta( j, source->fileset() ,
// 						      buf, align, k, idx, 
// 						      constant_length ? length : buf.gmeta(k).len(j) );
// 		    }
		  
// 		}
// 	      }
// 	    }


      
      // ... not yet implemented...
      
    }

  
  
  //
  // Apply any genotype masks? 
  //
  
  const int n_var = target->calls.size();

  if ( mask && mask->assuming_null_is_reference() )
    {
      std::vector<int> t(3);
      t[0] = 0; t[1] = 255; t[2] = 255;
      for ( unsigned int i = 0 ; i < n_var ; i++ )
	{
	  if ( target->calls.genotype(i).null() ) 
	    {
	      target->calls.genotype(i).set_alternate_allele_count(0);
	      
	      // also set PL to indicate this, with slight chance of HET
	      target->calls.genotype(i).meta.set( PLINKSeq::META_GENO_PHRED() , t );
	    }
	}
    }

  if ( mask && mask->genotype_mask() ) 
    {
      for ( unsigned int i = 0 ; i < n_var ; i++ )
	{
	  if ( ! mask->eval( target->calls.genotype(i) ) )
	    target->calls.genotype(i).failed( true );
	}
    }

  
  //
  // Call variant meta-filter (include="expr") that depends on g-functions
  // e.g.  include="g(DP>10) > 0.8 && DB" 
  // note: expressions that did not include any g() would have been evaluated 
  // earlier, so no need to do so again. The mask object let's us know whether or 
  // not the expression requires genotypes
  //

  // Commented out, as this is done in vardb.cpp now
  
  //   if ( mask && mask->filter_expression() && mask->filter_expression_requires_genotypes() ) 
  //     {
  //       if ( ! mask->calc_filter_expression( parent->consensus , *target ) ) return false;
  //     }


  return true;
  
}



bool Variant::frequency_filter( Mask * mask )
{
  
  //
  // Now we have a variant, and fully constructed genotypes
  // (potentially allowing for genotype-level filters, do we have any
  // variant-frequency-based filters?
  //
  
  if ( mask )
    {

      if ( mask->count_filter() )
	{
	  int m = 0; // minor allele
	  int c = 0; // total counts
	  bool altmin = n_minor_allele( m , c );	  	  
	  if ( ! mask->count_filter( m ) ) return false;	       
	}      

      if ( mask->frequency_filter() )
	{
	  int m = 0; // minor allele
	  int c = 0; // total counts	  
	  bool altmin = n_minor_allele( m , c );
	  if ( ! mask->frequency_filter( (double)m/(double)c ) ) return false;
	}
      
    }
  
  return true;

}


bool Variant::null_filter( Mask * mask )
{
  if ( !mask ) return true;
  if ( ! mask->null_filter() ) return true;
  int n = n_null();
  return mask->null_filter( n );
}


bool Variant::case_control_filter( Mask * mask )
{
  if ( ! mask ) return true;
  if ( ! mask->case_control_filter() ) return true;

  int alta = 0 , tota = 0;
  int altu = 0 , totu = 0;
  int dira = n_minor_allele( alta , tota , CASE );
  int diru = n_minor_allele( altu , totu , CONTROL );

  return mask->case_control_filter( alta , altu );
}



//
// Primary access functions
//

// Number of individuals

int Variant::size() const
{
  return align ?  align->size() : 0 ;
}


Individual * Variant::ind(const int i) const
{
  return align->ind(i);
}

Individual * Variant::ind(const std::string & s) const
{
  return align->ind(s);
}

int Variant::ind_n( const std::string & id ) const
{
  return align->ind_n(id);
}

Genotype * Variant::genotype(const int i)
{
  if ( i < 0 || i >= size() ) return NULL;
  return &consensus(i);  
}

const Genotype * Variant::genotype(const int i) const
{
  if ( i < 0 || i >= size() ) return NULL;
  return &(consensus(i));  
}

const Genotype * Variant::genotype( const SampleVariant * svar , const int j) const
{
  // for sample f, get genotype the j'th individual
  // the function below will handle corrupt cases (-1)  
  // Under a flat alignment, no data would have been 
  // stored in the sample-variants 
  // otherwise, get from original slot
  
  if ( align->flat() ) return genotype( align->get_slot( svar->fileset() , j ) ); 
  else return &(svar->calls.genotype(j));            
 
}



std::map<int, const Genotype *> Variant::all_genotype(const int i) const 
{

  // return a map of SV-ID to genotypes for this individual
  
  std::map<int, const Genotype*> g;
  if ( align->unique( i ) )
    {
      
      // this returns a { file-ID, slot-ID } pair
      int2 k = align->unique_mapping(i);
      
      // convert file-ID to a list of SV-IDs (i.e. as an individual
      // could feature >1 time in the same file for the same variant)
      
      std::map<int,std::vector<int> >::const_iterator j = ftosv.find( k.p1 );
      
      if ( j != ftosv.end() )
	{
	  const std::vector<int> & ii = j->second;
	  for (int jj = 0 ; jj < ii.size(); jj++ )
	    g.insert( make_pair( ii[jj], &(svar[ ii[jj] ]( k.p2 ) ) ) );
	}
      
    }
  else
    {
      std::set<int2> k = align->multiple_mapping(i);
      std::set<int2>::iterator ki = k.begin();
      while ( ki != k.end() )
	{	  
	  std::map<int,std::vector<int> >::const_iterator j= ftosv.find( ki->p1 ); 
	  if ( j !=ftosv.end() )
	    {
	      const std::vector<int> & ii= j->second;
	      for (int jj =0 ; jj < ii.size(); jj++ )
		g[ ii[jj] ] = &(svar[ ii[jj] ]( ki->p2 ) );
	    }
	  ++ki;
	}      
    }
  return g;
}



std::map<int, Genotype *> Variant::all_genotype(int i)
{

  std::map<int, Genotype*> g;
  if ( align->unique(i) )
    {
      int2 k = align->unique_mapping(i);

      std::map<int,std::vector<int> >::iterator j = ftosv.find( k.p1 );
      
      if ( j != ftosv.end() )
	{
	  std::vector<int> & ii = j->second;
	  for (int jj = 0 ; jj < ii.size(); jj++ )
	    g.insert( make_pair( ii[jj], &(svar[ ii[jj] ]( k.p2 ) ) ) );
	}
    }
  else
    {
      std::set<int2> k = align->multiple_mapping(i);
      std::set<int2>::iterator ki = k.begin();
      while ( ki != k.end() )
	{

	  std::map<int,std::vector<int> >::iterator j= ftosv.find( ki->p1 );	  
	  if ( j !=ftosv.end() )
	    {
	      std::vector<int> & ii= j->second;
	      for (int jj =0 ; jj < ii.size(); jj++ )
		g[ ii[jj] ] = &(svar[ ii[jj] ]( ki->p2 ) );
	    }	  
	  ++ki;
	}      
    }
  return g;
}


//
// Sample-level access functions
//

std::vector<int> Variant::samples() const 
{
  return svtof;
}


// for a given file #, return the SVAR slot
// or -1 if not present, or that file appears more than once

int Variant::unique_svar_slot( int f ) const
{
  std::map<int,std::vector<int> >::const_iterator i = ftosv.find( f );
  if ( i == ftosv.end() ) return -1;
  if ( i->second.size() != 1 ) return -1;
  return *(i->second.begin());
}

bool Variant::file_present( const int f ) const
{
  std::map<int,std::vector<int> >::const_iterator i = ftosv.find( f );
  if ( i == ftosv.end() ) return false;
  return i->second.size() != 0;
}


std::set<int> Variant::unique_files() const 
{
  std::set<int> f;
  std::map<int,std::vector<int> >::const_iterator i = ftosv.begin();
  while ( i != ftosv.end() )
    {
      f.insert( i->first );
      ++i;
    }
  return f;
}

void Variant::set_first_sample() 
{
  si = 0;
}

SampleVariant & Variant::sample()
{
  return svar.size() == 0 ? consensus : svar[si];
}

bool Variant::next_sample()
{
  return ++si < svar.size();
}

int Variant::sample_n()
{
  return si;
}

std::string Variant::print_meta(const std::string & key , const std::string & delim) const
{
  // static meta-attributes --> in Variant
  // single sample          --> in Consensus SampleVariant
  // otherwise,             --> need to same individual SampleVariants

  if ( MetaMeta::static_variant( key ) ) 
    return meta.as_string( key , "," );
  
  if ( align->single_sample() )
    return consensus.meta.as_string( key , "," );

  std::string r = "";
  for (int i=0; i<svar.size(); i++)
    {
      if ( i != 0 ) r += delim;	    
      r += svar[i].meta.as_string(key, "," ); 
    }
  return r;

}

std::string Variant::print_meta_filter(const std::string & delim) const
{  
  if ( align->single_sample() ) return consensus.filter();  
  std::string r = "";
  for (int i=0; i < svar.size(); i++ )
    {
      std::string s = svar[i].filter();
      if ( s == "" ) s = ".";
      if ( i != 0 ) r += delim;
      r += s;	
    }
  return r;
}

std::string Variant::print_samples( const std::string & delim ) const
{
  std::stringstream s;
  for (int i=0; i< svar.size(); i++) 
    {
      if ( i != 0 ) s << delim;
      s << GP->vardb.file_tag( svar[i].fileset() );      
    }
  return s.str();
}


void SampleVariant::info( const std::string & s , VarDBase * vardb , int file_id ) 
{

  // store original text 
  other_info = s; 

  // parse semi-colon delimited list
  std::vector<std::string> f = Helper::parse(s,";");

  for (int i=0; i<f.size(); i++)
    {

      std::vector<std::string> k = Helper::char_split( f[i] , '=' );
      mType mt = MetaInformation<VarMeta>::type( k[0] );
      
      if ( mt == META_UNDEFINED ) 
	{	  
	  MetaInformation<VarMeta>::field( k[0] , k.size() > 1 ? META_TEXT : META_FLAG , 1 , "undeclared tag" );	  
 	  if ( vardb ) 
 	    vardb->insert_metatype( file_id , k[0] , k.size() > 1 ? META_TEXT : META_FLAG , 1 , META_GROUP_VAR , "undeclared tag" );
	  plog.warn("undefined INFO field (absent in VCF header)", k[0] );
	}
    }

  // parse on ; separators into MetaInformation<VarMeta>   
  // 3rd arg indicates that unknown fields should be 
  // accepted (as a flag, or string value)

  //  -- although, given the above fix, this should no longer
  //     be necessary

  meta.parse(s,";",true); 

}



void SampleVariant::filter( const std::string & s , VarDBase * vardb , int file_id )
{ 
  filter_info = ""; 
  
  // parse semi-colon delimited list
  std::vector<std::string> f = Helper::parse(s,";");
  
  for (int i=0; i<f.size(); i++)
    {
      if ( f[i] == "." || f[i] == "0" ) f[i] = PLINKSeq::PASS_FILTER();      
      
      mType mt = MetaInformation<VarFilterMeta>::type( f[i] );      
      if ( mt == META_UNDEFINED ) 
	{	  
	  MetaInformation<VarFilterMeta>::field( f[i] , META_FLAG , 1 , "undeclared filter tag" );	  
 	  if ( vardb ) 
 	    vardb->insert_metatype( file_id , f[i] , META_FLAG , 1 , META_GROUP_FILTER , "undeclared filter tag" );
	  plog.warn("undefined FILTER field (absent in VCF header)", f[i] );
	}

      meta_filter.set( f[i] );      
      filter_info += f[i];
      if ( i < f.size()-1 ) filter_info += ";";
    }
}



void SampleVariant::set_allelic_encoding()
{
  VariantSpec * ps = SampleVariant::decoder.decode( "GT " + ref + " " + alt );

  VariantSpec::set_format( 0 , NULL ); // 0 just means GT is first field (e.g. if GT:DP:GL etc)
  
  specification( ps ); 

  // TODO: revisit this, but for now a "safe" way to determine this is whether 
  //       alt has a comma in it (i.e. >1 alternate allele specified)

  simple = alt.find(",") == std::string::npos;    

}


std::string SampleVariant::label( const Genotype & g , bool unphased ) const
{
  
  std::stringstream s; 

  const std::string allele_delim = 
    unphased ? "/" : ( g.phased() ? "|" : ( g.pswitch() ? "\\" : "/" ) ) ; 

  if( g.null() )
    {
      s << ( g.haploid() ? "." : "." + allele_delim + "." );
    }  
  else if ( simple ) 
    {
      if ( g.haploid() )
	{
	  s << ( g.pat() ? alt : ref );
	}
      else
	{
	  std::string a1 = g.pat() ? alt : ref;
	  std::string a2 = g.mat() ? alt : ref;
	  if ( unphased && a1 > a2 ) { std::string t=a1;a1=a2; a2=t; }		  
	  s << a1 << allele_delim << a2;		  
	}    
    }
  else
    {      
      if( g.null() )
	s << ( g.haploid() ? "." : "." + allele_delim + "." );      
      else if ( g.more() )	
	{
	  s << spec->printGenotype( g.code() , unphased ) ;
	}
      else
	{	  
	  if ( g.haploid() )
	    {
	      s << ( g.pat() ? alt : ref );
	    }
	  else
	    {
	      std::string a1 = g.pat() ? alleles[1].name() : ref;
	      std::string a2 = g.mat() ? alleles[1].name() : ref;
	      if ( unphased && a1 > a2 ) { std::string t=a1;a1=a2; a2=t; }
	      s << a1 << allele_delim << a2;
	    }
	}
    }

  return s.str();
}


std::string SampleVariant::num_label( const Genotype & g ) const
{
  
  std::stringstream s; 

  if( g.null() )
    s << ( g.haploid() ? "." : "./." );  
  else if ( simple ) 
    {      
      s << ( g.pat() ? "1" : "0" );
      if ( !g.haploid() )
	{	      
	  if ( g.pswitch() ) s << "\\";
	  else if ( g.phased() ) s << "|";
	  else s << "/";	
	  s << ( g.mat() ? "1" : "0" );
	}
    }      
  else
    {      

      if( g.null() )
	s << ( g.haploid() ? "." : "./." );
      else if ( g.more() )	
	s << spec->num_printGenotype( g.code() ) ;
      else
	{
	  if ( g.pat() ) s << "1"; else s << "0";
	  if ( !g.haploid() )
	    {
	      if ( g.pswitch() ) s << "\\";
	      else if ( g.phased() ) s << "|";
	      else s << "/";	
	      if ( g.mat() ) s << "1"; else s << "0";
	    }
	}
    }           
    
  return s.str();
}


bool Variant::concordant( int s1, const Genotype * g1, int s2, const Genotype * g2 ) const
{
  SampleVariant * p1 = psample(s1);
  if ( ! p1 ) return true;
  SampleVariant * p2 = psample(s2);
  if ( ! p2 ) return true;
  return concordant( p1, g1, p2, g2 );
}

bool Variant::concordant( SampleVariant * s1, const Genotype * g1, SampleVariant * s2, const Genotype * g2 ) const
{
  if ( g1->null() || g2->null() ) return true;
  if ( *g1 == *g2 ) return true;
  std::string a1 = s1->label( *g1 , true );
  std::string a2 = s2->label( *g2 , true );
  return a1 == a2;
}


std::map<std::string,int> SampleVariant::allele_count(const int i ) const
{

  const Genotype & g = calls.genotype(i);

  std::map<std::string,int> a;
  
  if( g.null() ) return a;
  
  if ( simple ) 
    {      
      if ( g.haploid() )
	{
	  a[ g.pat() ? alt : ref ]++;
	}
      else
	{
	  a[ g.pat() ? alt : ref ]++;
	  a[ g.mat() ? alt : ref ]++;
	}
    }
  else
    {      
      if ( g.more() )		
	{
	  std::map<std::string,int> b = spec->allele_counts( g.code() );
	  std::map<std::string,int>::iterator i = b.begin();
	  while ( i != b.end() )
	    {
	      a[ i->first ] += i->second;
	      ++i;
	    }
	}
      else
	{
	  if ( g.haploid() )
	    {
	      a[ g.pat() ? alleles[1].name() : ref ]++;
	    }
	  else
	    {
	      a[ g.pat() ? alleles[1].name() : ref ]++;
	      a[ g.mat() ? alleles[1].name() : ref ]++;
	    }
	}
    }
  return a;
}


bool SampleVariant::align_reference_alleles( SampleVariant & s1 , SampleVariant & s2 , bool alt )
{

  // get min and max posisitions on scale where 0 is Variant bp1

  int s1_1  = s1.offset;
  int s1_2  = s1.offset + s1.ref.size() - 1;
  int s2_1  = s2.offset;
  int s2_2  = s2.offset + s2.ref.size() - 1;
  int min_len = s1_1 < s2_1 ? s1_1 : s2_1 ;
  int max_len = s1_2 > s2_2 ? s1_2 : s2_2 ;
  int len = max_len - min_len + 1;

  std::vector<char> conref;
  int c1 = 0, c2 = 0;
  for (int i=min_len; i<=max_len; i++)
    {
      bool in1 = i >= s1_1 && i <= s1_2;
      bool in2 = i >= s2_1 && i <= s2_2;
      if ( in1 && in2 ) 
	{
	  if ( s1.ref[c1] != s2.ref[c2] ) return false;
	  else conref.push_back( s1.ref[c1] );
	  ++c1; ++c2;
	}
      else if ( in1 ) { conref.push_back( s1.ref[c1] ); ++c1; }
      else if ( in2 ) { conref.push_back( s2.ref[c2] ); ++c2; }
      else Helper::halt("internal error in align_ref_alleles()");
    }

  // if we've made it here, then we must have constructed a unifying
  // reference sequence that is okay for both SVs.  In this instance,
  // swap in the new allele coding, along with updating the alternate
  // alleles for each SV

  // padding: left s1_1 - min_len  on LEFT 
  //               s1_2 - max_len  on RIGHT 
  
  // reference 1

  int left1 = s1_1 - min_len;
  int right1 = max_len - s1_2;
  std::string pleft1 = "";
  std::string pright1 = "";
  for (int i=0; i<left1; i++) pleft1 += conref[i];
  for (int i=0; i<right1; i++) pright1 += conref[s1_2+i+1];
  
  int left2 = s2_1 - min_len;
  int right2 = max_len - s2_2;
  std::string pleft2 = "";
  std::string pright2 = "";
  for (int i=0; i<left2; i++) pleft2 += conref[i];
  for (int i=0; i<right2; i++) pright2 += conref[s2_2+i+1];
  
  // update ref alleles
  
  s1.ref = pleft1 + s1.ref + pright1;
  s2.ref = pleft2 + s2.ref + pright2;

  // now for each alternate allele, remake 

  std::vector<std::string> alts = Helper::char_split( s1.alt , ',' );
  s1.alt = "";
  for (int a=0; a<alts.size(); a++)
    {
      if ( a != 0 ) s1.alt += ",";
      s1.alt += pleft1 + alts[a] + pright1;
    }

  alts = Helper::char_split( s2.alt , ',' );
  s2.alt = "";
  for (int a=0; a<alts.size(); a++)
    {
      if ( a != 0 ) s2.alt += ",";
      s2.alt += pleft2 + alts[a] + pright2;
    }
  
  if ( s1.ref != s2.ref ) Helper::halt( "internal error (2) in align_ref_alleles()" );

  // the alt alleles may still be internally inconsistent, but this will be caught in the next section
  
  // this is wasteful / plenty of room for optimisation, but (hopefully) this scenario should not happen v often.

  s1.parse_alleles();
  s2.parse_alleles();

  // align offsets also, now that alleles are changed
  s1.offset = s2.offset = s1.offset < s2.offset ? s1.offset : s2.offset; 
  
  return true;
  
}


void SampleVariant::collapse_alternates( int altcode )
{
  
  // if altcode == 0, then collapse all alt alleles into a single class
  // if altcode > 1, then collapse all non-altcode alt alleles to either missing or reference

  // (For the consensus) recode as a biallelic variant
  
  if ( alleles.size() < 3 ) return;
  if ( altcode > alleles.size() -1 ) return;
  if ( altcode == 0 ) 
    {
      // all ALT become ONE
      alt = alleles[1].name();
      for (int a = 2 ; a < alleles.size(); a++)
	alt += "_" + alleles[a].name();
    }
  else // make other alleles as part of REF
    {
      alt = alleles[altcode].name();
      for (int a = 1 ; a < alleles.size(); a++)
	if ( a != altcode ) ref += "_" + alleles[a].name();
    }
  
  parse_alleles();
  
  for (int i=0; i< calls.size(); i++)
    {
      Genotype & g = calls.genotype(i);
      if ( ! g.null() )
	{
	  if ( altcode )
	    {
	      int ac = g.allele_count( altcode );
	      g.set_alternate_allele_count( ac );
	    }
	  else // merge all alts into ONE
	    {
	      int ac = g.allele_count( ) ;
	      g.set_alternate_allele_count( ac );
	    }
	}
      g.more( false );      
    }

  // as this is a simple biallelic variant, set the following
  simple = true;
  spec = NULL;
  
}


bool SampleVariant::vcf_expand_buffer( Variant * parent )
{

  // If not a valid variant, do not try to expand genotypes
  if ( ! parent->valid() ) return false;
  
  // Call genotypes, add to variant   
  for ( int i=9; i < vcf_direct_buffer.size(); i++)
    {
      Genotype g = spec->callGenotype( vcf_direct_buffer[i] , parent ); 
      calls.add( g );  
    }
  vcf_direct_buffer.clear();
  return true;
}


bool SampleVariant::has_nonreference( const bool also_poly ) const
{
  bool nonref = false;
  std::set<int> npoly;
  for (int i=0; i<calls.size(); i++)
    {
      if ( calls.genotype(i).nonreference() ) 
	{
	  // leave when hit first non-reference call
	  if ( ! also_poly ) return true;
	  else 
	    {
	      std::vector<int> ac = calls.genotype(i).allele_list();
	      for (int i=0; i<ac.size(); i++) npoly.insert( ac[i] );
	      if ( npoly.size() > 1 ) return true;
	    }
	}
    }
  return false; // we never found a non-ref site (and also ref, if also_poly=T) in this file
}
