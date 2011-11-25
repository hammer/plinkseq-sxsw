#include <memory>

#include "variant.h"
#include "genotype.h"
#include "helper.h"
#include "mask.h"
#include "bcf.h"

#include "gstore.h"

extern GStore * GP;


Variant::Variant( bool b ) 
{ 
  init(); 
  is_valid = b; 
}
  

Variant::Variant( const std::string & n, int c, int b )
{		
  init();      
  vname = n;
  chr = c;
  bp = bp2 = b;      
}

void Variant::init()
{
  chr = bp = bp2 = 0;    
  vname = ".";
  meta.clear();
  align = NULL;
  is_valid = true;
  is_multi_sample = false;
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

  // Note -- we may have a unflat, single-sample alignment, if we
  // are letting the different SVARs be compiled within the same
  // VCF as a single variant.  In that case, the force_unflat()
  // flag will be set and we can use that.
  
  if ( ( ! align->multi_sample() ) && align->flat() ) 
    {
      
      int n_alleles = consensus.parse_alleles();

      // for biallelic markers, we can leave now
      if ( n_alleles == 2 ) 
	{	  
	  return true;
	}
      
      consensus.set_allelic_encoding();
            
      return true;
    }
  
  
  for (int i = 0 ; i < svar.size(); i++ )
    svar[i].parse_alleles();



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

	      //need_to_resolve = true;
	      //obs_alleles.insert( sv.ref );

	      // should probably just halt?

	      plog.warn(  " **serious** incompatible REF sequences " , coordinate() ); 

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
  
  if ( align->flat() && ! need_to_resolve ) 
    {
      return true;
    }
  

  //
  // Handle flat case, where ALT alleles are different between samples -- will need to 
  // recode the genotypes of some individuals
  //

  
  if ( align->flat() )
    {
      const int n = align->size();
      for (int i=0; i<n; i++)
	{
	  
	  int2 j = align->unique_mapping(i);      
	  
	  // Does this individual feature in multiple samples? 
	  
	  bool multiple_samples = j.p1 != -1 ? ftosv[ j.p1 ].size() > 1 : true;
	  
	  if ( ! multiple_samples )
	    { 
	      
	      // not sure this would ever fail, but just in case (perhaps
	      // also need to check ftosv/j.p1 above for range?)
	      
	      if ( ftosv[ j.p1 ].size() == 1 ) 
		{
		  SampleVariant * p = psample( ftosv[ j.p1 ][0] );	      
		  if ( p )  
		    {		  
		      Genotype & g = consensus(i);

// 		      // Get ACGT label
// 		      std::string label = p->label( g ); // get ACGT label

		      // Recall with consensus encoding
		      consensus.recall( g , p );
		      
		    }	 
		}
	    }
	  else
	    Helper::halt("internal error in Variant::make_consensus()");
	}
      return true;
    }

      
  //
  // Genotypes in non-flat alignments, (which might need to be recoded)
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
	  
	  // not sure this would ever fail, but just in case (perhaps
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
		      //std::string label = p->label( g );		  
		      //g = consensus.spec->callGenotype( label , this , true );  //T=ACGT mode
		      consensus.recall( g , p );		      
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
// 		      std::string label = p->label( g );
// 		      g = consensus.spec->callGenotype( label , this , true ); // T=ACGT mode
		      consensus.recall( g , p );
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
		  
		  if ( ! Genotype::equivalent ( g , consensus(i) ) )
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

  std::stringstream s; 
  
  if ( ! biallelic() ) return "0" + delim + "0";
  
  if ( g.more() ) // do not display multi-allelic variants (but should not happen, given above)
    s << "0" << delim << "0";
  else
    {
      if( g.null() )
	s << "0" << delim << "0";      
      else
	{
	  
	  if ( g.acode1() ) 
	    s << consensus.alt; else s << consensus.ref;
	  
	  s << delim;
	  
	  if ( g.haploid() ) // Haploid --> homozygote in PED files...
	    {
	      if ( g.acode1() ) 
		s << consensus.alt; else s << consensus.ref;		  
	    }
	  else	    
	    {
	      if ( g.acode2() ) 
		s << consensus.alt; else s << consensus.ref;
	    }
	}
    }
  
  return s.str();
}



std::string Variant::geno_label( const Genotype & g ) const
{
  return consensus.label( g );
}

std::string Variant::geno_label( const int s , const Genotype & g ) const
{
  return svar[s].label( g );
}


std::string Variant::phased_geno_label( const Genotype & g ) const
{
  return consensus.label( g , true );
}

std::string Variant::phased_geno_label( const int s , const Genotype & g ) const
{
  return svar[s].label( g , true );
}


std::string Variant::label( const int i  , const std::string & delim ) const
{   
  
  // get basic textual representation of genotype from SampleVariant function
  // (w/ phase shown)

  std::string s = consensus.label( consensus(i) , true );
  
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
	  const SampleVariant * svar = psample( j->first );
	  if ( svar ) 
	    s +=  ( j != gm.begin() ? delim : "" ) + svar->label( *(j->second) , true );	    
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

  if ( flat() && ! infile_overlap() ) return ss.str();

  std::map<int, const Genotype *> gm = all_genotype(i);	  
  std::map<int, const Genotype *>::iterator j = gm.begin();
  if ( gm.size() > 1 )
    {
      ss << " {";
      j = gm.begin();	  
      while ( j != gm.end() )
	{
	  const SampleVariant * svar = psample( j->first );
	  if ( svar ) ss << ( j != gm.begin() ? delim : "" ) << j->second->meta ;
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



int Variant::n_alleles() const
{
  return consensus.alleles.size();    
}

const Allele & Variant::allele(const int n) const
{
  return consensus.alleles[n];
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
    if ( ! genotype(i)->null() ) ++n;
  return n;
}


bool Variant::n_minor_allele( int * m_ , int * n_ , double * maf_ , const affType & aff ) const
{
  
  // returns counts of any non-reference allele
  // optionally, if aff is not UNKNOWN_PHE , split this by case/control
  
  int m = 0;
  int n = 0;
  
  if ( aff == UNKNOWN_PHE )
    {

      for (int i=0; i< size() ; i++) 
	{
	  const Genotype * g = genotype(i);	  
	  if ( ! g->null() )
	    {
	      m += g->allele_count();
	      n += g->copy_number();	    
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
	      if ( ! g->null() )
		{
		  m += g->allele_count();
		  n += g->copy_number();
		}		      
	    }
	}      
    }


  // If reference allele is minor allele, swap the minor allele count,
  // and let the caller know about this

  double maf = (double)m / (double)n ;
  
  bool altmin = true;
  
  if ( maf > 0.5 ) 
    {
      m = n - m;
      maf = 1 - maf;
      altmin = false;
    }
  
  // Update external values
  
  if ( m_ ) *m_ = m;
  if ( n_ ) *n_ = n;
  if ( maf_ ) *maf_ = maf; 

  // Alternate allele is minor allele
  return altmin;  
}



std::map<std::string,int> Variant::allele_counts( const affType & aff ) const 
{
  return consensus.allele_counts( aff, this );
}


std::map<std::string,int> Variant::genotype_counts( const SampleVariant & svar , const affType & aff , bool unphased ) const
{
  return sample_genotypes( svar ).genotype_counts( aff , this , unphased );
}

std::map<std::string,int> Variant::genotype_counts( const int si , const affType & aff , bool unphased ) const
{  
  return sample_genotypes( si ).genotype_counts( aff , this , unphased );
}

std::map<std::string,int> Variant::genotype_counts( const affType & aff , bool unphased ) const
{  
  return consensus.genotype_counts( aff , this , unphased );
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
  int n=0,m=0;
  n_minor_allele(&n,&m);
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
  if ( consensus.alt == "A" && ( consensus.ref == "C" || consensus.ref == "T" ) ) return true;
  if ( consensus.alt == "C" && ( consensus.ref == "A" || consensus.ref == "G" ) ) return true;
  if ( consensus.alt == "G" && ( consensus.ref == "C" || consensus.ref == "T" ) ) return true;
  if ( consensus.alt == "T" && ( consensus.ref == "A" || consensus.ref == "G" ) ) return true;

  // this means if we have complex, collapsed variants, we get ti == tv == false, i.e. not defined  
  return false;  
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


std::string Variant::VCF()
{
    
    // Construct a string that is a VCF format entry
    
    std::ostringstream s;
    
    // VCF is tab-delimited 
    
    s << Helper::chrCode( chr ) << "\t"
      << bp << "\t"
      << vname << "\t"
      << consensus.ref << "\t"
      << consensus.alt << "\t";
    
    if ( consensus.qual < 0 )
	s << "." << "\t";
    else
	s << consensus.qual << "\t";
    
    s << consensus.filter_info << "\t"
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
      for (unsigned int j=0; j<keys.size(); j++) if ( MetaMeta::display( keys[j] ) ) allkeys.insert( keys[j] );
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
	  bool altmin = n_minor_allele( &m , &c );	  	  
	  if ( ! mask->count_filter( m ) ) return false;	       
	}      

      if ( mask->frequency_filter() )
	{
	  double maf;
	  bool altmin = n_minor_allele( NULL , NULL , &maf );
	  if ( ! mask->frequency_filter( maf ) ) return false;
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
  int dira = n_minor_allele( &alta , &tota , NULL , CASE );
  int diru = n_minor_allele( &altu , &totu , NULL , CONTROL );

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

int Variant::size( const int svar_id ) const
{
  return svar_id == -1 ? size() : ( align ? align->size( svar[svar_id].fset ) : 0 ) ;
}

void Variant::resize(const int n ) 
{
  consensus.calls.size(n);
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
	    g.insert( std::make_pair( ii[jj], &(svar[ ii[jj] ]( k.p2 ) ) ) );
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
	    g.insert( std::make_pair( ii[jj], &(svar[ ii[jj] ]( k.p2 ) ) ) );
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

int Variant::unique_svar_slot( int file_id ) const
{
  std::map<int,std::vector<int> >::const_iterator i = ftosv.find( file_id );
  if ( i == ftosv.end() ) return -1;
  if ( i->second.size() != 1 ) return -1;
  return *(i->second.begin());
}


bool Variant::file_present( const int file_id ) const
{
  std::map<int,std::vector<int> >::const_iterator i = ftosv.find( file_id );
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

std::set<std::string> Variant::meta_filter( ) const
{  
  std::set<std::string> r;
  
  if ( align->single_sample() ) 
    {
      std::vector<std::string> f = consensus.filters();
      for (int i=0; i<f.size(); i++) r.insert(f[i]);
      return r;
    }
 
 for (int i=0; i < svar.size(); i++ )
    {
      std::vector<std::string> f = svar[i].filters();
      for (int i=0; i<f.size(); i++) r.insert(f[i]);           
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




bool Variant::concordant( int s1, 
			  const Genotype * g1, 
			  int s2, 
			  const Genotype * g2 ) const
{
  
  const SampleVariant * p1 = psample(s1);
  if ( ! p1 ) return true;
  
  const SampleVariant * p2 = psample(s2);
  if ( ! p2 ) return true;
  
  return concordant( p1, g1, p2, g2 );
}


bool Variant::concordant( const SampleVariant * s1,
			  const Genotype * g1, 
			  const SampleVariant * s2, 
			  const Genotype * g2 ) const
{
  if ( g1->null() || g2->null() ) return true;
  if ( *g1 == *g2 ) return true;
  std::string a1 = s1->label( *g1 , true );
  std::string a2 = s2->label( *g2 , true );
  return a1 == a2;
}


std::map<std::string,int> SampleVariant::allele_count( const int i ) const
{  
  const Genotype & g = calls.genotype(i);  
  std::map<std::string,int> a;  
  if ( g.null() ) return a;  
  a[ alleles[ g.acode1() ].name() ]++;
  if ( ! g.haploid() ) a[ alleles[ g.acode2() ].name() ]++;
  return a;
}



bool SampleVariant::align_reference_alleles( SampleVariant & s1 , 
					     SampleVariant & s2 , 
					     bool alt )
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


bool Variant::has_nonreference( const SampleVariant & svar ) const
{

  // Is this is consensus, just return
  if ( svar.fset == 0 ) return svar.has_nonreference( false );

  // If this is a specific sample variant, we might need to extract the subset of individuals
  // from the consensus.

  // 1) Where are genotypes stored for this sample? 
  
  SampleVariant & sv_geno = sample_genotypes( svar );
  
  // If in the SV itself (a non-flat alignment)  
  if ( sv_geno.fset ) 
    return sv_geno.has_nonreference( false );
  else // otherwise, if looking in the consensus, need a i-mask
    return sv_geno.has_nonreference( false , align->file2consensus( svar.fset ) );
}

std::vector<int> Variant::indiv_mask( const int file_id ) const
{
  const std::vector<int> * i = align->file2consensus( file_id );
  if ( i ) return *i;
  std::vector<int> j;
  return j;
}



bool Variant::has_nonreference_by_file( const int file_id ) const
{
  std::vector<const SampleVariant *> sv = fsample( file_id );
  for ( int s = 0 ; s < sv.size() ; s++ )
    if ( has_nonreference( *sv[s] ) ) return true;    
  return false;
}

std::string Variant::pp_reference() const
{ 
  if ( consensus.ref.size() < 10 ) return consensus.ref; 
  return consensus.ref.substr(0,5) + "...(" + Helper::int2str( consensus.ref.size() ) + "bp)";
}

std::string Variant::pp_alternate() const 
{ 
  if ( consensus.alt.size() < 10 ) return consensus.alt; 
  return consensus.alt.substr(0,5) + "...(" + Helper::int2str( consensus.alt.size() ) + "bp)";
}


