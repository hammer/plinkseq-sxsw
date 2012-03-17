
#include "svar.h"

#include "mask.h"
#include "genotype.h"
#include "vardb.h"
#include "gstore.h"

extern GStore * GP;

std::ostream & operator<<( std::ostream & out, const SampleVariant & v)
{   
  out    << v.ref << "/" << v.alt;  
  return out;
}


std::string SampleVariant::file_name() const 
{
  return GP->vardb.attached() ? GP->vardb.file_tag( fset ) : ".";
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


std::map<std::string,int> SampleVariant::genotype_counts( const affType & aff , const Variant * parent , bool phased ) const 
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
	c[ parent->consensus.label( *(parent->genotype( this , i ) ) , phased ) ] += count ? 1 : 0 ;
      else
	c[ label( *(parent->genotype( this , i ) ) , phased ) ] += count ? 1 : 0 ;

    }

  return c;
}



blob SampleVariant::encode_var_BLOB() const
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
  
  // TODO: Just add the filter (string) as a single element (for now)
  //  actually, we should change filter_info to be a parsed list

  pbVar.add_filter( filter_info );

  std::string s;
  pbVar.SerializeToString(&s);
  return blob(s);

}

   

blob SampleVariant::encode_vmeta_BLOB() const
{

  //
  // Set variant meta-data
  //
  

  VariantMetaBuffer pbVarMeta;
  
  std::vector<std::string> keys = meta.keys();
  
  for (unsigned int k=0; k<keys.size(); k++) 
    {
      
      VariantMetaUnit * m = pbVarMeta.add_vmeta();
      
      meta_index_t midx = MetaInformation<VarMeta>::field( keys[k] );
      
      m->set_name( keys[k] );
      
      // Number of elements in meta-value?
      
      int num = meta.size( keys[k] );
      
      switch ( midx.mt ) {
      case META_INT :
	{
	  m->set_type( VariantMetaUnit::INT );
	  const std::vector<int> * v = meta.ptr_int( midx.key );
	  for (int j=0; j<num; j++)
	    m->add_int_value( (*v)[j] );
	  break;
	}
      case META_FLOAT :
	{
	  m->set_type( VariantMetaUnit::FLOAT );
	  const std::vector<double> * v = meta.ptr_double( midx.key );
	  for (int j=0; j<num; j++)
	    m->add_double_value( (*v)[j] );
	  break;
	}
      case META_TEXT :
	{
	  m->set_type( VariantMetaUnit::TEXT );
	  const std::vector<std::string> * v = meta.ptr_string( midx.key );
	  for (int j=0; j<num; j++)
	    m->add_string_value( (*v)[j] );
	  break;
	}
      case META_BOOL :
	{
	  m->set_type( VariantMetaUnit::BOOL );
	  const std::vector<bool> * v = meta.ptr_bool( midx.key );
	  for (int j=0; j<num; j++)
	    m->add_bool_value( (*v)[j]  );
	  break;
	}
      case META_FLAG :
	{
	  m->set_type( VariantMetaUnit::BOOL );
	  // TODO -- shoudn't something else be set here??
	  m->add_bool_value( true );
	  break;
	}
      default :
	{
	  m->set_type( VariantMetaUnit::TEXT );
	  const std::vector<std::string> * v = meta.ptr_string( midx.key );
	  for (int j=0; j<num; j++)
	    m->add_string_value( (*v)[j] );
	  break;
	}
      }
      
    }
  

  // Return a BLOB 
  
  std::string s;
  pbVarMeta.SerializeToString(&s);
  return blob(s);

}



blob SampleVariant::encode_geno_BLOB() const
{ 
  
  //
  // Set all genotypes 
  //
  
  GenotypeBuffer pbGeno;
  unsigned int n = calls.size();  
  for ( unsigned int i = 0 ; i < n; i++)
    pbGeno.add_geno( calls.genotype(i).pack() );	

  // Return a BLOB   

  std::string s;
  pbGeno.SerializeToString(&s);
  return blob(s);

}


blob SampleVariant::encode_gmeta_BLOB() const
{

  unsigned int n = calls.size();
  
  //
  // Genotype meta-data
  //

  GenotypeMetaBuffer pbGMeta;
  

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
      
      GenotypeMetaUnit * m = pbGMeta.add_gmeta();

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
		    m->set_type( GenotypeMetaUnit::INT );
		    if ( ! constant_length )
		      m->add_len( meta.size( k->first ) );
		    const std::vector<int> * d = meta.ptr_int( mk_key );
		    for (unsigned int j = 0 ; j < d->size(); j++) 
		      m->add_int_value( (*d)[j] );
		    break;
		  }
		case META_FLOAT :
		  {
		    m->set_type( GenotypeMetaUnit::FLOAT );
		    if ( ! constant_length )
		      m->add_len( meta.size( k->first) );
		    const std::vector<double> * d = meta.ptr_double( mk_key );
		    for (unsigned int j = 0 ; j < d->size(); j++) 
		      m->add_double_value( (*d)[j] );
		    break;
		  }
		case META_TEXT :
		  {
		    m->set_type( GenotypeMetaUnit::TEXT );
		    if ( ! constant_length )
		      m->add_len( meta.size( k->first ) );
		    const std::vector<std::string> * d  = meta.ptr_string( mk_key );
		    for (unsigned int j = 0 ; j < d->size(); j++) 
		      m->add_string_value( (*d)[j] );
		    break;
		  }
		case META_BOOL :
		  {
		    m->set_type( GenotypeMetaUnit::BOOL );
		    if ( ! constant_length )
		      m->add_len( meta.size( k->first ) );
		    const std::vector<bool> * d = meta.ptr_bool( mk_key );
		    for (unsigned int j = 0 ; j < d->size(); j++) 
		      m->add_bool_value( (*d)[j] );
		    break;
		  }
		default : // add as text if unsure
		  {
		    m->set_type( GenotypeMetaUnit::TEXT );
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
  
  std::string s;  
  pbGMeta.SerializeToString(&s);    
  return blob(s);
  
}



void SampleVariant::store_BLOBs( blob * var_blob ,
				 blob * vmeta_blob ,
				 blob * geno_blob ,
				 blob * gmeta_blob )
{
  var_buf.ParseFromString( var_blob->get_string() );
  if ( vmeta_blob ) vmeta_buf.ParseFromString( vmeta_blob->get_string() );
  if ( geno_blob ) geno_buf.ParseFromString( geno_blob->get_string() );
  if ( gmeta_blob ) gmeta_buf.ParseFromString( gmeta_blob->get_string() );
}


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


  if ( svar->bcf || vcf_direct ) return true; 


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
  
  svar->alt = var_buf.alt();
  svar->ref = var_buf.ref();
  svar->qual = var_buf.quality();
  
  std::string my_filter = "";
  unsigned int num = var_buf.filter().size();
  for (int i=0;i<num; i++)    
    my_filter += var_buf.filter(i);  
  

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



void SampleVariant::decode_BLOB_alleles()
{ 
  
  // this function sometimes called when deciding to merge two SampleVariants into one Variant
  // (i.e. based on MERGE rules).  Thus we may need to peek at the REF and ALT prior to the 
  // decode_BLOB_basic() being called.
  
  if ( bcf || vcf_direct ) {
    return; // should already have been done
  }

  alt = var_buf.alt();
  ref = var_buf.ref();
}



bool SampleVariant::decode_BLOB_vmeta( Mask * mask, Variant * parent , SampleVariant * sample )
{    
  
  
  // SampleVariant meta-information
  
  //  1) to extract from a PB 
  
  //  2) to apply some variant-level masks 
  //     (if fail here, can save time by skipping genotype information extraction)
  
  
  // Possible that we do not need to look at any variant
  // meta-information, in which case, simply return now

  if ( mask && ! mask->load_variant_meta() ) return true;
  
  
  // For BCF-derived SVs, or those read direct from a VCF, skip the first step
  // i.e. genotypes already extracted
  
  if ( ! ( sample->bcf || sample->vcf_direct ) ) 
    {

      unsigned int num = vmeta_buf.vmeta().size();
  
      for (unsigned int k=0; k<num; k++)
	{	
	  
	  // Do we have permission to skip this field? 
	  
	  if ( mask && ! mask->load_variant_meta( vmeta_buf.vmeta(k).name() ) ) continue;
	  
	  // Any fields that are flagged in MetaMeta::pop_static() are 
	  // also copied to the consensus variant:
	  
	  bool incon = parent && MetaMeta::static_variant( vmeta_buf.vmeta(k).name() ) ; 
	  
	  switch ( vmeta_buf.vmeta(k).type() ) {
	    
	  case VariantMetaUnit::INT :
	    {
	      std::vector<int> d( vmeta_buf.vmeta(k).int_value().size() );
	      for ( unsigned int j = 0 ; j < d.size(); j++ ) 
		d[j] = vmeta_buf.vmeta(k).int_value(j);
	      sample->meta.set( vmeta_buf.vmeta(k).name() , d );
	      if ( incon ) parent->meta.set( vmeta_buf.vmeta(k).name() , d );
	      break;
	    }
	  case VariantMetaUnit::FLOAT :
	    {
	      std::vector<double> d( vmeta_buf.vmeta(k).double_value().size() );
	      for ( unsigned int j = 0 ; j < d.size(); j++ ) 
		d[j] = vmeta_buf.vmeta(k).double_value(j);
	      sample->meta.set( vmeta_buf.vmeta(k).name() , d );
	      if ( incon ) parent->meta.set( vmeta_buf.vmeta(k).name() , d );
	      break;
	    }
	  case VariantMetaUnit::BOOL :
	    {
	      // Stored in PB as BOOL, but this could either be a flag or a bool
	      
	      meta_index_t midx = MetaInformation<VarMeta>::field( vmeta_buf.vmeta(k).name() );
	      if ( midx.mt == META_BOOL ) 
		{
		  std::vector<bool> d( vmeta_buf.vmeta(k).bool_value().size() );
		  for ( unsigned int j = 0 ; j < d.size(); j++ ) 
		    d[j] = vmeta_buf.vmeta(k).bool_value(j);
		  sample->meta.set( vmeta_buf.vmeta(k).name() , d );
		  if ( incon ) parent->meta.set( vmeta_buf.vmeta(k).name() , d );
		}
	      else if ( midx.mt == META_FLAG ) 
		{
		  sample->meta.set( vmeta_buf.vmeta(k).name() );
		  if ( incon ) parent->meta.set( vmeta_buf.vmeta(k).name() );	      
		}
	      break;
	    }
	  default :
	    {
	      
	      //
	      // NOTE: legacy -- we now save FLAGs as BOOL in PB, so will only need the above
	      //
	      
	      meta_index_t midx = MetaInformation<VarMeta>::field( vmeta_buf.vmeta(k).name() );
	      if ( midx.mt == META_FLAG ) 
		{	      
		  sample->meta.set( vmeta_buf.vmeta(k).name() );
		  if ( incon ) parent->meta.set( vmeta_buf.vmeta(k).name() );	      
		}
	      else
		{
		  std::vector<std::string> d( vmeta_buf.vmeta(k).string_value().size() );
		  for ( unsigned int j = 0 ; j < d.size(); j++ ) 
		    d[j] = vmeta_buf.vmeta(k).string_value(j);
		  sample->meta.set( vmeta_buf.vmeta(k).name() , d );
		  if ( incon ) 
		    parent->meta.set( vmeta_buf.vmeta(k).name() , d );
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
      if ( ! mask->calc_filter_expression( *parent , *sample ) ) return false; 
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
  
    
  // Possible that we do not require any individual/genotypic level
  // information at all, in which case, return now.
    
  if ( mask && ! mask->load_genotype_data() ) 
    {
      target->calls.size( align ? align->size() : 0 );
      return true;
    }

  
  // Sample A : 1 2 3
  // Sample B : 4 5 6
  // Sample C : 3 4 

  // Consensus
  // Uniq      1(A1)  2(A2)                       5(B2)  6(B3)
  // Mult                    3(A3,C1)  4(B1,C2)
  //

  
  //
  // Decode genotype information, for 0+ individuals. If an optional
  // alignment is specified, use that to only extract that subset
  //
  
  if ( ! ( target->bcf || vcf_direct ) ) 
    {
      
      if ( mask == NULL || ( mask && mask->load_genotype_data() ) ) 
	{
	  
	  // Number of individuals in PBuffer
	  unsigned int n_buffer = geno_buf.geno().size();
	  
	  // Number of individuals we actually want
	  unsigned int n_variant = align ? align->size() : n_buffer ;
	  
	  //
	  // Allocate space as needed
	  //
	  
	  target->calls.size( n_variant );
      
  
	  //
	  // Add basic genotype information
	  //
	  
	  std::set<int> complex_genotype;

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

	      // Basic genotype
	      
	      if ( slot != -1 )        	      
		{	      
		  Genotype g;
		  if ( ! g.unpack( geno_buf.geno(i) ) )  complex_genotype.insert( i );
		  target->calls.add( g, slot );
		}	
	    }

	}

    
      //
      // Append genotype meta-information
      //
      
      unsigned int m = gmeta_buf.gmeta().size();
      
      for ( unsigned int k = 0 ; k < m ; k++ )
	{
	  
	  // Can we skip this?

	  if ( mask && ! mask->load_genotype_meta( gmeta_buf.gmeta(k).name() ) )
	    continue;
	  

	  // Does this have a set length, how are missing individuals
	  // handlded?
	  
	  // Mode: 
	  //  0) constant_length or no? 
	  
	  //  1) all_nonmissing -- no index, just list all values
	  //  2) indiv_index    -- specify IDs for people w/ data
	  //  3) missing_index  -- specify IDs for people w/out data
	  
	  bool constant_length = gmeta_buf.gmeta(k).has_fixed_len();
	  int  length = 0;
	  if ( constant_length ) 
	    length = gmeta_buf.gmeta(k).fixed_len();
	  
	  // Should only get one of these:
	  bool missing_index = gmeta_buf.gmeta(k).missing_index().size() > 0;
	  bool indiv_index = gmeta_buf.gmeta(k).indiv_index().size() > 0;
	  bool all_nonmissing = ! ( missing_index || indiv_index );
	  	  

	  // Set buffer sizes

	  // (N from genotype, not gmeta, buffer)
	  unsigned int n_buffer = geno_buf.geno().size();
	  
	  int nlen = indiv_index ?
	    (int)gmeta_buf.gmeta(k).indiv_index().size() :
	    n_buffer;
	  
	  int idx = 0;
	  
	  if ( all_nonmissing ) 
	    {
	      
	      switch (  gmeta_buf.gmeta(k).type() ) {
		
	      case GenotypeMetaUnit::INT :
		{
	      
		  for (int j=0;j<nlen; j++)
		    {
		      
		      idx = target->addIntGenMeta( j , source->fileset() , 
						   gmeta_buf, align, k, idx, 
						   constant_length ? length : gmeta_buf.gmeta(k).len(j) );		  
		    }
		  break;
		}
	      case GenotypeMetaUnit::FLOAT :
		{
		  for (int j=0;j<nlen; j++)
		    {
		      idx = target->addFloatGenMeta( j,  source->fileset() ,
						     gmeta_buf, align, k, idx, 
						     constant_length ? length : gmeta_buf.gmeta(k).len(j) );
		    }
		  break;
		}
	      case GenotypeMetaUnit::BOOL :
		{
		  for (int j=0;j<nlen; j++)
		    {
		      idx = target->addBoolGenMeta( j, source->fileset() ,
						    gmeta_buf, align, k, idx, 
						    constant_length ? length : gmeta_buf.gmeta(k).len(j) );
		    }		   
		  break;
		}
	      default :
		{
		  for ( int j=0;j<nlen; j++)
		    {
		      idx = target->addStringGenMeta( j, source->fileset() ,
						      gmeta_buf, align, k, idx, 
						      constant_length ? length : gmeta_buf.gmeta(k).len(j) );
		    }
		  
		}
	      }
	    }
	  
	  //
	  // Or using an individual-index?
	  //
	  
	  else if ( indiv_index )
	    {
	      switch (  gmeta_buf.gmeta(k).type() ) {
		
	      case GenotypeMetaUnit::INT :
		{
		  for ( int j=0;j<nlen; j++)
		    {
		      idx = target->addIntGenMeta( gmeta_buf.gmeta(k).indiv_index(j) , source->fileset() , 
						   gmeta_buf, align, k, idx, 
						   constant_length ? length : gmeta_buf.gmeta(k).len(j) );
		    }
		  break;
		}
	      case GenotypeMetaUnit::FLOAT :
		{
		  for ( int j=0;j<nlen; j++)
		    {
		      idx = target->addFloatGenMeta( gmeta_buf.gmeta(k).indiv_index(j) , source->fileset() , 
						     gmeta_buf, align, k, idx, 
						     constant_length ? length : gmeta_buf.gmeta(k).len(j) );
		    }
		  break;
		}
	      case GenotypeMetaUnit::BOOL :
		{
		  for ( int j=0;j<nlen; j++)
		    {
		      idx = target->addBoolGenMeta( gmeta_buf.gmeta(k).indiv_index(j) ,  source->fileset() ,
						    gmeta_buf, align, k, idx, 
						    constant_length ? length : gmeta_buf.gmeta(k).len(j) );
		    }
		  break;
		}
	      default :
		{		
		  for ( int j=0;j<nlen; j++)
		    {
		      idx = target->addStringGenMeta( gmeta_buf.gmeta(k).indiv_index(j),  source->fileset() ,
						      gmeta_buf, align, k, idx, 
						      constant_length ? length : gmeta_buf.gmeta(k).len(j) );
		    }
		}
		
	      }
	      
	    }
	  
	  
	  else // assume missing index
	    {
	      
	      // Missing index 
	      int skip = gmeta_buf.gmeta(k).missing_index(0);
	      int scnt = 0;
	      int cnt = 0;
	      
	      switch (  gmeta_buf.gmeta(k).type() ) {
		
	      case GenotypeMetaUnit::INT :
		{		    
		  for (int j=0;j<nlen; j++)
		    {
		      if ( j == skip )
			{
			  ++scnt;
			  skip = scnt < gmeta_buf.gmeta(k).missing_index().size() ? 
			    gmeta_buf.gmeta(k).missing_index(scnt) : -1 ;
			}
		      else
			{
			  idx = target->addIntGenMeta( j, source->fileset() , 
						       gmeta_buf, align, k, idx, 
						       constant_length ? length : gmeta_buf.gmeta(k).len(cnt++) );
			}
		    }
		  break;
		}
	      case GenotypeMetaUnit::FLOAT :
		{
		  for (int j=0;j<nlen; j++)
		    {
		      if ( j == skip )
			{
			  ++scnt;
			  skip = scnt < gmeta_buf.gmeta(k).missing_index().size() ? 
			    gmeta_buf.gmeta(k).missing_index(scnt) : -1 ;
			}
		      else
			{
			  idx = target->addFloatGenMeta( j,  source->fileset() ,
							 gmeta_buf, align, k, idx, 
							 constant_length ? length : gmeta_buf.gmeta(k).len(cnt++) );
			}
		    }
		  break;
		}
	      case GenotypeMetaUnit::BOOL :
		{
		  
		  for (int j=0;j<nlen; j++)
		    {
		      if ( j == skip )
			{
			  ++scnt;
			  skip = scnt < gmeta_buf.gmeta(k).missing_index().size() ? 
			    gmeta_buf.gmeta(k).missing_index(scnt) : -1 ;
			}
		      else
			{
			  idx = target->addBoolGenMeta( j,  source->fileset() ,
							gmeta_buf, align, k, idx, 
							constant_length ? length : gmeta_buf.gmeta(k).len(cnt++) );
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
			  skip = scnt < gmeta_buf.gmeta(k).missing_index().size() ? 
			    gmeta_buf.gmeta(k).missing_index(scnt) : -1 ;
			}
		      else
			{
			  idx = target->addStringGenMeta( j ,  source->fileset() ,
							  gmeta_buf, align, k, idx, 
							  constant_length ? length : gmeta_buf.gmeta(k).len(cnt++) );
			  
			}
		    }
		  
		}
		
	      }
	      
	    }
	  	  
	} // Next individual/genotype
                
    } // end of extracting all genotypic (meta) information from PB/BLOB



  //
  // Extract genotype information from a BCF-encoded buffer
  //

  else if ( target->bcf ) 
    {     
      
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

      // TODO: fix this really slow way to get allele count, that duplicates previous effort...
      
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
	  
	  const std::string & tag = format[t];
	  
	  // this should exist in the BCF, or else we would have received an error by now
	  
	  BCF::bcf_meta_t bt = target->bcf->bcftype[ tag ];
	  
	  int nalt = alt.size() + 1;
	  int ngen = (int) (nalt * (nalt+1) * 0.5);
	  if ( bt.len == -1 ) bt.len = nalt - 1;
	  else if ( bt.len == -2 ) bt.len = nalt;
	  else if ( bt.len == -3 ) bt.len = ngen;
	  
	  
	  // Unpack
	  
	  if ( bt.type == BCF::BCF_genotype )
	    {

	      for ( int i=0; i < n_buffer ; i++ )
		{
		  if ( s2t[i] != -1 )
		    {			
		      target->calls.genotype( s2t[i] ).bcf( target->bcf_genotype_buf[ p ] );		      
		    }
		  ++p;
		}
	    }

	  
	  else if ( bt.type == BCF::BCF_int32 )
	    {

	      for ( int i=0; i < n_buffer ; i++ )
		{
		  if ( s2t[i] != -1 ) 
		    {
		      std::vector<int> tmp( bt.len );
		      for (int j=0;j< bt.len; j++)
			tmp[j] = target->bcf_genotype_buf[ p + j * sizeof(uint32_t) ];
		      target->calls.genotype( s2t[i] ).meta.set( tag , tmp );		      
		    }
		  p += bt.len * sizeof(uint32_t);
		}
	    }

	  else if ( bt.type == BCF::BCF_uint8 )
	    {
	      for ( int i=0; i < n_buffer ; i++ )
		{
		  if ( s2t[i] != -1 ) 
		    {
		      std::vector<int> tmp( bt.len );
		      for (int j=0;j< bt.len; j++)
			tmp[j] = target->bcf_genotype_buf[ p + j * sizeof(uint8_t) ];
		      target->calls.genotype( s2t[i] ).meta.set( tag , tmp );
		    }
		  p += bt.len * sizeof(uint8_t);
		}
	    }
	  
	  else if ( bt.type == BCF::BCF_uint16 )
	    {
	      for ( int i=0; i < n_buffer ; i++ )
		{
		  if ( s2t[i] != -1 ) 
		    {
		      std::vector<int> tmp( bt.len );
		      for (int j=0;j< bt.len; j++)
			tmp[j] = target->bcf_genotype_buf[ p + j * sizeof(uint16_t) ];
		      target->calls.genotype( s2t[i] ).meta.set( tag , tmp );
		    }
		  p += bt.len * sizeof(uint16_t);
		}
	    }
	  

	  else if ( bt.type == BCF::BCF_double )
	    {
	      for ( int i=0; i < n_buffer ; i++ )
		{
		  if ( s2t[i] != -1 ) 
		    {
		      std::vector<double> tmp( bt.len );
		      for (int j=0;j< bt.len; j++)
			tmp[j] = target->bcf_genotype_buf[ p + j * sizeof(double) ];
		      target->calls.genotype( s2t[i] ).meta.set( tag , tmp );
		    }
		  p += bt.len * sizeof(double);
		}
	    }


	  else if ( bt.type == BCF::BCF_float )
	    {
	      for ( int i=0; i < n_buffer ; i++ )
		{
		  if ( s2t[i] != -1 ) 
		    {
		      std::vector<double> tmp( bt.len );
		      for (int j=0;j< bt.len; j++)
			tmp[j] = target->bcf_genotype_buf[ p + j * sizeof(float) ];
		      target->calls.genotype( s2t[i] ).meta.set( tag , tmp );
		    }
		  p += bt.len * sizeof(float);
		}
	    }

	  else if ( bt.type == BCF::BCF_string )
	    {
	      Helper::halt("BCF_string parsing not implemented yet for FORMAT");

	      for ( int i=0; i < n_buffer ; i++ )
		{
		  if ( s2t[i] != -1 ) 
		    {
		      std::vector<std::string> tmp( bt.len );
		      for (int j=0;j< bt.len; j++)
			tmp[j] = target->bcf_genotype_buf[ p + j * tmp[j].size() ];
		      target->calls.genotype( s2t[i] ).meta.set( tag , tmp );
		    }
		  p += bt.len * sizeof(float);
		}
	    }
	  
	} // next tag
      
    }
  
  
  //
  // Reading directly from a VCF buffer
  //

  else if ( vcf_direct )    
    {
      
      // If not a valid variant, do not try to expand genotypes
      if ( ! parent->valid() ) return false;
      
      // Individual counter
      int j = 0;
      

      unsigned int n_variant = align ? align->size() : vcf_direct_buffer.size()-9 ;

      // Target will always be consensus when reading a single VCF from the command line
      //  -- but not necessarily if reading BGZF-compressed VCFs that are indexed in thre VARDB
      // In any case, target *should* be set already to the appropriate place.

      target->calls.size( n_variant );
      
      // Call genotypes, add to variant   

      // check # of allowable alleles
      // TODO: needless re-parsing of the string, can speed this up if needed
      
      int na = 1;
      if ( vtarget->alternate() != "." ) 
	{ Helper::char_tok( vtarget->alternate() , &na , ',' ); ++na; } 
      
      
      for ( int i=9; i < vcf_direct_buffer.size(); i++)
	{

	  int slot = j;

	  if ( align )
	    {
	      // is this needed?? -- don't think both are needed 
	      slot = align->sample_remapping( source->fileset()  , j ) ;
	      if ( align->flat() ) slot = align->get_slot( source->fileset()  , slot );
	    }	  
	  
	  if ( slot != -1 )
	    {	      	    

	      Genotype g( vcf_direct_buffer(i) , vcf_gt_field , *vcf_formats , na );
	      	     
	      if ( mask && mask->fixxy() ) 
		{
		  
		    // Is this a flagged chromomse?
		    ploidy_t p = mask->ploidy( Helper::chrCode( parent->chromosome() ) );
		    
		    // AA --> A  | AB --> .  | ./. --> .
		    
		    if ( p == PLOIDY_HAPLOID )
		    {
			g.make_haploid();
		    }
		    else if ( p == PLOIDY_X  
			      && align->ind( slot )->sex() != FEMALE 			    
			      && ! mask->pseudo_autosomal( *parent ) )
		    {
			g.make_haploid();
		    }
		    else if ( p == PLOIDY_Y )
		    {
			if ( align->ind( slot )->sex() != MALE )
			    g.null( true );
			else if ( ! mask->pseudo_autosomal( *parent ) ) 
			    g.make_haploid();
		    }	  
		}
		
		target->calls.add( g , slot );
	    }
	  
	  ++j;
	}
      
      // Clear buffer
      vcf_direct_buffer.clear();
      
    }


  

  //
  // Now that we've extracted all the genotypic information and 
  // meta-information: do we want to apply any genotype masks? 
  //
  
  const int n_var = target->calls.size();

  
  // 1) mask 'assume-ref'.  Assume missing genotypes are
  // obligatorarily homozygous for the reference allele; also populate
  // PL fields to indicate this.
  
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
  
  
  //
  // 2) mask any per-individual segments
  //

  if ( mask && mask->genotype_segmask() ) 
    {
      
      Region vreg( *parent );
      
      for ( unsigned int i = 0 ; i < n_var ; i++ )
	{
	  if ( ! mask->eval_segmask( i , vreg ) )
	    {
	      target->calls.genotype(i).null( true );
	    }
// 	  else
// 	    {
// 	      // std::cout << "DEBUG: keeping " << i << " " << *parent << " " << GP->indmap(i)->id() << "\n";
// 	    }
	}
    }
  
  //
  // 3). mask 'geno' conditions -- i.e. determineing whether to zero-out specific genotypes
  // because they meet/fail to meet certain criteria, e.g. geno=DP:ge:10
  //
  
  if ( mask && mask->genotype_mask() ) 
    {
      for ( unsigned int i = 0 ; i < n_var ; i++ )
	{
	  if ( ! mask->eval( target->calls.genotype(i) ) )
	    target->calls.genotype(i).null( true );
	}
    }
  
  return true;
  
}


void SampleVariant::info( const std::string & s , VarDBase * vardb , int file_id , Variant * parent ) 
{
  

  // store original text 
  other_info = s; 
  
  if ( s == "." ) return;
  
  // parse semi-colon delimited list; true means escape quotes
  int ntok;
  Helper::char_tok f( s , &ntok , ';' , true );

  std::vector<Helper::char_tok*> ptoks;

  for (int i=0; i<f.size(); i++)
  {      
      
    int ntok2;
    Helper::char_tok * k = new Helper::char_tok( f(i) , &ntok2 , '=' , true );
    ptoks.push_back( k );
    
    mType mt = MetaInformation<VarMeta>::type( (*k)(0) );
      
      if ( mt == META_UNDEFINED ) 
      {	  
	  MetaInformation<VarMeta>::field( (*k)(0) , k->size() > 1 ? META_TEXT : META_FLAG , 1 , "undeclared tag" );	  
	  if ( vardb ) 
	      vardb->insert_metatype( file_id , (*k)(0) , k->size() > 1 ? META_TEXT : META_FLAG , 1 , META_GROUP_VAR , "undeclared tag" );
	  plog.warn("undefined INFO field (absent in VCF header)", (*k)(0) );
      }
            
  }

  // parse on ; separators into MetaInformation<VarMeta>   
  // 3rd arg indicates that unknown fields should be 
  // accepted (as a flag, or string value)
  
  //  -- although, given the above fix, this should no longer
  //     be necessary
  
  if ( MetaMeta::force_consensus() )
      parent->consensus.meta.parse( s , ';' , true );
  else
      meta.parse( s , ';' , true ); 
  
  // Now we've added to svar, if the parent pointer is non-null, check whether
  // or not any tags are static, in which case we also want to add this value
  // to the parent
  
  // add any STATIC tags to parent 
  // or apply the 'force consensus' mode
  
  if ( parent ) 
  {
      for (int i=0; i<ptoks.size(); i++)
      {	  
	  if ( MetaMeta::static_variant( (*ptoks[i])(0) ) )
	  {	      	      
	      const char * str = (*ptoks[i])(0);
	      
	      mType mt = MetaInformation<VarMeta>::type( str );

	      if      ( mt == META_UNDEFINED ) Helper::halt( "internal error in SampleVariant::info()" );

	      if      ( mt == META_FLAG )      parent->meta.set( str );
	      else if ( mt == META_INT )       parent->meta.set( str , meta.get_int( str ) );
	      else if ( mt == META_FLOAT )     parent->meta.set( str , meta.get_double( str ) );
	      else if ( mt == META_TEXT )      parent->meta.set( str , meta.get_string( str ) );
	      else if ( mt == META_BOOL )      parent->meta.set( str , meta.get_bool( str ) );
	  }
	  
	  delete ptoks[i];
      }
  }

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
  
  // not sure what this is doing now.. where/when does break-down of 
  //  'parse_alleles()'  get done??
  
//   VariantSpec * ps = SampleVariant::decoder.decode( "GT " + ref + " " + alt );
//   VariantSpec::set_format( 0 , NULL ); // 0 just means GT is first field (e.g. if GT:DP:GL etc)
//   specification( ps ); 
  
  // TODO: revisit this, but for now a "safe" way to determine this is whether 
  //       alt has a comma in it (i.e. >1 alternate allele specified)
  
  // simple = alt.find(",") == std::string::npos;    

}


std::string SampleVariant::allele1_label( const Genotype & g ) const
{
  if ( g.null() ) return ".";
  return alleles[ g.acode1() ].name();
}

std::string SampleVariant::allele2_label( const Genotype & g ) const
{
  if ( g.null() ) return ".";
  if ( g.haploid() ) return ".";
  return alleles[ g.acode2() ].name();
}

std::string SampleVariant::label( const Genotype & g , bool phased ) const
{

  std::stringstream s; 
  
  const std::string allele_delim = phased ? ( g.phased() ? "|" : "/" ) : "/" ; 

  if ( g.null() )
    {
      s << ( g.haploid() ? "." : "." + allele_delim + "." );
    }  
  else 
    {
      if ( g.haploid() )
	{
	  s << alleles[ g.acode1() ].name() ;
	}
      else
	{
	  std::string a1 = alleles[ g.acode1() ].name();
	  std::string a2 = alleles[ g.acode2() ].name();
	  if ( (!phased) && a1 > a2 ) { std::string t=a1;a1=a2; a2=t; }		  
	  s << a1 << allele_delim << a2;	  
	}    
    }
  return s.str();
}


std::string SampleVariant::num_label( const Genotype & g ) const
{
  std::stringstream s; 
  s << g;
  return s.str();
}


void SampleVariant::collapse_alternates( const Variant * parent , int altcode )
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
	      int ac = g.allele_count( );
	      g.set_alternate_allele_count( ac );
	    }
	}
      
    }
      
}


bool SampleVariant::has_nonreference( const bool also_poly , const std::vector<int> * imask ) const
{
  
  // If remap is NON-NULL, this means we only want to look at these IDs 
  // i.e. if this is consensus in a flat alignment, but we want a file-specific summary
  
  bool nonref = false;
  
  std::set<int> npoly;

  const int n = imask ? imask->size() : calls.size();

  for ( int i = 0; i < n; i++)
    {      

      const int idx = imask ? (*imask)[i] : i ;
      
      if ( calls.genotype( idx ).nonreference() ) 
	{	  
	  
	  // leave when hit first non-reference call
	  if ( ! also_poly ) 
	    return true;
	  else 
	    {
	      std::vector<int> ac = calls.genotype( idx ).allele_list( alleles.size() );
	      for (int a=0; a<ac.size(); a++) npoly.insert( ac[a] );
	      if ( npoly.size() > 1 ) return true; 
	    }
	}
    }
  return false; // we never found a non-ref site (and also ref, if also_poly=T) in this file
}



void SampleVariant::recall( Genotype & g , SampleVariant * p )
{
  
  // null genotypes always remain null

  if ( g.null() ) return;

  std::string str1 = p->alleles[ g.acode1() ].name();
  
  int a0 = 0;
  
  for ( int a=0; a<alleles.size(); a++)
    {
      if ( alleles[a].name() == str1 ) 
	{
	  if ( g.haploid() ) 
	    {
	      g.genotype( a );
	      return;
	    }
	  a0 = a;
	  break;
	}
    }
  
  
  // Second swap allele

  std::string str2 = p->alleles[ g.acode2() ].name();
  
  for ( int a=0; a<alleles.size(); a++)
    {
      if ( alleles[a].name() == str2 ) 
	{
	  g.genotype( a0 , a );
	  break;
	}
    }
  
}



  
inline int SampleVariant::addIntGenMeta( int j , int f , 
					 const GenotypeMetaBuffer & v, 
					 IndividualMap * align, 
					 int k,   // meta-info slot
					 int idx, // current counter 
					 int l ) // length arg 
{
    
    if ( align )
    {
	j = align->sample_remapping( f , j);	  
	if ( align->flat() ) j = align->get_slot( f , j );
    }      
    
    if( j == -1 ) return idx + l;     
    
    MetaInformation<GenMeta> & gmeta = calls.genotype( j ).meta ;
    
    if ( l == 1 ) 
	gmeta.set( v.gmeta(k).name() , v.gmeta(k).int_value( idx++ ) );
    else
    {
	std::vector<int> t(l);
	for ( int i = 0 ; i < l; i++) t[i] = v.gmeta(k).int_value( idx++ );
	gmeta.set( v.gmeta(k).name() , t );
    }      
    return idx;      
}



inline int SampleVariant::addFloatGenMeta( int j , int f , 
					   const GenotypeMetaBuffer & v, 
					   IndividualMap * align, 
					   int k,   // meta-info slot
					   int idx, // current counter 		     
					   int l )  // length arg
{
  
  if ( align )
    {
      j = align->sample_remapping( f , j);	  
      if ( align->flat() ) j = align->get_slot( f , j );
    }
  
  if( j == -1 ) return idx + l;
  
  MetaInformation<GenMeta> & gmeta = calls.genotype( j ).meta ;
  
  if ( l == 1 ) 
    gmeta.set( v.gmeta(k).name() , v.gmeta(k).double_value( idx++ ) );
  else
    {
      std::vector<double> t(l);
      for ( int i = 0 ; i < l; i++) t[i] = v.gmeta(k).double_value( idx++ );
      gmeta.set( v.gmeta(k).name() , t );
    }
  return idx;      
}


inline int SampleVariant::addStringGenMeta( int j , int f , 
					    const GenotypeMetaBuffer & v, 
					    IndividualMap * align, 
					    int k,   // meta-info slot
					    int idx, // current counter 		     
					    int l )  // length arg
{      
  if ( align ) 
    {
      j = align->sample_remapping( f , j);
      if ( align->flat() ) j = align->get_slot( f , j );
    }
  
  if( j == -1 ) return idx + l;
  
  MetaInformation<GenMeta> & gmeta = calls.genotype( j ).meta ;      
  
  if ( l == 1 )  
    gmeta.set( v.gmeta(k).name() , v.gmeta(k).string_value( idx++ ) );
  else
    {
      std::vector<std::string> t(l);
      for ( int i = 0 ; i < l; i++) t[i] = v.gmeta(k).string_value( idx++ );
      gmeta.set( v.gmeta(k).name() , t );
    }  
  
  return idx;  
}


inline int SampleVariant::addBoolGenMeta( int j , int f , 
					  const GenotypeMetaBuffer & v, 
					  IndividualMap * align, 
					  int k,   // meta-info slot
					  int idx, // current counter 		     
					  int l )  // length arg
{
  
  if ( align )
    {
      j = align->sample_remapping( f , j);	  
      if ( align->flat() ) j = align->get_slot( f , j );
    }

  if( j == -1 ) return idx + l;

  MetaInformation<GenMeta> & gmeta = calls.genotype( j ).meta ;
 
 if ( l == 1 ) 
   gmeta.set( v.gmeta(k).name() , v.gmeta(k).bool_value( idx++ ) );
 else
   {
     std::vector<bool> t(l);
     for ( int i = 0 ; i < l; i++) t[i] = v.gmeta(k).bool_value( idx++ );
     gmeta.set( v.gmeta(k).name() , t );
   }
 
 return idx;
}  


  /// pretty print versions of the above
std::string SampleVariant::pp_reference() const
{
  if ( ref.size() < 10 ) return ref;
  return ref.substr(0,5) + "...(" + Helper::int2str( ref.size() ) + "bp)";
}


std::string SampleVariant::pp_alternate() const
{
  if ( alt.size() < 10 ) return alt;
  return alt.substr(0,5) + "...(" + Helper::int2str( alt.size() ) + "bp)";
}
