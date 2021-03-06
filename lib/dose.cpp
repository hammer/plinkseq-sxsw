#include "plinkseq/dose.h"

#include "plinkseq/gstore.h"
#include "plinkseq/helper.h"
#include "plinkseq/vardb.h"
#include "plinkseq/seqdb.h"

#include <fstream>
#include <cmath>

extern GStore g;


int DoseReader::read_fam( const int file_id ) 
{
  
  // expect simple format of just an ID list here
  
  if ( ! Helper::fileExists( fam_filename ) ) return 0;
  
  std::ifstream FAM( fam_filename.c_str() , std::ios::in );
  
  int ni = 0;
  while ( ! FAM.eof() )
    {
      std::string t;
      FAM >> t;
      if ( t == "" ) continue;
      ++ni;      
      // Register this individual in VARDB
      Individual ind(t);           
      vardb->insert( file_id , ind );      
    }

  if ( ni == 0 ) 
    {
      FAM.close();
      plog.warn("no individuals specified in " + fam_filename );
      return 0;
    }
  
  return ni;
}

void DoseReader::read_meta( const int file_id )
{

  if ( ! has_meta ) return;

  // read header -- expecting #meta1 meta2 meta2
  //                rsID m1 m2 m3
  
  InFile meta(meta_filename);
  std::string l = meta.readLine();
  if ( l.size() == 0 ) Helper::halt( "could not open meta-information file [" + meta_filename + "]" );
  if ( l[0] != '#' ) Helper::halt("expecting header row to start with '#' for " + meta_filename );
  metas = Helper::parse( l.substr(1) , "\t " );
  
  // register all as floats for scalar Variant-meta tags
  for (int m=0;m<metas.size();m++)
    {
      registerMetatype( metas[m] , META_FLOAT , 1, META_GROUP_VAR , "." );
      vardb->insert_metatype( file_id , metas[m] , META_FLOAT, 1, META_GROUP_VAR , "." );
    }
  meta.close();  
}



bool DoseReader::read_dose( const std::string & f )
{
  


  const bool store_as_dosage = storage_format == "as-dosage" ;
  const bool store_as_prob3  = storage_format == "as-posteriors" ; 

  if ( store_as_dosage ) 
    MetaInformation<GenMeta>::field( tagname , META_FLOAT , 1 , "ALT dosage (0..2)" );
  else if ( store_as_prob3) 
    MetaInformation<GenMeta>::field( tagname , META_FLOAT , 3 , "Genotype probs (REF,HET,HOM)" );
  
  //
  // Open dosage file 
  //

  dose_filename = f;
  
  if ( ! Helper::fileExists( dose_filename ) ) return false;
  
  InFile inp_dosage( dose_filename );
  



  //
  // Map file?
  //

  InFile inp_map;
  if ( separate_map )
    {      
      Helper::checkFileExists( map_filename );
      inp_map.open( map_filename.c_str() );
    }

  //
  // If no IDs in the main doseage file, we require a separate map
  //

  if ( ! has_rsid )
    {
      if ( ! separate_map ) 
	Helper::halt( "must specify a separate map if --no-ids in primary dosage file" );
    }


  //
  // Get slot in VARDB, and track file-source, if we have not already seen this fileset
  //
  

  uint64_t file_id = vardb->lookup_file_id( filetag );   
  
  bool already_loaded_header = file_id != 0 ;
  
  if ( ! already_loaded_header )
    {
      file_id = vardb->insert( dose_filename , filetag );      
      vardb->insert_header( file_id , "format", "dosage-data" );
    }
  

  
  //
  // Do we need a separate meta-file?
  //
  
  InFile inp_meta;

  if ( has_meta )
    {
      read_meta( file_id );
      Helper::checkFileExists( meta_filename );
      inp_meta.open( meta_filename.c_str() );
      std::string skip_hdr = inp_meta.readLine();
    }

  //
  // Get individual map
  //
  
  int ni = 0;
  
  if ( separate_fam ) 
    {
      ni = read_fam( file_id );

      // and also skip any header row in the dosage file
      if ( skip_header ) 
	{
	  std::string ignored = inp_dosage.readLine();
	}
    }
  else // read from header
    {

      std::string line = inp_dosage.readLine();

      // expecting format #ID1 ID2 ID3
      
      if ( line.size() == 0 || 
	   line[0] != '#' ) 
	Helper::halt("expecting first row in dosage file to start '#' followed by ID list" );      
      
      std::vector<std::string> ids = Helper::parse( line.substr(1) , "\t " );
      ni = ids.size();      
      
      // insert people, if not already done
      if ( ! already_loaded_header ) 
	for (int i=0;i<ni;i++)
	  vardb->insert( file_id , Individual( ids[i] ) );      	        
    }
  
  if ( ni == 0 ) Helper::halt( "no inividuals specified in either header or --indiv-file" );


  //
  // Expected number of fields in the dosage file
  //
  
  int per_term = 1;
  if      ( format_prob2 ) per_term = 2;
  else if ( format_prob3 ) per_term = 3;
  
  int expected_fields = ni * per_term ;
  if ( has_rsid ) ++expected_fields;
  if ( separate_map ) 
    {
      if ( map_pos_only ) expected_fields += 2;    // implies A1 and A2 in main dosage file
      if ( map_allele_only ) expected_fields += 2;  // implies chr/bp in main dosage file      
    }
  else expected_fields += 4; // implies chp, bp, A1 and A2 all in dosage file (in this order)

  const int start_pos = expected_fields - ni * per_term;
    
  
  //
  // Read each variant, one at a time, and save into VARDB
  //

  
  long int cnt = 0;

  while (1) 
    {
      
      // Read dosage line
      
      std::string dosage_line = inp_dosage.readLine();

      if ( inp_dosage.eof() ) break;
      if ( dosage_line == "" ) continue;

      //
      // Two options for parsing: use the quick char_tok(), although currently this 
      // only allows one single char as a delimiter, i.e. will choke on files that 
      // have tab and space, or multiple spaces, etc. 
      //
      
      //  default =                       tab   with char_tok
      //  --format space-delimited        space with char_tok
      //  --format whitespace-delimited   whitespace parse()
      //
      
      int ntok;
      
      
      Helper::char_tok toks( &(dosage_line[0]) , dosage_line.size() , &ntok , spaced ? ' ' : '\t' );


      //std::vector<std::string> toks = Helper::parse( dosage_line , " \t" ) ; 

      
      ntok = toks.size();

      // expecting rsID + 2(pos) + 2(alleles) + n_genotypes
            
      if ( ntok < expected_fields ) 
	Helper::halt( "invalid dosage line entry, wrong number of rows (expecting " 
		      + Helper::int2str( expected_fields ) + " but found " 
		      + Helper::int2str( ntok ) + " ):\n" + dosage_line );
      
      
      // Create a new variant

      Variant v;

     
            
      //
      // Assign MAP information
      //
      
      std::string id = "", chr = "", a1 = "" , a2 = "";
      int bp = 0;
      

      // are we reading from a separate MAP file?
      
      if ( separate_map ) 
	{
	  std::vector<std::string> mtok = inp_map.tokenizeLine("\t ");

	  if ( mtok.size() == 0 ) continue;
	  
	  if ( map_pos_only ) 
	    {

	      if ( mtok.size() != 3 ) 
		Helper::halt( "problem with map, expecting 3 fields:\n[" + Helper::stringize( mtok , "\t") + "]" );	      

	      // missing chr/bp code implies skip
	      if ( mtok[1] == "." || mtok[2] == "." ) continue;

	      id = mtok[0];
	      chr = mtok[1];
	      if ( ! Helper::str2int( mtok[2] , bp ) ) 
		Helper::halt( "problem with map line base position: " + mtok[2] );	      
	    }
	  else if ( map_allele_only ) 
	    {
	      if ( mtok.size() != 3 ) Helper::halt( "problem with map line, expecting 3 fields:\n[" + Helper::stringize( mtok , "\t" ) +"]" );
	      id = mtok[0];
	      a1 = mtok[1];
	      a2 = mtok[2];
	    }
	  else
	    {
	      if ( mtok.size() != 5 ) Helper::halt( "problem with map line, expecting five fields\n[" + Helper::stringize( mtok , "\t" ) + "]" );
	      // missing chr/bp code implies skip
	      if ( mtok[1] == "." || mtok[2] == "." ) continue;

	      id = mtok[0];
	      chr = mtok[1];
	      if ( ! Helper::str2int( mtok[2] , bp ) ) 
		Helper::halt( "problem with map line base position: " + mtok[2] );	      
	      a1 = mtok[3];
	      a2 = mtok[4];
	    }	  
	}
      

      // if no ID in the main file, skip that column;

      const int id_offset = has_rsid ? 0 : -1;      
      
      if ( separate_map )
	{
	  if ( map_pos_only ) // .dose contains a1,a2
	    {
	      a1 = toks[ 1 + id_offset ];
	      a2 = toks[ 2 + id_offset ];
	    }
	  else if ( map_allele_only ) // .dose contains chr/bp
	    {
	      // missing chr/bp code implies skip
	      if ( toks[ 1 + id_offset ] == "." || toks[ 2 + id_offset ] == "." ) continue;

	      chr = toks[ 1 + id_offset ];	      	      
	      if ( ! Helper::str2int( toks[ 2 + id_offset ] , bp ) ) Helper::halt( "invalid base-position: " + std::string( toks[ 2 + id_offset ] ) );	      
	    }
	}
      else
	{
	  // missing chr/bp code implies skip
	  if ( toks[ 1 + id_offset ] == "." || toks[ 2 + id_offset ] == "." ) continue;

	  chr = toks[ 1 + id_offset ];
	  if ( ! Helper::str2int( toks[ 2 + id_offset ] , bp ) ) Helper::halt( "invalid base-position: " + std::string( toks[ 2 + id_offset ] ) );
	  a1 = toks[ 3 + id_offset ];
	  a2 = toks[ 4 + id_offset ];	  
	}
      

      //
      // check that ID matches
      //

      if ( has_rsid ) 
	{
	  if ( separate_map )
	    {
	      if ( id != toks[0] ) 
		{
		  std::string tid = toks[0];
		  Helper::halt( "mismatch of ID codes: [" + tid + "] in dosage, but [" + id + "] in map-file" );
		}
	    }
	  else
	    {
	      // need to set ID from first column (required)
	      id = toks[0];
	    }	  
	}




      
      //
      // We seem to be all set
      //
      
      v.chromosome( Helper::chrCode( chr ) ) ;
      v.position( bp );
      v.name( id );
      
      // Use SEQDB to look up reference allele. Implies that all
      // alleles are known to be on the positive strand;
      
      bool ref1 = true;  // A1 is reference allele
      
      if ( a1 == "." || a1 == "0" ||  a1 == "X" )
	{
	  std::string t = a2;
	  a2 = a1;
	  a1 = t;
	}
      
      // TODO: this currently only works for SNPs
      // TODO: handle symbolic alleles (e.g. R, I, D)
      
      if ( seqdb )
	{

	  std::string ref = seqdb->lookup( v.chromosome() , bp );
	  
	  Helper::str2upper(ref);

	  if ( ref != "N" ) 
	    {
	      
	      if ( ref == a2 ) ref1 = false;
	      
	      else if ( ref != a1 ) 
		{
		  
		  // either attempt to fix strand, or give error
		  // only for biallelic, single nucleotide markers 
		  if ( 1 ) 
		    {
		      		      
 		      if      ( a1 == "A" ) a1 = "T";
 		      else if ( a1 == "C" ) a1 = "G";
 		      else if ( a1 == "G" ) a1 = "C";
 		      else if ( a1 == "T" ) a1 = "A";		      
		      
 		      if      ( a2 == "A" ) a2 = "T";
 		      else if ( a2 == "C" ) a2 = "G";
 		      else if ( a2 == "G" ) a2 = "C";
 		      else if ( a2 == "T" ) a2 = "A";

		      bool flipped = false;

 		      if      ( ref == a1 ) { flipped = true; ref1 = true; }
 		      else if ( ref == a2 ) { flipped = true; ref1 = false; }
 		      else 
 			plog.warn( "attempted to fix, but mismatching reference allele versus SEQDB" , 
 				   v.displaycore() + " " + a1+"/"+a2 + " vs " + ref );

		      if ( flipped )
			plog.warn( "mismatching reference allele versus SEQDB, but resolved by strand-flip : " , 
				   v.displaycore() + " " + a1+"/"+a2 + " vs " + ref );
		      
		      
 		    }
		}		
	    }
	}
      

      if ( ref1 ) 
	{	  
	  v.consensus.reference( a1 );
	  v.consensus.alternate( a2 );	  
	}
      else
	{	  
	  v.consensus.reference( a2 );
	  v.consensus.alternate( a1 );	  
	}
      

      if ( cnt % 1000 == 0  )
	plog.counter1( "inserted " + Helper::int2str( cnt ) + " variants" );

      
      //
      // Get genotypes
      //

      v.resize( ni );
      //      v.reserve( ni );

      
      // read 'start_pos' slots first      
	
      int p = start_pos;

      //
      // Re-use the same genotype
      //

      Genotype g;


      // dosages will always be present, so just create once upfront
      
      if ( store_as_dosage ) 
	g.meta.set( tagname , (double)2.0 );	 // note -- must explicitly be a double 
      else if ( store_as_prob3 )	
	g.meta.set( tagname , std::vector<double>(3 , 0.0 ) );
      
      // get a pointer to this dosage store
      meta_key_t mk = MetaInformation<GenMeta>::field( tagname ).key;
      double * dstore = g.meta.mutptr1_double( mk );

      for (int i=0; i<ni; i++)
	{
	  
	  //
	  // Get dosage
	  //
	  
	  int geno = 0; // 1,2,3 = ref,het,hom


	  //
	  // By default, if ref1==T we assume 
	  //   A1  A2    dosage is for 'A2' allele
	  //             posterior is 11 12 22
	  
	  // But always store as dosage of 'alternate' allele when possible
	  
	  //   dosage = EC = alternate allele (2)
	  //   PP     = 11 / 12 / 22  ( REF / HET / HOM )
	  
	  // Thus by default, we assume that A1 is reference
	  // If A2 is reference, then swap codes: (reverse PP / calculate 2 - dosage )
	  
	  //
	  // assume just biallelic for now.
	  //

	  if ( format_prob2 ) 
	    {
	      
	      std::vector<double> pp(3);
	      
	      if ( Helper::str2dbl( toks[p++] , pp[0] ) ) 
		{
		  if ( Helper::str2dbl( toks[p++] , pp[1] ) )
		    {		  		      
		      
		      pp[2] = 1 - pp[0] - pp[1];
		      
		      if ( ! ref1 ) 
			{
			  double tmp = pp[2];
			  pp[2] = pp[0];
			  pp[0] = tmp;
			}
		      
		      if ( store_as_dosage ) 
			{
			  // 0..2 scale, of '2/alternate' allele
			  g.meta.set( tagname , pp[1] + 2 * pp[2] );
			}
		      else if ( store_as_prob3 )
			{
			  g.meta.set( tagname , pp );
			}
		      
		      if ( make_hard_call ) 
			{
			  if      ( pp[0] >= hard_call_prob_threshold ) geno = 1;
			  else if ( pp[1] >= hard_call_prob_threshold ) geno = 2;
			  else if ( pp[2] >= hard_call_prob_threshold ) geno = 3;
			}
		    }
		}
	    } 
	  else if ( format_prob3 ) 
	    {
	      
	      std::vector<double> pp(3);
	      
	      if ( Helper::str2dbl( toks[p++] , pp[0] ) && 
		   Helper::str2dbl( toks[p++] , pp[1] ) && 
		   Helper::str2dbl( toks[p++] , pp[2] ) )
		{		  
		  
		  if ( ! ref1 )
		    {
		      double tmp = pp[2];
		      pp[2] = pp[0];
		      pp[0] = tmp;
		    }		  
		  
		  if ( store_as_dosage ) 
		    {
		      // 0..2 scale, of 'B' allele
		      g.meta.set( tagname , pp[1] + 2 * pp[2] );
		    }
		  else if ( store_as_prob3 )
		    {		      
		      g.meta.set( tagname , pp );
		    }
		  
		  if ( make_hard_call ) 
		    {
		      if      ( pp[0] >= hard_call_prob_threshold ) geno = 1;
		      else if ( pp[1] >= hard_call_prob_threshold ) geno = 2;
		      else if ( pp[2] >= hard_call_prob_threshold ) geno = 3;
		    }
		}
	    }
	  else if ( format_dose2 ) 
	    {

	      double d = atof( toks[p++] );
	      
	      if ( d >= 0 && d <= 2 )
		//	      if ( Helper::str2dbl( toks[p++] , d ) ) 
		{
		  
		  // assume diploid...
		  if ( ! ref1 ) d = 2 - d;
		  
		  if ( store_as_dosage ) 
		    {
		      *dstore = d; 		      
		    }
		  else if ( store_as_prob3 )
		    {
		      
		      std::vector<double> pp(3);

		      // 'impute' probabilities assuming HWE

		      if ( d > 1 ) 
			{
			  pp[0] = 0;
			  pp[1] = 2 - d;
			  pp[2] = d - 1;
			}
		      else
			{
			  pp[0] = 1 - d;
			  pp[1] = d;
			  pp[2] = 0;
			}
		      
		      *dstore     = pp[0];
		      *(dstore+1) = pp[1];
		      *(dstore+2) = pp[2];

		      // g.meta.set( tagname , pp ) ;
		    }
		 
		  if ( make_hard_call ) 
		    {
		      if      ( d <= hard_call_dosage_threshold ) geno = 1;
		      else if ( d >= 2.0 - hard_call_dosage_threshold ) geno = 3;
		      else if ( fabs( d - 1.0 ) <= hard_call_dosage_threshold ) geno = 2;
		    }
		}
	    }
	  else if ( format_dose1 ) 
	    {	      
	      double d;
	      
	      if ( Helper::str2dbl( toks[p++] , d ) ) 
		{		  

		  // internally, use 0..2 encoding of dosage always

		  d *= 2;
		  
		  if ( !ref1 ) d = 1 - d;
		  
		  if ( store_as_dosage ) 
		    {
		      g.meta.set( tagname , d  ) ;
		    }
		  
		  else if ( store_as_prob3 )
		    {
		      
		      std::vector<double> pp(3);
		      // fudge, but probably not too bad??
		      // no other choice...
		      if ( d > 1 ) 
			{
			  pp[0] = 0;
			  pp[1] = 2 - d;
			  pp[2] = d - 1;
			}
		      else
			{
			  pp[0] = 1 - d;
			  pp[1] = d;
			  pp[2] = 0;
			}
		      
		      g.meta.set( tagname , pp ) ;
		    }
		  
		  if ( make_hard_call ) 
		    {
		      if      ( d <= hard_call_dosage_threshold ) geno = 1;
		      else if ( d >= 2.0 - hard_call_dosage_threshold ) geno = 3;
		      else if ( fabs( d - 1.0 ) <= hard_call_dosage_threshold ) geno = 2;		      
		    }
		  
		}
	    }
	
	  
	  //
	  // Make a hard-call?
	  //
	  
	  if ( geno == 0 ) g.null( true );
	  else g.set_alternate_allele_count( geno - 1 );
	  
	  //
	  // Add to variant
	  //
	  
	  //	  v.add( g );
	  v.add( g , i );
	  
	}  // add next individual/genotype
    
      
      //
      // Separate meta-information?
      //

      if ( has_meta ) 
	{
	  
	  std::string id_meta;
	  inp_meta >> id_meta;
	  
	  if ( id_meta != "." && id != id_meta ) 
	    Helper::halt( "mismatching rows in meta-information file: " + id_meta + " vs " + id );

	  for (int m=0;m<metas.size();m++)
	    {
	      
	      double d;
	      std::string x = "";
	      inp_meta >> x;
	      
	      if ( Helper::str2dbl( x , d ) ) 
		{
		  v.consensus.meta.set( metas[m] , d );
		}

	    }
	}
           
      //
      // Put this variant into the database
      //
      
      vardb->insert_consensus( file_id , v );      
      
      
      //
      // Keep track of # of variants added
      //

      ++cnt;

      
    } // next SNP

  
  //
  // Ensure variant index is set on VARDB
  //
  
  int2 niv = vardb->make_summary( file_id ) ;

  //
  // Clean-up and leave
  //

  inp_dosage.close();
  if ( separate_map ) inp_map.close();
  if ( has_meta ) inp_meta.close();

  plog << "inserted " 
       << cnt << " variants, for " 
       << ni << " individuals from "
       << dose_filename << "\n";
  
  return cnt;

}



double DoseReader::Rsq( const Variant & var , double * aaf )
{
  
  
  // Calculate INFO score given alternate allele frequency
  
  const int n = var.size();
  
  bool dosage = false;
  bool prob   = false;

  double frq = 0;
  int cnt = 0;

  for (int i=0; i<n; i++)
    if ( ! var(i).null() )
      {
	
	if ( ! var(i).meta.has_field( Genotype::soft_call_label ) ) continue;
	
	const std::vector<double> & x = var(i).meta.get_double( Genotype::soft_call_label );	
	
	if ( x.size() == 1 ) 
	  {
	    dosage = true;
	    frq += x[0]; // assumes 0..2 encoding
	    ++cnt;
	  }
	else if ( x.size() == 3 ) 
	  {
	    prob = true;
	    frq += x[1] + 2 * x[2];  // assumes biallelic, or REF vs (ALTs)
	    ++cnt;
	  }
      }

  //
  // Sanity check: should be either *all* dosage or *all* probabilities
  //
  
  if ( ( dosage && prob )  || cnt == 0 ) 
    {
      if ( aaf ) *aaf = -1;
      return -1;
    }

  
  //
  // AAF as estimated from the dosage data (2*cnt asusmes diploid for all)
  //

  frq /= (double) 2 * cnt;
  if ( aaf ) *aaf = frq;
  

  //
  // Given HWE, theoretical variance of dosage
  //

  double theoreticalVariance = frq * ( 1 - frq );

  //
  // Get empirical variance of dosage (on 0..1 scaling)
  //

  double dosageSSQ = 0;
  for (int i=0; i<n; i++)
    if ( ! var(i).null() )
      {
	if ( ! var(i).meta.has_field( Genotype::soft_call_label ) ) continue;
	const std::vector<double> & x = var(i).meta.get_double( Genotype::soft_call_label );	
	
	if ( dosage )
	  {
	    double t = x[0]/2.0 - frq;
	    dosageSSQ += t*t;
	  }
	else // implies prob
	  {
	    double t = (x[1]+2*x[2])/2.0-frq ;
	    dosageSSQ += t*t; 
	  }
      }

  double empiricalVariance  = 2 * ( dosageSSQ / (double)cnt);
  
  double rsq = theoreticalVariance > 0 ? empiricalVariance / theoreticalVariance : 0 ;
  
  return rsq;

}

