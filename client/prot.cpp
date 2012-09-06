#include "func.h"
#include "plinkseq/protdb.h"
#include "plinkseq/output.h"

extern GStore g;

struct aux_protview { 
  bool loc_info;
  Region gene;
  std::set<Feature> features;
};

void f_view_variants_prot_loc_annot( Variant & v , void * p )
{
  std::cout << v << "\t";  
  return;
}

bool Pseq::ProtDB::lookup( const std::string & db , 
			   const std::string & gene , 
			   const std::string & lgroup , 
			   Mask * mask )
{

  Out & pout = Out::stream( "prot" );

  bool query_locdb = g.locdb.attached();
  
  std::set<std::string> transcripts;
  
  std::string group = lgroup;
  
  if ( query_locdb )
    {
      
      // swap in, e.g., 'refseq' if not specified
      if ( group == "" ) group = PLINKSeq::DEFAULT_LOC_GROUP();
      
      // See if the 'gene' maps to multiple transcripts
      transcripts = g.locdb.targetted_lookup_alias( gene ,
						    PLINKSeq::DEFAULT_GENE_SYMBOL() , 
						    group );      
    }
  
  if ( transcripts.size() == 0 ) transcripts.insert( gene );


  //
  // Attach a PROTDB
  //

  ProtDBase protdb;
  
  protdb.attach( db );
  
  if ( ! protdb.attached() ) Helper::halt( "could not attach PROTDB" );  


  //
  // Repeat query to PROTDB for each unique transcript
  //

  std::set<std::string>::iterator tt = transcripts.begin();
  while ( tt != transcripts.end() )
    {
      
      //
      // collect info
      //
      
      aux_protview aux;
      
      
      //
      // Get protein features/domains
      //
      
      aux.features = protdb.fetch( *tt );
      
      
      //
      // Query VARDB
      //
      
      if ( mask )
	{

	  Mask m2 = *mask;

      	  m2.include_loc( group );
	  m2.subset_loc( group , *tt );
	  
	  g.vardb.iterate( f_view_variants_prot_loc_annot , &aux , m2 );

	}
          
      
      // 
      // Write integrated map, with variants, exomes and variants
      //
      
      pout << "\n*** " << *tt << "\n\n";
      
      std::set<Feature>::iterator ii = aux.features.begin();
      while ( ii != aux.features.end() )
	{
	  pout << *ii << "\n";
	  ++ii;
	}


      ++tt;

    } // next transcript
  
  return true;
}

