#include "pp.h"

#include "pp.pb.h"
#include "filemap.h"
#include "regions.h"
#include "locdb.h"
#include "annot.h"
#include "gstore.h"
#include "sqlz.h"

#include <iostream>

extern GStore * GP;

const std::string PPH2DBase::transcript_set_name = "refseq";

bool PPH2DBase::score( const Variant & v , double & score , int & prediction )
{
  
  score = 0;  
  prediction = 0;
  
  bool observed = false;
  bool nonsynon = false;
 
  // obtain all overlapping RefSeq transcripts


  std::set<Region> regions = locdb->get_regions( transcript_set_name , v );

  std::set<Region>::const_iterator i = regions.begin();
  
  while ( i != regions.end() )
    {

      //
      // Look-up gene in PPH2 database
      //

      PPH2Set * s = lookup( i->name );

      if ( ! s ) { ++i; continue; } 
      
      observed = true;

      //
      // Given genomic position of this variant, get annotation using
      // PSEQ databases 
      //
      
      std::set<SeqInfo> seqann = Annotate::annotate( v.chromosome() , 
						     v.position() , 
						     v.alternate(), 
						     v.reference() ,
						     *i );

      
      std::set<SeqInfo>::iterator j = seqann.begin();
      while ( j != seqann.end() )
	{
	  
	  //
	  // is this a synomous change? 
	  //
	  
	  if ( j->ref_aa == j->alt_aa ) 
	    {
	      ++j;
	      continue;
	    }
	  
	  nonsynon = true;
	  
	  //
	  // See if we have a match for this position and reference/alternate AA PPH2
	  //
	  
	  const PPH2Position * p = s->position( j->ppos1 , j->ref_aa , j->alt_aa );
	  
	  if ( ! p ) 
	    {	      
	      ++j;
	      continue;
	    }
	  
	  
	  //
	  // Some temporary output
	  //
	  
	  //        plog << " transcript  = " << j->transcript << "\n"
	  //  		<< " aa position = " << j->ppos1 << " " << j->ppos2 << "\n"
	  //  		<< " aa change   = " << j->ref_aa << ">" << j->alt_aa << "\n"
	  //  		<< " status      = " << j->status() << "\n";
	  
	  std::string ref_aa = Annotate::translate_reference( *i );
	  
	  //        plog << "PSEQ entry contains " << ref_aa.size() << " positions\n";
	  //        plog << "PPH2 entry contains " << s->max_position() << " AAs\n";
	  
	  //
	  // does the implied length in PPH2 match PSEQ?
	  //
	  
	  // if ( ref_aa.size() != s->max_position() )
	  //   plog << "flagging " << s->protein_name << " " << s->transcript_name << " with discrepant lengths\n";
	  
	  
	  //        plog << "_found PPH2: " << p->reference << ">" << p->alternate 
	  // 		 << " = " << p->prediction << " " << p->score << "\n";
	  
	  
	  if ( p->prediction > 0 ) 
	    {
	      
	      observed = true;
	      
	      if ( p->prediction > prediction ) 
		prediction = p->prediction;
	      
	      if ( p->score > score ) 
		score = p->score;
	      
	    }
	 
	  ++j;
	}
      
      ++i;
    }
  
  if ( ! observed ) prediction = nonsynon ? -1 : 0 ;
  
  return observed;
}


bool PPH2DBase::attach( const std::string & name ) 
{

  sql.open(name);

//   sqlite3_create_function( sql.pointer(), "mycompress", 1, SQLITE_UTF8, 0, &compressFunc, 0, 0);
//   sqlite3_create_function( sql.pointer(), "myuncompress", 1, SQLITE_UTF8, 0, &uncompressFunc, 0, 0);
  
  sql.synchronous(false);
  
  sql.query(" CREATE TABLE IF NOT EXISTS main("
	    "  prot_id INTEGER PRIMARY KEY , "
	    "  data   BLOB ); " );

  sql.query(" CREATE TABLE IF NOT EXISTS genes("
	    "  gene_name  VARCHAR(20) NOT NULL, "
            "  prot_name  VARCHAR(20) NOT NULL, "
	    "  prot_id    INTEGER PRIMARY KEY ); " );
  
  init();

  return true;
}
  


void PPH2DBase::init() 
{

  stmt_insert_pph2_scores = 
    sql.prepare( " INSERT OR IGNORE INTO main ( prot_id, data ) values( :prot_id, :data ); " );

  stmt_fetch_pph2_scores = 
    sql.prepare( " SELECT data FROM main WHERE prot_id == :prot_id ; " );

  stmt_insert_pph2_id = 
    sql.prepare( " INSERT OR IGNORE INTO genes ( gene_name, prot_name ) values( :gene_name, :prot_name ); " );

  stmt_fetch_pph2_id = 
    sql.prepare( " SELECT prot_id, prot_name FROM genes WHERE gene_name == :gene_name ; " );

}

void PPH2DBase::release()
{
  sql.finalise( stmt_insert_pph2_scores );
  sql.finalise( stmt_fetch_pph2_scores );
  sql.finalise( stmt_insert_pph2_id );
  sql.finalise( stmt_fetch_pph2_id );
}

void PPH2DBase::index() 
{
  sql.query( "CREATE INDEX IF NOT EXISTS i1 ON genes(gene_name);" );
}

void PPH2DBase::drop_index() 
{
  sql.query( "DROP INDEX IF EXISTS i1;" );
}

void PPH2DBase::load( const std::string & filename )
{

  drop_index();
  
  InFile F( filename );
  
  PPH2Set pset;

  // Assume file is sorted by transcript

  std::string curr = "";

  while ( ! F.eof() ) 
    {


      // format: protein transcript pos AA1 AA2 0/1/2/3 prob-score
      // where 0/1/2/3 = ?/benign/possibly/probably damaging
      
      std::vector<std::string> tok = F.tokenizeLine( );
      
      if ( tok.size() == 0 ) continue;

      if ( tok.size() != 7 ) 
	{
	  plog.warn("found input row with wrong # of columns");
	  plog << tok.size() << " : ";
	  for (int i=0; i<tok.size(); i++)
	    plog << tok[i] << " ";
	  plog << "\n";
	  continue;
	}

      if ( pset.transcript_name != tok[1] )
	{

	  if ( pset.transcript_name != "" ) 
	    {
	      insert(pset);
	      pset.reset();	      
	    }
	  pset.protein_name = tok[0];
	  pset.transcript_name = tok[1];
	}
 
      accumulate( pset , tok );
      
    }
  F.close();
  
  index();
}



void PPH2DBase::accumulate( PPH2Set & pset , 
			    const std::vector<std::string> & tok )
{

  PPH2Position pos;
  int p = 0;
  if ( ! Helper::str2int( tok[2] , p ) )
    Helper::halt("bad format for position");
  pos.reference = tok[3];
  pos.alternate = tok[4];

  if ( ! Helper::str2dbl( tok[6] , pos.score ) )
    pos.score = 0;
  
  if ( ! Helper::str2int( tok[5] , pos.prediction ) )
    pos.prediction = 0;
  
  // add to gene in the appropriate slot for amino-acod 'p'
  pset.scores[p].insert(make_pair( pos.reference+pos.alternate, pos) );
  
}

  

void PPH2DBase::insert( const PPH2Set & pset )
{
  
// encode a set of PPH2 scores for a gene, and stick in database as BLOB
  
  PolyPhen2Buffer b;

  b.set_protein_name( pset.protein_name );
  b.set_transcript_name( pset.transcript_name );

  std::map<int,std::map<std::string,PPH2Position> >::const_iterator i = pset.scores.begin();
  
  while ( i != pset.scores.end() )
    {
      std::map<std::string,PPH2Position>::const_iterator j = i->second.begin();
      while ( j != i->second.end() )
	{

	  b.add_position( i->first );
	  b.add_reference( j->second.reference );
	  b.add_alternate( j->second.alternate );
	  b.add_score( j->second.score ); 
	    
	  if ( j->second.prediction == 0 ) 
	    b.add_prediction( PolyPhen2Buffer::UNKNOWN );
	  else if ( j->second.prediction == 1 ) 
	    b.add_prediction( PolyPhen2Buffer::BENIGN );
	  else if ( j->second.prediction == 2 ) 
	    b.add_prediction( PolyPhen2Buffer::POSS );
	  else if ( j->second.prediction == 3 ) 
	    b.add_prediction( PolyPhen2Buffer::PROB );
	  
	  ++j;
	}
      ++i;
    }
  
  // encode as BLOB

  std::string s;  
  b.SerializeToString(&s);
  blob encoded(s);

  // store gene/protein names
  sql.bind_text( stmt_insert_pph2_id , ":gene_name" , pset.transcript_name );
  sql.bind_text( stmt_insert_pph2_id , ":prot_name" , pset.protein_name );
  sql.step( stmt_insert_pph2_id );
  sql.reset( stmt_insert_pph2_id );

  // retrieve prot_id
  sql.bind_text( stmt_fetch_pph2_id , ":gene_name" , pset.transcript_name );
  sql.step( stmt_fetch_pph2_id );
  uint64_t prot_id = sql.get_int64( stmt_fetch_pph2_id , 0 );
  sql.reset( stmt_fetch_pph2_id );

  // insert actual data
  sql.bind_int64( stmt_insert_pph2_scores , ":prot_id" , prot_id );
  sql.bind_blob( stmt_insert_pph2_scores , ":data" , encoded );
  sql.step( stmt_insert_pph2_scores );
  sql.reset( stmt_insert_pph2_scores );

}
  
bool PPH2DBase::present( const std::string & gene )
{
  PPH2Set * s = lookup( gene );
  return s != NULL;
}

PPH2Set * PPH2DBase::lookup( const std::string & gene )
{
  
  // Lookup (from cache first) based on transcript name (RefSeq)
  
  std::map<std::string, PPH2Set>::iterator i = cache.find( gene );
  
  if ( i != cache.end() ) 
    return &(i->second);
  
  // do we need to clear cache? 
  
  if ( cache.size() > 10 ) 
    cache.clear();

  sql.bind_text( stmt_fetch_pph2_id , ":gene_name" , gene );
  sql.step( stmt_fetch_pph2_id );
  uint64_t prot_id = sql.get_int64( stmt_fetch_pph2_id , 0 );
  sql.reset( stmt_fetch_pph2_id );
  if ( prot_id == 0 ) return NULL;

  sql.bind_int64( stmt_fetch_pph2_scores , ":prot_id" , prot_id );

  if ( sql.step( stmt_fetch_pph2_scores ) )
    {

      blob b = sql.get_blob( stmt_fetch_pph2_scores , 0 );
      
      PolyPhen2Buffer pb;
      
      pb.ParseFromString( b.get_string() );
      
      //plog << pb.DebugString() << "\n\n";

      PPH2Set ps;
      
      ps.protein_name = pb.protein_name();

      ps.transcript_name = pb.transcript_name();

      int n = pb.score().size();
      
      for (int i=0; i<n; i++)
	{
	  PPH2Position pp;
	  
	  int position = pb.position(i);
	  
	  pp.reference = pb.reference(i);
	  pp.alternate = pb.alternate(i);
	  pp.score = pb.score(i);
	  pp.prediction = pb.prediction(i);

	  ps.scores[ position ].insert( make_pair( pp.reference + pp.alternate , pp ) );
	  
	}
       
      // and store
      cache[ gene ] = ps;
    }

 
  sql.reset( stmt_fetch_pph2_scores );

  return cache.find( gene ) == cache.end() ? NULL : &( cache[gene] );

}

