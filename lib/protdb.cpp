#include "plinkseq/protdb.h"
#include "plinkseq/locdb.h"
#include "plinkseq/output.h"

#include <iostream>

std::set<Feature> ProtDBase::fetch( const std::string & transcript )
{

  // simple dump of all features for a gene

  ProtFeatureSet s;

  sql.bind_text( stmt_fetch_given_transcript , ":transcript_id" , transcript );
  while ( sql.step( stmt_fetch_given_transcript ) )
    {
      Feature f;
      f.source_id = sql.get_text( stmt_fetch_given_transcript , 0 ) ;
      f.feature_id = sql.get_text( stmt_fetch_given_transcript , 1 ) ;
      f.feature_name = sql.get_text( stmt_fetch_given_transcript , 2 ) ;
      f.protein_id = sql.get_text( stmt_fetch_given_transcript , 3 ) ;
      f.start = sql.get_int( stmt_fetch_given_transcript , 4 ) ;
      f.stop = sql.get_int( stmt_fetch_given_transcript , 5 ) ;
      f.mstr = sql.get_text( stmt_fetch_given_transcript , 6 ) ;
      s.add( transcript , f );
    }
	  
  sql.reset( stmt_fetch_given_transcript );
        
  return s.get( transcript );
}
  

/// Map a PROTDB source/feature to a LOCDB (i.e. to make a genomic region list)

int ProtDBase::map_to_genomic( LocDBase * locdb , 
			       const std::string & lgroup , 
			       const std::string & source , 
			       const std::string & feature , 
			       std::map<std::string,Region> * plmap )
{
  
  if ( (!locdb) || ! locdb->attached() ) Helper::halt( "no attached LOCDB" );
  if ( ! sql.is_open() ) Helper::halt( "no attached PROTDB" );

  // use this LOCDB group to connect transcripts in PROTDB, and map to genomic co-ordinates
  int loc_id = locdb->lookup_group_id( lgroup );
  if ( loc_id == 0 ) Helper::halt( "could not find group " + lgroup + " in the LOCDB" );
  

  bool select_feature = feature != ".";

  // use the 'source' as the new loc-group name 
  int new_id = locdb->set_group_id( ( select_feature ? source + "::" + feature : source ) ) ;


  //
  // Output stream
  //

  Out & pout = Out::stream( "protmap" );  

  //
  // Use existing 'Annotate' class to store table of transcripts for quick matching
  //

  bool need_to_delete = false;

  if ( plmap == NULL ) 
    {
      need_to_delete = true;
      plog << "loading " << lgroup << " transcripts... ";
      plmap = new std::map<std::string,Region>;
      std::set<Region> loci = locdb->get_regions( lgroup );
      std::set<Region>::iterator ii = loci.begin();
      while ( ii != loci.end() )
	{
	  (*plmap)[ ii->name ] = *ii;
	  ++ii;
	}
      loci.clear();
    }
  
  plog << " using " << plmap->size() << "\n";



  
  //
  // Iterate through PROTDB
  //

  int mapped = 0;

  sqlite3_stmt * s = select_feature ? stmt_fetch_given_annot_feature : stmt_fetch_given_annot ;


  sql.bind_text( s , ":source_id" , source );
  
  if ( select_feature )
    sql.bind_text( s , ":feature_id" , feature );
  
  while ( sql.step( s ) )
    {
      
      std::string transcript = sql.get_text( s , 0 );
      std::string protein_id = sql.get_text( s , 1 );
      int start = sql.get_int( s , 2 );
      int stop  = sql.get_int( s , 3 );
      std::string mstr = sql.get_text( s , 4 );
      
      std::string feature_id = ".";
      std::string feature_name = sql.get_text( s , 5 );
      
      if ( ! select_feature )
	{
	  feature_id = sql.get_text( s , 6 );	  
	}

      
      //
      // Look up genomic co-ordinates, given transcript and AA positions 
      //
      
      // convert AA (1-based) to BP (1-based) 
      start = 1 + ( start - 1 ) * 3;
      stop  = 1 + ( stop -  1 ) * 3;
      
      std::map<std::string,Region>::iterator ii = plmap->find( transcript );
      if ( ii == plmap->end() )
	{
	  pout << "#no region for " << transcript << "\n";
	}
      else
	{
	  const Region & region = ii->second;
	  // determine the genomic position for the start/stop pairs
	  
	  int lwr = 1;
	  int ns = region.subregion.size();
	  int genomic_start = 0;
	  int genomic_stop  = 0;
	  int chr = region.chromosome();

	  for (int e = 0; e < ns ; e++ )
	    {
	      // only consider actual CDS regions
	      if ( region.subregion[e].CDS() )
		{
		  int sz = region.subregion[e].stop.position() - region.subregion[e].start.position() + 1;
		  int upr = lwr + sz;
		  if ( start >= lwr && start <= upr ) genomic_start = ( start - lwr ) + region.subregion[e].start.position();
		  if ( stop  >= lwr && stop  <= upr ) genomic_stop = ( stop - lwr ) + region.subregion[e].start.position();
		  // advance to first position in next 
		  lwr += sz + 1 ;
		}
	      if ( genomic_stop != 0 ) break;
	    }
	  
	  if ( genomic_start == 0 || genomic_stop == 0 ) 
	    pout << "#trouble mapping to transcript " << transcript << " for " << start << ".." << stop << "\n";

	  ++mapped;

	  pout << transcript << "\t"
	       << protein_id << "\t"
	       << start << ".."
	       << stop << "\t"
	       << Helper::chrCode( chr ) << ":"
	       << genomic_start << ".."
	       << genomic_stop << "\t"
	       << source << "\t"
	       << feature << "\t"
	       << feature_name << "\t"

	       << mstr << "\n";
	  
	}

      // next transcript
    }

  sql.reset( s );
  
  if ( need_to_delete ) delete plmap;

  return mapped;
  
}




bool ProtDBase::attach( const std::string & name ) 
{
  
  sql.open(name);
  
  sql.synchronous( false );

  // just have a single flat table; will not need to be queried often,
  // so no need to make in normal form etc

  sql.query(" CREATE TABLE IF NOT EXISTS main("
	    "  pk_id          INTEGER PRIMARY KEY , "
	    "  transcript_id  VARCHAR(2) NOT NULL , "
	    "  protein_id     VARCHAR(20) , "
	    "  source_id      VARCHAR(20) NOT NULL , "
	    "  feature_id     VARCHAR(20) NOT NULL , "
	    "  feature_name   VARCHAR(20) , "
	    "  start          INTEGER NOT NULL , "
	    "  stop           INTEGER NOT NULL , "
	    "  meta           VARCHAR(20) ); " );

  sql.query(" CREATE TABLE IF NOT EXISTS mapping("
	    "  pk_id          INTEGER NOT NULL , "
	    "  chr            VARCHAR(20) NOT NULL , "
	    "  start          INTEGER NOT NULL , "
	    "  stop           INTEGER NOT NULL ) ;  " );


  init();

  return true;
}
  


void ProtDBase::init() 
{
  
  stmt_insert = 
    sql.prepare( " INSERT OR IGNORE INTO main ( transcript_id , protein_id , source_id , feature_id , feature_name , start , stop , meta ) "
		 " values( :transcript_id , :protein_id , :source_id , :feature_id , :feature_name , :start , :stop , :meta ) ; " );
  
  stmt_fetch_given_transcript = 
    sql.prepare( " SELECT source_id , feature_id , feature_name , protein_id , start , stop , meta FROM main WHERE transcript_id == :transcript_id ; " );

  stmt_fetch_given_annot = 
    sql.prepare( " SELECT transcript_id , protein_id , start , stop , meta , feature_name , feature_id FROM main WHERE source_id == :source_id ; " );

  stmt_fetch_given_annot_feature = 
    sql.prepare( " SELECT transcript_id , protein_id , start , stop , meta , feature_name FROM main WHERE source_id == :source_id AND feature_id == :feature_id ; " );

  stmt_insert_mapping = 
    sql.prepare( " INSERT OR IGNORE INTO mapping ( pk_id , chr , start , stop ) values( :pk_id , :chr , :start , :stop ) ; " );

  stmt_fetch_mapping = 
    sql.prepare( " SELECT chr , start , stop FROM mapping WHERE pk_id == :pk_id ; " );

}

void ProtDBase::release()
{
  sql.finalise( stmt_insert );
  sql.finalise( stmt_fetch_given_transcript );
  sql.finalise( stmt_fetch_given_annot );
  sql.finalise( stmt_fetch_given_annot_feature );
}

void ProtDBase::index() 
{
  sql.query( " CREATE INDEX IF NOT EXISTS transIdx ON main( transcript_id ) ;" );
  sql.query( " CREATE INDEX IF NOT EXISTS transAnnotIdx ON main( transcript_id , source_id ) ;" );
}

void ProtDBase::drop_index() 
{
  sql.query( "DROP INDEX IF EXISTS transIdx;" );
  sql.query( "DROP INDEX IF EXISTS transAnnotIdx;" );
}


void ProtDBase::load( const std::string & filename )
{
  
  // do not clear out table if already exists

  drop_index();
  
  Helper::checkFileExists( filename );

  InFile F( filename );

  int inserted = 0;

  while ( ! F.eof() ) 
    {
      
      // expecting tab-delimited format: 
      // #transcript_id  protein_id source_id  feature_id  feature_name begin  end  meta 
      
      std::string l = F.readLine();
      std::cout << "proc [ " << l << "]\n";

      if ( l.size() == 0 ) continue;
      if ( l[0] == '#' ) continue; // skip header/comments
      
      int n = 0 ;
      Helper::char_tok tok( l , &n , '\t' );

      if ( tok.size() != 8 ) 
	Helper::halt( "expecting 8 tab-separated columns in input:\n" + l + "\n" );

      ProtFeatureSet s;
      Feature f;
      
      std::string transcript_id = tok(0);
      
      f.protein_id      = tok( 1 );
      f.source_id       = tok( 2 );
      f.feature_id      = tok( 3 );
      f.feature_name    = tok( 4 );

      if ( ! Helper::str2int( tok( 5 ) , f.start ) ) f.start = 0;
      if ( ! Helper::str2int( tok( 6 ) , f.stop ) ) f.stop = 0;

      // any meta-information to add?  keep as strings
      if ( tok( 7 ) != "" && tok( 7 ) != "." )
	{
	  f.mstr = tok( 7 );
	  //f.meta.set( tok( 7 ) );
	}

      // create a convenient item to hold this row,
      s.add( transcript_id , f );

      // and insert into DB
      insert( s );
      ++inserted;
      
    }

  plog << "inserted " << inserted << " features into PROTDB\n";

  F.close();
  
  plog << "adding index now...\n";
  index();
}



void ProtDBase::insert( const ProtFeatureSet & fset )
{ 
  std::map<std::string,std::set<Feature> >::const_iterator ii = fset.feat.begin();
  while ( ii != fset.feat.end() )
    {
      std::set<Feature>::const_iterator jj = ii->second.begin();
      while ( jj != ii->second.end() )
	{
	  
	  sql.bind_text( stmt_insert , ":transcript_id" , ii->first );
	  sql.bind_text( stmt_insert , ":protein_id" , jj->protein_id );
	  sql.bind_text( stmt_insert , ":source_id" , jj->source_id );
	  sql.bind_text( stmt_insert , ":feature_id" , jj->feature_id );
	  sql.bind_text( stmt_insert , ":feature_name" , jj->feature_name );
	  sql.bind_int( stmt_insert , ":start" , jj->start );
	  sql.bind_int( stmt_insert , ":stop" , jj->stop );
	  sql.bind_text( stmt_insert , ":meta" , jj->mstr );
	  sql.step( stmt_insert );
	  sql.reset( stmt_insert );
	  ++jj;
	}
      ++ii;
    }
}
  

