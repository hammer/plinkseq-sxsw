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
      f.pstart = sql.get_int( stmt_fetch_given_transcript , 4 ) ;
      f.pstop = sql.get_int( stmt_fetch_given_transcript , 5 ) ;
      f.mstr = sql.get_text( stmt_fetch_given_transcript , 6 ) ;
      
      f.chr = sql.get_text( stmt_fetch_given_transcript , 7 ) ;
      f.gstart = sql.get_int( stmt_fetch_given_transcript , 8 ) ;
      f.gstop = sql.get_int( stmt_fetch_given_transcript , 9 ) ;
      
      s.add( transcript , f );
    }
	  
  sql.reset( stmt_fetch_given_transcript );
        
  return s.get( transcript );
}
  

/// Lookup all transcripts and protein domains given a variant
ProtFeatureSet ProtDBase::lookup( const Variant & v )
{
  
  // simple dump of all features for a gene
  
  ProtFeatureSet s;

  std::string chrc = Helper::chrCode( v.chromosome() ) ;
  sql.bind_text( stmt_fetch_given_genomic_coord , ":chr" , chrc ); 
  sql.bind_int( stmt_fetch_given_genomic_coord , ":pos" , v.position() ); 
  
  while ( sql.step( stmt_fetch_given_genomic_coord ) )
    {
      
      Feature f;
      
      f.source_id = sql.get_text( stmt_fetch_given_genomic_coord , 1 );
      
      f.protein_id = sql.get_text( stmt_fetch_given_genomic_coord , 2 ) ;
      f.pstart = sql.get_int( stmt_fetch_given_genomic_coord , 3 ) ;
      f.pstop = sql.get_int( stmt_fetch_given_genomic_coord , 4 ) ;

      f.mstr = sql.get_text( stmt_fetch_given_genomic_coord , 5 ) ;
      
      f.feature_name = sql.get_text( stmt_fetch_given_genomic_coord , 6 ) ;
      f.feature_id = sql.get_text( stmt_fetch_given_genomic_coord , 7 );
             
      f.chr = sql.get_text( stmt_fetch_given_genomic_coord , 8 ) ;
      f.gstart = sql.get_int( stmt_fetch_given_genomic_coord , 9 ) ;
      f.gstop = sql.get_int( stmt_fetch_given_genomic_coord , 10 ) ;

      s.add( sql.get_text( stmt_fetch_given_genomic_coord , 0 ) , f );
    }
	  
  sql.reset( stmt_fetch_given_genomic_coord );
        
  return s;

}



bool ProtDBase::attach( const std::string & name ) 
{
  
  sql.open(name);
  
  sql.synchronous( false );

  sql.query("PRAGMA journal_mode = OFF;");
  sql.query("PRAGMA encoding='UTF-8';");

  // just have a single flat table; will not need to be queried often,
  // so no need to make in normal form etc

  sql.query(" CREATE TABLE IF NOT EXISTS main("
	    "  pk_id          INTEGER PRIMARY KEY , "
	    "  transcript_id  VARCHAR(20) NOT NULL , "
	    "  protein_id     VARCHAR(20) , "
	    "  source_id      VARCHAR(20) NOT NULL , "
	    "  feature_id     VARCHAR(20) NOT NULL , "
	    "  feature_name   VARCHAR(20) , "
	    "  pstart         INTEGER NOT NULL , "
	    "  pstop          INTEGER NOT NULL , "
	    "  meta           VARCHAR(20) , "
	    "  chr            VARCHAR(20) NOT NULL , "
	    "  gstart         INTEGER NOT NULL , "
	    "  gstop          INTEGER NOT NULL , " 
            "  strand         INTEGER NOT NULL ) ;  " );
  
  sql.query(" CREATE TABLE IF NOT EXISTS group_id( group_id VARCHAR(20) NOT NULL ) ; " );
  
  sql.query(" CREATE TABLE IF NOT EXISTS sources("
	    "  source_id VARCHAR(20) , CONSTRAINT c1 UNIQUE( source_id)  ); " );
  
  init();

  return true;
}
  


void ProtDBase::init() 
{

  stmt_dump_types = 
    sql.prepare( "SELECT source_id FROM sources ORDER BY source_id; ");

  stmt_insert = 
    sql.prepare( " INSERT OR IGNORE INTO main ( transcript_id , protein_id , source_id , feature_id , feature_name , pstart , pstop , meta , chr , gstart , gstop , strand  ) "
		 " values( :transcript_id , :protein_id , :source_id , :feature_id , :feature_name , :pstart , :pstop , :meta , :chr , :gstart , :gstop , :strand ) ; " );
  
  stmt_fetch_given_transcript = 
    sql.prepare( " SELECT source_id , feature_id , feature_name , protein_id , pstart , pstop , meta , chr , gstart , gstop FROM main WHERE transcript_id == :transcript_id ; " );

  stmt_fetch_given_annot = 
    sql.prepare( " SELECT transcript_id , protein_id , pstart , pstop , meta , feature_name , feature_id , chr , gstart , gstop FROM main WHERE source_id == :source_id ; " );

  stmt_fetch_given_annot_feature = 
    sql.prepare( " SELECT transcript_id , protein_id , pstart , pstop , meta , feature_name , chr , gstart , gstop FROM main WHERE source_id == :source_id AND feature_id == :feature_id ; " );

  stmt_fetch_given_genomic_coord = 
    sql.prepare( " SELECT transcript_id , source_id , protein_id , pstart , pstop , meta , feature_name , feature_id , chr , gstart , gstop FROM main WHERE chr == :chr AND gstart <= :pos AND gstop >= :pos ; " );

  stmt_fetch_given_genomic_coord_and_source = 
    sql.prepare( " SELECT transcript_id , source_id , protein_id , pstart , pstop , meta , feature_name , feature_id , chr , gstart , gstop FROM main WHERE source_id == :source_id AND chr == :chr AND gstart <= :pos AND gstop >= :pos ; " );


}

void ProtDBase::release()
{
  sql.finalise( stmt_insert );
  sql.finalise( stmt_dump_types );
  sql.finalise( stmt_fetch_given_transcript );
  sql.finalise( stmt_fetch_given_annot );
  sql.finalise( stmt_fetch_given_annot_feature );
  sql.finalise( stmt_fetch_given_genomic_coord );
  sql.finalise( stmt_fetch_given_genomic_coord_and_source );

}

void ProtDBase::index() 
{
  sql.query( " CREATE INDEX IF NOT EXISTS transIdx ON main( transcript_id ) ;" );
  sql.query( " CREATE INDEX IF NOT EXISTS transAnnotIdx ON main( transcript_id , source_id ) ;" );
  sql.query( " CREATE INDEX IF NOT EXISTS gpositionIdx1 ON main( chr , gstart , gstop ) ; " );
  sql.query( " CREATE INDEX IF NOT EXISTS gpositionIdx2 ON main( source_id , chr , gstart , gstop ) ; " );
}

void ProtDBase::drop_index() 
{
  sql.query( "DROP INDEX IF EXISTS transIdx;" );
  sql.query( "DROP INDEX IF EXISTS transAnnotIdx;" );
  sql.query( "DROP INDEX IF EXISTS gpositionIdx1;" );
  sql.query( "DROP INDEX IF EXISTS gpositionIdx2;" );
}


void ProtDBase::load( const std::string & filename , 
		      LocDBase * locdb ,
		      const std::string & lgroup )
		      
{
  
  // ensure output
  bool orig_mode = plog.silent();
  bool orig_mode2 = plog.silent_except_errors();
  plog.silent( false );
  plog.silent_except_errors( false );

  if ( (!locdb) || ! locdb->attached() ) Helper::halt( "no attached LOCDB" );
  if ( ! sql.is_open() ) Helper::halt( "no attached PROTDB" );
  

  // use this LOCDB group to connect transcripts in PROTDB, and map to genomic co-ordinates                                                                                                                     
  int loc_id = locdb->lookup_group_id( lgroup );
  if ( loc_id == 0 ) Helper::halt( "could not find group " + lgroup + " in the LOCDB" );

  // keep track of nominal locus group applicable to this protdb

  sql.query( "DELETE FROM group_id; ");
  sql.query( "INSERT OR REPLACE INTO group_id ( group_id ) values( '" + lgroup + "' );" );

   
  // load transcripts, make name key
  std::map<std::string,Region> plmap;    
  plog << "loading " << lgroup << " transcripts... ";  
  std::set<Region> loci = locdb->get_regions( lgroup );
  std::set<Region>::iterator ii = loci.begin();
  while ( ii != loci.end() )
    {
      plmap[ ii->name ] = *ii;
      ++ii;
    }
  loci.clear();  
  plog << " using " << plmap.size() << "\n";


  // track sources
  std::map<std::string,int> sources;

  // start transaction for entire load0
  sql.begin();
  
  // do not clear out table if already exists

  drop_index();
  
  Helper::checkFileExists( filename );

  InFile F( filename );

  int inserted = 0;
  int failed = 0;

  while ( ! F.eof() ) 
    {
      
      // expecting tab-delimited format: 
      // #transcript_id  protein_id source_id  feature_id  feature_name begin  end  meta 
      
      std::string l = F.readLine();

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
     
      
      if ( ! Helper::str2int( tok( 5 ) , f.pstart ) ) f.pstart = 0;
      if ( ! Helper::str2int( tok( 6 ) , f.pstop ) ) f.pstop = 0;

      // any meta-information to add?  keep as strings
      if ( tok( 7 ) != "" && tok( 7 ) != "." )
	{
	  f.mstr = tok( 7 );	  
	}

      // create a convenient item to hold this row,
      s.add( transcript_id , f );

      // and insert into DB

      if ( insert( s , &plmap ) )
	{	  
	  ++inserted;
	  
	  if ( inserted % 1000 == 0 ) 
	    plog.counter( "parsed " + Helper::int2str( inserted ) + " rows" );
	  
	  sources[ f.source_id ]++;
	  
	  }
      else ++failed;
    }
  
  plog.counter("\n");  
  plog << "inserted " << inserted << " features into PROTDB, failed for " << failed << " rows of input\n";
  
  
  F.close();
  
  plog << "adding index...\n";
  
  index();
  
  // commit
  sql.commit();

  std::map<std::string,int>::iterator jj = sources.begin();
  while ( jj != sources.end() ) 
    {
      sql.query( "INSERT OR IGNORE INTO sources ( source_id ) values( '" + jj->first + "'  ); " );
      plog << jj->second << "\t"
	   << jj->first << " entries\n";
      ++jj;
    }

  plog.silent( orig_mode );
  plog.silent_except_errors( orig_mode2 );
}


std::set<std::string> ProtDBase::get_sources()
{
  std::set<std::string> s;
  if ( ! attached() ) return s;
  while ( sql.step( stmt_dump_types ) ) 
    s.insert( sql.get_text( stmt_dump_types , 0 ) );
  return s;
}



bool ProtDBase::insert( const ProtFeatureSet & fset , std::map<std::string,Region> * plmap )
{ 

  std::map<std::string,std::set<Feature> >::const_iterator ii = fset.feat.begin();
  while ( ii != fset.feat.end() )
    {
      
      bool negative_strand = false;
      bool positive_strand = false;
      int bp_size = 0;

      std::set<Feature>::const_iterator jj = ii->second.begin();
      while ( jj != ii->second.end() )
	{

	  // first, calculate genomic position
	  
	  int start = 1 + ( jj->pstart - 1 ) * 3;
	  int stop  = 1 + ( jj->pstop -  1 ) * 3;

	  
	  std::map<std::string,Region>::iterator kk = plmap->find( ii->first );
	  if ( kk == plmap->end() )
	    {
	      plog.warn( "could not find in LOCDB transcript" , ii->first );
	      return false;
	    }
	  else
	    {
	      
	      const Region & region = kk->second;

	      // strand and total size
	      

	      int ns = region.subregion.size();

              for (int e = 0; e < ns ; e++ )
		{
		  if ( region.subregion[e].CDS() )
		    {
		      bp_size += region.subregion[e].stop.position() - region.subregion[e].start.position() + 1;
		      int s = region.subregion[e].meta.get1_int( PLINKSeq::TRANSCRIPT_STRAND() );
		      if ( s == 0 ) continue;
		      negative_strand = s < 0 ;
		      positive_strand = s > 0 ;			  
		    }		  
		}
	      
	      if ( positive_strand == negative_strand ) 
		{
		  plog.warn( "strand ambiguous for transcript" , ii->first );
		  return false;
		}

	      if ( bp_size % 3 != 0 ) 
		{
		  plog.warn( "CDS size not mod 3 transcript" , ii->first );
		  return false;
		}

	      // determine the genomic position for the start/stop pairs                                                                                                                                            
	      int lwr = 1;
	      int genomic_start = 0;
	      int genomic_stop  = 0;
	      int chr = region.chromosome();
	      
	      if ( negative_strand ) 
		{
                  for (int e = ns-1; e >= 0 ; e-- )
                    {
                      if ( region.subregion[e].CDS() )
                        {
                          int sz = region.subregion[e].stop.position() - region.subregion[e].start.position() + 1;
                          int upr = lwr + sz - 1;
			  
			  //std::cout << ii->first << "\t" << start << " " << stop << " " << lwr << " " << upr << "\n";
			  
                          if ( start >= lwr && start <= upr ) 
			    {			      
			      genomic_start = region.subregion[e].stop.position() - ( start - lwr ) ;
			      //std::cout << "flagging start = " << genomic_start << "\n";
			    }
			  
                          if ( stop  >= lwr && stop  <= upr ) genomic_stop = region.subregion[e].stop.position() - ( stop - lwr ) ;
                          
			  // advance to first position in next                                                                                                   
			  
                          lwr += sz ;
                        }
                      if ( genomic_stop != 0 ) break;
                    }

		}
	      else
		{

		  for (int e = 0; e < ns ; e++ )
		    {
		      // only consider actual CDS regions                                                                                                                                                               
		      if ( region.subregion[e].CDS() )
			{
			  int sz = region.subregion[e].stop.position() - region.subregion[e].start.position() + 1;
			  int upr = lwr + sz - 1;
			  if ( start >= lwr && start <= upr ) genomic_start = ( start - lwr ) + region.subregion[e].start.position();
			  if ( stop  >= lwr && stop  <= upr ) genomic_stop = ( stop - lwr ) + region.subregion[e].start.position();
			  // advance to first position in next                                                                                                                                                          
			  lwr += sz ;
			}
		      if ( genomic_stop != 0 ) break;
		    }

		}

	      if ( negative_strand ) 
		{
		  genomic_stop -= 2;
		}
	      else
		{
		  genomic_stop += 2;
		}
	      

	      if ( genomic_start == 0 || genomic_stop == 0 )
		{
		  plog.warn( "trouble mapping to length of LOCDB transcript " , ii->first ) ;
		  //std::cout << "prob -- " << ii->first << " " << start <<" " << stop << " " << bp_size << " " << genomic_start << " " << genomic_stop << " " << (negative_strand ? "-1" : "+1") << "\n";
		  return false;
		}
	      
	      
	      // enter into protdb
	      std::string chrc = Helper::chrCode( chr );
	      
	      sql.bind_text( stmt_insert , ":transcript_id" , ii->first );
	      sql.bind_text( stmt_insert , ":protein_id"    , jj->protein_id );
	      sql.bind_text( stmt_insert , ":source_id"     , jj->source_id );
	      sql.bind_text( stmt_insert , ":feature_id"    , jj->feature_id );
	      sql.bind_text( stmt_insert , ":feature_name"  , jj->feature_name );
	      sql.bind_int(  stmt_insert , ":pstart" , jj->pstart );
	      sql.bind_int(  stmt_insert , ":pstop"  , jj->pstop );
	      sql.bind_text( stmt_insert , ":meta"   , jj->mstr );
	      sql.bind_text( stmt_insert , ":chr"    , chrc );
	      sql.bind_int(  stmt_insert , ":gstart" , genomic_start );
	      sql.bind_int(  stmt_insert , ":gstop"  , genomic_stop );
	      sql.bind_int(  stmt_insert , ":strand" , negative_strand ? -1 : +1 );
	      sql.step( stmt_insert );
	      sql.reset( stmt_insert );
	    }
	  
	  ++jj;
	}
      ++ii;
    }  

  return true;
}

  
void ProtDBase::dump( Out & pout )
{

  sqlite3_stmt * stmt_dump = 
    sql.prepare( "SELECT * FROM main ORDER BY chr , gstart , gstop ; " );
  
  while ( stmt_dump ) 
    {
      
      Feature f;

      std::string transcript_id = sql.get_text( stmt_dump , 1 );      
      f.protein_id = sql.get_text( stmt_dump , 2 ) ;
      f.source_id = sql.get_text( stmt_dump , 3 ) ;
      f.feature_id = sql.get_text( stmt_dump , 4 ) ;
      f.feature_name = sql.get_text( stmt_dump , 5 ) ;

      f.pstart = sql.get_int( stmt_dump , 6 ) ;
      f.pstop = sql.get_int( stmt_dump , 7 ) ;
      f.mstr = sql.get_text( stmt_dump , 8 ) ;
      
      f.chr = sql.get_text( stmt_dump , 9 ) ;
      f.gstart = sql.get_int( stmt_dump , 10 ) ;
      f.gstop = sql.get_int( stmt_dump , 11 ) ;
      int strand = sql.get_int( stmt_dump , 12 ) ;
      
      pout << transcript_id << ( strand == 1 ? "\t+\t" : "\t-\t" ) << f << "\n";
      
    }

  sql.reset( stmt_dump );
  sql.finalise( stmt_dump );

  return;
}

