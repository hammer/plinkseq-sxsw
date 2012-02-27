
#include "seqdb.h"
#include "filemap.h"
#include "gstore.h"
#include "sqlz.h"

#include <algorithm>
#include <cmath>

extern GStore * GP;

using namespace std;
using namespace Helper;


// IUPAC codes 
// A    A
// C    C
// G    G
// T    T

// M    A or C
// R    A or G
// W    A or T
// S    C or G
// Y    C or T
// K    G or T

// V    A or C or G
// H    A or C or T
// D    A or G or T
// B    C or G or T

// N    G or A or T or C


bool SeqDBase::iupac( const std::string & c , const std::string & b )
{

  // is a given base 'b' consistent with the IUPAC code 'c'

  if ( b == "A" ) return c == "A" || c == "N" || c == "M" || c == "R" || c == "W" || c == "V" || c == "H" || c == "D";
  if ( b == "C" ) return c == "C" || c == "N" || c == "M" || c == "S" || c == "Y" || c == "V" || c == "H" || c == "B";
  if ( b == "G" ) return c == "G" || c == "N" || c == "R" || c == "S" || c == "K" || c == "V" || c == "D" || c == "B";
  if ( b == "T" ) return c == "T" || c == "N" || c == "W" || c == "Y" || c == "K" || c == "H" || c == "D" || c == "B";
  return false;

}
    
bool SeqDBase::create( const std::string & name )
{
  if ( Helper::fileExists(name) ) return false;
  sql.open( name );
  sql.close();
  attach( name );
  return true;
}

bool SeqDBase::attach( const std::string & name )
{
  
  if ( name == "-" || name == "." ) { dettach(); return false; } 

  // Only open existing databases  
  if ( ! Helper::fileExists( name ) ) { dettach(); return false; } 
  
  sql.open( name );
  
  // register compression functions 
  
  sqlite3_create_function( sql.pointer(), "mycompress", 1, SQLITE_UTF8, 0, &compressFunc, 0, 0);
  sqlite3_create_function( sql.pointer(), "myuncompress", 1, SQLITE_UTF8, 0, &uncompressFunc, 0, 0);
  
  sql.synchronous(false);
  
  sql.query(" CREATE TABLE IF NOT EXISTS refseq("
            "   chr       INTEGER NOT NULL , "
            "   bp1       INTEGER NOT NULL , "
            "   bp2       INTEGER NOT NULL , "
	    "   seq       BLOB ); ");

  sql.query( " CREATE TABLE IF NOT EXISTS meta("
	     "   key    VARCHAR(20) NOT NULL , "
	     "   value  VARCHAR(20) NOT NULL ) ; " );
	       
  init();
  setMinMax();
  
  // set some core values
  lookup_meta();
  
}


void SeqDBase::index()
{
    sql.query( "DROP INDEX IF EXISTS ind1;");
    sql.query( "CREATE INDEX ind1 ON refseq(chr,bp1,bp2);");
    release();
    init();
}

bool SeqDBase::init()
{    
    stmt_insert =
	sql.prepare( "INSERT OR REPLACE INTO refseq ( chr, bp1, bp2, seq ) "
		     " values( :chr, :bp1, :bp2, mycompress( :seq ) ); " );
    stmt_lookup = 
	sql.prepare( " SELECT bp1,bp2,myuncompress( seq ) FROM refseq "
		     " WHERE chr == :chr AND :bp1 <= bp2 AND :bp2 >= bp1 ORDER BY chr,bp1 ; ");

    stmt_getmeta = 
      sql.prepare( "SELECT key,value FROM meta ; " );

    stmt_putmeta = 
      sql.prepare( "INSERT OR REPLACE INTO meta ( key , value ) values (:key,:value) ; " );
    

}


void SeqDBase::insert_meta( const std::map< std::string , std::string > & m ) 
{

  // store internally in class
  meta = m;

  // and in DB
  std::map< std::string , std::string >::const_iterator i = m.begin();
  sql.begin();
  while ( i != m.end() )
    {
      sql.bind_text( stmt_putmeta , ":key" , i->first );
      sql.bind_text( stmt_putmeta , ":value" , i->second );
      sql.step( stmt_putmeta );
      sql.reset(stmt_putmeta);
    ++i;
    }  
  sql.commit();
}

std::map< std::string , std::string > SeqDBase::lookup_meta()
{

  while ( sql.step( stmt_getmeta ) )
    {
      std::string k = sql.get_text( stmt_getmeta , 0 ) ;
      std::string v = sql.get_text( stmt_getmeta , 1 ) ;

      meta[k] = v;
      
      if ( k == PLINKSeq::SEQDB_GENOME_BUILD_KEY() ) 
	genome_build = v;
      else if ( k == PLINKSeq::SEQDB_REPEAT_MODE_KEY() )
	{
	  if ( v == "none" ) rpt_mode = SEQ_RPT_NONE ;
	  else if ( v == "N" ) rpt_mode = SEQ_RPT_N;
	  else if ( v == "lower" ) rpt_mode = SEQ_RPT_LOWER;
	  else rpt_mode = SEQ_RPT_UNKNOWN;
	}
      else if ( k == PLINKSeq::SEQDB_NAME_KEY() ) name = v;
      else if ( k == PLINKSeq::SEQDB_DESCRIPTION_KEY() ) desc = v;
      else if ( k == PLINKSeq::SEQDB_IUPAC_KEY() ) use_iupac = v == "1";
    }
  sql.reset(stmt_getmeta);  

  return meta;
}

bool SeqDBase::release()
{    
  sql.finalise( stmt_lookup );
  sql.finalise( stmt_insert );
}

bool SeqDBase::dettach()
{
  if ( attached() ) 
    {
      release();
      sql.close();  
    }
}


void  SeqDBase::loadFASTA( const std::string & filename ,
			   const std::map< std::string , std::string > & mymeta )
{

  if ( ! attached() ) return;

  
  InFile f( filename );
  
  int c = 0;
  int p1 = 0;
  std::string s = "";
  
  sql.begin();
  
  while ( ! f.eof() )
    {
      std::string l = f.readLine();

      if ( l == "" ) continue;
      
      if ( l.substr(0,1) == ">" )
	{
	  // Do we need to flush the buffer?
	  if ( s.size()>0 )
	    {
	      insert( c, p1, p1 + s.size() -1 , s );			
	    }
	  
	  if ( l.find("_") == std::string::npos ) 
	    {
	      // change "1" to "chr1", etc
	      std::string chr_str = Helper::defaultChrPrefix(l.substr( 1 ));
	      c = chrCode( chr_str );
	      
	      if ( c > 0 ) 
		{
		  plog << "Reading chromosome " << Helper::chrCode( c ) << "\n";
		  p1 = 1;
		  s = "";
		}
	    }
	  else
	    {
	      c = 0;
	      plog << "Skipping " << l << "\n";
	    }
	}
      else
	{
	  if ( c > 0 ) 
	    {
	      s += l;
	      if ( s.size() > BATCH_SIZE ) 
		{
		  insert( c, p1, p1 + s.size() -1 , s );					  
		  p1 += s.size();
		  s="";				
		}
	    }
	}
      
    }
  

  // 
  // Flush last batch
  //

  if ( s.size()>0 )
    {
      insert( c, p1, p1 + s.size() -1 , s );			
      s="";
    }
  
 
  sql.commit();
 
   
  //
  // Index, and set min/max for each chromosome
  //
  
  index();
  
  setMinMax();


  //
  // Parse meta-information
  //

  
  insert_meta( mymeta );
  
  
  //
  // Clean-up
  //
  
  f.close();
    
}


void  SeqDBase::insert(const int chr, const int bp1, const int bp2, const std::string & s)
{
    sql.bind_int( stmt_insert , ":chr" , chr );
    sql.bind_int( stmt_insert , ":bp1" , bp1  );
    sql.bind_int( stmt_insert , ":bp2" , bp2  );
    sql.bind_text( stmt_insert , ":seq" , s  );
    sql.step( stmt_insert );
    sql.reset( stmt_insert );
}
        

std::string SeqDBase::lookup( const Variant & v, int window )
{
  return lookup( v.chromosome(), v.position() - window , v.position() + window ); 
}

std::string SeqDBase::lookup( const Region & region ) 
{
  return lookup( region.chromosome() , region.start.position() , region.stop.position() ) ;
}



std::string SeqDBase::lookup( int chr , int bp1 , int bp2 )
{
    
    if ( bp2 == 0 ) bp2 = bp1;
    else if ( bp2 < bp1 ) return "";

    // Fully contained within cache?
    
    if ( chr == cache_chr && 
	 bp1 >= cache_bp1 && 
	 bp2 <= cache_bp2 ) 
    {
	return chunk.substr( bp1 - cache_bp1 , bp2-bp1+1 );
    }
    
    // Is this query within bounds
    
    std::map<int,int2>::iterator i = chrminmax.find( chr );
    if ( i == chrminmax.end() ) return "";
    
    if ( bp1 < i->second.p1 ) bp1 = i->second.p1;
    if ( bp2 < i->second.p1 ) bp2 = i->second.p1;

    if ( bp1 > i->second.p2 ) bp1 = i->second.p2;
    if ( bp2 > i->second.p2 ) bp2 = i->second.p2;

    
    // Otherwise, clear cache

    cache_chr = chr;
    
    sql.bind_int( stmt_lookup , ":chr" , chr );
    sql.bind_int( stmt_lookup , ":bp1" , bp1 );
    sql.bind_int( stmt_lookup , ":bp2" , bp2 );
  
    bool first = true; 
    

    while ( sql.step( stmt_lookup ) )
    {
       int p1 = sql.get_int( stmt_lookup , 0 );
       int p2 = sql.get_int( stmt_lookup , 1 );

       if ( first ) 
       { 
	   cache_bp1 = p1;       
	   cache_bp2 = p2;
	   chunk = sql.get_text( stmt_lookup , 2 ); 
	   first = false;
       }
       else
       {  
	   // A standard, sequential chunk? 
	   if ( p1 == cache_bp2 + 1 ) 
	   {
	       std::string s = sql.get_text( stmt_lookup , 2 );	       
	       chunk += s;
	       cache_bp2 = p2;
	   }
	   else Helper::halt("Not implemented yet -- only sequential ref chunks");
	   
       }
    }
    
    sql.reset( stmt_lookup );	

    return chunk.substr( bp1 - cache_bp1 , bp2-bp1+1 );
}

bool SeqDBase::N( const Region & region , int & n , int & tot ) 
{
  if ( rpt_mode == SEQ_RPT_NONE || rpt_mode == SEQ_RPT_UNKNOWN ) return false;

  std::string s = lookup( region );
  n = 0;
  tot = s.size();
  if ( tot == 0 ) return false;
  if ( rpt_mode == SEQ_RPT_N )
    {
      for (int i=0; i<tot; i++)
	if ( s[i] == 'N' ) ++n;
    }
  else if ( rpt_mode == SEQ_RPT_LOWER )
    {
      for (int i=0; i<tot; i++)
	if ( s[i] == 'a' || s[i] == 'c' || s[i] == 'g' || s[i] == 't' ) ++n;
    }
  
  return true;
}

bool SeqDBase::GC( const Region & region , int & gc , int & tot ) 
{
  std::string s = lookup( region );
  gc = 0;
  tot = s.size();
  if ( tot == 0 ) return false;
  for (int i=0; i<tot; i++)
    if ( s[i] == 'G' || s[i] == 'C' ) ++gc;
  return true;
}

bool SeqDBase::ACGT( const Region & region , int & a , int & c , int & g , int & t , int & n ) 
{
  std::string s = lookup( region );
  a = c = g = t = n = 0;
  int tot = s.size();
  if ( tot == 0 ) return false;
  for (int i=0; i<tot; i++)
    if ( s[i] == 'A' ) ++a;
    else if ( s[i] == 'C' ) ++c;
    else if ( s[i] == 'G' ) ++g;
    else if ( s[i] == 'T' ) ++t;
    else ++n;
  return true;
}

void SeqDBase::setMinMax()
{

  chrminmax.clear();
  sqlite3_stmt * s = 
    sql.prepare(" SELECT chr,min(bp1),max(bp2) FROM refseq GROUP BY chr; ");    
  while ( sql.step( s ) )
    {
      chrminmax.insert(make_pair( sql.get_int(s,0) , 
				  int2( sql.get_int(s,1), sql.get_int(s,2) )));
    }
  sql.finalise( s );    
}    



std::string SeqDBase::summary( bool ugly )
{

  std::stringstream ss;


  if ( ! ugly ) ss << "---Sequence DB summary---\n\n";

  std::map<int,int2>::iterator i = chrminmax.begin();
  while ( i != chrminmax.end() )
    {
      int2 d = i->second;
      
      if ( i->first > 0 ) 
	{
	  if ( ugly ) ss << "SEQDB\t"
			 << "REGION=" << chrCode( i->first ) << ":" 
			 << d.p1 << ".." << d.p2 << "\t"
			 << "MB=" << (d.p2-d.p1)/1000000 
			 << "\n";
	  else
	    ss << chrCode( i->first ) << ":" 
	       << d.p1 << ".." << d.p2 << "\t"
	       << "MB=" << (d.p2-d.p1)/1000000 
	       << "\n";
	}
	
	++i;
    }

  if ( ! ugly ) ss << "\n";
  
  std::map<string,std::string>::iterator j = meta.begin();
  while ( j != meta.end() )
    {
      if ( ugly )
	ss << "SEQDB\t" 
	   << j->first << "\t"
	   << j->second << "\n";
      else
	ss << "SEQDB meta-information: " 
	   << j->first << " = "
	   << j->second << "\n";
      
      ++j;
    }
  
    return ss.str();
}


void SeqDBase::dump( const Region & region , bool compact )
{
  
  std::string s = lookup( region );
  
  std::string chr = Helper::chrCode( region.chromosome() );
  
  int bp = region.start.position();
  
  if ( compact ) 
    {
      for (int i=0; i<s.size(); i++)
	plog << s[i];
      plog << "\n";
    }
  else
    for (int i=0; i<s.size(); i++)
      {
	plog << chr << "\t"
	     << bp++ << "\t"
	     << s[i] << "\n";
      }
  
  return;
}
