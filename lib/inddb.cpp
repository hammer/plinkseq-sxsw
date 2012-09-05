
#include "plinkseq/inddb.h"
#include "plinkseq/filemap.h"
#include "plinkseq/gstore.h"

extern GStore * GP;

void IndDBase::wipe( const std::string & n )
{
  if ( Helper::fileExists(n) ) Helper::remove_file( n );  
}

bool IndDBase::new_db( const std::string & n )
{
  
  if ( Helper::fileExists(n) ) Helper::remove_file( n );
  
  sql.open(n);
  sql.synchronous(false);
  sql.query("PRAGMA encoding='UTF-8'");
    
       
  // Individual store
    
  sql.query( " CREATE TABLE IF NOT EXISTS individuals("
	     "   indiv_id INTEGER PRIMARY KEY , "
	     "   name     VARCHAR(20) NOT NULL , "
	     "   fam_id   VARCHAR(20) , "
	     "   ind_id   VARCHAR(20) , "
	     "   pat_id   VARCHAR(20) , "
	     "   mat_id   VARCHAR(20) , "
	     "   sex      CHAR(1) ); " );
  
  //	       "   CONSTRAINT uniqID UNIQUE ( name ) ); " );
  
  // Individual/phenotypic meta-information
  
  sql.query( " CREATE TABLE IF NOT EXISTS phenotypes("
	     "    indiv_id  INTEGER NOT NULL , "
	     "    pheno_id  INTEGER NOT NULL , "
	     "    value     REAL  ) ; " );
  
  //               " CONSTRAINT uniqval UNIQUE ( indiv_id,pheno_id ) ); ");
  
    
  // Phenotype meta-information
  
  sql.query( " CREATE TABLE IF NOT EXISTS metaphenotypes("
	     "    pheno_id    INTEGER PRIMARY KEY,"
	     "    type        CHAR(1) , "
	     "    name        VARCHAR(40) NOT NULL,"
	     "    missing     VARCHAR(1) , "
	     "    description TEXT , "
  	     " CONSTRAINT uniqphe UNIQUE ( name ) ); " );
  
  // Prepate some key queries
  
  index();
    
  init();
  
  //plog << "Created individual-database: " << n << "\n";
  
  return true;
  
}


bool IndDBase::attach( const std::string & n )
{

  if ( n == "-" || n == "." ) { dettach(); return false; } 

  if ( ! Helper::fileExists(n) )
    {
      new_db(n);
      return true;
    }
  
  sql.open(n);
  
  index();
  
  init();
  
  //
  // Load new meta-information table
  //
  
  set_metatypes( );
  
  return true;
}


void IndDBase::set_metatypes(bool clear )
{

  if ( clear ) 
    {
      MetaInformation<IndivMeta>::reset();    
    }
  
  std::map<std::string, std::vector<std::string> > pi = fetch_phenotype_info();
  std::map<std::string, std::vector<std::string> >::iterator i = pi.begin();
  
  while ( i != pi.end() )
    {
      std::string name = i->first;
      std::string type = i->second[0];
      std::string desc = i->second[1];
      
      if ( Helper::is_int( type ) ) registerMetatype( name,  META_INT , 1, META_GROUP_INDIV, desc );
      else if ( Helper::is_float( type ) ) registerMetatype( name,  META_FLOAT , 1, META_GROUP_INDIV, desc );
      else registerMetatype( name,  META_TEXT , 1, META_GROUP_INDIV, desc );	
      ++i;
    }  

}


bool IndDBase::dettach()
{
  if ( attached() )
    {
      release();
      sql.close();
    }
  return true;
}


bool IndDBase::init()
{
  
    // 
    // Insertions
    // 
    
    stmt_insert_individual = 
	sql.prepare(" INSERT OR REPLACE INTO individuals "
		    "          ( name, fam_id, ind_id, pat_id, mat_id, sex ) "
		    "   values ( :name, :fid, :iid, :pat, :mat, :sex ); " );

    stmt_update_individual = 
	sql.prepare(" UPDATE individuals "
		    " SET fam_id = :fid, ind_id = :iid, pat_id = :pat, mat_id = :mat , sex = :sex "
		    " WHERE name == :name ; ");

    stmt_replace_individual_id = 
      sql.prepare( "UPDATE individuals SET name = :new_id WHERE name == :old_id; ");
    
    stmt_insert_phenotype = 
	sql.prepare(" INSERT OR REPLACE INTO phenotypes ( indiv_id , pheno_id , value ) "
		    " values ( :indiv_id , :pheno_id , :value ); " );

    stmt_insert_metaphenotype =
	sql.prepare( " INSERT OR REPLACE INTO metaphenotypes ( pheno_id , type , name , missing, description ) "
		     " values ( :pheno_id , :type , :name , :missing, :description ) ; ");    
    
    
    //
    // Queries
    //
    
    stmt_fetch_individuals = 
	sql.prepare(" SELECT indiv_id , name FROM individuals ORDER BY indiv_id; ");

    stmt_lookup_id = 
	sql.prepare(" SELECT indiv_id FROM individuals WHERE name == :name; ");

    stmt_lookup_pheno_id = 
	sql.prepare(" SELECT pheno_id FROM metaphenotypes WHERE name == :name; ");

    stmt_fetch_sex = 
	sql.prepare(" SELECT sex FROM individuals WHERE name == :name; ");

    stmt_fetch_phenotype_list = 
	sql.prepare(" SELECT pheno_id,name,type,missing,description "
		    " FROM metaphenotypes; ");
    
    stmt_fetch_individual = 
	sql.prepare(" SELECT * FROM individuals WHERE indiv_id == :indiv_id ; ");
    
    stmt_fetch_phenotype_values = 
	sql.prepare(" SELECT mp.name,mp.type,p.value FROM phenotypes AS p , metaphenotypes AS mp "
		    " WHERE p.pheno_id == mp.pheno_id "
		    "   AND indiv_id == :indiv_id ; ");
    
    
    return true;
}

bool IndDBase::release() 
{
    
  sql.finalise( stmt_insert_individual  ); 
  sql.finalise( stmt_update_individual  ); 
  sql.finalise( stmt_replace_individual_id );

  sql.finalise( stmt_insert_phenotype  ); 
  sql.finalise( stmt_insert_metaphenotype  );
  
  sql.finalise( stmt_fetch_individuals  ); 
  sql.finalise( stmt_lookup_id  ); 
  sql.finalise( stmt_lookup_pheno_id  ); 
  sql.finalise( stmt_fetch_individual  ); 
  
  sql.finalise( stmt_fetch_phenotype_values );
  sql.finalise( stmt_fetch_phenotype_list );
  
  return true;
}

int IndDBase::fetch_id(const std::string & name)
{
  if ( ! attached() ) return 0;
  int code = 0;
  sql.bind_text( stmt_lookup_id , ":name" , name );
  if ( sql.step( stmt_lookup_id ) )
    {
      code = sql.get_int( stmt_lookup_id , 0 );
    }
  sql.reset( stmt_lookup_id );
  return code;
}

int IndDBase::fetch_pheno_id( const std::string & name)
{
  if ( ! attached() ) return 0;
  int code = 0;
  sql.bind_text( stmt_lookup_pheno_id , ":name" , name );
  if ( sql.step( stmt_lookup_pheno_id ) )
    code = sql.get_int( stmt_lookup_pheno_id , 0 );
  sql.reset( stmt_lookup_pheno_id );
  return code;
}

bool IndDBase::index()
{
  sql.query("CREATE INDEX IF NOT EXISTS ind1 ON individuals(name); ");
  sql.query("CREATE INDEX IF NOT EXISTS phe1 ON phenotypes(indiv_id); ");
  sql.query("CREATE INDEX IF NOT EXISTS phe2 ON phenotypes(pheno_id); ");
}

bool IndDBase::drop_index()
{
  sql.query("DROP INDEX IF EXISTS ind1;");
  sql.query("DROP INDEX IF EXISTS phe1;");
  sql.query("DROP INDEX IF EXISTS phe2;");
}


uint64_t IndDBase::insert_phenotype( const std::string & name , 
				     const std::string & type , 
				     const std::string & miss ,
				     const std::string & desc )
{
  
  sql.bind_text( stmt_insert_metaphenotype , ":name" , name );
  sql.bind_text( stmt_insert_metaphenotype , ":type" , type );
  sql.bind_text( stmt_insert_metaphenotype , ":missing" , miss );
  sql.bind_text( stmt_insert_metaphenotype , ":description" , desc );
  
  sql.step( stmt_insert_metaphenotype );	    
  sql.reset( stmt_insert_metaphenotype ); 
  
  if ( Helper::is_int( type ) ) registerMetatype( name,  META_INT , 1, META_GROUP_INDIV, desc );
  else if ( Helper::is_float( type ) ) registerMetatype( name,  META_FLOAT , 1, META_GROUP_INDIV, desc );
  else registerMetatype( name,  META_TEXT , 1, META_GROUP_INDIV, desc );	
      
  return fetch_pheno_id( name );
  
}


uint64_t IndDBase::insert( const Individual & ind , bool * newone )
{

  int id = fetch_id( ind.id() );
  
  if ( id == 0 ) 
    {
      
      // This is a new record
      
      sql.bind_text( stmt_insert_individual , ":name" , ind.id() );
      sql.bind_text( stmt_insert_individual , ":fid" , ind.fid() );
      sql.bind_text( stmt_insert_individual , ":iid" , ind.iid() );
      sql.bind_text( stmt_insert_individual , ":pat" , ind.father() );
      sql.bind_text( stmt_insert_individual , ":mat" , ind.mother() );
      sql.bind_text( stmt_insert_individual , ":sex" , Helper::int2str( ind.sex() ) ); 
      
      sql.step( stmt_insert_individual );
      sql.reset( stmt_insert_individual );
      
      if ( newone ) *newone = true;      
      id = fetch_id( ind.id() );
    }
  else
    {
      // ..append an existing record
      
      sql.bind_text( stmt_update_individual , ":name" , ind.id() );
      sql.bind_text( stmt_update_individual , ":fid" , ind.fid() );
      sql.bind_text( stmt_update_individual , ":iid" , ind.iid() );
      sql.bind_text( stmt_update_individual , ":pat" , ind.father() );
      sql.bind_text( stmt_update_individual , ":mat" , ind.mother() );
      sql.bind_text( stmt_update_individual , ":sex" , Helper::int2str( ind.sex() ) );
      
      sql.step( stmt_update_individual );
      sql.reset( stmt_update_individual );
      if( newone ) *newone = false;      
    }  
  return id;
}


void IndDBase::insert( const uint64_t i, const uint64_t p , const int x )
{  
  sql.bind_int( stmt_insert_phenotype , ":indiv_id" , i );
  sql.bind_int( stmt_insert_phenotype , ":pheno_id" , p );
  sql.bind_int( stmt_insert_phenotype , ":value" , x );
  sql.step( stmt_insert_phenotype );
  sql.reset( stmt_insert_phenotype );
}
  
void IndDBase::insert( const uint64_t i, const uint64_t p , const double x )
{  
  sql.bind_int( stmt_insert_phenotype , ":indiv_id" , i );
  sql.bind_int( stmt_insert_phenotype , ":pheno_id" , p );
  sql.bind_double( stmt_insert_phenotype , ":value" , x );
  sql.step( stmt_insert_phenotype );
  sql.reset( stmt_insert_phenotype );
}

void IndDBase::insert( const uint64_t i, const uint64_t p , const std::string & x )
{  
  sql.bind_int( stmt_insert_phenotype , ":indiv_id" , i );
  sql.bind_int( stmt_insert_phenotype , ":pheno_id" , p );
  sql.bind_text( stmt_insert_phenotype , ":value" , x );
  sql.step( stmt_insert_phenotype );
  sql.reset( stmt_insert_phenotype );
}


bool IndDBase::load_ped_info( const std::string & filename )
{

  if ( ! attached() ) Helper::halt( "no attached INDDB" );
  
  if ( ! Helper::fileExists(filename) ) 
    {
      plog.warn("could not find pedigree file " + filename );
      return false;
    }

  // Expect this fixed format:
  // A "." indicates that no FID/IID is specified
  // "." or "0" similar
  
  // ID  FID  IID  PAT  MAT  SEX 
  
  InFile f( filename );
  
  int cnt1 = 0 , cnt2 = 0;
  
  // place entire load in a single transaction
  sql.begin();
  
  while ( ! f.eof() )
    {
      
      std::vector<std::string> tok = f.tokenizeLine();

      const int sz = tok.size();
      if ( sz == 0 ) continue;
      
  	  const std::string& firstChar = tok[0].substr(0,1);
  	  if (firstChar == "#") continue; // ignore comment lines
      
      if (tok.size() < 6) { // allow for extra meta-information fields
    	  plog.warn("found line in pedigree file with less than 6 tab-delimited fields");
    	  continue;
      }
      
      Individual ind( tok[0] );
      
      ind.fid( tok[1] );
      ind.iid( tok[2] );
      ind.pat( tok[3] );
      ind.mat( tok[4] );
      ind.sex( tok[5] );

      bool newone;

      insert( ind , &newone );

      if ( newone ) ++cnt1;
      else ++cnt2;
      
    }
  
  f.close();

  sql.commit();
  
  plog << "Inserted " << cnt1 
       << " new individuals, updated " << cnt2 
       << " existing individuals\n";

  if ( cnt1 + cnt2 ) 
    {
      if ( GP && GP->has_project_file() )
	GP->fIndex.append_to_projectfile( Helper::fullpath(filename) , "PED" );      
    }

  return true;
  
}


bool IndDBase::load_phenotypes( const std::string & filename )
{

  if ( ! attached() ) Helper::halt( "no attached INDDB" );

  if ( ! Helper::fileExists(filename) ) 
    {
      plog.warn("could not find phenotype file " + filename );
      return false;
    }

  // Expect a variant of FAM file format
  
  InFile f( filename );

  std::map<std::string,int> phe_codes1;
  std::map<int,int> phe_codes2;
  std::vector<std::string> phe_codes;
  std::map<std::string,mType> type_codes;
  
  std::map<std::string,std::string> mis_codes;
  
  int expected_col_count = -1;
  int inserted = 0;
  
  sql.begin();

  drop_index();
 
  while ( ! f.eof() )
    {
      
      std::string s = f. readLine();
      
      if ( s == "" ) continue;
      
      // Meta-information? 
      
      if ( s.size() > 2 && s.substr(0,2) == "##" )
	{
	  
	  std::vector<std::string> tok = Helper::quoted_parse( s.substr(2) );
	  
	  // need at least a name
	  if ( tok.size() < 1 ) continue;	  
	  std::string name = tok[0];

	  // defaults
	  std::string type = "Integer"; 
	  std::string miss = "0";
	  std::string desc = "Dichotomous phenotype";	  
	  if ( tok.size() >= 2 ) type = tok[1];
	  if ( tok.size() >= 3 ) miss = tok[2];
	  if ( tok.size() >= 4 ) desc = tok[3];
	  
	  int code = insert_phenotype( name, type, miss, desc );

	  phe_codes1[ name ] = code;
	  mis_codes[ name ] = miss;
	  
	  if ( Helper::is_int( type ) ) type_codes[ name ] = META_INT;
	  else if ( Helper::is_float( type ) ) type_codes[ name ] = META_FLOAT;
	  else type_codes[ name ] = META_TEXT;

	}

      // Or header line?
      
      else if ( s.substr(0,1) == "#" )
	{
	  
	  // #ID phe1 phe2 phe3
	  
	  std::vector<std::string> tok = Helper::parse( s , " \t");
	  
	  if ( tok.size() < 2 ) { plog.warn("malformed phenotype file"); continue; } 
	  if ( tok[0] != "#ID" ) { plog.warn("malformed phenotype file"); continue; } 
	  
	  for ( int i = 1 ; i < tok.size(); i++ )
	    {
	      std::map<std::string,int>::iterator k = phe_codes1.find( tok[i] );	      

	      if ( k == phe_codes1.end() ) 
		Helper::halt( tok[i] + " in header of phenotype file but not defined" );
	      
	      phe_codes2[ i ] = k->second;
	      phe_codes.push_back( tok[i] );
	      
	    }
	  expected_col_count = tok.size();
	}
      
      // Or data ? 
      
      else 
	{
	  
	  // Skip, if we haven't seen a header
	  if ( expected_col_count == -1 ) continue;
	 
	  std::vector<std::string> tok = Helper::parse( s , " \t");
	  
	  if ( tok.size() != expected_col_count ) 
	    {
	      plog.warn("row in phenotype file with wrong number of fields");
	      continue;
	    }
	  
	  int indiv_id = fetch_id( tok[0] );
	  
	  // if individual does not exist, create

	  if ( indiv_id == 0 ) 
	    {
	      std::string period = ".";

	      sql.bind_text( stmt_insert_individual , ":name" , tok[0] );
	      sql.bind_text( stmt_insert_individual , ":fid" , period );
	      sql.bind_text( stmt_insert_individual , ":fid" , period );
	      sql.bind_text( stmt_insert_individual , ":iid" , period );
	      sql.bind_text( stmt_insert_individual , ":pat" , period );
	      sql.bind_text( stmt_insert_individual , ":mat" , period );
	      sql.bind_text( stmt_insert_individual , ":sex" , period );
	      
	      sql.step( stmt_insert_individual );
	      sql.reset( stmt_insert_individual );
	      
	      // and grab the ID
	      indiv_id = fetch_id( tok[0] );

	    }
	  
	  //
	  // Insert actual phenotypes
	  //
	  
	  for ( int i = 1; i < tok.size(); i++ )
	    {

	      // skip undefined phenotypes
	      if ( phe_codes2[i] == 0 ) continue;

	      // skip missing values
	      if ( tok[i] == mis_codes[ phe_codes[i-1] ] ) continue;
	      
	      mType mt = type_codes[ phe_codes[ i-1 ] ];
	      
	      // skip invalid values for numerics (as MT will be registered)
	      
 	      if ( mt == META_INT )
 		{    
 		  int x;
 		  if ( Helper::str2int( tok[i] , x ) )
 		    insert( indiv_id , phe_codes2[i] , x );
 		}
 	      else if ( mt == META_FLOAT )
 		{
 		  double x;
 		  if ( Helper::str2dbl( tok[i] , x ) )
 		    insert( indiv_id , phe_codes2[i] , x );
 		}
 	      else 
 		{
 		  insert( indiv_id , phe_codes2[i] , tok[i] );
 		}


	    }
	  ++inserted;	    
	}
    }
  
  f.close();

  index();
  
  sql.commit();

  plog << "Processed " << inserted << " rows\n";
  
  if ( inserted && GP && GP->has_project_file() ) 
    GP->fIndex.append_to_projectfile( Helper::fullpath( filename ) , "PHE" );
    
  return true;
}


//
// Fetch all individuals
//

std::vector<Individual> IndDBase::fetch( )
{
  std::vector<Individual> inds;
  if ( ! attached() ) return inds;   
  while ( sql.step( stmt_fetch_individuals ) )
    {	
      inds.push_back( fetch( sql.get_int64( stmt_fetch_individuals , 0 ) ) );
    }
  sql.reset( stmt_fetch_individuals );
  return inds;
}

bool IndDBase::fetch( Individual * person )
{

  if ( ! attached() ) return false;

  // Look-up this individual based upon their ID. 
  // Attach all phenotype and meta-information we find
  // Return true if find somebody matching that description
  
  uint64_t idx = fetch_id( person->id() );
 
  if ( idx == 0 ) return false;
  
  if ( fetch( person , fetch_id( person->id() ) ) )
    {
      person->missing( false );
      return true;
    }
  else
    {
      person->missing( true );
      return false;
    }

}

Individual IndDBase::fetch( uint64_t indiv_id )
{
  Individual ind;
  if ( ! fetch( &ind, indiv_id ) ) ind.missing( true );
  return ind;
}


bool IndDBase::fetch( Individual * person , uint64_t indiv_id )
{

  sql.bind_int64( stmt_fetch_individual , ":indiv_id" , indiv_id );
  
  bool obs = false;
  
  if ( sql.step( stmt_fetch_individual ) )
    {
      obs = true;
      person->idx( sql.get_int(  stmt_fetch_individual , 0 ) );
      person->id( sql.get_text(   stmt_fetch_individual , 1 ) );
      person->fid( sql.get_text(  stmt_fetch_individual , 2 ) );
      person->iid( sql.get_text(  stmt_fetch_individual , 3 ) );
      person->pat( sql.get_text(  stmt_fetch_individual , 4 ) );
      person->mat( sql.get_text(  stmt_fetch_individual , 5 ) );

      int s = sql.get_int(  stmt_fetch_individual , 6 ) ;
      if ( s == 1 ) 
	person->sex( MALE );
      else if ( s == 2 ) 
	person->sex( FEMALE );
      else 
	person->sex( UNKNOWN_SEX );      
    }
  
  sql.reset( stmt_fetch_individual );
  
  // Get any phenotype information
  
  sql.bind_int64( stmt_fetch_phenotype_values , ":indiv_id" , indiv_id );
  while ( sql.step( stmt_fetch_phenotype_values ) )
    {
      obs = true;

      std::string pheno = sql.get_text( stmt_fetch_phenotype_values , 0 );
      std::string type = sql.get_text( stmt_fetch_phenotype_values , 1 );
      
      if ( Helper::is_int( type ) ) person->meta.set( pheno , sql.get_int( stmt_fetch_phenotype_values , 2 ) );
      else if ( Helper::is_float( type ) ) person->meta.set( pheno , sql.get_double( stmt_fetch_phenotype_values , 2 ) );
      else person->meta.set( pheno , sql.get_text( stmt_fetch_phenotype_values , 2 ) );

    }
  sql.reset( stmt_fetch_phenotype_values );
  
  return obs;
}



void IndDBase::load_meta( std::vector<Individual> & ind, const std::string & name )
{

  // this function used by rint.cpp
  
  // we should probably scrap this and get a single framework for
  // appending phenotypic info.  Might be useful to be able to not
  // append *all* phenotypic info. available however (i.e. if v. large
  // phenotype datasets attached.

  int pcode = fetch_pheno_id( name );
  if ( pcode == 0 ) return;
  
  sql.begin();
  
  for ( int i = 0 ; i < ind.size(); i++)
    {
      
      sql.bind_int64( stmt_fetch_phenotype_values , ":indiv_id" , ind[i].idx() );
      
      while ( sql.step( stmt_fetch_phenotype_values ) )
	{
	  std::string pheno = sql.get_text( stmt_fetch_phenotype_values , 0 );
	  std::string type = sql.get_text( stmt_fetch_phenotype_values , 1 );
	  
	  if ( Helper::is_int( type ) )
	    ind[i].meta.set( pheno , sql.get_int( stmt_fetch_phenotype_values , 2 ) );	    
	  else if ( Helper::is_float( type ) )
	    ind[i].meta.set( pheno , sql.get_double( stmt_fetch_phenotype_values , 2 ) );
	  else
	    ind[i].meta.set( pheno , sql.get_text( stmt_fetch_phenotype_values , 2 ) );	  
	}
      sql.reset( stmt_fetch_phenotype_values );
    }  
    sql.commit();        
}


int IndDBase::size()
{
  std::vector<int> r = 
    sql.intTable( "SELECT count(*) FROM individuals;" , 1 );
  if ( r.size() == 1 ) 
    return r[0];
  return -1;    
}

std::string IndDBase::summary( bool ugly )
{
  std::stringstream ss;

  if ( ugly ) 
    ss << "INDDB\t"
       << "N=" << size() << "\n";
  else
    {
      ss << "---Individual DB summary---\n\n";

      ss << size() << " unique individuals\n";
      
      std::map<std::string,std::vector<std::string> > t = fetch_phenotype_info();
      
      std::map<std::string,std::vector<std::string> >::iterator p = t.begin();

      while ( p != t.end() )
	{
	  ss << "Phenotype : " << p->first << " "
	     << "(" << p->second[0]  << ") "
	     << p->second[1] << "\n";
	  ++p;
	}      
    }

  return ss.str();
}
  


std::map<std::string,std::vector<std::string> > IndDBase::fetch_phenotype_info()
{
  std::map<std::string,std::vector<std::string> > m;
  while ( sql.step( stmt_fetch_phenotype_list ) )
    {
      std::string name = sql.get_text(  stmt_fetch_phenotype_list , 1 ) ;      
      std::vector<std::string> type_desc;
      type_desc.push_back( sql.get_text(  stmt_fetch_phenotype_list , 2 ) );      
      // 3 = missing data code
      type_desc.push_back( sql.get_text(  stmt_fetch_phenotype_list , 4 ) );      
      m[name] = type_desc;
    } 
  sql.reset( stmt_fetch_phenotype_list );
  return m;
}


sType IndDBase::sex( const std::string & id ) 
{  
  sql.bind_text( stmt_fetch_sex , ":name" , id );
  if ( sql.step( stmt_fetch_sex ) )
    {
      int s = sql.get_int( stmt_fetch_sex , 0 );
      sql.reset( stmt_fetch_sex );
      return s == 2 ? FEMALE : ( s == 1 ? MALE : UNKNOWN_SEX ) ;
    }
  sql.reset( stmt_fetch_sex );
  return UNKNOWN_SEX;
}


bool IndDBase::replace_individual_id( const std::string & old_id , const std::string & new_id )
{
  sql.bind_text( stmt_replace_individual_id , ":old_id" , old_id );
  sql.bind_text( stmt_replace_individual_id , ":new_id" , new_id );
  sql.step( stmt_replace_individual_id );
  sql.reset( stmt_replace_individual_id );
  return true;
}
