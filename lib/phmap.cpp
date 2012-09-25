
#include "plinkseq/phmap.h"
#include "plinkseq/indmap.h"
#include "plinkseq/inddb.h"
#include "plinkseq/filemap.h"
#include "plinkseq/gstore.h"

extern GStore * GP;

PhenotypeMap::PhenotypeMap(IndDBase * pinddb) 
{   
  phenotype_name = ".";
  phenotype_type = PHE_NONE;
  use_strata = false;
  strata_name = ".";
  inddb = pinddb;
}


PhenotypeMap::~PhenotypeMap()
{
  reset();
}

void PhenotypeMap::begin()
{ 
  if ( inddb && inddb->attached() ) { inddb->begin(); } 
} 

void PhenotypeMap::commit()
{ 
  if ( inddb && inddb->attached() ) { inddb->commit(); }
} 

void PhenotypeMap::reset()
{
  // Separate copies of individuals will have been made, so free those now  
  std::map<std::string,Individual*>::const_iterator i = pmap.begin();
  while ( i != pmap.end() )
    {
      delete i->second;
      ++i;
    }  
  pmap.clear();
  phenotype_name = ".";
  phenotype_type = PHE_NONE;  
  use_strata = false;
  strata_name = ".";
}

int PhenotypeMap::set_strata( const std::string & s )
{

  strata_name = ".";
  use_strata = false;
  
  if ( s == "" || s == "." ) 
    {
      plog.warn("no valid strata specified");
      return 0;
    }

  // accumlate strata code across 1+ vars
  std::map<std::string,std::string> smap;
  
  std::vector<std::string> svec = Helper::char_split( s , ',' );
  for (int j = 0 ; j < svec.size(); j++ )
    {
      mType mt = MetaInformation<IndivMeta>::type( svec[j] ) ;
  
      if ( ! ( mt == META_TEXT || mt == META_INT || mt == META_CHAR ) ) 
	{
	  plog.warn("strata arg(s) must be text or integer");
	  return 0;
	}
      
      strata_name = svec[j];
      
      std::map< std::string, Individual* >::const_iterator i = pmap.begin();
      while( i != pmap.end() )
	{
	  
	  Individual * person = i->second;
	  
	  if ( ! person->meta.has_field( strata_name ) ) 
	    {
	      if ( smap[ i->first ] == "" ) 
		smap[ i->first ] = ".";
	    }
	  else if ( mt == META_TEXT || mt == META_CHAR ) 
	    {
	      if ( smap[ i->first ]  != "" ) 
		smap[ i->first ] += ",";
	      smap[ i->first ] += person->meta.get1_string( strata_name );
	    }
	  else if ( mt == META_INT )
	    {
	      if ( smap[ i->first ]  != "" ) 
		smap[ i->first ] += ",";
	      smap[ i->first ] += Helper::int2str( person->meta.get1_int( strata_name ) );
	    }
	  
	  ++i;
	}
    }


  //
  // Set actual groupings
  //

  int nonmissing = 0;
  std::map< std::string, Individual* >::const_iterator i = pmap.begin();
  while( i != pmap.end() )
    {
      
      Individual * person = i->second;
      person->group( smap[ i->first ] );
      if ( smap[ i->first ] != "." ) 
	++nonmissing;
      ++i;
    }
  
  if ( nonmissing ) 
    {
      use_strata = true;
      strata_name = s; // reset to whole name, if comma-delim
    }
  else 
    {
      use_strata = false;
      strata_name = ".";
    }
  return nonmissing;

}

void PhenotypeMap::align( const std::set<std::string> & ids )
{
  
  // Given a set of IDs from an IndividualMap, ensure that we align 

  // By definition, everybody in the indiv-map should already be in 
  // the pheno-map, but just in case...

  std::set<std::string>::iterator i = ids.begin();
  while ( i != ids.end() )
    {      
      new_individual( *i );  // does nothing if person already exists
      ++i;
    }
  
  // Now, in the other direction
  std::map< std::string, Individual *>::iterator j = pmap.begin();
  while ( j != pmap.end() )
    {
      if ( ids.find( j->first ) == ids.end() )
	{
	  delete j->second;
	  pmap.erase(j++);
	}
      else
	++j;
    }
}

bool PhenotypeMap::phenotype_exists( const std::string & phenotype ) const
{
  // catch case of on-the-fly created variable:
  if ( phenotype == phenotype_name ) return true;
  mType mt = MetaInformation<IndivMeta>::type( phenotype ) ;
  return mt != META_UNDEFINED;  
}

pType PhenotypeMap::type( const std::string & phenotype ) const 
{

  // catch case of on-the-fly created variable:
  // (i.e. residualized phenotypes are stored here
  if ( phenotype == phenotype_name ) return phenotype_type;

  // otherwise, look up from meta-info
  mType mt = MetaInformation<IndivMeta>::type( phenotype ) ;
  if ( mt == META_INT ) return PHE_DICHOT;
  if ( mt == META_FLOAT ) return PHE_QT;
  if ( mt == META_TEXT ) return PHE_FACTOR;
  return PHE_NONE;  
}

int PhenotypeMap::set_phenotype( const std::string & phenotype )
{
  
  int nonmissing = 0;
  
  phenotype_name = phenotype;
  
  mType mt = MetaInformation<IndivMeta>::type( phenotype ) ;
  
  if ( mt == META_UNDEFINED ) 
    Helper::halt( "could not find phenotype " + phenotype );

  if ( mt == META_INT )     
    phenotype_type = PHE_DICHOT; 
  else if ( mt == META_FLOAT ) 
    phenotype_type = PHE_QT;
  else 
    phenotype_type = PHE_FACTOR;
    
  std::map< std::string, Individual* >::const_iterator i = pmap.begin();
  while( i != pmap.end() )
    {
      
      Individual * person = i->second;
      
      if ( phenotype_type == PHE_DICHOT ) 
	{
	  if ( ! person->meta.has_field( phenotype ) )
	    {
	      person->missing( true );
	      person->affected( UNKNOWN_PHE );
	    }
	  else if ( person->meta.get1_int( phenotype ) == 2 )
	    {
	      person->affected( CASE );
	      person->missing( false );
	      nonmissing++;
	    }
	  else if ( person->meta.get1_int( phenotype ) == 1 )
	    {
	      person->affected( CONTROL );
	      person->missing( false );
	      nonmissing++;
	    }
	  else
	    {
	      person->affected( UNKNOWN_PHE );
	      person->missing( true );
	    }
	}
      
      else if ( phenotype_type == PHE_QT )
	{
	  if ( ! person->meta.has_field( phenotype ) )
	    person->missing( true );
	  else
	    {	      
	      person->missing( false );
	      person->qt( person->meta.get1_double( phenotype ) );
	      nonmissing++;				  
	    }
	}
      else 
	{
	  if ( ! person->meta.has_field( phenotype ) )
	    {
	      person->missing( true );
	      person->group(0);
	    }
	  else
	    {
	      person->missing( false );
	      person->group( person->meta.get1_string( phenotype ) );
	      nonmissing++;				  
	    }
	}
      
      ++i;
    }
  
  phenotype_name = phenotype;

  return nonmissing;
}

int PhenotypeMap::attach_dichot_phenotype( const std::string & pname , const std::vector<int> & phe , const IndividualMap & imap )
{

  int nonmissing = 0;
  phenotype_name = pname;
  phenotype_type = PHE_DICHOT;

  // register as a type if does not exist
  MetaInformation<IndivMeta>::field( pname , META_INT , 1 , "." );
  
  const int n = imap.size();
    for (int i = 0 ; i < n ; i++ ) 
    {
      
      Individual * person = imap(i);
      
      // Is this person a case? 
      if ( phe[i] == 2 ) 
	{
	  person->affected( CASE );
	  person->meta.set( pname , 2 );
	  ++nonmissing;
	}
      else if ( phe[i] == 1 ) 
	{
	  // A control
	  person->affected( CONTROL );
	  person->meta.set( pname , 1 );
	  ++nonmissing;
	}
      else
	{
	  person->affected( UNKNOWN_PHE );		  
	  person->meta.set( pname , 0 );
	}
    }
  return nonmissing;  
}


int PhenotypeMap::attach_qt_phenotype( const std::string & pname , 
				       const std::vector<bool> & missing , 
				       const std::vector<double> & phe , 
				       const IndividualMap & imap )
{
  
  int nonmissing = 0;
  phenotype_name = pname;
  phenotype_type = PHE_QT;

  // register as a type if does not exist
  MetaInformation<IndivMeta>::field( pname , META_FLOAT , 1 , "." );
  
  const int n = imap.size();

  for (int i = 0 ; i < n ; i++ ) 
    {
      
      Individual * person = imap(i);
      
      // Is this person a case? 
      if ( ! missing[i] ) 
	{
	  person->qt( phe[i] );
	  person->meta.set( pname , phe[i] );
	  person->missing( false );
	  ++nonmissing;
	}
      else
	{
	  person->missing( true );	  
	}
    }
  return nonmissing;  
}



int PhenotypeMap::make_phenotype( const std::string & make_phenotype )
{
  
  int nonmissing = 0;

  // expect in format: GRP=L1,L2
  //                   GRP=L1,L2:L3
  // where GRP should be FACTOR or INT

  std::vector<std::string> p = Helper::char_split( make_phenotype , '=' );
  
  // Well-formed specification?
  if ( p.size() != 2 ) 
    {
      plog.warn("make-phenotype arg not well formed (" , make_phenotype );
      return 0;
    }

  // Can we find a phenotype?
  if ( set_phenotype( p[0] ) == 0 ) 
    {
      plog.warn("could not find phenotype values for", p[0] );
      return 0;
    }

  
  // Is this a factor (or a dichot, for a simple swap)?
  
  if ( type() != PHE_FACTOR && type() != PHE_DICHOT ) 
    {
      plog.warn("make-phenotype arg must be a factor");
      return 0;
    }


  std::vector<std::string> p2 = Helper::char_split( p[1] , ':' );
  
  if ( p2.size() != 1 && p2.size() !=2 ) 
    {
      plog.warn("make-phenotype arg not well formed");
      return 0;
    }

  bool explicit_missing = p2.size() == 2 ;
   
  std::set<std::string> grp1;
  std::set<std::string> grp2;
  
  std::vector<std::string> t = Helper::char_split( p2[0] , ',' );
  for (int i=0; i<t.size(); i++) grp1.insert( t[i] );
  

  if ( explicit_missing ) 
    {
      std::vector<std::string> t = Helper::char_split( p2[1] , ',' );
      for (int i=0; i< t.size(); i++) grp2.insert( t[i] );
    }


  //
  // We seem all okay now, so let's make the phenotype
  //

  pType original_type = phenotype_type;
  phenotype_name = make_phenotype;
  phenotype_type = PHE_DICHOT;
  
  std::map< std::string , Individual*>::iterator i = pmap.begin();
  
  while ( i != pmap.end() )
    {

      Individual * person = i->second;
      
      std::string label = original_type == PHE_DICHOT ?
	( person->affected() == CASE ? "2" : ( person->affected() == CONTROL ? "1" : "." ) ) 
	: 
	person->group_label() ;
      
      if ( person->missing() ) 
	{
	  if ( !explicit_missing )  
	    {
	      person->affected( CONTROL );
	
	      ++nonmissing;
	    }
	  else
	    person->affected( UNKNOWN_PHE );
	}
      else
	{

	  // Is this person a case? 
	  if ( grp1.find( label ) != grp1.end() )
	    {
	      person->affected( CASE );
	      ++nonmissing;
	    }
	  else
	    {
	      // A control
	      if ( ! explicit_missing )
		{
		  person->affected( CONTROL );
		  ++nonmissing;
		}
	      else
		{
		  if ( grp2.find( label ) != grp2.end() )
		    {
		      person->affected( CONTROL );
		      ++nonmissing;
		    }
		  else
		    person->affected( UNKNOWN_PHE );		  
		}
	    }	    
	}      
      
      ++i;
    }
  
  return nonmissing;
}



std::map<std::string,int> PhenotypeMap::summarise_phenotype( const std::string & phenotype )
{
  set_phenotype( phenotype );
  return summarise_phenotype();
}

std::map<std::string,int> PhenotypeMap::summarise_phenotype()
{
  
  mType mt = MetaInformation<IndivMeta>::type( phenotype_name ) ;
  
  // Int   -- interpreted as case/control : return CASE/CONTROL/MISSING
  // Float -- QT                          : return NON-MISSING/MISSING
  // Other -- FACTOR                      : return N at each level, MISSING
  
  int missing = 0;
  
  std::map<std::string,int> r;
  std::set<Individual*> seen;

  r["NA"] = 0;
  //  r["OBS"] = 0;
  
  std::map< std::string, Individual* >::const_iterator i = pmap.begin();
  while( i != pmap.end() )
    {
            
      Individual * person = i->second;

      if ( phenotype_type == PHE_DICHOT  )
	{
	  if ( person->affected( ) == CASE )
	    { 
	      r["CASE"]++;
	      //r["OBS"]++;
	    }
	  else if ( person->affected() == CONTROL ) 
	    {
	      r["CONTROL"]++;
	      //r["OBS"]++;
	    }
	  else 
	    r["NA"]++;
	  
	}
      else if ( phenotype_type == PHE_QT )
	{
	  if ( person->missing() ) 
	    r["NA"]++;
	  else
	    r["OBS"]++;
	}
      else if ( phenotype_type == PHE_FACTOR )
	{
	  if ( person->missing() ) 
	    r["NA"]++;
	  else
	    {
	      //r["OBS"]++;
	      r[ person->group_label() ]++;
	    }
	}
      ++i;
    }
  return r;
}


Individual * PhenotypeMap::new_individual( const std::string & id )
{
  
  // Already exists?

  Individual * person = ind(id);

  if ( person ) return person;

  
  // Otherwise, create a new individual

  person = new Individual( id );
  
  
  // Track in the phenotype map
  
  pmap[ id ] = person ;


  // Lookup in INDDB, attaching any phenotypic information that exists
  
  if ( inddb ) inddb->fetch( person );

  
  // And return a pointer to this new person
  
  return person;
    
}

Data::Vector<double> PhenotypeMap::get_pheno( const std::string & p , const IndividualMap & indmap ) const
{

  // this assumes that the phenotype will already been attached as meta-information

  // --> this is indeed current practice, although likely to change when INDDB's are
  //     used to store very large amounts of phenotype data. In this case, these functions
  //     should secondarily perform a direct lookup from the INDDB.

  const int n = indmap.size();

  Data::Vector<double> d( n );
  
  for (int r=0; r<n; r++)
    {
      Individual * person = indmap(r);      
      if ( person->meta.has_field( p ) )
	{
	  mType mt = MetaInformation<IndivMeta>::type( p );
	  if      ( mt == META_INT ) d(r) = person->meta.get1_int( p );
	  else if ( mt == META_FLOAT ) d(r) = person->meta.get1_double( p );
	  else if ( mt == META_BOOL ) d(r) = person->meta.get1_bool( p );
	  else d.set_elem_mask( r );
	}
      else 
	d.set_elem_mask( r );      
    }
  return d;
}


std::string PhenotypeMap::phenotype(const std::string & p , const int i ) const
{
  // return for a single individual a printable version of the phenotype, or '.' if it does
  // not exist

  const int n = GP->indmap.size();
  if ( i < 0 || i >= n ) return ".";
  Individual * person = GP->indmap(i);
  if ( person->meta.has_field( p ) )
    {
      mType mt = MetaInformation<IndivMeta>::type( p );      
      if      ( mt == META_INT ) return Helper::int2str( person->meta.get1_int( p ) );
      else if ( mt == META_FLOAT ) return Helper::dbl2str( person->meta.get1_double( p ) );
      else if ( mt == META_BOOL ) return person->meta.get1_bool( p ) ? "T" : "F" ;
      else if ( mt == META_TEXT ) return person->meta.get1_string( p );
      else return ".";
    }
  return ".";
}

Data::Matrix<double> PhenotypeMap::covariates( const std::vector<std::string> & c , const IndividualMap & indmap , 
					       VarDBase * vardb )
{

  // Create a matrix of covariate values
  // The order of rows of this matrix corresponds to the indmap given
  
  // To add -- function to automatically downcode factors?
  // Return a matrix of covariate values

  // If w cannot find a covariate by name, see if it exists as a SNP/region in the VARDB, if a VARDB is passed in

  const int n = indmap.size();

  Data::Matrix<double> d( n , c.size() );

  for (int p=0; p<c.size(); p++)
    {	  
      
      // Does this variable exist in the INDDB at all?
      
      bool covariate_exists = phenotype_exists( c[p] );
      
      if ( !covariate_exists )
	{

	  if ( ! vardb ) 
	    Helper::halt("could not find covariate " + c[p] );

	  // If in form rs12345:AC
	  //  means use AC tag from rs12345
	  //  if so, determine type of tag

	  std::vector<std::string> t = Helper::parse( c[p] , ":" );
	  std::string gettag = "";
	  if ( t.size() == 2 ) gettag = t[1];
	  else if ( t.size() == 3 ) 
	    {
	      t[0] = t[0] + ":" + t[1];  // assume this was in form chr1:12345:TAG
	      gettag = t[2];
	    }

	  bool usetag = gettag != "";
	  
	  const std::string & cp = t[0];
	  
	  bool is_float = false , is_int = false;
	  if ( usetag ) 
	    {
	      mType mt = MetaInformation<GenMeta>::type( gettag );
	      if      ( mt == META_INT ) is_int = true;
	      else if ( mt == META_FLOAT ) is_float = true;
	      else Helper::halt("could not find valid (float/int) genotype tag " + gettag );
	    }
	  
	  // search for in the VARDB: first by rsID name, then by region
	  bool is_region = false;
	  Region region( cp , is_region );
	  
	  // is it instead an ID
	  if ( ! is_region ) 
	    {
	      region = vardb->get_position_from_id( cp , cp );
	      is_region = region.start.chromosome() != 0;
	    }
	  
	  if ( ! is_region ) 
	    Helper::halt("could not find a single covariate value to attach for " + cp ); 
	  
	  std::set<Variant> s = vardb->fetch( region );
	  
	  if ( s.size() != 1 ) 
	    Helper::halt("could not find a single covariate value to attach for for " + cp ); 
	  
	  const Variant & var = *s.begin();

	  // extract tag, or genotype value
	  for (int r=0;r<n;r++)
	    {
	      if ( usetag )
		{
		  if ( var(r).meta.has_field( gettag ) ) 
		    {
		      // has to be one of these two, we checked above
		      if      ( is_float ) d(r,p) = var(r).meta.get1_double( gettag );
		      else if ( is_int ) d(r,p) = var(r).meta.get1_int( gettag );		      		      
		    }
		  else
		    d.set_row_mask( r );
		}
	      else
		{
		  // use standard genotype (alt-allele count) as covariate value
		  if ( var(r).null() ) d.set_row_mask( r ) ;
		  else d(r,p) = var(r).minor_allele_count( true ) ; // assume REF is MAJOR allele here...
		}
	    } 
	}
      
      // If we are here, it must be a normal covariate that exists in the INDDB
      
      if ( covariate_exists ) 
	{
	  for (int r=0; r<n; r++)
	    {
	      
	      Individual * person = indmap(r);
	      
	      if ( person->meta.has_field( c[p] ) )
		{
		  mType mt = MetaInformation<IndivMeta>::type( c[p] );
		  if ( mt == META_INT ) d(r,p) = person->meta.get1_int( c[p] );
		  else if ( mt == META_FLOAT ) d(r,p) = person->meta.get1_double( c[p] );
		  else if ( mt == META_BOOL ) d(r,p) = person->meta.get1_bool( c[p] );
		  else d.set_row_mask( r );	      
		}
	      else // for now, require completely non-missing data
		d.set_row_mask( r );
	    }
	}

    } // next covariate

  return d;
}


void PhenotypeMap::direct_load( const std::string & filename , const std::string & label )
{
  
  Helper::fileExists( filename );
  InFile f( filename );
    
  int expected_col_count = -1;
  
  // Details for the single phenotype to load
  mType mt = META_UNDEFINED;
  int to_load = 0;
  std::string mis_code = ".";
  
  while ( ! f.eof() )
    {
      
      std::string s = f. readLine();
      if ( s == "" ) continue;
      
      // Meta-information? 
      
      if ( s.size() > 2 && s.substr(0,2) == "##" )
	{
	  
	  std::vector<std::string> tok = Helper::quoted_parse( s.substr(2) );
	  
	  if ( tok.size() != 4 ) continue;
	  std::string name = tok[0];	  
	  if ( name != label ) continue;

	  std::string type = tok[1];
	  std::string miss = tok[2];
	  std::string desc = tok[3];
	  
	  mis_code = tok[2];
	  
	  if ( Helper::is_int( type ) ) 
	    {
	      MetaInformation<IndivMeta>::field( name , META_INT , 1 , desc );
	      mt = META_INT;
	      phenotype_type = PHE_DICHOT;
	      phenotype_name = label;
	    }
	  else if ( Helper::is_float( type ) ) 
	    {
	      MetaInformation<IndivMeta>::field( name , META_FLOAT , 1 , desc );
	      mt = META_FLOAT;
	      phenotype_type = PHE_QT;
	      phenotype_name = label;
	    }
	  else 
	    {
	      MetaInformation<IndivMeta>::field( name , META_TEXT , 1 , desc );
	      mt = META_TEXT;
	      phenotype_type = PHE_FACTOR;
	      phenotype_name = label;
	    }
	  
	}

      // Or header line?      
      else if ( s.substr(0,1) == "#" )
	{
	  
	  if ( mt == META_UNDEFINED ) 
	    Helper::halt( "phenotype " + label + " not defined in header of " + filename );

	  // #ID phe1 phe2 phe3	  
	  std::vector<std::string> tok = Helper::parse( s , " \t");	  
	  if ( tok.size() < 2 ) { plog.warn("malformed phenotype file"); continue; } 
	  if ( tok[0] != "#ID" ) { plog.warn("malformed phenotype file"); continue; } 
	  
	  for ( int i = 1 ; i < tok.size(); i++ )
	    {
	      if ( tok[i] != label ) continue;
	      to_load = i;
	    }

	  if ( to_load == 0 ) 
	    Helper::halt( "could not find phenotype " + label + " in " + filename );
	  
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
	  	  
	  
	  Individual * person = ind( tok[0] );
	  
	  
	  //
	  // If person not in indmap/VCF header, skip
	  //

	  if ( ! person ) continue;
	  

	  //
	  // Insert actual phenotypes
	  //

	  // skip 
	  if ( to_load == 0 ) Helper::halt( "phenotype " + label + " not found" );

	  // skip missing values
	  if ( tok[to_load] == mis_code ) continue;
	  
	  // skip invalid values for numerics (as MT will be registered)  
	  if ( mt == META_INT )
	    {    
	      int x;
	      if ( Helper::str2int( tok[to_load] , x ) )
		person->meta.set( label , x );
	    }
	  else if ( mt == META_FLOAT )
	    {
	      double x;
	      if ( Helper::str2dbl( tok[to_load] , x ) )
		person->meta.set( label , x );		    
	    }
	  else 
	    {
	      person->meta.set( label , tok[to_load]  );		    		  		  
	    }
	  
	}
    }
  
  f.close();
  
  set_phenotype( label );

}
