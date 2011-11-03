
#include "filemap.h"
#include "helper.h"
#include "defs.h"
#include "options.h"
#include "bcf.h"

#include <iostream>
#include <fstream>

using namespace std;

using namespace Helper;

extern Log plog;


std::map<string,fType> FileMap::fTypeMap;

void FileMap::setTypes()
{

  // project file 
  fTypeMap["FIDX"] = FIDX;
     
  // Some special project folders
  fTypeMap["OUTPUT"] = OUTPUT;
  fTypeMap["RESOURCES"] = RESOURCES;

  // Some special system files 
  fTypeMap["LOG"] = LOG;
  fTypeMap["TEMP"] = TEMP;
  fTypeMap["METAMETA"] = METAMETA;
  
  // Core databases
  fTypeMap["INDDB"] = INDDB;
  fTypeMap["VARDB"] = VARDB;
  fTypeMap["SEGDB"] = SEGDB;
  fTypeMap["LOCDB"] = LOCDB;
  fTypeMap["REFDB"] = REFDB;
  fTypeMap["SEQDB"] = SEQDB;
  fTypeMap["NETDB"] = NETDB;
  fTypeMap["WGTDB"] = WGTDB;

  // main input files
  fTypeMap["VCF"] = VCF;
  fTypeMap["BCF"] = BCF_FILE;
  fTypeMap["GTF"] = GTF;
  
  fTypeMap["PHE"] = PHE;
  fTypeMap["IND"] = IND;
  
  // misc. utility 'types'
  fTypeMap["INVALID"] = INVALID;
  fTypeMap["UNKNOWN"] = UNKNOWN;
  
}


std::string File::typeName() const
{
  std::map<std::string,fType>::iterator i = FileMap::fTypeMap.begin();
  while ( i != FileMap::fTypeMap.end() )
    {
      if ( i->second == ft )
	return i->first;
      ++i;
    }
  return "INVALID";
}


std::string FileMap::typeName(const fType & t)
{
  std::map<std::string,fType>::iterator i = FileMap::fTypeMap.begin();
  while ( i != FileMap::fTypeMap.end() )
    {
      if ( i->second == t )
	return i->first;
      ++i;
    }
  return "INVALID";
}

std::string InFile::readLine()
{
  std::string sline;
  gzifstream & gg = *this;
  std::getline( gg ,  sline );
  return sline;
}

std::vector< std::string > InFile::tokenizeLine( const std::string & delim )
{
  gzifstream & gg = *this;
  std::vector<std::string> tokens;
  std::string sline;
  std::getline( gg ,  sline );
  return Helper::parse( sline , delim );
}


void FileMap::setCoreFiles( const std::string & f )
{
  
  //
  // Read a user-specified main file index  
  //

  reset();
  
  checkFileExists( f );
  
  addSpecial( FIDX , f );
  
  InFile fidx( f.c_str(), std::ios::in );
  
  while ( ! fidx.eof() )
    {
      
      std::vector<std::string> names = fidx.tokenizeLine( );
      
      if ( names.size() == 0 ) continue;
      
      if ( parse_for_variable( names[0] ) ) continue;
      else names[0] = replace_variable( names[0] );

      if ( names.size() < 2 ) 
	{
	  plog.warn("invalid row in project file (should be tab-delimited)",
		    Helper::print( names ) );
	  continue;

	}
      fType ft = FileMap::type( names[1] );
      
      if ( ft == INVALID ) continue;
      
      if ( ft == OUTPUT ) 
	{
	  // ensure specifies a directory w/ trailing "/"
	  std::string folder = names[0];

	  if ( folder.substr( folder.size()-1,1 ) != "/" )
	    folder += "/";
	  addSpecial( OUTPUT , folder );
	}
      else if ( ft == RESOURCES )
	{
	  // ensure specifies a directory w/ trailing "/"
	  std::string folder = names[0];
	  if ( folder.substr( folder.size()-1,1 ) != "/" )
	    folder += "/";
	  addSpecial( RESOURCES , folder );
	}
      else if ( ft == METAMETA )
	{
	  MetaMeta::load( names[0] );
	}
      else if ( ft == TEMP )
	{	    
	  PLINKSeq::SQLITE_SCRATCH_FOLDER() = names[0];
	}
      else if ( ft == LOG ) addSpecial( LOG , names[0] );
      else if ( ft == VARDB ) addSpecial( VARDB , names[0] );
      else if ( ft == INDDB ) addSpecial( INDDB , names[0] );
      else if ( ft == SEGDB ) addSpecial( SEGDB , names[0] );
      else if ( ft == LOCDB ) addSpecial( LOCDB , names[0] );
      else if ( ft == NETDB ) addSpecial( NETDB , names[0] );
      else if ( ft == WGTDB ) addSpecial( WGTDB , names[0] );
      else if ( ft == REFDB ) addSpecial( REFDB , names[0] );
      else if ( ft == SEQDB ) addSpecial( SEQDB , names[0] );
      else if ( ft == BCF_FILE )   add_BCF( names[0] );
    }    
  
  fidx.close();
  
  
  ///////////////////////////////////////////
  // Check that some key defaults are given
  
  if ( ! file( OUTPUT ) ) 
    addSpecial( OUTPUT , f + "_out/" );
  
  if ( ! file( RESOURCES ) )
    addSpecial( RESOURCES , f + "_res/" );
  
  if ( ! file( LOG ) )
    addSpecial( LOG , file( OUTPUT )->name() + "log.txt" );
  
  //
  // project-specific data
  //
  
  if ( ! file( VARDB ) )
    addSpecial( VARDB , file( OUTPUT )->name() + "vardb" );
  
  if ( ! file( INDDB ) )
    addSpecial( INDDB , file( OUTPUT )->name() + "inddb" );
  
  
  //
  // (shared) resources
  //
  
  if ( ! file( LOCDB ) )
    addSpecial( LOCDB , file( RESOURCES )->name() + "locdb" );
  
  if ( ! file( REFDB ) )
    addSpecial( REFDB , file( RESOURCES )->name() + "refdb" );
  
  if ( ! file( SEQDB ) )
    addSpecial( SEQDB , file( RESOURCES )->name() + "seqdb" );
  
  
    
  //
  // Create main folders
  //
  
  make_dir( file( OUTPUT )->name() );
  
  make_dir( file( RESOURCES )->name() );
  
}


bool FileMap::readFileIndex( const std::string & f )
{

  // FORMAT: uncompressed, plain text
  //         1 line per file; 2 fields
  //         fullpath/name    filetype    comments
  
  InFile fidx( f , std::ios::in );

  
  while ( ! fidx.eof() ) 
    {
      
      std::vector<std::string> names = fidx.tokenizeLine( );

      if ( names.size() == 0 ) continue;      
      
      std::string filename = names[0];
      
      if ( parse_for_variable( filename ) )
	  continue;
      else filename = replace_variable( filename );
      
      // Ignore badly formated line
      
      if ( names.size() < 2 )
	continue;
      

      
      //
      // Check a legal type has been specified
      //
      
      fType ft = FileMap::type( names[1] );

      if ( ft == INVALID ) continue;

      //
      // Replace variables
      //

      // 0          1       2      3,4,...
      // file-name  TYPE    TAG    COMMENT... 
      
      std::string file_tag;
      
      if ( names.size() > 2 ) 
	file_tag = names[2];

      // Compile any comment
      
      std::string comment = "";
      
      for (unsigned int i=3; i<names.size(); i++) 
	comment += " " + names[i]; 
      
      // Only insert non-core file
      
      if ( special_files.find( ft ) == special_files.end() )
	add( names[0] , ft , file_tag , comment );
      
    }  
  
  fidx.close();
  
  return true;
  
}


void FileMap::reset()
{ 
  std::map<std::string,File*>::iterator i = fmap.begin();
  while ( i != fmap.end() )
    {
      if ( i->second )
	delete i->second;
      ++i;
    }
  fmap.clear();
  special_files.clear();
}


File * FileMap::file( const fType & t ) const
{
  std::map<fType,File*>::const_iterator i = special_files.find(t);
  return i != special_files.end() ? i->second : NULL ;
}


bool FileMap::exists( const fType & t ) const
{
  std::map<fType, File*>::const_iterator i = special_files.find(t);
  return i != special_files.end();     
}

void FileMap::addSpecial(fType t, std::string n)
{
  special_files[ t ] = new File(n,t);
}


bool FileMap::exists( const std::string & f ) const
{
  std::map<std::string,File*>::const_iterator i = fmap.find( f );
  return i != fmap.end();
}


std::string FileMap::tilde_expansion( const std::string & f )
{
  wordexp_t exp_result;
  wordexp( f.c_str() , &exp_result, 0);
  std::string nf = exp_result.we_wordv[0];
  wordfree(&exp_result);
  return nf;
}

File * FileMap::add( const std::string & n, 
		     fType t, 
		     const std::string & tag ,
		     const std::string & comment )
{
  
  // return pointer 
  std::map< std::string, File*>::iterator i = fmap.find( n );
  if ( i != fmap.end() ) return fmap.find(n)->second;
  
  File * f = new File(n,t);  
  f->included( fileExists(n) );  
  f->comment( comment );
  f->tag( tag );
  
  // add to map
  fmap.insert(make_pair( f->name() , f )) ;
  
  return f;
}



std::string FileMap::summary() const
{
  std::stringstream ss;

  ss << "FILE_INDEX" << "\t"
     << "TYPE=INDEX" << "\t"
     << "NAME=" << special_files.find( FIDX )->second->name() << "\n";
  
  std::map< std::string, File* >::const_iterator f = fmap.begin();
  while ( f != fmap.end() )
    {
      ss << "FILE_INDEX" << "\t"
	 << "TYPE=" << FileMap::typeName( f->second->type() ) << "\t"
	 << "NAME=" << f->second->name() << "\n";
      ++f;
    }
  
  std::map<fType,File*>::const_iterator i = special_files.begin();
  while ( i != special_files.end() )
    {
      if ( i->first != FIDX ) // already listed this first
	ss << "FILE_INDEX" << "\t" 
	   << "TYPE=" << FileMap::typeName( i->first ) << "\t"
	   << "NAME=" << i->second->name() << "\n";
      ++i;
    }

  return ss.str();
}


File * FileMap::firstFile() 
{       
  if ( fmap.size() == 0 ) 
    {
      currf = fmap.end();
      return NULL;
    }
  else
    {
      currf = fmap.begin();
      return currf->second;
    }
}

File * FileMap::nextFile() 
{
  ++currf;
  return currf == fmap.end() ? NULL : currf->second ;
}


std::set<File*> FileMap::get( const fType t ) const
{
  std::set<File*> s;
  std::map< std::string, File* >::const_iterator f = fmap.begin();
  while ( f != fmap.end() )
    {
      if ( f->second->included() && f->second->type() == t )
	s.insert( f->second );
      ++f;
    }
  return s;  
}


File * FileMap::file( const std::string & t ) const
{
  std::map<std::string,File*>::const_iterator i = fmap.find(t);
  return i != fmap.end() ? i->second : NULL ;
}


bool FileMap::parse_for_variable( const std::string & t )
{
  if ( t.substr(0,1) == "#" ) 
    {
      std::string s = t.substr(1) ;
      // in format  #var=val
      
      if ( s.find("=") != std::string::npos ) 
	{
	  vmap[ "${" + s.substr( 0 , s.find("=") ) + "}" ] 
	    = s.substr( s.find("=") + 1 );
	}
      return true;
    }
  return false;
}


std::string FileMap::replace_variable( std::string & s )
{
  std::map<std::string,std::string>::iterator i = vmap.begin();
  while ( i != vmap.end() )
    {
      if ( s.find( i->first ) != std::string::npos )
	{
	  size_t sz = s.find( i->first );	    
	  s.replace( sz , sz + (i->first).size() , i->second );	    
	}
      ++i;
    }
  return s;
}


BCF * FileMap::bcf( const std::string & filename )
{
  return bcf_map[ filename ]; // NULL if not in map
}


BCF * FileMap::add_BCF( const std::string & f )
{
  BCF * bcf = new BCF( f );
  if ( bcf ) 
    {
      bcf_map[ f ] = bcf;
      // also add to normal filemap
      add( f , BCF_FILE , "" , "BCF" );
    }
  return bcf;
}


bool FileMap::append_to_projectfile( const std::string & s , const std::string & t )
{
  
  if ( exists( s ) ) return false; // already present, nothing to do
  
  std::string projectfile = special_files.find( FIDX )->second->name();
  
  if ( projectfile == "." ) return false;
  
  if ( ! Helper::fileExists( projectfile ) )
    {
      plog.warn("could not find projectfile",projectfile);
      return false;
    }
  
  // open in append-to-end mode
  std::ofstream O1( projectfile.c_str() , std::ios::out | std::ios::app );
  O1 << s << "\t" << t << "\n";
  O1.close();
  
  // add to internal map
  add( s , type(t) , "" , "" );
}


bool FileMap::remove_from_projectfile( const std::string & s )
{
  std::string projectfile = special_files.find( FIDX )->second->name();
  if ( projectfile == "." ) return false;
  if ( ! Helper::fileExists( projectfile ) )
    {
      plog.warn("could not find projectfile",projectfile);
      return false;
    }

  // open in append-to-end mode
  InFile O1( projectfile );
  std::vector<std::string> lines;
  while ( ! O1.eof() )
    {
      std::string l = O1.readLine();
      if ( l == "" ) continue;
      std::vector<std::string> tok = Helper::char_split(l,'\t');
      if ( ! ( tok[0] == s || ( tok.size()>1 && tok[1] == s ) )  ) lines.push_back(l);
    }
  O1.close();
  
  // write back out, 
  std::ofstream O2( projectfile.c_str() );
  for (int l=0;l<lines.size();l++)
    O2 << lines[l] << "\n";    
  O2.close();

}

bool FileMap::write_new_projectfile( const std::string & filename ) 
{
  
  std::ofstream O2( filename.c_str() );
  
  std::map< std::string, File* >::const_iterator f = fmap.begin();
  while ( f != fmap.end() )
    {
      O2 << f->second->name() << "\t" 
	 << FileMap::typeName( f->second->type() ) << "\n";
      ++f;
    }
  
  std::map<fType,File*>::const_iterator i = special_files.begin();
  while ( i != special_files.end() )
    {
      if ( i->first != FIDX ) // not needed
	O2 << i->second->name() << "\t"
	   << FileMap::typeName( i->first )  << "\n";
      ++i;
    }
  
  O2.close();
  return true;
}
