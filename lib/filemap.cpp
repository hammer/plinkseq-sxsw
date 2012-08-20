
#include "plinkseq/filemap.h"
#include "plinkseq/helper.h"
#include "plinkseq/defs.h"
#include "plinkseq/bcf.h"
#include "plinkseq/vcfz.h"
#include "plinkseq/gstore.h"

#include <iostream>
#include <fstream>

using namespace std;

using namespace Helper;

extern Log plog;
extern GStore * GP;

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
  fTypeMap["VCFZ"] = BGZF_VCF;
  fTypeMap["GTF"] = GTF;
  
  fTypeMap["PHE"] = PHE;
  fTypeMap["IND"] = IND;
  
  // misc. utility 'types'
  fTypeMap["INVALID"] = INVALID;
  fTypeMap["UNKNOWN"] = UNKNOWN;

  // project password
  fTypeMap["PASSWD"] = PWD;

  fTypeMap["PARAM"] = PARAM;

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
  
  Helper::checkFileExists( f );
  
  addSpecial( FIDX , f );
  
  InFile fidx( f.c_str(), std::ios::in );

  // Old, default mode is in format   {file-name}  {type}   {description}

  // New (proj 1+) is                 {key}       {value}   {ignored comments}
  // and has first line               PROJN       1  
  
  int project_file_mode = 0;

  bool firstline = true;

  while ( ! fidx.eof() )
    {
      
      std::vector<std::string> names = fidx.tokenizeLine( );
            
      if ( names.size() == 0 ) continue;
      
      if ( firstline ) 
	{
	  if ( names[0] == "PROJN" )
	    {
	      if ( names.size() < 2 ) Helper::halt("invalid PROJN line in project file");
	      if ( ! Helper::str2int( names[1] , project_file_mode ) )
		Helper::halt("invalid PROJN line in project file");
	    }
	  firstline = false;
	}

      // Swap names around from old mode
      if ( project_file_mode == 0 )
	{
	  if ( names.size() > 1 ) 
	    {
	      std::string tmp = names[0];
	      names[0] = names[1];
	      names[1] = tmp;
	    }	  
	}

      std::string & value = names[1];
      std::string & key   = names[0];

      if ( parse_for_variable( value ) ) continue;
      else names[1] = replace_variable( value );
      
      if ( names.size() < 2 ) 
	{
	  plog.warn("invalid row in project file (should be tab-delimited)",
		    Helper::print( names ) );
	  continue;

	}
      fType ft = FileMap::type( key );
      
      if ( ft == INVALID ) continue;
      
      if ( ft == OUTPUT ) 
	{
	  // ensure specifies a directory w/ trailing "/"
	  std::string folder = value;
	  
	  if ( folder.substr( folder.size()-1,1 ) != "/" )
	    folder += "/";
	  addSpecial( OUTPUT , folder );
	}
      else if ( ft == RESOURCES )
	{
	  // ensure specifies a directory w/ trailing "/"
	  std::string folder = value;
	  if ( folder.substr( folder.size()-1,1 ) != "/" )
	    folder += "/";
	  addSpecial( RESOURCES , folder );
	}
      else if ( ft == METAMETA )
	{
	  MetaMeta::load( value );
	}
      else if ( ft == TEMP )
	{	    
	  PLINKSeq::SQLITE_SCRATCH_FOLDER() = value;
	}
      else if ( ft == LOG ) addSpecial( LOG , value );
      else if ( ft == VARDB ) addSpecial( VARDB , value );
      else if ( ft == INDDB ) addSpecial( INDDB , value );
      else if ( ft == SEGDB ) addSpecial( SEGDB , value );
      else if ( ft == LOCDB ) addSpecial( LOCDB , value );
      else if ( ft == NETDB ) addSpecial( NETDB , value );
      else if ( ft == WGTDB ) addSpecial( WGTDB , value );
      else if ( ft == REFDB ) addSpecial( REFDB , value );
      else if ( ft == SEQDB ) addSpecial( SEQDB , value );
      else if ( ft == BCF_FILE )   add_BCF( value );
      else if ( ft == BGZF_VCF )   add_VCFZ( value );

      // non-file operations
      else if ( ft == PWD ) GP->set_pwd( value ); 
      else if ( ft == PARAM ) GP->set_param( value );
      
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

  int project_file_mode = 0;
  bool firstline = true;
  
  while ( ! fidx.eof() ) 
    {
      
      std::vector<std::string> names = fidx.tokenizeLine( );

      if ( names.size() == 0 ) continue;      
      
      if ( firstline ) 
	{
	  if ( names[0] == "PROJN" )
	    {
	      if ( names.size() < 2 ) Helper::halt("invalid PROJN line in project file");
	      if ( ! Helper::str2int( names[1] , project_file_mode ) )
		Helper::halt("invalid PROJN line in project file");
	    }
	  firstline = false;
	}

      // Swap names around from old mode
      if ( project_file_mode == 0 )
	{
	  if ( names.size() > 1 ) 
	    {
	      std::string tmp = names[0];
	      names[0] = names[1];
	      names[1] = tmp;
	    }	  
	}

      std::string & value = names[1];
      std::string & key   = names[0];

      
      std::string filename = value;
      
      if ( parse_for_variable( filename ) )
	  continue;
      else filename = replace_variable( filename );
      
      // Ignore badly formated line
      
      if ( names.size() < 2 )
	continue;
      

      
      //
      // Check a legal type has been specified
      //
      
      fType ft = FileMap::type( key );

      if ( ft == INVALID ) continue;

      //
      // Replace variables
      //

      std::string file_tag;
      
      if ( names.size() > 2 ) file_tag = names[2];
      
      // Compile any comment
      
      std::string comment = "";
      
      for (unsigned int i=3; i<names.size(); i++) 
	comment += " " + names[i]; 
      
      // Only insert non-core file
      
      if ( special_files.find( ft ) == special_files.end() )
	add( value , ft , file_tag , comment );
      
    }  
  
  fidx.close();
  
  return true;
  
}


void FileMap::reset()
{ 
    std::map<std::string,File*>::iterator i = fmap.begin();
    while ( i != fmap.end() )
    {
	if ( i->second ) delete i->second;
	i->second = NULL;
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



std::string FileMap::summary( bool ugly ) const
{

  std::stringstream ss;

  if ( ugly ) 
    {
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
    }
  else
    {
      ss << "---File-index summary---\n\n";

      ss << "Core project specification index : " << special_files.find( FIDX )->second->name()  << "\n"; 
      
      std::map<fType,File*>::const_iterator i = special_files.begin();
      while ( i != special_files.end() )
	{
	  if ( i->first != FIDX ) // already listed this first
	    ss << "Core " << FileMap::typeName( i->second->type() ) << " file : " << i->second->name() << "\n";
	  ++i;
	}
      
      std::map< std::string, File* >::const_iterator f = fmap.begin();
      while ( f != fmap.end() )
	{
	  ss << "Added " << FileMap::typeName( f->second->type() ) << " : " << f->second->name() << "\n";
	  ++f;
	}
        
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


VCFZ * FileMap::vcfz( const std::string & filename )
{
    return vcfz_map[ filename ]; // NULL if not in map
}

VCFZ * FileMap::add_VCFZ( const std::string & f )
{
    
  VCFZ * vcfz = new VCFZ( f );
  
  if ( vcfz ) 
    {
      vcfz_map[ f ] = vcfz;
      
      // also add to normal filemap
      add( f , BGZF_VCF , "" , "VCFZ" );
    }
  
  return vcfz;
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
  O1 << t << "\t" << s << "\n";
  O1.close();
  
  // add to internal map
  add( s , type(t) , "" , "" );
}


bool FileMap::remove_from_projectfile( const std::string & s )
{
  // this function okay for either v1 or v2 format (value/key vs key/value)
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

  O2 << "PROJN\t2\n";

  std::map< std::string, File* >::const_iterator f = fmap.begin();
  while ( f != fmap.end() )
    {
      O2 << FileMap::typeName( f->second->type() ) << "\t" 
	 << f->second->name() << "\n";
      ++f;
    }
  
  std::map<fType,File*>::const_iterator i = special_files.begin();
  while ( i != special_files.end() )
    {
      if ( i->first != FIDX ) // not needed
	O2 << FileMap::typeName( i->first ) << "\t"
	   << i->second->name() << "\n";
      ++i;
    }
  
  O2.close();
  return true;
}
