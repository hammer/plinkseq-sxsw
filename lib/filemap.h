#ifndef __FILEMAP_H__
#define __FILEMAP_H__

#include <iostream>
#include <iomanip>
#include <fstream>
#include <vector>
#include <set>


#ifdef __WIN
#include <direct>
#else
#include <sys/stat.h>
#endif

#include <sys/stat.h>
#include <sys/types.h>
#include <cerrno>

#include <cstdio>
#ifdef WINDOWS
#include <direct.h>
#define GetCurrentDir _getcwd
#else
#include <unistd.h>
#define GetCurrentDir getcwd
#endif

#include "helper.h"
#include "defs.h"
#include "meta.h"
#include "zfstream.h"
#include "bcf.h"

extern Log plog;


// wrappers for input and output streams, allowing future 
// addition on content at this level

class InFile : public gzifstream {
 public:

  InFile() { } 
  
 InFile( const std::string & n , 
	 ios_base::openmode mode = ios_base::in ) 
   : gzifstream(n.c_str(), mode) { Helper::checkFileExists(n); } 
  
  std::string readLine();
  std::vector< std::string > tokenizeLine( const std::string & delim = PLINKSeq::DELIM() );
};

class OutFile : public gzofstream {
 public:
 OutFile( const std::string & n , 
	  ios_base::openmode mode = ios_base::out ) 
   : gzofstream(n.c_str(), mode) { }   
};


// class that represents some of the key information
// with a "file" (that might or might not yet exist)

class File {

 private:

  std::string fname;
  std::string ftag;
  std::string vnote;
  bool include;
  fType ft;

 public:
  
 File( const std::string & n , const fType & t )
   : fname(n) , ft(t) 
   { 
    include = true;
    vnote = "";
    ftag = "";
  } 
 
  MetaInformation<FileMeta>   meta;

  void name(const std::string & n) { fname = n; }
  std::string name() const { return fname; }

  std::string tag() const { return ftag; }
  void tag(std::string s)  { ftag=s; }
 
  std::string comment() const { return vnote; }
  void comment(std::string s)  { vnote=s; }
  
  bool included() const { return include; }
  void included(bool i) { include = i; }
  
  void type(const fType & t) { ft = t; }
  fType type() const { return ft; }
  
  std::string typeName() const;

};


// organise the Files into one place, i.e. that corresponds to the
// project-specification file

class FileMap {

 private:
    
    bool parse_for_variable( const std::string & );
    std::string replace_variable( std::string & );
    std::map<std::string,std::string> vmap;

 public:

  //
  // File-type information and functions
  //

  static void setTypes();
  
  static std::map<std::string,fType> fTypeMap;
  
  
  
  static std::string typeName(const fType&);
  
  static fType type(std::string n)
  {
    std::map<std::string,fType>::iterator i = FileMap::fTypeMap.find(n);
    if ( i == FileMap::fTypeMap.end() ) 
      return INVALID;
    else 
      return i->second;
  }

  
    
  // Main file MAP
  
  std::map<std::string,File*> fmap;

  // Special instance for tracking BCF files
  std::map<std::string,BCF*> bcf_map;
  
  // Special/core files, of which we expect exactly one
  //  e.g. VARDB
  
  std::map<fType,File*> special_files;  
  

  
  //
  // Current file
  //
  
  std::map< std::string, File*>::iterator currf;
  

  //
  // Some key values, that need defaults
  //




 public:


  FileMap() { }  
  
  ~FileMap()
    {
      std::map<std::string,BCF*>::iterator i = bcf_map.begin();
      while ( i != bcf_map.end() )
	{
	  if ( i->second ) delete i->second;
	  ++i;
	}
    }

  void reset();



  //
  // Add to index
  //

  File * add( const std::string & , 
	      fType, 
	      const std::string & , 
	      const std::string & );
  
  void addSpecial( fType, std::string );

  
  //
  // Populate index, from file or database
  //
  
  void setCoreFiles( const std::string & );

  bool readFileIndex( const std::string & );

  //
  // Project file direct manipulation
  //
  
  bool append_to_projectfile( const std::string & , const std::string & );
  bool remove_from_projectfile( const std::string & );
  bool write_new_projectfile( const std::string & );

  //
  // Lookup information on files
  //

  bool exists( const fType & ) const;
  bool exists( const std::string & ) const;
  
  File * file( const fType & ) const;
  File * file( const std::string &  ) const;
   
  std::string summary() const;

  //
  // Specific BCF indexing
  //
  
  BCF * bcf( const std::string & );

  BCF * add_BCF( const std::string & f );

  //
  // Group access functions ( i.e. get all valid VCF )
  //

  std::set<File*> get( const fType ) const;
  
  
  // Iterators over whole fileset

  File * firstFile();
  
  File * nextFile();


  //
  // Misc.
  //
  
  bool make_dir( const std::string & s ) 
  { 	
    
    if ( mkdir( s.c_str() , 0777) == -1 )
      {
	if ( errno == EEXIST ) return true;
	plog.warn("could not create directory: " + s ); 
	return false;
      }	   
    return true;
  }


  static std::string working_directory()
    {      
      char cCurrentPath[ FILENAME_MAX ];      
      if ( ! GetCurrentDir(cCurrentPath, sizeof(cCurrentPath) ) )
	Helper::halt("problem getting current working directory in FileMap()");      
      cCurrentPath[sizeof(cCurrentPath) - 1] = '\0'; // not really required
      return cCurrentPath;
    }
  
  
};



#endif
