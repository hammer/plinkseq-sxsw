#include <iostream>
#include <sstream>
#include <fstream>
#include <iomanip>

#include <cstdio>
#include <cerrno>

#include "helper.h"
#include "regions.h"
#include "variant.h"
#include "gstore.h"

#ifdef R_SHLIB
#include "rint.h"
#endif

extern GStore * GP;

using namespace std;
using namespace Helper;

//////////////////////////////////////////////////////////////////////
// Logging and error handling


void Helper::halt( const std::string & msg )
{
#ifdef R_SHLIB
 R_error( msg );
#else
 std::cerr << "pseq error : " << msg << "\n";
 if ( GP && GP->gseq_mode() )
   {
     std::ofstream O1( GP->gseq_history().c_str() , std::ios::out | std::ios::app );
     O1 << "_STATUS" << "\t"
	<< GP->gseq_job() << "\t"
	<< "error: "<< msg  << "\n";
     O1.close();
   } 
 exit(1);
#endif
}

void Log::warn(const std::string & msg, const std::vector<std::string> & spec )
{
  warn( msg , Helper::stringize( spec , " " ) );
}

void Log::warn(const std::string & msg , const std::string & spec ) 
{

  if ( ignore_warnings ) return;

  if ( ! silent_mode )
    {
#ifdef R_SHLIB
  R_warning( msg );
#else
  if ( warnings[ msg ] == 0 && early_warn )
    {
      std::cerr << "plinkseq warning: " << msg << " : " << spec << "\n";
      std::cerr.flush();      
    }
#endif
    }

  // keep track of how many times each warning issued
  warnings[ msg ]++;
  if ( spec != "" && warnings[ msg ] < 10 ) 
    warnings_specific[ msg ].push_back( spec );

}




void Log::print_warnings()
{
  std::map<std::string,int>::const_iterator i = warnings.begin();
  while ( i != warnings.end() )
    {

      std::vector<std::string> & arr = warnings_specific[ i->first ];            
      
      std::string msg = "";

      if ( arr.size() == 0 ) 
	{
	  if ( i->second > 1 ) 
	    msg += "plinkseq warning: " + i->first + " (repeated " + Helper::int2str( i->second ) + " times)\n";
	  else
	    msg += "plinkseq warning: " + i->first + "\n";
	}
      else
	{
	  for (int j=0; j<arr.size(); j++)
	    msg += "plinkseq warning: " + i->first + " : " + arr[j] + "\n";
	  if ( i->second > arr.size() )
	    msg += "plinkseq warning: " + i->first + " (repeated " + Helper::int2str( i->second ) + " times)\n";
	}	  
      
#ifdef R_SHLIB
      // R should keep track of warning count, so no need for further action here?
#else
      if ( ! silent_mode ) std::cerr << msg ;
#endif

      if ( output_file ) file << msg ;
      ++i;
    }
}



void Helper::NoMem()
{
  std::cerr << "*****************************************************\n"
	    << "* FATAL ERROR    Exhausted system memory            *\n"
	    << "*****************************************************\n\n";

  if ( GP && GP->gseq_mode() )
    {
      std::ofstream O1( GP->gseq_history().c_str() , std::ios::out | std::ios::app );
      O1 << "_STATUS" << "\t"
	 << GP->gseq_job() << "\t"
	 << "failed: out of memory" << "\n";
      O1.close();
    }
    
  exit(1);
}



//////////////////////////////////////////////////////////////////////
// Conversion functions


std::string Helper::int2str(int n)
{
  std::ostringstream s2( std::stringstream::out );
  s2 << n;
  return s2.str();
}

std::string Helper::uint64_t2str(uint64_t n)
{
  std::ostringstream s2( std::stringstream::out );
  s2 << n;
  return s2.str();
}

std::string Helper::longint2str(long int n)
{
  std::ostringstream s2( std::stringstream::out );
  s2 << n;
  return s2.str();
}

std::string Helper::dbl2str_fixed(double n, int prc)
{
  std::ostringstream s2;
  s2 << setiosflags( ios::fixed );
  if ( prc > 0 )
    s2.precision(prc);
  s2 << n;
  return s2.str();
}

std::string Helper::dbl2str(double n, int prc)
{
  std::ostringstream s2;
  if ( prc > 0 )
    s2.precision(prc);
  s2 << n;
  return s2.str();
}

std::string Helper::flt2str(const float n, const int prc)
{
  std::ostringstream s2;
  if ( prc > 0 )
    s2.precision(prc);
  s2 << n;
  return s2.str();
}

bool Helper::valid_name( const std::string & s )
{
  // a valid name starts with a letter, contains only A-Z, a-z, digits or _
  if ( s.size() == 0 ) return false;
  char c = s[0];
  if ( c >= '0' && c <= '9' ) return false;
  for (int i=0; i<s.size(); i++)
    {
      char c = s[i];
      bool okay = false;
      if ( c >= '0' && c <= '9' ) okay = true;
      else if ( c >= 'a' && c <= 'z' ) okay = true;
      else if ( c >= 'A' && c <= 'Z' ) okay = true;
      else if ( c == '_' ) okay = true;      
      if ( ! okay ) return false;
    }
  return true;
}

// C++ style conversion functions

bool Helper::str2bool( const std::string & s , const std::string & miss )
{
  // returns a positive flag
  if ( s == "0" || s == "F" || s == "f" || s == miss ) return false;
  return true;
}

bool Helper::str2int(const std::string & s , int & i)
{
  return from_string<int>(i,s,std::dec);
}

int Helper::str2int(const std::string & s)
{
    int i = 0;
    from_string<int>(i,s,std::dec);
    return i;
}

bool Helper::str2uint64_t(string s , uint64_t & i)
{
  return from_string<uint64_t>(i,s,std::dec);
}

bool Helper::str2dbl(const std::string & s , double & i)
{
  return from_string<double>(i,s,std::dec);
}

double Helper::str2dbl(const std::string & s)
{
    double d = 0;
    from_string<double>(d,s,std::dec);
    return d;
}

// C-style conversion functions

bool Helper::str2uint64_t(const char * c , uint64_t & i)
{
    
}

bool Helper::str2int(const char * c , int & i)
{
    errno = 0;
    char * endptr;
    i = strtol( c , &endptr , 10 );
    if ( *endptr == '\0' ) return true;
    i = 0; 
    return false;
}

bool Helper::str2dbl(const char * c , double & i)
{
    errno = 0;
    char * endptr;
    i = strtod( c , &endptr );
    if ( *endptr == '\0' ) return true;
    i = 0; 
    return false;
}

int Helper::str2int(const char * c )
{
    char * endptr;
    int i = strtol( c , &endptr , 10 );
    if ( *endptr == '\0' ) return i;
    plog.warn("problem converting string to integer");
    return 0;    
}

double Helper::str2dbl(const char * c)
{
    char * endptr;
    double i = strtod( c , &endptr  );
    if ( *endptr == '\0' ) return i;
    plog.warn("problem converting string to integer");
    return 0;        
}

bool Helper::str2bool( const char * c )
{
    // match on single, first character
    if ( *c == '0' || *c == 'F' || *c == 'f' || *c == '.' ) return false;
    return true;
}



std::string Helper::search_replace( std::string & str , const std::string & search , const std::string & replace )
{
  string::size_type pos = 0;
  while ( (pos = str.find(search, pos)) != string::npos ) 
    {
      str.replace( pos, search.size(), replace );
      pos++;
    }
  return str;
}

std::string Helper::header( const std::string & s , const int len , const std::string & rep ) 
{
  if ( s.size() <= len ) return s;
  return s.substr(0,len-3) + rep + rep + rep;
}


bool Helper::inCommaList( const std::string & lst, const std::string & term)
{
  std::vector<std::string> tok = Helper::char_split( lst , ',' );
  std::vector<std::string>::iterator tok_iter = tok.begin();
  while ( tok_iter != tok.end() )
    {
      if ( *tok_iter == term )
	return true;
      ++tok_iter;
    }
  return false;
}

std::vector<std::string> Helper::char_split( const std::string & s , const char c , bool empty )
{

  std::vector<std::string> strs;  
  if ( s.size() == 0 ) return strs;
  int p=0;

  for (int j=0; j<s.size(); j++)
    {	        
      if ( s[j] == c ) 
	{ 	      
	  if ( j == p ) // empty slot?
	    {
	      if ( empty ) strs.push_back( "." );
	      ++p;
	    }
	  else
	    {
	      strs.push_back(s.substr(p,j-p)); 
	      p=j+1; 
	    }
	}	  
    }
  
  if ( empty && p == s.size() ) 
    strs.push_back( "." );
  else if ( p < s.size() )
    strs.push_back( s.substr(p) );
  
  return strs;
}

std::vector<std::string> Helper::char_split( const std::string & s , const char c , const char c2 , bool empty )
{
  std::vector<std::string> strs;  
  if ( s.size() == 0 ) return strs;
  int p=0;

  for (int j=0; j<s.size(); j++)
    {	        
      if ( s[j] == c || s[j] == c2 ) 
	{ 	      
	  if ( j == p ) // empty slot?
	    {
	      if (empty) strs.push_back( "." );
	      ++p;
	    }
	  else
	    {
	      strs.push_back(s.substr(p,j-p)); 
	      p=j+1; 
	    }
	}	  
    }
  
  if ( empty && p == s.size() ) 
    strs.push_back( "." );
  else if ( p < s.size() )
    strs.push_back( s.substr(p) );
  
  return strs;
}


std::vector<std::string> Helper::char_split( const std::string & s , const char c , const char c2 , const char c3 , bool empty )
{
  std::vector<std::string> strs;  
  if ( s.size() == 0 ) return strs;
  int p=0;

  for (int j=0; j<s.size(); j++)
    {	        
      if ( s[j] == c || s[j] == c2 || s[j] == c3 ) 
	{ 	      
	  if ( j == p ) // empty slot?
	    {
	      if ( empty ) strs.push_back( "." );
	      ++p;
	    }
	  else
	    {
	      strs.push_back(s.substr(p,j-p)); 
	      p=j+1; 
	    }
	}	  
    }
  
  if ( empty && p == s.size() ) 
    strs.push_back( "." );
  else if ( p < s.size() )
    strs.push_back( s.substr(p) );
  
  return strs;
}


std::vector<std::string> Helper::quoted_char_split( const std::string & s , const char c , bool empty )
{

  std::vector<std::string> strs;  
  if ( s.size() == 0 ) return strs;
  int p=0;

  bool in_quote = false;

  for (int j=0; j<s.size(); j++)
    {	        

      if ( s[j] == '"' ) in_quote = ! in_quote;

      if ( (!in_quote) && s[j] == c ) 
	{ 	      
	  if ( j == p ) // empty slot?
	    {
	      if ( empty ) strs.push_back( "." );
	      ++p;
	    }
	  else
	    {
	      strs.push_back(s.substr(p,j-p)); 
	      p=j+1; 
	    }
	}	  
    }
  
  if ( empty && p == s.size() ) 
    strs.push_back( "." );
  else if ( p < s.size() )
    strs.push_back( s.substr(p) );
  
  return strs;
}


std::vector<std::string> Helper::quoted_char_split( const std::string & s , const char c , const char c2 , bool empty )
{

  std::vector<std::string> strs;  
  if ( s.size() == 0 ) return strs;
  int p=0;

  bool in_quote = false;

  for (int j=0; j<s.size(); j++)
    {	        

      if ( s[j] == '"' ) in_quote = ! in_quote;

      if ( (!in_quote) && ( s[j] == c || s[j] == c2 ) ) 
	{ 	      
	  if ( j == p ) // empty slot?
	    {
	      if ( empty ) strs.push_back( "." );
	      ++p;
	    }
	  else
	    {
	      strs.push_back(s.substr(p,j-p)); 
	      p=j+1; 
	    }
	}	  
    }
  
  if ( empty && p == s.size() ) 
    strs.push_back( "." );
  else if ( p < s.size() )
    strs.push_back( s.substr(p) );
  
  return strs;
}


std::vector<std::string> Helper::quoted_char_split( const std::string & s , const char c , const char c2 , const char c3 , bool empty )
{

  std::vector<std::string> strs;  
  if ( s.size() == 0 ) return strs;
  int p=0;

  bool in_quote = false;

  for (int j=0; j<s.size(); j++)
    {	        

      if ( s[j] == '"' ) in_quote = ! in_quote;

      if ( (!in_quote) && ( s[j] == c || s[j] == c2 || s[j] == c3 ) ) 
	{ 	      
	  if ( j == p ) // empty slot?
	    {
	      if (empty) strs.push_back( "." );
	      ++p;
	    }
	  else
	    {
	      strs.push_back(s.substr(p,j-p)); 
	      p=j+1; 
	    }
	}	  
    }
  
  if ( empty && p == s.size() ) 
    strs.push_back( "." );
  else if ( p < s.size() )
    strs.push_back( s.substr(p) );
  
  return strs;
}




std::vector<std::string> Helper::tokenizeLine(std::ifstream & F1)
{
  std::string sline;
  getline(F1,sline);
  std::string buf; 
  std::stringstream ss(sline); 
  std::vector<std::string> tokens; 
  while (ss >> buf)
    tokens.push_back(buf);
  return tokens;
}

std::string Helper::stringize( const std::set<std::string> & s , const std::string & delim  )
{
  std::string r = "";
  std::set<std::string>::iterator i = s.begin();
  while ( i != s.end() )
    {
      if ( i != s.begin() ) r += delim;
      r += *i;
      ++i;
    }
  return r;
}

std::string Helper::stringize( const std::vector<std::string> & s , const std::string & delim  )
{
  std::string r = "";
  std::vector<std::string>::const_iterator i = s.begin();
  while ( i != s.end() )
    {
      if ( i != s.begin() ) r += delim;
      r += *i;
      ++i;
    }
  return r;
}

std::string Helper::stringizeKeyPairList( const std::map<std::string,std::string> & m , bool show_keys )
{
  std::string s;
  if ( show_keys ) 
    {     
      std::map<std::string,std::string>::const_iterator i = m.begin();
      while ( i != m.end() )
	{
	  std::string v = i->second == "" ? "." : i->second;
	  if ( i != m.begin() ) s += ",";
	  s += i->first + "=" + v;	  
	  ++i;
	}
    }
  else
    {
      std::map<std::string,std::string>::const_iterator i = m.begin();
      while ( i != m.end() )
	{
	  if ( i->second == "" ) { ++i; continue; }
	  if ( s != "" ) s += ",";
	  s += i->second;	
	  ++i;
	}
      if ( s == "" ) s = ".";
    }
  return s;
}

std::map<std::string,std::string> Helper::quoted_comma_keypair_split( const std::string & item )
{  

  std::map<std::string,std::string> pairs;  
  std::vector<std::string> tok = Helper::quoted_char_split( item , ',' );
  for (int i=0; i<tok.size(); i++)
    {
      const std::string & p = tok[i];
      int ppos = p.find("=");
      if ( ppos != std::string::npos ) 
	pairs[ p.substr( 0 , ppos ) ] = p.substr( ppos+1 );
    }
  return pairs; 
}


set<std::string> Helper::parseCommaList(const std::string & lst)
{
  std::set<std::string> results;
  std::vector<std::string> tok = Helper::char_split( lst , ',' );
  std::vector<std::string>::iterator tok_iter = tok.begin();
  while ( tok_iter != tok.end() )
    {
      results.insert( *tok_iter );
      ++tok_iter;
    }
  return results;
}

bool Helper::ends_with( const std::string & v , const std::string & e )
{
  if ( v.size() < e.size() ) return false;
  return v.substr( v.size() - e.size() ) == e;
}

std::vector<std::string> Helper::whitespace(const std::string & sline)
{
  std::string buf;
  std::stringstream ss(sline);
  std::vector<std::string> tokens;
  while (ss >> buf)
    tokens.push_back(buf);
  return tokens;
}


std::string Helper::filelist2commalist( const std::string & f )
{

  // If a string starts with a @file, then replace "@file" with contents of file
  // can be a list of real values and files;
  
  //   @file1,val1,val2,@file2
  //   @file2,#value2    --> extract value file2 item from list
  //   @file1,#@file2    --> extract contents of file2 form file1
  
  // return a simple string of comma-separated values
  // final +3 means look up from col 3 (col 1 by default

  //  +file1+3,+file+3
  
  std::vector<std::string> o = Helper::quoted_parse( f , " \t\n" );

  // return string

  std::string rstr = "";

  for ( int i = 0 ; i < o.size() ; i++ )
    {
      
      if ( i > 0 ) rstr += " ";

      // this could be a key=value thing.
      std::vector<std::string> tok = Helper::quoted_parse( o[i] , "=" );
      if ( tok.size() < 1 || tok.size() > 2 ) continue;

      bool keyval = tok.size() == 2;

      std::string key = keyval ? (tok[0]+"=\"") : "";

      std::string & s2 = keyval ? tok[1] : tok[0];

      std::vector< std::string > l = Helper::parse( s2 , "," );
      
      std::set< std::string > keeps;
      std::set< std::string > removes;
      std::vector< bool > special( l.size(), false );
      
      for ( int j = 0 ; j < l.size(); j++ )
	{
	  std::string js = Helper::unquote( l[j] );

	  // Special character?
	  if ( js.substr(0,1) == "@" )
	    {
	      special[j] = true;
	      
	      if ( js.size() < 2 ) 
		{
		  continue;
		}
	      
	      // treat "++symbol"  as value "+symbol"
	      
	      if ( js.size() > 2 && js.substr(1,1) == "@" )
		{
		  keeps.insert( js.substr(1) );		  
		  continue;
		}
	      
	      // else insert contents of file
	      inserter( keeps , js.substr(1) );
	    }
	  else if ( js.substr(0,1) == "#" )
	    {
	      special[j] = true;
	      
	      // remove a file or value?
	      if ( js.size() >= 3 && js.substr(1,1) == "@" ) 
		{	      
		  inserter( removes , js.substr(2) );
		}
	      else
		{
		  removes.insert( js.substr(1) );
		}
	    }
	  
	}
      
      // construct string list
      
      rstr += key;
      bool extra = false;
      bool added_includes = false;
      for ( int j = 0 ; j < l.size(); j++ )
	{
	  if ( special[j] )
	    {
	      if ( ! added_includes ) 
		{
		  added_includes = true;
		  std::set< std::string>::iterator k = keeps.begin();
		  while ( k != keeps.end() )
		    {
		      if ( removes.find( *k ) == removes.end() )
			{
			  if ( extra ) rstr += ",";
			  extra = true;
			  rstr += *k;
			}
		      ++k;
		    }
		}
	    }
	  else 
	    {
	      // a normal entry

	      if ( extra ) rstr += ",";
	      extra = true;
	      rstr += Helper::unquote( l[j] );
	    }
	}

      // place quote back around whole thing
      if ( keyval) rstr += "\"";
    }

  return rstr;
}


void Helper::inserter( std::set< std::string > & strset , const std::string & filespec )
{
  int col = 0;
  int comma = filespec.find("#");
  std::string filename = filespec;
  if ( comma != string::npos )
    {
      filename = filespec.substr( 0 , comma ) ;
      if ( ! str2int( filespec.substr( comma+1 ) , col ) )
	{
	  plog.warn( "trouble with: " + filespec );
	  return;
	}

      // user gives us 1-based
      --col;
      if ( col < 0 ) 
	{
	  plog.warn( "trouble with column value: " + col );
	  return;
	}
    }
  
  if ( ! Helper::fileExists( filename ) ) 
    {
      plog.warn("could not find " + filename );
      return;
    }
  
  InFile IN1( filename );
  while ( ! IN1.eof() )
    {
      std::string sline;
      getline(IN1,sline);
      std::vector<std::string> cols = parse( sline, "\t" );
      if ( cols.size() <= col ) 
	{
	  if ( cols.size() != 0 )
	    plog.warn( filename + " row with " + Helper::int2str( cols.size() ) + " fields when field " + Helper::int2str(col+1) + " requested" ); 
	  continue;
	}
      strset.insert( cols[col] );
    }  
  IN1.close(); 
}


std::vector<std::string> Helper::parse(const string & item, const string & s , bool empty )
{  
  if ( s.size() == 1 ) return Helper::char_split( item , s[0] , empty ); 
  if ( s.size() == 2 ) return Helper::char_split( item , s[0] , s[1] , empty ); 
  if ( s.size() == 3 ) return Helper::char_split( item , s[0] , s[1] , s[2] , empty ); 
  Helper::halt("silly internal error in parse/char_split");
}  

std::vector<std::string> Helper::quoted_parse(const string & item , const std::string & s , bool empty )
{
  if ( s.size() == 1 ) return Helper::quoted_char_split( item , s[0] , empty ); 
  if ( s.size() == 2 ) return Helper::quoted_char_split( item , s[0] , s[1] , empty ); 
  if ( s.size() == 3 ) return Helper::quoted_char_split( item , s[0] , s[1] , s[2] , empty ); 
  Helper::halt("silly internal error in parse/char_split");
}


std::string Helper::remove_tags( const std::string & s)
{
  if ( s=="" ) return "";
  // convert "<text>" to "text"
  bool at_start = s.substr(0,1) == "<";
  bool at_end = s.substr( s.size()-1,1) == ">"; 
  if ( ! ( at_start || at_end ) ) return s;
  int len = s.size();
  if ( at_start ) --len;
  if ( at_end ) --len;
  return s.substr( at_start ? 1 : 0 , len );
}

std::string add_tags( const std::string & s)
{
  // convert "text" to "<text>"  
  return "<" + s + ">";
}


std::string Helper::coordinate( const int chr, 
				const int bp1 , 
				const int bp2  )
{
  std::string r = Helper::chrCode( chr );
  if ( bp1 > 0 ) r += ":" + Helper::int2str( bp1 );
  if ( bp2 > bp1 ) r += ".." + Helper::int2str( bp2 );
  return r;  
}


string Helper::chrCode( int c , bool prefix )
{
  // if available, use VARDB table
  if ( GP && GP->vardb.attached() ) return GP->vardb.chr_name(c);

  // otherwise, this will work for the standard human chromsomes 1-22 and X,Y,M (codes 23-25)
  if ( c == 23 ) return  prefix ? "chrX" : "X";
  else if ( c == 24 ) return prefix ? "chrY" : "Y";
  else if ( c == 25 ) return prefix ? "chrM" : "M";
  return prefix ? "chr"+ int2str(c) : int2str(c);
}

bool Helper::chr_known( const std::string & c )
{
  if ( GP && GP->vardb.attached() ) return GP->vardb.chr_known(c);
  return chrCode(c) != 0;
}


int Helper::chrCode(const std::string & c)
{

  // if available, use VARDB table
  if ( GP && GP->vardb.attached() ) return GP->vardb.chr_code(c);
  
  // Otherwise, use hard-coded human-specific values
  int cn;  
  if ( Helper::str2int( c , cn ) ) return cn;
  
  // Ignore chr1_random, etc  
  if( c.size() > 5 ) return 0;
  
  // Remove 'chr' prefix
  std::string c2 ="";
  if ( c.size() > 3 && c.substr(0,3) == "chr" )
    c2 = c.substr(3);
  
  if ( Helper::str2int( c2 , cn ) ) return cn;    
  if ( c2 == "X" ) return 23;
  if ( c2 == "Y" ) return 24;
  if ( c2 == "M" ) return 25;

  return 0;
}



//////////////////////////////////////////////////////////////////////
// File handling

void Helper::ensure_folder( std::string & f) 
{
  if ( f.substr( f.size() - 1 , 1 ) != "/" ) 
    f += "/";
}

std::string Helper::fullpath( const std::string & f)
{

#ifdef WINDOWS
  return f;
#endif
  
  if ( f == "" ) Helper::halt("missing filename in fullpath()");  
  if ( f.substr(0,1) != "/" ) return FileMap::working_directory() + "/" + f;
  return f;
}


bool Helper::fileExists( const std::string & f )
{
  std::ifstream inp;
  inp.open(f.c_str(), ifstream::in);
  if(inp.fail())
    {
      inp.clear(ios::failbit);
      inp.close();
      return false;
    }
  inp.close();
  return true;
}

bool Helper::fileExists(File * f)
{
  return fileExists( f->name() );
}

bool Helper::checkFileExists(File * f)
{
  return checkFileExists( f->name() );
}

bool Helper::checkFileExists( const std::string & f )
{
  ifstream inp;
  inp.open(f.c_str(), ifstream::in);
  if(inp.fail())
    {
      inp.clear(ios::failbit);
      inp.close();
      string msg = "No file [ " + f + " ] exists.";
      halt(msg);
    }
  inp.close();
  return true;
}

bool Helper::checkFileExists( const std::vector<std::string> & f )
{
  for (int k=0; k<f.size(); k++) checkFileExists( f[k] );
  return true;
}



//////////////////////////
// Pretty printing

std::string Helper::sw( const std::string & s , int n)
{
  int l = n - s.size();
  if ( l < 1 ) return " " + s;
  std::string t = s;
  t.insert( t.begin(), l , ' ' );
  return t;
}

std::string Helper::sw(double d , int n)
{
  std::string s = realnum(d) ? dbl2str(d) : "NA";
  int l = n - s.size();
  if ( l < 1 ) return " " + s;
  s.insert(s.begin(), l , ' ' );
  return s;
}

std::string Helper::sw(double d , int f, int n)
{
  std::string s = realnum(d) ? ( f < 0 ? dbl2str(d,-f) : dbl2str_fixed(d,f) ) : "NA";
  int l = n - s.size();
  if ( l < 1 ) return " " + s;
  s.insert(s.begin(), l , ' ' );
  return s;
}

std::string Helper::sw(int i , int n)
{
  std::string s = realnum(i) ? int2str(i) : "NA";
  int l = n - s.size();
  if ( l < 1 ) return " " + s;
  s.insert(s.begin(), l , ' ' );
  return s;
}

std::string Helper::sw(uint64_t i , int n)
{
  std::string s = uint64_t2str(i);
  int l = n - s.size();
  if ( l < 1 ) return " " + s;
  s.insert(s.begin(), l , ' ' );
  return s;
}

bool Helper::realnum(double d)
{
  double zero = 0;
  if (d != d || d == 1/zero || d == -1/zero)
    return false;
  else
    return true;
}



bool Helper::is_int(const std::string & s )
{
  if ( s == "Integer" ) return true;
  if ( s.size() < 3 ) return false;
  std::string t = s;
  Helper::str2upper(t);
  return t.substr(0,3) == "INT";
}

bool Helper::is_long(const std::string & s )
{
  if ( s == "Long" ) return true;
  if ( s == "Integer64" ) return true;
  std::string t = s;
  Helper::str2upper(t);
  if ( t == "LONG" ) return true;
  if ( t == "INT64" || t == "INTEGER64" ) return true;
  return false;
}

bool Helper::is_float(const std::string & s )
{
  if ( s == "Float" ) return true;
  std::string t = s;
  Helper::str2upper(t);
  return t == "FLOAT";
}

bool Helper::is_string(const std::string & s )
{
  if ( s == "String" ) return true;
  if ( s == "Text" ) return true;
  std::string t = s;
  Helper::str2upper(t);
  if ( t == "TEXT" ) return true;
  return t.substr(0,3) == "STR";
}

bool Helper::is_char(const std::string & s )
{
  if ( s.size() < 4 ) return false;
  std::string t = s;
  Helper::str2upper(t);
  return t.substr(0,4) == "CHAR";
}

bool Helper::is_text(const std::string & s )
{
  return is_string(s) || is_char(s);
}

bool Helper::is_bool(const std::string & s )
{
  if ( s == "Bool" ) return true;
  std::string t = s;
  Helper::str2upper(t);
  return t == "BOOL";
}

bool Helper::is_flag( const std::string & s )
{
  if ( s == "Flag" ) return true;
  std::string t = s;
  Helper::str2upper(t);
  return t == "FLAG";
}

double Helper::chi2x2(int a, int b, int c, int d)
{

  double row1 = a + b;
  double row2 = c + d;
  double col1 = a + c;
  double col2 = b + d;

  if ( row1 * row2 * col1 * col2 == 0 )
    return 0;

  double total = col1 + col2;

  double E_a = ( row1 * col1 ) / total;
  double E_b = ( row1 * col2 ) / total;
  double E_c = ( row2 * col1 ) / total;
  double E_d = ( row2 * col2 ) / total;

  return ((a-E_a)*(a-E_a) ) / E_a +
    ((b-E_b)*(b-E_b))/E_b +
    ((c-E_c)*(c-E_c))/E_c +
    ((d-E_d)*(d-E_d))/E_d ;

}


vector<string> Helper::load_string_list(const string & n)
{
  vector<string> r;
  if ( ! fileExists(n) ) return r ;
  std::ifstream f;
  f.open( n.c_str() , ios::in );
  while ( ! f.eof() )
    {
      string line;
      cin >> line;
      if ( line == "" ) continue;
      r.push_back(line);
    }
  f.close();
  return r;
}

/*
// This function implements an exact SNP test of Hardy-Weinberg
// Equilibrium as described in Wigginton, JE, Cutler, DJ, and
// Abecasis, GR (2005) A Note on Exact Tests of Hardy-Weinberg
// Equilibrium. American Journal of Human Genetics. 76: 000 - 000
//
// Written by Jan Wigginton
*/

double Helper::hwe( const Variant & v , int * phom1 , int * phets , int * phom2 )
{

  if ( ! v.biallelic() ) 
    {
      if ( phom1 ) *phom1 = 0;
      if ( phets ) *phets = 0;
      if ( phom2 ) *phom2 = 0;
      return 1;
    }
  int hets = 0, hom1 = 0 , hom2 = 0;
  for (int i=0; i<v.size(); i++)
    if ( ! v(i).null() )
      {
	int ac = v(i).allele_count( );
	if ( ac == 0 ) ++hom1;
	else if ( ac == 1 ) ++hets;
	else if ( ac == 2 ) ++hom2;
      }
  if ( phom1 ) *phom1 = hom1;
  if ( phets ) *phets = hets;
  if ( phom2 ) *phom2 = hom2;
  return Helper::SNPHWE( hets,hom1,hom2 );
}

double Helper::SNPHWE(int obs_hets, int obs_hom1, int obs_hom2)
{
  
  if (obs_hom1 + obs_hom2 + obs_hets == 0 ) return 1;
  
  if (obs_hom1 < 0 || obs_hom2 < 0 || obs_hets < 0) 
    {
      Helper::halt("Internal error: negative count in HWE test: "
		   + Helper::int2str(obs_hets)+" "
		   + Helper::int2str(obs_hom1)+" "
		   + Helper::int2str(obs_hom2));
    }

  int obs_homc = obs_hom1 < obs_hom2 ? obs_hom2 : obs_hom1;
  int obs_homr = obs_hom1 < obs_hom2 ? obs_hom1 : obs_hom2;

  int rare_copies = 2 * obs_homr + obs_hets;
  int genotypes   = obs_hets + obs_homc + obs_homr;

  double * het_probs = (double *) malloc((size_t) (rare_copies + 1) * sizeof(double));
  if (het_probs == NULL) 
    Helper::halt("Internal error: SNP-HWE: Unable to allocate array" );
  
  int i;
  for (i = 0; i <= rare_copies; i++)
    het_probs[i] = 0.0;
  
  /* start at midpoint */
  int mid = rare_copies * (2 * genotypes - rare_copies) / (2 * genotypes);
  
  /* check to ensure that midpoint and rare alleles have same parity */
  if ((rare_copies & 1) ^ (mid & 1))
    mid++;
  
  int curr_hets = mid;
  int curr_homr = (rare_copies - mid) / 2;
  int curr_homc = genotypes - curr_hets - curr_homr;

  het_probs[mid] = 1.0;
  double sum = het_probs[mid];
  for (curr_hets = mid; curr_hets > 1; curr_hets -= 2)
    {
      het_probs[curr_hets - 2] = het_probs[curr_hets] * curr_hets * (curr_hets - 1.0)
	/ (4.0 * (curr_homr + 1.0) * (curr_homc + 1.0));
      sum += het_probs[curr_hets - 2];

      /* 2 fewer heterozygotes for next iteration -> add one rare, one common homozygote */
      curr_homr++;
      curr_homc++;
    }

  curr_hets = mid;
  curr_homr = (rare_copies - mid) / 2;
  curr_homc = genotypes - curr_hets - curr_homr;
  for (curr_hets = mid; curr_hets <= rare_copies - 2; curr_hets += 2)
    {
      het_probs[curr_hets + 2] = het_probs[curr_hets] * 4.0 * curr_homr * curr_homc
	/((curr_hets + 2.0) * (curr_hets + 1.0));
      sum += het_probs[curr_hets + 2];
      
      /* add 2 heterozygotes for next iteration -> subtract one rare, one common homozygote */
      curr_homr--;
      curr_homc--;
    }
  
  for (i = 0; i <= rare_copies; i++)
    het_probs[i] /= sum;

  /* alternate p-value calculation for p_hi/p_lo
   double p_hi = het_probs[obs_hets];
   for (i = obs_hets + 1; i <= rare_copies; i++)
     p_hi += het_probs[i];
   
   double p_lo = het_probs[obs_hets];
   for (i = obs_hets - 1; i >= 0; i--)
      p_lo += het_probs[i];
   
   double p_hi_lo = p_hi < p_lo ? 2.0 * p_hi : 2.0 * p_lo;
  */

  double p_hwe = 0.0;
  /*  p-value calculation for p_hwe  */
  for (i = 0; i <= rare_copies; i++)
    {
      if (het_probs[i] > het_probs[obs_hets])
	continue;
      p_hwe += het_probs[i];
    }
   
  p_hwe = p_hwe > 1.0 ? 1.0 : p_hwe;

  free(het_probs);

  return p_hwe;
}


int_range::int_range( const std::string & s , const int smode )
{
  set(s, smode );
}

void int_range::reset( )
{
  lwr = upr = -1;
  has_lwr = has_upr = false;
}

void int_range::set( const std::string & s , const int smode )
{
  
  // constrained to positive values
  
  // expect 1-2
  // or     -2      means 0-2 
  //        2-            2- inf

  //        2    by default -->  -2  


  reset();

  //
  // explicitly set no filter?
  //

  if ( s == "-" || s == "" ) return;
  
  if ( s == "*" || s == "." )
    {
      // i.e. everything
      has_lwr = has_upr = false;
    }

  //
  // Parse, expecting at least one "-"
  //
  
  std::vector<std::string> tok = Helper::char_split( s , ':' , false );
  if ( tok.size() != 2 ) 
    {
      std::vector<std::string> tok2 = Helper::char_split( s , '-' , false );
      if ( tok2.size() == 2 ) tok = tok2;
    }
  
  if ( tok.size() == 2 )
    {
      has_lwr = Helper::str2int( tok[0] , lwr );
      has_upr = Helper::str2int( tok[1] , upr );
      if ( lwr > upr ) 
	{
	  int t = lwr;
	  lwr = upr;
	  upr = t;
	}
    }
  else if ( tok.size() == 1 ) // single char means upper bound (inc. below)
    {
      
      // 2- means lower bound, otherwise upper:
      
      if ( s.substr( s.size()-1 , 1 )  == "-" || s.substr( s.size()-1,1) == ":" )
	has_lwr = Helper::str2int( tok[0] , lwr );
      else if ( s.substr(0,1) == ":" || s.substr(0,1) == "-" )
	has_upr = Helper::str2int( tok[0] , upr );
      else
	{
	  if ( smode == -1 ) 
	    {
	      has_lwr = false; 
	      has_upr = Helper::str2int( tok[0] , upr ); 	  	      
	    }
	  else if ( smode == 1 ) 
	    {
	      has_lwr = Helper::str2int( tok[0] , lwr ); 
	      has_upr = false;
	    }
	  else
	    {
	      has_lwr = Helper::str2int( tok[0] , lwr ); 
	      has_upr = Helper::str2int( tok[0] , upr ); 	  
	    }
	}    
    }
}


bool int_range::in( const int i ) const
{
  if ( i < 0 ) return false;
  if ( has_lwr && i < lwr ) return false;
  if ( has_upr && i > upr ) return false;
  return true;
}


//
// prop_range
//

dbl_range::dbl_range( const std::string & s, const int smode )
{
  set(s , smode );
}

void dbl_range::reset( )
{
  lwr = 0; upr = 0;
  has_lwr = has_upr = false;
}
 
void dbl_range::set( const std::string & s , const int smode )
{

  // try comma first, as this will allow -ve numbers
  // -2,2

  bool colon = true;
  reset();

  if ( s == "-" || s == "" )
    {
      return;
    }

  if ( s == "*" || s == "." )
    {
      // i.e. everything
      has_lwr = has_upr = false;
      return;
    }
  
  std::vector<std::string> tok = Helper::char_split( s , ':' , false );
  if ( tok.size() != 2 )
    {
      std::vector<std::string> tok2 = Helper::char_split( s , '-' , false );
      if ( tok2.size() == 2 ) tok = tok2;
      colon = false;
    }
  
  if ( tok.size() == 2 )
    {
      has_lwr = Helper::str2dbl( tok[0] , lwr );
      has_upr = Helper::str2dbl( tok[1] , upr );
      if ( lwr > upr ) 
	{
	  double t = lwr;
	  lwr = upr;
	  upr = t;
	}
    }
  else if ( tok.size() == 1 ) // single char means upper bound (inc. below)
    {

      // ,2   means  * -- 2
      // 2,   means  2 -- *

      // 2    means  2 -- 2

      // -2   means  -2 -- -2
      // 2-   means  2  -- *
      
      // 2- means lower bound, otherwise upper:
      
      if ( ( colon && s.substr( s.size()-1 , 1 )  == ":" ) 
	   || ( (!colon) && s.substr( s.size()-1,1) == "-" ) )
	{
	  has_lwr = Helper::str2dbl( tok[0] , lwr );
	}
      else if ( s.substr(0,1) == ":" )
	has_upr = Helper::str2dbl( tok[0] , upr ); 
      else 
	{

	  if ( smode == -1 ) 
	    {
	      has_lwr = false; 
	      has_upr = Helper::str2dbl( tok[0] , upr ); 	  	      
	    }
	  else if ( smode == 1 ) 
	    {
	      has_lwr = Helper::str2dbl( tok[0] , lwr ); 
	      has_upr = false;
	    }
	  else
	    {
	      has_lwr = Helper::str2dbl( tok[0] , lwr ); 
	      has_upr = Helper::str2dbl( tok[0] , upr ); 	  
	    }
	}	
    }
  
}

bool dbl_range::in( const double d ) const
{
  //  if ( d < 0 || d > 1 ) return false;
  if ( has_lwr && d < lwr ) return false;
  if ( has_upr && d > upr ) return false;
  return true;
}



//
// delete a file
//

bool Helper::remove_file( const std::string & f )
{

  if( remove( f.c_str() ) == -1 ) 
    {
      plog.warn("could not delete file: " + f );
      return false; 
    }
  
  return true;
}


std::string Helper::quote_value( const std::string & s )
{
  // find first quote
  size_t p = s.find("=");  
  if ( p == std::string::npos ) return s;
  bool start_quote = s.substr(p+1,1) == "\"" ;
  bool end_quote = s.substr(s.size()-1) == "\"" ;
  return s.substr(0,p+1) + ( start_quote ? "" : "\"" ) + s.substr(p+1) + ( end_quote ? "" : "\"" );
}

std::string Helper::unquote( const std::string & s )
{
  if ( s == "" ) return s;
  int s1 = s.substr(0,1) == "\"" ? 1 : 0;
  int s2 = s.substr(s.size()-1) == "\"" ? s.size()-1-s1 : s.size()-s1;
  return s.substr(s1,s2);
}

void Helper::str2upper( std::string & s )
{
  std::string::iterator i = s.begin();
  std::string::iterator end = s.end();  
  while (i != end) {
    *i = std::toupper((unsigned char)*i);
    ++i;
  }
}


//
// C-style tokenizer for parsing VCFs
//

Helper::char_tok::char_tok() 
{ 
    s = NULL; clear(); 
}

Helper::char_tok::char_tok( const char_tok & rhs ) 
    : len(rhs.len) , d(rhs.d) , p(rhs.p) , s(NULL) , escape_quotes(rhs.escape_quotes) 
{
    if ( rhs.s )
    {	
	s = new char[ rhs.len + 1 ];    
	memcpy( s , rhs.s , rhs.len + 1 );
    }
} 
	  
Helper::char_tok & Helper::char_tok::operator= ( const Helper::char_tok &rhs ) 
{
    // clean up first any existing data
    if ( s ) delete [] s;
    s = NULL;
    len = rhs.len;
    d = rhs.d;
    p = rhs.p;
    escape_quotes = rhs.escape_quotes;
    if ( rhs.s ) 
    {	
	// use memcpy, as now we have many \0 delimiters within the string
	s = new char[ rhs.len + 1 ];	
	memcpy( s , rhs.s , rhs.len + 1 );
    } 
    return *this;
}


Helper::char_tok::char_tok( const std::string & istr , int * ps , const char d , bool eq ) 
    : d(d) , escape_quotes(eq)
{
    len = istr.size();
    init( istr.c_str(), ps );
}

Helper::char_tok::char_tok( const char * istr , int l , int * ps , const char d , bool eq ) 
    : d(d) , len(l) , escape_quotes(eq)
{     
    // if len is positive, assume that is length of 's'
    if ( ! len ) len = strlen(istr);    
    init( istr , ps );
}

Helper::char_tok::~char_tok()
{
    if ( s ) delete [] s;
}


void Helper::char_tok::init( const char * istr , int * ps )
{
    
    if ( ! istr ) { s = NULL; return; }

    // make a copy we can modify
    s = new char[ len+1 ] ;    
    strcpy( s, istr );
    
    p.clear();
    p.push_back(0); // first token starts at position 0	    
    
    // tokenize new string (+1 means position of next token)
    // consecutive tokens result in empty strings, as desired
    
    if ( escape_quotes )
    {
	bool inqt = false;
	for (int i=0; i < len; i++) 
	{
		if ( s[i] == '\"' ) inqt = ! inqt;
		if ( (!inqt) && s[i] == d ) { s[i] = '\0'; p.push_back(i+1); }	   
	    }
    }
    else // do not escape delimiters in quotes
    {
	for (int i=0; i < len; i++) 
	    if ( s[i] == d ) { s[i] = '\0'; p.push_back(i+1); }	   
    }

    // return number of tokens ( == # of delimiters + 1 )
    *ps = p.size(); 
    
}

void Helper::char_tok::clear()
{
    if ( s ) delete [] s;
    s = NULL;
    d = '\t';
    p.clear();
    len = 0;
}

