#include "helper.h"

#include <fstream>
#include <sstream>

bool fileExists(const std::string & f )
{
  std::ifstream inp;
  inp.open( f.c_str(), std::ifstream::in );
  if( inp.fail() )
    {
      inp.clear( std::ios::failbit );
      inp.close();
      return false;
    }
  inp.close();
  return true;
}


std::string int2str(int i)
{
  std::stringstream ss;
  ss << i;
  return ss.str();
}


std::vector<std::string> char_split( const std::string & s , const char c , bool empty )
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


std::vector<std::string> quoted_char_split( const std::string & s , const char c , bool empty )
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
