#include <iostream>
#include <cstdlib>

#include <vector>
#include <string>
#include <sstream>

#include <unistd.h>
#include <getopt.h>
#include <cstdlib>

const int MAXBUF = 10000;

std::vector<std::string> char_split( const std::string & s , const char c , bool empty = true )
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


bool tokenize( std::vector<std::string> & data , const int n )
{  
  char line[MAXBUF];
  std::cin.getline( line, MAXBUF, '\n' );  
  std::string sline = line;  
  data = char_split( sline , '\t' );  
  return data.size() == n; 
}

int main(int argc , char ** argv )
{
  
  int col = 0;
  int row = 0;

  std::vector<std::string> headers;
  std::vector<std::string> data;

  while ( ! std::cin.eof() ) 
    {
      // parse a line?
      if ( col ) 
	{
	  
	  if ( ! tokenize( data , col ) ) 
	    {
	      if ( data.size() > 0 ) std::cout << "read bad line...\n";
	      continue;
	    }
	  
	  ++row;
	  
	  for (int i=0; i<col; i++ ) 
	    std::cout << row << "\t" << i+1 << "\t" << headers[i] << "\t" << data[i] << "\n";
	  
	  std::cout << "\n";
	}
      else // this is first row -- read as headers
	{
	  char line[MAXBUF];
	  std::cin.getline( line, MAXBUF, '\n' );
	  std::string buf;
	  std::string sline = line;
	  headers = char_split( sline , '\t' );
	  col = headers.size();
	  data.resize( col );
	}
	
    }
  exit(0);
}

