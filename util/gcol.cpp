#include <iostream>
#include <cstdlib>

#include <vector>
#include <set>
#include <map>
#include <string>
#include <sstream>

#include <unistd.h>
#include <getopt.h>
#include <cstdlib>


const int MAXBUF = 50000;

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
  
  std::vector<std::string> ordhdr;
    
  if ( argc < 2 ) 
    {
      std::cerr << "gcol: no columns specified\n";
      exit(1);
    }

  for ( int i = 1 ; i < argc ; i++)
      ordhdr.push_back( argv[i] );

  std::vector<int> cols;
  std::vector<std::string> data;
  int col = 0;

  while ( ! std::cin.eof() ) 
    {
      
      // parse a line?
      if ( col ) 
	{
	  if ( ! tokenize( data , col ) ) 
	    {
	      if ( data.size() > 0 ) std::cout << "bad line...\n";
	      continue;
	    }
	  
	  for (int i=0; i<cols.size(); i++ ) 
	    {
	      if ( i ) std::cout << "\t";
	      if ( cols[i] == -1 ) 
		std::cout << ordhdr[i];
	      else
		std::cout << data[ cols[i] ] ;
	    }
	  std::cout << "\n";
	}
      else // this is first row -- read as headers
	{
	  char line[MAXBUF];
	  std::cin.getline( line, MAXBUF, '\n' );
	  std::string buf;
	  std::string sline = line;
	  std::vector<std::string> headers = char_split( sline , '\t' );
	  col = headers.size();
	  data.resize( col );
	  
	  bool first = true;
	  
	  std::map<std::string,int> hset;
	  for (int i=0; i<headers.size(); i++)
	    hset.insert( make_pair( headers[i] , i )) ;
	  
	  for (int i=0; i<ordhdr.size(); i++)
	    {
	      if ( hset.find( ordhdr[i] ) != hset.end() ) 
		{
		  if ( ! first ) std::cout << "\t";
		  first = false;
		  std::cout << ordhdr[i] ;
		  cols.push_back( hset[ ordhdr[i] ] );
		}
	      else  // insert string literal
		{
		  if ( ! first ) std::cout << "\t";
		  first = false;
		  std::cout << ordhdr[i] ;
		  cols.push_back( -1 );
		}
	    }
	
	  if ( cols.size() == 0 ) 
	    {
	      std::cerr << "gcol: no matching columns\n";
	      exit(1);
	    }
	  
	  std::cout << "\n";
	}
	
    }
  exit(0);
}

