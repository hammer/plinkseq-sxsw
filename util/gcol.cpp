#include <iostream>
#include <cstdlib>
#include <vector>
#include <map>
#include <string>
#include <sstream>
#include <unistd.h>
#include <getopt.h>


const int MAXBUF = 50000;


class char_tok
{
  
  // takes a single char delimiter (assumed tab, or :)
  // always returns empty tokens
  // no allowance for quotes 
  
  // Usage:
  
  // int n;
  // char_tok tok( s , &n , ' ' );
  // for (int i = 0 ; i < n; i++) fs( tok(i) , cnt );
  //
  //  tok(i)  returns char* to \0-terminated C-string
  
  
public:  
  char_tok();
  char_tok( const char_tok & rhs );
  char_tok& operator= ( const char_tok &rhs );
  char_tok( const std::string & istr , int * ps , const char d = '\t' , bool eq = false );
  char_tok( const char * istr , int len , int * ps , const char d = '\t' , bool eq = false );
  ~char_tok();
  void init( const char * istr , int * ps );
  // return the i'th token as a C-style, \0 terminated
  const char * operator() (const int i) const { return s + p[i]; }
  const char * operator[] (const int i) const { return s + p[i]; }
  int size() const { return p.size(); } 
  void clear();
  
private:
  
  char * s;    
  int len;
  char d; 
  std::vector<int> p;
  bool escape_quotes;
};


char_tok::char_tok() 
{ 
  s = NULL; clear(); 
}

char_tok::char_tok( const char_tok & rhs ) 
  : len(rhs.len) , d(rhs.d) , p(rhs.p) , s(NULL) , escape_quotes(rhs.escape_quotes) 
{
  if ( rhs.s )
    {
      s = new char[ rhs.len + 1 ];    
      memcpy( s , rhs.s , rhs.len + 1 );
    }
} 
  
char_tok & char_tok::operator= ( const char_tok &rhs ) 
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


char_tok::char_tok( const std::string & istr , int * ps , const char d , bool eq ) 
  : d(d) , escape_quotes(eq)
{
  len = istr.size();
  init( istr.c_str(), ps );
}

char_tok::char_tok( const char * istr , int l , int * ps , const char d , bool eq ) 
  : d(d) , len(l) , escape_quotes(eq)
{     
  // if len is positive, assume that is length of 's'
  if ( ! len ) len = strlen(istr);    
  init( istr , ps );
}

char_tok::~char_tok()
{
  if ( s ) delete [] s;
}


void char_tok::init( const char * istr , int * ps )
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

void char_tok::clear()
{
  if ( s ) delete [] s;
  s = NULL;
  d = '\t';
  p.clear();
  len = 0;
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
  int col = 0;    // # cols in file
  int prtcol = 0; // # of cols to print

  while ( ! std::cin.eof() ) 
    {
      
      // parse a line?
      if ( col ) 
	{
	  
	  char line[MAXBUF];
	  std::cin.getline( line , MAXBUF , '\n' );
	  int n;
	  char_tok tok( line , &n , '\t' );  
	  
	  if ( n != col ) 
	    {
	      if ( n  > 1 ) std::cerr << "*** skipping bad line (" << n << " fields, not " << col << "): " << line << "\n";
	      continue;
	    }
	  
	  for (int i=0; i<prtcol; i++ ) 
	    {	      
	      if ( i ) std::cout << "\t";
	      if ( cols[i] == -1 ) std::cout << ordhdr[i];
	      else std::cout << tok( cols[i] ) ;
	    }
	  std::cout << "\n";
	}
      else // this is first row -- read as headers
	{
	  
	  std::string line;
	  std::getline( std::cin , line );
	  char_tok headers( line , &col , '\t' );	  
	  data.resize( col );
	  
	  bool first = true;
	  
	  std::map<std::string,int> hset;
	  for (int i=0; i<headers.size(); i++)
	    hset.insert( std::make_pair( headers[i] , i )) ;
	  
	  for (int i=0; i<ordhdr.size(); i++)
	    {
	      if ( hset.find( ordhdr[i] ) != hset.end() ) 
		{
		  if ( ! first ) std::cout << "\t";
		  first = false;
		  std::cout << ordhdr[i] ;
		  cols.push_back( hset[ ordhdr[i] ] );
		  std::cout << "adding col " << ordhdr[i] << "\n";
		}
	      else  // insert string literal
		{
		  if ( ! first ) std::cout << "\t";
		  first = false;
		  std::cout << ordhdr[i] ;
		  cols.push_back( -1 );
		  std::cout << "added str lit " << ordhdr[i] << "\n";
		}
	    }
	  
	  if ( cols.size() == 0 ) 
	    {
	      std::cerr << "gcol: no matching columns\n";
	      exit(1);
	    }

	  prtcol = cols.size();

	  std::cout << "\n";
	}
	
    }
  exit(0);
}

