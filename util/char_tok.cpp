
#include "char_tok.h"
#include <cstring>

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

