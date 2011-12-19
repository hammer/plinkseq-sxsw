#ifndef __CHAR_TOK_H__
#define __CHAR_TOK_H__

#include <vector>

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


#endif
