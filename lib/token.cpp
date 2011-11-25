
#include "token.h"

#include <sstream>
#include <cmath>
#include <algorithm>

// TODO

// for now, vector x scalar comparisons do not allow for type conversion

void Token::set()
{
    ttype = UNDEF;
}

void Token::set( const std::string & s )
{
    ttype = STRING;
    sval = s;
}

void Token::set( const double d )
{
    ttype = FLOAT;
    fval = d;
}

void Token::set( const int i )
{
    ttype = INT;
    ival = i;
}

void Token::set( const bool b )
{
    ttype = BOOL;
    bval = b;
}


void Token::set( const std::vector<std::string> & s )
{
  if ( s.size() == 1 ) 
    set( s[0] );
  else
    {
      ttype = STRING_VECTOR;
      svec = s;
    }
}

void Token::set( const std::vector<double> & d )
{
  if ( d.size() == 1 ) 
    set( d[0] );
  else 
    {
      ttype = FLOAT_VECTOR;
      fvec = d;
    }
}

void Token::set( const std::vector<int> & i )
{
  if ( i.size() == 1 ) 
    set( i[0] ); 
  else 
    {
      ttype = INT_VECTOR;
      ivec = i;
    }
}

void Token::set( const std::vector<bool> & b )
{
  if (b.size()==1) 
    set(b[0]);
  else 
    { 
      ttype = BOOL_VECTOR;
      bvec = b;
    }
}


void Token::function( const std::string & fn )
{
  ttype = FUNCTION;
  tname = fn;
}

void Token::oper( Token::tok_type t )
{
  ttype = t;    
}

void Token::variable( const std::string & mf )
{
  ttype = VARIABLE;
  tname = mf;
}

void Token::init()
{

    tok_map[ "*" ]  = MULTIPLY_OPERATOR;
    tok_map[ "^" ]  = POWER_OPERATOR;
    tok_map[ "/" ]  = DIVIDE_OPERATOR;
    tok_map[ "%" ]  = MOD_OPERATOR;
    tok_map[ "%%" ] = MOD_OPERATOR;
    tok_map[ "+" ]  = ADD_OPERATOR;
    tok_map[ "-" ]  = SUBTRACT_OPERATOR;
    tok_map[ "&&" ] = AND_OPERATOR;
    tok_map[ "&" ]  = AND_OPERATOR;
    tok_map[ "||" ] = OR_OPERATOR;
    tok_map[ "|" ]  = OR_OPERATOR;
    tok_map[ "=" ]  = ASSIGNMENT_OPERATOR;
    tok_map[ "==" ] = EQUAL_OPERATOR;
    tok_map[ "!=" ] = UNEQUAL_OPERATOR;
    tok_map[ "!" ]  = NOT_OPERATOR;
    tok_map[ ">" ]  = GREATER_THAN_OPERATOR;
    tok_map[ ">=" ] = GREATER_THAN_OR_EQUAL_OPERATOR;
    tok_map[ "<" ]  = LESS_THAN_OPERATOR;
    tok_map[ "<=" ] = LESS_THAN_OR_EQUAL_OPERATOR;


    //
    // Reverse mapping
    //  
    
    std::map<std::string,Token::tok_type>::iterator i = tok_map.begin();
    while ( i != tok_map.end() ) 
      {
	tok_unmap[ i->second ] = i->first;
	++i;
      }
    

    //
    // Token function map
    //
    
    fn_map[ "set"  ]   = 1;  // number of args
    fn_map[ "sqrt" ]   = 1;  // square-root
    fn_map[ "sqr"  ]   = 1;  // X^2
    fn_map[ "pow"  ]   = 2;  // X^N    
    fn_map[ "ifelse" ] = 3;  // ifelse( cond , T , F )
    fn_map[ "n" ]      = 0;  // number of people in file
    
    // vector functions

    fn_map[ "element" ] = 2;  // element(Y,i)      extract element 'i' from vector 'Y'
    fn_map[ "length" ]  = 1;  // length(Y)   length of Y
    fn_map[ "min" ]     = 1;  // min(Y)      minimum element in Y
    fn_map[ "max" ]     = 1;  // max(Y)      max elenebt in Y
    fn_map[ "sum" ]     = 1;  // sum(Y)      sum of elements in Y
    fn_map[ "mean" ]    = 1;  // mean(Y)     mean of elements in Y
    fn_map[ "sort" ]    = 1;  // sort(Y)     returns sorted vector Y (asc.)

    // vector creation 
    fn_map[ "vec" ]    = 1;  // vec( '1,0,1' ) -- floating point vector
    fn_map[ "int" ]    = 1;  // vec( '1,0,1' )  ints
    fn_map[ "str" ]    = 1;  // vec( '1,0,1' )  strings
    fn_map[ "bool" ]   = 1;  // vec( '1,0,1' )  bools

    
    // here 'x' can be a complex expression, where has( Y < 2 ) 
    fn_map[ "any" ]    = 2;   // any( expr )    T/F if Y contains 1+ element matching 'x'
    fn_map[ "count" ]  = 2;   // count( expr )  return int of # of matches
    
}



bool Token::is_bool(bool * b ) const
{ 
    if ( ttype == BOOL ) 
    {
	if ( b ) *b = bval;
	return true;
    }
    return false;
}


bool Token::is_string( std::string * s ) const 
{ 
    if ( ttype == STRING ) 
    {
	if ( s ) *s = sval;
	return true;
    }
    return false;
}


bool Token::is_float( double * f ) const
{ 
    if ( ttype == FLOAT ) 
    {
	if ( f ) *f = fval;
	return true;
    }
    return false;
}


bool Token::is_int( int * i ) const
{ 
    if ( ttype == INT ) 
    {
	if ( i ) *i = ival;
	return true;
    }
    return false;
}



bool Token::is_bool_vector( std::vector<bool> * b ) const
{ 
    if ( ttype == BOOL_VECTOR ) 
    {
      if ( b ) *b = bvec;
      return true;
    }
    return false;
}

bool Token::is_string_vector( std::vector<std::string> * s ) const 
{ 
  if ( ttype == STRING_VECTOR ) 
    {
      if ( s ) *s = svec;
      return true;
    }
  return false;
}


bool Token::is_float_vector( std::vector<double> * f ) const
{ 
  if ( ttype == FLOAT_VECTOR ) 
    {
      if ( f ) *f = fvec;
      return true;
    }
  return false;
}


bool Token::is_int_vector( std::vector<int> * i ) const
{ 
  if ( ttype == INT_VECTOR ) 
    {
      if ( i ) *i = ivec;
      return true;
    }
  return false;
}



bool Token::is_operator() const
{
  // note -- treat parenthesis separately

  return ttype == EQUAL_OPERATOR ||
    ttype == UNEQUAL_OPERATOR ||
    ttype == ASSIGNMENT_OPERATOR ||
    ttype == NOT_OPERATOR ||
    ttype == AND_OPERATOR ||
    ttype == OR_OPERATOR ||     
    ttype == GREATER_THAN_OPERATOR ||
    ttype == GREATER_THAN_OR_EQUAL_OPERATOR ||
    ttype == LESS_THAN_OPERATOR ||
    ttype == LESS_THAN_OR_EQUAL_OPERATOR ||
    ttype == MOD_OPERATOR ||
    ttype == MULTIPLY_OPERATOR || 
    ttype == DIVIDE_OPERATOR ||
    ttype == ADD_OPERATOR ||
    ttype == SUBTRACT_OPERATOR;
}


bool Token::is_scalar() const
{
  return ttype == INT || ttype == FLOAT || ttype == STRING || ttype == BOOL;
}

bool Token::is_vector() const
{
  return ttype == INT_VECTOR || ttype == FLOAT_VECTOR || ttype == STRING_VECTOR || ttype == BOOL_VECTOR;
}


bool Token::is_function() const
{
  return ttype == FUNCTION;
}

bool Token::is_ident() const
{
  return ! ( is_operator() || is_function() || is_left_paren() || is_right_paren() || is_separator() );
}

bool Token::is_variable() const
{
  return ttype == VARIABLE;
}

Token::Token( const std::string & s )
{
  ttype = STRING;
  sval = s;
  init();
}

Token::Token( const double d )
{
  ttype = FLOAT;
  fval = d;
  init();
}

Token::Token( const int i )
{
  ttype = INT;
  ival = i;
  init();
}

Token::Token( const bool b )
{
  ttype = BOOL;
  bval = b;
  init();
}


Token::Token( const std::vector<std::string> & s )
{
  ttype = STRING_VECTOR;
  svec = s;
  init();
}

Token::Token( const std::vector<double> & d )
{
  ttype = FLOAT_VECTOR;
  fvec = d;
  init();
}

Token::Token( const std::vector<int> & i )
{
  ttype = INT_VECTOR;
  ivec = i;
  init();
}

Token::Token( const std::vector<bool> & b )
{
  ttype = BOOL_VECTOR;
  bvec = b;
  init();
}



Token::Token( const Token & rhs )
{
    *this = rhs;
}


Token Token::operator!() const
{
  // handles only bools and ints
  // both scalars and vectors
  
  if ( is_bool() ) return Token( ! bval ); 
  else if ( is_int() ) return Token( ival == 0 );
  else if ( is_bool_vector() )
    {
      std::vector<bool> ans( bvec.size() );
      for ( int i=0; i<bvec.size(); i++) ans[i] = ! bvec[i];
      return Token( ans );
    }  
  else if ( is_int_vector() )
    {
      std::vector<bool> ans( ivec.size() );
      for ( int i=0; i<ivec.size(); i++) ans[i] = ! ivec[i];
      return Token( ans );
    }  
  else return Token();
}


Token & Token::operator=(const Token & rhs)
{
  Token ret;
  
  ttype = rhs.ttype;
  tname = rhs.tname;
  
  ival = rhs.ival;
  sval = rhs.sval;
  fval = rhs.fval;
  bval = rhs.bval;
  
  ivec = rhs.ivec;
  svec = rhs.svec;
  fvec = rhs.fvec;
  bvec = rhs.bvec;
  
  return *this;
}


int Token::size() const
{
  if      ( is_scalar() ) return 1;
  else if ( ! is_vector() ) return 0;
  else if ( ttype == INT_VECTOR ) return ivec.size();
  else if ( ttype == FLOAT_VECTOR ) return fvec.size();
  else if ( ttype == STRING_VECTOR ) return svec.size();
  else if ( ttype == BOOL_VECTOR ) return bvec.size();
  return 0;
}


Token Token::operator!=(const Token & rhs ) const
{
  // vector x vector comparison not defined

  if ( is_vector() && rhs.is_vector() ) return Token();

  // vector != scalar 

  if ( is_vector() )
    {
      const int sz = size();
      if ( sz == 0 ) return Token();
      std::vector<bool> ans( sz );      
      if      ( is_int_vector() )    for (int i=0; i<sz; i++) ans[i] = ivec[i] != rhs.ival;
      else if ( is_float_vector() )  for (int i=0; i<sz; i++) ans[i] = fvec[i] != rhs.fval; 
      else if ( is_string_vector() ) for (int i=0; i<sz; i++) ans[i] = svec[i] != rhs.sval;  
      else if ( is_bool_vector() )   for (int i=0; i<sz; i++) ans[i] = bvec[i] != rhs.bval; 
      return Token( ans );
    }
  
  // scalar != vector

  if ( rhs.is_vector() )
    {
      const int sz = rhs.size();
      if ( sz == 0 ) return Token();
      std::vector<bool> ans( sz );      
      if      ( rhs.is_int_vector() )    for (int i=0; i<sz; i++) ans[i] = rhs.ivec[i] != ival;
      else if ( rhs.is_float_vector() )  for (int i=0; i<sz; i++) ans[i] = rhs.fvec[i] != fval;
      else if ( rhs.is_string_vector() ) for (int i=0; i<sz; i++) ans[i] = rhs.svec[i] != sval;
      else if ( rhs.is_bool_vector() )   for (int i=0; i<sz; i++) ans[i] = rhs.bvec[i] != bval;
      return Token( ans );      
    }

      
  // scalar x scalar
  if ( is_bool()   && rhs.is_bool()   )  return Token( bval != rhs.bval ); 
  if ( is_int()    && rhs.is_int()    )  return Token( ival != rhs.ival ); 
  if ( is_float()  && rhs.is_float()  )  return Token( fval != rhs.fval ); 
  if ( is_string() && rhs.is_string() )  return Token( sval != rhs.sval ); 

  // also allow int / bool and  int / float comparisons
  if ( is_int()   && rhs.is_bool()    )  return Token( ival != rhs.bval );
  if ( is_bool()  && rhs.is_int()     )  return Token( bval != rhs.ival );

  if ( is_float() && rhs.is_int()     )  return Token( fval != rhs.ival );
  if ( is_int()   && rhs.is_float()   )  return Token( ival != rhs.fval );

  return Token();
}


Token Token::operator==(const Token & rhs) const
{

  // vector x vector comparison not defined
  if ( is_vector() && rhs.is_vector() ) return Token();

  // vector != scalar 

  if ( is_vector() )
    {
      const int sz = size();
      if ( sz == 0 ) return Token();
      std::vector<bool> ans( sz );      
      if      ( is_int_vector() )    for (int i=0; i<sz; i++) ans[i] = ivec[i] == rhs.ival;
      else if ( is_float_vector() )  for (int i=0; i<sz; i++) ans[i] = fvec[i] == rhs.fval; 
      else if ( is_string_vector() ) for (int i=0; i<sz; i++) ans[i] = svec[i] == rhs.sval;  
      else if ( is_bool_vector() )   for (int i=0; i<sz; i++) ans[i] = bvec[i] == rhs.bval; 
      return Token( ans );
    }
  
  // scalar != vector

  if ( rhs.is_vector() )
    {
      const int sz = rhs.size();
      if ( sz == 0 ) return Token();
      std::vector<bool> ans( sz );      
      if      ( rhs.is_int_vector() )    for (int i=0; i<sz; i++) ans[i] = rhs.ivec[i] == ival;
      else if ( rhs.is_float_vector() )  for (int i=0; i<sz; i++) ans[i] = rhs.fvec[i] == fval;
      else if ( rhs.is_string_vector() ) for (int i=0; i<sz; i++) ans[i] = rhs.svec[i] == sval;
      else if ( rhs.is_bool_vector() )   for (int i=0; i<sz; i++) ans[i] = rhs.bvec[i] == bval;
      return Token( ans );      
    }

      
  // scalar x scalar
  if ( is_bool()   && rhs.is_bool()   )  return Token( bval == rhs.bval ); 
  if ( is_int()    && rhs.is_int()    )  return Token( ival == rhs.ival ); 
  if ( is_float()  && rhs.is_float()  )  return Token( fval == rhs.fval ); 
  if ( is_string() && rhs.is_string() )  return Token( sval == rhs.sval ); 

  // also allow int / bool and  int / float comparisons
  if ( is_int()   && rhs.is_bool()    )  return Token( ival == rhs.bval );
  if ( is_bool()  && rhs.is_int()     )  return Token( bval == rhs.ival );

  if ( is_float() && rhs.is_int()     )  return Token( fval == rhs.ival );
  if ( is_int()   && rhs.is_float()   )  return Token( ival == rhs.fval );

  return Token();

}

Token Token::operator+(const Token & rhs) const
{

  // vector x vector comparison not defined
  if ( is_vector() && rhs.is_vector() ) return Token();

  // concatenate strings ( in which case, A+B != B+A, so do before next step)
  if ( is_string_vector() && rhs.is_string() ) 
    {
      const int sz = size();
      if ( sz == 0 ) return Token();
      std::vector<std::string> ans( sz );      
      for (int i=0; i<sz; i++) ans[i] = svec[i] + rhs.sval;
      return Token( ans );            
    }
  
  if ( is_string() && rhs.is_string_vector() )
    {
      const int sz = rhs.size();
      if ( sz == 0 ) return Token();
      std::vector<std::string> ans( sz );      
      for (int i=0; i<sz; i++) ans[i] = sval + rhs.svec[i];
      return Token( ans );            
    }
  
  if ( is_int_vector() ) 
    {
      const int sz = size();
      if ( sz == 0 ) return Token();
      std::vector<int> ans( sz );      
      if      ( rhs.is_int() )   for (int i=0; i<sz; i++) ans[i] = ivec[i] + rhs.ival;
      else if ( rhs.is_float() ) for (int i=0; i<sz; i++) ans[i] = ivec[i] + (int)rhs.fval;
      else if ( rhs.is_bool() )  for (int i=0; i<sz; i++) ans[i] = ivec[i] + (int)rhs.bval;
      return Token( ans );            
    }
  
  if ( rhs.is_int_vector() ) 
    {
      const int sz = rhs.size();
      if ( sz == 0 ) return Token();
      std::vector<int> ans( sz );      
      if      ( is_int() )   for (int i=0; i<sz; i++) ans[i] = ival      + rhs.ivec[i];
      else if ( is_float() ) for (int i=0; i<sz; i++) ans[i] = (int)fval + rhs.ivec[i];
      else if ( is_bool() )  for (int i=0; i<sz; i++) ans[i] = (int)bval + rhs.ivec[i];
      return Token( ans );            
    }

  else if ( is_float_vector() )
    {
      const int sz = size();
      if ( sz == 0 ) return Token();
      std::vector<double> ans( sz );      
      if ( rhs.is_int() )        for (int i=0; i<sz; i++) ans[i] = fvec[i] + rhs.ival;
      else if ( rhs.is_float() ) for (int i=0; i<sz; i++) ans[i] = fvec[i] + rhs.fval;
      else if ( rhs.is_bool() )  for (int i=0; i<sz; i++) ans[i] = fvec[i] + (double)rhs.bval;
      return Token( ans );                  
    }

  else if ( rhs.is_float_vector() ) 
    {
      const int sz = rhs.size();
      if ( sz == 0 ) return Token();
      std::vector<double> ans( sz );      
      if      ( is_int() )   for (int i=0; i<sz; i++) ans[i] = ival         + rhs.fvec[i];
      else if ( is_float() ) for (int i=0; i<sz; i++) ans[i] = fval         + rhs.fvec[i];
      else if ( is_bool() )  for (int i=0; i<sz; i++) ans[i] = (double)bval + rhs.fvec[i];
      return Token( ans );            
    }

  else if ( is_bool_vector() )
    {
      const int sz = size();
      if ( sz == 0 ) return Token();
      std::vector<double> ans( sz );      
      if ( rhs.is_int() )        for (int i=0; i<sz; i++) ans[i] = bvec[i] + rhs.ival;
      else if ( rhs.is_float() ) for (int i=0; i<sz; i++) ans[i] = bvec[i] + rhs.fval;
      else if ( rhs.is_bool() )  for (int i=0; i<sz; i++) ans[i] = bvec[i] + (double)rhs.bval;
      return Token( ans );                  
    }

  else if ( rhs.is_bool_vector() ) 
    {
      const int sz = rhs.size();
      if ( sz == 0 ) return Token();
      std::vector<double> ans( sz );      
      if      ( is_int() )   for (int i=0; i<sz; i++) ans[i] = ival         + rhs.bvec[i];
      else if ( is_float() ) for (int i=0; i<sz; i++) ans[i] = fval         + rhs.bvec[i];
      else if ( is_bool() )  for (int i=0; i<sz; i++) ans[i] = (double)bval + rhs.bvec[i];
      return Token( ans );            
    }


  // scalar + scalar

  if ( is_int() ) 
    {
      if ( rhs.is_int() ) return Token( ival + rhs.ival );
      if ( rhs.is_float() ) return Token( ival + rhs.fval );
    }

  if ( is_float() ) 
    {
      if ( rhs.is_int() ) return Token( fval + rhs.ival );
      if ( rhs.is_float() ) return Token( fval + rhs.fval );
    }

  if ( is_string() ) // concatenate
    {
      if ( rhs.is_string() ) return Token( sval + rhs.sval );
    }

  return Token();
}

Token Token::operator-(const Token & rhs) const
{

  // vector x vector comparison not defined
  if ( is_vector() && rhs.is_vector() ) return Token();

  // vector x scalar ops
  if ( is_int_vector() ) 
    {
      const int sz = size();
      if ( sz == 0 ) return Token();
      std::vector<int> ans( sz );      
      if      ( rhs.is_int() )   for (int i=0; i<sz; i++) ans[i] = ivec[i] - rhs.ival;
      else if ( rhs.is_float() ) for (int i=0; i<sz; i++) ans[i] = ivec[i] - (int)rhs.fval;
      else if ( rhs.is_bool() )  for (int i=0; i<sz; i++) ans[i] = ivec[i] - (int)rhs.bval;
      return Token( ans );            
    }
  
  if ( rhs.is_int_vector() ) 
    {
      const int sz = rhs.size();
      if ( sz == 0 ) return Token();
      std::vector<int> ans( sz );      
      if      ( is_int() )   for (int i=0; i<sz; i++) ans[i] = ival      - rhs.ivec[i];
      else if ( is_float() ) for (int i=0; i<sz; i++) ans[i] = (int)fval - rhs.ivec[i];
      else if ( is_bool() )  for (int i=0; i<sz; i++) ans[i] = (int)bval - rhs.ivec[i];
      return Token( ans );            
    }

  else if ( is_float_vector() )
    {
      const int sz = size();
      if ( sz == 0 ) return Token();
      std::vector<double> ans( sz );      
      if ( rhs.is_int() )        for (int i=0; i<sz; i++) ans[i] = fvec[i] - rhs.ival;
      else if ( rhs.is_float() ) for (int i=0; i<sz; i++) ans[i] = fvec[i] - rhs.fval;
      else if ( rhs.is_bool() )  for (int i=0; i<sz; i++) ans[i] = fvec[i] - (double)rhs.bval;
      return Token( ans );                  
    }

  else if ( rhs.is_float_vector() ) 
    {
      const int sz = rhs.size();
      if ( sz == 0 ) return Token();
      std::vector<double> ans( sz );      
      if      ( is_int() )   for (int i=0; i<sz; i++) ans[i] = ival         - rhs.fvec[i];
      else if ( is_float() ) for (int i=0; i<sz; i++) ans[i] = fval         - rhs.fvec[i];
      else if ( is_bool() )  for (int i=0; i<sz; i++) ans[i] = (double)bval - rhs.fvec[i];
      return Token( ans );            
    }


  // scalar + scalar

  if ( is_int() ) 
    {
      if ( rhs.is_int()   ) return Token( ival - rhs.ival );
      if ( rhs.is_float() ) return Token( ival - rhs.fval );
    }

  if ( is_float() ) 
    {
      if ( rhs.is_int()   ) return Token( fval - rhs.ival );
      if ( rhs.is_float() ) return Token( fval - rhs.fval );
    }

  return Token();

}

Token Token::operator*(const Token & rhs) const
{

  // vector x vector comparison not defined
  if ( is_vector() && rhs.is_vector() ) return Token();

  // vector x scalar ops
  if ( is_int_vector() ) 
    {
      const int sz = size();
      if ( sz == 0 ) return Token();
      std::vector<int> ans( sz );      
      if      ( rhs.is_int() )   for (int i=0; i<sz; i++) ans[i] = ivec[i] * rhs.ival;
      else if ( rhs.is_float() ) for (int i=0; i<sz; i++) ans[i] = ivec[i] * (int)rhs.fval;
      else if ( rhs.is_bool() )  for (int i=0; i<sz; i++) ans[i] = ivec[i] * (int)rhs.bval;
      return Token( ans );            
    }
  
  if ( rhs.is_int_vector() ) 
    {
      const int sz = rhs.size();
      if ( sz == 0 ) return Token();
      std::vector<int> ans( sz );      
      if      ( is_int() )   for (int i=0; i<sz; i++) ans[i] = ival      * rhs.ivec[i];
      else if ( is_float() ) for (int i=0; i<sz; i++) ans[i] = (int)fval * rhs.ivec[i];
      else if ( is_bool() )  for (int i=0; i<sz; i++) ans[i] = (int)bval * rhs.ivec[i];
      return Token( ans );            
    }

  else if ( is_float_vector() )
    {
      const int sz = size();
      if ( sz == 0 ) return Token();
      std::vector<double> ans( sz );      
      if ( rhs.is_int() )        for (int i=0; i<sz; i++) ans[i] = fvec[i] * rhs.ival;
      else if ( rhs.is_float() ) for (int i=0; i<sz; i++) ans[i] = fvec[i] * rhs.fval;
      else if ( rhs.is_bool() )  for (int i=0; i<sz; i++) ans[i] = fvec[i] * (double)rhs.bval;
      return Token( ans );                  
    }

  else if ( rhs.is_float_vector() ) 
    {
      const int sz = rhs.size();
      if ( sz == 0 ) return Token();
      std::vector<double> ans( sz );      
      if      ( is_int() )   for (int i=0; i<sz; i++) ans[i] = ival         * rhs.fvec[i];
      else if ( is_float() ) for (int i=0; i<sz; i++) ans[i] = fval         * rhs.fvec[i];
      else if ( is_bool() )  for (int i=0; i<sz; i++) ans[i] = (double)bval * rhs.fvec[i];
      return Token( ans );            
    }

  else if ( is_bool_vector() )
    {
      const int sz = size();
      if ( sz == 0 ) return Token();
      std::vector<double> ans( sz );      
      if      ( rhs.is_int() )   for (int i=0; i<sz; i++) ans[i] = bvec[i] * rhs.ival;
      else if ( rhs.is_float() ) for (int i=0; i<sz; i++) ans[i] = bvec[i] * rhs.fval;
      else if ( rhs.is_bool() )  for (int i=0; i<sz; i++) ans[i] = bvec[i] * (double)rhs.bval;
      return Token( ans );                  
    }

  else if ( rhs.is_bool_vector() ) 
    {
      const int sz = rhs.size();
      if ( sz == 0 ) return Token();
      std::vector<double> ans( sz );      
      if      ( is_int() )   for (int i=0; i<sz; i++) ans[i] = ival         * rhs.bvec[i];
      else if ( is_float() ) for (int i=0; i<sz; i++) ans[i] = fval         * rhs.bvec[i];
      else if ( is_bool() )  for (int i=0; i<sz; i++) ans[i] = (double)bval * rhs.bvec[i];
      return Token( ans );            
    }

  // scalar + scalar

  if ( is_int() ) 
    {
      if ( rhs.is_int()   ) return Token( ival * rhs.ival );
      if ( rhs.is_float() ) return Token( ival * rhs.fval );
      if ( rhs.is_bool()  ) return Token( ival * rhs.bval );
    }

  if ( is_float() ) 
    {
      if ( rhs.is_int()   ) return Token( fval * rhs.ival );
      if ( rhs.is_float() ) return Token( fval * rhs.fval );
      if ( rhs.is_bool()  ) return Token( fval * rhs.bval );
    }

  if ( is_bool() ) 
    {
      if ( rhs.is_int()   ) return Token( bval * rhs.ival );
      if ( rhs.is_float() ) return Token( bval * rhs.fval );
      if ( rhs.is_bool()  ) return Token( bval * rhs.bval );
    }

  return Token();
}

Token Token::operator^(const Token & rhs) const
{
  if ( rhs.is_vector() ) Helper::halt( "not allowed vector expression 'x' ^ vector" );

  // vector ^ scalar
  // scalar ^ scalar

  if ( is_int_vector() ) 
    {
      const int sz = size();
      if ( sz == 0 ) return Token();
      std::vector<double> ans( sz );      
      if      ( rhs.is_int() )   for (int i=0; i<sz; i++) ans[i] = pow( ivec[i] , rhs.ival );
      else if ( rhs.is_float() ) for (int i=0; i<sz; i++) ans[i] = pow( ivec[i] , rhs.fval );
      return Token( ans );
    }

  if ( is_float_vector() ) 
    {
      const int sz = size();
      if ( sz == 0 ) return Token();
      std::vector<double> ans( sz );      
      if      ( rhs.is_int() )   for (int i=0; i<sz; i++) ans[i] = pow( fvec[i] , rhs.ival );
      else if ( rhs.is_float() ) for (int i=0; i<sz; i++) ans[i] = pow( fvec[i] , rhs.fval );
      return Token( ans );
    }

  if ( is_int() ) 
    {
      if ( rhs.is_int() ) return Token( pow( ival , rhs.ival ) );
      if ( rhs.is_float() ) return Token( pow( ival , rhs.fval ) );
    }
  if ( is_float() ) 
    {
      if ( rhs.is_int() ) return Token( pow( fval , rhs.ival ) );
      if ( rhs.is_float() ) return Token( pow( fval , rhs.fval ) );
    }
  return Token();
}

Token Token::operator/(const Token & rhs) const
{
  // vector x vector comparison not defined
  if ( is_vector() && rhs.is_vector() ) return Token();

  // vector x scalar ops
  if ( is_int_vector() ) 
    {
      const int sz = size();
      if ( sz == 0 ) return Token();
      std::vector<double> ans( sz );      
      if      ( rhs.is_int() )   for (int i=0; i<sz; i++) ans[i] = ivec[i] / (double)rhs.ival;
      else if ( rhs.is_float() ) for (int i=0; i<sz; i++) ans[i] = ivec[i] / rhs.fval;
      else if ( rhs.is_bool() )  for (int i=0; i<sz; i++) ans[i] = ivec[i] / (double)rhs.bval;
      return Token( ans );            
    }
  
  if ( rhs.is_int_vector() ) 
    {
      const int sz = rhs.size();
      if ( sz == 0 ) return Token();
      std::vector<double> ans( sz );      
      if      ( is_int() )   for (int i=0; i<sz; i++) ans[i] = ival / (double)rhs.ivec[i];
      else if ( is_float() ) for (int i=0; i<sz; i++) ans[i] = fval / (double)rhs.ivec[i];
      else if ( is_bool() )  for (int i=0; i<sz; i++) ans[i] = bval / (double)rhs.ivec[i];
      return Token( ans );            
    }

  else if ( is_float_vector() )
    {
      const int sz = size();
      if ( sz == 0 ) return Token();
      std::vector<double> ans( sz );      
      if ( rhs.is_int() )        for (int i=0; i<sz; i++) ans[i] = fvec[i] / rhs.ival;
      else if ( rhs.is_float() ) for (int i=0; i<sz; i++) ans[i] = fvec[i] / rhs.fval;
      else if ( rhs.is_bool() )  for (int i=0; i<sz; i++) ans[i] = fvec[i] / (double)rhs.bval;
      return Token( ans );                  
    }

  else if ( rhs.is_float_vector() ) 
    {
      const int sz = rhs.size();
      if ( sz == 0 ) return Token();
      std::vector<double> ans( sz );      
      if      ( is_int() )   for (int i=0; i<sz; i++) ans[i] = ival         / rhs.fvec[i];
      else if ( is_float() ) for (int i=0; i<sz; i++) ans[i] = fval         / rhs.fvec[i];
      else if ( is_bool() )  for (int i=0; i<sz; i++) ans[i] = (double)bval / rhs.fvec[i];
      return Token( ans );            
    }


  // scalar / scalar  (always return a double)

  if ( is_int() ) 
    {
      if ( rhs.is_int()   ) return Token( ival / (double)rhs.ival );
      if ( rhs.is_float() ) return Token( ival / rhs.fval );
    }

  if ( is_float() ) 
    {
      if ( rhs.is_int()   ) return Token( fval / (double)rhs.ival );
      if ( rhs.is_float() ) return Token( fval / rhs.fval );
    }

  return Token();

}

Token Token::operator%(const Token & rhs) const
{

  if ( rhs.is_vector() ) Helper::halt( "not allowed vector expression 'x' % vector" );
  if ( ! rhs.is_int() ) return Token();

  // scalar % scalar
  // vector % scalar
  
  if ( is_int_vector() ) 
    {
      const int sz = size();
      if ( sz == 0 ) return Token();
      std::vector<int> ans( sz );      
      for (int i=0; i<sz; i++) ans[i] = (int)(ival % rhs.ival );
      return Token(ans);
    }
  
  if ( is_int() ) return Token( (int)(ival % rhs.ival) );

  return Token();
}


Token Token::operator<(const Token & rhs) const
{

  // vector x vector comparison not defined
  if ( is_vector() && rhs.is_vector() ) return Token();
  
  // vector x scalar ops
  if ( is_int_vector() ) 
    {
      const int sz = size();
      if ( sz == 0 ) return Token();
      std::vector<bool> ans( sz );      
      if      ( rhs.is_int() )   for (int i=0; i<sz; i++) ans[i] = ivec[i] < rhs.ival;
      else if ( rhs.is_float() ) for (int i=0; i<sz; i++) ans[i] = ivec[i] < rhs.fval;
      return Token( ans );            
    }
  
  else if ( rhs.is_int_vector() ) 
    {
      const int sz = rhs.size();
      if ( sz == 0 ) return Token();
      std::vector<bool> ans( sz );
      if      ( is_int() )   for (int i=0; i<sz; i++) ans[i] = ival < rhs.ivec[i];
      else if ( is_float() ) for (int i=0; i<sz; i++) ans[i] = fval < rhs.ivec[i];
      return Token( ans ); 
    }

  else if ( is_float_vector() )
    {
      const int sz = size();
      if ( sz == 0 ) return Token();
      std::vector<bool> ans( sz );      
      if ( rhs.is_int() )        for (int i=0; i<sz; i++) ans[i] = fvec[i] < rhs.ival;
      else if ( rhs.is_float() ) for (int i=0; i<sz; i++) ans[i] = fvec[i] < rhs.fval;
      return Token( ans );                  
    }

  else if ( rhs.is_float_vector() ) 
    {
      const int sz = rhs.size();
      if ( sz == 0 ) return Token();
      std::vector<bool> ans( sz );      
      if      ( is_int() )   for (int i=0; i<sz; i++) ans[i] = ival         < rhs.fvec[i];
      else if ( is_float() ) for (int i=0; i<sz; i++) ans[i] = fval         < rhs.fvec[i];
      return Token( ans );            
    }


  else if ( is_string_vector() && rhs.is_string() ) 
    {
      const int sz = size();
      if ( sz == 0 ) return Token();
      std::vector<bool> ans( sz );      
      for (int i=0; i<sz; i++) ans[i] = svec[i] < rhs.sval;
      return Token( ans );            
    }
  
  else if ( is_string() && rhs.is_string_vector() ) 
    {
      const int sz = size();
      if ( sz == 0 ) return Token();
      std::vector<bool> ans( sz );      
      for (int i=0; i<sz; i++) ans[i] = svec[i] < rhs.sval;
      return Token( ans );            
    }

  
  // scalar < scalar  
  else if ( is_int() ) 
    {
      if ( rhs.is_int() )  return Token( ival < rhs.ival );
      if ( rhs.is_float() ) return Token( ival < rhs.fval );
    }
  
  else if ( is_float() ) 
    {
      if ( rhs.is_int() )  return Token( fval < rhs.ival );
      if ( rhs.is_float() ) return Token( fval < rhs.fval );
    }
  
  else if ( is_string() ) 
    {
      if ( rhs.is_string() )  return Token( sval < rhs.sval );
    }
  
  return Token();
}


Token Token::operator>(const Token & rhs) const
{

  // vector x vector comparison not defined
  if ( is_vector() && rhs.is_vector() ) return Token();
  
  // vector x scalar ops
  if ( is_int_vector() ) 
    {
      const int sz = size();
      if ( sz == 0 ) return Token();
      std::vector<bool> ans( sz );      
      if      ( rhs.is_int() )   for (int i=0; i<sz; i++) ans[i] = ivec[i] > rhs.ival;
      else if ( rhs.is_float() ) for (int i=0; i<sz; i++) ans[i] = ivec[i] > rhs.fval;
      return Token( ans );            
    }
  
  else if ( rhs.is_int_vector() ) 
    {
      const int sz = rhs.size();
      if ( sz == 0 ) return Token();
      std::vector<bool> ans( sz );      
      if      ( is_int() )   for (int i=0; i<sz; i++) ans[i] = ival      > rhs.ivec[i];
      else if ( is_float() ) for (int i=0; i<sz; i++) ans[i] =      fval > rhs.ivec[i];
      return Token( ans );            
    }

  else if ( is_float_vector() )
    {
      const int sz = size();
      if ( sz == 0 ) return Token();
      std::vector<bool> ans( sz );      
      if ( rhs.is_int() )        for (int i=0; i<sz; i++) ans[i] = fvec[i] > rhs.ival;
      else if ( rhs.is_float() ) for (int i=0; i<sz; i++) ans[i] = fvec[i] > rhs.fval;
      return Token( ans );                  
    }

  else if ( rhs.is_float_vector() ) 
    {
      const int sz = rhs.size();
      if ( sz == 0 ) return Token();
      std::vector<bool> ans( sz );      
      if      ( is_int() )   for (int i=0; i<sz; i++) ans[i] = ival         > rhs.fvec[i];
      else if ( is_float() ) for (int i=0; i<sz; i++) ans[i] = fval         > rhs.fvec[i];
      return Token( ans );            
    }

  else if ( is_string_vector() && rhs.is_string() ) 
    {
      const int sz = size();
      if ( sz == 0 ) return Token();
      std::vector<bool> ans( sz );      
      for (int i=0; i<sz; i++) ans[i] = svec[i] > rhs.sval;
      return Token( ans );            
    }
  
  else if ( is_string() && rhs.is_string_vector() ) 
    {
      const int sz = size();
      if ( sz == 0 ) return Token();
      std::vector<bool> ans( sz );      
      for (int i=0; i<sz; i++) ans[i] = svec[i] > rhs.sval;
      return Token( ans );            
    }



  if ( is_int() ) 
    {
      if ( rhs.is_int() )  return Token( ival > rhs.ival );
      if ( rhs.is_float() ) return Token( ival > rhs.fval );
    }
  
  if ( is_float() ) 
    {
      if ( rhs.is_int() )  return Token( fval > rhs.ival );
      if ( rhs.is_float() ) return Token( fval > rhs.fval );
    }
  
  if ( is_string() ) 
    {
      if ( rhs.is_string() )  return Token( sval > rhs.sval );
    }
  
  return Token();
}

Token Token::operator>=(const Token & rhs) const
{
  return ! ( *this < rhs );
}

Token Token::operator<=(const Token & rhs) const
{
  return ! ( *this > rhs );
}


Token Token::operator&&(const Token & rhs) const
{
  // vector x vector comparison not defined
  if ( is_vector() && rhs.is_vector() ) return Token();


  // lazy evaluation of RHS
  if ( is_bool() && !bval ) return Token( false );
  if ( is_int() && !ival ) return Token( false );
  
  if ( is_bool() ) 
    {
      if ( rhs.is_bool() ) return Token( bval && rhs.bval );
      if ( rhs.is_int() ) return Token( bval && rhs.ival );
    }
  
  if ( is_int() ) 
    {
      if ( rhs.is_bool() ) return Token( ival && rhs.bval );
      if ( rhs.is_int() ) return Token( ival && rhs.ival );
    }

  return Token();
}


Token Token::operator||(const Token & rhs) const
{

  // vector x vector comparison not defined
  if ( is_vector() && rhs.is_vector() ) return Token();
  
  // lazy evaluation of RHS
  if ( is_bool() && bval ) return Token( true );
  if ( is_int() && ival ) return Token( true );
  
  if ( is_bool() ) 
    {
      if ( rhs.is_bool() ) return Token( bval || rhs.bval );
      if ( rhs.is_int() ) return Token( bval || rhs.ival );
    }
  
  if ( is_int() )
    {
      if ( rhs.is_bool() ) return Token( ival || rhs.bval );
      if ( rhs.is_int() ) return Token( ival || rhs.ival );
    }
  return Token();    
}

Token Token::operands( Token & t)
{
  if ( ttype == NOT_OPERATOR ) return !t;
  else return Token();
}

Token Token::operands( Token & right, Token & left )
{
  switch( ttype )
    {
    case ASSIGNMENT_OPERATOR : return right;
    case ADD_OPERATOR : return left + right ;      
    case SUBTRACT_OPERATOR : return left - right ;
    case MULTIPLY_OPERATOR : return left * right ;
    case DIVIDE_OPERATOR : return left / right ;
    case MOD_OPERATOR : return left % right ;
    case AND_OPERATOR : return left && right ;
    case OR_OPERATOR : return left || right ;
    case LESS_THAN_OPERATOR : return left < right ;
    case LESS_THAN_OR_EQUAL_OPERATOR : return left <= right ;
    case GREATER_THAN_OPERATOR : return left > right ;
    case GREATER_THAN_OR_EQUAL_OPERATOR : return left >= right ;
    case EQUAL_OPERATOR : return left == right ;
    case UNEQUAL_OPERATOR : return left != right ;
    }
  return Token();
}



int Token::as_int() const
{
  switch ( ttype ) 
    {
    case INT : return ival;
    case FLOAT : return (int)fval;
    case BOOL : return bval ? 1 : 0;
    }
  if ( ttype == STRING )
    {
      int i;
      if ( Helper::from_string<int>( i , sval, std::dec ) ) 
	return i;
    }
  return 0; // UNDEF
}


double Token::as_float() const
{
  switch ( ttype ) 
    {
    case INT : return (double)ival;
    case FLOAT : return fval;
    case BOOL : return bval ? 1.0 : 0.0;
    }
  if ( ttype == STRING )
    {
      double d;
      if ( Helper::from_string<double>( d , sval, std::dec ) ) 
	return d;
    }
  return 0.0; // UNDEF   
}


std::string Token::as_string() const
{
  // unlike the above scalar as_type() functions, 
  // this automatically converts vectors to strings, comma-delimited

  if ( ttype == STRING ) return sval;

  std::stringstream ss;

  if ( ttype == INT ) ss << ival;
  else if ( ttype == FLOAT ) ss << fval;
  else if ( ttype == BOOL ) ss << ( bval ? "T" : "F" );
  else if ( ttype == STRING_VECTOR )
    {
      for (int i=0;i<svec.size();i++) 
	ss << ( i ? "," : "" ) << svec[i]; 	
    }
  else if ( ttype == INT_VECTOR ) 
    {
      for (int i=0;i<ivec.size();i++) 
	ss << ( i ? "," : "" ) << ivec[i]; 	
    }
  else if ( ttype == FLOAT_VECTOR ) 
    {
      for (int i=0;i<fvec.size();i++) 
	ss << ( i ? "," : "" ) << fvec[i]; 	
    }
  else if ( ttype == BOOL_VECTOR )
    {
      for (int i=0;i<bvec.size();i++) 
	ss << ( i ? "," : "" ) << ( bvec[i] ? "T" : "F" );
    }
  else 
    ss << ".";

  return ss.str();
}




bool Token::as_bool() const
{

  // unlike as_int and as_float, as_bool automatically converts boolean vectors
  // to a scalar vector, based on any(T)

  if      ( ttype == BOOL )   return bval;    
  else if ( ttype == INT )    return ival;
  else if ( ttype == FLOAT )  return fval;
  else if ( ttype == STRING ) return !( sval == "" || sval == "." ) ;

  else if ( ttype == BOOL_VECTOR )  for (int i=0;i<bvec.size(); i++) { if ( bvec[i] ) return true; } 
  else if ( ttype == INT_VECTOR )   for (int i=0;i<ivec.size(); i++) { if ( ivec[i] ) return true; }
  else if ( ttype == FLOAT_VECTOR ) for (int i=0;i<fvec.size(); i++) { if ( fvec[i] ) return true; }
  else if ( ttype == STRING_VECTOR ) 
    for (int i=0;i<svec.size(); i++) { if ( ! ( svec[i] == "." || svec[i] == "" ) ) return true; }

  return false;
}



int Token::int_element(const int i) const
{
  if ( i < 0 || i > size() ) return 0;
  if ( ttype == INT_VECTOR ) return ivec[i];
  if ( ttype == INT ) return ival;
  return 0;      
}

double Token::float_element(const int i) const
{
  if ( i < 0 || i > size() ) return 0;
  if ( ttype == FLOAT_VECTOR ) return fvec[i];
  if ( ttype == FLOAT ) return fval;
  return 0;      
}


std::string Token::string_element(const int i) const
{
  if ( i < 0 || i > size() ) return ".";
  if ( ttype == STRING_VECTOR ) return svec[i];
  if ( ttype == STRING ) return sval;
  return ".";      
}

bool Token::bool_element(const int i) const
{
  if ( i < 0 || i > size() ) return false;
  if ( ttype == BOOL_VECTOR ) return bvec[i];
  if ( ttype == BOOL ) return bval;
  return false;  
}


std::vector<int> Token::as_int_vector() const
{

  if ( ttype == INT_VECTOR ) return ivec;

  std::vector<int> ans( size() );
  
  if ( ttype == FLOAT_VECTOR ) 
    {
      for (int i=0; i<fvec.size(); i++) ans[i] = (int)fvec[i];
      return ans;
    }
  
  if ( ttype == BOOL_VECTOR ) 
    {
      for (int i=0; i<bvec.size(); i++) ans[i] = bvec[i] ? 1 : 0 ;
      return ans;
    }
  
  if ( ttype == STRING_VECTOR ) 
    {
      for (int i=0; i<svec.size(); i++) 
	if ( ! Helper::from_string<int>( ans[i] , svec[i], std::dec ) ) ans[i] = 0;	
      return ans;
    }

  switch ( ttype ) 
    {
    case INT   : ans[0] = ival; return ans;
    case FLOAT : ans[0] = (int)fval; return ans;
    case BOOL  : ans[0] = bval ? 1 : 0; return ans;
    }
  
  if ( ttype == STRING )
    {      
      if ( ! Helper::from_string<int>( ans[0] , sval, std::dec ) ) ans[0] = 0;
      return ans;
    }
  
  // should be size 0 vector
  return ans; // UNDEF
  
}

std::vector<double> Token::as_float_vector() const
{
  if ( ttype == FLOAT_VECTOR ) return fvec;

  std::vector<double> ans( size() );
  
  if ( ttype == INT_VECTOR ) 
    {
      for (int i=0; i<ivec.size(); i++) ans[i] = ivec[i];
      return ans;
    }
  
  if ( ttype == BOOL_VECTOR ) 
    {
      for (int i=0; i<bvec.size(); i++) ans[i] = bvec[i] ? 1 : 0 ;
      return ans;
    }
  
  if ( ttype == STRING_VECTOR ) 
    {
      for (int i=0; i<svec.size(); i++) 
	if ( ! Helper::from_string<double>( ans[i] , svec[i], std::dec ) ) ans[i] = 0;	
      return ans;
    }
  
  switch ( ttype ) 
    {
    case INT   : ans[0] = ival; return ans;
    case FLOAT : ans[0] = fval; return ans;
    case BOOL  : ans[0] = bval ? 1 : 0; return ans;
    }
  
  if ( ttype == STRING )
    {      
      if ( ! Helper::from_string<double>( ans[0] , sval, std::dec ) ) ans[0] = 0;
      return ans;
    }
  
  // should be size 0 vector
  return ans; // UNDEF

}


std::vector<std::string> Token::as_string_vector() const
{
  if ( ttype == STRING_VECTOR ) return svec;
  std::vector<std::string> ans;
  Helper::halt( "as_string_vector() automatic type conversion not defined" );
  return ans;
}


std::vector<bool> Token::as_bool_vector() const
{
  if ( ttype == BOOL_VECTOR ) return bvec;

  std::vector<bool> ans;
  if ( is_scalar() ) { ans.push_back( as_bool() ); return ans; }
  
  ans.resize( size() );
  if      ( ttype == INT_VECTOR )   { for (int i=0; i<ivec.size(); i++) ans[i] = ivec[i]; }
  else if ( ttype == FLOAT_VECTOR ) { for (int i=0; i<fvec.size(); i++) ans[i] = fvec[i]; }
  else if ( ttype == STRING_VECTOR ) { for (int i=0; i<svec.size(); i++) ans[i] = true; }
  
  return ans;
}



//
// Token functions
//

Token TokenFunctions::fn_assign( Token & lhs , const Token & rhs )
{
  
  // Currently, this does not allow ASSIGNs to GenMeta, or any other
  // kind of meta.  This is fine -- just a bit messy as we use the
  // Eval class for gfuncs also... okay to keep this as is, VarMeta
  // specific, but try to clean up at some point.

  if ( ! meta ) return Token();
  
  bool b;
  if ( rhs.is_bool(&b) )
    {
      
      MetaInformation<VarMeta>::field( lhs.name() , META_BOOL , 1 , "" );      
      meta->set( lhs.name() , b );
      lhs.set( b );
      return Token( true );
    }

  int i;
  if ( rhs.is_int(&i) )
    {      
      MetaInformation<VarMeta>::field( lhs.name() , META_INT , 1 , "" );
      meta->set( lhs.name() , i );
      lhs.set( i );
      return Token( true );
    }

  double f;
  if ( rhs.is_float(&f) ) 
    {
      MetaInformation<VarMeta>::field( lhs.name() , META_FLOAT , 1 , "" );
      meta->set( lhs.name() , f );
      lhs.set( f );
      return Token( true );
    }

  std::string s;
  if ( rhs.is_string(&s) )
    {
      MetaInformation<VarMeta>::field( lhs.name() , META_TEXT , 1 , "" );
      meta->set( lhs.name() , s );
      lhs.set( s );
      return Token( true );
    }


  // now the same for vectors

  std::vector<double> fv;
  if ( rhs.is_float_vector(&fv) ) 
    {
      MetaInformation<VarMeta>::field( lhs.name() , META_FLOAT , -1 , "" );
      meta->set( lhs.name() , fv );
      lhs.set( fv );
      return Token( true );
    }

  std::vector<bool> bv;
  if ( rhs.is_bool_vector(&bv) )
    {      
      MetaInformation<VarMeta>::field( lhs.name() , META_BOOL , -1 , "" );      
      meta->set( lhs.name() , bv );
      lhs.set( bv );
      return Token( true );
    }

  std::vector<int> iv;
  if ( rhs.is_int_vector(&iv) )
    {                  
      MetaInformation<VarMeta>::field( lhs.name() , META_INT , -1 , "" );
      meta->set( lhs.name() , iv );
      lhs.set( iv );
      return Token( true );
    }


  std::vector<std::string> sv;
  if ( rhs.is_string_vector(&sv) )
    {
      MetaInformation<VarMeta>::field( lhs.name() , META_TEXT , -1 , "" );
      meta->set( lhs.name() , sv );
      lhs.set( sv );
      return Token( true );
    }
    
  return Token( true );  
} 

Token TokenFunctions::fn_set( const Token & tok ) const
{
    return tok.is_set();
} 

Token TokenFunctions::fn_sqrt( const Token & tok ) const
{
  if ( tok.is_int() ) return Token( sqrt( tok.as_int() ) );
  if ( tok.is_float() ) return Token( sqrt( tok.as_float() ) );
  if ( tok.is_int_vector() || tok.is_float_vector() ) 
    {
      std::vector<double> ans = tok.as_float_vector();
      for (int i=0; i<ans.size(); i++) ans[i] = sqrt( ans[i] );
      return Token( ans );
    }
  return Token();
}

Token TokenFunctions::fn_pow( const Token & tok , const Token & tok2 ) const
{
  bool is_int_type = tok.is_int() || tok.is_int_vector();
  bool is_float_type = tok.is_float() || tok.is_float_vector();
  if ( ! ( is_int_type || is_float_type ) ) return Token();
  if ( ! ( tok2.is_int() || tok2.is_float() ) ) return Token();
  
  if ( is_int_type && tok2.is_int() ) 
    {
      if ( tok.is_scalar() ) return Token( (int)(pow( tok.as_int() , tok2.as_int() )) );
      std::vector<int> ans = tok.as_int_vector();
      const int ev = tok2.as_int();
      for (int i=0; i<ans.size(); i++) ans[i] = (int)(pow( ans[i] , ev ) );
      return Token( ans );
    }
  
  double ev = tok2.as_float();

  if ( tok.is_int() || tok.is_float() ) return Token( pow( tok.as_float() , ev ) );
  else if ( tok.is_int_vector() || tok.is_float_vector() )
    {
      std::vector<double> ans = tok.as_float_vector();
      for (int i=0; i<ans.size(); i++) ans[i] = pow( ans[i] , ev ) ;
      return Token( ans );
    }

  return Token();
}


// CHECK -- *** We can remove the 3 functions below *** 

Token TokenFunctions::fn_n() 
{
  return 0;//Token( var->size() );
}


Token TokenFunctions::fn_g( const Token & cond )
{

  // hmm -- these will be harder, if we are to allow arbitrary
  // expressions that apply per-genotype...  leave for now..
  return Token();

}


Token TokenFunctions::fn_gmean( const Token & field , const Token & cond )
{
  // which genotype meta-field do we want to loot at?
  std::string key;
  if ( ! field.is_string( &key ) ) return Token();

  return Token();  
}


Token TokenFunctions::fn_ifelse( const Token & cond , const Token & opt1 , const Token & opt2 ) const
{
  
  // cond ? opt1 : opt2

  // condition must evaluate to a scalar boolean, or be easily converted
  // no automatic conversion of vector conditions
  
  bool b;
  
  if ( ! cond.is_bool(&b) ) 
    {
      // allow implicit conversion of int -> bool
      if ( cond.is_int() ) b = cond.as_bool();
      else return Token();
    }  

  // opt1 and opt2 must have the same type, or be easily converted
  
  if ( opt1.type() == opt2.type() ) 
    {
      return b ? opt1 : opt2 ;
    }    
  
  // If opt1 and opt2 are not the same, attempt to convert in some cases

  // Safe automatic upcasts
  // F I --> F F
  // F B --> F F
  // I B --> I I (0/1)
  // do not allow any STRINGs to be converted
  
  Token tmp1 = opt1;
  Token tmp2 = opt2;

  Token::tok_type t1 = tmp1.type();
  Token::tok_type t2 = tmp2.type();

  if ( t1 == Token::UNDEF || t2 == Token::UNDEF ) return Token();
  
  if ( t1 == Token::STRING || t2 == Token::STRING ) 
    Helper::halt("ifelse(?,T,F) cannot specify incompatible return types");
  
  if      ( t1 == Token::FLOAT )  tmp2 = Token( tmp2.as_float() );
  else if ( t2 == Token::FLOAT )  tmp1 = Token( tmp1.as_float() );
  else if ( t1 == Token::INT   )  tmp2 = Token( tmp2.as_int()   ); 
  else if ( t2 == Token::INT   )  tmp1 = Token( tmp1.as_int()   );  
  else return Token();  // all other comparisons not defined

  // Now evaluate the original expression
  return b ? tmp1 : tmp2 ;
  
}


Token TokenFunctions::fn_vec_length( const Token & tok) const 
{
  return Token( tok.size() ) ;
}


Token TokenFunctions::fn_vec_any( const Token & tok1 , const Token & tok2 ) const
{
  return Token( fn_vec_count( tok1 , tok2 ) > 0 );  
}


Token TokenFunctions::fn_vec_count( const Token & tok1 , const Token & tok2 ) const
{
  plog.warn( "fn_vec_count() not yet implemented" );
  return Token();  
}


Token TokenFunctions::fn_vec_sum( const Token & tok ) const
{
  
  Token::tok_type ttype = tok.type();
  
  if ( tok.is_scalar() ) return tok;
  
  if ( ttype == Token::INT_VECTOR ) 
    {
      int sm = 0;
      std::vector<int> tmp = tok.as_int_vector();
      for (int i=0;i<tmp.size(); i++) sm += tmp[i];      
      return Token( sm );
    }
  
  if ( ttype == Token::FLOAT_VECTOR ) 
    {
      double sm = 0;
      std::vector<double> tmp = tok.as_float_vector();
      for (int i=0;i<tmp.size(); i++) sm += tmp[i];
      return Token( sm );
    }

  // convert to INT
  if ( ttype == Token::BOOL_VECTOR )
    {
      int sm = 0;
      std::vector<bool> tmp = tok.as_bool_vector();
      for (int i=0;i<tmp.size(); i++) sm += tmp[i];
      return Token( sm );
    }

  // undefined
  return Token();
       
}

Token TokenFunctions::fn_vec_mean( const Token & tok1 ) const
{
  return fn_vec_sum( tok1 ) / fn_vec_length( tok1 ); 
}


Token TokenFunctions::fn_vec_sort( const Token & tok ) const
{

  if ( ! tok.is_vector() ) return tok;  

  Token::tok_type ttype = tok.type();        

  if ( ttype == Token::INT_VECTOR ) 
    {
      std::vector<int> t = tok.as_int_vector();
      std::sort( t.begin() , t.end() );
      return Token( t );
    }
  else if ( ttype == Token::FLOAT_VECTOR ) 
    {
      std::vector<double> t = tok.as_float_vector();
      std::sort( t.begin() , t.end() );
      return Token( t );
    }
  else if  ( ttype == Token::STRING_VECTOR ) 
    {
      std::vector<std::string> t = tok.as_string_vector();
      std::sort( t.begin() , t.end() );
      return Token( t );  
    }
  else if ( ttype == Token::BOOL_VECTOR ) 
    {
      std::vector<bool> t = tok.as_bool_vector();
      std::sort( t.begin() , t.end() );
      return Token( t );
    }
  return Token();
}


Token TokenFunctions::fn_vec_extract( const Token & tok , const Token & idx ) const
{

  if ( ! ( idx.is_int() || idx.is_int_vector() || idx.is_bool_vector() ) ) 
    Helper::halt( "index for vector subscripting is not an integer value, integer vector or boolean vector" );
  
  // extract out a single element and return as a scalar
  
  if ( idx.is_int() )
    {
      
      int i = idx.as_int();
      
      // subscript out of range? (1-based)
      if ( i < 1 || i > tok.size() ) return Token();
      
      // if scalar, return whole thing (i.e. could only have been x[0] so pointless)
      if ( ! tok.is_vector() ) return tok;
      
      Token::tok_type ttype = tok.type();      
      
      if      ( ttype == Token::INT_VECTOR )    return Token( tok.int_element(i-1) ); 
      else if ( ttype == Token::FLOAT_VECTOR )  return Token( tok.float_element(i-1) );
      else if ( ttype == Token::STRING_VECTOR ) return Token( tok.string_element(i-1) );
      else if ( ttype == Token::BOOL_VECTOR )   return Token( tok.bool_element(i-1) );
      
      return Token();

    }

  // extract out a list 
  
  if ( idx.is_int_vector() )
    {
      
      Token::tok_type ttype = tok.type();

      if ( ttype == Token::INT_VECTOR ) 
	{
	  std::vector<int> ans;
	  for (int i=0;i<idx.size();i++) ans.push_back( tok.int_element( idx.int_element(i) - 1 ) );
	  return Token(ans);
	}
      else if ( ttype == Token::FLOAT_VECTOR )
	{
	  std::vector<double> ans;
	  for (int i=0;i<idx.size();i++) ans.push_back( tok.float_element( idx.int_element(i) - 1 ) );
	  return Token(ans);
	}
      else if ( ttype == Token::STRING_VECTOR )
	{
	  std::vector<std::string> ans;
	  for (int i=0;i<idx.size();i++) ans.push_back( tok.string_element( idx.int_element(i) - 1 ) );
	  return Token(ans);
	}
      else if ( ttype == Token::BOOL_VECTOR )
	{
	  std::vector<bool> ans;
	  for (int i=0;i<idx.size();i++) ans.push_back( tok.bool_element( idx.int_element(i) - 1 ) );
	  return Token(ans);
	}
      return Token();
    }

  
  // extract from a boolean vector
  
  if ( idx.is_bool_vector() )
    {
      // here require that length of index matches (no cycling, as per R)
      if ( idx.size() != tok.size() ) return Token();

      Token::tok_type ttype = tok.type();
      
      if ( ttype == Token::INT_VECTOR ) 
	{
	  std::vector<int> ans;
	  for (int i=0;i<idx.size();i++) if ( idx.bool_element(i-1) ) ans.push_back( tok.int_element(i-1) );
	  return Token(ans);
	}
      else if ( ttype == Token::FLOAT_VECTOR )
	{
	  std::vector<double> ans;
	  for (int i=0;i<idx.size();i++) if ( idx.bool_element(i-1) ) ans.push_back( tok.float_element(i-1) );
	  return Token(ans);
	}
      else if ( ttype == Token::STRING_VECTOR )
	{
	  std::vector<std::string> ans;
	  for (int i=0;i<idx.size();i++) if ( idx.bool_element(i-1) ) ans.push_back( tok.string_element(i-1) );
	  return Token(ans);
	}
      else if ( ttype == Token::BOOL_VECTOR )
	{
	  std::vector<bool> ans;
	  for (int i=0;i<idx.size();i++) if ( idx.bool_element(i-1) ) ans.push_back( tok.bool_element(i-1) );
	  return Token(ans);
	}
      return Token();
    
    }

}


Token TokenFunctions::fn_vec_min( const Token & tok ) const
{
  if ( ! tok.is_vector() ) return tok;
  Token s = fn_vec_sort(tok);
  Token::tok_type ttype = tok.type();  
  if      ( ttype == Token::INT_VECTOR ) return Token( s.int_element(0) );
  else if ( ttype == Token::FLOAT_VECTOR ) return Token( s.float_element(0) );
  else if ( ttype == Token::BOOL_VECTOR ) return Token( s.bool_element(0) );
  else if ( ttype == Token::STRING_VECTOR ) return Token( s.string_element(0) );
  return Token();
}


Token TokenFunctions::fn_vec_maj( const Token & tok ) const
{
  if ( ! tok.is_vector() ) return tok;
  Token s = fn_vec_sort(tok);
  Token::tok_type ttype = tok.type();  
  const int sz = tok.size() - 1;
  if      ( ttype == Token::INT_VECTOR ) return Token( s.int_element(sz) );
  else if ( ttype == Token::FLOAT_VECTOR ) return Token( s.float_element(sz) );
  else if ( ttype == Token::BOOL_VECTOR ) return Token( s.bool_element(sz) );
  else if ( ttype == Token::STRING_VECTOR ) return Token( s.string_element(sz) );
  return Token();
}


Token TokenFunctions::fn_vec_new_float( const Token & tok ) const
{
  // create a new float vector from a comma-delimited string literal, '1,2,3'
  if ( tok.type() != Token::STRING ) return Token();
  std::vector<std::string> t = Helper::char_split( tok.as_string() , ',' , false );
  std::vector<double> d( t.size() );
  for (int i=0; i<d.size(); i++)
    if ( ! Helper::from_string<double>( d[i] , t[i] , std::dec ) ) return Token();
  return Token( d );
}


Token TokenFunctions::fn_vec_new_int( const Token & tok ) const
{
  // create a new float vector from a comma-delimited string literal, '1,2,3'
  if ( tok.type() != Token::STRING ) return Token();
  std::vector<std::string> t = Helper::char_split( tok.as_string() , ',' , false );
  std::vector<int> d( t.size() );
  for (int i=0; i<d.size(); i++)
    if ( ! Helper::from_string<int>( d[i] , t[i] , std::dec ) ) return Token();
  return Token( d );
}

Token TokenFunctions::fn_vec_new_str( const Token & tok ) const
{
  // create a new float vector from a comma-delimited string literal, '1,2,3'
  if ( tok.type() != Token::STRING ) return Token();
  return Token( Helper::char_split( tok.as_string() , ',' , false ) );
}


Token TokenFunctions::fn_vec_new_bool( const Token & tok ) const
{
  // create a new float vector from a comma-delimited string literal, '1,2,3'
  if ( tok.type() != Token::STRING ) return Token();
  std::vector<std::string> t = Helper::char_split( tok.as_string() , ',' , false );
  std::vector<bool> d( t.size() );
  for (int i=0; i<d.size(); i++)
    { 
      if( t[i] == "T" || t[i] == "true" || t[i] == "1" ) d[i] = true;
      else if ( t[i] == "F" || t[i] == "false" || t[i] == "0" ) d[i] = false;
      else return Token();
    }
   return Token( d );
}

