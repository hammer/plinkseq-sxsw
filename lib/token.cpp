
#include "token.h"

#include <sstream>
#include <cmath>

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

    tok_map[ "*" ] = MULTIPLY_OPERATOR;
    tok_map[ "^" ] = POWER_OPERATOR;
    tok_map[ "/" ] = DIVIDE_OPERATOR;
    tok_map[ "%" ] = MOD_OPERATOR;
    tok_map[ "%%" ] = MOD_OPERATOR;
    tok_map[ "+" ] = ADD_OPERATOR;
    tok_map[ "-" ] = SUBTRACT_OPERATOR;
    tok_map[ "&&" ] = AND_OPERATOR;
    tok_map[ "&" ] = AND_OPERATOR;
    tok_map[ "||" ] = OR_OPERATOR;
    tok_map[ "|" ] = OR_OPERATOR;
    tok_map[ "=" ] = ASSIGNMENT_OPERATOR;
    tok_map[ "==" ] = EQUAL_OPERATOR;
    tok_map[ "!=" ] = UNEQUAL_OPERATOR;
    tok_map[ "!" ] = NOT_OPERATOR;
    tok_map[ ">" ] = GREATER_THAN_OPERATOR;
    tok_map[ ">=" ] = GREATER_THAN_OR_EQUAL_OPERATOR;
    tok_map[ "<" ] = LESS_THAN_OPERATOR;
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

    fn_map[ "set"  ] = 1; // number of args
    fn_map[ "sqrt" ] = 1; // square-root
    fn_map[ "sqr"  ] = 1; // X^2
    fn_map[ "pow"  ] = 2; // X^N    
    fn_map[ "ifelse" ] = 3;  // ifelse( cond , T , F )
    fn_map[ "n" ] = 0;    // number of people in file
    
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

Token::Token( const Token & rhs )
{
    *this = rhs;
}

Token Token::operator!()
{
  if ( is_bool() ) return Token( ! bval ); 
  else if ( is_int() ) return Token( ival == 0 );
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
    return *this;
}

Token Token::operator!=(const Token & rhs )
{
    if ( is_bool() && rhs.is_bool() ) return Token( bval != rhs.bval ); 
    if ( is_int() && rhs.is_int() ) return Token( ival != rhs.ival ); 
    if ( is_float() && rhs.is_float() ) return Token( fval != rhs.fval ); 
    if ( is_string() && rhs.is_string() ) return Token( sval != rhs.sval ); 
    if ( is_int() && rhs.is_bool() ) return Token( ival != rhs.bval );
    if ( is_bool() && rhs.is_int() ) return Token( bval != rhs.ival );
    if ( is_float() && rhs.is_int() ) return Token( fval != rhs.ival );
    if ( is_int() && rhs.is_float() ) return Token( ival != rhs.fval );
    return Token();
}

Token Token::operator==(const Token & rhs)
{
  if ( is_bool() && rhs.is_bool() ) return Token( bval == rhs.bval ); 
  if ( is_int() && rhs.is_int() ) return Token( ival == rhs.ival ); 
  if ( is_float() && rhs.is_float() ) return Token( fval == rhs.fval ); 
  if ( is_string() && rhs.is_string() ) return Token( sval == rhs.sval ); 
  if ( is_int() && rhs.is_bool() ) return Token( ival == rhs.bval );
  if ( is_bool() && rhs.is_int() ) return Token( bval == rhs.ival );
  if ( is_float() && rhs.is_int() ) return Token( fval == rhs.ival );
  if ( is_int() && rhs.is_float() ) return Token( ival == rhs.fval );
  return Token();
}

Token Token::operator+(const Token & rhs)
{
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

Token Token::operator-(const Token & rhs)
{
    if ( is_int() ) 
    {
	if ( rhs.is_int() ) return Token( ival - rhs.ival );
	if ( rhs.is_float() ) return Token( ival - rhs.fval );
    }
    if ( is_float() ) 
    {
	if ( rhs.is_int() ) return Token( fval - rhs.ival );
	if ( rhs.is_float() ) return Token( fval - rhs.fval );
    }
    return Token();
}

Token Token::operator*(const Token & rhs)
{
    if ( is_int() ) 
    {
	if ( rhs.is_int() ) return Token( ival * rhs.ival );
	if ( rhs.is_float() ) return Token( ival * rhs.fval );
    }
    if ( is_float() ) 
    {
	if ( rhs.is_int() ) return Token( fval * rhs.ival );
	if ( rhs.is_float() ) return Token( fval * rhs.fval );
    }
    return Token();
}

Token Token::operator^(const Token & rhs)
{
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

Token Token::operator/(const Token & rhs)
{
    if ( is_int() ) 
    {
	if ( rhs.is_int() )  return Token( (double)ival / (double)rhs.ival );
	if ( rhs.is_float() ) return Token( (double)ival / rhs.fval );
    }
    if ( is_float() ) 
    {
	if ( rhs.is_int() ) return Token( fval / (double)rhs.ival );
	if ( rhs.is_float() ) return Token( fval / rhs.fval );
    }
    return Token();
    return Token();
}

Token Token::operator%(const Token & rhs)
{
    if ( is_int() ) 
    {
	if ( rhs.is_int() ) return Token( (int)(ival % rhs.ival) );
    }
    return Token();
}


Token Token::operator<(const Token & rhs)
{
    if ( is_int() ) 
    {
	if ( rhs.is_int() )  return Token( ival < rhs.ival );
	if ( rhs.is_float() ) return Token( ival < rhs.fval );
    }

    if ( is_float() ) 
    {
	if ( rhs.is_int() )  return Token( fval < rhs.ival );
	if ( rhs.is_float() ) return Token( fval < rhs.fval );
    }

    if ( is_string() ) 
    {
	if ( rhs.is_string() )  return Token( sval < rhs.sval );
    }

    return Token();
}

Token Token::operator>(const Token & rhs)
{
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

Token Token::operator>=(const Token & rhs)
{
    if ( is_int() ) 
    {
	if ( rhs.is_int() )  return Token( ival >= rhs.ival );
	if ( rhs.is_float() ) return Token( ival >= rhs.fval );
    }

    if ( is_float() ) 
    {
	if ( rhs.is_int() )  return Token( fval >= rhs.ival );
	if ( rhs.is_float() ) return Token( fval >= rhs.fval );
    }

    if ( is_string() ) 
    {
	if ( rhs.is_string() )  return Token( sval >= rhs.sval );
    }

    return Token();
}

Token Token::operator<=(const Token & rhs)
{
    if ( is_int() ) 
    {
	if ( rhs.is_int() )  return Token( ival <= rhs.ival );
	if ( rhs.is_float() ) return Token( ival <= rhs.fval );
    }

    if ( is_float() ) 
    {
	if ( rhs.is_int() )  return Token( fval <= rhs.ival );
	if ( rhs.is_float() ) return Token( fval <= rhs.fval );
    }

    if ( is_string() ) 
    {
	if ( rhs.is_string() )  return Token( sval <= rhs.sval );
    }

    return Token();
}

Token Token::operator&&(const Token & rhs)
{
  // lazy evaluation 
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

Token Token::operator||(const Token & rhs)
{
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
    switch ( ttype ) {
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
    switch ( ttype ) {
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
  if ( ttype == STRING ) return sval;
  std::stringstream ss;
  if ( ttype == INT ) ss << ival;
  else if ( ttype == FLOAT ) ss << fval;
  else if ( ttype == BOOL ) ss << ( bval ? "T" : "F" );
  return ss.str();
}


bool Token::as_bool() const
{
    if ( ttype == BOOL ) return bval;    
    if ( ttype == INT ) return ival;
    if ( ttype == FLOAT ) return fval;
    if ( ttype == STRING ) return !( sval == "" || sval == "." ) ;
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

  if ( !meta ) return Token();
  
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
      registerMetatype( lhs.name() , META_INT, 1 , META_GROUP_VAR , "" );
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
  return Token( true );  
} 

Token TokenFunctions::fn_set( const Token & tok )
{
    return tok.is_set();
} 

Token TokenFunctions::fn_sqrt( const Token & tok )    
{
    if ( tok.is_int() ) return Token( sqrt( tok.as_int() ) );
    if ( tok.is_float() ) return Token( sqrt( tok.as_float() ) );
    return Token();
}

Token TokenFunctions::fn_pow( const Token & tok , const Token & tok2 )
{
    if ( tok.is_int() && tok2.is_int() ) return Token( (int)(pow( tok.as_int() , tok2.as_int() )) );
    if ( ( tok.is_int() || tok.is_float() ) && ( tok2.is_int() || tok2.is_float() ) )
	return Token( pow( tok.as_float() , tok2.as_float() ) );
    return Token();
}


// *** We can remove the 3 functions below *** 

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


Token TokenFunctions::fn_ifelse( const Token & cond , const Token & opt1 , const Token & opt2 )
{

  // cond ? opt1 : opt2

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
  
  // Safe automatic upcasts
  // F I --> F F
  // F B --> F F
  // I B --> I I (0/1)
  // do not allow any STRINGs to be converted
      
  return b ? tmp1 : tmp2 ;
  
}
