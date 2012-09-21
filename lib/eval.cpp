#include <string>
#include <cstdio>
#include <iostream>
#include <vector>
#include <cmath>
#include <algorithm>

#include "plinkseq/eval.h"
#include "plinkseq/token.h"
#include "plinkseq/meta.h"

// initialise static members
std::map<std::string,Token::tok_type> Token::tok_map; 
std::map<Token::tok_type,std::string> Token::tok_unmap; 
std::map<std::string,int> Token::fn_map; 


void Eval::init()
{
  is_valid = false;
  
  errs = "";

  genmeta_mode = false;


  // G-functions (# of args, including expr) 

  // return vector evaluated on genotypes
  gdef[ "g" ] = 1;      // expression

  // set genotype to null of expression evaluates to false
  gdef[ "gf" ]   = 1;

  // set genotype to VCF-formatted value (or # of alternate alleles)
  gdef[ "gs" ]   = 1;

  // phenotype extraction/assignment
  gdef[ "p" ] = 1;      // expression

  // vector creation functions
//   gdef[ "int" ] = 1;
//   gdef[ "bool" ] = 1;
//   gdef[ "vec" ] = 1;
//   gdef[ "str" ] = 1;


}

bool Eval::get_token( std::string & input ,  Token & tok )
{

  // no more input 

  if ( input.size() == 0 ) return false;

  std::string c = input.substr(0,1);
  
  // To interpret whether a + or - sign, or '.' is the start of a number 
  // rather than a binary operator
  // e.g 
  // (C - 2)  vs ( C < -2 )
  
  // eat leading space
  
  while ( 1 ) 
    {
      if ( c == " " ) input = input.substr(1);
      else break;
      if ( input.size() == 0 ) return false;	    
      c = input.substr(0,1);
    }
  


  // token is either func()
  //                 number
  //                 variable
  //                 operator 
  


  std::map<std::string,Token::tok_type>::iterator i = 
    Token::tok_map.find( c ) ;

  // Numeric, allow to start with + or - or .
  
  if ( ( c >= "0" && c <= "9" ) || c == "." || ( (!previous_value) && ( c=="-" || c=="+" ) ) )
    {
      int p = 1;

      bool need_leading_zero = c == ".";
      
      while ( 1 ) 
	{
	  
	  if ( p == input.size() ) break;
	  char d = input[p];
	  if ( d >= '0' && d <= '9' ) c += input.substr(p,1);
	  else if ( d == '.' ) c += input.substr(p,1); 
	  else if ( d == 'e' || d == 'E' ) c += input.substr(p,1);
	  else if ( d == '+' || d == '-' ) {
	    // + and - okay, if previous char was 'e' or 'E'
	    if ( input.substr(p-1,1) == "E" || input.substr(p-1,1) == "e" )
	      c += input.substr(p,1);
	    else break;
	  }
	  else break;
	  ++p;
	}

      // a valid int or float? 
      int i;
      double d;
      bool err = false;

      if ( ! Helper::from_string<int>( i , need_leading_zero ? "0"+c : c , std::dec ) ) err = true;
      if ( ! Helper::from_string<double>( d , need_leading_zero ? "0"+c : c , std::dec ) ) err = true;

      if ( ! err ) 
	{
	  if ( i == d ) tok.set(i);
	  else tok.set(d);
	  previous_value = true;
	}
      
    }
    
  else if ( i != Token::tok_map.end() )
    {

      // looks like an operator
      // if potentially a two-char operator, look for 2nd char
      
      if ( c == "%" && input.substr(1,1) == "%" ) c = "%%";
      if ( c == "<" && input.substr(1,1) == "=" ) c = "<=";
      if ( c == ">" && input.substr(1,1) == "=" ) c = ">=";
      if ( c == "&" && input.substr(1,1) == "&" ) c = "&&";
      if ( c == "|" && input.substr(1,1) == "|" ) c = "||";
      if ( c == "!" && input.substr(1,1) == "=" ) c = "!=";
      if ( c == "=" && input.substr(1,1) == "=" ) c = "==";
      
      tok.oper( Token::tok_map[ c ] ); 
      
      previous_value = false;

      // eat input name 
      input = input.substr( c.size() );
      return true;
      
    }
  

  // If not an operator or a function, must be 
  // a variable (meta-field of a literal)
  // literal strings should be in 'single quotes'
  
  
  // Handle some other special case syntax first: ( ) and ,
  
  if ( c == "(" ) 
    {
      tok.oper( Token::LEFT_PARENTHESIS );
      previous_value = false;
    }
  
  else if ( c == ")" )
    {
      tok.oper( Token::RIGHT_PARENTHESIS );
      previous_value = true;
    }
  
  else if ( c == "," )
    {
      tok.oper( Token::ARG_SEPARATOR );
      previous_value = false;      
    }

  // a literal string
  else if ( c == "'" )
    {
      int p = 1;
      bool found = false;
      while ( 1 ) 
	{
	  if ( p == input.size() ) break;
	  if ( input.substr(p,1) == "'" ) { c+= "'"; found = true; break; }
	  c += input.substr(p,1);
	  ++p;
	}
      
      tok.set( c.substr(1,c.size()-2) ); 
      previous_value = true;
    }
  
  // a literal string with specific start/stop quotes, allowing nesting
  else if ( c == "{" )
    {
      int cnt = 1;
      int p = 1;
      bool found = false;
      while ( 1 ) 
	{
	  
	  if      ( p == input.size() ) break;
	  else if ( input.substr(p,1) == "{" ) ++cnt;
	  else if ( input.substr(p,1) == "}" ) 
	    { 
	      --cnt;
	      if ( cnt == 0 ) 
		{
		  c+= "}"; 
		  found = true; break; 
		}
	    }
	  
	  c += input.substr(p,1);
	  ++p;
	}      

      tok.set( c.substr(1,c.size()-2) ); 

      previous_value = true;
    }

  // otherwise, assume a variable that start: a-z, A-Z or _
  // or a variable name is ends in "("
  
  else if ( ( c >= "a" && c <= "z" ) || 
	    ( c >= "A" && c <= "Z" ) ||
	    c == "_" ) 
    {
      
      // Is this a function?
      bool isfn = false;

      // read until space or next operator char, or comma

      int p = 1;
      while ( 1 ) 
	{
	  if ( p == input.size() ) break;
	  
	  if ( input.substr(p,1) == "(" ) 
	    { isfn=true; break; }
	  
	  // e.g. end of set(AB) for 'AB'	  
	  if ( input.substr(p,1) == ")" ) break; 
	  
	  if ( input.substr(p,1) == "," ) break;

	  if ( input.substr(p,1) == " " ) break;

	  if ( Token::tok_map.find( input.substr(p,1) ) 
	       != Token::tok_map.end() ) 
	    break;
	  
	  c += input.substr(p,1);
	  ++p;
	}
      
      if ( isfn )
	{
	  
	  // does this look like a valid function name?
	  std::map<std::string,int>::iterator f = 
	    Token::fn_map.find( c );
	  	  
	  if ( f == Token::fn_map.end() ) 
	    return false;
	  
	  // store function token
	  tok.function( c );	    
	}
      else if ( c == "T" || c == "true" )
	{
	  tok.set( true );
	}
      else if ( c == "F" || c == "false" )
	{
	  tok.set( false );
	}
      else // add token as variable
	{	  
	  tok.variable( c );
	}
      previous_value = true;
    }
  
  // remove this input

  input = input.substr( c.size() );

  return true;
  
}


int Eval::op_preced( const Token & c )
{
  
  Token::tok_type t = c.type();
  
  switch(t) 
    {
    case Token::NOT_OPERATOR :  return 9;
      
    case Token::MULTIPLY_OPERATOR :  
    case Token::DIVIDE_OPERATOR :  
    case Token::MOD_OPERATOR :  return 8;
      
    case Token::ADD_OPERATOR :
    case Token::SUBTRACT_OPERATOR : return 7;
      
    case Token::LESS_THAN_OPERATOR : 
    case Token::LESS_THAN_OR_EQUAL_OPERATOR : 
    case Token::GREATER_THAN_OPERATOR : 
    case Token::GREATER_THAN_OR_EQUAL_OPERATOR : return 6;
      
    case Token::EQUAL_OPERATOR : 
    case Token::UNEQUAL_OPERATOR : return 5;

    case Token::AND_OPERATOR : return 4;

    case Token::OR_OPERATOR : return 3;
      
    case Token::LEFT_PARENTHESIS : 
    case Token::RIGHT_PARENTHESIS : return 2;
      
    case Token::ASSIGNMENT_OPERATOR : return 1;
      
    case Token::ARG_SEPARATOR : return 0;
    }
  return 0;
}

bool Eval::op_left_assoc(const Token & tok )
{
  
  Token::tok_type t = tok.type();
  
  switch(t)    
    {
    case Token::MULTIPLY_OPERATOR :
    case Token::DIVIDE_OPERATOR :
    case Token::MOD_OPERATOR :
    case Token::ADD_OPERATOR :
    case Token::EQUAL_OPERATOR :
    case Token::AND_OPERATOR :
    case Token::OR_OPERATOR :
    case Token::ARG_SEPARATOR : 
    case Token::SUBTRACT_OPERATOR : return true;
      
    case Token::ASSIGNMENT_OPERATOR :
    case Token::NOT_OPERATOR : return false;
    }
  return false;
}

unsigned int Eval::op_arg_count( const Token & tok )
{

  Token::tok_type t = tok.type();
  
  switch( t )  
    {	      
    case Token::NOT_OPERATOR :  return 1;
      
    case Token::ASSIGNMENT_OPERATOR :       
    case Token::MULTIPLY_OPERATOR :  
    case Token::DIVIDE_OPERATOR :  
    case Token::MOD_OPERATOR : 
    case Token::ADD_OPERATOR :
    case Token::SUBTRACT_OPERATOR : 
    case Token::LESS_THAN_OPERATOR : 
    case Token::LESS_THAN_OR_EQUAL_OPERATOR : 
    case Token::GREATER_THAN_OPERATOR : 
    case Token::GREATER_THAN_OR_EQUAL_OPERATOR : 
    case Token::EQUAL_OPERATOR : 
    case Token::UNEQUAL_OPERATOR : 
    case Token::AND_OPERATOR : 
    case Token::OR_OPERATOR : return 2;
            
    case Token::FUNCTION : return Token::fn_map[ tok.name() ];
    }
  return 0;
}


bool Eval::shunting_yard( const std::string & oinput, 
			  std::vector<Token> & output )
{
  
  std::string input = oinput;

  // collect output tokens here
  output.resize( input.size() );
  output.clear();
  
  Token sc;
  
  std::vector<Token> stack;
  
  // redundant, but use 'sl' to keep track of stack size also
  
  unsigned int sl = 0;
  
  previous_value = false;

  // Parse input string token by token
  
  while( 1 )
    {
      Token c;
      
      // get next token -- if end of input break
      
      if ( ! get_token( input , c ) ) break;
      
      // If the token is a number (identifier), then add it to output queue.
      
      if( c.is_ident() )  
	{
	  output.push_back(c);
	}
      
      // If the token is a function token, then push it onto the stack.
      
      else if ( c.is_function() )   
	{
	  stack.push_back( c );
	  ++sl;
	}
	
      // If the token is a function argument separator (e.g., a comma):
      
      else if ( c.is_separator() )
	{
	  bool pe = false;
	  
	  while( sl > 0 )   	    
	    {	
	      sc = stack.back();

	      if( sc.is_left_paren() )
		{		    
		  pe = true;
		  break;
		}
	      else  
		{
		  // Until the token at the top of the stack is a left
		  // parenthesis, pop operators off the stack onto the
		  // output queue.
		  
		  output.push_back( sc );
		  stack.pop_back();
		  sl--;
		}
	    }
	  
	  // If no left parentheses are encountered, either the
	  // separator was misplaced
	  // or parentheses were mismatched.
	  
	  if( !pe )   
	    {
	      errmsg( "separator or parentheses mismatched" );
	      return false;
	    }
	}
      

      // If the token is an operator, op1, then:
      
      else if ( c.is_operator() )  
	{
	  
	  while ( sl > 0 )    
	    {
	      
	      sc = stack.back();
	      
	      // While there is an operator token, o2, at the top of
	      // the stack op1 is left-associative and its precedence
	      // is less than or equal to that of op2, or op1 is
	      // right-associative and its precedence is less than
	      // that of op2,
	      
	      if( sc.is_operator() &&
		  ( ( op_left_assoc(c) && ( op_preced(c) <= op_preced(sc) ) ) ||
		    ( !op_left_assoc(c) && ( op_preced(c) < op_preced(sc) ) )))   {
		
		// Pop o2 off the stack, onto the output queue;
		
		output.push_back( sc );
		stack.pop_back();
		sl--;
	      }
	      else
		{
		  break;
		}
	    }
	  
	  // push op1 onto the stack.
	  
	  stack.push_back( c );
	  ++sl;
	}
      
      
      // If the token is a left parenthesis, then push it onto the stack.
      
      else if ( c.is_left_paren() ) 
	{	    
	  stack.push_back( c );
	  ++sl;
	}
      
      
      // If the token is a right parenthesis:
      
	else if ( c.is_right_paren() )
	  {
	    
	    bool pe = false;
	    
	    // Until the token at the top of the stack is a left
	    // parenthesis, pop operators off the stack onto the
	    // output queue
	    
	    while( sl > 0 )     
	      {
		sc = stack.back();
		
		if( sc.is_left_paren() )
		  {
		    pe = true;
		    break;
		  }
		else  
		  {
		    output.push_back( sc );
		    stack.pop_back();
		    sl--;
		  }
	      }
	    
	    // If the stack runs out without finding a left
	    // parenthesis, then there are mismatched parentheses.
	    
	    if( ! pe )  
	      {
		errmsg( "parentheses mismatched" );
		return false;
	      }
	    
	    
	    // Pop the left parenthesis from the stack, but not onto
	    // the output queue.
	    
	    stack.pop_back();
	    sl--;
	    

	    // If the token at the top of the stack is a function
	    // token, pop it onto the output queue.

	    if( sl > 0 )   
	      {
		sc = stack.back();
		
		if( sc.is_function() )
		  {
		    output.push_back( sc );
		    stack.pop_back();
		    sl--;
		  }
	    }
	  }
	else  
	  {	    
	    errmsg( "unknown token" ) ;
	    return false; 
	  }
      
    }

  
  
  // When there are no more tokens to read: While there are still
  // operator tokens in the stack:
  
  while ( sl > 0 )  
    {
      sc = stack.back();
      
      if ( sc.is_left_paren() || sc.is_right_paren() )
	{
	  errmsg( "parentheses mismatched" );
	  return false;
        }
      
      output.push_back( sc );
      stack.pop_back();
      --sl;
    }
  
  return true;
}


bool Eval::execute( const std::vector<Token> & input )
{

  Token sc;
  
  std::vector<Token> stack;
  
  // redundant, but keep track of stack size here also
  unsigned int sl = 0;
  
  
  // While there are input tokens left
  
  for (int i = 0 ; i < input.size() ; i++ )
    {
      
      // Read the next token from input.
      
      Token c = input[i];
      
      // If the token is a value or identifier
      
      if ( c.is_ident() )
	{
	  
	  // Push it onto the stack.
	  
	  stack.push_back( c );
	  ++sl;
        }
      
      // Otherwise, the token is an operator
      // (operator here includes both operators, and functions).
	
      else if ( c.is_operator() || c.is_function() )
	{

	  // It is known a priori that the operator takes n arguments.
	  
	  int nargs = op_arg_count(c);
	  
	  // If there are fewer than n values on the stack
	  
	  if ( sl < nargs && nargs != -1 ) 
	    {
	      errmsg( "not enough arguments for " + c.name() ) ;
	      return false;
	    }
	  
	  //
	  // Else, Pop the top n values from the stack.
	  // Evaluate the operator, with the values as arguments.
	  //
	  
	  Token res;

	  if ( c.is_function() ) 
	    {
	      
	      std::vector<Token> args;
	      
	      // for fixed 'nargs' functions, just count them off 
	      
	      if ( nargs == -1 ) 
		{
		  nargs = stack.back().as_int();
		  stack.pop_back(); 
		  sl--;
		}
	      
	      while ( nargs > 0 ) 
		{		  
		  sc = stack.back();
		  stack.pop_back(); 
		  sl--;
		  args.push_back(sc);
		  --nargs;
		}
	      

	      //
	      // Perform function calls:
	      //
	      
	      // correct # of arguments? 
	      // -1 means variable length
	      
	      const int nargs = Token::fn_map[ c.name() ]; 
	      if ( args.size() != nargs && nargs != -1 ) 
		{
		  errmsg( "wrong number of arguments for " + c.name() );
		  return false;
		}

	      // note: args are in reverse order here
	      
	      if      ( c.name() == "set" )  res = func.fn_set( args[0] );
	      
	      else if ( c.name() == "sqrt" ) res = func.fn_sqrt( args[0] );
	      else if ( c.name() == "sqr" )  res = func.fn_sqr( args[0] );
	      else if ( c.name() == "pow" )  res = func.fn_pow( args[1] , args[0] );
	      
	      else if ( c.name() == "exp" ) res = func.fn_exp( args[0] );
	      else if ( c.name() == "log" ) res = func.fn_log( args[0] );
	      else if ( c.name() == "log10" ) res = func.fn_log10( args[0] );

	      else if ( c.name() == "ifelse" ) res = func.fn_ifelse( args[2], args[1], args[0] );
	      
	      // vector functions
	      else if ( c.name() == "element" ) res = func.fn_vec_extract( args[1] , args[0] );
	      
	      else if ( c.name() == "length" )  res = func.fn_vec_length( args[0] );	      
	      
	      else if ( c.name() == "min" )     res = func.fn_vec_min( args[0] );	      
	      else if ( c.name() == "max" )     res = func.fn_vec_maj( args[0] );
	      
	      else if ( c.name() == "sum" )     res = func.fn_vec_sum( args[0] );
	      else if ( c.name() == "mean" )    res = func.fn_vec_mean( args[0] );
	      else if ( c.name() == "sort" )    res = func.fn_vec_sort( args[0] );

	      else if ( c.name() == "vec_func" )  res = func.fn_vec_new_float( args );
	      else if ( c.name() == "int_func" )   res = func.fn_vec_new_int( args );
	      else if ( c.name() == "str_func" )   res = func.fn_vec_new_str( args );     
	      else if ( c.name() == "bool_func" )  res = func.fn_vec_new_bool( args ); 
	      
	      else if ( c.name() == "any" )     res = func.fn_vec_any( args[1] , args[0] );	      
	      else if ( c.name() == "count" )   res = func.fn_vec_count( args[1] , args[0] );	      

	      else if ( c.name() == "g_func" )   res = func.fn_vec_g( args[0] , this );
	      else if ( c.name() == "gf_func" )  res = func.fn_vec_gnull( args[0] , this );
	      else if ( c.name() == "gs_func" )  res = func.fn_vec_gset( args[0] , this );
	      else if ( c.name() == "n" )        res = Token( gvar ? gvar->calls.size() : 0 );

	      else if ( c.name() == "p_func" )   res = genmeta_mode ? func.fn_vec_1pheno( args[0] , indiv ) : func.fn_vec_pheno( args[0] );

	    }
	  else
	    {
	      
	      if ( nargs == 1 )  // unary operator
		{
		  sc = stack.back();
		  stack.pop_back(); sl--;		                        
		  res = c.operands( sc );		    
		}
	      else // binary operator
		{
		  Token t0 = stack.back();
		  stack.pop_back(); sl--;
		  
		  sc = stack.back();
		  stack.pop_back(); sl--;		    
		  
		  // For assignment (that impacts meta-information)
		  // need to call a function also

		  if ( c.is_assignment() ) 
		    {
		      if ( sc.name() == "p_func" )
			res = func.fn_assign_pheno( sc , t0 ); 
		      else if ( genmeta_mode )
			res = func.fn_assign_gen( sc, t0 );  // res==T
		      else
			res = func.fn_assign_var( sc, t0 );  // res==T
		      
		      // and bind new value if needed
		      bind( &sc );
		      
		    }		  
		  else
		    res = c.operands( t0 , sc );		    
		  
		}
	    }
	
	  // Push the returned results, if any, back onto the stack.
	  stack.push_back( res );
	  ++sl;
        }
      
    }

 

  // If there is only one value in the stack
  // That value is the result of the calculation.
  
  if ( sl != 1 || stack.size() != 1 ) 
    {
      errmsg( "badly formed expression" );
      return false;
    }
  
  if ( sl == 1 ) 
    {
      sc = stack.back(); 
      stack.pop_back();
      sl--;
      
      // store result in primary slot, e
      e = sc;
      
      return true;
    }
  
  
  // If there are more values in the stack
  // (Error) The user input has too many values.
  
  errmsg( "badly formed expression: too many values" );
  return false;
}

bool Eval::parse( const std::string & input )
{
  
  gvar = NULL; 

  delete_symbols();    
  
  std::string input2 = input;
  
  if ( ! expand_indices( &input2 ) ) return false;

  if ( ! expand_vargs( &input2 ) ) return false;

  // this may contain several statements, delimited by ";"
  // evaluate each sequential, to perform any assignments into
  // meta data, but only the last will be reflected in the 'e' 
  // endpoint

  
  std::vector<std::string> etok = Helper::parse( input2, ";" );

  // set number of evals we need to do

  neval = etok.size();
  
  output.resize( neval );
  
  is_valid = true;


  for (int i=0; i<etok.size(); i++)
    {      
      output[i].clear();    

      errs = "";

      if ( ! extract_gfunc( &(etok[i]) ) )
	is_valid = false;

      if ( ! shunting_yard( etok[i], output[i] ) )
	is_valid = false;
      
    }
  
  // set pointers to all variables now construction of tokens is complete
  for (int i=0; i<etok.size(); i++)  
    locate_symbols( output[i] ); 
  
  return is_valid;

}

bool Eval::extract_gfunc( std::string * s )
{

  // replace g( DP > 10 && PL[1] < 0.2 ) with 
  //  with,  g( 'DP > 10 && PL[1] < 0.2' ) 
  
  // NOTE -- not allowed to have nested g() functions
  
  // Also set p(X) -> p('X')
  
  while ( 1 ) 
    {

      bool found = false;
      
      std::map<std::string,int>::iterator i = gdef.begin();
      
      while ( i != gdef.end() )
	{
	  
	  size_t p = s->find( i->first + "(" );
	  
	  if ( p != std::string::npos ) 
	    {	      

	      // this will also match 
	      //       g(X)
	      //     log(X)
	      // so we need to check that an actual g() or p() has been found
	      
	      if ( p > 1 )
		{
		  if ( ( (*s)[p-1] >= 'A' && (*s)[p-1] <= 'Z' ) || 
		       ( (*s)[p-1] >= 'a' && (*s)[p-1] <= 'z' ) || 
		       ( (*s)[p-1] >= '0' && (*s)[p-1] <= '9' ) || 
		       (*s)[p-1] >= '_' ) { p = std::string::npos; } 
		}
	    }
	  

	  if ( p != std::string::npos )
	    {
	      
	      found = true;
	      
	      std::vector<std::string> arglist;
	      arglist.push_back( i->first ); // track fn name: g() or gnull()

	      // look for closing brace
	      
	      int bc = 0;
	      int q = p;
	      int pp = p;
	      
	      while ( ++q )
		{
		  
		  // gone past end of string?
		  if ( q == s->size() ) return false;
		  
		  char c = s->substr(q,1)[0];
		  
		  if ( c == '(' ) 
		    {
		      ++bc;  //includes first paran
		      if ( bc == 1 ) pp = q+1; // track position of first argument
		    }
		  else if ( c == ')' )
		    {
		      --bc;
		      if ( bc == 0 ) 
			{
			  std::string str = s->substr(pp,q-pp);
			  str.erase(remove_if(str.begin(), str.end(), isspace), str.end());
			  arglist.push_back( str );
			  break;
			}
		    }		  
		}
	      	      

	      // gfunc() spans positions p to q

	      // correct number of arguments? (typically 1 or 2)
	      // allow that first 'arglist' is function name
	      // if one short, assume we skipped the filter -- make to include all

	      std::string label = arglist[0] + "_func({" + arglist[1] + "})";
	      
	      // store internally too, if a genotype-function
	      if ( arglist[0] == "g" || arglist[0] == "gf" )
		gfunc[ label ] = arglist;
	
      
	      // swap out entire expression
	      // replace gfunc() with a variable

	      s->replace( p , (q-p+1) , label );
	      
	      // also replace any ';' with ':' so that we do not confuse the parser; these will be swapped back
	      // when this expression is parsed
	      
	      for (int i=0;i<s->size();i++)
		if ( (*s)[i] == ';' ) (*s)[i] = ':';

	    }
	  
	  ++i;  // next possible gfunc()
	}
      
      // are we sure there are no more gfunc()s in input?
      if ( ! found ) break;
      
    }
  
  return true;

}

void Eval::bind( SampleVariant & ref_svar , SampleVariant & ref_gvar , bool reset )   
{    
  

  // We have svar and gvar, as the gvar might be the consensus SampleVariant (e.g. 
  // from a flat alignment) whereas the variant meta-information will still be 
  // in the original SampleVariant.

  if ( reset ) reset_symbols();
  
  // Bind any gfunc()s  

  gvar = &ref_gvar;
  
  // Standard meta-information  
  bind( ref_svar.meta , false );
  
  // For assignment ( --> to meta )  
  func.attach( ref_gvar.meta );
}

void Eval::bind( SampleVariant & ref_svar , bool reset )   
{     

  if ( reset ) reset_symbols();

  gvar = &ref_svar;

  // For lookup ( --> from meta )
  bind( ref_svar.meta , false );

  // For assignment ( --> to meta )  
  func.attach( ref_svar.meta );

}


void Eval::bind( Variant & var , bool reset )   
{     

  if ( reset ) reset_symbols();
  
  // genotypes in consensus SampleVariant
  gvar = &var.consensus;
  
  // Standard meta-information  
  bind( var.meta , false );

  // For assignment ( --> to meta )  
  func.attach( var.meta );
  
}


template<class T> void Eval::assign_to( MetaInformation<T> & m )
{
  func.attach( m );
}

template<class T> void Eval::bind( MetaInformation<T> & m , bool reset )   
{
  
  //
  // Two scenarios: a variable might already exist in the MetaField (in which case
  // we use that type and value to specify the token;  alternatively, we may be 
  // expecting to populate the meta-field with a created token (i.e. assignment). At 
  // this stage we do not know the type of to-be-assiged meta-fields
  //
  
  // 
  // In GenMeta mode, we also have some special values: 
  //

  // GT    means genotype as VCF-style string  (0/1)
  // GT_A  means genotype as allele count (# non-ref alleles)
  
  std::map<std::string,std::set<Token*> >::iterator i = vartb.begin();
  while ( i != vartb.end() )
    { 
      
      std::set<Token*>::iterator tok = i->second.begin();
      
      while ( tok != i->second.end() )
	{
	  
	  // In GenMeta mode, we allow direct reading of the GT in a few ways
	  // GT GT_A GT_NR GT_NULL

	  bool set_genotype = false;
	  
	  if ( genmeta_mode ) 
	    {
	      
	      set_genotype = true;
	      
	      if ( i->first == PLINKSeq::VCF_GENOTYPE() )
		(*tok)->set( gvar->num_label( gvar->calls.genotype( indiv ) ) );
	      else if ( i->first == PLINKSeq::VCF_GENOTYPE_AC() )
		(*tok)->set( gvar->calls.genotype( indiv ).allele_count() );
	      else if ( i->first == PLINKSeq::VCF_GENOTYPE_NONREF() )
		(*tok)->set( gvar->calls.genotype( indiv ).allele_count() != 0 );
	      else if ( i->first == PLINKSeq::VCF_GENOTYPE_NULL() )
		(*tok)->set( gvar->calls.genotype( indiv ).null() );
	      else
		{
		  set_genotype = false;
		}	      
	    }

	  if ( ! set_genotype ) // variable will be represnted as a standard VarMeta
	    {

	      mType mt = MetaInformation<T>::type( i->first );

	      if ( mt != META_UNDEFINED ) 
		{
		  
		  meta_index_t midx = MetaInformation<T>::field( i->first );
		  
		  if ( midx.mt == META_FLAG ) 
		    {
		      (*tok)->set( m.has_field( i->first ) ) ; 
		    }
		  else if ( m.has_field( i->first ) )
		    {
		      
		      if ( midx.len == 0 || midx.len == 1 ) // scalars -- flags len == 0 (?check)
			{
			  if      ( midx.mt == META_INT   )  { (*tok)->set( m.get1_int( i->first ) );  }
			  else if ( midx.mt == META_FLOAT )  { (*tok)->set( m.get1_double( i->first ) ); }
			  else if ( midx.mt == META_TEXT  )  { (*tok)->set( m.get1_string( i->first ) ); }
			  else if ( midx.mt == META_BOOL  )  { (*tok)->set( m.get1_bool( i->first ) );   }	      
			}
		      else //vectors
			{		      
			  if      ( midx.mt == META_INT   )  { (*tok)->set( m.get_int( i->first ) );    }
			  else if ( midx.mt == META_FLOAT )  { (*tok)->set( m.get_double( i->first ) ); }
			  else if ( midx.mt == META_TEXT  )  { (*tok)->set( m.get_string( i->first ) ); }
			  else if ( midx.mt == META_BOOL  )  { (*tok)->set( m.get_bool( i->first ) );   }	      
			}
		    }
		  else
		    (*tok)->set(); // UNDEFINED
		}
	      else 
		{ 	      
		  (*tok)->set(); // UNDEFINED
		}
	    }	
	  ++tok;
	}
      ++i;
    }    
}



Token Eval::eval_gfunc( const std::string & expr , int gmode )
{
  
  //               gmode
  // gnull() --->  0
  // g()     --->  1 
  // gset()  --->  2


  // 'gnull' function modifies actual data
  
  bool gnull_mode = gmode == 0;
  bool gset_mode  = gmode == 2;

  //
  // Parse genotype-specific expression
  //

  if ( gvar == NULL ) return Token();

  Eval e;

  // replace any delimiter fields that were previously set to ':' instead of ';'

  std::string e2 = expr;

  for (int i=0;i<e2.size();i++)
    if ( e2[i] == ':' ) e2[i] = ';';

  if ( ! e.parse( e2 ) ) return Token();


  // set gen-meta mode

  e.genmeta( true );
  e.gvar = gvar;      
  
  //
  // For each individual, evaluate e and store resulting token
  //
  
  const int n = gvar->calls.size();
  
  std::vector<Token> rval(n);

  std::set<Token::tok_type> types;

  for (int j=0; j < n; j++)
    {

      // attach genotype meta-information
      
      e.indiv_pointer( j );
      e.bind( gvar->calls.genotype(j).meta );
      e.assign_to( gvar->calls.genotype(j).meta );

      // evaluate g(expr)

      e.evaluate();	  
      
      // gnull mode -- set genotype to null if does not pass
      
      if ( gnull_mode )
	{
	  bool passed = false;	  
	  bool valid = e.value( passed );	  
	  if ( (!valid) || (!passed) )
	    gvar->calls.genotype(j).null( true );
	}
      else if ( gset_mode ) 
	{	  
	  Helper::halt( "gset() not yet implemented" );
	}
      else // mode to return a vector of expressions
	{
	  // needs to eval to a scalar value for each individual
	  
	  Token tok = e.value();
	  Token::tok_type t = tok.type();
	  if ( ! tok.is_scalar() ) t = Token::UNDEF;	      

	  if ( t != Token::UNDEF ) 
	    {
	      types.insert( t );		  
	      rval[j] = tok;
	    }
	  else
	    rval[j] = Token();
	}
      
      // next individual
    }
  
  
  // 
  // Return something
  //
  
  if ( gnull_mode ) return Token( true );

  
  // Determine appropriate type of return vector
  // bool < int < float < string
  
  Token::tok_type t = Token::UNDEF;
  
  if      ( types.find( Token::STRING ) != types.end() ) t = Token::STRING_VECTOR;
  else if ( types.find( Token::FLOAT  ) != types.end() ) t = Token::FLOAT_VECTOR;
  else if ( types.find( Token::INT    ) != types.end() ) t = Token::INT_VECTOR;
  else if ( types.find( Token::BOOL   ) != types.end() ) t = Token::BOOL_VECTOR;
      
  Token tok;	  
  
  if ( t == Token::BOOL_VECTOR )	    
    {	      
      std::vector<bool> x( n );
      for (int i=0;i<n;i++) x[i] = rval[i].as_bool();
      tok.set( x );
    }
  else if ( t == Token::FLOAT_VECTOR ) 
    {
      std::vector<double> x( n );
      for (int i=0;i<n;i++) x[i] = rval[i].as_float();
      tok.set( x );
    }
  else if ( t == Token::INT_VECTOR ) 
    {	    
      std::vector<int> x( n );
      for (int i=0;i<n;i++) x[i] = rval[i].as_int();	      
      tok.set( x );
    }
  else if ( t == Token::STRING_VECTOR )
    {
      std::vector<std::string> x( n );
      for (int i=0;i<n;i++) x[i] = rval[i].as_string();	      
      tok.set( x );
    }

  return tok;

}
  


void Eval::bind( const Token * ntok )
{

  // we have have just created a meta-field from an assignment --
  // check whether this value now needs to be bound to any other
  // variable tokens (i.e. in next expression). ie. this function is
  // only called after an assignment is made and a meta-variable newly
  // created
  
  std::map<std::string,std::set<Token*> >::iterator i = vartb.find( ntok->name() );
  if ( i != vartb.end() )
    { 
      std::set<Token*>::iterator tok = i->second.begin();
      
      while ( tok != i->second.end() )
	{	  

	  if ( ntok == *tok ) { ++tok; continue; } 
	  
	  // we assume this new token has a well-defined type
	  // copy it over
	  
	  *(*tok) = *ntok;
	  
	  ++tok;	  
	}
    }
}

bool Eval::evaluate()
{
  for (int i=0; i<neval; i++)
    if ( is_valid ) is_valid = execute( output[i] );
  return is_valid;
}

bool Eval::valid() const 
{
  return is_valid;
}

std::string Eval::errmsg() const
{
  return errs;
}

void Eval::errmsg( const std::string & e )
{
  errs += e + "\n";
}


Token::tok_type Eval::rtype() const
{
  return e.type();
}

bool Eval::value(bool & b)
{

  if ( e.is_bool(&b) ) return true;

  int i;
  if ( e.is_int(&i) ) { b=i; return true; }

  // for vectors, evaluate to T is at least one T

  std::vector<bool> bv;
  if ( e.is_bool_vector(&bv) ) 
    {
      b = false;
      for (int i=0;i<bv.size();i++) 
	if ( bv[i] ) { b=true; break; }
      return true;
    }
  
  std::vector<int> iv;
  if ( e.is_int_vector(&iv) ) 
    {
      b = false;
      for (int i=0;i<iv.size();i++) 
	if ( iv[i] ) { b=true; break; }
      return true;
    }

  return false;
}

Token Eval::value() const
{
  return e;
}

bool Eval::value(int & i)
{
  if ( e.is_int(&i) ) return true;
  bool b;
  if ( e.is_bool(&b) ) { i=b; return true; }
  return false;
}

bool Eval::value(double & d)
{
  if ( e.is_float(&d) ) return true;
  int i;
  if ( e.is_int(&i) ) { d=i; return true; }
  bool b;
  if ( e.is_bool(&b) ) { d=b; return true; }
  return false;
}

bool Eval::value(std::string & s)
{
  if ( e.is_string(&s) ) return true;
  return false;
}


std::string Eval::result() const
{
  std::stringstream ss;
  ss << e;
  return ss.str();
}


void Eval::reset_symbols()
{
  // clear all variables -- put Tokens in the UNDEFINED state
  std::map<std::string,std::set<Token*> >::iterator i = vartb.begin();
  while ( i != vartb.end() )
    {
      std::set<Token*>::iterator k = i->second.begin();
      while ( k != i->second.end() )
	{
	  (*k)->set();
	  ++k;
	}
      ++i;
    }
  e.set(); // also clear end expression
}


void Eval::delete_symbols() 
{
  vartb.clear();
}


bool Eval::expand_vargs( std::string * s )
{
  // change vec(1,2,3) to vec(22,23,24,3)
  //        vec()      to vec(0)
  //        vec(22)    to vec(22,1)
  // etc i.e. put VA at end (as will come off stack first)


  std::vector<std::string> fname;
  fname.push_back("vec(");
  fname.push_back("int(");
  fname.push_back("str(");
  fname.push_back("bool(");
  
  for (int f = 0 ; f < fname.size() ; f++)
    while ( 1 ) 
      {
	
	size_t p = s->find( fname[f] );
	
	if ( p != std::string::npos ) 
	  {	      
	    
	    // this will also match 
	    //     blahint(X)
	    
	    if ( p > 1 )
	      {
		if ( ( (*s)[p-1] >= 'A' && (*s)[p-1] <= 'Z' ) || 
		     ( (*s)[p-1] >= 'a' && (*s)[p-1] <= 'z' ) || 
		     ( (*s)[p-1] >= '0' && (*s)[p-1] <= '9' ) || 
		     (*s)[p-1] >= '_' ) { p = std::string::npos; } 
	      }

	  }
	
	bool found = false;

	if ( p != std::string::npos )
	  {
	    
	    found = true;
	    
	    // count # of upper level ARGs, by counting commas until the next matching ")"
	    
	    //  vec( 1 , 2 , 3 )     -->  2 commas , 3 args
	    //  int( pow(2,3) , 7 )  ---> 1 comma  , 2 args
	    
	    
	    // look for closing brace
	    
	    int bc = 0;
	    int q = p;
	    int pp = p;
	    int commacnt = 0;

	    while ( ++q )
	      {
		
		// gone past end of string?
		if ( q == s->size() ) return false;
		
		char c = s->substr(q,1)[0];
		
		if ( c == '(' ) 
		  {
		    ++bc;  //includes first paran		    
		  }
		else if ( c == ')' )
		  {
		    --bc;
		    
		    // end of expression
		    if ( bc == 0 ) 
		      {
			break;
		      }
		  }
		else if ( bc == 1 && c == ',' ) ++commacnt;
	      }
	    
	    int len = q-p+1;	    
	    std::string orig = s->substr(p,q-p+1);    
	    orig = fname[f].substr( 0 , fname[f].size() - 1 ) + "_func("  + orig.substr( fname[f].size() ) ;
	    orig = orig.substr(0,orig.size()-1);
	    orig += "," + Helper::int2str( commacnt+1 ) + ")";
	    
	    s->replace( p , len , orig );
	    	    
	  }
	
	// are we sure there are no more gfunc()s in input?
	if ( ! found ) break;
	
      }	  

  return true;
}


bool Eval::expand_indices( std::string * s )
{

  // simply convert  X[i]  to  xfn(X,i)
  
  while ( 1 ) 
    {
      
      // search for opening index
      // assume it comes after a variable
      //  x[y]
      //  backtrack to get 'x'
      //  then forwards to get 'y'
      //  then make extract(x,y)
      // valid delimiters backwards are ' ' \n \t , % < > & | ! ~ = * + - / ; { } (

      // but if we hit ")" before anything else (except whitespace) then we need to jump to the next 
      // "(" and then look for an identifier 
      
      
      int p = s->find( "[" );
      if ( p == std::string::npos ) return true;
      int q = p;
      bool anything = false;

      while ( --q )
	{
	  if ( q < 0 ) return false;
	  if ( q == 0 ) { ++q; break; }
	  char c = s->substr(q,1)[0];
	  
	  if ( c == ')' )
	    {
	      // jump to next '(', allowing nesting
	      int nest = 1;
	      while (1) 
		{
		  --q;
		  if ( s->substr(q,1) == ")" ) ++nest;
		  else if ( s->substr(q,1) == "(" ) --nest;
		  if ( nest == 0 ) break; 
		}
	    }

	  if ( c == ',' || c == '&' || c == '%' || c == '>' || c == '<' || c == '|' || c == '(' ||
	       c == '!' || c == '~' || c == '^' || c == '=' || c == '*' || c == '+' || c == '-' || 
	       c == '/' || c == ';' || c == ':'  )
	    { ++q; break; }	  
	  
	  // whitespace is a delimiter (outside of parens) 
	  
	  if ( c == ' ' || c == '\n' || c == '\t' )
	    {
	      if ( anything ) { ++q; break; } 
	    }
	  else 
	    anything = true; // i.e. we must have encountered an identifier

	}
      
      std::string vec_idx = s->substr(q,p-q) ;
  
      std::string  arg_idx;
      int r = p;      
      while ( ++p )
	{  
	  if ( p == s->size() ) return false;
	  char c = s->substr(p,1)[0];
	  // do not allow nested indexing
	  if ( c == '[' ) return false;	  
	  if ( c == ']' )
	    {
	      arg_idx = s->substr(r+1,p-r-1);
	      break;
	    }  
	}
      
      std::string label = "element(" + vec_idx + "," + arg_idx + ")";
      s->replace( q , (p-q+1) , label );      

    } // search for next []  

  return true;    
}



void Eval::locate_symbols( std::vector<Token> & tok ) 
{
  
  // set pointers to output table for the updating, as this 
  // does not shift between runs
  
  for (int i=0; i<tok.size(); i++)
    if ( tok[i].is_variable() )
      {
	vartb[ tok[i].name() ].insert( &(tok[i]) );
      }
  reset_symbols(); // symbol table
}

