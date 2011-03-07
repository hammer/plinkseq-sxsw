#ifndef __PSEQ_TOKEN_H__
#define __PSEQ_TOKEN_H__

#include "variant.h"

#include <iostream>
#include <ostream>
#include <map>

class Token {
  
  friend std::ostream & operator<<( std::ostream & out , const Token & tok )
  {
    if ( tok.is_bool() ) out << ( tok.bval ? "T" : "F" );
    else if ( tok.is_int() ) out << tok.ival << "i";
    else if ( tok.is_float() ) out << tok.fval << "f";
    else if ( tok.is_string() ) out << tok.sval;
    else if ( tok.is_function() ) out << "fn(" << tok.name() << ")";
    else if ( tok.is_variable() ) out << "var(" << tok.name() << ")";
    else if ( tok.is_operator() ) out << Token::tok_unmap[ tok.type() ];
    else if ( tok.is_left_paren() ) out << "(";
    else if ( tok.is_right_paren() ) out << ")";
    else if ( tok.is_separator() ) out << ",";
    else out << ".";
    return out;
  }
  

 public:
  
  enum tok_type { UNDEF = 0 , 
		  INT , 
		  FLOAT , 
		  STRING , 
		  BOOL , 
		  ARG_SEPARATOR , 
		  FUNCTION , 
		  VARIABLE , 
		  MULTIPLY_OPERATOR , 
		  POWER_OPERATOR , 
		  DIVIDE_OPERATOR , 
		  MOD_OPERATOR, 
		  ADD_OPERATOR, 
		  SUBTRACT_OPERATOR ,
		  AND_OPERATOR , 
		  OR_OPERATOR , 
		  NOT_OPERATOR , 
		  EQUAL_OPERATOR , 
		  UNEQUAL_OPERATOR , 
		  GREATER_THAN_OPERATOR , 
		  LESS_THAN_OPERATOR ,
		  GREATER_THAN_OR_EQUAL_OPERATOR , 
		  LESS_THAN_OR_EQUAL_OPERATOR ,
		  ASSIGNMENT_OPERATOR , 
		  LEFT_PARENTHESIS , 
		  RIGHT_PARENTHESIS };

  
  
  // Constructors
  
  Token() { ttype = UNDEF; init(); }
  Token( const std::string & s );
  Token( const double d );
  Token( const int i );
  Token( const bool b );
  Token( const Token & );
  
  static void init();

  // Setters    
  void set( const std::string & s );
  void set( const double d );
  void set( const int i );
  void set( const bool b );
  void set();
  void function( const std::string & fn );
  void oper( Token::tok_type );
  void variable( const std::string & mf );


  // Evaluate: unary, binary and ternary operators
  Token operands( Token & t);
  Token operands( Token & t1, Token & t2);
  
  // Overload token operators; we'll overload '^' to mean power
  // function, even though it has wrong precedence (bitwise OR),
  // because in this context it is only used in simple paiwise ops

  Token operator!();
  Token & operator=(const Token & rhs);
  Token operator==(const Token & rhs);
  Token operator!=(const Token & rhs);
  Token operator+(const Token & rhs);
  Token operator^(const Token & rhs);
  Token operator-(const Token & rhs);
  Token operator*(const Token & rhs);
  Token operator/(const Token & rhs);
  Token operator%(const Token & rhs);
  Token operator<(const Token & rhs);
  Token operator>(const Token & rhs);
  Token operator>=(const Token & rhs);
  Token operator<=(const Token & rhs);
  Token operator&&(const Token & rhs);
  Token operator||(const Token & rhs);
    
    
  // Queries
  bool is_bool(bool * b = NULL) const ;
  bool is_string(std::string * s = NULL) const;
  bool is_float(double * f = NULL) const;
  bool is_int(int * i = NULL) const;
  bool is_operator() const;
  bool is_assignment() const { return ttype == ASSIGNMENT_OPERATOR; }
  bool is_function() const;
  bool is_ident() const;
  bool is_variable() const;    
  bool is_separator() const { return ttype == ARG_SEPARATOR; }
  bool is_left_paren() const { return ttype == LEFT_PARENTHESIS; }
  bool is_right_paren() const { return ttype == RIGHT_PARENTHESIS; }
  bool is_set() const { return ttype != UNDEF; } 

  
  // fetch actual data
  int as_int() const;
  double as_float() const;
  std::string as_string() const;
  bool as_bool() const;
  
  
  // fetch meta-data
  std::string name() const { return tname; }
  tok_type type() const { return ttype; }
  
  
  // name maps
  static std::map<std::string,tok_type> tok_map;   
  static std::map<tok_type,std::string> tok_unmap; 
  static std::map<std::string,int> fn_map;   
  
 private:
  
  // type and name (if function or variable)
  tok_type ttype;
  std::string tname;  
  
  // actual storage slots
  int ival;
  double fval;
  std::string sval;
  bool bval;

  
};


class TokenFunctions{ 
  
 public:
  
  void attach( MetaInformation<VarMeta> & m ) { meta = &m; }
  void attach( MetaInformation<GenMeta> & m ) { genmeta = &m; }
  
  Token fn_assign( Token & lhs , const Token & rhs );
  Token fn_set( const Token & tok );    
  Token fn_sqrt( const Token & tok );    
  Token fn_pow( const Token & tok , const Token & tok2 );
  Token fn_sqr( const Token & tok ) { return fn_pow( tok , Token(2) ); }
  Token fn_ifelse( const Token & cond , const Token & left , const Token & right );

  Token fn_n();
  Token fn_g( const Token & cond );
  Token fn_gmean( const Token & field , const Token & cond );

 private:  
  
  MetaInformation<VarMeta> * meta;
  MetaInformation<GenMeta> * genmeta;  // ignore this basically for now...
  
};

#endif
