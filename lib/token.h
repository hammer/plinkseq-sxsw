#ifndef __PSEQ_TOKEN_H__
#define __PSEQ_TOKEN_H__

#include "variant.h"

#include <iostream>
#include <ostream>
#include <map>

class Eval;

class Token {
  
  friend std::ostream & operator<<( std::ostream & out , const Token & tok )
  {
    if ( tok.is_vector() )
      {
	int l = tok.size() > 5 ? 5 : tok.size() ;

	out << "[";
	
	for ( int i=0; i<l; i++ ) 
	  {
	    if ( i ) out << "," ; 
	    if      ( tok.is_bool_vector() ) out << ( tok.bvec[i] ? "T" : "F" );
	    else if ( tok.is_int_vector() ) out << tok.ivec[i];
	    else if ( tok.is_float_vector() ) out << tok.fvec[i];
	    else if ( tok.is_string_vector() ) out << tok.svec[i];
	  }

	if ( tok.size() > l ) out << "... ("<< tok.size()<< " elements) ";
	
	if      ( tok.is_bool_vector() ) out << "]b";
	else if ( tok.is_int_vector() ) out << "]i";
	else if ( tok.is_float_vector() ) out << "]f";
	else if ( tok.is_string_vector() ) out << "]s";
      }    
    else if ( tok.is_bool() ) out << ( tok.bval ? "T" : "F" );     
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
		  INT_VECTOR , 
		  FLOAT_VECTOR , 
		  STRING_VECTOR , 
		  BOOL_VECTOR ,		  
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

  Token( const std::vector<std::string> & s );
  Token( const std::vector<double> & d );
  Token( const std::vector<int> & i );
  Token( const std::vector<bool> & b );

  Token( const Token & );
  
  static void init();

  // Setters    

  void set( const std::string & s );
  void set( const double d );
  void set( const int i );
  void set( const bool b );

  void set( const std::vector<std::string> & s );
  void set( const std::vector<double> & d );
  void set( const std::vector<int> & i );
  void set( const std::vector<bool> & b );

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

  Token operator!() const;
  Token & operator=(const Token & rhs);
  Token operator==(const Token & rhs) const;
  Token operator!=(const Token & rhs) const;
  Token operator+(const Token & rhs) const;
  Token operator^(const Token & rhs) const;
  Token operator-(const Token & rhs) const;
  Token operator*(const Token & rhs) const;
  Token operator/(const Token & rhs) const;
  Token operator%(const Token & rhs) const;
  Token operator<(const Token & rhs) const;
  Token operator>(const Token & rhs) const;
  Token operator>=(const Token & rhs) const;
  Token operator<=(const Token & rhs) const;
  Token operator&&(const Token & rhs) const;
  Token operator||(const Token & rhs) const;
    
    
  // Queries

  int size() const;

  bool is_bool(bool * b = NULL) const ;
  bool is_string(std::string * s = NULL) const;
  bool is_float(double * f = NULL) const;
  bool is_int(int * i = NULL) const;

  bool is_scalar() const;
  bool is_vector() const;

  bool is_bool_vector( std::vector<bool> * b = NULL ) const;
  bool is_string_vector( std::vector<std::string> * s = NULL ) const;
  bool is_float_vector( std::vector<double> * f = NULL ) const;
  bool is_int_vector( std::vector<int> * f = NULL ) const;

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
  int         as_int()     const;
  double      as_float()   const;
  std::string as_string()  const;
  bool        as_bool()    const;

  // a specific element from a vector, no type-conversion
  int         int_element(const int ) const;
  double      float_element(const int) const;
  std::string string_element(const int) const;
  bool        bool_element(const int) const;

  std::vector<int>         as_int_vector()    const;
  std::vector<double>      as_float_vector()  const;
  std::vector<std::string> as_string_vector() const;
  std::vector<bool>        as_bool_vector()   const;

  
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
  int           ival;
  double        fval;
  std::string   sval;
  bool          bval;

  std::vector<int>          ivec;
  std::vector<double>       fvec;
  std::vector<std::string>  svec;
  std::vector<bool>         bvec;
  
};


class TokenFunctions{ 
  
 public:
  

  Token fn_set( const Token & tok ) const;    
  Token fn_sqrt( const Token & tok ) const;    
  Token fn_pow( const Token & tok , const Token & tok2 ) const;
  Token fn_sqr( const Token & tok ) const { return fn_pow( tok , Token(2) ); }
  Token fn_ifelse( const Token & cond , const Token & left , const Token & right ) const;

  // vector functions
  Token fn_vec_length( const Token & tok) const;
  Token fn_vec_extract( const Token & tok , const Token & idx ) const;
  Token fn_vec_min( const Token & tok ) const;
  Token fn_vec_maj( const Token & tok ) const;
  Token fn_vec_sum( const Token & tok ) const;
  Token fn_vec_mean( const Token & tok ) const;
  Token fn_vec_sort( const Token & tok ) const;

  // genotype extraction
  Token fn_vec_g( const Token & tok , Eval * e ) const;
  Token fn_vec_gnull( const Token & tok , Eval * e ) const;
  Token fn_vec_gset( const Token & tok , Eval * e ) const;
  
  // phenotype extraction/assignment
  Token fn_vec_pheno( const Token & tok ) const;
  Token fn_vec_1pheno( const Token & rhs , int ) const;
  Token fn_assign_pheno( Token & lhs , const Token & rhs );


  Token fn_vec_new_float( const Token & tok ) const;
  Token fn_vec_new_int( const Token & tok ) const;
  Token fn_vec_new_str( const Token & tok ) const;
  Token fn_vec_new_bool( const Token & tok ) const;
  
  Token fn_vec_any( const Token & tok1 , const Token & tok2 ) const;
  Token fn_vec_count( const Token & tok1 , const Token & tok2 ) const;

  void attach( MetaInformation<VarMeta> & m );
  void attach( MetaInformation<GenMeta> & m );
  
  Token fn_assign_var( Token & lhs , const Token & rhs );
  Token fn_assign_gen( Token & lhs , const Token & rhs );
  
 private:  
  
  MetaInformation<VarMeta> * meta;
  MetaInformation<GenMeta> * genmeta;
  
};

#endif
