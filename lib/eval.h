#ifndef __PSEQ_EVAL_H__
#define __PSEQ_EVAL_H__

#include "token.h"

#include <cstring>
#include <ostream>
#include <map>
#include <set>
#include <vector>


class Eval {
    
  friend std::ostream & operator<<( std::ostream & out , Eval & rhs )
  {
    out << rhs.e;
    return out;
  }

  
 public:    
  
  // constructors

  Eval() 
    { init(); }

  Eval( const std::string & input ) 
    { init(); parse(input); }
  
  void init(); 

  // primary functions

  bool parse( const std::string & input );

  void bind( SampleVariant & svar , bool reset = true );
  void bind( SampleVariant & svar , SampleVariant & gvar , bool reset = true );
  void bind( Variant & var , bool reset = true );

  template<class T> void assign_to( MetaInformation<T> & m );

  template<class T> void bind( MetaInformation<T> & m , bool reset = true );

  void bind( const Token * );

  void eval_gfunc( SampleVariant & svar );

  bool evaluate();
    
  bool valid() const;
  
  std::string errmsg() const;

  std::string result() const;

  bool requires_genotypes() const
  { return gfunc.size() > 0; } 

  // queries into value of expression
  bool value(bool & b);
  bool value(int & );
  bool value(double &);
  bool value(std::string &);
  
 private:

  //work horses
  bool get_token( std::string & input ,  Token & );
  bool previous_value;
  bool execute( const std::vector<Token> & );
  bool shunting_yard( const std::string & input, std::vector<Token> & );

  // helpers
  int op_preced(const Token & tok );
  bool op_left_assoc(const Token & tok );
  unsigned int op_arg_count(const Token & tok );

  bool variant_mode; // as opposed to SampleVariant

  // expression in RPN notation (post parse) (for each eval)
  std::vector< std::vector<Token> >output;   
  
  // keep track of state (errors?)
  bool is_valid;
  std::string errs;
  
  // slot for final expression value
  Token e;
  
  // Symbol table for variables, and related functions
  bool extract_gfunc( std::string * );
  std::map<std::string,std::set<Token*> > vartb;
  void delete_symbols();
  void reset_symbols();
  void locate_symbols( std::vector<Token> & );

  // G-functions
  std::map<std::string,int> gdef; 
  std::map<std::string,std::vector<std::string> > gfunc;

  // Keep token-functions in one place
  TokenFunctions func;

  // the number of evals
  int neval;

};




#endif

