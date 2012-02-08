#ifndef __HELPER_H__
#define __HELPER_H__

#include <string>
#include <vector>
#include <iostream>
#include <sstream>
#include <fstream>
#include <set>
#include <map>
#include <inttypes.h>


class Log;

extern Log plog;

class File;

class Variant;

// forward declaration, as used in Log
namespace Helper { std::string int2str(int); } 

class int2 {
 public:
    int p1;
    int p2;
    int2() { p1=p2=0; }
    int2(int p1, int p2) : p1(p1) , p2(p2) { }
    bool operator< (const int2 & b) const
        {
            return (p1 < b.p1 || (p1 == b.p1 && p2 < b.p2) );
        }

    bool operator== (const int2 & b) const
    {
      return p1 == b.p1 && p2 == b.p2;
    }

    bool operator!= (const int2 & b) const
    {
      return p1 != b.p1 || p2 != b.p2;
    }

};



////////////////////////////
// Logging and error handling

class Log {

  bool     silent_mode;   // write to STOUT?
  bool     output_file;   // write main output to file?
  bool     prolix_mode;   // write any prolix output to file?

  // keep track of warnings issued (and # of times)
  std::map<std::string,int> warnings;
  std::map<std::string,std::vector<std::string> > warnings_specific;

  // major output files
  std::ofstream file;
  std::ofstream prolix_file;

  // In R mode, use this as buffer 
  std::stringstream rstream;
  
  // data-dumper modes
  std::string  igrp;          // track i-group IDs
  std::string  igrp_header;   // used in tabular-mode
  int          coln;          // and col max.
  std::vector<std::string> cols;  // track header columns
  bool mode_tabular;      // mode (long or tabular)
  bool mode_header;       // header/numbers in long mode
  bool ignore_warnings; 
  bool early_warn;

 public:
  
  Log( bool s ,  // show? i.e. debug is silent by default
       std::string filename = "", 
       std::string prolix_filename = "" )
    { 	        
      
      silent_mode = !s;
      output_file = false;
      prolix_mode = false;      
      ignore_warnings = false;
      early_warn = false;

      mode_tabular = true;
      mode_header = false;
      igrp = igrp_header = "";
      
      if ( filename != "" ) logfile( filename );     
      if ( prolix_filename != "" ) prolix_logfile( prolix_filename );
   }
  
  ~Log()
    {
      if ( output_file ) file.close();
      if ( prolix_mode ) prolix_file.close();
    }


  //
  // Long-format/tabular helper functions for data rows
  //
  
  void tabular() { mode_tabular = true; }
  void longmode() { mode_tabular = false; }
  
  void header() { mode_header = true; } 
  void noheader() { mode_header = false; } 
  
  std::map<std::string,std::string> databuf0;
  std::map<std::string,std::map<std::string,std::string> > databufk;
  
  void data_reset() 
  { 
    databuf0.clear();
    databufk.clear();
    igrp = "";
    coln = 0; 
    cols.clear(); 
  }

  void data_group_header( const std::string & s )
  {
    igrp_header = s;    
    if ( mode_tabular ) 
      {
	// this should really be the first column, 
	// and always be present
	(*this) << s;
      }   
  }
  
  void data_header( const std::string & s ) 
  { 
    ++coln;    
    if ( mode_tabular ) 
      {
	//we assume an igrp_header has already been written
	//so always put in a tab here
	(*this) << "\t" << s;
      }    
    cols.push_back(s);    
  } 

  
  void data_header_done()
  {
    if ( mode_tabular ) (*this) << "\n"; 
    coln = 0;
  }

  template<class T> bool data( const T & v )
    {
      return data( v , "0" );
    }

  template<class T> bool data( const T & v , int k )
    {
      return data( v , Helper::int2str(k) );
    }

  template<class T> bool data( const T & v , const std::string & k )
    {
      // start, coln will be 0; use 1-based codes for headers
      if ( coln == cols.size() ) return false;
      return data( v , cols[ coln++ ] , k );
    }
  
  template<class T> bool data( const T & v , const std::string & j, int k )
    {
      return data( v , j , Helper::int2str(k) );
    }
  
  template<class T> bool data( const T & v , const std::string & j , const std::string & k )
    {
           
      if ( mode_tabular ) 
	{
	  std::stringstream ss;
	  ss << v;
	  if ( k == "0" ) databuf0[ j ] = ss.str();
	  else databufk[ k ][ j ] = ss.str();
	}
      else
	{
	  (*this) << igrp << "\t"
		  << j 
		  << ( k == "0" ? "" : "{" + k + "}" ) 
		  << "\t" 
		  << v << "\n";
	}
      
      return true;
    }
  
  
  template<class T> void data_group( const T & i )
  {
    std::stringstream ss;
    ss << i;
    igrp = ss.str();
  }

  void print_data_group()
  { 

    if ( mode_tabular )
      {
	
	if ( databufk.size() == 0 )
	  {
	    (*this) << igrp;
	    for (int i=0; i<cols.size(); i++)
	      (*this) << "\t" << databuf0[ cols[i] ];
	    (*this) << "\n";
	  }
	else
	  {
	    
	    std::map<std::string,std::map<std::string,std::string> >::iterator ik = databufk.begin();
	    while ( ik != databufk.end() )
	      {

		(*this) << igrp;
				
		for (int i=0; i<cols.size(); i++)
		  {
		    std::map<std::string,std::string>::iterator f0 = databuf0.find( cols[i] );
		    if ( f0 != databuf0.end() )
		      {
			(*this) << "\t" << f0->second;
		      }
		    else
		      {
			std::map<std::string,std::string>::iterator f1 = ik->second.find( cols[i] );
			if ( f1 != ik->second.end() )
			  {
			    (*this) << "\t" << f1->second;
			  }
			else
			  (*this) << "\t.";
		      }
		  }
		
		(*this) << "\n";
		
		++ik;
	      }
	  }
      }

    // clear stores
    databuf0.clear();
    databufk.clear();
    
    // reset row counter
    coln = 0;
  }
  
  // General mode switches
  void silent(const bool b) { silent_mode = b; }
  bool silent() const { return silent_mode; }

  void precision( const int p )
  {
    std::cout.precision(p);
  }
  
  bool logfile() const { return output_file; } 
  
  void logfile( const std::string & f )     
  {
    file.open( f.c_str() , std::ios::out );
    output_file = true;
  }
  
  void close_logfile()     
  {
    if ( output_file ) file.close();
    output_file = false;
  }
  
  bool prolix_logfile() const { return prolix_mode; }

  void prolix_logfile(const std::string f)
    {
      prolix_file.open( f.c_str() , std::ios::out );
      prolix_mode = true;
    }
  
  void close_prolix_logfile()
  {
    if ( prolix_mode ) prolix_file.close();
    prolix_mode = false;
  }

  void counter( const std::string & msg )
  {
#ifndef R_SHLIB
    if ( silent_mode ) return;
    std::cerr << msg << "        \r";
    std::cerr.flush();
#endif
  }  

  
  template<class T>
  Log & operator<<(const T & msg) 
  {
    // Output to file
    if ( output_file )
      {
	file << msg;
	file.flush();
      }

    if ( silent_mode ) return *this;
    
#ifdef R_SHLIB 
    rstream << msg;
#else
    std::cout << msg;
#endif

    return *this;
  }
  

  void flush()
  {    
#ifdef R_SHLIB 
    // nothing
#else
    if ( ! silent_mode ) std::cout.flush();
#endif    
  }

  std::string R_flush()
    {
      std::string r = rstream.str();
      rstream.str( std::string() );
      return r;
    }

  //
  // Warnings and errors, to be handled separately
  //

  void warn(const std::string & msg, const std::string & spec = "" );

  void warn(const std::string & msg, const std::vector<std::string> & spec );
  
  void show_warnings(const bool b ) { ignore_warnings = !b; } 

  void early_warnings(const bool b ) { early_warn = b; }
  
  void print_warnings();



  //
  // Verbose output information
  //

  template<class T>
  void prolix(const T & msg)
  {
    if ( silent_mode ) return;
    if ( prolix_mode ) 
      {
	prolix_file << msg;
	prolix_file.flush();
      }
  }
  
 
};


namespace Helper 
{
    
  void halt( const std::string & msg);
  void NoMem();

  ////////////////////////////
  // File handling
  
  bool remove_file( const std::string & );

  ////////////////////////////
  // Conversion functions
  

  void str2upper( std::string & );
  std::string search_replace( std::string & , const std::string & , const std::string & );
  std::string header( const std::string & s , const int len = 60 , const std::string & rep = "." ); 

  std::string int2str(int);
  std::string longint2str(long int);
  
  std::string uint64_t2str(uint64_t);
  std::string lng2str(long int);
  std::string dbl2str(double,int prc = -1);
  std::string dbl2str_fixed(double, int prc = -1);
  
  std::string flt2str(const float, const int prc = -1);
  
  bool valid_name( const std::string & s );

  template <class T>
    bool from_string(T& t,
		     const std::string& s,
		     std::ios_base& (*f)(std::ios_base&));

  
  template <class T> 
    T lexical_cast( const std::string & s )
    {
      T t;
      bool okay = from_string<T>( t , s , std::dec ); 
      if ( ! okay ) plog.warn("problem converting " + s );
      return t;
    }
  
  bool str2uint64_t(const std::string s , uint64_t & i);
  bool str2int(const std::string & s , int & i);
  bool str2dbl(const std::string & s , double & i);  
  int str2int(const std::string & s);
  double str2dbl(const std::string & s);  
  bool str2bool( const std::string & s , const std::string & miss = "." );

  // C-style versions

  bool str2uint64_t(const char * c , uint64_t & i);
  bool str2int(const char * c , int & i);
  bool str2dbl(const char * c , double & i);  
  int str2int(const char * c );
  double str2dbl(const char * c);
  bool str2bool( const char * c );


  template <class T>
    bool from_string(T& t,
		     const std::string& s,
		     std::ios_base& (*f)(std::ios_base&))
    {
      std::istringstream iss(s);
      return !(iss >> f >> t).fail();
    }

  template<class T>
  void append( std::vector<T> & x , const std::vector<T> & y )
    {
      for (int i=0;i<y.size();i++)
	x.push_back( y[i] );
    }

  template <class T>
    std::string print( const std::vector<T> & x , const bool columns = false , const bool nums = false , const std::string & sep = " " )
    {
      std::stringstream ss;
      for (int i=0; i<x.size(); i++) 
	{
	  if ( columns ) 
	    {
	      if ( nums ) ss << i << "\t";
	      ss << x[i] << "\n";
	    }
	  else
	    {
	      if ( i != 0 ) ss << sep;
	      ss << x[i];
	    }	  
	}
	return ss.str();
    }

  
  std::string quote_value( const std::string & );
  std::string unquote( const std::string & );

  bool inCommaList( const std::string & , const std::string & );
  std::set<std::string> parseCommaList(const std::string & lst);
  std::vector<std::string> parse(const std::string &, const std::string &, bool empty = false );
  std::vector<std::string> quoted_parse(const std::string &, const std::string & s = ",", bool empty = false );
  std::vector<std::string> whitespace(const std::string & sline);
  std::map<std::string,std::string> quoted_comma_keypair_split( const std::string & item );
  std::string stringizeKeyPairList( const std::map<std::string,std::string> & , bool show_keys = true );
  std::string stringize( const std::set<std::string> & , const std::string & delim = "," );
  std::string stringize( const std::vector<std::string> & , const std::string & delim = "," );

  std::vector<std::string> char_split( const std::string & , const char , bool empty = true );
  std::vector<std::string> char_split( const std::string & , const char , const char , bool empty = true );
  std::vector<std::string> char_split( const std::string & , const char , const char , const char , bool empty = true );

  std::vector<std::string> quoted_char_split( const std::string & s , const char c , bool empty = true );
  std::vector<std::string> quoted_char_split( const std::string & s , const char c , const char , bool empty = true );
  std::vector<std::string> quoted_char_split( const std::string & s , const char c , const char , const char , bool empty = true );

  std::vector<std::string> tokenizeLine(std::ifstream &);

  //
  // faster C-string, in place tokenizer for VCF genotype entries
  //
  
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
  

  
  //
  // String list helper functions
  //

  
  std::string filelist2commalist( const std::string & f );
  void inserter( std::set< std::string > & strset , const std::string & filespec );

  bool ends_with( const std::string & , const std::string & );
  
  std::string remove_tags( const std::string & );
  std::string add_tags( const std::string & );

  std::string chrCode(int c, bool prefix=true);

  int chrCode(const std::string & c);

  std::string coordinate( const int chr, 
			  const int bp1 = 0 , 
			  const int bp2 = 0 );

  std::string metatype_summary( const bool pretty = true );

  double hwe( const Variant & , int * phom1 = NULL , int * phet = NULL , int * phom2 = NULL );

  double SNPHWE(int obs_hets, int obs_hom1, int obs_hom2);  


  ////////////////////////////
  // File handling

  void ensure_folder( std::string & );
  bool fileExists( const std::string & );
  bool fileExists( File* );
  bool checkFileExists( const std::string &);
  bool checkFileExists( const std::vector<std::string> &);
  bool checkFileExists(File *);
  bool checkFileRectangular(File *);
  std::string fullpath( const std::string & );

  ////////////////////////////
  // Pretty-printing
  
  std::string sw(const std::string & s , int n);
  std::string sw(double d , int n);
  std::string sw(double d , int f, int n);
  std::string sw(int i , int n);
  std::string sw(uint64_t , int n);
  bool realnum(double d);
    
  /////////////////////////////  
  // Type specificiation

  bool is_int(const std::string & s );
  bool is_long(const std::string & s );
  bool is_float(const std::string & s );
  bool is_string(const std::string & s );
  bool is_char(const std::string & s );
  bool is_text(const std::string & s );
  bool is_bool(const std::string & s );
  bool is_flag(const std::string & s );
  


  ///////////////////////////
  // Statistics
  
  double chi2x2(int,int,int,int);
  bool realnum(double);
  
    ///////////////////////////
    // Loaders

    std::vector<std::string> load_string_list(const std::string &); 

}



////////////////////////////
// SQL BLOB object

class blob {
  
 public:
  
  std::string s;
  
  const char * p;
  
  int l;
  
  blob() { }
  
  blob(const std::string & t) 
    {
      s = t;
      p = s.data();
      l = s.size();
    }
  
  void set_string(const std::string & tmp)
  {
    s = tmp;
    p = s.data();
    l = s.size();
  }
  
  std::string get_string()
    {
      return std::string((const char *)p,l);
    }
  
};

class GroupInfo {

 public:

    GroupInfo() 
	{
	    idx = 0;
	    name = "";
	    temp = false;
	    description = "";
	}

    uint64_t idx;
    std::string name;
    bool temp;
    std::string description;

    bool operator<( const GroupInfo & b ) const 
	{
	  if ( idx == b.idx ) return name < b.name;
	  return idx < b.idx;
	}

    friend std::ostream & operator<<( std::ostream & out , const GroupInfo & g)
	{
	    out << g.idx << " "
		<< g.name << " "
		<< g.temp << " "
		<< g.description;
	    return out;
	}
};


struct int_string_pair {
    int_string_pair(int i, std::string s) : i(i) , s(s) { } 
    int i;
    std::string s;  
};


class int_range {

  int lwr, upr;
  bool has_lwr, has_upr;  

 public:

  friend std::ostream & operator<<( std::ostream & out , const int_range & r )
  {
    if ( r.has_lwr ) out << r.lwr ; else out << "*";
    if ( r.has_upr ) 
      {
	if ( r.upr != r.lwr || ! r.has_lwr ) out << ":" << r.upr ; 
      }
    else 
      if ( r.has_lwr ) out << ":*";
  }

  int_range() { reset(); }

  // smode =  0  2 means  2:2
  //         -1  2 means  *:2 
  //         +1  2 means  2:*

  int_range( const std::string & , const int smode = 0 ); 
  void set( const std::string & , const int smode = 0 );
  void reset();
  void set( const int a, const int b ) { lwr=a; upr=b; has_lwr=has_upr=true; }
  bool in( const int i ) const; 
  int lower() const { return has_lwr ? lwr : -1; }
  int upper() const { return has_upr ? upr : -1; }  
  bool has_lower() const { return has_lwr; } 
  bool has_upper() const { return has_upr; } 

  bool operator<( const int_range & rhs ) const
  {
    if ( rhs.has_lwr && ! has_lwr ) return true;
    if ( has_lwr && ! rhs.has_lwr ) return false;
    if ( has_lwr ) 
      {
	if ( lwr < rhs.lwr ) return true;
	if ( lwr > rhs.lwr ) return false;
      }    
    if ( rhs.has_upr & ! has_upr ) return false;
    if (     has_upr & ! rhs.has_upr ) return true;    
    return upr < rhs.upr;   
  }
};

class dbl_range {
  double lwr, upr;
  bool has_lwr, has_upr;  
 public:

  friend std::ostream & operator<<( std::ostream & out , const dbl_range & r )
  {
    if ( r.has_lwr ) out << r.lwr ; else out << "*";
    if ( r.has_upr ) 
      {
	if ( r.upr != r.lwr || ! r.has_lwr ) out << ":" << r.upr ; 
      }
    else 
      if ( r.has_lwr ) out << ":*";
    return out;
  }
  dbl_range() { reset(); }
  dbl_range( const std::string & , const int smode = 0 );
  void set( const std::string & , const int smode = 0 );
  void reset();
  void set( const double a, const double b ) { lwr=a; upr=b; has_lwr=has_upr=true; }
  bool in( const double d ) const; 
  double lower() const { return has_lwr ? lwr : 0; }
  double upper() const { return has_upr ? upr : 0; }  
  bool has_lower() const { return has_lwr; } 
  bool has_upper() const { return has_upr; } 
  bool operator<( const dbl_range & rhs ) const
  {
    if ( rhs.has_lwr && ! has_lwr ) return true;
    if ( has_lwr && ! rhs.has_lwr ) return false;
    if ( has_lwr ) 
      {
	if ( lwr < rhs.lwr ) return true;
	if ( lwr > rhs.lwr ) return false;
      }    
    if ( rhs.has_upr & ! has_upr ) return false;
    if (     has_upr & ! rhs.has_upr ) return true;    
    return upr < rhs.upr;   
  }
};


#endif
