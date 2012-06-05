
#ifndef __META_H__
#define __META_H__

#include <map>
#include <set>
#include <string>
#include <vector>
#include <iostream>
#include <sstream>

#include "defs.h"
#include "helper.h"

extern Log plog;

typedef std::string                        meta_name_t;
typedef int                                meta_key_t;
typedef std::pair<mType,meta_key_t>        meta_typed_key_t;

typedef std::vector<std::string> meta_string_set_t;
typedef std::vector<int>         meta_int_set_t;
typedef std::vector<double>      meta_double_set_t;
typedef std::vector<bool>        meta_bool_set_t;

typedef std::map<meta_key_t, meta_string_set_t>  meta_string_t;
typedef std::map<meta_key_t, meta_int_set_t>     meta_int_t;
typedef std::map<meta_key_t, meta_double_set_t>  meta_double_t;
typedef std::map<meta_key_t, meta_bool_set_t>    meta_bool_t;
typedef std::set<meta_key_t>                     meta_flag_t;

struct meta_index_t { 
  meta_key_t   key;  
  meta_name_t  name;
  mType        mt;
  int          len;
  std::string  description;
  
  bool operator<(const meta_index_t & b) const
  {
    if ( mt < b.mt ) return true;
    if ( mt > b.mt ) return false;
    return key < b.key;
  }
  
};

// also see enum definition in defs.h

struct MiscMeta      { static const int n = 0; };
struct VarMeta       { static const int n = 1; };
struct GenMeta       { static const int n = 2; };
struct LocMeta       { static const int n = 3; };
struct RefMeta       { static const int n = 4; };
struct FileMeta      { static const int n = 5; };
struct IndivMeta     { static const int n = 6; };
struct AlleleMeta    { static const int n = 7; };
struct VarFilterMeta { static const int n = 8; };

class MetaMeta {

 public:
  
  static void load( const std::string & );
  static void clear();
  
  /// indicator to show that meta-field show be consider not sample-specific
  static void set_static( const std::string & s) { pop_static.insert(s); } 
  static bool static_variant( const std::string & );
  
  // indicator to show whether all meta-fields should be forced into the consensus

  static void set_force_consensus( const bool b ) { force_consensus_mode = b; }
  static bool force_consensus() { return force_consensus_mode; }

  /// Hide/show specific signals?
  static void hide( const std::string & );
  static void show( const std::string & );
  
  // Internal variant (that will never be shown, unless subsequently explicitly declared)
  static void is_internal( const std::string & );
  static void is_external( const std::string & attrib )
    {
      if ( internal_mask.find( attrib ) != internal_mask.end() )
	{
	  internal_mask.erase( internal_mask.find( attrib ) );
	}  
      if ( internal_mask.size() == 0 ) masked_internal = false;
    }
  

  // Should I display this variant?
  static bool display( const std::string & );
  

 private:

  static std::set<std::string> pop_static;
  static bool force_consensus_mode;

  static bool masked_hide;
  static bool masked_show;
  static bool masked_internal;

  static std::set<std::string> show_mask;
  static std::set<std::string> hide_mask;
  static std::set<std::string> internal_mask;

};

void registerMetatype( const std::string & name, 
		       mType mt, 
		       int num , 
		       int grp , 
		       const std::string & desc );


template<class T> 
class MetaInformation {
  
 public:
  
  void parse( const std::string & str, char delim = ';' , bool autoadd = false , const std::string * prefix = NULL )
  {
    
    // Take delimited string and parse
    // Format is  " ;" ->  X=Y  ->  X single      
    
    // X=1;Y=2,3;Z;A=House;
    // X="1 ;23";Y=22;
    // but allow quotes
    
    const bool escape_quotes = true;
	    
    int ntok;
    Helper::char_tok tokens( str , &ntok , delim , escape_quotes ); 
    
    //std::vector<std::string> tokens = Helper::quoted_parse( str , delim , empty );
    
    const int sz = tokens.size();
    
    for ( int i = 0 ; i < sz ; i++ )
      {
	
	// Expecting a key=value pairing; allow quotes
	
	int ntok2;
	Helper::char_tok tokens_b( tokens(i) , &ntok2 , '=' , escape_quotes );
	
	// skip handle weird things like trying to parse '==' 
	if ( tokens_b.size() == 0 ) continue;
	
	// get (with optional prefix -- for RefVariant meta)
	
	const std::string key = prefix ? *prefix + "_" + tokens_b[0] : tokens_b[0];	      	      
	
	// a key/value pair
	
	if ( tokens_b.size() == 2 ) // a key = value pair
	  {	    
	    if ( autoadd ) 
	      {
		if ( ! MetaInformation::exists( key ) ) 
			MetaInformation::field( key , META_TEXT );
	      }
	    
	    // will find the appropriate type and add value to key
	    parse_set( key, Helper::unquote( tokens_b[1] ) );
	  } 	
	else 
	  { 
	    if ( autoadd ) MetaInformation::field( key , META_FLAG );	    
	    set( key ); // set as META_FLAG 
	  } 
	
      }
    
  }
  
    
    
    //
    // Set value functions
    //    
    

    // 1) FLAGs, i.e no value, expect string, but needs to 
    //    have been pre-registered
    
  void set( const meta_name_t & name )
  {	  	    
    if ( ! MetaInformation::exists( name ) ) return;
    meta_index_t midx = MetaInformation::field(name,META_FLAG);      
    if ( midx.mt != META_FLAG ) return; 
    m_flag.insert( midx.key );
  }

  // A key is specified directly, so assume it exists
  void set( const meta_key_t & key )
  {	  	    
    m_flag.insert( key );
  }
    
  // 2) Key/value pairs
    
  void parse_set(const meta_name_t & name , const std::string & str )
  {
    
    std::vector<std::string> parsed = Helper::quoted_parse( str );
    
    mType mt = MetaInformation::type( name ) ;
    
	if ( mt == META_INT )
	  {
	    std::vector<int> value;
	    for (unsigned int i=0; i<parsed.size(); i++)
	      {
		bool okay = true;
		int v;
		try { v = Helper::lexical_cast<int>(parsed[i]); }
		catch ( std::exception& e) { okay = false; }
		if ( okay ) value.push_back( v );
	      }
	    set( name , value );
	  }
	else if ( mt == META_FLOAT )
	  {
	    std::vector<double> value;
	    for (unsigned int i=0; i<parsed.size(); i++)
	      {
		bool okay = true;
		double v;
		try { v =  Helper::lexical_cast<double>(parsed[i]); }
		catch ( std::exception& e) { okay = false; }
		if ( okay ) value.push_back( v );
	      }
	    set( name , value );
	  }
	else if ( mt == META_BOOL )
	  {
	    std::vector<bool> value;
	    for (unsigned int i=0; i<parsed.size(); i++)
	      {
		bool okay = true;
		bool v;
		try { v =  Helper::lexical_cast<bool>(parsed[i]); }
		catch ( std::exception & e) { okay = false; }
		if ( okay ) value.push_back( v );
	      }
	    set( name , value );
	  }
	else if ( mt == META_TEXT )
	  {
	    set( name , parsed );
	  }

      }
    
    
    //
    //  Given NAME, add a std::vector of typed items
    //
    
    void set(const meta_name_t & name , const std::vector<int> & value )
      {	  
	set( MetaInformation::field(name).key , value );
      }
    
    void set(const meta_name_t & name , const std::vector<double> & value )
      {	  
	set( MetaInformation::field(name).key , value );
      }
    
    void set(const meta_name_t & name , const std::vector<std::string> & value )
      {	  
	set( MetaInformation::field(name).key , value );
      }
    
    void set(const meta_name_t & name , const Helper::char_tok & value )
	{
	    set( MetaInformation::field(name).key , value );
	}
    
    void set(const meta_name_t & name , const std::vector<bool> & value )
      {	  
	set( MetaInformation::field(name).key , value );
      }
   
    
    //
    // Actual functions that do the insertion
    //

    void set(const meta_key_t key , const std::vector<int> & value)
    {
      m_int[key] = value;
    }
    
    void set(const meta_key_t key , const std::vector<double> & value)
    {
      m_double[key] = value;
    }
    
    void set(const meta_key_t key , const std::vector<std::string> & value)
    {
      m_string[key] = value;
    }
    
    void set(const meta_key_t key , const Helper::char_tok & value)
    {
	std::vector<std::string> & t = m_string[key];
	t.resize( value.size() , "." );
	for (int i=0; i<value.size(); i++) t[i] = value[i];
    }

    void set(const meta_key_t key , const std::vector<bool> & value)
    {
      m_bool[key] = value;
    }
    
    
    //
    // Convenience functions, to set specific values
    //      
    

    void set(const meta_name_t & name , const int value ) 
    {
      set( MetaInformation::field(name).key , value );
    }
    
    void set(const meta_name_t & name , const double value ) 
    {
      set( MetaInformation::field(name).key , value );
    }
    
    void set(const meta_name_t & name , const bool value ) 
    {
      set( MetaInformation::field(name).key , value );
    }
    
    void set(const meta_name_t & name , const std::string & value ) 
    {
      set( MetaInformation::field(name).key , value );
    }
    
    void set(const meta_name_t & name , const char * value )
    {
      set( MetaInformation::field(name).key , std::string(value) );
    }
    
    
    //
    // Add a std::vector, given a single value
    //

    void set(const meta_key_t key , const int value )
	{
	    std::vector<int> t;
	    t.push_back(value);
	    m_int[key] = t;
	}

    void set(const meta_key_t key , const double value )
	{
	    std::vector<double> t;
	    t.push_back(value);
	    m_double[key] = t;
	}

    void set(const meta_key_t key , const std::string & value )
	{
	    std::vector<std::string> t;
	    t.push_back(value);
	    m_string[key] = t;
	}

    void set(const meta_key_t key , const bool value )
	{
	    std::vector<bool> t;
	    t.push_back(value);
	    m_bool[key] = t;
	}
    

    //
    // Append a value onto a variable-length list
    //

    bool add_if_unique(const meta_name_t & name , const int value )
      { 
	meta_index_t midx = MetaInformation::field(name);
	if ( midx.len == -1 ) return add( midx.key , value , true );
	return false;
      }
    
    bool add_if_unique(const meta_name_t & name , const double value )
      { 
	meta_index_t midx = MetaInformation::field(name);
	if ( midx.len == -1 ) return add( midx.key , value , true );
	return false;
      }

    bool add_if_unique(const meta_name_t & name , const std::string & value )
      { 
	meta_index_t midx = MetaInformation::field(name);
	if ( midx.len == -1 ) return add( midx.key , value , true );
	return false;
      }

    bool add_if_unique(const meta_name_t & name , const bool value )
      { 
	meta_index_t midx = MetaInformation::field(name);
	if ( midx.len == -1 ) return add( midx.key , value , true );
	return false;
      }


    bool add(const meta_name_t & name , const int value , bool uniq = false ) 
      {
	meta_index_t midx = MetaInformation::field(name);
	if ( midx.len == -1 ) return add( midx.key , value , uniq );
	return false;
      }
    
    bool add(const meta_name_t & name , const double value , bool uniq = false  ) 
      {
	meta_index_t midx = MetaInformation::field(name);
	if ( midx.len == -1 ) return add( midx.key , value , uniq );
	return false;
      }
    
    bool add(const meta_name_t & name , const bool value , bool uniq = false  ) 
      {
	meta_index_t midx = MetaInformation::field(name);
	if ( midx.len == -1 ) return add( midx.key , value , uniq );
	return false;
      }
    
    bool add(const meta_name_t & name , const std::string & value , bool uniq = false  ) 
      {
	meta_index_t midx = MetaInformation::field(name);
	if ( midx.len == -1 ) return add( midx.key , value , uniq );
	return false;
      }
    
    //
    // Perform the actual appends, optional checking that this 
    // exact value does not already exist
    //
    
    bool add(const meta_key_t key , const int value , bool uniq )
      {
	if ( ! uniq ) 
	  {
	    m_int[key].push_back( value );
	    return true;
	  }

	if ( m_int.find( key ) != m_int.end() )
	  {
	    std::vector<int> & v = m_int[key];
	    for (unsigned int i=0; i<v.size(); i++) 
	      if ( v[i] == value ) return false;	    
	  }
	
	// okay to add
	m_int[key].push_back( value );
	return true;
      }

    
    bool add(const meta_key_t key , const double value , bool uniq  )
      {
	if ( ! uniq ) 
	  {
	    m_double[key].push_back( value );
	    return true;
	  }

	if ( m_double.find( key ) != m_double.end() )
	  {
	    std::vector<double> & v = m_double[key];
	    for (unsigned int i=0; i<v.size(); i++) 
	      if ( v[i] == value ) return false;	    
	  }
	
	// okay to add
	m_double[key].push_back( value );
	return true;
      }
    
    bool add(const meta_key_t key , const std::string & value , bool uniq  )
      {
	if ( ! uniq ) 
	  {
	    m_string[key].push_back( value );
 	    return true;
	  }

 	if ( m_string.find( key ) != m_string.end() )
	  {
	    std::vector<std::string> & v = m_string[key];
	    for (unsigned int i=0; i<v.size(); i++) 
	      if ( v[i] == value ) return false;
	  }
	
	// okay to add
	m_string[key].push_back( value );
	return true;
     }
    
    bool add(const meta_key_t key , const bool value , bool uniq  )
      {
	if ( ! uniq ) 
	  {
	    m_bool[key].push_back( value );
	    return true;
	  }

 	if ( m_bool.find( key ) != m_bool.end() )
	  {
	    std::vector<bool> & v = m_bool[key];
	    for (unsigned int i=0; i<v.size(); i++) 
	      if ( v[i] == value ) return false;	    
	  }
	
	// okay to add
	m_bool[key].push_back( value );
	return true;
      }
    

    //
    // Efficient set, with precalculated index list
    //
    
    bool set( const Helper::char_tok & tok , 
	      const std::vector<meta_index_t*> * midx )
	{
	    
	    // If last values missing, will be set to missing implicitly,
	    // but do not allow a longer-than-expected list
	    
	    if ( tok.size() > midx->size() ) return false;
	    
	    for ( int i = 0 ; i < tok.size(); i++ )
	    {
		
		if ( (*midx)[i] == NULL ) continue;
		
		meta_index_t * m = (*midx)[i];
		
		if ( m->len == 1 )  // single value fields
		{
		    switch( m->mt ) 
		    { 		  
			case META_INT   : 
			    int x;
			    if ( Helper::str2int( tok[i] , x ) ) set( m->key , x );
			    break;  		  
			case META_FLOAT : 
			    double d;
			    if ( Helper::str2dbl( tok[i] ,  d) ) set( m->key , d );
			    break;  		  
			case META_TEXT  : 
			    if ( *tok[i] != '.' ) set( m->key , tok[i] );
			    break;  		  
			case META_BOOL  : 
			    set( m->key , *tok[i] != '0' && *tok[i] != '.' && *tok[i] != 'F' );
			    break;  		  
		case META_FLAG  : 
		  set( m->key );
		  break;  
		}
		}
	  else
	    {
		// Assume comma-delimited vector
		int tok2size = 0;
		
		Helper::char_tok tok2( tok[i] , 0 , &tok2size , ',' );
		
		// If an exact length given, check we match
		if ( m->len > 1 && m->len != tok2.size() ) continue;
		// this implicitly handles the case in which the whole 
		// vector is specified as missing 
		
		// i.e. "."  instead of "0.0,0.1,0.9" for example (and not ".,.,.")
		
		switch( m->mt )
                {
		    case META_INT   :
		    {
			std::vector<int> x( tok2size );
			bool okay = true;
			for (int i=0;i<tok2size;i++) if ( ! Helper::str2int( tok2[i] , x[i] ) ) okay = false;
			if ( okay ) set( m->key , x );
			break;
		  }
                case META_FLOAT :
		  {
		      std::vector<double> x( tok2size );
		      bool okay = true;
		      for (int i=0;i<tok2size;i++) if ( ! Helper::str2dbl( tok2[i] , x[i] ) ) okay = false;
		      if ( okay ) set( m->key , x );
		      break;
		  }
		    case META_TEXT  :
		    {
			set( m->key , tok2 );
			break;
		    }
		    case META_BOOL  :
		    {
		    std::vector<bool> x( tok2size );
		    bool okay = true;
		    for (int i=0;i<tok2size;i++) 
		      {
			  if ( *tok2[i] == '0' || *tok2[i] == 'F' || *tok2[i] == 'f' || *tok2[i] == '.' ) x[i] = false;
			  else if ( *tok2[i] == '1' || *tok2[i] == 'T' || *tok2[i] == 't' ) x[i] = true;
			  else okay = false;
		      }
		    if ( okay ) set( m->key , x );
		    break;
		    }
                }	      
	    }
	}
    }
    
	      
	      
    //
    // Queries
    //
    
    bool hasField(const std::string & k) const
    { return has_field(k); }
    
    bool has_field(const std::string & k) const
    { 
      	meta_index_t midx = MetaInformation::field(k);
	return has_field( std::make_pair( midx.mt, midx.key ) );
    }

    bool has_field(const meta_typed_key_t & k) const
    {
                  
      if ( k.first == META_INT )
	return m_int.find( k.second ) != m_int.end();

      if ( k.first == META_FLOAT )
	return m_double.find( k.second ) != m_double.end();
      
      if ( k.first == META_TEXT )
	return m_string.find( k.second ) != m_string.end();
      
      if ( k.first == META_BOOL )
	return m_bool.find( k.second ) != m_bool.end();

      if ( k.first == META_FLAG )
	return m_flag.find( k.second ) != m_flag.end();
      
      return false;
    }
    
    
    //
    // Return number of elements in meta-field
    //
    

    int size( const std::string & k ) const
    {      
      meta_index_t midx = MetaInformation::field(k);
      return size( std::make_pair( midx.mt , midx.key ) );
    }
    
    
    int size( const meta_typed_key_t & k ) const
    {
      
      if ( k.first == META_FLAG ) return 0;
      
      if ( k.first == META_INT )
	{
	  meta_int_t::const_iterator i = m_int.find( k.second );
	  if ( i == m_int.end() ) return 0;
	  return i->second.size();
	}
      
      if ( k.first == META_FLOAT )
	{
	  meta_double_t::const_iterator i = m_double.find( k.second );
	  if ( i == m_double.end() ) return 0;
	  return i->second.size();
	}
      
      if ( k.first == META_TEXT )
	{
	  meta_string_t::const_iterator i = m_string.find( k.second );
	  if ( i == m_string.end() ) return 0;
	  return i->second.size();
	}
      
      if ( k.first == META_BOOL )
	{
	  meta_bool_t::const_iterator i = m_bool.find( k.second );
	  if ( i == m_bool.end() ) return 0;
	  return i->second.size();
	}
      
      return 0;
      
    }
    
    
    //
    // Get std::vector of elements associated to a key; type-strict
    //

    std::vector<std::string> get_flags() const
      {
	// Return all flags
	std::vector<std::string> t;
	meta_flag_t::const_iterator i1 = m_flag.begin();
	while ( i1 != m_flag.end() ) { t.push_back( field( META_FLAG, *i1 ) ); ++i1; }
	return t;
      }
    
    std::vector<std::string> get_string(meta_key_t key) const
      {
	    std::vector<std::string> t;
	    meta_string_t::const_iterator i = m_string.find( key );
	    return i == m_string.end() ? t : i->second; 
	}
    
    std::vector<int> get_int(meta_key_t key) const
	{
	    std::vector<int> t;
	    meta_int_t::const_iterator i = m_int.find( key );
	    return i == m_int.end() ? t : i->second; 
	}
    
    std::vector<double> get_double(meta_key_t key) const
	{
	    std::vector<double> t;
	    meta_double_t::const_iterator i = m_double.find( key );
	    return i == m_double.end() ? t : i->second; 
	}

    std::vector<bool> get_bool(meta_key_t key) const
	{
	    std::vector<bool> t;
	    meta_bool_t::const_iterator i = m_bool.find( key );
	    return i == m_bool.end() ? t : i->second; 
	}



    //
    // More efficient, pass pack only references to the originals
    //

    const std::vector<std::string> * ptr_string(meta_key_t key) const
      {
	meta_string_t::const_iterator i = m_string.find( key );
	return i == m_string.end() ? NULL : &(i->second);
      }
    
    const std::vector<int> * ptr_int(meta_key_t key) const
	{
	  meta_int_t::const_iterator i = m_int.find( key );
	  return i == m_int.end() ? NULL : &(i->second);
	}
    
    const std::vector<double> * ptr_double(meta_key_t key) const
      {
	meta_double_t::const_iterator i = m_double.find( key );
	return i == m_double.end() ? NULL : &(i->second); 
      }
    
    const std::vector<bool> * ptr_bool(meta_key_t key) const
      {
	meta_bool_t::const_iterator i = m_bool.find( key );
	return i == m_bool.end() ? NULL : &(i->second); 
      }
    



    //
    // Type specific queries, given name
    // 

    
    std::vector<std::string> get_string(const meta_name_t & name) const
	{
	    return get_string( MetaInformation::field( name ).key );
	}
    
    std::vector<int> get_int(const meta_name_t & name) const
	{
	    return get_int( MetaInformation::field( name ).key );
	}

    std::vector<double> get_double(const meta_name_t & name) const
	{
	    return get_double( MetaInformation::field( name ).key );
	}

    std::vector<bool> get_bool(const meta_name_t & name) const
	{
	    return get_bool( MetaInformation::field( name ).key );
	}


    //
    // Temporary measure: obtain first item from list
    //    i.e. assumes 1 and only 1 value

    std::string get1_string(const meta_name_t & name ) const 
	{
	    std::vector<std::string> t = get_string( MetaInformation::field( name ).key );
	    return t.size() > 0 ? t[0] : "";
	}
    
    int get1_int(const meta_name_t & name ) const 
	{
	    std::vector<int> t = get_int( MetaInformation::field( name ).key );
	    return t.size() > 0 ? t[0] : -1;
	}

    double get1_double(const meta_name_t & name ) const 
	{
	    std::vector<double> t = get_double( MetaInformation::field( name ).key );
	    return t.size() > 0 ? t[0] : -1;
	}
    
    bool get1_bool(const meta_name_t & name ) const 
	{
	    std::vector<bool> t = get_bool( MetaInformation::field( name ).key );
	    return t.size() > 0 ? t[0] : false;
	}

    
    //
    // Does a flag exist?
    //
    
    bool flag( const std::string & name ) const
    {
      // flag has not been registered
      if ( ! has_field( name ) ) return false;
      // does this particular obj possess flag?
      return flag( MetaInformation::field( name ).key ); 
    }
    

    bool flag( const meta_key_t k ) const 
    {
      meta_flag_t::const_iterator i1 = m_flag.begin();
      while ( i1 != m_flag.end() ) { if ( *i1 == k ) return true; ++i1; } 
      return false;
    }
    

    //
    // Convenience function, to return as a string list
    //
    
    std::string as_string(const meta_name_t & name , const std::string & delim = " " ) const 
      {

	meta_index_t midx = MetaInformation::field( name );

	if ( midx.mt == META_FLAG ) return flag( midx.key ) ? "set" : "";
	
	if ( midx.mt == META_INT )
	  {
	    std::stringstream s;
	    std::vector<int> t = get_int( midx.key );
	    for (int i=0; i<t.size(); i++)
	      {
		if ( i > 0 ) s << delim;
		s << t[i];
	      }		
	    return s.str();
	  }

	
	if ( midx.mt == META_FLOAT )
	  {
	    std::stringstream s;
	    std::vector<double> t = get_double( midx.key );
	    for (int i=0; i<t.size(); i++)
	      {
		if ( i > 0 ) s << delim;
		s << t[i];
	      }		
	    return s.str();

	  }

	if ( midx.mt == META_TEXT )
	  {
	    std::stringstream s;
	    std::vector<std::string> t = get_string( midx.key );
	    for (int i=0; i<t.size(); i++)
	      {
		if ( i > 0 ) s << delim;
		s << t[i];
	      }		
	    return s.str();

	  }

	if ( midx.mt == META_BOOL )
	  {
	    std::stringstream s;
	    std::vector<bool> t = get_bool( midx.key );
	    for (int i=0; i<t.size(); i++)
	      {
		if ( i > 0 ) s << delim;
		s << t[i];
	      }		
	    return s.str();

	  }

	return "";

      }

    
    //
    // Number of elements in total
    //

    int size() const 
    {
      return m_string.size() + m_int.size() + m_double.size() + m_bool.size() + m_flag.size();
    }
    

    //
    // Get a list of all keys
    //

    
    std::vector<std::string> keys() const
	{
	    std::vector<std::string> t;

	    meta_flag_t::const_iterator i1 = m_flag.begin();
	    while ( i1 != m_flag.end() ) { t.push_back( field( META_FLAG, *i1 ) ); ++i1; }
	    
	    meta_bool_t::const_iterator i2 = m_bool.begin();
	    while ( i2 != m_bool.end() ) { t.push_back( field( META_BOOL, i2->first ) ); ++i2; }
	    
	    meta_string_t::const_iterator i3 = m_string.begin();
	    while ( i3 != m_string.end() ) { t.push_back( field( META_TEXT, i3->first ) ); ++i3; }

	    meta_int_t::const_iterator i4 = m_int.begin();
	    while ( i4 != m_int.end() ) { t.push_back( field( META_INT, i4->first ) ); ++i4; }

	    meta_double_t::const_iterator i5 = m_double.begin();
	    while ( i5 != m_double.end() ) { t.push_back( field( META_FLOAT, i5->first ) ); ++i5; }

	    return t;
	}
    

    //
    // List all keys, but in meta_typed_key_t form
    //

    std::vector<meta_typed_key_t> typed_keys() const
	{
	  std::vector<meta_typed_key_t> t;

	  meta_flag_t::const_iterator i1 = m_flag.begin();
	  while ( i1 != m_flag.end() ) { t.push_back( std::make_pair( META_FLAG, *i1 ) ); ++i1; }
	  
	  meta_bool_t::const_iterator i2 = m_bool.begin();
	  while ( i2 != m_bool.end() ) { t.push_back( std::make_pair( META_BOOL, i2->first ) ); ++i2; }
	  
	  meta_string_t::const_iterator i3 = m_string.begin();
	  while ( i3 != m_string.end() ) { t.push_back( std::make_pair( META_TEXT, i3->first ) ); ++i3; }
	  
	  meta_int_t::const_iterator i4 = m_int.begin();
	  while ( i4 != m_int.end() ) { t.push_back( std::make_pair( META_INT, i4->first ) ); ++i4; }
	  
	  meta_double_t::const_iterator i5 = m_double.begin();
	  while ( i5 != m_double.end() ) { t.push_back( std::make_pair( META_FLOAT, i5->first ) ); ++i5; }
	  
	  return t;
	}
    
    

    template<class U> std::string print( const std::vector<U> & val) const 
	{
	    std::stringstream s;
	    for ( unsigned int i=0; i<val.size(); i++) 
		{
		    s << val[i];
		    if ( i != val.size() - 1 ) s << ",";
		}
	    return s.str();
	}
    

    std::string printValues( std::set<std::string> & k , std::string sep = ";" ) const 
      {
	std::stringstream out;	
	std::set<std::string>::iterator i = k.begin();
	while ( i != k.end() )
	  {
	    
	    if ( i != k.begin() ) out << sep;       
	    
	    const meta_index_t & midx = field( *i );
	    
	    if ( ! hasField( midx.name ) ) out << ".";
	    else if      ( midx.mt == META_INT )
	      out << print( m_int.find( midx.key )->second );
	    else if ( midx.mt == META_FLOAT )
	      out << print( m_double.find( midx.key )->second );
	    else if ( midx.mt == META_TEXT )
	      out << print( m_string.find( midx.key )->second );
	    else if ( midx.mt == META_BOOL )
	      out << print( m_bool.find( midx.key )->second );
	    else out << ".";	    

	    ++i;

	  }	
	return out.str();
      }

    std::string printValues(std::string sep = ";" ) const 
      {	    
	std::stringstream out;	
	for (unsigned int i = 0 ; i < ordered.size(); i++) 
	  {
	    const meta_index_t & midx = ordered[i];	    
	    if ( ! hasField( midx.name ) ) continue;	    
	    if      ( midx.mt == META_INT    && hasField( midx.name ) ) 
	      out << print( m_int.find( midx.key )->second );
	    else if ( midx.mt == META_FLOAT  && hasField( midx.name ) ) 
	      out << print( m_double.find( midx.key )->second );
	    else if ( midx.mt == META_TEXT   && hasField( midx.name ) ) 
	      out << print( m_string.find( midx.key )->second );
	    else if ( midx.mt == META_BOOL   && hasField( midx.name ) ) 
	      out << print( m_bool.find( midx.key )->second );
	    else out << ".";	    
	    if ( i != ordered.size() - 1 ) out << sep;       
	  }	
	return out.str();
      }
    
    
    //
    // Display this instance of meta-information
    //

    friend std::ostream & operator<<(std::ostream & out, const MetaInformation  & m ) 
	{ 

	  bool first = true;
	  
	  for (unsigned int i=0; i < MetaInformation::ordered.size(); i++)
	    {		
	      
	      meta_index_t midx = MetaInformation::ordered[i];
	      
	      if ( ! m.hasField( midx.name ) ) continue;
	      	      
	      if ( ! MetaMeta::display( midx.name ) ) continue;

	      if ( ! first ) out << ";";
	      first = false;
	      
	      out << midx.name;
	      
	      if      ( midx.mt == META_INT   ) out << "=" << m.print( m.m_int.find( midx.key )->second );
	      else if ( midx.mt == META_FLOAT ) out << "=" << m.print( m.m_double.find( midx.key )->second );
	      else if ( midx.mt == META_TEXT  ) out << "=" << m.print( m.m_string.find( midx.key )->second );
	      else if ( midx.mt == META_BOOL  ) out << "=" << m.print( m.m_bool.find( midx.key )->second );
	      //else if ( midx.mt != META_FLAG  ) out << "=.";
	      
	    }
	  
	  if ( first ) out << ".";

	  return out;
	}
    
    std::string print( meta_name_t & name ) 
      {
	meta_index_t midx = MetaInformation::field(name);
	if ( ! hasField( midx.name) ) return ".";
	if ( midx.mt == META_INT ) return print( m_int[ midx.key ] );
	if ( midx.mt == META_FLOAT ) return print( m_double[ midx.key ] );
	if ( midx.mt == META_TEXT ) return print( m_string[ midx.key ] );
	if ( midx.mt == META_BOOL ) return print( m_bool[ midx.key ] );
	if ( midx.mt == META_FLAG ) return "flag set";
	return ".";
      }

    std::string display(std::string indent = "\t") const
      { 

	std::string s;	
	for (unsigned int i=0; i< MetaInformation::ordered.size(); i++)
	  {

	    meta_index_t midx = MetaInformation::ordered[i];
	    if ( ! MetaMeta::display( midx.name ) ) continue;
	    if ( ! hasField( midx.name ) ) continue;

	    s += indent + midx.name;
	    if      ( midx.mt == META_INT   ) s += " = " + print( m_int.find( midx.key )->second ) + "\n";
	    else if ( midx.mt == META_FLOAT ) s += " = " + print( m_double.find( midx.key )->second ) + "\n";
	    else if ( midx.mt == META_TEXT  ) s += " = " + print( m_string.find( midx.key )->second ) + "\n";
	    else if ( midx.mt == META_BOOL  ) s += " = " + print( m_bool.find( midx.key )->second ) + "\n";
	    else if ( midx.mt == META_FLAG  ) s += " flag set\n";
	  }
	return s;
      }


    static std::string headers( int grp = 1 )
      {
	
	// filters are always displayed (as they will be in the VCF field)
	// as the straight string is printed

	const bool as_filter = grp == META_GROUP_FILTER ;
	const bool as_format = grp == META_GROUP_GEN ;
	
	std::stringstream s;	
	for (unsigned int i=0; i< MetaInformation::ordered.size(); i++)
	  {
	    meta_index_t midx = MetaInformation::ordered[i];

	    if ( grp != META_GROUP_FILTER && ! MetaMeta::display( midx.name ) ) continue;
	    
	    // VCF 4.0 headers
	    // ##INFO=<ID=ID,Number=number,Type=type,Description=description>
	    
	    // Or, if "as_filter=T" only ID and Description
	    
	    if ( as_format ) 
	      s << "##FORMAT=<ID=" << midx.name;
	    else if ( as_filter )
	      s << "##FILTER=<ID=" << midx.name;
	    else
	      s << "##INFO=<ID=" << midx.name;
	    
	    if ( ! as_filter ) 
	      {
		s << ",Number=" << midx.len 
		  << ",Type=";
		if ( midx.mt == META_INT ) s << "Integer";
		else if ( midx.mt == META_FLOAT ) s << "Float";
		else if ( midx.mt == META_TEXT ) s << "String";
		else if ( midx.mt == META_BOOL ) s << "Bool";
		else if ( midx.mt == META_FLAG ) s << "Flag";
		else s << "String";
	      }

	    s << ",Description=\"" << Helper::unquote( midx.description ) << "\">"
	      << "\n";
	  }
	return s.str();
      }    


    static std::string display_header(const std::string & prefix = "") 
      { 
	std::string s;	
	bool first = true;
	for (unsigned int i=0; i< MetaInformation::ordered.size(); i++)
	  {
	    meta_index_t midx = MetaInformation::ordered[i];
	    if ( ! MetaMeta::display( midx.name ) ) continue;
	    
	    if ( first ) { first = false; } 
	    else { s += "\t"; }

	    s += prefix + midx.name;

	  }
	return s;
      }
    
    static std::string display_header_static(const std::string & prefix = "") 
      { 
	std::string s;	
	bool first = true;
	for (unsigned int i=0; i< MetaInformation::ordered.size(); i++)
	  {
	    meta_index_t midx = MetaInformation::ordered[i];
	    if ( ! MetaMeta::display( midx.name ) ) continue;
	    if ( ! MetaMeta::static_variant( midx.name ) ) continue;
	    if ( first ) { first = false; } 
	    else { s += "\t"; }
	    s += prefix + midx.name;
	  }
	return s;
      }

    static std::string display_header_nonstatic(const std::string & prefix = "") 
      { 
	std::string s;	
	bool first = true;
	for (unsigned int i=0; i< MetaInformation::ordered.size(); i++)
	  {
	    meta_index_t midx = MetaInformation::ordered[i];
	    if ( ! MetaMeta::display( midx.name ) ) continue;
	    if ( MetaMeta::static_variant( midx.name ) ) continue;
	    if ( first ) { first = false; } 
	    else { s += "\t"; }
	    s += prefix + midx.name;
	  }
	return s;
      }

    std::string display_row( const std::string & missing = "." ) const
      { 
	std::string s;	
	bool first = true;
	for (unsigned int i=0; i< MetaInformation::ordered.size(); i++)
	  {		
	    meta_index_t midx = MetaInformation::ordered[i];
	    if ( ! MetaMeta::display( midx.name ) ) continue;	    

	    if ( first ) { first = false; } 
	    else { s += "\t"; }

	    if ( midx.mt == META_FLAG  ) 
	      {
		if ( hasField( midx.name ) ) 
		  s += "1"; 
		else 
		  s += "0";     	    
	      }
	    else if ( ! hasField( midx.name ) ) s += missing;
	    else if ( midx.mt == META_INT   ) s += print( m_int.find( midx.key )->second );
	    else if ( midx.mt == META_FLOAT ) s += print( m_double.find( midx.key )->second );
	    else if ( midx.mt == META_TEXT  ) s += print( m_string.find(  midx.key )->second );
	    else if ( midx.mt == META_BOOL  ) s += print( m_bool.find(  midx.key )->second );
	    
	  }	  
	return s;
      }

    // specialisation for static VarMeta attribs
    std::string display_row_static( const std::string & missing = "." ) const
      { 
	std::string s;	
	bool first = true;
	for (unsigned int i=0; i< MetaInformation::ordered.size(); i++)
	  {		
	    meta_index_t midx = MetaInformation::ordered[i];
	    if ( ! MetaMeta::display( midx.name ) ) continue;	    
	    if ( ! MetaMeta::static_variant( midx.name ) ) continue;	    
	    if ( first ) { first = false; } 
	    else { s += "\t"; }

	    if ( midx.mt == META_FLAG  ) 
	      {
		if ( hasField( midx.name ) ) 
		  s += "1"; 
		else 
		  s += "0";     	    
	      }
	    else if ( ! hasField( midx.name ) ) s += missing;
	    else if ( midx.mt == META_INT   ) s += print( m_int.find( midx.key )->second );
	    else if ( midx.mt == META_FLOAT ) s += print( m_double.find( midx.key )->second );
	    else if ( midx.mt == META_TEXT  ) s += print( m_string.find(  midx.key )->second );
	    else if ( midx.mt == META_BOOL  ) s += print( m_bool.find(  midx.key )->second );
	    
	  }	  
	return s;
      }

    // specialisation for static VarMeta attribs
    std::string display_row_nonstatic( const std::string & missing = "." ) const
      { 
	std::string s;	
	bool first = true;
	for (unsigned int i=0; i< MetaInformation::ordered.size(); i++)
	  {		
	    meta_index_t midx = MetaInformation::ordered[i];
	    if ( ! MetaMeta::display( midx.name ) ) continue;	    
	    if ( MetaMeta::static_variant( midx.name ) ) continue;	    
	    if ( first ) { first = false; } 
	    else { s += "\t"; }

	    if ( midx.mt == META_FLAG  ) 
	      {
		if ( hasField( midx.name ) ) 
		  s += "1"; 
		else 
		  s += "0";     	    
	      }
	    else if ( ! hasField( midx.name ) ) s += missing;
	    else if ( midx.mt == META_INT   ) s += print( m_int.find( midx.key )->second );
	    else if ( midx.mt == META_FLOAT ) s += print( m_double.find( midx.key )->second );
	    else if ( midx.mt == META_TEXT  ) s += print( m_string.find(  midx.key )->second );
	    else if ( midx.mt == META_BOOL  ) s += print( m_bool.find(  midx.key )->second );
	    
	  }	  
	return s;
      }


    
    // Append meta-information from one thing to annother, with an optional identifier
    template<class V>
      void append( V & m , const std::string & prefix = "" )
      {
	
	bool wp = prefix != "";
	
	std::vector<std::string> keys = m.keys();
	for (unsigned int i = 0; i < keys.size(); i++)
	  {

	    meta_index_t midx = V::field( keys[i] );
	    std::string name = wp ? prefix + "_" + midx.name : midx.name ;
	    
	    if ( ! exists( name ) )
	      {		
		field( name , midx.mt , midx.len , midx.description );		
	      }

	    // Assume same type otherwise....
	    // still uncertain how FLAGs should be handled...
	    
	    if ( midx.mt == META_TEXT ) 
	      set( name , m.get_string( keys[i] ) );
	    else if ( midx.mt == META_INT ) 
	      set( name , m.get_int( keys[i] ) );
	    else if ( midx.mt == META_FLOAT ) 
	      set( name , m.get_double( keys[i] ) );
	    else if ( midx.mt == META_BOOL ) 
	      set( name , m.get_bool( keys[i] ) );
	    else if ( midx.mt == META_FLAG ) 
	      {
		set( name );
	      }
	    else
	      set( name , m.get_string( keys[i] ) );
	    
	  }
      }


    
    // Reset just this instance of meta-information
    
    void clear()
	{
	    m_string.clear();
	    m_int.clear();
	    m_double.clear();
	    m_bool.clear();
	}
       

    
    //
    // Static functions
    //

    static std::vector<std::string> field_names()
    { 
      std::vector<std::string> s;	
      for (unsigned int i=0; i< MetaInformation::ordered.size(); i++)
	{
	  meta_index_t midx = MetaInformation::ordered[i];
	  s.push_back( midx.name );	  
	}
      return s;
    }

    static int n_keys()
    {
      return ordered.size();
    }
    
    static int n_visible_keys()
      {
	int n = 0;
	for ( unsigned int i=0; i < ordered.size(); i++)
	  if ( MetaMeta::display( ordered[i].name ) ) ++n;
	return n;	
      }

    // specialization for Variant VarMeta
    static int n_visible_keys_static()
      {
	int n = 0;
	for ( unsigned int i=0; i < ordered.size(); i++)
	  if ( MetaMeta::display( ordered[i].name ) && MetaMeta::static_variant( ordered[i].name) ) ++n;
	return n;	
      }

    // specialization for SampleVariant/Consensus VarMeta
    static int n_visible_keys_nonstatic()
      {
	int n = 0;
	for ( unsigned int i=0; i < ordered.size(); i++)
	  if ( MetaMeta::display( ordered[i].name ) && ! MetaMeta::static_variant( ordered[i].name ) ) ++n;
	return n;
      }

    static mType type(const meta_name_t & s) 
      {
	std::map<meta_name_t,meta_index_t>::iterator i = nameMap.find( s );
	return i == nameMap.end() ? META_UNDEFINED : i->second.mt;
      }
    
    static meta_index_t * index(const meta_name_t & s)
      {
	std::map<meta_name_t,meta_index_t>::iterator i = nameMap.find( s );
	return i == nameMap.end() ? NULL : &(i->second);
      }
    
    static meta_name_t field( const mType mt , const meta_key_t & key )
      {
	meta_index_t midx;
	midx.mt = mt;
	midx.key = key;
	std::set<meta_index_t>::iterator i = indexSet.find( midx );
	return i == indexSet.end() ? "" : i->name ;
      }
    
    static bool exists( const meta_name_t & s )
      { 
	std::map<meta_name_t,meta_index_t>::iterator i = nameMap.find( s );
	return i != nameMap.end();
      }


    static std::map<meta_name_t, meta_index_t> dump_types()
	{
	    return nameMap;
	}



    static meta_index_t field(const meta_name_t & s,
			      const mType t = META_UNDEFINED ,
			      const int num = -1 ,
			      const std::string desc = "" )
      {

	// this may be a system(hidden) flag, but that gets declared after
	// this initial declaration, so treat everything as external for now
	
	MetaMeta::is_external( s );


	std::map<meta_name_t,meta_index_t>::iterator i = nameMap.find( s );

	if ( i == nameMap.end() )
	  {

	    meta_index_t midx; 

	    midx.mt = t;

	    if      ( t == META_INT   ) midx.key = cnt_int++; 
	    else if ( t == META_FLOAT ) midx.key = cnt_double++;
	    else if ( t == META_BOOL  ) midx.key = cnt_bool++;
	    else if ( t == META_FLAG  ) midx.key = cnt_flag++;
	    else  // default to string if none of above
	      {
		midx.mt = META_TEXT;
		midx.key = cnt_string++;
	      }
	    
	    midx.name = s;
	    midx.description = desc;
	    midx.len = num;

	    nameMap[ midx.name ] = midx;
	    indexSet.insert( midx );	    
	    ordered.push_back( midx );
	    	    
	    return midx;
	  }
	else
	  {
	    return i->second;
	  }
      }
    
    
    static std::string list_fields( const std::string & t )
    {

      std::stringstream ss;
      
      std::map< meta_name_t , meta_index_t >::iterator i = nameMap.begin();
      
      while ( i != nameMap.end() )
	{
	  
	  meta_index_t & midx = i->second;
	  
	  ss << t << "\t" 
	     << "NAME=" << midx.name << "\t";
	  
	  if ( MetaMeta::display( midx.name ) ) 
	    ss << "DISPLAY=Y\t";
	  else
	    ss << "DISPLAY=N\t";
	  
	  switch ( midx.mt ) {
	  case META_INT:
	    ss << "TYPE=Integer\t";
	    break;
	  case META_FLOAT :
	    ss << "TYPE=Float\t";
	    break;					 
	  case META_TEXT :
	    ss << "TYPE=String\t";
	    break;
	  case META_BOOL :
	    ss << "TYPE=Bool\t";
	    break;
	  case META_FLAG :
	    ss << "TYPE=Flag\t";
	    break;
	  default :
	    ss << "TYPE=Undefined\t";
	  }		
	  
	  ss << "LEN=" << midx.len << "\t"
	     << "DESC=" << midx.description << "\n";
	  
	  ++i;
	}
      return ss.str();
    }
    

    static std::string pretty_list_fields( const std::string & t )
    {
	
	std::stringstream ss;
	
	std::map< meta_name_t , meta_index_t >::iterator i = nameMap.begin();
	
	while ( i != nameMap.end() )
	{
	    
	    meta_index_t & midx = i->second;
	    
	    if ( MetaMeta::display( midx.name ) ) 
	    {
		
		ss << midx.name << " : "
		   << midx.description << " (" 
		   << t << ", ";
		
		switch ( midx.mt ) {
		    case META_INT:
			ss << "Integer";
			break;
		    case META_FLOAT :
			ss << "Float";
			break;					 
		    case META_TEXT :
			ss << "String";
		      break;
		    case META_BOOL :
			ss << "Bool";
			break;
		    case META_FLAG :
			ss << "Flag";
			break;
		    default :
			ss << "Undef.";
		}
		
		if ( midx.len > 1 ) 
		    ss << " x " << midx.len;
		else if ( midx.len == -1 ) 
		    ss << " variable-length vector";
		
		ss << ")\n"; 
	    }
	    
	    ++i;
	}
	return ss.str();
    }

    
    static void reset()
	{
	    nameMap.clear();
	ordered.clear();
	cnt_int = cnt_double = cnt_string = cnt_bool = cnt_flag = 0;
      }
    
    
 private:
    
    //
    // Members
    //

    meta_string_t     m_string;
    meta_int_t        m_int;
    meta_double_t     m_double;
    meta_bool_t       m_bool;
    meta_flag_t       m_flag;

    //
    // Class of meta-information
    //

    static T c;
    

    //
    // Static members (track type information)
    //
     
    static std::map<meta_name_t, meta_index_t> nameMap;

    static std::set<meta_index_t> indexSet;

    static std::vector<meta_index_t> ordered;

    static int cnt_int;
    static int cnt_double;
    static int cnt_string;
    static int cnt_bool;
    static int cnt_flag;

};

template<class T>  std::map<meta_name_t,meta_index_t> MetaInformation<T>::nameMap;
template<class T>  std::set<meta_index_t>             MetaInformation<T>::indexSet;
template<class T>  std::vector<meta_index_t>          MetaInformation<T>::ordered;

template<class T>  int MetaInformation<T>::cnt_int    = 0;
template<class T>  int MetaInformation<T>::cnt_bool   = 0;
template<class T>  int MetaInformation<T>::cnt_string = 0;
template<class T>  int MetaInformation<T>::cnt_double = 0;
template<class T>  int MetaInformation<T>::cnt_flag   = 0;
template<class T>  T   MetaInformation<T>::c;


#endif
