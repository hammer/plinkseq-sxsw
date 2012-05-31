#ifndef __REFERENCE_H__
#define __REFERENCE_H__

#include <string>
#include <vector>

#include "sqlwrap.h"
#include "regions.h"

#include "helper.h"
#include "variant.h"

class Region;

struct RefVarInfo {
  
  RefVarInfo()
  {
    build = date = version = name = description = "N/A";
  }
  
  std::string name;
  std::string file;
  std::string date;
  std::string version;
  std::string description;
  std::string build;
  
  friend std::ostream & operator<<( std::ostream & out , const RefVarInfo & r)
  {
    out << "  NAME:  " << r.name << "\n"
	<< "  FILE: " << r.file << "\n"
	<< "  DATE: " << r.date << "\n"
	<< "   VER: " << r.version << "\n"
	<< "  DESC: " << r.description << "\n"
	<< " BUILD: " << r.build << "\n";
    return out;
  }
  
};


class RefVariant {
  
 public:
  
  RefVariant()
    {
      grp_id = 0;
      rname = "";
      chr = 0;
      bp1 = bp2 = 0;
      ref = alt = ".";
      rvalue = "";
      obs = false;
    }
  
  RefVariant(int grp_id, 
	     const std::string & name, 
	     int chr, int bp1, int bp2, 
	     const std::string & ref , 
	     const std::string & alt , 
	     const std::string & value)
    : grp_id(grp_id), rname(name), chr(chr) , 
    bp1(bp1) , bp2(bp2), ref(ref) , alt(alt) , 
    rvalue(value) 
    {	    
      obs = true;
    }
    
    bool operator<( const RefVariant & rhs ) const 
    {
      if ( chr < rhs.chr ) return true;
      if ( chr > rhs.chr ) return false;
      
      if ( bp1 < rhs.bp1 ) return true;
      if ( bp1 > rhs.bp1 ) return false;
      
      if ( bp2 < rhs.bp2 ) return true;
      if ( bp2 > rhs.bp2 ) return false;
      
      return alt < rhs.alt;

    }
    
  std::string name() const { return rname; }
  void name(const std::string n) { rname = n; }
  
  std::string value() const { return rvalue; }
  void value(const std::string & s) { rvalue = s; }
  
  bool observed() const { return obs; }
  void observed(const bool b) { obs=b; }
  
  int group() const { return grp_id; }
  void group(const int g) { grp_id = g; }
  
  int chromosome() const { return chr; }
  void chromosome(const int c) { chr = c; }

  int start() const { return bp1; }
  void start(const int p) { bp1 = p; }
  
  int stop() const { return bp2; }
  void stop(const int p) { bp2 = p; }
  
  std::string reference() const { return ref; }
  void reference(const std::string & r ) { ref = r; } 

  std::string alternate() const { return alt; }
  void alternate(const std::string & a ) { alt = a; } 
  
  // But also allow more extensible meta-information
  
  MetaInformation<RefMeta> meta;
  
  friend std::ostream & operator<<( std::ostream & out , const RefVariant & rv)
    {
      if ( rv.obs )
	{
	  out << Helper::chrCode( rv.chr ) << ":"
	      << rv.bp1;
	  if ( rv.bp2 > rv.bp1 ) 
	    out << ".." << rv.bp2;	  
	  out << ":" << rv.ref << "/" << rv.alt << ":" << rv.rname;
	}	  
      else
	out << ".";
      return out;
    }
  
 private:
  
    int grp_id;
    std::string rname;
    int chr;
    int bp1, bp2;
    std::string ref;
    std::string alt;
    std::string rvalue;
    bool obs;
    
};



class RefDBase { 
    
 public:
    
    RefDBase()
	{
	  values = false;
	  vstring = "";
	}
    
    ~RefDBase()
      {
	dettach();
	sql.close();
      }
    
    
    //
    // Core functions to connect, close to database
    //
    
    bool attach(std::string name);
    bool init();
    bool release();
    bool dettach();
    
    void index();
    void drop_index();

    bool attached() { return sql.is_open(); }
    
    std::string summary(bool);
    
    std::vector<std::string> fetch_groups();
    
    void set_metatypes( bool clear = false);

    // Main insertion functions

    bool load_VCF( const std::string & file,
		   const std::string & tag ,
		   const std::string & comment,
		   const std::set<std::string> & includes,
		   const std::set<std::string> & excludes,
		   const std::set<Region> * pfilter );
    
    // these are now redundant...
    void set_values( const std::string & s ) { values = true; vstring = s; } 
    bool values;
    std::string vstring;

    std::map<std::string,mType> populate_metatypes( std::map<std::string,int> * meta , const int );

    std::map< std::string, std::map<std::string,std::string> > get_metatypes();

    void attach_metainformation( RefVariant & , const Variant & );
				
    void attach_metainformation( RefVariant & rv , 
				 std::map<std::string,int> * meta , 
				 std::set<int> * flags , 
				 std::set<int> * skip , 
				 std::map<std::string,mType> & , 
				 std::vector<std::string> & );
    
    uint64_t loadRefVariants(const std::string & filename , 
			     const std::string & grp, 
			     int col_name,
			     int col_chr,
			     int col_bp1,
			     std::map<std::string,int> * meta = NULL , 
			     std::set<int> * flags = NULL , 
			     std::set<int> * skip = NULL , 
			     int col_bp2 = -1 ,
			     bool zero1 = false , 
			     bool zero2 = false );

    // Core inserts
    
    uint64_t insert( const std::string & filename , const std::string & grp );

    void insert_header( uint64_t file_id , const std::string & name , const std::string & value );
    
    // old void insert_metatype( const std::string & filename , const std::string & grp );
    void insert_metatype( uint64_t file_id , const std::string & name ,  mType mt , int num, int mgrp,  const std::string & description );

    // Insert from Variant (->RefVariant)
    bool     insert( int , const Variant & );

    // Directly insert RefVariant
    bool     refInsertion(const RefVariant &);
        
    
    // Core queries
    
    RefVariant lookup( const Variant & , const int grp_id );
    RefVariant lookup( const Variant & , const std::string &  );

    std::set<RefVariant> lookup( const Region & , const int , const int limit = 0 );
    std::set<RefVariant> lookup( const Region & , const std::string & , const int limit = 0 );

    bool annotate( Variant & , const std::string & );
    bool annotate( Variant & , const int grp_id );
    
    int count( const Region & , const std::string & );

    // Data-dumper
    
    void dump( const std::string & grp , bool with_meta );

    // Groups
    
    uint64_t set_group_id(const std::string grp, 
			  const bool temp = false , 
			  const std::string desc = "");
    
    uint64_t lookup_group_id(const std::string grp);
    
    void flush(const uint64_t id = 0);
    
    // Iterator

    bool init_iterate( const std::string & );
    bool iterate( RefVariant * );

    // Misc
      
    std::string filename() const { return fname; }


 private:
    
    SQL sql;
    std::string fname;

    RefVariant construct( sqlite3_stmt * );
    void construct_inplace( sqlite3_stmt * s , RefVariant * rv );
    
    void check_version();
    
    std::map<std::string,meta_key_t> mtmap;
    std::map<int,std::string> grpmap;
      
    sqlite3_stmt * stmt_dump;
    sqlite3_stmt * stmt_lookup;
    sqlite3_stmt * stmt_lookup_range;
    sqlite3_stmt * stmt_lookup_dbsnp;
    sqlite3_stmt * stmt_lookup_range_count;

    sqlite3_stmt * stmt_insert;
    sqlite3_stmt * stmt_insert_group;

    sqlite3_stmt * stmt_meta_dump;
    sqlite3_stmt * stmt_meta_insert_prep;
    sqlite3_stmt * stmt_meta_insert_prep2;
    sqlite3_stmt * stmt_meta_insert;
    sqlite3_stmt * stmt_get_meta;
    sqlite3_stmt * stmt_fetch_metatypes;
    sqlite3_stmt * stmt_fetch_groups;
    
    sqlite3_stmt * stmt_lookup_group_name;
    sqlite3_stmt * stmt_insert_group_name;
    sqlite3_stmt * stmt_remove_group1;
    sqlite3_stmt * stmt_remove_group2;
    sqlite3_stmt * stmt_update_group_count;
    sqlite3_stmt * stmt_fetch_group_count;

};


#endif

