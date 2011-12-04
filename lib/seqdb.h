#ifndef __PLINKSEQ_SEQDB_H__
#define __PLINKSEQ_SEQDB_H__

#include <string>
#include <vector>

#include "sqlwrap.h"
#include "regions.h"

#include "helper.h"
#include "locdb.h"

enum seqdb_repeat_mode_t {
  SEQ_RPT_UNKNOWN = 0 , 
  SEQ_RPT_NONE = 1 , 
  SEQ_RPT_LOWER = 2 , 
  SEQ_RPT_N = 3
};

class SeqDBase { 
    
 public:
    
    SeqDBase() 
	{
	    cache_chr = 0;
	    cache_bp1 = cache_bp2 = 0;
	    
	    rpt_mode = SEQ_RPT_UNKNOWN;
	    use_iupac = false;
	    genome_build = "";
	    name = desc = "";
	}
    
    ~SeqDBase()
	{
	    dettach();
	    sql.close();
	}
    
    
    //
    // Core functions to connect, close to database
    //
    
    bool create( const std::string & );
    bool attach( const std::string & name);
    bool init();
    bool release();
    bool dettach();

    void index();

    bool attached() { return sql.is_open(); }
    
    std::string summary(bool);


    void setMinMax();
    
    std::map<int,int2> getMinMax() const { return chrminmax; }

    //
    // Main insertion functions
    //
    
    void  loadFASTA( const std::string & filename , const std::map< std::string , std::string > & );

    void  insert(const int,const int,const int,const std::string &);
    
    //
    // Core queries
    //

    std::string lookup( const Variant & , int );
    std::string lookup( const Region & region );
    std::string lookup( int chr , int bp1 , int bp2 = 0 );
 
    bool N( const Region & region , int & n , int & tot );   
    bool GC( const Region & region , int & gc , int & tot );
    bool ACGT( const Region & region , int & a , int & c , int & g , int & t , int & n );

    void dump( const Region & , bool compact = false );
    
    //
    // Helper functions
    //
        
    static bool iupac( const std::string & , const std::string & );

    std::map< std::string , std::string > lookup_meta();

    void insert_meta( const std::map< std::string , std::string > & );

 private:

    //
    // Main datastore
    //

    SQL sql;


    //
    // Meta variables
    //

    seqdb_repeat_mode_t rpt_mode;

    std::string name ;

    std::string desc ;

    std::string  genome_build;

    bool use_iupac;

    std::map<int,int2> chrminmax;
    
    std::map<std::string,std::string> meta;

    //
    // Cache
    //

    static const unsigned int BATCH_SIZE = 500000;
    
    std::string chunk;

    int cache_chr, cache_bp1, cache_bp2;


    //
    // SQL prepared queries
    //
    
    sqlite3_stmt * stmt_insert;
    sqlite3_stmt * stmt_lookup;

    sqlite3_stmt * stmt_getmeta;
    sqlite3_stmt * stmt_putmeta;
        
};



#endif

