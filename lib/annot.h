#ifndef __PLINKSEQ_ANNOT_H__
#define __PLINKSEQ_ANNOT_H__

#include <string>
#include <vector>

#include "regions.h"
#include "options.h"
#include "helper.h"
#include "variant.h"
#include "locdb.h"


class SeqDBase;
class RefDBase;


enum seq_annot_t { UNDEF   =  0 ,     // could not annotate 
                   MONO    =  1 ,     // monomorphic site 

		   IGR     =  2 ,     // intergenic region
		   INTRON  =  3 ,     // intronic		   
		   UTR5    =  4 ,     // 5' UTR allele -- not used
		   UTR3    =  5 ,     // 3' UTR allele -- not used		   

 		   SYN     =  10 ,    // synonymous allele 		   

		   MIS     =  20 ,    // missense allele
 		   PART    =  21 ,    // partial codon  -- not used
		   SPLICE  =  22 ,    // split-site allele -- not used
 		   NON     =  23 ,    // nonsense allele		   		  
 		   FS      =  24 ,    // frameshift  -- not used
		   RT      =  25 };   // readthrough


struct SeqInfo { 

  static std::map< seq_annot_t , std::string> types;
  
  SeqInfo( seq_annot_t t ) : type(t) { } 
  
  SeqInfo( const std::string & transcript, 
	   const seq_annot_t & type , 
	   const int cpos1 = 0, 
	   const std::string & ref_seq = "",
	   const std::string & alt_seq = "",
	   const int ppos1 = 0, 
	   const std::string & ref_aa = "",
	   const std::string & alt_aa = "" )
    : transcript(transcript), type(type), 
      ref_seq(ref_seq), ref_aa(ref_aa), 
      alt_seq(alt_seq), alt_aa(alt_aa),
      cpos1(cpos1), cpos2(cpos1), 
      ppos1(ppos1), ppos2(ppos1) 
  { } 
  
  
  bool operator<( const SeqInfo & rhs ) const
  {
    if ( transcript < rhs.transcript ) return true;
    if ( transcript > rhs.transcript ) return false;
    return type < rhs.type;
  }
  
  
  seq_annot_t type;
  
  std::string transcript;

  int cpos1;  // position in DNA
  int cpos2;

  int ppos1;  // position in protein
  int ppos2;

  std::string ref_seq;
  std::string ref_aa;
  
  std::string alt_seq;
  std::string alt_aa;  
  
  // note -- these function depend on exact coding of seq_annot_t (see above)
  bool coding() const { return type > 19 ; } 
  bool exonic() const { return type > 9 ; }   
  
  std::string codon() const;
  
  std::string protein() const;
  
  std::string status() const 
  { 
    std::map<seq_annot_t,std::string>::iterator i = types.find(type);
    if ( i == types.end() ) return ".";
    return i->second;
  }

  std::string summary() const 
  {

  }

  std::string gene_name() const { return transcript; }
  
};



class Annotate {
    
    // Helper functions

    static std::string getrc(const std::string &);

    static std::string translate(std::string &, int, std::vector<std::string> &);
    
    // DNA base --> AA
    static std::map<std::string,std::string> t;


    // Pointers to data-base functions

    static SeqDBase * seqdb;

    static LocDBase * db;

    // Hold region map, populated by fetch_transcripts()

    static std::map<uint64_t,Region> rmap;
    static uint64_t transcript_group_id;

    static bool load_transcripts( uint64_t id );

    static bool load_transcripts( uint64_t id , const std::set<Region> & );
    
 public:
    
    // AA names
    static std::map<std::string,std::string> aa;

    static void setDB( const fType t );

    static void setDB( LocDBase * , SeqDBase * );

    static void init();
    
    static void clear()
    {
      transcript_group_id = 0;
      rmap.clear();
    }
    
    
    // Add as meta-information any interesting annotation, returning 'T' if
    // anything was added

    static bool load_transcripts( fType d , const std::string & );

    
    static bool annotate(Variant & var);

    static std::set<SeqInfo> lookup(Variant & var);
    
    static std::set<SeqInfo> annotate( int chr, 
				       int bp1 , 
				       const std::string & alternate , 
				       const std::string & reference , 
				       const Region & r );

    static std::set<SeqInfo> annotate( int,int,
				       const std::string & alt,
				       const std::string & ref = "", 
				       const std::set<Region> * pregion = NULL );

   
    static std::string translate_reference( const Region & region , bool verbose = false);
    
};

#endif

