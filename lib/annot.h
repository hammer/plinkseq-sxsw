#ifndef __PLINKSEQ_ANNOT_H__
#define __PLINKSEQ_ANNOT_H__

#include <string>
#include <vector>

#include "regions.h"
#include "helper.h"
#include "variant.h"
#include "locdb.h"


class SeqDBase;
class RefDBase;


enum seq_annot_t { UNDEF   =  0 ,     // could not annotate 
                   MONO    =  1 ,     // monomorphic site 

		   IGR     =  2 ,     // intergenic region

		   // near gene?
		   INTRON  =  3 ,     // intronic		   
		   UTR5    =  4 ,     // 5' UTR allele -- not used
		   UTR3    =  5 ,     // 3' UTR allele -- not used		   
		   
		   // exonic
		   SYN      =  10 ,    // synonymous allele 		   		   
		   
		   // non-synon coding
		   MIS      =  20 ,    // missense allele
 		   PART     =  21 ,    // partial codon  -- not used
		   SPLICE5  =  22 ,    // 5' splice-site
		   SPLICE3  =  23 ,    // 3' splice-site 
		   ESPLICE5 =  27 ,    // Essential 5' splice-site
                   ESPLICE3 =  28 ,    // Essential 3' splice-site
 		   NON      =  24 ,    // nonsense allele		   		  
 		   FS       =  25 ,    // frameshift 
		   RT       =  26 };   // readthrough

struct SeqInfo { 
  
  // note -- these function depend on exact coding of seq_annot_t (see above)

  bool missense() const { return type == 20 ; } 
  bool nonsense() const { return type == 24 ; }
  bool readthrough() const { return type == 26 ; }
  bool frameshift() const { return type == 25 ; }    
  bool splice() const { return type == 22 || type == 23; }   
  bool esplice() const { return type == 27 || type == 28; }

  bool coding() const { return type > 9 ; } 
  bool synon() const { return type == 10 ; } 
  bool nonsyn() const { return type > 19 ; }   
  bool intergenic() const { return type == 2 ; }
  bool intronic() const { return type == 3 ; }
  bool invalid() const { return type < 2 ; }
  
  static std::map< seq_annot_t , std::string> types;
  
  SeqInfo( seq_annot_t t ) : type(t) 
  { 
    genomic_ref = genomic_alt = "";
    cpos1 = 0;
    ref_seq = alt_seq = "";
    ppos1 = 0;
    ref_aa = alt_aa = "";
    transcript = "";
  } 
  
  SeqInfo( const std::string & transcript, 
	   const seq_annot_t & type , 
	   const std::string & genomic_ref = "", 
	   const std::string & genomic_alt = "" ,
	   const int cpos1 = 0, 
	   const std::string & ref_seq = "",
	   const std::string & alt_seq = "",
	   const int ppos1 = 0, 
	   const std::string & ref_aa = "",
	   const std::string & alt_aa = "" )
    : transcript(transcript), type(type), 
      genomic_ref(genomic_ref) , genomic_alt(genomic_alt),
      ref_seq(ref_seq), ref_aa(ref_aa), 
      alt_seq(alt_seq), alt_aa(alt_aa),
      cpos1(cpos1), cpos2(cpos1), 
      ppos1(ppos1), ppos2(ppos1) 
  {
    splicedist = 0;
  } 
    
  bool operator<( const SeqInfo & rhs ) const
  {
    if ( transcript < rhs.transcript ) return true;
    if ( transcript > rhs.transcript ) return false;
    if ( type < rhs.type ) return true;
    if ( type > rhs.type ) return false;
    return genomic_alt < rhs.genomic_alt;
  }
  
  
  seq_annot_t type;
  
  std::string transcript;

  int splicedist; // for splice-sites only
  
  int cpos1;  // position in DNA
  int cpos2;

  int ppos1;  // position in protein
  int ppos2;

  std::string genomic_ref;
  std::string genomic_alt;

  std::string ref_seq;
  std::string ref_aa;
  
  std::string alt_seq;
  std::string alt_aa;  
  

  std::string genomic() const;

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
    return transcript + ":" + status() + ":" + codon() + ":" + protein(); 
  }
  
  std::string gene_name() const { return transcript == "" ? "." : transcript ; }
  
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

    static bool annotate(Variant & var , Region * pregion = NULL );

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

