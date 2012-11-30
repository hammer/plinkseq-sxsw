#ifndef __PLINKSEQ_ANNOT_H__
#define __PLINKSEQ_ANNOT_H__

#include <string>
#include <vector>

#include "regions.h"
#include "plinkseq/helper.h"
#include "plinkseq/variant.h"
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
		   INDEL    = 11 , // any indel
		   // non-synon coding
		   MIS      =  20 ,    // missense allele
 		   PART     =  21 ,    // partial codon  -- not used
 		   CODONINSERTION		=  22 ,	   // codon insertion
 		   CODONDELETION		=  23 ,    // codon deletion
		   STOPINSERTION               =  34 ,    // stop insertion                                                                                                                 
                   STOPDELETION                =  35 ,    // stop deletion
		   OOFCODONINSERTION               =  36 ,    // out of frame codon insertion
		   OOFCODONDELETION               =  37 , // out of frame codon deletion
 		   SPLICE 	= 24	, // general splice +/- 5bp
 		   EXONIC_UNKNOWN = 38, // overlaps an exon, but since the ALT is 'N', cannot know its exact coding impact
 		   
 		   // Special class of splice variants : Faustino and Cooper. Pre-mrna splicing and human disease. AG|G   AG|GTNAG. This is consistent with splicing motif measures.
		   DONORIN2  =  25 ,       // donor splice-site |[GT]
		   DONOREX2AG = 26 ,       // donor splice-site ex2ag [AG]|
		   DONORIN45AG = 27 ,      // donor splice-site in45ag |GTN[AG]
		   ACCEPTOREX1G = 28 ,     // acceptor splice-site ex1g |[G]
		   ACCEPTORIN2  =  29 ,    // 3' splice-site [AG]|
		   
		   // Additional LoF annotations 
		   SL = 30 ,	  // Start Loss
 		   NON      =  31 ,    // nonsense allele
 		   FRAMESHIFT       =  32 ,    // frameshift
		   RT       =  33 };   // readthrough

struct SeqInfo { 
  
  // note -- these function depend on exact coding of seq_annot_t (see above)

  bool missense() const { return type == 20 ; } 
  bool nonsense() const { return type == 31 ; }
  bool startlost() const { return type == 30 ; }
  bool readthrough() const { return type == 33 ; }
  bool frameshift() const { return type == 32 ; }
  bool codondeletion() const { return type == 23 || type == 35 || type == 37; }
  bool codoninsertion() const { return type == 22 || type == 34 || type == 36;}
  bool splice() const { return type == 24 ;}
  bool csplice() const { return type == 26 || type == 27 || type == 28 ; }
  bool esplice() const { return type == 25 || type == 29; }

  bool coding() const { return type > 9 ; } 
  bool synon() const { return type == 10 ; } 
  bool indel() const { return type == 11 ; }
  bool nonsyn() const { return type > 19 ; }   
  bool intergenic() const { return type == 2 ; }
  bool intronic() const { return type == 3 ; }
  bool exonic_unknown() const { return type == EXONIC_UNKNOWN; }
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
    fs_stop = 0;
    origpepsize = 0;
    newpepsize = 0;
    exon=0;
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
	   const std::string & alt_aa = "" ,
	   const int fs_stop = 0 ,
	   const int origpepsize = 0,
	   const int newpepsize = 0, 
	   const int exon = 0 )
    : transcript(transcript), type(type), 
      genomic_ref(genomic_ref) , genomic_alt(genomic_alt),
      ref_seq(ref_seq), ref_aa(ref_aa), 
      alt_seq(alt_seq), alt_aa(alt_aa),
      cpos1(cpos1), cpos2(cpos1), 
      ppos1(ppos1), ppos2(ppos1),
      fs_stop(fs_stop),
      origpepsize(origpepsize),
      newpepsize(newpepsize), 
      exon(exon)
  {
    splicedist = 0;
    nmd = 0;
    ofptv = 0;
    exin = 0;
    iseq = "";
    eseq = "";
    alt = "";
    splice_type = "";
  }
    
  bool operator<( const SeqInfo & rhs ) const
  {
    if ( transcript != rhs.transcript ) return transcript < rhs.transcript;
    if ( cpos1 != rhs.cpos1 ) return cpos1 < rhs.cpos1;
    if ( cpos2 != rhs.cpos2 ) return cpos2 < rhs.cpos2;
    if ( genomic_alt != rhs.genomic_alt ) return genomic_alt < rhs.genomic_alt;
    if ( ppos1 != rhs.ppos1 ) return ppos1 < rhs.ppos1;
    if ( ppos2 != rhs.ppos2 ) return ppos2 < rhs.ppos2;
    /* if ( type != rhs.type ) */ return type < rhs.type;
  }
  
  
  seq_annot_t type;
  
  std::string transcript;

  int splicedist; // for splice-sites only
  int ofptv; // for splice-sites only at the moment. If out of frame protein truncating variant. 
  int nmd; // for LoF variants 
  int exin; // for splice variants -- what exon is closest 
  int cpos1;  // position in DNA
  int cpos2;

  std::string iseq; // intronic splice sequence
  std::string eseq; // exonic splice sequence
  std::string splice_type; // splice donor or acceptor
  std::string alt; // alternate allele coded properly to match ref

  int ppos1;  // position in protein
  int ppos2;

  int fs_stop; // number of AA positions to new stop
  int origpepsize; //original pep size
  int newpepsize; // new pep size
  int exon; //exon number of mutation
  
  std::string genomic_ref;
  std::string genomic_alt;

  std::string ref_seq;
  std::string ref_aa;
  
  std::string alt_seq;
  std::string alt_aa;  
  

  std::string genomic() const;

  std::string codon() const;
  
  std::string protein() const;
    
  void details( Variant & ) const;

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
 public:
    static std::string getrc(const std::string &);
    static std::string getc(const std::string &);
    static std::string getr(const std::string &);

    static std::string translate(std::string &, int, std::vector<std::string> &);    
 private:

    // DNA base --> AA
    static std::map<std::string,std::string> t;


    // Pointers to data-base functions

    static SeqDBase * seqdb;

    static LocDBase * db;

    // Hold region map, populated by fetch_transcripts()

    static std::map<uint64_t,Region> rmap;
    
    static uint64_t transcript_group_id;


    /* static bool load_transcripts( uint64_t id ); */

    /* static bool load_transcripts( uint64_t id , const std::set<Region> & ); */

    // get transcript from LOCDB and/or cache
    bool fetch_transcript( const std::string & name );

    // primary annotation function -- here pregion can either be 
    static std::set<SeqInfo> annotate( int,int,
				       const std::string & alt,
				       const std::string & ref, 
				       const std::vector<uint64_t> & pregion );

    static bool annotate(Variant & var , const std::vector<uint64_t> & );

    // pull a region from the cache
    static Region * from_cache( uint64_t id );
    
    // add regions to the cache, if not already in there
    static void add_transcripts( const std::vector<uint64_t> & id );


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
    
    // Tell Annotate which group to look at for transcripts
    static bool set_transcript_group( const std::string & );
    
    // 3 main entry points to the Annotate::annotate() command; these all call
    // the same underlying code -- the private annotate(...) -- and return true 
    // if any 'interesting' annotation was added
    
    static bool annotate(Variant & var );
    static bool annotate(Variant & var , const Region & );
    static bool annotate(Variant & var , const std::set<Region> & );

    // other helper functions: pull a region from the store:
    static Region * pointer_to_region( const std::string & );
    
    static std::string translate_reference( const Region & region , bool verbose = false);
    
};

#endif

