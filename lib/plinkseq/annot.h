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


enum seq_annot_t { UNDEF,     // could not annotate
                   MONO,     // monomorphic site

		   IGR,     // intergenic region

		   // near gene?
		   INTRON,     // intronic

		   // non-protein-coding:
		   UTR5,     // 5' UTR
		   UTR3,     // 3' UTR
		   NPC_RNA, // variants in non-protein-coding exons

		   // protein-coding:
		   SYN,    // synonymous allele
		   INDEL, // any indel

		   // non-synon coding
		   MIS,    // missense allele
 		   PART,    // partial codon  -- not used
 		   CODONINSERTION,	   // codon insertion
 		   CODONDELETION,    // codon deletion
		   STOPINSERTION,    // stop insertion
           STOPDELETION,    // stop deletion
		   OOFCODONINSERTION,    // out of frame codon insertion
		   OOFCODONDELETION, // out of frame codon deletion
 		   SPLICE, // general splice +/- 5bp
 		   EXONIC_UNKNOWN, // overlaps an exon, but since the ALT is 'N', cannot know its exact coding impact
		   
 		   // Special class of splice variants : Faustino and Cooper. Pre-mrna splicing and human disease. AG|G   AG|GTNAG. This is consistent with splicing motif measures.
		   DONORIN2,       // donor splice-site |[GT]
		   DONOREX2AG,       // donor splice-site ex2ag [AG]|
		   DONORIN45AG,      // donor splice-site in45ag |GTN[AG]
		   ACCEPTOREX1G,     // acceptor splice-site ex1g |[G]
		   ACCEPTORIN2,    // 3' splice-site [AG]|
		   SPLICEDEL, // a deletion that includes both sides of splice junction

		   // Additional LoF annotations 
		   SL,	  // Start Loss
 		   NON,    // nonsense allele
 		   FRAMESHIFT,    // frameshift
		   RT // readthrough
};

struct SeqInfo { 
  bool missense() const { return type == MIS ; }
  bool nonsense() const { return type == NON || type == STOPINSERTION; }
  bool startlost() const { return type == SL ; }
  bool readthrough() const { return type == RT || type == STOPDELETION; }
  bool frameshift() const { return type == FRAMESHIFT ; }
  bool codondeletion() const { return type == CODONDELETION || type == OOFCODONDELETION; }
  bool codoninsertion() const { return type == CODONINSERTION || type == OOFCODONINSERTION;}
  bool splice() const { return type == SPLICE;}
  bool csplice() const { return type == DONOREX2AG || type == DONORIN45AG || type == ACCEPTOREX1G ; }
  bool esplice() const { return type == DONORIN2 || type == ACCEPTORIN2 || type == SPLICEDEL; }
  bool exonic_unknown() const { return type == EXONIC_UNKNOWN; }

  bool coding() const { return type != UNDEF && type != MONO && type != IGR && type != INTRON && type != UTR5 && type != UTR3 ; }
  bool nonsyn() const { return type == PART || missense() || nonsense() || startlost() || readthrough() || frameshift() || codondeletion() || codoninsertion() || splice() || csplice() || esplice(); }
  bool synon() const { return type == SYN ; }
  bool indel() const { return type == INDEL || frameshift() || codondeletion() || codoninsertion(); }
  bool intergenic() const { return type == IGR ; }
  bool intronic() const { return type == INTRON; }
  bool utr3() const { return type == UTR3; }
  bool utr5() const { return type == UTR5; }
  bool npcRNA() const { return type == NPC_RNA; }

  bool invalid() const { return type == UNDEF || type == MONO; }

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
	   const std::set<std::string> & aliases,
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
    : transcript(transcript), aliases(aliases), type(type),
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
  std::set<std::string> aliases;

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
  
  std::string transcript_name() const { return transcript == "" ? "." : transcript ; }
};



class Annotate {
    
    // Helper functions
 public:
    static std::string getrc(const std::string &);
    static std::string getc(const std::string &);
    static std::string getr(const std::string &);

    static std::string translate(std::string &, int, std::vector<std::string> &, unsigned int& missingBases);

 private:

    // DNA base --> AA
    static std::map<std::string,std::string> t;


    // Pointers to data-base functions

    static SeqDBase * seqdb;

    static LocDBase * db;

    // Hold region map, populated by fetch_transcripts()

    static std::map<uint64_t,Region> rmap;
    
    static uint64_t transcript_group_id;

    static std::set<uint64_t> alias_group_ids;

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

    typedef std::map<std::string, int> AnnotToCount;
    static std::string getWorstAnnotation(const AnnotToCount& potentialWorstAnnotations);

    // pull a region from the cache
    static Region * from_cache( uint64_t id );
    
    // add regions to the cache, if not already in there
    static void add_transcripts( const std::vector<uint64_t> & id );

    static const std::string DEFAULT_PRIORITIZED_WORST_ANNOTATIONS[];
    static std::list<std::string> PRIORITIZED_WORST_ANNOTATIONS;

 public:
    static void setWorstAnnotationPriorities(std::string prioritiesFile);
    
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
    
    // Tell Annotate which alias group to add in for transcripts
    static bool set_alias_groups( const std::set<std::string>& addAliases );

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

