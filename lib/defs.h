#ifndef __DEFS_H__
#define __DEFS_H__

#include <vector>
#include <string>
#include <map>

#include <inttypes.h>

struct PLINKSeq {    

  static void register_standard_metatypes();

  static std::string & VERSION_NUMBER();

  static int         & VARDB_VERSION_NUMBER();

  static std::string & SQLITE_SCRATCH_FOLDER();  
  static std::string & PASS_FILTER();  

  static std::string & TRANSCRIPT_FRAME();
  static std::string & TRANSCRIPT_STRAND();

  static std::string & META_LSET();
  static std::string & META_LGRP();
  static std::string & META_VSET();
  static std::string & META_VGRP();
  static std::string & META_LSSET();
  static std::string & META_LSGRP();
  static std::string & META_GROUP();

  static std::string & META_GENE();
  static std::string & META_ANNOT();
  static std::string & META_ANNOT_FLAG();

  static std::string & META_GENO_LIK();
  static std::string & META_GENO_PHRED();
  static std::string & META_GENO_POSTPROB();

  static std::string & SEQDB_GENOME_BUILD_KEY();
  static std::string & SEQDB_REPEAT_MODE_KEY();
  static std::string & SEQDB_NAME_KEY();
  static std::string & SEQDB_DESCRIPTION_KEY();
  static std::string & SEQDB_IUPAC_KEY();

  static std::string & DELIM();
  static std::string & VCF_MISSING_CHAR();
  static std::string & DEFAULT_LOC_GROUP();

  static std::string & ANNOT();
  static std::string & ANNOT_TYPE();
  static std::string & ANNOT_GENE();
  static std::string & ANNOT_CODING();
  static std::string & ANNOT_EXONIC();
  static std::string & ANNOT_CODON();
  static std::string & ANNOT_PROTEIN();

  static std::string & DEFAULT_PHENOTYPE();

  static std::string & META_DP();
  static std::string & META_GENO_DP();
  static std::string & META_GENO_GQ();

  static std::string & DEFAULT_G1K();
  static std::string & DEFAULT_DBSNP();

};



///////////////////////////////////////
// Type definitions

class Genotype;

typedef  std::vector<double>  vec_d;
typedef  std::vector<float>   vec_f;
typedef  std::vector<bool>    vec_b;
typedef  std::vector<int>     vec_i;
typedef  std::vector<std::string>  strList;

typedef  uint64_t         ID_t;

enum sType { MALE , 
	     FEMALE , 
	     UNKNOWN_SEX }; 

enum mType { META_FLAG      = 0 ,   //  <none> 
	     META_UNDEFINED = 1 ,   //  string
	     META_TEXT      = 2 ,   //  string
	     META_INT       = 3 ,   //  int32
	     META_FLOAT     = 4 ,   //  double
	     META_BOOL      = 5 ,   //  bool 
	     META_CHAR      = 6 };  //  string, 1char limit

// Also see meta.h
enum mGroup { META_GROUP_MISC = 0 , 
	      META_GROUP_VAR = 1 , 
	      META_GROUP_GEN = 2 , 
	      META_GROUP_LOC = 3 ,  // == SEG 
	      META_GROUP_REF = 4 , 
	      META_GROUP_FILE = 5 , 
	      META_GROUP_INDIV = 6 ,
	      META_GROUP_ALLELE = 7 ,
	      META_GROUP_FILTER = 8 };  


enum downcode_mode_t { 
  DOWNCODE_MODE_NONE ,
  DOWNCODE_MODE_ALL_ALT ,  
  DOWNCODE_MODE_EACH_ALT
}; 

  
enum fType { INVALID,
	     UNKNOWN,
	     OUTPUT,     // Folder for output
	     METAMETA ,  // Meta-information on meta-information
	     RESOURCES,  // Folder for core databases, libraries, etc
	     TEMP,       // Folder for SQLite temp files
	     VCF , 
	     BCF_FILE , 
	     GTF , 
	     PHE , 
	     IND ,
	     LOCDB ,    
	     SEGDB , 
	     INDDB ,
	     VARDB , 
	     LOG ,     // Output log file
	     FIDX ,    // Main file index (input)
	     SEQDB,    // Human genome reference database
	     REFDB } ;

	         

//
// Disease state variable
//

enum affType { UNKNOWN_PHE = 0 , 
	       CONTROL     = 1 , 
	       CASE        = 2 };

enum pType { PHE_NONE    = 0 ,
	     PHE_DICHOT  = 1 , 
	     PHE_QT      = 2 ,
	     PHE_FACTOR  = 3 };
	     


//
//  FileModes
//

enum fMode { READ = 1 , 
	     WRITE = 2 }; 


//
// Variant type (not currently used / likely to change )
//

enum vType { STANDARD = 1,      // S  STANDARD SNP					    
	     XCHR = 2,          // X  STANDARD X CHR (i.e. either diploid or haploid)    
	     HAPLOID = 3,       // H  HAPLOID					    
	     EXT_HARD = 4,      // G  GENERAL FORMAT (HARD CALL)			    
	     EXT_SOFT = 5 } ;   // Q  GENERAL FORMAT (SOFT CALLS)                        



#endif
