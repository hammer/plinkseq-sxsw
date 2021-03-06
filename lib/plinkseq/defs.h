#ifndef __DEFS_H__
#define __DEFS_H__

#include <vector>
#include <string>
#include <map>

#include <inttypes.h>

struct PLINKSeq {    

  static void register_standard_metatypes();

  static std::string & VERSION_NUMBER();

  static int         & PROJECT_VERSION_NUMBER();
  static int         & VARDB_VERSION_NUMBER();
  static int         & LOCDB_VERSION_NUMBER();
  static int         & REFDB_VERSION_NUMBER();

  static std::string & SQLITE_SCRATCH_FOLDER();  
  static std::string & PASS_FILTER();  

  static std::string & CURRENT_VCF_VERSION();
  
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
  static std::string & META_GENO_ALT_DOSAGE();

  static std::string & SEQDB_GENOME_BUILD_KEY();
  static std::string & SEQDB_REPEAT_MODE_KEY();
  static std::string & SEQDB_NAME_KEY();
  static std::string & SEQDB_DESCRIPTION_KEY();
  static std::string & SEQDB_IUPAC_KEY();

  static std::string & DELIM();
  static std::string & VCF_MISSING_CHAR();
  static std::string & DEFAULT_LOC_GROUP();
  static std::string & DEFAULT_GENE_SYMBOL();

  static std::string & VCF_END_FIELD();
  static std::string & VCF_GENOTYPE();
  static std::string & VCF_GENOTYPE_ACGT();
  static std::string & VCF_GENOTYPE_AC();
  static std::string & VCF_GENOTYPE_NONREF();
  static std::string & VCF_GENOTYPE_NULL();

  static std::string & ANNOT();
  static std::string & ANNOT_TYPE();
  static std::string & ANNOT_GENE();
  static std::string & ANNOT_ALIAS_GROUPS();
  static std::string & ANNOT_ALIAS_GROUP_WORST();
  static std::string & ANNOT_CODING();
  static std::string & ANNOT_EXONIC();
  static std::string & ANNOT_CHANGE();
  static std::string & ANNOT_CODON();
  static std::string & ANNOT_PROTEIN();
  static std::string & ANNOT_DETAILS();
  static std::string & ANNOT_SUMMARY();

  static std::string & DEFAULT_PHENOTYPE();

  static std::string & META_DP();
  static std::string & META_GENO_DP();
  static std::string & META_GENO_GQ();
 

  static double & DEFAULT_PREV();
  static std::string & DEFAULT_AD();
  static std::string & DEFAULT_TRANS();
  static std::string & DEFAULT_FUNC();
  static double & DEFAULT_AB_HETMIN();
  static double & DEFAULT_AB_HETMAX();
  static double & DEFAULT_AB_HOMMAX();



  static std::string & DEFAULT_G1K();
  static std::string & DEFAULT_DBSNP();

  static long unsigned & DEFAULT_RNG_SEED();

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

enum sType { UNKNOWN_SEX , 
	     MALE , 
	     FEMALE };
	     

enum ploidy_t { PLOIDY_UNKNOWN = 0 , 
		PLOIDY_HAPLOID = 1 , 
		PLOIDY_AUTOSOMAL = 2 , 
		PLOIDY_X = 3 , 
		PLOIDY_Y = 4 , 
		PLOIDY_OTHER = 5 };

enum genotype_model_t { GENOTYPE_MODEL_ALLELIC = 0 , 
			GENOTYPE_MODEL_ALLELIC2 ,
			GENOTYPE_MODEL_ALLELIC3 ,
			GENOTYPE_MODEL_DOM ,
			GENOTYPE_MODEL_REC ,
			GENOTYPE_MODEL_REC2 , 
			GENOTYPE_MODEL_CN , 
			GENOTYPE_MODEL_NULL ,
			GENOTYPE_MODEL_DOSAGE , 
			GENOTYPE_MODEL_PROB_REF , 
			GENOTYPE_MODEL_PROB_HET , 
			GENOTYPE_MODEL_PROB_HOM ,
			GENOTYPE_MODEL_UNSPEC };

enum mType { META_FLAG      = 0 ,   //  <none> 
	     META_UNDEFINED = 1 ,   //  string
	     META_TEXT      = 2 ,   //  string
	     META_INT       = 3 ,   //  int32
	     META_FLOAT     = 4 ,   //  double
	     META_BOOL      = 5 ,   //  bool 
	     META_CHAR      = 6 };  //  string, 1char limit

// Also see plinkseq/meta.h
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

enum merge_mode_t { 
    MERGE_MODE_NONE,       // never merge if different ALT alleles
    MERGE_MODE_EXACT,      // default; merge only if ALT (and REF) have similar sizes
    MERGE_MODE_ANY_OVERLAP // try to merge any overlapping variants
};
  
enum fType { INVALID,
	     UNKNOWN,
	     OUTPUT,     // Folder for output
	     METAMETA ,  // Meta-information on meta-information
	     RESOURCES,  // Folder for core databases, libraries, etc
	     TEMP,       // Folder for SQLite temp files
	     VCF , 
	     BCF_FILE , 
	     BGZF_VCF ,
	     GTF , 
	     PHE , 
	     IND ,
	     LOCDB ,    
	     SEGDB , 
	     INDDB ,
	     VARDB , 
	     NETDB , 
	     WGTDB ,
	     LOG ,     // Output log file
	     FIDX ,    // Main file index (input)
	     SEQDB,    // Human genome reference database
	     REFDB,
	     PWD ,      // not a file, special code for project pasword
             PARAM } ;  // not a file, special code for any parameter setting

	         

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
