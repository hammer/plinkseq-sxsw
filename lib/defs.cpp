
#include "defs.h"
#include "meta.h"

std::string & PLINKSeq::VERSION_NUMBER() { static std::string s = "0.09(2-Apr-2012)"; return s; }

int & PLINKSeq::PROJECT_VERSION_NUMBER() { static int i = 2; return i; }
int & PLINKSeq::VARDB_VERSION_NUMBER() { static int i = 5; return i; }

// v2 added exon/CDS/start/stop features encoded in strand column
int & PLINKSeq::LOCDB_VERSION_NUMBER() { static int i = 2; return i; } 
int & PLINKSeq::REFDB_VERSION_NUMBER() { static int i = 2; return i; } 

std::string & PLINKSeq::SQLITE_SCRATCH_FOLDER() { static std::string s = ""; return s; }

std::string & PLINKSeq::CURRENT_VCF_VERSION() { static std::string s = "VCFv4.1"; return s; } 

// internals
std::string & PLINKSeq::TRANSCRIPT_FRAME() { static std::string s = "_FRAME"; return s; }
std::string & PLINKSeq::TRANSCRIPT_STRAND() { static std::string s = "_STRAND"; return s; }
std::string & PLINKSeq::META_LSET() { static std::string s = "_LSET"; return s; }
std::string & PLINKSeq::META_LGRP() { static std::string s = "_LGRP"; return s; }
std::string & PLINKSeq::META_VSET() { static std::string s = "_VSET"; return s; }
std::string & PLINKSeq::META_VGRP() { static std::string s = "_VGRP"; return s; }
std::string & PLINKSeq::META_LSSET() { static std::string s = "_LSSET"; return s; }
std::string & PLINKSeq::META_LSGRP() { static std::string s = "_LSGRP"; return s; }
std::string & PLINKSeq::META_GROUP() { static std::string s = "_GROUP"; return s; }

// user

std::string & PLINKSeq::PASS_FILTER() { static std::string s = "PASS"; return s; }
std::string & PLINKSeq::DELIM() { static std::string s = "\t"; return s; }
std::string & PLINKSeq::VCF_MISSING_CHAR() { static std::string s = "."; return s; }

std::string & PLINKSeq::VCF_END_FIELD() { static std::string s = "END"; return s; }
std::string & PLINKSeq::VCF_GENOTYPE() { static std::string s = "GT"; return s; }
std::string & PLINKSeq::VCF_GENOTYPE_ACGT() { static std::string s = "GT_S"; return s; }
std::string & PLINKSeq::VCF_GENOTYPE_AC() { static std::string s = "GT_A"; return s; }
std::string & PLINKSeq::VCF_GENOTYPE_NONREF() { static std::string s = "GT_NR"; return s; }
std::string & PLINKSeq::VCF_GENOTYPE_NULL() { static std::string s = "GT_NULL"; return s; }

std::string & PLINKSeq::DEFAULT_LOC_GROUP() { static std::string s = "refseq"; return s; }
std::string & PLINKSeq::DEFAULT_GENE_SYMBOL() { static std::string s = "symbol"; return s; }
std::string & PLINKSeq::DEFAULT_PHENOTYPE() { static std::string s = ""; return s; } 

std::string & PLINKSeq::META_GENE() { static std::string s = "GENE"; return s; }
std::string & PLINKSeq::META_ANNOT() { static std::string s = "TYPE"; return s; }
std::string & PLINKSeq::META_ANNOT_FLAG() { static std::string s = "Nonsense,Missense"; return s; }

std::string & PLINKSeq::META_GENO_LIK() { static std::string s = "GL"; return s; }
std::string & PLINKSeq::META_GENO_PHRED() { static std::string s = "PL"; return s; }
std::string & PLINKSeq::META_GENO_POSTPROB() { static std::string s = "PP"; return s; }
std::string & PLINKSeq::META_GENO_ALT_DOSAGE() { static std::string s = "EC"; return s; }

std::string & PLINKSeq::META_GENO_DP() { static std::string s = "DP"; return s; }
std::string & PLINKSeq::META_GENO_GQ() { static std::string s = "GQ"; return s; }
std::string & PLINKSeq::META_DP() { static std::string s = "DP"; return s; }

std::string & PLINKSeq::DEFAULT_DBSNP() { static std::string s = "dbsnp"; return s; }
std::string & PLINKSeq::DEFAULT_G1K() { static std::string s = "g1k"; return s; }

long unsigned & PLINKSeq::DEFAULT_RNG_SEED() { static long unsigned i = time(0); return i; } 

// misc sys.
std::string & PLINKSeq::SEQDB_GENOME_BUILD_KEY() { static std::string s = "BUILD"; return s; }
std::string & PLINKSeq::SEQDB_REPEAT_MODE_KEY() { static std::string s = "REPEATMODE"; return s; }
std::string & PLINKSeq::SEQDB_NAME_KEY() { static std::string s = "NAME"; return s; }
std::string & PLINKSeq::SEQDB_DESCRIPTION_KEY() { static std::string s = "DESC"; return s; }
std::string & PLINKSeq::SEQDB_IUPAC_KEY() { static std::string s = "IUPAC"; return s; }


std::string & PLINKSeq::ANNOT()         { static std::string s = "_ANNOT"; return s; }
std::string & PLINKSeq::ANNOT_TYPE()    { static std::string s = "_ANNOT_TYPE"; return s; }
std::string & PLINKSeq::ANNOT_GENE()    { static std::string s = "_ANNOT_GENE"; return s; }
std::string & PLINKSeq::ANNOT_CODING()  { static std::string s = "_ANNOT_CODING"; return s; }
std::string & PLINKSeq::ANNOT_EXONIC()  { static std::string s = "_ANNOT_EXONIC"; return s; }
std::string & PLINKSeq::ANNOT_CHANGE()  { static std::string s = "_ANNOT_CHANGE"; return s; }
std::string & PLINKSeq::ANNOT_CODON()   { static std::string s = "_ANNOT_CODON"; return s; }
std::string & PLINKSeq::ANNOT_PROTEIN() { static std::string s = "_ANNOT_PROTEIN"; return s; }
std::string & PLINKSeq::ANNOT_SUMMARY() { static std::string s = "_ANNOT_SUMMARY"; return s; }


void PLINKSeq::register_standard_metatypes()
{

  // ensure the following are always defined, no matter whether they
  // appear in the VCF or not

  //
  // INFO fields
  // 

  MetaInformation<VarMeta>::field( "AA" , META_TEXT   , 1  , "Ancestral allele" );
  MetaInformation<VarMeta>::field( "AC" , META_INT    , -1 , "Allele count for each alternate allele" );
  MetaInformation<VarMeta>::field( "AF" , META_FLOAT  , -1 , "Frequency for each alternate allele (estimated from primary data, not called genotypes)" );
  MetaInformation<VarMeta>::field( "AN" , META_INT    , 1 , "Total number of alleles in called genotypes" );
  MetaInformation<VarMeta>::field( "BQ" , META_FLOAT  , 1 , "RMS base quality at this position" );
  MetaInformation<VarMeta>::field( "CIGAR" , META_TEXT  , 1 , "CIGAR string describing how to align an alternate allele to the reference allele" );
  MetaInformation<VarMeta>::field( "DB" , META_INT    , 1 , "dbSNP membership" );
  MetaInformation<VarMeta>::field( "DP" , META_INT    , 1 , "Combined depth across samples" );
  MetaInformation<VarMeta>::field( "END" , META_INT    , 1 , "End position of the variant" );
  MetaInformation<VarMeta>::field( "H2" , META_FLAG    , 0 , "Membership in HapMap2" );
  MetaInformation<VarMeta>::field( "H3" , META_FLAG    , 0 , "Membership in HapMap3" );
  MetaInformation<VarMeta>::field( "MQ" , META_FLOAT   , 1 , "RMS mapping quality" );
  MetaInformation<VarMeta>::field( "MQ0" , META_INT    , 1 , "Number of MQ == 0 reads covering this record" );
  MetaInformation<VarMeta>::field( "NS" , META_INT    , 1 , "Number of samples with data" );
  MetaInformation<VarMeta>::field( "SB" , META_FLOAT  , 1 , "Strand bias" );
  MetaInformation<VarMeta>::field( "SOMATIC" , META_FLAG  , 0 , "Somatic mutation" );
  MetaInformation<VarMeta>::field( "VALIDATED" , META_FLAG  , 0 , "Validated by follow-up experiment" );

  //
  // FORMAT 
  //

  MetaInformation<GenMeta>::field( "DP" , META_INT   , 1 , "Read depth" );
  MetaInformation<GenMeta>::field( "FT" , META_TEXT  , -1 , "Sample genotype filter(s)" );
  MetaInformation<GenMeta>::field( "GL" , META_FLOAT  , -1 , "log10-scaled likelihoods; for AA,AB,BB where A=ref, B=alt" );
  MetaInformation<GenMeta>::field( "PL" , META_INT  , -1 , "phred-scaled genotype likelihoods; for AA,AB,BB where A=ref, B=alt" );
  MetaInformation<GenMeta>::field( "PP" , META_FLOAT  , -1 , "P(genotype | data)" );
  MetaInformation<GenMeta>::field( "EC" , META_FLOAT  , -1 , "Expected count (dosage) of alternate alleles");
  MetaInformation<GenMeta>::field( "GQ" , META_INT  , 1 , "phred-scaled genotype quality, -10log_10p(genotype call is wrong)" );
  MetaInformation<GenMeta>::field( "HQ" , META_INT  , 2 , "haplotype qualities, two phred qualities" );
  

  //
  // FILTER
  //

  // TODO: I needed to change this so we can compile statically...
  MetaInformation<VarFilterMeta>::field( "PASS" , META_FLAG  , 1 , "Passed variant FILTERs" );
  //MetaInformation<VarFilterMeta>::field( PLINKSeq::PASS_FILTER() , META_FLAG  , 1 , "Passed variant FILTERs" );

  registerMetatype( PLINKSeq::TRANSCRIPT_FRAME() , 
		    META_INT, 1 , META_GROUP_LOC , "CDS Frame" );
  
  registerMetatype( PLINKSeq::TRANSCRIPT_STRAND() , 
		    META_INT, 1 , META_GROUP_LOC , "CDS Strand" );


  //
  // Unless these are explicitly encountered, make invisible
  //

  MetaMeta::is_internal("AA");
  MetaMeta::is_internal("AC");
  MetaMeta::is_internal("AF");
  MetaMeta::is_internal("AN");
  MetaMeta::is_internal("BQ");
  MetaMeta::is_internal("CIGAR");
  MetaMeta::is_internal("DB");
  MetaMeta::is_internal("DP");
  MetaMeta::is_internal("END");
  MetaMeta::is_internal("H2");
  MetaMeta::is_internal("H3");
  MetaMeta::is_internal("MQ");
  MetaMeta::is_internal("MQ0");
  MetaMeta::is_internal("NS");
  MetaMeta::is_internal("SB");
  MetaMeta::is_internal("SOMATIC");
  MetaMeta::is_internal("VALIDATED");
  
  MetaMeta::is_internal("HQ");
  MetaMeta::is_internal("DP");
  MetaMeta::is_internal("GL");
  MetaMeta::is_internal("PL");
  MetaMeta::is_internal("PP");
  MetaMeta::is_internal("EC");
  MetaMeta::is_internal("GQ");
  MetaMeta::is_internal("FT");

  MetaMeta::is_internal( PLINKSeq::TRANSCRIPT_STRAND() );
  MetaMeta::is_internal( PLINKSeq::TRANSCRIPT_FRAME() );

  MetaMeta::is_internal("PASS");
  
}
