#ifndef __PSEQ_BED_H__
#define __PSEQ_BED_H__

#include <string>
#include <vector>

class ifstream;
class SeqDBase;
class IndDBase;
class VarDBase;

struct BEDLocus {
  int chr;
  std::string name;
  double pos;
  int bp;
  std::string allele1;
  std::string allele2;
};

class BEDReader {
  
 public:
  
  BEDReader( VarDBase * v )
    {
      vardb = v;     
      inddb = NULL;
      seqdb = NULL;
      file_id = 0;
    }
  

  void set_root( std::string n )
  {
    bed_filename = n + ".bed";
    bim_filename = n + ".bim";
    fam_filename = n + ".fam";
  }
  
  void use_seqdb( SeqDBase * s ) 
  {
    seqdb = s;
  }
  
  void set_tag( const std::string & t ) 
  {
    ftag = t;
  }

  void store_phenotypes( IndDBase * i, std::string p ) 
  {
    inddb = i;
    phenotype_name = p;
  }

  bool read_bed();
  

 private:
  
  bool read_header( std::ifstream & B );

  int read_fam();

  int read_bim();

  
  // ID handling: default is to use FID only unless not unique, in which case, use 
  // look to IID, otherwise, use FID_IID

  int id_mode;

  std::string fam_filename;
  std::string bim_filename;
  std::string bed_filename;

  std::string ftag;
  
  std::string phenotype_name;
  
  int file_id;

  // MAP information
  std::vector<BEDLocus> locus;
  
  // DB pointers

  VarDBase * vardb;
  IndDBase * inddb;
  SeqDBase * seqdb;
  
};


#endif
