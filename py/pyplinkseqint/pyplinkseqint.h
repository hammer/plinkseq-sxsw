#include "plinkseq.h"

#include <string>
#include <vector>

// TODO(mrivas): Add support to access statistics. 
// TODO(mrivas): Add support for netdb.
// TODO(hammer): Add flag to ensure initialization, project loaded
// TODO(hammer): Complete Py_variant Done M.A. Rivas
// TODO(mrivas): Complete Py_variantgroup
// TODO(hammer): Handle multiply sampled variants
// TODO(hammer): Allow phenotype list in Py_ind_list - Done M.A. Rivas


struct Py_Feature{
  
  std::string SOURCE_ID;
  std::string FEATURE_ID;
  std::string FEATURE_NAME;
  std::string PROTEIN_ID;
  int PSTART;
  int PSTOP;
  std::string MSTR;
  std::string CHR;
  int GSTART;
  int GSTOP;

  bool operator<( const Py_Feature & rhs ) const
  {
    if ( PSTART < rhs.PSTART ) return true;
    if ( PSTART > rhs.PSTART ) return false;

    if ( PSTOP < rhs.PSTOP ) return true;
    if ( PSTOP > rhs.PSTOP ) return false;
    
    if ( SOURCE_ID < rhs.SOURCE_ID ) return true;
    if ( SOURCE_ID > rhs.SOURCE_ID ) return false;

    return FEATURE_ID < rhs.FEATURE_ID;

  }

};

struct Py_ProtFeatureSet {
  Py_ProtFeatureSet()
  {
    feat.clear();
  };

  void add( const std::string & t , const Py_Feature & f ){
    feat[ t ].insert( f );
  }

  std::set<Py_Feature> get( std::string & t ){
    return feat[t];
  }
  void clear(){
    feat.clear();
  }
  std::set<Py_Feature> * features( std::string & t )
  {
    std::map<std::string,std::set<Py_Feature> >::iterator ii = feat.find( t );
    if ( ii == feat.end() ) return NULL;
    return &(ii->second);
  };
  std::map<std::string,std::set<Py_Feature> > feat;

};

struct Py_Phenotype { 
  std::string LABELS;
  std::vector<std::string> PHENOTYPE;
};

struct Py_individual_map {
  std::vector<std::string> ID;
  Py_Phenotype PHENO;
  
};

struct Py_genotype {
  std::vector<int> GT;
};

struct Py_sample_variant {
  int FSET;
  std::string REF;
  std::string ALT;
  float QUAL;
  Py_genotype GENO;
};



struct Py_variant {
  int CHR;
  int BP1;
  int BP2;
  std::string ID;
  int NS;
  Py_sample_variant CON;
};

struct Py_variantGroup{
 int NV;
 int SIZE;
 std::string COORD;
 int MIDPOS;
 std::string NAME;
 int NIND;
 int SPAN;
std::vector<Py_variant> VARS;


};

struct Py_locGroup{
 std::string NAME;
 std::string COORD;
 std::string ALIAS;

};

void Py_init_Pyplinkseq();
std::string Py_gstore_version();
void Py_set_project(std::string);
std::string Py_summary();
std::vector<Py_variant> Py_iterate(std::string, int);
std::vector<Py_variantGroup> Py_iterateGroup(std::string mask);
std::vector<Py_locGroup> Py_locview(std::string group);
Py_individual_map Py_ind_list(std::string,std::string);
void Py_seqdbattach(std::string);
void Py_refdbattach(std::string);
void Py_locdbattach(std::string);
void Py_protdbattach(std::string);
void Py_annotate_load(std::string);
void Py_locdb_load_gtf(std::string,std::string);
void Py_locdb_collapse_subregions(std::string,std::string);
std::string Py_seqdb_annotate(int,int,std::string,std::string,std::string,std::string);
std::set<Py_Feature> Py_protdb_fetch(std::string, std::string);
//Py_ProtFeatureSet Py_protdb_lookup(int,int,std::string,std::string,std::string);
// Database Summaries 
std::string Py_locdb_summary();
std::string Py_seqdb_summary();
std::string Py_refdb_summary();
std::string Py_protdb_summary();
