#include "plinkseq.h"

#include <string>
#include <vector>

// TODO(hammer): Add flag to ensure initialization, project loaded
// TODO(hammer): Complete Py_variant
// TODO(hammer): Handle multiply sampled variants
// TODO(hammer): Allow phenotype list in Py_ind_list
struct Py_individual_map {
  std::vector<std::string> ID;
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


void Py_init_Pyplinkseq();
std::string Py_gstore_version();
void Py_set_project(std::string);
std::string Py_summary();
std::vector<Py_variant> Py_iterate(std::string, int);
Py_individual_map Py_ind_list(std::string);
