#include "pyplinkseqint.h"

#include <map>

GStore *g;

void Py_init_Pyplinkseq() {
  g = new GStore;
  return;
}

std::string Py_gstore_version() {
  std::map<std::string, std::string> g_version = g->version();
  std::map<std::string, std::string>::iterator iter;
  std::string strToReturn;

  for (iter = g_version.begin(); iter != g_version.end(); ++iter) {
    strToReturn.append(iter->first);
    strToReturn.append("=");
    strToReturn.append(iter->second);
    strToReturn.append("\n");
  }

  return strToReturn.c_str();
}

void Py_set_project(std::string project) {
  std::string project_str(project);
  g->set_project(project_str);
}

std::string Py_summary() {
  return g->summary(false);
}

void accumulate_func(Variant &v, void *p) {
  // Variant
  Py_variant var;
  var.CHR = v.chromosome();
  var.BP1 = v.position();
  var.BP2 = v.stop();
  var.ID = v.name();
  var.NS = v.n_samples();

  // Sample variant
  SampleVariant &sv = v.sample(-1);
  Py_sample_variant svar;
  svar.FSET = sv.fileset();
  svar.REF = sv.reference();
  svar.ALT = sv.alternate();
  svar.QUAL = sv.quality();

  // Genotype
  Py_genotype py_geno;
  std::vector<int> gt;
  int nind = v.size(-1);
  for (int i = 0; i < nind; i++) {
    Genotype* geno = v.genotype(-1, i);
    if (geno->more() || geno->null()) {
      gt.push_back(-1); // Code NA as -1
    } else {
      gt.push_back(geno->allele_count());
    }
  }
  py_geno.GT = gt;
  svar.GENO = py_geno;
  var.CON = svar;

  // Add new variant to the vector of variants
  std::vector<Py_variant>* d = (std::vector<Py_variant>*)p;
  d->push_back(var);
}

std::vector<Py_variant> Py_iterate(std::string mask, int limit) {
  Mask m(mask, "", true, false);
  g->register_mask(m);
  m.limit(limit);
  std::vector<Py_variant> vars;
  g->vardb.iterate(accumulate_func, &vars, m);
  return vars;
}

Py_individual_map Py_ind_list(std::string mask) {
  Mask m(mask, "", true, false);
  g->register_mask(m);
  g->indmap.populate(g->vardb, g->phmap, m);
  const int n = g->indmap.size();
  std::vector<std::string> ids;
  for (int i = 0; i < n; i++) {
    ids.push_back(g->indmap.ind(i)->id());
  }
  Py_individual_map ind_map;
  ind_map.ID = ids;
  return ind_map;
}
