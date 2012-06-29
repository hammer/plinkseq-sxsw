from libcpp.vector cimport vector
from libcpp.string cimport string
from cython.operator cimport dereference as deref, preincrement as inc

cdef extern from "pyplinkseqint.h":
  struct Py_individual_map:
    vector[string] ID

  struct Py_genotype:
    vector[int] GT

  struct Py_sample_variant:
    int FSET
    string REF
    string ALT
    float QUAL
    Py_genotype GENO

  struct Py_variant:
    int CHR
    int BP1
    int BP2
    string ID
    int NS
    Py_sample_variant CON

  cdef void Py_init_Pyplinkseq()
  cdef string Py_gstore_version()
  cdef void Py_set_project(string)
  cdef string Py_summary()
  cdef vector[Py_variant] Py_iterate(string, int)
  cdef Py_individual_map Py_ind_list(string)


class IndividualMap:
  def __init__(self, ID):
    self.ID = ID

class Genotype:
  def __init__(self, GT):
    self.GT = GT

class SampleVariant:
  def __init__(self, FSET, REF, ALT, QUAL, GENO):
    self.FSET = FSET
    self.REF = REF
    self.ALT = ALT
    self.QUAL = QUAL
    self.GENO = GENO

class Variant:
  def __init__(self, CHR, BP1, BP2, ID, NS, CON):
    self.CHR = CHR
    self.BP1 = BP1
    self.BP2 = BP2
    self.ID = ID
    self.NS = NS
    self.CON = CON

def init_Pyplinkseq():
  Py_init_Pyplinkseq()

def gstore_version():
  return Py_gstore_version().c_str()

def set_project(char* project):
  Py_set_project(string(project))

def summary():
  return Py_summary().c_str()

def var_fetch(char* mask="", int limit=1000):
  vars_list = []

  cdef Py_variant var
  cdef Py_sample_variant svar
  cdef Py_genotype geno

  cdef vector[int] ind_gts_vec
  cdef vector[int].iterator nested_it

  # Unpack vector of variants
  cdef vector[Py_variant] vars_vec = Py_iterate(string(mask), limit)
  cdef vector[Py_variant].iterator it = vars_vec.begin()
  while it != vars_vec.end():
    var = deref(it)

    # Unpack vector of individual genotypes
    ind_gts = []
    svar = var.CON
    geno = svar.GENO
    ind_gts_vec = geno.GT
    nested_it = ind_gts_vec.begin()
    while nested_it != ind_gts_vec.end():
      ind_gts.append(deref(nested_it))
      inc(nested_it)

    # Build up Variant object
    genotype = Genotype(ind_gts)
    sample_variant = SampleVariant(svar.FSET,
                                   svar.REF.c_str(),
                                   svar.ALT.c_str(),
                                   svar.QUAL, genotype)
    vars_list.append(Variant(var.CHR,
                             var.BP1,
                             var.BP2,
                             var.ID.c_str(),
                             var.NS,
                             sample_variant))

    inc(it)

  return vars_list

def ind_fetch(char* mask=""):
  # Unpack vector of IDs
  cdef Py_individual_map py_ind_map = Py_ind_list(string(mask))
  cdef vector[string] id_vec = py_ind_map.ID
  cdef vector[string].iterator it = id_vec.begin()
  id_list = []
  while it != id_vec.end():
    id_list.append(deref(it).c_str())
    inc(it)
  return IndividualMap(id_list)


# Initialize the module
Py_init_Pyplinkseq()