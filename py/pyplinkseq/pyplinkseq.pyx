from libcpp.vector cimport vector
from libcpp.string cimport string
from libcpp.set cimport set as cpp_set

from cython.operator cimport dereference as deref, preincrement as inc

cdef extern from "pyplinkseqint.h":     
  cdef string Py_locdb_summary()
  cdef string Py_refdb_summary()
  cdef string Py_seqdb_summary()
  cdef string Py_protdb_summary()
  cdef void Py_annotate_load( string )
  cdef void Py_locdb_collapse_subregions( string , string)
  cdef void Py_locdb_load_gtf( string , string )
  cdef void Py_locdbattach( string )
  cdef void Py_seqdbattach( string )
  cdef void Py_refdbattach( string )
  cdef void Py_protdbattach( string )
  struct Py_Phenotype:
    string LABELS
    vector[string] PHENOTYPE

  struct Py_Feature:
    string SOURCE_ID 
    string FEATURE_ID
    string FEATURE_NAME
    string PROTEIN_ID
    int PSTART
    int PSTOP
    string MSTR
    string CHR
    int GSTART
    int GSTOP
    

  struct Py_individual_map:
    vector[string] ID
    Py_Phenotype PHENO

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
   
  struct Py_variantGroup:
    int NV
    int SIZE
    string COORD
    int SPAN
    int MIDPOS
    string NAME
    int NIND
    vector[Py_variant] VARS

  struct Py_locGroup:
    string NAME
    string COORD
    string ALIAS
      
  cdef void Py_init_Pyplinkseq()
  cdef string Py_gstore_version()
  cdef void Py_set_project(string)
  cdef string Py_summary()
  cdef cpp_set[Py_Feature] Py_protdb_fetch(string,string)
  cdef vector[Py_variant] Py_iterate(string, int)
  cdef vector[Py_variantGroup] Py_iterateGroup(string)
  cdef vector[Py_locGroup] Py_locview(string)
  cdef Py_individual_map Py_ind_list(string, string)
  cdef string Py_seqdb_annotate(int, int, string, string, string, string)

class ProtFeature:
  def __init__(self, SOURCE_ID, FEATURE_ID, FEATURE_NAME, PROTEIN_ID, PSTART, PSTOP, MSTR, CHR, GSTART, GSTOP):
    self.SOURCE_ID = SOURCE_ID
    self.FEATURE_ID = FEATURE_ID
    self.FEATURE_NAME = FEATURE_NAME
    self.PROTEIN_ID = PROTEIN_ID
    self.PSTART = PSTART
    self.PSTOP = PSTOP
    self.MSTR = MSTR
    self.CHR = CHR
    self.GSTART = GSTART
    self.GSTOP = GSTOP

class Phenotype:
  def __init__(self, LABELS, PHENOTYPE):
    self.LABELS = LABELS
    self.PHENOTYPE = PHENOTYPE

class IndividualMap:
  def __init__(self, ID , PHENO):
    self.ID = ID
    self.PHENO	= PHENO

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

class VariantGroup:
  def __init__(self, NV, SIZE, COORD, SPAN, MIDPOS, NAME, NIND, CON):
    self.NV = NV
    self.SIZE = SIZE
    self.COORD = COORD
    self.SPAN = SPAN
    self.MIDPOS = MIDPOS
    self.NAME = NAME
    self.NIND = NIND
    self.CON = CON

class GeneData:
  def __init__(self, NAME, COORD, DET, CON):
    self.NAME = NAME
    self.COORD = COORD
    self.DET = DET
    self.CON = CON

class GeneTest:
  def __init__(self, NAMETEST, PVAL, I):
    self.NAMETEST = NAMETEST
    self.PVAL = PVAL
    self.I = I
  
class LocGroup:
  def __init__(self, NAME, COORD, ALIAS):
    self.NAME = NAME
    self.COORD = COORD
    self.ALIAS = ALIAS 
    
def init_Pyplinkseq():
  Py_init_Pyplinkseq()

def gstore_version():
  return Py_gstore_version().c_str()

def set_project(char* project):
  Py_set_project(string(project))

def refdbattach(char* refdb):
  Py_refdbattach(string(refdb))

def seqdbattach(char* seqdb):
  Py_seqdbattach(string(seqdb))

def protdbattach( char* protdb):
  Py_protdbattach(string(protdb))

def locdbattach(char* locdb):
  Py_locdbattach(string(locdb))

def locdb_load_gtf(char* filename, char* name):
  Py_locdb_load_gtf( string(filename) , string(name))

def locdb_collapse_subregions( char* id , char* name):
  Py_locdb_collapse_subregions( string(id), string(name))

def protdb_summary():
  return Py_protdb_summary().c_str()

def locdb_summary():
  return Py_locdb_summary().c_str()

def seqdb_summary():
  return Py_seqdb_summary().c_str()

def refdb_summary():
  return Py_refdb_summary().c_str()

def annotate_load(char* name):
  Py_annotate_load(string(name))

def annotate(int chr, int bp, char* ref, char* alt , char* info , char* transcript):
   refallele = string(ref)
   altallele = string(alt)   
   if len(info) > 0:
      infodat = string(info)

   else:
      infodat = "annotfull"    
   trnsc = string(transcript)
   return Py_seqdb_annotate( chr , bp , refallele , altallele , infodat, trnsc)

def summary():
  return Py_summary().c_str()

def loc_view(char* group):
  loc_list = []
  #Unpack vector of locus
  cdef Py_locGroup loc
  cdef vector[Py_locGroup] loc_vec = Py_locview(string(group))
  cdef vector[Py_locGroup].iterator it = loc_vec.begin()
  while it != loc_vec.end():
    loc = deref(it)
    name = loc.NAME
    coord = loc.COORD
    alias = loc.ALIAS
    loc_list.append(LocGroup(name,coord,alias))
  return loc_list 

def var_groupfetch(char* mask=""):
  group_list = []
  
  cdef Py_variantGroup group
  cdef Py_variant var
  cdef Py_sample_variant svar
  cdef Py_genotype geno
  cdef vector[Py_variant] vars_vec
  cdef vector[Py_variant].iterator vit
  cdef vector[int] ind_gts_vec
  cdef vector[int].iterator nested_it
  #Unpack vector of groups
  cdef vector[Py_variantGroup] group_vec = Py_iterateGroup(string(mask))
  cdef vector[Py_variantGroup].iterator it = group_vec.begin()
  while it != group_vec.end():
    group = deref(it) 
    nvar = group.NV
    size = group.SIZE
    coord = group.COORD
    span = group.SPAN
    midpos = group.MIDPOS
    name = group.NAME
    nind = group.NIND
    vars_list = []
    #Unpack vector of variants 
    vars_vec = group.VARS
    vit = vars_vec.begin()
    while vit != vars_vec.end():
      var = deref(vit)
      #Unpack vector of individual genotypes
      ind_gts = []
      svar = var.CON
      geno = svar.GENO
      ind_gts_vec = geno.GT
      nested_it = ind_gts_vec.begin()
      while nested_it != ind_gts_vec.end():
        ind_gts.append(deref(nested_it))
        inc(nested_it)
      #Build up Variant object
      genotype = Genotype(ind_gts)
      sample_variant = SampleVariant(svar.FSET,
      								 svar.REF.c_str(),
      								 svar.ALT.c_str(),
      								 svar.QUAL,genotype)
      vars_list.append(Variant(var.CHR,
      		                   var.BP1,
      		                   var.BP2,
      		                   var.ID.c_str(),
      		                   var.NS,
      		                   sample_variant))
      inc(vit)
    group_list.append(VariantGroup(nvar,
    							   size,
    							   coord,
    							   span,
    							   midpos,
    							   name,
    							   nind,
    							   vars_list))
    inc(it)
  return group_list 

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

def ind_fetch(char* mask ,char* phenotype):
  indt_list = []
  cdef Py_Phenotype phenodat 
  cdef Py_individual_map py_ind_map = Py_ind_list(string(mask),string(phenotype))  
  cdef vector[string] id_vec = py_ind_map.ID
  cdef vector[string].iterator it = id_vec.begin()
  id_list = []
  while it != id_vec.end():
    id_list.append(deref(it).c_str())
    inc(it)
#  cdef string label_str = "test"
  phenodat = py_ind_map.PHENO
  cdef string label_str = phenodat.LABELS
  cdef vector[string] pheno_vec = phenodat.PHENOTYPE
  cdef vector[string].iterator phen_it = pheno_vec.begin()
  phenos = []
  while phen_it != pheno_vec.end():
    phenos.append(deref(phen_it).c_str())
    inc(phen_it)
#  phenos = ["test1"]
  return IndividualMap(id_list,Phenotype(label_str,phenos))  

def protdbfetch( char* db , char* transcript ):
  cdef cpp_set[Py_Feature] protfeat = Py_protdb_fetch(string(db), string(transcript))
  cdef cpp_set[Py_Feature].iterator prot_it = protfeat.begin()
  protdbs = []
  while prot_it != protfeat.end():
     protint = deref(prot_it)
     source_id = protint.SOURCE_ID
     feature_id = protint.FEATURE_ID
     feature_name = protint.FEATURE_NAME
     protein_id = protint.PROTEIN_ID
     pstart = protint.PSTART
     pstop = protint.PSTOP
     mstr = protint.MSTR
     chr = protint.CHR
     gstart = protint.GSTART
     gstop = protint.GSTOP     
     protdbs.append(ProtFeature(source_id,feature_id,feature_name,protein_id,pstart,pstop,mstr,chr,gstart,gstop))
     inc(prot_it)
  return protdbs

# Initialize the module
Py_init_Pyplinkseq()