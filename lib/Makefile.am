# TODO: should this be a shared library or a convenience library?
# I tried to make this a "convenience" library with the noinst_ prefix, but it didn't work
# TODO: set up a pkginclude

lib_LTLIBRARIES = libplinkseq.la

libplinkseq_la_SOURCES = \
allele.h \
annot.h \
bcf.h \
bed.h \
bgzf.h \
cdflib.h \
config.h \
crandom.h \
dcdflib.h \
defs.h \
em.h \
eval.h \
filemap.h \
fisher.h \
genotype.h \
glm.h \
gstore.h \
helper.h \
inddb.h \
individual.h \
indmap.h \
ipmpar.h \
iterrep.h \
khash.h \
knetfile.h \
locdb.h \
mask.h \
matrix.h \
meta.h \
netdb.h \
options.h \
permute.h \
phmap.h \
pp.h \
pp.pb.h \
pseq.h \
refdb.h \
regions.h \
segments.h \
seqdb.h \
sqlite3.h \
sqlwrap.h \
sqlz.h \
statistics.h \
svar.h \
token.h \
vardb.h \
varfunc.h \
variant.h \
variant.pb.h \
vcf.h \
vgroup.h \
zfstream.h \
gstore.cpp \
globals.cpp \
variant.cpp \
svar.cpp \
genotype.cpp \
meta.cpp \
regions.cpp \
filemap.cpp \
vcf.cpp \
bcf.cpp \
zfstream.cpp \
sqlwrap.cpp \
helper.cpp \
options.cpp \
locdb.cpp \
vardb.cpp \
segments.cpp \
varfunc.cpp \
variant.pb.cpp \
refdb.cpp \
seqdb.cpp \
annot.cpp \
inddb.cpp \
mask.cpp \
reffuncs.cpp \
vgroup.cpp \
allele.cpp \
individual.cpp \
indmap.cpp \
phmap.cpp \
permute.cpp \
crandom.cpp \
defs.cpp \
em.cpp \
bed.cpp \
varmeta.cpp \
pp.cpp \
pp.pb.cpp \
eval.cpp \
token.cpp \
statistics.cpp \
dcdflib.cpp \
bgzf.cpp \
fisher.cpp \
glm.cpp \
matrix.cpp \
iterate.cpp \
vcfiterate.cpp \
netdb.cpp \
sqlite3.c

libplinkseq_la_CPPFLAGS = $(PROTOBUF_CFLAGS)
libplinkseq_la_LIBADD = $(PROTOBUF_LIBS)

ACLOCAL_AMFLAGS = -I m4
