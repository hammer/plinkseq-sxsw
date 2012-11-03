include ../common.am

# TODO: should this be a shared library or a convenience library?
# I tried to make this a "convenience" library with the noinst_ prefix, but it didn't work

lib_LTLIBRARIES = libplinkseq.la

# Install header files to plinkseq directory
nobase_include_HEADERS = \
plinkseq/eval.h \
plinkseq/allele.h \
plinkseq/annot.h \
plinkseq/bcf.h \
plinkseq/vcfz.h \
plinkseq/bed.h \
plinkseq/dose.h \
plinkseq/bgzf.h \
plinkseq/crandom.h \
plinkseq/defs.h \
plinkseq/em.h \
plinkseq/filemap.h \
plinkseq/fisher.h \
plinkseq/genotype.h \
plinkseq/glm.h \
plinkseq/gstore.h \
plinkseq/helper.h \
plinkseq/output.h \
plinkseq/inddb.h \
plinkseq/individual.h \
plinkseq/indmap.h \
plinkseq/ipmpar.h \
plinkseq/iterrep.h \
plinkseq/prob.h \
plinkseq/locdb.h \
plinkseq/mask.h \
plinkseq/matrix.h \
plinkseq/meta.h \
plinkseq/netdb.h \
plinkseq/permute.h \
plinkseq/phmap.h \
plinkseq/protdb.h \
plinkseq/pp.h \
plinkseq.h \
plinkseq/refdb.h \
plinkseq/regions.h \
plinkseq/segments.h \
plinkseq/seqdb.h \
plinkseq/sqlite3.h \
plinkseq/sqlwrap.h \
plinkseq/sqlz.h \
plinkseq/statistics.h \
plinkseq/svar.h \
plinkseq/token.h \
plinkseq/vardb.h \
plinkseq/varfunc.h \
plinkseq/variant.h \
plinkseq/variant.pb.h \
plinkseq/vcf.h \
plinkseq/vgroup.h \
plinkseq/zfstream.h

libplinkseq_la_SOURCES = \
cdflib.h \
dcdflib.h \
khash.h \
pp.pb.h \
eval.cpp \
token.cpp \
vcfz.cpp \
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
output.cpp \
helper.cpp \
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
dose.cpp \
varmeta.cpp \
pp.cpp \
pp.pb.cpp \
prob.cpp \
dcdflib.cpp \
bgzf.cpp \
fisher.cpp \
glm.cpp \
iterate.cpp \
vcfiterate.cpp \
netdb.cpp \
protdb.cpp \
statistics.cpp \
matrix.cpp \
sqlite3.c

# right now, --with-debug and --with-static don't change anything here

libplinkseq_la_CPPFLAGS =
libplinkseq_la_LIBADD = -lz

if WITH_PROTOBUF
libplinkseq_la_CPPFLAGS += -I $(PROTOBUF_PREFIX)/include
libplinkseq_la_LIBADD += -L$(PROTOBUF_PREFIX)/lib
endif

libplinkseq_la_LIBADD += -lprotobuf

ACLOCAL_AMFLAGS = -I m4
