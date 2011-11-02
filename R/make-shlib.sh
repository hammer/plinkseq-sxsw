
REXE=~/R

rm Rplinkseq.so

# compile the R interface code
g++ -shared -fPIC -O2 -pg -c rint.cpp -I../lib -I/psych/genetics/pseq/share/include  -I../share/lib64/R/include/

# also, compile a different version of helper, that handles errors differently (i.e. via the R error())
g++ -shared -fPIC -O2 -c ../lib/helper.cpp -DR_SHLIB=1 -I. -I../lib -I/psych/genetics/pseq/share/include  -I../share/lib64/R/include/

$REXE CMD SHLIB \
 -o Rplinkseq.so \
 ./rint.o \
 ./helper.o \
../lib/allele.o \
../lib/crandom.o \
../lib/filemap.o \
../lib/gstore.o \
../lib/iterate.o \
../lib/meta.o \
../lib/pp.pb.o \
../lib/seqdb.o \
../lib/vardb.o \
../lib/vcfiterate.o \
../lib/annot.o \
../lib/dcdflib.o \
../lib/fisher.o \
../lib/knetfile.o \
../lib/options.o \
../lib/refdb.o \
../lib/varfunc.o \
../lib/vcf.o \
../lib/bcf.o \
../lib/defs.o \
../lib/genotype.o \
../lib/inddb.o \
../lib/locdb.o \
../lib/permute.o \
../lib/reffuncs.o \
../lib/sqlwrap.o \
../lib/variant.o \
../lib/svar.o \
../lib/vgroup.o \
../lib/bed.o \
../lib/em.o \
../lib/glm.o \
../lib/individual.o \
../lib/mask.o \
../lib/phmap.o \
../lib/regions.o \
../lib/statistics.o \
../lib/variant.pb.o \
../lib/zfstream.o \
../lib/bgzf.o \
../lib/eval.o \
../lib/globals.o \
../lib/indmap.o \
../lib/matrix.o \
../lib/pp.o \
../lib/segments.o \
../lib/token.o \
../lib/varmeta.o \
../lib/sqlite3.o \
 -L/psych/genetics/pseq/share/lib \
 -L../lib \
 -lprotobuf \
 -lz  \
 -lpthread




