
CXX = g++
CXXFLAGS = -O0 -pg -shared -fPIC  -fno-common

INC = -I../share/include/ 

R_SHLIB = 

GSTORE = libpseq.so

LIBHDR = gstore.h variant.h genotype.h meta.h regions.h filemap.h	\
defs.h options.h helper.h zfstream.h vcf.h bcf.h spec.h sqlwrap.h	\
locdb.h vardb.h segments.h varfunc.h variant.pb.h refdb.h seqdb.h	\
annot.h inddb.h individual.h indmap.h phmap.h mask.h vgroup.h allele.h	\
permute.h crandom.h iterrep.h em.h bed.h pp.h pp.pb.h eval.h token.h	\
statistics.h matrix.h cdflib.h dcdflib.h ipmpar.h bcf.h knetfile.h	\
khash.h bgzf.h fisher.h glm.h

LIBSRC = gstore.cpp globals.cpp variant.cpp genotype.cpp meta.cpp	\
regions.cpp filemap.cpp vcf.cpp bcf.cpp spec.cpp zfstream.cpp		\
sqlwrap.cpp helper.cpp options.cpp locdb.cpp vardb.cpp segments.cpp	\
varfunc.cpp variant.pb.cpp refdb.cpp seqdb.cpp annot.cpp inddb.cpp	\
mask.cpp reffuncs.cpp vgroup.cpp allele.cpp individual.cpp indmap.cpp	\
phmap.cpp permute.cpp crandom.cpp defs.cpp em.cpp bed.cpp varmeta.cpp	\
pp.cpp pp.pb.cpp eval.cpp token.cpp statistics.cpp dcdflib.cpp		\
bgzf.cpp knetfile.cpp bcf.cpp fisher.cpp glm.cpp matrix.cpp		\
iterate.cpp vcfiterate.cpp

LIBOBJ = $(LIBSRC:.cpp=.o)

$(GSTORE) : $(LIBOBJ)
	ar -crs $@ $(LIBOBJ)

$(LIBOBJ) : $(LIBHDR) 

.cpp.o : 
	$(CXX) $(CXXFLAGS) -c $*.cpp $(INC)

.SUFFIXES : .cpp .o $(SUFFIXES)

clean :
	rm *.o *.so *~
