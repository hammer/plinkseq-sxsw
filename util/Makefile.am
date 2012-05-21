bin_PROGRAMS = tab2vcf gcol behead smp

tab2vcf_SOURCES = \
char_tok.h \
tab2vcf.cpp \
char_tok.cpp 

gcol_SOURCES = gcol.cpp

behead_SOURCES = behead.cpp

smp_SOURCES = char_tok.h \
char_tok.cpp \
smp.cpp

smp_LDADD = ../lib/libplinkseq.la -ldl -lpthread -lz

smp_CPPFLAGS = -I../lib
