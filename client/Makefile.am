bin_PROGRAMS = pseq

pseq_SOURCES = \
args.h \
assoc.h \
clusters.h \
compare.h \
config.h \
extra.h \
func.h \
genic.h \
ibd.h \
ibs.h \
loaders.h \
netassoc.h \
summaries.h \
util.h \
views.h \
util.cpp \
ibs.cpp \
main.cpp \
func.cpp \
summaries.cpp \
views.cpp \
assoc.cpp \
genic.cpp \
compare.cpp \
ibd.cpp \
loaders.cpp \
clusters.cpp \
varfunc.cpp \
glm-tests.cpp \
stepup.cpp \
kbac.cpp \
netassoc.cpp \
gseq.cpp \
extra.cpp

# note that we don't need to link to any shared libs here
# they are all linked transitively from libplinkseq
# we still need to include headers until we create public headers for the lib
pseq_LDADD = ../lib/libplinkseq.la

pseq_CPPFLAGS = -I../lib
pseq_LDFLAGS = -lz -ldl

if WITH_PROTOBUF
pseq_CPPFLAGS += -I $(PROTOBUF_PREFIX)/include
endif

if DEBUG
pseq_LDFLAGS += -no-install all-static
endif

ACLOCAL_AMFLAGS = -I m4
