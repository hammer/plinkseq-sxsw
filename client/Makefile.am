# TODO: enable -static flag (page 137)

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
views.h \
main.cpp \
func.cpp \
summaries.cpp \
views.cpp \
assoc.cpp \
genic.cpp \
compare.cpp \
ibd.cpp \
ibs.cpp \
loaders.cpp \
clusters.cpp \
varfunc.cpp \
util.cpp \
glm-tests.cpp \
stepup.cpp \
kbac.cpp \
netassoc.cpp \
gseq.cpp \
extra.cpp

pseq_CPPFLAGS = -I../lib $(PROTOBUF_CFLAGS)


# note that we don't need to link to any shared libs here
# they are all linked transitively from libplinkseq
pseq_LDADD = ../lib/libplinkseq.la

ACLOCAL_AMFLAGS = -I m4
