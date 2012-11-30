include ../common.am

bin_PROGRAMS = pdas

pdas_SOURCES = \
das.cpp \
das.h \
socks.cpp \
socks.h

pdas_LDADD = ../lib/libplinkseq.la -lpthread -lz -ldl 

pdas_CPPFLAGS = -I../lib

if WITH_PROTOBUF
pdas_CPPFLAGS += -I $(PROTOBUF_PREFIX)/include
endif

ACLOCAL_AMFLAGS = -I m4
