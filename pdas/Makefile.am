bin_PROGRAMS = pdas

pdas_SOURCES = \
das.cpp \
das.h \
socks.cpp \
socks. h

pdas_LDADD = ../lib/libplinkseq.la -ldl -lpthread -lz

pdas_CPPFLAGS = -I../lib

if WITH_PROTOBUF
pdas_CPPFLAGS += -I $(PROTOBUF_PREFIX)/include
endif

ACLOCAL_AMFLAGS = -I m4
