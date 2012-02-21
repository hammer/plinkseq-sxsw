bin_PROGRAMS = index.cgi
bin_SCRIPTS = index

index_cgi_SOURCES = \
index.cpp \
index.h \
autohtml.h \
cgi.h

index_cgi_LDADD = ../lib/libplinkseq.la -ldl -lpthread -lz

index_cgi_CPPFLAGS = -I../lib

if WITH_PROTOBUF
index_cgi_CPPFLAGS += -I $(PROTOBUF_PREFIX)/include
endif


ACLOCAL_AMFLAGS = -I m4
