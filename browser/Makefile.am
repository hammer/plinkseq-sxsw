bin_PROGRAMS = pbrowse.cgi mongoose
bin_SCRIPTS = pbrowse

mongoose_SOURCES = \
mongoose.h \
mongoose.c \
main.c

pbrowse_cgi_SOURCES = \
psb.cpp \
psb.h \
cgi.h

pbrowse_cgi_LDADD = ../lib/libplinkseq.la -ldl -lpthread -lz

pbrowse_cgi_CPPFLAGS = -I../lib

if WITH_PROTOBUF
pbrowse_cgi_CPPFLAGS += -I $(PROTOBUF_PREFIX)/include
endif


mongoose_CPPFLAGS = 

mongoose_LDFLAGS = -ldl -lpthread

ACLOCAL_AMFLAGS = -I m4
