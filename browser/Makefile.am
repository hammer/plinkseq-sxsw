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

pbrowse_cgi_LDADD = ../lib/libplinkseq.la

pbrowse_cgi_CPPFLAGS = -I../lib

mongoose_CPPFLAGS = 

mongoose_LDFLAGS = 

ACLOCAL_AMFLAGS = -I m4
