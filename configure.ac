#                                               -*- Autoconf -*-
# Process this file with autoconf to produce a configure script.

# TODO: add --with-gseq and --with-acdb options
# TODO: should we use AC_CHECK_TYPES? (are we worried about int64_t cross platform issues?)
# todo: should the path to protobuf be a "precious variable"?  (page 206)

AC_PREREQ(2.61)
AC_INIT([plinkseq], [0.07], [bthomas@atgu.mgh.harvard.edu])

# foreign allows automake to run without the AUTHORS, COPYING, etc. files
# silent-rules allows you to add --enable-silent-rules to ./configure for less disastrous output
AM_INIT_AUTOMAKE([-Wall -Werror silent-rules foreign])

LT_PREREQ([2.2])
LT_INIT([])

AC_CONFIG_MACRO_DIR([m4])

# note that util is not a subdir, it is part of the root build
AC_CONFIG_SUBDIRS([lib client])
AC_CONFIG_FILES([Makefile util/Makefile])

# Checks for programs.
AC_PROG_CXX
AC_PROG_CC

AC_OUTPUT

# This is just a fun message to tell the user what's up after configure is compete

echo \
"-------------------------------------------------

 PLINK/SEQ Version ${PACKAGE_VERSION}

 Prefix: '${prefix}'.
 C Compiler: '${CC} ${CFLAGS} ${CPPFLAGS}'
 C++ Compiler: '${CXX} ${CXXFLAGS} ${CPPFLAGS}'

 Package features:
   Gseq: Yes
   ACDB: Yes

 Now type 'make @<:@<target>@:>@'
   where the optional <target> is:
     all                - build all binaries
     install            - install everything

--------------------------------------------------"