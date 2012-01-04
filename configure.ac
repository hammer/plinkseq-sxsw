#                                               -*- Autoconf -*-
# Process this file with autoconf to produce a configure script.

# TODO: add --with-gseq and --with-acdb options
# TODO: should we use AC_CHECK_TYPES? (are we worried about int64_t cross platform issues?)
# todo: should the path to protobuf be a "precious variable"?  (page 206)

AC_PREREQ(2.61)
AC_INIT([plinkseq], [0.08], [bthomas@atgu.mgh.harvard.edu])

# foreign allows automake to run without the AUTHORS, COPYING, etc. files
AM_INIT_AUTOMAKE([-Wall -Werror foreign])

LT_PREREQ([2.2])
LT_INIT([pic-only disable-static])

AC_CONFIG_MACRO_DIR([m4])
AC_CONFIG_SRCDIR([lib/pseq.h])
AC_CONFIG_HEADERS([config.h])

# --enable-debug
AC_ARG_ENABLE([debug],
AS_HELP_STRING([--enable-debug], [[Enable debugging build, which builds faster. \
"make install" doesn't work - pseq must be run from within the build directory. ]]), 
[debug=true], 
[debug=false]
)
AM_CONDITIONAL([DEBUG], [test x$debug = xtrue])

# --enable-static
AC_ARG_ENABLE([static],
AS_HELP_STRING([--enable-static], [[Static build, which creates a single \
statically linked pseq binary. Will not work on some platforms, most notably Mac OS X. ]]),
[static=true],
[static=false]
)
AM_CONDITIONAL([STATIC], [test x$static = xtrue])

#
# Protobuf dependency
#

# did user specify --with-protobuf? propagate to automake if she did. 
AC_ARG_WITH([protobuf],
        [AS_HELP_STRING(
                [--with-protobuf=/path/to/protobuf],
                [Specify path to protobuf]
                )],
        [AC_SUBST([PROTOBUF_PREFIX], [${withval}])],
        [protobuf=false]
)

# set WITH_PROTOBUF variable so the Makefile knows whether protobuf was specified on command line
AM_CONDITIONAL([WITH_PROTOBUF], [test "x$protobuf" != xfalse])

# if --with-protobuf wasn't given, try to build with machine installation of protobuf
AS_IF(
[test "x$protobuf" = "xfalse"],
AC_CHECK_LIB(protobuf, main, 
[AC_SUBST([HAVE_LIBPROTOBUF], [1])],
[AC_MSG_ERROR([[

Cannot find protobuf, which is required to build PLINK/SEQ. 
Please install protobuf (2.4.1) on your system. 
If you do not have root access, you can install protobuf to 
a local directory and link to it with:

     ./configure --with-protobuf=/path/to/protobuf/prefix

]])
])
)

AC_CONFIG_FILES([
Makefile 
lib/Makefile 
client/Makefile 
browser/Makefile 
browser/pbrowse
util/Makefile
])

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

 Now type 'make @<:@<target>@:>@'
   where the optional <target> is:
     all                - build all binaries
     install            - install everything

--------------------------------------------------"
