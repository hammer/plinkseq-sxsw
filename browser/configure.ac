#                                               -*- Autoconf -*-
# Process this file with autoconf to produce a configure script.

AC_PREREQ(2.61)
AC_INIT([plinkseq-browser], [0.07], [brettpthomas@gmail.com])
AM_INIT_AUTOMAKE([-Wall foreign])
LT_PREREQ([2.2])
LT_INIT([])

# added by autoscan because there are C++ source files
AC_LANG([C++])

AC_CONFIG_MACRO_DIR([m4])
AC_CONFIG_SRCDIR([psb.h])
AC_CONFIG_HEADER([config.h])

# Checks for programs.
# (compilers)
AC_PROG_CXX
AC_PROG_CC

# --enable-debug
AC_ARG_ENABLE([debug],
AS_HELP_STRING([--enable-debug], [[Enable debugging build, which builds faster.	\
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

# Checks for header files.
# TODO: need to run autoscan and fill this out

# Checks for typedefs, structures, and compiler characteristics.
AC_HEADER_STDBOOL
AC_C_CONST
AC_C_RESTRICT
AC_TYPE_UINT64_T

# Checks for library functions.
AC_HEADER_STDC
AC_CHECK_FUNCS([pow sqrt])

AC_CONFIG_FILES([Makefile pbrowse])
AC_OUTPUT
