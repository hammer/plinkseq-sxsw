#                                               -*- Autoconf -*-
# Process this file with autoconf to produce a configure script.

AC_PREREQ(2.61)
AC_INIT([plinkseq-client], [0.07], [brettpthomas@gmail.com])
AM_INIT_AUTOMAKE([-Wall -Werror silent-rules foreign])
LT_PREREQ([2.2])
LT_INIT([])

# added by autoscan because there are C++ source files
AC_LANG([C++])

AC_CONFIG_MACRO_DIR([m4])
AC_CONFIG_SRCDIR([args.h])
AC_CONFIG_HEADER([config.h])

# Checks for programs.
# (compilers)
AC_PROG_CXX
AC_PROG_CC

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

AC_CONFIG_FILES([Makefile])
AC_OUTPUT
