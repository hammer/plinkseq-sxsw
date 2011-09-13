#                                               -*- Autoconf -*-
# Process this file with autoconf to produce a configure script.

AC_PREREQ(2.61)
AC_INIT([plinkseq-lib], [0.07], [brettpthomas@gmail.com])

AM_INIT_AUTOMAKE([-Wall -Werror silent-rules foreign])

LT_PREREQ([2.2])
LT_INIT([])

# TODO: should we use AC_SYS_LARGEFILE? it breaks if included with AC_LANG([C++])...
AC_LANG([C++])

AC_CONFIG_MACRO_DIR([m4])
AC_CONFIG_SRCDIR([gstore.h])
AC_CONFIG_HEADERS([config.h])

# Checks for programs.
AC_PROG_CXX
AC_PROG_CC

#
# Set up protobuf dependency
# TODO: should we use pkg-config for this instead?
#
# AC_ARG_VAR([PROTOBUF_PREFIX], [Path where we can find protobuf. If you installed protobuf by source, this was specified by the --prefix option. 
# PROTOBUF_PREFIX/lib should contain libprotobuf.so, and PROTOBUF_PREFIX/include should contain a google/ directory. 
# ])

# if test -n "$PROTOBUF_PREFIX"; then
# 	CPPFLAGS="-I$PROTOBUF_PREFIX/include $CPPFLAGS "
# 	LDFLAGS="-L$PROTOBUF_PREFIX/lib $LDFLAGS"
# 	AC_MSG_WARN([Just set the make vars...........
# 	
# 	])
# fi
# 
# AC_CHECK_LIB([protobuf], [main])
# 
# have_protobuf=no
# AC_CHECK_LIB([protobuf], [main], [have_protobuf=yes])
# 
# if test "x${have_protobuf}" = xno; then 
# 	AC_MSG_ERROR([
# 	-------------------------
# 	Cannot find protobuf. Protobuf version 2.3.0 is required (both header and library files). 
# 	
# 	You can install protobuf via your system's package manager, or build and install it from source. 
# 	If you build from source, configure PLINK/SEQ with the PROTOBUF_PREFIX option set to 
# 	the directory where protobuf was installed, like this: 
# 		
# 		./configure PROTOBUF_PREFIX=/path/to/protobuf/prefix/
# 	
# 	-------------------------
# 	])
# fi
PKG_PROG_PKG_CONFIG
PKG_CHECK_MODULES([protobuf], [protobuf])

# Check for other libraries
AC_CHECK_LIB([dl], [main])
AC_CHECK_LIB([pthread], [main])
AC_CHECK_LIB([z], [main])

# Checks for header files.
AC_HEADER_STDC
AC_CHECK_HEADERS([arpa/inet.h fcntl.h inttypes.h limits.h netdb.h stddef.h stdint.h stdlib.h string.h sys/file.h sys/ioctl.h sys/mount.h sys/param.h sys/socket.h sys/time.h unistd.h utime.h])

# Checks for typedefs, structures, and compiler characteristics.
AC_HEADER_STDBOOL
AC_C_CONST
AC_C_INLINE
AC_TYPE_INT32_T
AC_TYPE_INT64_T
AC_TYPE_INT8_T
AC_TYPE_MODE_T
AC_TYPE_OFF_T
AC_TYPE_SIZE_T
AC_TYPE_SSIZE_T
AC_CHECK_MEMBERS([struct stat.st_blksize])
AC_HEADER_TIME
AC_STRUCT_TM
AC_TYPE_UINT16_T
AC_TYPE_UINT32_T
AC_TYPE_UINT64_T
AC_TYPE_UINT8_T
AC_C_VOLATILE

# Checks for library functions.
AC_FUNC_FSEEKO
AC_FUNC_MALLOC
AC_FUNC_MEMCMP
AC_FUNC_MMAP
AC_FUNC_REALLOC
AC_FUNC_SELECT_ARGTYPES
AC_FUNC_STAT
AC_FUNC_STRERROR_R
AC_FUNC_STRFTIME
AC_FUNC_UTIME_NULL
AC_FUNC_VPRINTF
AC_CHECK_FUNCS([clock_gettime fdatasync floor gethostbyaddr gethostbyname gettimeofday localtime_r memmove memset mkdir munmap pow select socket sqrt strchr strdup strerror strstr strtol utime])

AC_CONFIG_FILES([Makefile])
AC_OUTPUT
