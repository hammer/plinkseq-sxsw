#                                               -*- Autoconf -*-
# Process this file with autoconf to produce a configure script.

AC_PREREQ(2.61)
AC_INIT([plinkseq-R], [0.08], [brettpthomas@gmail.com])

AC_CONFIG_MACRO_DIR([m4])
AC_CONFIG_HEADERS([config.h])


# Checks for programs.
AC_PROG_CXX
AC_PROG_CC

#
# Protobuf dependency
#

# did user specify --with-protobuf? propagate to automake if she did. 
AC_ARG_WITH([protobuf],
        [AS_HELP_STRING(
                [--with-protobuf=/path/to/protobuf],
                [Specify path to protobuf]
                )],
        [protobuf_path=${withval}],
        [protobuf_path=""]
)
AC_ARG_WITH([static-protobuf],
						[AS_HELP_STRING(
							[--with-static-protobuf=/path/to/libprotobuf.a],
							[Specify path to protobuf static library]
						)],
						[protobuf_archive=${withval}],
						[protobuf_archive=""]
						)


if test [ -n "${protobuf_path}" ] ; then
	PROTOBUF_LIBS="-L${protobuf_path}/lib -R${protobuf_path}/lib -Wl,-rpath ${protobuf_path}/lib -lprotobuf"
	PROTOBUF_INCLUDE="-I${protobuf_path}/include"
elif test [ -n "${protobuf_archive}" ] ; then
	PROTOBUF_LIBS="${protobuf_archive}"
else
	PROTOBUF_LIBS="-lprotobuf"
fi

AC_SUBST(PROTOBUF_LIBS)
AC_SUBST(PROTOBUF_INCLUDE)

# Check for other libraries
# TODO: do we really need dl and pthread? 
# I need to explore how these are really used in our code
AC_CHECK_LIB([dl], [main])
AC_CHECK_LIB([pthread], [main])
AC_CHECK_LIB([z], [main])

OTHER_LIBS="-ldl -lpthread -lz"
AC_SUBST(OTHER_LIBS)

# the manual made me do it
: ${R_HOME=`R RHOME`}
if test -z "${R_HOME}"; then
   echo "could not determine R_HOME"
   exit 1
fi
CC=`"${R_HOME}/bin/R" CMD config CC`
CFLAGS=`"${R_HOME}/bin/R" CMD config CFLAGS`
CPPFLAGS=`"${R_HOME}/bin/R" CMD config CPPFLAGS`

AC_CONFIG_FILES([src/Makevars])
AC_OUTPUT
