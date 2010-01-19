#                                               -*- Autoconf -*-
# Process this file with autoconf to produce a configure script.
AC_PREREQ([2.63])
AC_INIT([haploid], [svn], [haploid-bugs@unc.edu])
AM_INIT_AUTOMAKE(haploid,svn)
AC_CONFIG_MACRO_DIR([m4])
AC_CONFIG_SRCDIR([src/recomb.c])
AC_CONFIG_HEADERS([config.h])

# Checks for programs.
# 
# don't need C++ right now
# AC_PROG_CXX
AC_PROG_AWK
AC_PROG_CC
AC_PROG_CPP
AC_PROG_INSTALL
AC_PROG_LN_S
AC_PROG_MAKE_SET
AC_PROG_RANLIB
AC_PROG_LIBTOOL

# Checks for libraries.
# FIXME: Replace `main' with a function in `-lm':
AC_CHECK_LIB([m], [pow])

# Checks for header files.
AC_CHECK_HEADERS([fenv.h stdlib.h string.h unistd.h])

# Checks for typedefs, structures, and compiler characteristics.
AC_C_INLINE

# Checks for library functions.
AC_FUNC_MALLOC
AC_CHECK_FUNCS([pow])

# system extensions
AC_USE_SYSTEM_EXTENSIONS(_GNU_SOURCE)

AC_CONFIG_FILES([Makefile])
AC_OUTPUT
