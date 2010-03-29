#                                               -*- Autoconf -*-
# Process this file with autoconf to produce a configure script.
AC_PREREQ([2.63])
AC_INIT([haploid], [svn], [haploid-bugs@unc.edu])
AM_INIT_AUTOMAKE(haploid,svn)
AC_CONFIG_MACRO_DIR([m4])
AC_CONFIG_SRCDIR([src/rec.c])
AC_CONFIG_HEADERS([config.h])
AC_PREFIX_DEFAULT(/usr)

# system extensions
AC_USE_SYSTEM_EXTENSIONS(_GNU_SOURCE)

# Checks for programs.
# 
AC_PROG_AWK
AC_PROG_CC
AC_PROG_CPP
AC_PROG_INSTALL
AC_PROG_LN_S
AC_PROG_MAKE_SET
LT_INIT

# Checks for libraries.
# FIXME: Replace `main' with a function in `-lm':
AC_CHECK_LIB([m], [pow])

# Checks for header files.
AC_CHECK_HEADERS([stdlib.h math.h limits.h string.h limits.h error.h])

# OpenMP:
AC_OPENMP

# Checks for typedefs, structures, and compiler characteristics.
AC_C_INLINE

# Checks for library functions.
AC_FUNC_MALLOC
AC_CHECK_FUNCS([pow] [error])

AC_CONFIG_FILES([Makefile])
AC_OUTPUT
