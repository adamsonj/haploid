# Copyright 2009 Joel J. Adamson
# $Id: Makefile 497 2009-06-20 00:23:00Z joel $: Makefile 468 2009-06-16 12:28:46Z joel $
# Joel J. Adamson	-- http://www.unc.edu/~adamsonj
# University of North Carolina at Chapel Hill
# CB #3280, Coker Hall
# Chapel Hill, NC 27599-3280
# <adamsonj@email.unc.edu>
#
# This file is part of haploid

# haploid is free software: you can redistribute it and/or modify it
# under the terms of the GNU General Public License as published by
# the Free Software Foundation, either version 3 of the License, or
# (at your option) any later version.

# haploid is distributed in the hope that it will be useful, but
# WITHOUT ANY WARRANTY; without even the implied warranty of
# MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the GNU
# General Public License for more details.

# You should have received a copy of the GNU General Public License
# along with haploid.  If not, see <http://www.gnu.org/licenses/>.

lib_LTLIBRARIES = libhaploid.la
libhaploid_la_SOURCES = src/recomb.c src/spec_func.c \
	src/mating.c src/geno_func.c src/debug.c
include_HEADERS = src/haploid.h src/bithacks.h
noinst_HEADERS = src/haploidpriv.h
noinst_PROGRAMS = rec_ck
rec_ck_SOURCES = tests/rec_ck.c

AM_CFLAGS = -std=gnu99 $(OPENMP_CFLAGS) -O3 -Wall -pedantic
ACLOCAL_AMFLAGS = -I m4 

# Documentation
info_TEXINFOS = doc/haploid.texi
haploid_TEXINFOS = doc/fdl-1.3.texi

# bin program
bin_PROGRAMS = rec_ck
rec_ck_SOURCES = rec_ck.c

# Tests and checks: each check must be a standalone program 
LDADD = -lm libhaploid.la
check_PROGRAMS = sim_stop tlta pop_ck 
sim_stop_SOURCES = tests/sim_stop.c
tlta_SOURCES = tests/tlta.c
pop_ck_SOURCES = tests/pop_ck.c

TESTS = sim_stop tlta pop_ck 
