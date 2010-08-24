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
libhaploid_la_SOURCES = src/rec.c src/spec_func.c \
	src/mating.c src/geno_func.c src/bits.c src/sparse.c
include_HEADERS = src/haploid.h 
noinst_HEADERS = src/sparse.h

ACLOCAL_AMFLAGS = -I m4 

# Documentation
info_TEXINFOS = doc/haploid.texi
haploid_TEXINFOS = doc/fdl-1.3.texi

# Tests and examples: each is a standalone program
LDADD = -lm libhaploid.la
check_PROGRAMS = sim_stop pop_ck sparse_test rec_test  diseq 
noinst_PROGRAMS = nrm rm_tlta tlta
rec_test_SOURCES = tests/rec_test.c tests/prtable.c
sim_stop_SOURCES = tests/sim_stop.c
pop_ck_SOURCES = tests/pop_ck.c
sparse_test_SOURCES = tests/sparse_test.c
diseq_SOURCES = tests/diseq.c
nrm_SOURCES = examples/nrm.c
nrm_CFLAGS = $(AM_CFLAGS) $(OPENMP_CFLAGS)
rm_tlta_SOURCES = examples/rm_tlta.c
rm_tlta_CFLAGS = $(AM_CFLAGS) $(OPENMP_CFLAGS)
tlta_SOURCES = examples/tlta.c tests/prtable.c
tlta_CFLAGS = $(AM_CFLAGS) $(OPENMP_CFLAGS)

TESTS = sim_stop pop_ck sparse_test rec_test diseq

# distribution:
sig: dist
	gpg -b --use-agent $(PACKAGE)-$(VERSION).tar.gz
	gpg --verify $(PACKAGE)-$(VERSION).tar.gz.sig
disthost: sig
	for file in  $(PACKAGE)-$(VERSION).tar.gz* ; \
		do rsync -avz $$file $(RMTHOST) ; done
