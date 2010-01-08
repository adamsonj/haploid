/* Haploid Selection models

   Copyright 2009 Joel J. Adamson 
   $Id$

   Joel J. Adamson	-- http://www.unc.edu/~adamsonj
   University of North Carolina at Chapel Hill
   CB #3280, Coker Hall
   Chapel Hill, NC 27599-3280
   <adamsonj@email.unc.edu>
   
   This file contains application-specific code for a two-sex haploid
   genetics model

   This file is part of haploid

   haploid is free software: you can redistribute it and/or modify it
   under the terms of the GNU General Public License as published by
   the Free Software Foundation, either version 3 of the License, or
   (at your option) any later version.

   haploid is distributed in the hope that it will be useful, but
   WITHOUT ANY WARRANTY; without even the implied warranty of
   MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the GNU
   General Public License for more details.

   You should have received a copy of the GNU General Public License
   along with haploid.  If not, see <http://www.gnu.org/licenses/>.
   
*/

#define IOMODE		1
#define ONLYEIGEN       0
#define ONLYSIM         1
#define BOTH            2
#define CALCMODE	-1
#define ECALCMODE	1240
#define ENOPARAM        666
#define LINEWID         72
#define MAXTRIES        10000
#define TOLEPS          1e-10
#define PRINTINIT       0

/* the derivative of f_i with respect to x_i */
#define DFIDXI(  x, w, wbar, eta, r, dDdxi ) \
  (1 - (x)  *   ((w) / (wbar)))		     \
	- (eta) * (r) * (dDdxi)

/* the derivative of f_i with respect to x_j */
#define DFIDXJ( x, wi, wj, wbar, eta, r, dDdxj )   \
     	(1.0/2.0) * (x) * ((wj) * (wi)/(wbar) - (eta) * (r) * (dDdxj)

/* prototypes for io.c */
/* NOT static */
int io_init (int argc, char ** argv);
void Xst_print (Xstate_t * X, int Xst_init_print);
void Xvec_print (double * vector, char * message);
void Xmat_print (double * matrix, char * message);

/* calculation prototypes */
int finalize_Xstate (Xstate_t * X);
int calc_results (Xstate_t * X, int cmode);
void update_X (Xstate_t * X);

/* eigen.c */
int get_eigen (Xstate_t * X);

/* prototypes for simulations */
int sim (Xstate_t * X, int max_tries);
