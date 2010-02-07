/*

  rec_ck.c: inspecting recombination tables 
  Copyright 2009 Joel J. Adamson 

  $Id$

  Joel J. Adamson	-- http://www.unc.edu/~adamsonj
  University of North Carolina at Chapel Hill
  CB #3280, Coker Hall
  Chapel Hill, NC 27599-3280
  <adamsonj@email.unc.edu>

  This file is part of haploid

  haploid is free software: you can redistribute it and/or modify it
  under the terms of the GNU General Public License as published by the
  Free Software Foundation, either version 3 of the License, or (at your
  option) any later version.

  haploid is distributed in the hope that it will be useful, but WITHOUT
  ANY WARRANTY; without even the implied warranty of MERCHANTABILITY or
  FITNESS FOR A PARTICULAR PURPOSE.  See the GNU General Public License
  for more details.

  You should have received a copy of the GNU General Public License
  along with haploid.  If not, see <http://www.gnu.org/licenses/>.

*/

#include "../src/haploid.h"
#include "../src/bithacks.h"
#include <errno.h>
#include <argp.h>

/* some program info to access with argp () */
const char * argp_program_version = "svn";
const char * argp_program_bug_address = "haploid-bugs@unc.edu";

static struct argp_option rec_ck_options[] = {
  { "numloci", 'n', "N", 0, "Number of loci" },
  { "rvector", 'r', "vector", 0, "Linkage map (length n-1)" },
  { "output", 'o', "FILE", 0, "Output file (default stdout)" },
  { "quiet", 'q', 0, 0, "Quietly examine for consistency (default is to print incorrect entries)" },
  { "verbose", 'v', 0, 0, "Verbose output (i.e. print entire table)" },
  { "keep-going", 'k', 0, 0, "Produce all errors (default behavior exits after first error)" },
  { 0 }
};

const char rec_ck_doc[] = "Program to test recombination algorithm";
const char rec_ck_args_doc[] = "rec_ck";

struct rec_ck_args {
  char * n;			/* number of loci */
  char * rvec;			/* linkage map (r-vector) */
  char * outfile;		/* output file name */
};

struct rec_ck_data
{
  unsigned long nloci;		/* number of loci */
  unsigned long geno;		/* number of genotypes */
  double * rvector;		/* linkage map */
  FILE * outfile;		/* output file pointer */
};

struct rec_ck_err_data {
  int err;			/* error code */
  double * err_row;		/* pointer to row (for printing with
				   debug_print_array_double ()) */
  int mom;
  int dad;
};

static int keep_going_p = 0;
static int run_quiet_p = 0;
static int verbose_p = 0;

static error_t
rec_ck_parser (int key, char * arg, struct argp_state * state)
{

  struct rec_ck_args * rec_ck_args = state->input;
  /* if there are no command line arguments, print usage and exit */
  if (state->argc == 1)
    argp_usage (state);

  switch (key)
    {
    case 'n':
      /* set the value of n to the adress of an integer obtained from
	 strtoul () */
      rec_ck_args->n = arg;
      break;
    case 'r':
      /* get the rvector as a string (must be quoted on
	 command-line) */
      rec_ck_args->rvec = arg;
      break;
    case 'o':
      /* store the file name in our struct */
      rec_ck_args->outfile = arg;
    case 'q':
      run_quiet_p = 1;
      break;
    case 'v':
      verbose_p = 1;
      break;
    case 'k':
      keep_going_p = 1;
      break;
    default:
      return ARGP_ERR_UNKNOWN;
    }
  return 0;      
}

static struct argp rec_ck_argp = 
  { rec_ck_options, rec_ck_parser, 0, rec_ck_doc };

int
rec_ck_rt (struct rec_ck_data * data);


int
rec_ck_set_data (struct rec_ck_args * rec_ck_args,
		 struct rec_ck_data * data);


int
rec_ck_print_row (unsigned long geno,
		  int mom, int dad, double * row);


int
rec_ck_find_err (unsigned long int geno,
		 double rec_table[geno][geno][geno],
		 struct rec_ck_err_data err_data[geno * geno]);

inline void
rec_ck_log_err (struct rec_ck_err_data * err_data,
		int i, int j, double * row);

int
main (int argc, char ** argv)
{

  /* initialize argument structure, and initialize defaults (why can't
     I do this in the definition?) */
  struct rec_ck_args rec_ck_args;
 
  /* parse arguments: */
  argp_parse (&rec_ck_argp, argc, argv,
	      0, 0, &rec_ck_args);
  
  /* set up the data, then pass the complete data struct and the
     outfile, quiet/verbose and keep-going args to the driver
     procedure */
  struct rec_ck_data rec_ck_data;
  int rcsd_err = rec_ck_set_data (&rec_ck_args, &rec_ck_data);
  if (rcsd_err)
    return rcsd_err;

  /* print and check the values set by set_rec_table (): */
  int rc_err = rec_ck_rt(&rec_ck_data);
  if (rc_err)
    return rc_err;
  
  fclose (stdout);
  return 0;
}

int
rec_ck_rt (struct rec_ck_data * rec_ck_data)
{

  /* set up recombination table: */
  int geno = rec_ck_data->geno;
  double rec_table[geno][geno][geno];
  int srt_ck =  set_rec_table (rec_ck_data->nloci,
			       geno,
			       rec_table, rec_ck_data->rvector);
  /* an often used message: */
  char * erange_msg = "Value out of range:\n";
  /* if set_rec_table detected problems, return now! */
  if (srt_ck)
    {
      fprintf (stderr, erange_msg, srt_ck);
      return ERANGE;
    }

  /* array of error data structures: */
  struct rec_ck_err_data err_data[geno * geno];

  int j = 0;
  int i;
  /* rec_ck_find_err returns the total number of errors, giving a
     bound for the do-loop */
  int num_err = rec_ck_find_err (geno, rec_table, err_data);
  if (num_err > 0 && verbose_p == 0)
    {
      do
	{
	  if (run_quiet_p == 0)
	    {
	      fprintf (stderr, erange_msg);
	      fflush (stderr);
	      rec_ck_print_row (geno,
				err_data[j].mom,
				err_data[j].dad,
				err_data[j].err_row);
	    }
	  j++;
	}
      while (keep_going_p == 1 && j < num_err);
    }
  else if (verbose_p == 1)
    for (i = 0; i < geno; i++)
      for (j = 0; j < geno; j++)
	rec_ck_print_row (geno, i, j, rec_table[i][j]);
  /* otherwise we're all just fine and don't need to print anything */
  return 0;
}

int
rec_ck_print_row (unsigned long geno, int mom, int dad, double * row)
{
  char * str_array[geno];
  for (int i = 0; i < geno; i++)
    {
      str_array[i] = malloc (sizeof (char) * 9 + 1);
      if (str_array[i] == NULL)
	{
	  fprintf (stderr, "Null pointer");
	  return ENOMEM;
	}
    }
  /* obtain a string of values: */
  debug_print_array_double (geno, row, "%-9.8fL ", str_array);
  /* print the array to the output file: */
  fprintf (stdout, "%s X %s:",
	       debug_printbits (mom),
	       debug_printbits (dad));
      
  for (int i = 0; i < geno; i++)
    /* future versions should allow the user to modify the
       precision */
    fprintf (stdout, "%s  ", str_array[i]);
  
  fprintf (stdout, "\n");
  return 0;
}

int
rec_ck_find_err (unsigned long int geno, double rec_table[geno][geno][geno],
		 struct rec_ck_err_data err_data[geno * geno])
{
  /* count the number of errors */
  int num_errs = 0;
  _Bool err_p = 0;
  int i, j, k;
  double sum, val;
  /*
    going row-by-row, check for
    1. values greater than one
    2. NANs
    3. sums greater than one
  */
  /* #pragma omp parallel for shared (rec_table) private (i,j) */
  /* parallelize this loop, sharing rec_table, making i,j private */
  /* we only need to wait for the join if keep_going_p == 1; that will
     have to use the omp API to specify */
  /* the specific error-checking tasks can be split into sections; if
     keep_going_p == 0 then the first one that finds an error should
     set the first element of the err_data array and then return; this
     could lead to a race condition, of course, so we'd have to block
     the first member of err_data */
  for (i = 0; i < geno; i++)
    {
      for (j = 0; j < geno; j++)
	{
	  sum = 0.0;
	  /* the return value (num_errs) should reflect the number of
	     problematic rows, not the total number of errors */
	  /* in all cases return as fast as possible: either condition
	     #1 or condition #2 supercedes condition #3, so this makes
	     sense */
	  for (k = 0; k < geno; k++)
	    {
	      val = rec_table[i][j][k];
	      /* conditions #1 and #2: */
	      err_p = (isnan (val) || val > 1.0);
	      if (err_p)
		{
		  num_errs++;			
		  rec_ck_log_err (err_data, i, j, rec_table[i][j]);
		  err_data++;
		}
	      else
		sum += val;
	      
	      if (num_errs > 0 && keep_going_p == 0)
		return num_errs;
	    }

	  if (sum > 1.0)
	    {
	      num_errs++;			
	      rec_ck_log_err (err_data, i, j, rec_table[i][j]);
	      err_data++;
	    }
	  if (num_errs > 0 && keep_going_p == 0)
	    return num_errs;
	}
    }
  
  return num_errs;
}

void
rec_ck_log_err (struct rec_ck_err_data * err_data,
		int i, int j, double * row)
{
  err_data->err = ERANGE;		
  err_data->mom = i;			
  err_data->dad = j;			
  err_data->err_row = row;
}
  
int
rec_ck_set_data (struct rec_ck_args * rec_ck_args,
		 struct rec_ck_data * rec_ck_data)
{
  /* --------------------begin ll region-------------------- */

  /* set outfile */
  if (rec_ck_args->outfile == NULL)
    rec_ck_data->outfile = stdout;
  else if ((rec_ck_data->outfile = freopen(rec_ck_args->outfile, "a", stdout))
      == NULL)
    {
      fprintf (stderr, "Unspecified output stream\n");
      return EACCES;
    }

  /* set numloci */
  /* get values from command-line */
  const char errmsg[] = "Representation error\n";
  char * tail;
  rec_ck_data->nloci = strtoul(rec_ck_args->n, &tail, 10);
  if (errno == ERANGE || tail == rec_ck_args->n)
    {
      fprintf (stderr, errmsg);
      return errno;
    }

  /* set rec_ck_data->rvector */
  rec_ck_data->rvector = malloc (sizeof (double) * rec_ck_data->nloci);
  double val;			/* temporary home for converted value */
  int j = 0;
  while (j < rec_ck_data->nloci)
    {
      /* set the entry in rvec to the current token */
      val = strtod(rec_ck_args->rvec, &tail);
      if (errno == ERANGE && tail == rec_ck_args->n)
	{
	  /* there was a representation error, and no conversion was
	     performed by strtod */
	  fprintf (stderr, errmsg);
	  return errno;
	}
      else if (tail == rec_ck_args->rvec && val == 0)
	/* no conversion was performed and the string was empty */
	/* rec_ck_args->rvec should now be an array of as many values
	   as input by the user */
	break;
      else
	/* assign the value and increment j */
	(rec_ck_data->rvector)[j++] = val;
      /* allocate more space */
      /* make sure we have enough room for the next entry: */
      while (isspace (*rec_ck_args->rvec)) rec_ck_args->rvec++;
      if (*rec_ck_args->rvec == 0)
	break;
      
      /* clean up: */
      rec_ck_args->rvec = tail;
    }
  /* --------------------end ll region-------------------- */
  rec_ck_data->geno = 1 << (rec_ck_data->nloci);

  return 0;
}

