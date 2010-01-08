/*
  
  io.c: input-output functions
  
  Copyright 2009 Joel J. Adamson 
  $Id: io.c 494 2009-06-18 20:35:41Z joel $

  Joel J. Adamson	-- http://www.unc.edu/~adamsonj
  University of North Carolina at Chapel Hill
  CB #3280, Coker Hall
  Chapel Hill, NC 27599-3280
  <adamsonj@email.unc.edu>
  

  This file is part of haploid

  haploid is free software: you can redistribute it and/or modify it
  under the terms of the GNU General Public License as published by
  the Free Software Foundation, either version 3 of the License,	or
  (at your option) any later version.

  haploid is distributed in the hope that it will be useful, but
  WITHOUT ANY WARRANTY; without even the implied warranty of
  MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the GNU
  General Public License for more details.

  You should have received a copy of the GNU General Public License
  along with haploid.  If not, see <http:	//www.gnu.org/licenses/>.

  Goals for this module:

  1. Define a printer that will output the values of an Xstate
  2. define a reader that can take the contents of a buffer and build
  and Xstate structure
  3. The actual reader should make use of readline or another i/o library

  Design:

  The program should operate in one of two modes: (a) stability mode,
  wherein the user enters a set of equilibria and a fitness scheme	and
  the program calculates the Jacobians and eigenvalues of the
  equilibria (it should also test if they are, in fact, stable);	or
  (b) recursion mode, where the user enters an initial condition and a
  fitness scheme and the program spits out the paths from the set of
  initial conditions.  Output in (a) mode should take the form of the
  matrices printed in some human-readable form, with the eigenvalues
  and eigenvectors printed under them.

  The program should take both interactive input (for testing)	or
  non-interactive (i.e. it should read a file).

  The file for recursion mode (b) would look like this:

  x: 0.25 0.10 0.10 0.55
  w: 1.0  1.0  1.0  2.0

  This would be easy to write a program to produce (e.g. awk	or
  Python); multiple requests would be separated by an empty line or a
  %.  For example

  %
  x: 0.25 0.10 0.10 0.55
  w: 1.0  1.0  1.0  2.0
  
  %
  x: 0.25 0.10 0.10 0.55
  w: 1.0  1.0  1.0  0.2

  GNU Readline can handle the interactive input and regular glibc
  functions can handle the file input

  So I have three tasks to establish:

  1. Decide where to get the input -- see modes.c
  2. get the input,
  3. put it into a buffer

  *Done with input*
  
  4. construct the Xstates
  5. pass control to the calculating modules
  6. get the output
  7. Print it (to a file, or to the terminal?)

  See further notes in func.h under function prototypes

*/
#include "func.h"
#include "haploid.h"

static int acc_in (char *inbuf, Xstate_t * X);

/* do not inline these next two: they are not part of regular program
   operation */
static int help_mesg (void);
static int version_mesg (void);

static inline int strtoarr (char *string, double *array, int len);

/* prototypes for io.c */
int
io_init (int argc, char **argv)
{
  /*

     choose an input mode to operate under then pass control to the input
     functions; decide between interactive input and non-interactive
     input, terminal output (stdout) or file-directed output

     there will thus be four available CL options:

     * input/output mode *
     -i, --interactive      : user enters data line by line (default if
     no -f option)

     NOTE: this option is really a nicety, the
     default behavior without this option should
     be to accept stdin anyway, unless -f is
     specified, so omitting -f, one needn't
     specify -i; nevertheless, its best to be
     explicit

     -f FILE, --file   =        FILE   : get data from file FILE
     -o FILE, --output =        FILE : file for output, if empty

     * calculating mode *
     -e, --eigen            : eigenvalue mode --- assess stability of
     given points

     -s, --simulation       : simulation mode (DEFAULT) 

     both modes take the same kind of data: a fitness scheme (prefixed by
     w:) followed by a series of points in 4-space (preceded by p:      or
     x:), followed by a recombination rate.  We want to analyze multiple
     fitness schemes and recombination rates.  Since this is a two-sex
     model, we need two lines of fitness values, which the user can
     prefix with "wm:" or "wf:" or not.  

     pass control to input functions; pass a pointer to the input file
     handle

     if acc_in returns without error, then we need to direct output to
     the output file; the most efficient way is to pass the output file
     handle, the calculating mode and the data buffer to a function that
     will successively invoke the calculation, gather the output,       and
     invoke the output functions; should this function return, we can
     either exit with an error, or print a newline and return 0

     NOTE: the scheme used here means that main is void and control
     passes out of main completely and does not return up the call stack

     remember to allocate buf!  and free it when done


   */

  /* possible options */


  char *infile = "";
  char *outfile = "";		/* fileNAMES */

  /* flags */
  static int cmode = CALCMODE;
  static int int_flag = IOMODE;
  static int out_flag = 0;	/* output to stdout */

  /* short options */
  char *shortoptions = "if:o:es";

  /* long options */
  static struct option longopts[] = {
    {"interactive", no_argument, &int_flag, 1},
    {"file", required_argument, &int_flag, 0},
    {"output", required_argument, &out_flag, 1},
    {"eigen", no_argument, &cmode, ONLYEIGEN},
    {"simulation", no_argument, &cmode, ONLYSIM},
    /* some options every program should have */
    {"help", no_argument, 0, 'h'},
    {"version", no_argument, 0, 'v'},
    {0, 0, 0, 0}
  };

  /* process command line */
  int c;
  int option_index = 0;
  while (1)
    {

      c = getopt_long (argc, argv, shortoptions, longopts, &option_index);
      if (c == -1)
	break;

      switch (c)
	{
	case 0:
	  if (longopts[option_index].has_arg == required_argument)
	    /* this could contain the filename */
	    {
	      if (int_flag != IOMODE)
		{
		  /* if int_flag has been set, we need the input filename */
		  infile = optarg;
		  break;
		}
	      else if (out_flag != 0)
		{
		  /* if out_flag has been set, we need the output filename */
		  outfile = optarg;
		  break;
		}
	      /* no other options require arguments */
	      else
		break;
	    }
	  /* short options */

	case 'i':
	  {
	    int_flag = 1;
	    break;
	  }
	case 'f':
	  {
	    int_flag = 0;
	    infile = optarg;
	    break;
	  }
	case 'o':
	  {
	    out_flag = 1;
	    outfile = optarg;
	    break;
	  }
	case 'e':
	  {
	    if (cmode == ONLYSIM)
	      {
		cmode = BOTH;
		break;
	      }
	    /* do nothing */
	    else if (cmode == ONLYEIGEN)
	      break;
	    else if (cmode == BOTH)
	      break;
	    else
	      {
		cmode = ONLYEIGEN;
		break;
	      }
	  }
	case 's':
	  {
	    if (cmode == ONLYEIGEN)
	      {
		cmode = BOTH;
		break;
	      }
	    else if (cmode == ONLYSIM)
	      break;
	    else if (cmode == BOTH)
	      break;
	    else
	      {
		cmode = ONLYSIM;
		break;
	      }
	  }
	case 'h':
	  {
	    help_mesg ();
	    break;
	  }
	case 'v':
	  {
	    version_mesg ();
	    break;
	  }
	case '?':
	  /* no known option: */
	  fprintf (stderr, "What are you kidding?\n");
	}
    }

  /* decide how much to allocate for a line */
  /*

     I would really like this to access an environment variable, but
     the associated conversions are not working well

   */
  int linewid = LINEWID;
  /* allocate the input buffer based on above allocation */
  char *inbuf = malloc (linewid * sizeof (char));
  /*

     be sure to check the results of all these calls to malloc

     or use "xmalloc"

   */

  /* allocate an Xstate */
  Xstate_t *P = malloc (sizeof (Xstate_t));
  Xstate_t *pre_calc_P = malloc (sizeof (Xstate_t));

  /* proceed if we have enough memory for both inbuf and P */
  if (inbuf != NULL && P != NULL && pre_calc_P != NULL)
    {
      /* open the input file */
      if (int_flag == 0)
	{
	  FILE *INFILE = freopen (infile, "r", stdin);

	  if (INFILE == NULL)
	    {
	      perror ("Input does not exist");
	      return errno;
	    }
	}

      if (out_flag != 0)
	{
	  FILE *OUTFILE = freopen (outfile, "w", stdout);

	  if (OUTFILE == NULL)
	    {
	      perror ("Failed to redirect stdout");
	      return errno;
	    }

	}

      int in_ck;
      do
	{
	  /* copy P so that we can start over next time we read the file */

	  /*

	     we have to copy P's state in acc_in so that each time we
	     start reading the file we will update the finalized
	     state

	     we don't want to be updating the state at the end of the
	     calculation, but the one at the start of the calculation

	   */

	  /* 

	     at this point P is ready, but has not been modified during calculations

	     save it NOW before something happens to it!

	   */

	  P->ready_flag = -1;
	  in_ck = acc_in (inbuf, P);
	  *pre_calc_P = *P;

	  if (P->ready_flag == 1)
	    calc_results (P, cmode);
	  else
	    continue;
	  /*

	     set P back to its pre-calculated state

	   */
	  *P = *pre_calc_P;
	  P->ready_flag = -1;

	  /* 

	     next time around the loop P will be modified during the
	     call to acc_in at the top of this loop; if it needs
	     updating with new data the data file will provide it; any
	     other data will stay the same as the first block of input

	   */
	}
      while (in_ck != EOF);

      /* if the data file is still open then modify P and start over */
      /* for right now I will skip this feature and do some testing
         calculations */
      /* probably do this with a while loop */
      free (P);
      free (pre_calc_P);
      free (inbuf);
      fprintf (stdout, "You win!\n");
    }
  else
    {
      perror ("Null pointer \n");
      return ENOMEM;
    }
  return 0;
}

static int
help_mesg (void)
{
  printf ("Usage: haploid [options]\n"
	  "Options: \n"
	  "\t -i, --interactive \tAccept user input on STDIN (default)\n"
	  "\t -f FILE, --file=FILE \tAccept input from file FILE\n"
	  "\t -o FILE, --output=FILE \tPrint results in file FILE (default STDOUT)\n"
	  "Calculation modes:\n"
	  "\t -e, --eigen \tCalculate Jacobi Matrix and eigenvalues\n"
	  "\t -s, --simulation \tSimulation mode (default)\n"
	  "Help options:\n"
	  "\t -h, --help \tThis message\n"
	  "\t -v, --version \tVersion message\n");
  exit (0);
}

static int
version_mesg (void)
{
  time_t secs;
  int rn = time (&secs);
  printf ("Version rn_%i\n", rn);
  exit (0);
}

static int
acc_in (char *inbuf, Xstate_t * X)
{
  /*

     populate the input buffer *buf with contents of *infile or *stdin;
     update the buffer with these contents

     input should be read as the data needed to initialize one Xstate

     return 0 for normal exit status, error code for exit status

   */
  size_t len = 0;
  /*

     just process the data right here, don't bothe with trying to
     save it for later

   */

  /*

     algorithm:
     1. get the first (and second?) character of line input (this
     could be the argument to while())
     a. if it is 'p', 'x', or 'w' proceed, storing this value
     b. the second character must be a colon or a space, consume
     it and proceed to the next character
     c. if the third character is a space, proceed
     d. if next character is a number, put it back, then get the
     line */

  int c = 0;
  int i = 0;			/* counter for fitness matrix row */

  while ((c = getchar ()))
    {

      /*

         0kay, the problem with this is that when getc encounters a
         blank line, it consumes the whole line, then getline consumes
         the NEXT line; then the switch () branches according to the
         blank line, which takes us to the finalize_Xstate return to
         io_init.  This means that if there's only a single blank line
         between intended data points we never update the Xstate; the
         new data is gone!

         conditionally getting the line only when there are data points
         doesn't work because then we just go along the line instead of
         skipping to the next line

         What we need is a way to unconditionally consume the line, and
         use the data if we get a blank line

         for the bottom three options here (comments, blank lines and
         EOF), the data from getline is not used.  I need the middle
         ground!  And I want you to give it to me!!!!!

       */

      /* 

         the solution: set the ready flag to zero when it is in
         processing (when we're getting data); before that it will be
         -1

         then if we encounter a blank line *before* adding data, we
         will know to keep going; if we encounter a blank line *after*
         processing, we will finalize and return

       */

      if (c == '\n' && X->ready_flag == 0)
	{
	  /* have we started processing the data? */
	  /* if so, finalize and return */
	  return finalize_Xstate (X);
	}
      else if (c == '\n' && X->ready_flag == -1)
	/* if the ready_flag is -1 but we encounter a newline, we need
	   to go to the top of the loop and get the first character
	   of the NEXT line and process according to that */
	return 0;
      else
	{
	  /* if we don't have a blank line at all, then proceed as below */

	  getline (&inbuf, &len, stdin);
	  switch (c)
	    {

	    case 'w':
	      {

		if (i > 1)
		  break;
		else
		  strtoarr (inbuf, X->ws[i++], 4);
		X->ready_flag = 0;
		break;
	      }
	    case 'p':
	    case 'x':
	      {
		/*

		   both p and x are suitable prefixes for haplotype
		   frequencies (initial conditions)

		 */

		strtoarr (inbuf, X->p, 4);
		X->ready_flag = 0;
		break;
	      }
	    case 'r':
	      {

		strtoarr (inbuf, &X->r, 1);
		X->ready_flag = 0;
		break;
	      }
	    case '#':
	      /* comment, skip */
	      break;
	    case EOF:
	      {
		/*

		   close the input file

		   return -1 (EOF) so the main loop will stop

		 */
		if (X->ready_flag == -1);
		else
		  finalize_Xstate (X);

		return EOF;
	      }
	    default:
	      break;
	    }
	}

      /*

         I can accomplish all of this with a while loop


         2. tokenize the numbers in the line, transfering them to
         double floats
         3. call the appropriate "update_" function depending on the
         character (use a switch statement)
         4. peek at the next line; if it is empty, return; the Xstate
         will be modified and the next time around, we will modify the
         same Xstate address (i.e., it may have the same fitness
         scheme, and a different initial condition, or the same initial
         condition and a different fitness scheme); at that point, we
         will read the stream again

         In the caller: decide whether to re-use the existing Xstate

         alternatively, 2 & 3 could be replaced by a single function
         that does both

       */

    }
  return finalize_Xstate (X);
}

static inline int
strtoarr (char *string, double *array, int len)
{
  int i;

  char *delims = " :\n";

  char *token = strtok (string, delims);
  for (i = 0; i < len; i++)
    {
      array[i] = strtod (token, NULL);
      token = strtok (NULL, delims);
    }

  return 0;
}

void
Xst_print (Xstate_t * X, int Xst_init_print)
{
  if (Xst_init_print == 0)
    {
      /* print a header */
      fprintf (stdout, " %-11s %-11s %-11s %-11s %-11s\n", "x1", "x2", "x3",
	       "x4", "D");
    }
  int i;
  for (i = 0; i < 4; i++)
    fprintf (stdout, "%- 9.8f ", X->p[i]);
  fprintf (stdout, "%- 9.8f\n", X->D);
}

void
Xmat_print (double *matrix, char *message)
{
  /* print the matrix *matrix */

  if (matrix == NULL)
    {
      /* if you gave it a null pointer, it will be unhappy */
      perror ("Null pointer\n");
      exit (ENOMEM);
    }
  else
    {
      fprintf (stdout, "\n%s \n", message);
      int i, j;
      for (j = 0; j < 4; j++)
	{
	  for (i = 0; i < 4; i++)
	    {
	      fprintf (stdout, "%- 9.8g\t", *matrix);
	      matrix++;
	    }
	  fprintf (stdout, "\n");
	}
    }
}

void
Xvec_print (double *vector, char *message)
{
  /* print the vector *vector */

  if (vector == NULL)
    {
      /* if you gave it a null pointer, it will be unhappy */
      perror ("Null pointer\n");
      exit (ENOMEM);
    }
  else
    {
      fprintf (stdout, "\n%s\n", message);
      int i;
      for (i = 0; i < 4; i++)
	{
	  fprintf (stdout, "%- 9.8f \n", *vector);
	  vector++;
	}
    }
}
