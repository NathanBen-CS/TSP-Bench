#include <stdio.h>
#include <stdlib.h>

#include "Instance.h"
#include "Utilities.h"
#include "tour.h"

int *tour;      /* tour[i] = city at position i (see tour.h) */
int *pos;       /* pos[c]  = position of city c              */

/*******************************************************************/

void tourAlloc (void) {

  tour = (int *) calloc (ncities, sizeof(int));     /* one slot per position  */
  pos  = (int *) calloc (ncities+1, sizeof(int));   /* cities are 1..ncities  */
  if (tour == 0 || pos == 0) fatal ("Out of memory allocating the tour");
}

/*******************************************************************/

void tourSetPositions (void) {

  int i;
  for (i=0; i<ncities; ++i)
    pos[tour[i]] = i;                 /* city tour[i] sits at position i */
}

/*******************************************************************/

/* Reverse the path from city "from" forward to city "to".

   Example: tour 1 2 3 4 5 6, reversePath(2,4) gives 1 4 3 2 5 6.

   Reversing the *other* part of the cycle (succ(to) .. pred(from)) gives the
   same cycle, just traversed the other way round, so we reverse whichever
   part is shorter.  That halves the worst case and is a big speed-up. */

void reversePath (int from, int to) {

  int i, j, len, k, a, b;

  i = pos[from];                           /* left end of the path            */
  j = pos[to];                             /* right end of the path           */
  len = j - i;                             /* path may wrap past the array end */
  if (len < 0) len += ncities;
  len++;                                   /* number of cities on the path    */

  if (2*len > ncities) {                   /* more than half the tour?        */
    i = pos[succ(to)];                     /* reverse the complement instead  */
    j = pos[pred(from)];
    len = ncities - len;
  }

  for (k = len/2; k > 0; --k) {            /* swap the two ends, move inwards */
    a = tour[i];
    b = tour[j];
    tour[i] = b;  pos[b] = i;              /* keep pos[] in step with tour[]  */
    tour[j] = a;  pos[a] = j;
    if (++i == ncities) i = 0;             /* left end steps right (wrapping) */
    if (--j < 0) j = ncities - 1;          /* right end steps left (wrapping) */
  }
}

/*******************************************************************/

/* 2-opt exchange: remove edges (x1,x2) and (y1,y2), add (x1,y1) and (x2,y2).

	before:  x1 -> x2 ... y1 -> y2
	after:   x1 -> y1 ... x2 -> y2      (the middle part is reversed)

   The removed edges must point the same way round the cycle (x2 follows x1
   exactly when y2 follows y1); otherwise the result would be two cycles.
   Works whichever way the tour array currently runs. */

void exchange2 (int x1, int x2, int y1, int y2) {

  if (succ(x1) == x2) reversePath (x2, y1);   /* array runs x1 x2 .. y1 y2 */
  else                reversePath (x1, y2);   /* array runs x2 x1 .. y2 y1 */
}

/*******************************************************************/

double tourLength (void) {

  int i;
  double cost = 0;

  for (i=0; i<ncities-1; ++i)
    cost += distance (tour[i], tour[i+1]);      /* each step along the tour */
  cost += distance (tour[ncities-1], tour[0]);  /* plus the trip home       */
  return cost;
}

/*******************************************************************/

void writeTour (char *filename) {   /* DIMACS format: one city per line */

  int i;
  FILE *out = fopen (filename, "w");

  if (out == 0) fatal ("Can't open tour output file");
  for (i=0; i<ncities; ++i) fprintf (out, "%d\n", tour[i]);
  fclose (out);
}
