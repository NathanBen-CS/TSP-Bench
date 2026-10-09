#include <stdio.h>
#include <stdlib.h>
#include <time.h>

#include "Instance.h"
#include "Utilities.h"

/*********************************************************************

Starter template for your own TSP algorithm.

Usage (the benchmark runner calls every solver this way):

	nn instance_file tour_output_file [seed]

nn is deterministic, so the seed is accepted and ignored.  A randomized
solver should pass it to sprand() from Utilities.c.

Lines printed as key=value on stdout become extra columns in the
benchmark results (e.g. read_time, solve_time, start_length).

What you get from Instance.c after readInstance():

	ncities          number of cities, numbered 1..ncities (NOT 0-based)
	distance(i, j)   distance between city i and city j (works for every
	                 supported instance type: EUC_2D, CEIL_2D, matrix)
	xx[i], yy[i]     coordinates of city i (only if isEuclidean == 1)

What you must produce:

	tour[0..ncities-1], a permutation of 1..ncities.  The tour returns
	from the last city to the first automatically.

The algorithm below is nearest neighbour: start at city 1, always walk
to the closest city not yet visited.  Replace solve() with your own.

 *********************************************************************/

static int *tour;

/*******************************************************************/

static void solve (void) {

  int *visited;
  int i, j, curr, best, bestdist, d;

  visited = (int *) calloc (ncities+1, sizeof(int));

  curr = 1;
  tour[0] = curr;
  visited[curr] = 1;

  for (i=1; i<ncities; ++i) {
    best = 0;
    bestdist = 0x7fffffff;
    for (j=1; j<=ncities; ++j) {
      if (!visited[j] && (d = distance (curr, j)) < bestdist) {
	bestdist = d;
	best = j;
      }
    }
    tour[i] = best;
    visited[best] = 1;
    curr = best;
  }

  free (visited);
}

/*******************************************************************/

static double tourLength (void) {

  int i;
  double cost = 0;

  for (i=0; i<ncities-1; ++i) cost += distance (tour[i], tour[i+1]);
  cost += distance (tour[ncities-1], tour[0]);
  return cost;
}

/*******************************************************************/

static void writeTour (char *s) {

  int i;
  FILE *out = fopen (s, "w");

  if (out == 0) fatal ("Can't open tour output file");
  for (i=0; i<ncities; ++i) fprintf (out, "%d\n", tour[i]);
  fclose (out);
}

/*******************************************************************/

int main (int argc, char **argv) {

  clock_t start;

  if (argc < 3) fatal ("Usage:  nn instance_file tour_output_file [seed]");

  start = clock();
  readInstance (argv[1]);
  printf ("read_time=%.2f\n", (double)(clock() - start) / CLOCKS_PER_SEC);

  tour = (int *) calloc (ncities, sizeof(int));

  start = clock();
  solve();
  printf ("solve_time=%.2f\n", (double)(clock() - start) / CLOCKS_PER_SEC);
  printf ("tour_cost=%0.0f\n", tourLength());

  writeTour (argv[2]);
  return 0;
}
