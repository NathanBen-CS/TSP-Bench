#include <stdio.h>
#include <time.h>

#include "Instance.h"
#include "Utilities.h"
#include "tour.h"
#include "hilbert.h"

/*********************************************************************

Space-filling curve heuristic (see hilbert.c for how it works).

Usage (benchmark contract):

	spacefill instance_file tour_output_file [seed]

Deterministic: the seed is accepted and ignored.  Geometric instances only
(EUC_2D / CEIL_2D); matrix instances stop with an error, and the benchmark
config marks this solver "geometric_only" so they are skipped.

Prints key=value lines picked up by the benchmark runner.

 *********************************************************************/

static double secondsSince (clock_t start) {
  return (double)(clock() - start) / CLOCKS_PER_SEC;
}

int main (int argc, char **argv) {

  clock_t start;

  if (argc < 3) fatal ("Usage:  spacefill instance_file tour_output_file [seed]");

  /* 1. read the instance (codes/Instance.c) */
  start = clock();
  readInstance (argv[1]);
  printf ("read_time=%.2f\n", secondsSince (start));

  /* 2. order the cities along the Hilbert curve: that order is the tour */
  start = clock();
  tourAlloc();
  hilbertTour();
  printf ("solve_time=%.2f\n", secondsSince (start));

  /* 3. report the length and write the tour, one city per line */
  printf ("tour_cost=%0.0f\n", tourLength());
  writeTour (argv[2]);
  return 0;
}
