#include <stdio.h>
#include <time.h>

#include "Instance.h"
#include "Utilities.h"
#include "tour.h"
#include "neighbors.h"
#include "greedystart.h"
#include "localsearch.h"

/*********************************************************************

2-opt local search with neighbour lists and a work queue
(see localsearch.c for the moves and speed tricks).

Usage (benchmark contract):

	twoopt instance_file tour_output_file [seed]

Deterministic: the seed is accepted and ignored.

Steps:
  1. read the instance
  2. build each city's 16-nearest-neighbour candidate list
  3. build a greedy start tour from those lists (greedystart.c)
  4. apply improving 2-opt moves until none is left
  5. write the tour

Prints key=value lines picked up by the benchmark runner.

 *********************************************************************/

#define NEIGHBORS   16      /* candidate partners per city in local search  */
#define GREEDY_CAND 10      /* candidate edges per city for the greedy start */

static double secondsSince (clock_t start) {
  return (double)(clock() - start) / CLOCKS_PER_SEC;
}

int main (int argc, char **argv) {

  clock_t start;
  long long moves2, movesOr;

  if (argc < 3) fatal ("Usage:  twoopt instance_file tour_output_file [seed]");

  /* 1. read the instance (codes/Instance.c) */
  start = clock();
  readInstance (argv[1]);
  printf ("read_time=%.2f\n", secondsSince (start));

  /* 2. candidate lists: each city's NEIGHBORS nearest cities */
  start = clock();
  buildNeighbors (NEIGHBORS);
  printf ("neighbor_time=%.2f\n", secondsSince (start));

  /* 3. greedy start tour, built from the candidate lists */
  start = clock();
  tourAlloc();
  greedyTour (GREEDY_CAND);
  printf ("start_time=%.2f\n", secondsSince (start));
  printf ("start_length=%0.0f\n", tourLength());

  /* 4. 2-opt until no improving move is left (Or-opt off) */
  start = clock();
  localSearch (0, &moves2, &movesOr);
  printf ("search_time=%.2f\n", secondsSince (start));
  printf ("moves_2opt=%lld\n", moves2);

  /* 5. report the length and write the tour, one city per line */
  printf ("tour_cost=%0.0f\n", tourLength());
  writeTour (argv[2]);
  return 0;
}
