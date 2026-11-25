#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <time.h>

#include "Instance.h"
#include "Utilities.h"
#include "tour.h"
#include "neighbors.h"
#include "greedystart.h"
#include "localsearch.h"
#include "boxes.h"
#include "curves3d.h"
#include "kdboxes3d.h"
#include "merges.h"
#include "diag.h"

/*********************************************************************

Box-TSP solver and its comparison variants, 3D (the 2D program is
d2/boxtsp.c; options, stages and printed columns are the same).  A
Karp-style partitioning heuristic: no worst-case guarantee.

	boxtsp3d instance_file tour_output_file [seed] [options]

	--start   boxes | greedy | curve       start tour       (default boxes)
	--cleanup none | 2opt | oropt          improvement      (default oropt)
	--rounds  n | 1                        cleanup length   (default n)
	--curve   hilbert3 | morton3 | peano3 | snake3 | moore3 | sierpinski3
	--merge   M1 | M2 | M3 | M4 | M5                        (default M1)
	--split   tree | kd                    curve cells, or kd-median boxes
	--t T     max cities per box, 2..12                     (default 10)
	--k K     neighbour-list size for greedy, M3-M5 and cleanup (default 16)
	--diag    also print the box diagnostics (diag.c); not for timed runs

Cleanup, exactly the 2D stages:
	--rounds n   Rn: 2opt  = 2-opt rounds until a whole round moves nothing
	                 oropt = that, then 2-opt + Or-opt rounds until neither moves
	--rounds 1   R1: one round of the work queue, as twoopt (2opt) and
	                 twoopt_or (oropt: 2-opt and Or-opt mixed from the start)
R1 after a box start is new in 3D (in 2D only greedy had an R1 program).

What --k controls: the neighbour lists of the local search and of merges
M3-M5, and greedy's round-1 candidates, which are the first min(10, K)
entries of the same lists.  It does not change greedy's later join rounds
(8 nearest fragment ends), M1, M2 or the boxes.

Deterministic: the seed is accepted and ignored.

 *********************************************************************/

#define GREEDY_CAND 10           /* same greedy start as twoopt / twoopt_or */

static double secondsSince (clock_t start) {
  return (double)(clock() - start) / CLOCKS_PER_SEC;
}

static void usage (void) {
  fatal ("Usage:  boxtsp3d instance_file tour_output_file [seed] "
	 "[--start boxes|greedy|curve] [--cleanup none|2opt|oropt] [--rounds n|1] "
	 "[--curve NAME] [--merge M1..M5] [--split tree|kd] [--t T] [--k K] [--diag]");
}

static void build (int kd, int curve, int t, BoxSet *bs, BoxStats *st) {
  if (kd) boxBuildKd (curve, t, bs, st);
  else    boxBuild (curve, t, bs, st);
}

int main (int argc, char **argv) {

  const char *startMode = "boxes", *cleanup = "oropt", *split = "tree", *roundsOpt = "n";
  int t = 10, k = 16, i, rounds, curve = CURVE_HILBERT, merge = MERGE_M1, diag = 0, kd, boxes;
  long long m2, mo, added = 0;
  double startLength;
  BoxStats st;
  BoxSet bs;
  clock_t clk;

  /* ---- 1. arguments: two files, optional seed, then --options ---- */
  if (argc < 3) usage();
  for (i = 3; i < argc; ++i) {
    if      (!strcmp (argv[i], "--start")   && i + 1 < argc) startMode = argv[++i];
    else if (!strcmp (argv[i], "--cleanup") && i + 1 < argc) cleanup   = argv[++i];
    else if (!strcmp (argv[i], "--rounds")  && i + 1 < argc) roundsOpt = argv[++i];
    else if (!strcmp (argv[i], "--split")   && i + 1 < argc) split     = argv[++i];
    else if (!strcmp (argv[i], "--t")       && i + 1 < argc) t = atoi (argv[++i]);
    else if (!strcmp (argv[i], "--k")       && i + 1 < argc) k = atoi (argv[++i]);
    else if (!strcmp (argv[i], "--diag"))                    diag = 1;
    else if (!strcmp (argv[i], "--curve")   && i + 1 < argc) {
      if ((curve = curveByName (argv[++i])) < 0) usage();
    }
    else if (!strcmp (argv[i], "--merge")   && i + 1 < argc) {
      if ((merge = mergeByName (argv[++i])) < 0) usage();
    }
    else if (argv[i][0] == '-') usage();
    /* anything else is the seed: accepted and ignored (deterministic) */
  }
  if (t < 2 || t > BOX_MAXT) fatal ("--t must be between 2 and 12");
  if (k < 1) fatal ("--k must be at least 1");
  if (strcmp (startMode, "boxes") && strcmp (startMode, "greedy")
      && strcmp (startMode, "curve")) usage();
  if (strcmp (cleanup, "none") && strcmp (cleanup, "2opt")
      && strcmp (cleanup, "oropt")) usage();
  if (strcmp (roundsOpt, "n") && strcmp (roundsOpt, "1")) usage();
  if (strcmp (split, "tree") && strcmp (split, "kd")) usage();
  kd = !strcmp (split, "kd");
  boxes = !strcmp (startMode, "boxes");
  printf ("t=%d\nk=%d\ncurve=%s\nmerge=%s\nsplit=%s\nrounds=%s\n",
	  t, k, curveInfo[curve].name, mergeName[merge], split, roundsOpt);

  /* ---- 2. read the instance (d3/instance3d.c) ---- */
  clk = clock();
  readInstance (argv[1]);
  printf ("read_time=%.2f\n", secondsSince (clk));
  tourAlloc();

  /* ---- 3. neighbour lists: needed by the greedy start, merges M3-M5
	    and the cleanup ---- */
  if (!strcmp (startMode, "greedy") || strcmp (cleanup, "none") || diag
      || (boxes && mergeNeedsNeighbors (merge))) {
    clk = clock();
    buildNeighbors (k);
    printf ("neighbor_time=%.2f\n", secondsSince (clk));
  }

  /* ---- 4. start tour ---- */
  if (boxes) {
    build (kd, curve, t, &bs, &st);                /* sort, boxes, exact loops */
    clk = clock();
    mergeLoopsAdded (merge, &bs, &st, &added);     /* join them into one tour  */
    st.glueTime = secondsSince (clk);
    printf ("sort_time=%.2f\nsplit_time=%.2f\ndp_time=%.2f\nglue_time=%.2f\n",
	    st.sortTime, st.splitTime, st.dpTime, st.glueTime);
    printf ("boxes=%d\ncapped_boxes=%d\nlargest_box=%d\n",
	    st.leaves, st.capped, st.largest);
    printf ("loops_sum=%lld\nglue_G=%lld\n", st.loopsSum, st.glueG);
    printf ("join_cost=%lld\nmerge_added=%lld\n", st.joinCost, added);
    if (!diag) boxFree (&bs);
  }
  else if (!strcmp (startMode, "greedy")) {
    clk = clock();
    greedyTour (GREEDY_CAND);
    printf ("greedy_time=%.2f\n", secondsSince (clk));
  }
  else {                                           /* --curve's order, no boxes */
    clk = clock();
    curveSort (curve, tour, 0);
    tourSetPositions();
    printf ("sort_time=%.2f\n", secondsSince (clk));
  }
  startLength = tourLength();
  printf ("start_length=%0.0f\n", startLength);
  if (boxes) {
    printf ("glue_share=%.4f\n", st.glueG / startLength);    /* join edges / tour */
    printf ("join_share=%.4f\n", st.joinCost / startLength); /* net join cost     */
  }

  /* ---- 5. cleanup ---- */
  if (strcmp (cleanup, "none") && !strcmp (roundsOpt, "1")) {     /* R1 */
    clk = clock();
    localSearch (!strcmp (cleanup, "oropt"), &m2, &mo);
    printf ("search_time=%.2f\nmoves_2opt=%lld\nmoves_oropt=%lld\n", secondsSince (clk), m2, mo);
  }
  else if (strcmp (cleanup, "none")) {                            /* Rn */
    clk = clock();
    localSearchToOptimum (0, &m2, &mo, &rounds);          /* 2-opt only */
    printf ("twoopt_time=%.2f\ntwoopt_moves=%lld\ntwoopt_rounds=%d\n",
	    secondsSince (clk), m2, rounds);
    printf ("after_2opt_length=%0.0f\n", tourLength());
    if (!strcmp (cleanup, "oropt")) {
      clk = clock();
      localSearchToOptimum (1, &m2, &mo, &rounds);        /* 2-opt + Or-opt */
      printf ("oropt_time=%.2f\noropt_moves=%lld\noropt_2opt_moves=%lld\noropt_rounds=%d\n",
	      secondsSince (clk), mo, m2, rounds);
    }
  }

  /* ---- 6. report the length and write the tour, one city per line ---- */
  printf ("tour_cost=%0.0f\n", tourLength());
  writeTour (argv[2]);

  if (diag) {                                      /* boxes against the final tour */
    if (!boxes) build (kd, curve, t, &bs, &st);
    printDiagnostics (&bs, &st);
  }
  return 0;
}
