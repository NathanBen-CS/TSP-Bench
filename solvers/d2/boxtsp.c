#include <math.h>
#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <time.h>

#include "Instance.h"
#include "Utilities.h"
#include "tour.h"
#include "hilbert.h"
#include "neighbors.h"
#include "greedystart.h"
#include "localsearch.h"
#include "boxes.h"
#include "curves.h"
#include "merges.h"
#include "rings2d.h"

/*********************************************************************

Box-TSP solver and its comparison variants (one program, so every variant
runs exactly the same code).  A Karp-style partitioning heuristic: no
worst-case guarantee.

Usage (benchmark contract, options after the two file names):

	boxtsp instance_file tour_output_file [seed] [options]

	--start   boxes | greedy | hilbert | curve | hull   start tour (default boxes)
	--cleanup none | 2opt | oropt          improvement       (default oropt)
	--curve   hilbert | moore | sierpinski | knopp | peano | morton | strip
	          curve the boxes follow (default hilbert); with --start curve,
	          the curve's own order with no boxes
	--merge   M1 | M2 | M3 | M4 | M5 | H0  how box loops are joined
	          (default M1, see merges.h)
	--split   tree | rings   boxes = curve cells (default), or hull rings
	          (convex layers, d2/rings2d.c): no t, no Held-Karp
	--rounds  n | 1          cleanup until nothing improves (default), or
	          one round of the work queue (R1, as twoopt / twoopt_or)
	--check   debug: assert that every ring is convex (and, up to 20,000
	          cities, that later rings lie strictly inside it)
	--t T     max cities per box, 2..12                      (default 10)
	--k K     neighbour-list size for greedy, M3-M5 and cleanup (default 16)

The defaults (Hilbert, M1) give exactly the tours of the original version.
--start hull is hull order alone: the rings in order, each entered at its
city nearest to where the ring before ended, no merge optimisation.
The hull options print extra key=value lines; the existing options print
exactly what they printed before.

The plan's first variants:

	hilbert-only        --start hilbert --cleanup none
	boxes-only          --start boxes   --cleanup none
	boxes+2opt          --start boxes   --cleanup 2opt
	boxes+2opt+oropt    --start boxes   --cleanup oropt      (the full method)
	greedy+2opt+oropt   --start greedy  --cleanup oropt      (baseline)

Cleanup stages (identical for every start):
	2opt   2-opt until a whole round with every city re-examined finds no
	       improving move (localSearchToOptimum in localsearch.c)
	oropt  the 2opt stage, then 2-opt + Or-opt together until neither
	       improves
Both use the bench's neighbour lists, work queue and reversal cap.

Deterministic: the seed is accepted and ignored.  Geometric instances
only (the curves need coordinates).

Prints key=value lines that the benchmark runner stores as extra columns.

 *********************************************************************/

#define GREEDY_CAND 10           /* same greedy start as twoopt / twoopt_or */

static double secondsSince (clock_t start) {
  return (double)(clock() - start) / CLOCKS_PER_SEC;
}

static void usage (void) {
  fatal ("Usage:  boxtsp instance_file tour_output_file [seed] "
	 "[--start boxes|greedy|hilbert|curve] [--cleanup none|2opt|oropt] "
	 "[--curve NAME] [--merge M1..M5|H0] [--split tree|rings] [--t T] [--k K] "
	 "[--rounds n|1] [--check]");
}

int main (int argc, char **argv) {

  const char *startMode = "boxes", *cleanup = "oropt", *split = "tree", *roundsOpt = "n";
  int t = 10, k = 16, i, rounds, curve = CURVE_HILBERT, merge = MERGE_M1, check = 0, rings, hull;
  RingStats rs;
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
    else if (!strcmp (argv[i], "--t")       && i + 1 < argc) t = atoi (argv[++i]);
    else if (!strcmp (argv[i], "--k")       && i + 1 < argc) k = atoi (argv[++i]);
    else if (!strcmp (argv[i], "--split")   && i + 1 < argc) split = argv[++i];
    else if (!strcmp (argv[i], "--rounds")  && i + 1 < argc) roundsOpt = argv[++i];
    else if (!strcmp (argv[i], "--check")) check = 1;
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
      && strcmp (startMode, "hilbert") && strcmp (startMode, "curve")
      && strcmp (startMode, "hull")) usage();
  if (strcmp (split, "tree") && strcmp (split, "rings")) usage();
  if (strcmp (roundsOpt, "n") && strcmp (roundsOpt, "1")) usage();
  hull = !strcmp (startMode, "hull");
  rings = hull || (!strcmp (split, "rings") && !strcmp (startMode, "boxes"));
  if (strcmp (cleanup, "none") && strcmp (cleanup, "2opt")
      && strcmp (cleanup, "oropt")) usage();
  printf ("t=%d\nk=%d\ncurve=%s\nmerge=%s\n", t, k, curveInfo[curve].name, mergeName[merge]);
  if (rings) printf ("split=rings\n");                   /* new options only: the */
  if (!strcmp (roundsOpt, "1")) printf ("rounds=1\n");   /*   old output is unchanged */

  /* ---- 2. read the instance (codes/Instance.c) ---- */
  clk = clock();
  readInstance (argv[1]);
  printf ("read_time=%.2f\n", secondsSince (clk));
  if (!isEuclidean) fatal ("boxtsp needs city coordinates (matrix instance given)");
  tourAlloc();

  /* ---- 3. neighbour lists: needed by the greedy start, merges M3-M5
	    and the cleanup ---- */
  if (!strcmp (startMode, "greedy") || strcmp (cleanup, "none")
      || (!strcmp (startMode, "boxes") && mergeNeedsNeighbors (merge))) {
    clk = clock();
    buildNeighbors (k);
    printf ("neighbor_time=%.2f\n", secondsSince (clk));
  }

  /* ---- 4. start tour ---- */
  if (rings) {                                     /* hull rings as the boxes */
    ringBuild (&bs, &st, &rs);
    printf ("sort_time=%.2f\nring_time=%.2f\n", rs.sortTime, rs.peelTime);
    printf ("rings=%d\nring_mean=%.2f\nring_largest=%d\nring_small=%d\n",
	    rs.rings, (double) ncities / rs.rings, rs.largest, rs.small);
    printf ("rings_per_n23=%.4f\n", rs.rings / pow ((double) ncities, 2.0 / 3.0));
    printf ("ring_edge_mean=%.1f\n", rs.edges ? (double) st.loopsSum / rs.edges : 0.0);
    if (nbList) {                                  /* mean nearest-neighbour distance */
      double nn = 0;
      for (i = 1; i <= ncities; ++i) nn += distance (i, neighbors (i)[0]);
      printf ("nn_mean=%.1f\n", nn / ncities);
    }
    if (check) {
      int f = ringCheck (&bs, ncities <= 20000);
      printf ("ring_check_failures=%d\n", f);
      if (f) fatal ("ring check failed");
    }
    if (hull) hullOrderTour (&bs);
    else {
      clk = clock();
      mergeLoopsAdded (merge, &bs, &st, &added);
      st.glueTime = secondsSince (clk);
      printf ("glue_time=%.2f\nloops_sum=%lld\nglue_G=%lld\n", st.glueTime, st.loopsSum, st.glueG);
      printf ("join_cost=%lld\nmerge_added=%lld\n", st.joinCost, added);
    }
    boxFree (&bs);
  }
  else if (!strcmp (startMode, "boxes")) {
    boxBuild (curve, t, &bs, &st);                 /* sort, boxes, exact loops */
    clk = clock();
    mergeLoopsAdded (merge, &bs, &st, &added);     /* join them into one tour  */
    st.glueTime = secondsSince (clk);
    boxFree (&bs);
    printf ("sort_time=%.2f\nsplit_time=%.2f\ndp_time=%.2f\nglue_time=%.2f\n",
	    st.sortTime, st.splitTime, st.dpTime, st.glueTime);
    printf ("boxes=%d\ncapped_boxes=%d\nlargest_box=%d\n",
	    st.leaves, st.capped, st.largest);
    printf ("loops_sum=%lld\nglue_G=%lld\n", st.loopsSum, st.glueG);
    printf ("join_cost=%lld\nmerge_added=%lld\n", st.joinCost, added);
  }
  else if (!strcmp (startMode, "greedy")) {
    clk = clock();
    greedyTour (GREEDY_CAND);
    printf ("greedy_time=%.2f\n", secondsSince (clk));
  }
  else if (!strcmp (startMode, "hilbert")) {       /* Hilbert order, no boxes */
    clk = clock();
    hilbertSort (BOX_BITS, tour, 0);               /* same curve as the boxes */
    tourSetPositions();
    printf ("sort_time=%.2f\n", secondsSince (clk));
  }
  else {                                           /* --curve's order, no boxes */
    clk = clock();
    curveSort (curve, tour, 0);
    tourSetPositions();
    printf ("sort_time=%.2f\n", secondsSince (clk));
  }
  startLength = tourLength();
  printf ("start_length=%0.0f\n", startLength);
  if (!strcmp (startMode, "boxes")) {
    printf ("glue_share=%.4f\n", st.glueG / startLength);    /* join edges / tour */
    printf ("join_share=%.4f\n", st.joinCost / startLength); /* net join cost     */
  }

  /* ---- 5. cleanup ---- */
  if (strcmp (cleanup, "none") && !strcmp (roundsOpt, "1")) {     /* R1 */
    clk = clock();
    localSearch (!strcmp (cleanup, "oropt"), &m2, &mo);
    printf ("search_time=%.2f\nmoves_2opt=%lld\nmoves_oropt=%lld\n", secondsSince (clk), m2, mo);
  }
  else if (strcmp (cleanup, "none")) {
    clk = clock();
    localSearchToOptimum (0, &m2, &mo, &rounds);          /* 2-opt only */
    printf ("twoopt_time=%.2f\ntwoopt_moves=%lld\ntwoopt_rounds=%d\n",
	    secondsSince (clk), m2, rounds);
    printf ("after_2opt_length=%0.0f\n", tourLength());
  }
  if (!strcmp (cleanup, "oropt") && strcmp (roundsOpt, "1")) {
    clk = clock();
    localSearchToOptimum (1, &m2, &mo, &rounds);          /* 2-opt + Or-opt */
    printf ("oropt_time=%.2f\noropt_moves=%lld\noropt_2opt_moves=%lld\noropt_rounds=%d\n",
	    secondsSince (clk), mo, m2, rounds);
  }

  /* ---- 6. report the length and write the tour, one city per line ---- */
  printf ("tour_cost=%0.0f\n", tourLength());
  writeTour (argv[2]);
  return 0;
}
