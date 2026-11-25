#include "testlib3d.h"
#include "tour.h"
#include "neighbors.h"
#include "greedystart.h"
#include "localsearch.h"
#include "boxes.h"
#include "curves3d.h"
#include "kdboxes3d.h"
#include "merges.h"

/*********************************************************************

Tests for the box solver's pieces in 3D: the 3D version of boxtest2d.c.
(Reader, neighbour lists and the curves themselves: unit3d.c.)

	boxtest3d          run every test, print PASS/FAIL per test,
	                   exit code 0 only if all pass

 *********************************************************************/

/* ---------------------------------------------------------------- 1. boxes */

/* Every box is exactly one curve cell (all its keys share the cell's prefix,
   no other city falls in that cell), size <= t unless the depth cap was hit. */
static void testCurveBoxes (void) {

  int c, type, ti, ts[] = { 4, 10, 12 }, n = 5000, nb, j, k, t, ok, capped;
  int *order, *start, *level;
  unsigned long long *keys, size, cell;
  char name[160], msg[160];

  for (c = 0; c < NCURVES3; ++c)
    for (type = 0; type < NTYPES; ++type)
      for (ti = 0; ti < 3; ++ti) {
	t = ts[ti];
	makeInstance (n, type);
	order = (int *) malloc (n * sizeof(int));
	keys = (unsigned long long *) malloc (n * sizeof *keys);
	start = (int *) malloc ((n + 1) * sizeof(int));
	level = (int *) malloc ((n + 1) * sizeof(int));
	curveSort (c, order, keys);
	nb = boxSplitBase (keys, n, curveInfo[c].base, curveInfo[c].digits, t, start, level);
	ok = isPermutation (order, n) && start[0] == 0 && start[nb] == n;
	capped = 0;
	for (j = 0; j < nb && ok; ++j) {
	  if (start[j + 1] <= start[j]) ok = 0;
	  size = 1;
	  for (k = 0; k < curveInfo[c].digits - level[j]; ++k) size *= curveInfo[c].base;
	  cell = keys[start[j]] / size;
	  for (k = start[j]; k < start[j + 1]; ++k)
	    if (keys[k] / size != cell) ok = 0;
	  if (start[j] > 0 && keys[start[j] - 1] / size == cell) ok = 0;
	  if (start[j + 1] < n && keys[start[j + 1]] / size == cell) ok = 0;
	  if (start[j + 1] - start[j] > t) {
	    capped++;
	    if (level[j] != curveInfo[c].digits) ok = 0;   /* only at the depth cap */
	  }
	}
	sprintf (name, "%s boxes are whole curve cells, size <= t (%s, t=%d)", curveInfo[c].name, typeName[type], t);
	sprintf (msg, "%d boxes, %d at the depth cap", nb, capped);
	report (name, ok, msg);
	free (order); free (keys); free (start); free (level);
      }
}

/* A depth-capped box keeps the curve order of its cities. */
static void testCapped (void) {
  int n = 5000, j, k, ok = 1, capped = 0, curve;
  BoxSet bs;
  BoxStats st;
  char msg[120], name[120];
  for (curve = 0; curve < NCURVES3; ++curve) {
    makeInstance (n, DUPLICATES);
    boxBuild (curve, 10, &bs, &st);
    ok = 1;
    capped = 0;
    for (j = 0; j < bs.nb; ++j)
      if (bs.start[j + 1] - bs.start[j] > 10) {
	capped++;
	for (k = bs.start[j]; k < bs.start[j + 1]; ++k) if (bs.loops[k] != bs.order[k]) ok = 0;
      }
    sprintf (name, "%s: depth-capped boxes keep curve order (duplicates)", curveInfo[curve].name);
    sprintf (msg, "%d capped boxes, largest %d", capped, st.largest);
    report (name, ok && capped > 0 && capped == st.capped, msg);
    boxFree (&bs);
  }
}

/* kd boxes: every city in one box, every box but one holds exactly t. */
static void testKd (void) {
  int type, ti, ts[] = { 4, 10, 12 }, n = 5003, j, full, ok, s;
  BoxSet bs;
  BoxStats st;
  char name[160], msg[160];
  for (type = 0; type < NTYPES; ++type)
    for (ti = 0; ti < 3; ++ti) {
      makeInstance (n, type);
      boxBuildKd (CURVE_HILBERT3, ts[ti], &bs, &st);
      ok = isPermutation (bs.order, n) && isPermutation (bs.loops, n) && bs.start[0] == 0 && bs.start[bs.nb] == n;
      full = 0;
      for (j = 0; j < bs.nb; ++j) {
	s = bs.start[j + 1] - bs.start[j];
	if (s < 1 || s > ts[ti]) ok = 0;
	if (s == ts[ti]) full++;
      }
      sprintf (name, "kd boxes: all boxes but at most one hold exactly t (%s, t=%d)", typeName[type], ts[ti]);
      sprintf (msg, "%d boxes, %d full", bs.nb, full);
      report (name, ok && full >= bs.nb - 1 && st.capped == 0, msg);
      boxFree (&bs);
    }
}

/* ---------------------------------------------------------------- 2. DP */

static long long bfBest;
static void permute (int *p, int k, int s) {          /* brute force, city p[0] fixed */
  int i, tmp;
  long long L;
  if (k == s) {
    L = cycleLength (p, s);
    if (L < bfBest) bfBest = L;
    return;
  }
  for (i = k; i < s; ++i) {
    tmp = p[k]; p[k] = p[i]; p[i] = tmp;
    permute (p, k + 1, s);
    tmp = p[k]; p[k] = p[i]; p[i] = tmp;
  }
}

static void testDP (void) {

  int trial, s, c, cities[BOX_MAXT], loop[BOX_MAXT], p[BOX_MAXT], ok = 1, bad = 0;
  long long L;
  char msg[160];

  for (trial = 0; trial < 400; ++trial) {
    s = 1 + trial % 8;                            /* 1..8 cities */
    makeInstance (s, UNIFORM);
    for (c = 1; c <= s; ++c) {                    /* small coordinates, and */
      xx[c] = rnd() % 1000;                       /*   a duplicate now and  */
      yy[c] = rnd() % 1000;                       /*   then                 */
      zz[c] = rnd() % 1000;
    }
    if (trial % 7 == 0 && s > 1) { xx[2] = xx[1]; yy[2] = yy[1]; zz[2] = zz[1]; }
    for (c = 0; c < s; ++c) cities[c] = c + 1;
    L = heldKarpLoop (cities, s, loop);
    for (c = 0; c < s; ++c) p[c] = c + 1;
    bfBest = 0x7fffffffffffffffLL;                /* "infinity" */
    permute (p, 1, s);
    if (L != bfBest || !isPermutation (loop, s) || cycleLength (loop, s) != L) {
      ok = 0;
      if (bad++ < 3) printf ("      s=%d dp=%lld brute=%lld\n", s, L, bfBest);
    }
  }
  sprintf (msg, "400 random 3D sets of 1-8 cities, %d mismatches", bad);
  report ("Held-Karp DP = brute force, and its loop has that length", ok, msg);
}

/* ---------------------------------------------------------------- 3. merges */

/* Two unit cubes two apart.  A cube's best loop is 8 edges; joining them costs
   2 + 2 and saves 1 + 1, so the best merge gives 8 + 8 - 2 + 4 = 18. */
static void testTwoCubes (void) {

  int loop[8][3] = { {0,0,0},{0,1,0},{0,1,1},{0,0,1},{1,0,1},{1,1,1},{1,1,0},{1,0,0} };
  int order[16], start[3] = { 0, 8, 16 }, m, c;
  BoxSet bs;
  BoxStats st;
  long long added, L;
  char name[120], msg[160];

  makeInstance (16, UNIFORM);
  for (c = 0; c < 16; ++c) {
    xx[c+1] = 1000 * loop[c % 8][0] + (c >= 8 ? 3000 : 0);
    yy[c+1] = 1000 * loop[c % 8][1];
    zz[c+1] = 1000 * loop[c % 8][2];
    order[c] = c + 1;
  }
  bs.nb = 2; bs.order = order; bs.start = start; bs.loops = order;
  buildNeighbors (15);
  for (m = 0; m < NMERGES; ++m) {
    tourAlloc();
    st.loopsSum = 16000;
    mergeLoopsAdded (m, &bs, &st, &added);
    L = cycleLength (tour, 16);
    sprintf (name, "%s on two unit cubes two apart", mergeName[m]);
    sprintf (msg, "length %.0f (best merge 18), added %.0f, recomputed %.0f",
	     L / 1000.0, added / 1000.0, (L - 16000) / 1000.0);
    if (m == MERGE_M1)                                /* first cut is fixed */
      report (name, isPermutation (tour, 16) && added == L - 16000, msg);
    else
      report (name, isPermutation (tour, 16) && L == 18000 && added == L - 16000, msg);
    free (tour); free (pos);
  }
  free (nbList);
}

/* Every curve, merge and kind of data: one valid tour, and the extra length
   the method says it added equals the recomputed tour length minus the loops. */
static void testMergesValid (void) {

  int curve, type, m, n = 5000, kd;
  long long added, L, m1 = 0;
  BoxSet bs;
  BoxStats st;
  char name[160], msg[200];

  for (kd = 0; kd < 2; ++kd)
    for (curve = 0; curve < NCURVES3; ++curve) {
      if (kd && curve != CURVE_HILBERT3) continue;
      for (type = 0; type < NTYPES; ++type) {
	makeInstance (n, type);
	buildNeighbors (16);
	if (kd) boxBuildKd (curve, 10, &bs, &st); else boxBuild (curve, 10, &bs, &st);
	for (m = 0; m < NMERGES; ++m) {
	  tourAlloc();
	  mergeLoopsAdded (m, &bs, &st, &added);
	  L = cycleLength (tour, n);
	  if (m == MERGE_M1) m1 = L;
	  sprintf (name, "%s%s + %s (%s): valid tour, added cost = recomputed",
		   kd ? "kd boxes, " : "", curveInfo[curve].name, mergeName[m], typeName[type]);
	  sprintf (msg, "length %lld, added %lld, recomputed %lld", L, added, L - st.loopsSum);
	  report (name, isPermutation (tour, n) && added == L - st.loopsSum && L <= st.loopsSum + st.glueG
		  && L == (long long) tourLength() && (m != MERGE_M2 || L <= m1), msg);
	  free (tour); free (pos);
	}
	boxFree (&bs);
	free (nbList);
      }
    }
}

/* Tiny inputs: 2 to 30 cities, every curve, merge and several t, then the
   full cleanup.  Must give a valid tour and never crash. */
static void testTiny (void) {

  int curve, type, m, n, ti, ts[] = { 2, 4, 10 }, rounds, bad = 0, runs = 0;
  long long m2, mo, added;
  BoxSet bs;
  BoxStats st;
  char msg[160];

  for (n = 2; n <= 30; ++n)
    for (curve = 0; curve < NCURVES3; ++curve)
      for (type = 0; type < NTYPES; ++type)
	for (ti = 0; ti < 3; ++ti) {
	  makeInstance (n, type);
	  buildNeighbors (16);
	  boxBuild (curve, ts[ti], &bs, &st);
	  for (m = 0; m < NMERGES; ++m) {
	    tourAlloc();
	    mergeLoopsAdded (m, &bs, &st, &added);
	    if (!isPermutation (tour, n) || added != cycleLength (tour, n) - st.loopsSum) bad++;
	    localSearchToOptimum (0, &m2, &mo, &rounds);
	    localSearchToOptimum (1, &m2, &mo, &rounds);
	    if (!isPermutation (tour, n)) bad++;
	    runs++;
	    free (tour); free (pos);
	  }
	  boxFree (&bs);
	  free (nbList);
	}
  sprintf (msg, "%d runs, %d bad", runs, bad);
  report ("tiny inputs (2-30 cities), every curve, merge, t = 2, 4, 10, with cleanup", bad == 0, msg);
}

/* M2 is the best possible snake: brute force over every cut and direction. */
static long long bfSnake;
static int snakeR[64], snakeD[64];
static void snakeAll (const BoxSet *bs, int j, int *tmp) {
  int s, r, d, k, tn = 0, b, lo;
  if (j == bs->nb) {
    for (b = 0; b < bs->nb; ++b) {
      lo = bs->start[b]; s = bs->start[b + 1] - lo;
      for (k = 0; k < s; ++k) tmp[tn++] = bs->loops[lo + ((snakeR[b] + snakeD[b] * k) % s + s) % s];
    }
    if (cycleLength (tmp, tn) < bfSnake) bfSnake = cycleLength (tmp, tn);
    return;
  }
  s = bs->start[j + 1] - bs->start[j];
  for (r = 0; r < s; ++r)
    for (d = -1; d <= 1; d += 2) {
      snakeR[j] = r; snakeD[j] = d;
      snakeAll (bs, j + 1, tmp);
      if (s == 1) return;
    }
}

static void testM2Optimal (void) {
  int trial, n, bad = 0, done = 0, tmp[64];
  long long L;
  BoxSet bs;
  BoxStats st;
  char msg[160];
  for (trial = 0; trial < 80; ++trial) {
    n = 10 + trial % 7;
    makeInstance (n, trial % 2 ? CLUSTERED : UNIFORM);
    boxBuild (CURVE_HILBERT3, 6, &bs, &st);
    if (bs.nb > 6) { boxFree (&bs); continue; }          /* keep brute force small */
    tourAlloc();
    mergeLoops (MERGE_M2, &bs, &st);
    L = cycleLength (tour, n);
    bfSnake = 0x7fffffffffffffffLL;
    snakeAll (&bs, 0, tmp);
    if (L != bfSnake) bad++;
    done++;
    free (tour); free (pos);
    boxFree (&bs);
  }
  sprintf (msg, "%d small instances, %d where M2 missed the best snake", done, bad);
  report ("M2 = best snake by brute force over every cut and direction", bad == 0 && done > 0, msg);
}

/* ---------------------------------------------------------------- 4. greedy and cleanup */

static void testGreedy (void) {
  int type, n;
  char name[120];
  for (type = 0; type < NTYPES; ++type)
    for (n = 8; n <= 8000; n *= 10) {
      makeInstance (n, type);
      buildNeighbors (16);
      tourAlloc();
      greedyTour (10);
      sprintf (name, "greedy start: valid tour (%s, %d cities)", typeName[type], n);
      report (name, isPermutation (tour, n), "");
      free (tour); free (pos); free (nbList);
    }
}

static long long improving2opt (void) {        /* brute force over all edge pairs */
  int i, j, a, b, c, d;
  long long count = 0;
  for (i = 0; i < ncities; ++i) {
    a = tour[i]; b = tour[(i + 1) % ncities];
    for (j = i + 2; j < ncities; ++j) {
      c = tour[j]; d = tour[(j + 1) % ncities];
      if (d == a) continue;
      if (distance (a,b) + distance (c,d) - distance (a,c) - distance (b,d) > 0) count++;
    }
  }
  return count;
}

static void testLocalSearch (void) {

  int type, start, n = 300, rounds;
  long long L0, L1, L2, m2, mo, b2;
  BoxStats st;
  char name[120], msg[200];

  for (start = 0; start < 2; ++start)
    for (type = 0; type < NTYPES; ++type) {
      makeInstance (n, type);
      tourAlloc();
      buildNeighbors (n - 1);                 /* complete lists: nothing hidden */
      if (start) greedyTour (10); else boxTour (8, &st);
      L0 = cycleLength (tour, n);
      localSearchToOptimum (0, &m2, &mo, &rounds);       /* 2-opt */
      L1 = cycleLength (tour, n);
      b2 = improving2opt();
      sprintf (name, "2-opt from %s: never longer, brute force finds no improving move (%s)",
	       start ? "greedy" : "boxes", typeName[type]);
      sprintf (msg, "%lld -> %lld in %lld moves, %d rounds; %lld improving moves left", L0, L1, m2, rounds, b2);
      report (name, isPermutation (tour, n) && L1 <= L0 && b2 == 0, msg);
      localSearchToOptimum (1, &m2, &mo, &rounds);       /* 2-opt + Or-opt */
      L2 = cycleLength (tour, n);
      b2 = improving2opt();
      sprintf (name, "Or-opt from %s: never longer, no 2-opt move left (%s)", start ? "greedy" : "boxes", typeName[type]);
      sprintf (msg, "%lld -> %lld (%lld Or-opt moves)", L1, L2, mo);
      report (name, isPermutation (tour, n) && L2 <= L1 && b2 == 0, msg);
      free (tour); free (pos); free (nbList);
    }
}

static void testPipeline (void) {

  int type, n = 20000, rounds;
  long long L, m2, mo;
  BoxStats st;
  char name[120], msg[200];

  for (type = 0; type < NTYPES; ++type) {
    makeInstance (n, type);
    tourAlloc();
    boxTour (10, &st);
    buildNeighbors (16);
    localSearchToOptimum (0, &m2, &mo, &rounds);
    localSearchToOptimum (1, &m2, &mo, &rounds);
    L = cycleLength (tour, n);
    sprintf (name, "full method on 20,000 cities finishes with a valid tour (%s)", typeName[type]);
    sprintf (msg, "length %lld, %d boxes, %d capped, largest box %d", L, st.leaves, st.capped, st.largest);
    report (name, isPermutation (tour, n) && L == (long long) tourLength(), msg);
    free (tour); free (pos); free (nbList);
  }
}

int main (void) {
  testCurveBoxes();
  testCapped();
  testKd();
  testDP();
  testTwoCubes();
  testMergesValid();
  testTiny();
  testM2Optimal();
  testGreedy();
  testLocalSearch();
  testPipeline();
  return finish();
}
