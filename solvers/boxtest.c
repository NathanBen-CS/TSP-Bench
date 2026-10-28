#include <math.h>
#include <stdio.h>
#include <stdlib.h>
#include <string.h>

#include "Instance.h"
#include "Utilities.h"
#include "tour.h"
#include "hilbert.h"
#include "neighbors.h"
#include "localsearch.h"
#include "boxes.h"
#include "curves.h"
#include "merges.h"

/*********************************************************************

Tests for the box solver's pieces (Box-TSP build plan, Task 2).

	boxtest            run every test, print PASS/FAIL per test,
	                   exit code 0 only if all pass

Instances are generated in memory (uniform, clustered and duplicate-heavy)
and fed to the same code the solvers use, with the bench's rounding
(nearest integer, as EUC_2D).

 *********************************************************************/

static int failures = 0;

static void report (const char *name, int ok, const char *detail) {
  printf ("%-4s  %s%s%s\n", ok ? "PASS" : "FAIL", name, detail[0] ? " - " : "", detail);
  if (!ok) failures++;
}

/* ---------------------------------------------------------------- test instances */

static unsigned long long rs = 88172645463325252ULL;          /* xorshift64 */
static unsigned long long rnd (void) { rs ^= rs << 13; rs ^= rs >> 7; rs ^= rs << 17; return rs; }
static double urand (void) { return (rnd() >> 11) * (1.0 / 9007199254740992.0); }
static double nrand (void) {                                   /* standard normal */
  return sqrt (-2 * log (urand() + 1e-300)) * cos (6.283185307179586 * urand());
}

static int roundDist (int a, int b) {        /* same rule as EUC_2D in Instance.c */
  double dx = xx[a] - xx[b], dy = yy[a] - yy[b];
  return (int) (sqrt (dx*dx + dy*dy) + 0.5);
}

enum { UNIFORM, CLUSTERED, DUPLICATES };
static const char *typeName[] = { "uniform", "clustered", "duplicates" };

/* Make an n-city instance in the globals Instance.c normally fills. */
static void makeInstance (int n, int type) {

  int c, k, nc;
  double *cx, *cy, sd;

  free (xx); free (yy);
  ncities = n;
  xx = (double *) calloc (n + 1, sizeof(double));
  yy = (double *) calloc (n + 1, sizeof(double));
  isEuclidean = 1;
  distance = roundDist;

  if (type == UNIFORM) {                       /* like portgen */
    for (c = 1; c <= n; ++c) { xx[c] = rnd() % 1000000; yy[c] = rnd() % 1000000; }
  }
  else if (type == CLUSTERED) {                /* like portcgen: n/100 centres, */
    nc = n / 100 > 0 ? n / 100 : 1;            /* normal spread 1e6/sqrt(n)     */
    cx = (double *) malloc (nc * sizeof(double));
    cy = (double *) malloc (nc * sizeof(double));
    for (k = 0; k < nc; ++k) { cx[k] = urand() * 1e6; cy[k] = urand() * 1e6; }
    sd = 1e6 / sqrt ((double) n);
    for (c = 1; c <= n; ++c) {
      k = rnd() % nc;
      xx[c] = (int) (cx[k] + sd * nrand());
      yy[c] = (int) (cy[k] + sd * nrand());
    }
    free (cx); free (cy);
  }
  else {                                       /* many exact duplicates: */
    nc = n / 25 > 1 ? n / 25 : 2;              /*   ~25 cities per spot  */
    for (c = 1; c <= n; ++c) {
      k = rnd() % nc;
      xx[c] = (k * 7919) % 1000;
      yy[c] = (k * 104729) % 1000;
    }
  }
}

static int isPermutation (const int *a, int n) {   /* cities 1..n, each once */
  int i, ok = 1;
  char *seen = (char *) calloc (n + 1, 1);
  for (i = 0; i < n; ++i) {
    if (a[i] < 1 || a[i] > n || seen[a[i]]) { ok = 0; break; }
    seen[a[i]] = 1;
  }
  free (seen);
  return ok;
}

static long long cycleLength (const int *a, int n) {
  long long L = 0;
  int i;
  for (i = 0; i < n; ++i) L += distance (a[i], a[(i + 1) % n]);
  return L;
}

/* ---------------------------------------------------------------- 1. Hilbert */

static void testHilbert (void) {

  unsigned n = 16, x, y, d;
  int cx[256], cy[256], ok = 1;
  char msg[120];

  for (x = 0; x < n; ++x)                      /* cell at each curve position */
    for (y = 0; y < n; ++y) {
      d = (unsigned) hilbertIndex (x, y, n);
      cx[d] = x; cy[d] = y;
    }
  for (d = 1; d < n * n; ++d)                  /* consecutive -> neighbouring */
    if (abs (cx[d] - cx[d-1]) + abs (cy[d] - cy[d-1]) != 1) ok = 0;
  report ("Hilbert 16x16: consecutive positions are neighbouring cells", ok, "");

  sprintf (msg, "(0,0)=%llu (0,1)=%llu (1,1)=%llu (1,0)=%llu",
	   hilbertIndex (0,0,2), hilbertIndex (0,1,2), hilbertIndex (1,1,2), hilbertIndex (1,0,2));
  report ("Hilbert 2x2 visiting order", hilbertIndex (0,0,2) == 0 && hilbertIndex (0,1,2) == 1
	  && hilbertIndex (1,1,2) == 2 && hilbertIndex (1,0,2) == 3, msg);
}

/* ---------------------------------------------------------------- 2. boxes */

static void testBoxes (void) {

  int type, ti, ts[] = { 6, 8, 10, 12 }, n = 5000, nb, j, k, ok, capped, sorted;
  int *order, *start;
  unsigned long long *keys;
  char name[120], msg[160];

  for (type = 0; type < 3; ++type)
    for (ti = 0; ti < 4; ++ti) {
      makeInstance (n, type);
      order = (int *) malloc (n * sizeof(int));
      keys  = (unsigned long long *) malloc (n * sizeof(unsigned long long));
      start = (int *) malloc ((n + 1) * sizeof(int));
      hilbertSort (BOX_BITS, order, keys);
      nb = boxSplit (keys, n, BOX_BITS, ts[ti], start);

      ok = isPermutation (order, n) && start[0] == 0 && start[nb] == n;
      sorted = 1;
      for (k = 1; k < n; ++k) if (keys[k] < keys[k-1]) sorted = 0;
      capped = 0;
      for (j = 0; j < nb; ++j) {
	if (start[j + 1] <= start[j]) ok = 0;            /* empty box: wrong   */
	if (start[j + 1] - start[j] > ts[ti]) {          /* too big: only OK   */
	  capped++;                                      /*   if every city is */
	  for (k = start[j] + 1; k < start[j + 1]; ++k)  /*   in the same cell */
	    if (keys[k] != keys[start[j]]) ok = 0;
	}
      }
      sprintf (name, "boxes: every city in exactly one box, size <= t (%s, t=%d)",
	       typeName[type], ts[ti]);
      sprintf (msg, "%d boxes, %d over t (depth cap, same cell)", nb, capped);
      report (name, ok && sorted, msg);
      free (order); free (keys); free (start);
    }
}

/* ---------------------------------------------------------------- 3. DP */

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

  for (trial = 0; trial < 300; ++trial) {
    s = 4 + trial % 6;                            /* 4..9 cities */
    makeInstance (s, UNIFORM);
    for (c = 1; c <= s; ++c) {                    /* small coordinates, and */
      xx[c] = rnd() % 1000;                       /*   a duplicate now and  */
      yy[c] = rnd() % 1000;                       /*   then                 */
    }
    if (trial % 7 == 0) { xx[2] = xx[1]; yy[2] = yy[1]; }
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
  sprintf (msg, "300 random sets of 4-9 cities, %d mismatches", bad);
  report ("Held-Karp DP = brute force, and its loop has that length", ok, msg);
}

/* ---------------------------------------------------------------- 4. glue */

static void testGlue (void) {

  int type, ti, ts[] = { 6, 10, 12 }, n = 5000, ok;
  long long L;
  BoxStats st;
  char name[120], msg[200];

  for (type = 0; type < 3; ++type)
    for (ti = 0; ti < 3; ++ti) {
      makeInstance (n, type);
      tourAlloc();
      boxTour (ts[ti], &st);
      L = cycleLength (tour, n);                   /* recomputed independently */
      ok = isPermutation (tour, n) && L <= st.loopsSum + st.glueG;
      sprintf (name, "box tour is a valid tour, length <= loops + G (%s, t=%d)",
	       typeName[type], ts[ti]);
      sprintf (msg, "length %lld, loops %lld, G %lld, %d capped boxes",
	       L, st.loopsSum, st.glueG, st.capped);
      report (name, ok, msg);
      free (tour); free (pos);
    }
}

/* ---------------------------------------------------------------- 5. 2-opt and Or-opt */

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

/* Brute force over every Or-opt move (chains of 1-3 cities, every place,
   both orientations).  inScope counts the ones our search is designed to
   consider: a new edge from a chain end that is shorter than the saving g0
   from cutting the chain out (the gain criterion in localsearch.c). */
static void improvingOrOpt (long long *all, long long *inScope) {
  int i, len, dir, k, j, s1, s2, a, b, u, v, g0, inSeg, x;
  long long dSame, dRev, best, near;
  int N = ncities, seg[3];
  *all = *inScope = 0;
  for (i = 0; i < N; ++i)
    for (len = 1; len <= 3; ++len)
      for (dir = 1; dir >= -1; dir -= 2) {
	for (k = 0; k < len; ++k) seg[k] = tour[((i + dir * k) % N + N) % N];
	s1 = seg[0]; s2 = seg[len - 1];
	a = tour[((i - dir) % N + N) % N];
	b = tour[((i + dir * len) % N + N) % N];
	g0 = distance (a,s1) + distance (s2,b) - distance (a,b);
	if (g0 <= 0) continue;
	for (j = 0; j < N; ++j) {
	  u = tour[j]; v = tour[((j + dir) % N + N) % N];
	  inSeg = 0;
	  for (k = 0; k < len; ++k) if (u == seg[k] || v == seg[k]) inSeg = 1;
	  if (inSeg) continue;
	  dSame = g0 + distance (u,v) - distance (u,s1) - distance (s2,v);
	  dRev  = g0 + distance (u,v) - distance (u,s2) - distance (s1,v);
	  best = dSame > dRev ? dSame : dRev;
	  if (best <= 0) continue;
	  (*all)++;
	  near = distance (s1,u);
	  x = distance (s1,v); if (x < near) near = x;
	  if (len > 1) {
	    x = distance (s2,u); if (x < near) near = x;
	    x = distance (s2,v); if (x < near) near = x;
	  }
	  if (near < g0) (*inScope)++;
	}
      }
}

static void testLocalSearch (void) {

  int type, n = 300, rounds;
  long long L0, L1, L2, m2, mo, b2, orAll, orIn;
  BoxStats st;
  char name[120], msg[200];

  for (type = 0; type < 3; ++type) {
    makeInstance (n, type);
    tourAlloc();
    boxTour (8, &st);
    buildNeighbors (n - 1);                 /* complete lists: nothing hidden */
    L0 = cycleLength (tour, n);

    localSearchToOptimum (0, &m2, &mo, &rounds);       /* 2-opt */
    L1 = cycleLength (tour, n);
    b2 = improving2opt();
    sprintf (name, "2-opt: never longer, brute force finds no improving move (%s)",
	     typeName[type]);
    sprintf (msg, "%lld -> %lld in %lld moves, %d rounds; %lld improving moves left",
	     L0, L1, m2, rounds, b2);
    report (name, isPermutation (tour, n) && L1 <= L0 && b2 == 0, msg);

    localSearchToOptimum (1, &m2, &mo, &rounds);       /* 2-opt + Or-opt */
    L2 = cycleLength (tour, n);
    b2 = improving2opt();
    improvingOrOpt (&orAll, &orIn);
    sprintf (name, "Or-opt: never longer, no 2-opt or in-scope Or-opt move left (%s)",
	     typeName[type]);
    sprintf (msg, "%lld -> %lld (%lld Or-opt moves); left: %lld 2-opt, %lld Or-opt in "
	     "scope, %lld Or-opt outside the gain criterion", L1, L2, mo, b2, orIn, orAll);
    report (name, isPermutation (tour, n) && L2 <= L1 && b2 == 0 && orIn == 0, msg);

    free (tour); free (pos); free (nbList);
  }
}

/* ---------------------------------------------------------------- 6. whole pipeline */

static void testPipeline (void) {

  int type, n = 20000, rounds;
  long long L, m2, mo;
  BoxStats st;
  char name[120], msg[200];

  for (type = 0; type < 3; ++type) {
    makeInstance (n, type);
    tourAlloc();
    boxTour (10, &st);
    buildNeighbors (16);
    localSearchToOptimum (0, &m2, &mo, &rounds);
    localSearchToOptimum (1, &m2, &mo, &rounds);
    L = cycleLength (tour, n);
    sprintf (name, "full method on 20,000 cities finishes with a valid tour (%s)",
	     typeName[type]);
    sprintf (msg, "length %lld, %d boxes, %d capped, largest box %d",
	     L, st.leaves, st.capped, st.largest);
    report (name, isPermutation (tour, n) && L == (long long) tourLength(), msg);
    free (tour); free (pos); free (nbList);
  }
}

/* ================================================================ curves (2D extension) */

/* Square-cell curves on a small grid: every cell a different key, and
   consecutive keys are neighbouring cells (Morton: measured, not required). */
static void testGridCurves (void) {

  int curves[] = { CURVE_HILBERT, CURVE_MOORE, CURVE_PEANO, CURVE_MORTON, CURVE_STRIP };
  int ci, c, levels, ok, touching;
  unsigned long long side, n, k, *cx, *cy;
  unsigned x, y;
  char *seen, name[120], msg[200];

  for (ci = 0; ci < 5; ++ci) {
    c = curves[ci];
    levels = (c == CURVE_PEANO) ? 3 : 4;                 /* 27x27 or 16x16 */
    side = curveGridSide (c, levels);
    n = side * side;
    cx = (unsigned long long *) malloc (n * sizeof *cx);
    cy = (unsigned long long *) malloc (n * sizeof *cy);
    seen = (char *) calloc (n, 1);
    ok = 1;
    for (x = 0; x < side; ++x)                           /* cell at each key */
      for (y = 0; y < side; ++y) {
	k = curveGridKey (c, x, y, levels);
	if (k >= n || seen[k]) ok = 0;
	else { seen[k] = 1; cx[k] = x; cy[k] = y; }
      }
    sprintf (name, "%s: every cell of a %llux%llu grid gets its own position", curveInfo[c].name, side, side);
    report (name, ok, "");

    touching = 0;                                        /* consecutive cells touch? */
    if (ok)
      for (k = 1; k < n; ++k)
	if (llabs ((long long) cx[k] - (long long) cx[k-1]) + llabs ((long long) cy[k] - (long long) cy[k-1]) == 1)
	  touching++;
    sprintf (name, "%s: consecutive cells are neighbours", curveInfo[c].name);
    sprintf (msg, "%d of %llu steps (%.1f%%)", touching, n - 1, 100.0 * touching / (n - 1));
    if (c == CURVE_MORTON)                               /* expected to jump */
      printf ("NOTE  %s - %s; expected for Morton (negative control)\n", name, msg);
    else report (name, ok && (unsigned long long) touching == n - 1, msg);

    if (curveInfo[c].closed) {
      int closes = llabs ((long long) cx[n-1] - (long long) cx[0]) + llabs ((long long) cy[n-1] - (long long) cy[0]) == 1;
      sprintf (name, "%s: closed curve, last cell next to the first", curveInfo[c].name);
      report (name, closes, "");
    }
    free (cx); free (cy); free (seen);
  }
}

static int sharedCorners (double a[3][2], double b[3][2]) {
  int i, j, s = 0;
  for (i = 0; i < 3; ++i)
    for (j = 0; j < 3; ++j)
      if (a[i][0] == b[j][0] && a[i][1] == b[j][1]) s++;
  return s;
}

/* Triangle curves: the triangle with key k contains its own centre point
   (key of the centre = k), and consecutive triangles share a whole edge. */
static void testTriangleCurves (void) {

  int curves[] = { CURVE_SIERPINSKI, CURVE_KNOPP }, ci, c, levels = 10, ok, edges, shift;
  unsigned long long n = 1ULL << levels, k;
  double prev[3][2], cur[3][2], first[3][2], cxm, cym;
  char name[120], msg[200];

  for (ci = 0; ci < 2; ++ci) {
    c = curves[ci];
    shift = curveInfo[c].digits - levels;
    ok = 1;
    edges = 0;
    for (k = 0; k < n; ++k) {
      curveTriangle (c, k, levels, cur);
      cxm = (cur[0][0] + cur[1][0] + cur[2][0]) / 3;     /* centre of the triangle */
      cym = (cur[0][1] + cur[1][1] + cur[2][1]) / 3;
      if ((curveKey (c, cxm, cym) >> shift) != k) ok = 0;
      if (k == 0) memcpy (first, cur, sizeof cur);
      else if (sharedCorners (prev, cur) == 2) edges++;
      memcpy (prev, cur, sizeof cur);
    }
    sprintf (name, "%s: each of %llu triangles maps back to its own position", curveInfo[c].name, n);
    report (name, ok, "");
    sprintf (name, "%s: consecutive triangles share an edge", curveInfo[c].name);
    sprintf (msg, "%d of %llu steps", edges, n - 1);
    report (name, (unsigned long long) edges == n - 1, msg);
    if (curveInfo[c].closed) {
      sprintf (name, "%s: closed curve, last triangle shares an edge with the first", curveInfo[c].name);
      report (name, sharedCorners (prev, first) == 2, "");
    }
  }
}

/* Boxes for every curve: each box is exactly one curve cell (all its keys
   share the cell's prefix, no other city falls in that cell), size <= t
   unless the depth cap was hit. */
static void testCurveBoxes (void) {

  int c, type, n = 5000, t = 10, nb, j, k, ok, capped;
  int *order, *start, *level;
  unsigned long long *keys, size, cell, pw;
  char name[160], msg[160];

  for (c = 0; c < NCURVES; ++c)
    for (type = 0; type < 3; ++type) {
      makeInstance (n, type);
      order = (int *) malloc (n * sizeof(int));
      keys = (unsigned long long *) malloc (n * sizeof *keys);
      start = (int *) malloc ((n + 1) * sizeof(int));
      level = (int *) malloc ((n + 1) * sizeof(int));
      curveSort (c, order, keys);
      nb = boxSplitBase (keys, n, curveInfo[c].base, curveInfo[c].digits, t, start, level);
      ok = isPermutation (order, n) && start[0] == 0 && start[nb] == n;
      for (k = 1; k < n; ++k) if (keys[k] < keys[k-1]) ok = 0;
      capped = 0;
      for (j = 0; j < nb && ok; ++j) {
	pw = 1;
	for (k = 0; k < curveInfo[c].digits - level[j]; ++k) pw *= curveInfo[c].base;
	size = pw;                                    /* keys in this box's cell */
	cell = keys[start[j]] / size;                 /* the cell's prefix       */
	for (k = start[j]; k < start[j + 1]; ++k)
	  if (keys[k] / size != cell) ok = 0;         /* outside the cell        */
	if (start[j] > 0 && keys[start[j] - 1] / size == cell) ok = 0;      /* another */
	if (start[j + 1] < n && keys[start[j + 1]] / size == cell) ok = 0;  /* city in it */
	if (start[j + 1] - start[j] > t) {
	  capped++;
	  if (level[j] != curveInfo[c].digits) ok = 0;   /* only at the depth cap */
	}
      }
      sprintf (name, "%s boxes are whole curve cells, size <= t (%s)", curveInfo[c].name, typeName[type]);
      sprintf (msg, "%d boxes, %d at the depth cap", nb, capped);
      report (name, ok, msg);
      free (order); free (keys); free (start); free (level);
    }
}

/* ================================================================ merges (2D extension) */

/* Two unit squares two apart: the best merge gives the optimal tour, 10.
   Coordinates are x1000 because the bench rounds distances to integers. */
static void testTwoSquares (void) {

  double px[] = { 0, 1000, 1000, 0,   3000, 4000, 4000, 3000 };
  double py[] = { 0, 0,    1000, 1000, 0,   0,    1000, 1000 };
  int order[8] = { 1, 2, 3, 4, 5, 6, 7, 8 }, start[3] = { 0, 4, 8 };
  int loops[8] = { 1, 2, 3, 4, 5, 6, 7, 8 }, m, c;
  BoxSet bs;
  BoxStats st;
  long long added, L;
  char name[120], msg[160];

  makeInstance (8, UNIFORM);
  for (c = 1; c <= 8; ++c) { xx[c] = px[c-1]; yy[c] = py[c-1]; }
  bs.nb = 2; bs.order = order; bs.start = start; bs.loops = loops;
  buildNeighbors (7);
  for (m = 0; m < NMERGES; ++m) {
    tourAlloc();
    st.loopsSum = 8000;                               /* two loops of 4000 */
    mergeLoopsAdded (m, &bs, &st, &added);
    L = cycleLength (tour, 8);
    sprintf (name, "%s on two unit squares two apart", mergeName[m]);
    sprintf (msg, "length %.0f (optimum 10), reported added cost %.0f, recomputed %.0f",
	     L / 1000.0, added / 1000.0, (L - 8000) / 1000.0);
    if (m == MERGE_M1)                                /* first cut is fixed */
      report (name, isPermutation (tour, 8) && added == L - 8000, msg);
    else
      report (name, isPermutation (tour, 8) && L == 10000 && added == L - 8000, msg);
    free (tour); free (pos);
  }
  free (nbList);
}

/* Every merge on every data type: one valid tour, and the extra length the
   method says it added equals the recomputed tour length minus the loops. */
static void testMergesValid (void) {

  int curves[] = { CURVE_HILBERT, CURVE_SIERPINSKI, CURVE_PEANO };
  int ci, type, m, n = 5000;
  long long added, L, m1 = 0;
  BoxSet bs;
  BoxStats st;
  char name[160], msg[200];

  for (ci = 0; ci < 3; ++ci)
    for (type = 0; type < 3; ++type) {
      makeInstance (n, type);
      buildNeighbors (16);
      boxBuild (curves[ci], 10, &bs, &st);
      for (m = 0; m < NMERGES; ++m) {
	tourAlloc();
	mergeLoopsAdded (m, &bs, &st, &added);
	L = cycleLength (tour, n);
	if (m == MERGE_M1) m1 = L;
	sprintf (name, "%s + %s (%s): valid tour, added cost = recomputed",
		 curveInfo[curves[ci]].name, mergeName[m], typeName[type]);
	sprintf (msg, "length %lld, added %lld, recomputed %lld", L, added, L - st.loopsSum);
	report (name, isPermutation (tour, n) && added == L - st.loopsSum
		&& L == (long long) tourLength() && (m != MERGE_M2 || L <= m1), msg);
	free (tour); free (pos);
      }
      boxFree (&bs);
      free (nbList);
    }
}

/* M2 is the best possible snake: compare with every combination of cuts
   and directions on small instances (brute force). */
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
  for (trial = 0; trial < 40; ++trial) {
    n = 12 + trial % 7;                                  /* 12-18 cities */
    makeInstance (n, trial % 2 ? CLUSTERED : UNIFORM);
    boxBuild (CURVE_HILBERT, 4, &bs, &st);
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

int main (void) {

  testHilbert();
  testBoxes();
  testDP();
  testGlue();
  testLocalSearch();
  testPipeline();
  testGridCurves();
  testTriangleCurves();
  testCurveBoxes();
  testTwoSquares();
  testMergesValid();
  testM2Optimal();
  printf ("\n%s: %d failure(s)\n", failures ? "FAILED" : "ALL PASSED", failures);
  return failures ? 1 : 0;
}
