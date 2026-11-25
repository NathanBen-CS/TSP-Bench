#include <math.h>
#include <stdio.h>
#include <stdlib.h>
#include <string.h>

#include "Instance.h"
#include "Utilities.h"
#include "tour.h"
#include "neighbors.h"
#include "boxes.h"
#include "merges.h"
#include "rings2d.h"

/*********************************************************************

Tests for the hull rings (solvers/d2/rings2d.c) and merge H0.

	ringtest2d         run every test, print PASS/FAIL per test,
	                   exit code 0 only if all pass

Instances are made in memory with the bench's rounding (nearest integer,
as EUC_2D), like tests/boxtest2d.c.

 *********************************************************************/

static int failures = 0;

static void report (const char *name, int ok, const char *detail) {
  printf ("%-4s  %s%s%s\n", ok ? "PASS" : "FAIL", name, detail[0] ? " - " : "", detail);
  if (!ok) failures++;
}

static unsigned long long rs = 88172645463325252ULL;          /* xorshift64 */
static unsigned long long rnd (void) { rs ^= rs << 13; rs ^= rs >> 7; rs ^= rs << 17; return rs; }
static double urand (void) { return (rnd() >> 11) * (1.0 / 9007199254740992.0); }
static double nrand (void) {
  return sqrt (-2 * log (urand() + 1e-300)) * cos (6.283185307179586 * urand());
}

static int roundDist (int a, int b) {        /* same rule as EUC_2D in Instance.c */
  double dx = xx[a] - xx[b], dy = yy[a] - yy[b];
  return (int) (sqrt (dx*dx + dy*dy) + 0.5);
}

static void newInstance (int n) {
  free (xx); free (yy);
  ncities = n;
  xx = (double *) calloc (n + 1, sizeof(double));
  yy = (double *) calloc (n + 1, sizeof(double));
  isEuclidean = 1;
  distance = roundDist;
  tourAlloc();
}

enum { UNIFORM, CLUSTERED, DUPLICATES, CIRCLE };

static void makeInstance (int n, int type) {
  int c, k, nc;
  double cx[1000], cy[1000], sd, a;
  newInstance (n);
  if (type == UNIFORM)
    for (c = 1; c <= n; ++c) { xx[c] = rnd() % 1000000; yy[c] = rnd() % 1000000; }
  else if (type == CLUSTERED) {                /* like portcgen: n/100 centres */
    nc = n / 100 > 0 ? (n / 100 < 1000 ? n / 100 : 1000) : 1;
    for (k = 0; k < nc; ++k) { cx[k] = urand() * 1e6; cy[k] = urand() * 1e6; }
    sd = 1e6 / sqrt ((double) n);
    for (c = 1; c <= n; ++c) {
      k = rnd() % nc;
      xx[c] = (int) (cx[k] + sd * nrand());
      yy[c] = (int) (cy[k] + sd * nrand());
    }
  }
  else if (type == DUPLICATES) {               /* ~25 cities per spot */
    nc = n / 25 > 1 ? n / 25 : 2;
    for (c = 1; c <= n; ++c) {
      k = rnd() % nc;
      xx[c] = (k * 7919) % 1000;
      yy[c] = (k * 104729) % 1000;
    }
  }
  else                                         /* convex position: on a circle, one  */
    for (c = 1; c <= n; ++c) {                 /*   city per sector of 360/n degrees */
      a = 6.283185307179586 * (c - 1 + 0.8 * urand()) / n;   /* (fully random angles can put  */
      /* two cities so close that rounding pulls one inside: then they are not in convex position) */
      xx[c] = floor (500000 + 400000 * cos (a) + 0.5);
      yy[c] = floor (500000 + 400000 * sin (a) + 0.5);
    }
}

static int isPermutation (const int *a, int n) {
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

static void ringSizes (const BoxSet *bs, char *buf) {
  int j, o = 0;
  buf[0] = 0;
  for (j = 0; j < bs->nb && o < 200; ++j)
    o += sprintf (buf + o, "%s%d", j ? "," : "", bs->start[j + 1] - bs->start[j]);
}

/* ---------------------------------------------------------------- 1. convex position */

static void testConvexOptimal (void) {
  int trial, n, ok = 1, oneRing = 1, loop[BOX_MAXT];
  long long hk, ring;
  BoxSet bs; BoxStats st; RingStats r;
  char msg[160] = "";
  for (trial = 0; trial < 300; ++trial) {
    n = 3 + trial % 10;                          /* 3 .. 12 cities on a circle */
    makeInstance (n, CIRCLE);
    ringBuild (&bs, &st, &r);
    if (bs.nb != 1) { oneRing = 0; sprintf (msg, "trial %d: %d rings", trial, bs.nb); }
    ring = cycleLength (bs.loops, n);
    hk = heldKarpLoop (bs.order, n, loop);
    if (ring != hk || ring != st.loopsSum) {
      ok = 0;
      sprintf (msg, "trial %d, n=%d: ring %lld, Held-Karp %lld", trial, n, ring, hk);
    }
    boxFree (&bs);
  }
  report ("convex position: one ring of every city (300 sets, 3-12 cities)", oneRing, oneRing ? "" : msg);
  report ("convex position: the ring is exactly as long as the Held-Karp optimum", ok, ok ? "" : msg);
}

/* ---------------------------------------------------------------- 2. small random sets */

static void testSmallRandom (void) {
  int trial, n, j, lo, s, m, okRing = 1, okMerge = 1, okHull = 1, okCheck = 1, loop[BOX_MAXT];
  long long added;
  BoxSet bs; BoxStats st; RingStats r;
  char msg[200] = "", msg2[200] = "";
  for (trial = 0; trial < 300; ++trial) {
    n = 4 + trial % 9;                           /* 4 .. 12 cities */
    makeInstance (n, trial % 3 == 2 ? DUPLICATES : UNIFORM);
    buildNeighbors (n - 1 < 16 ? n - 1 : 16);
    ringBuild (&bs, &st, &r);
    if (ringCheck (&bs, 1)) okCheck = 0;
    for (j = 0; j < bs.nb; ++j) {                /* every ring: its hull order is optimal */
      lo = bs.start[j];
      s = bs.start[j + 1] - lo;
      if (s >= 4 && cycleLength (bs.loops + lo, s) != heldKarpLoop (bs.order + lo, s, loop)) {
	okRing = 0;
	sprintf (msg, "trial %d ring %d (%d cities)", trial, j, s);
      }
    }
    for (m = 0; m <= MERGE_H0; ++m) {           /* every merge: a valid tour of known length */
      mergeLoopsAdded (m, &bs, &st, &added);
      if (!isPermutation (tour, n) || cycleLength (tour, n) != st.loopsSum + added
	  || st.joinCost != added) {
	okMerge = 0;
	sprintf (msg2, "trial %d, %s", trial, mergeName[m]);
      }
    }
    hullOrderTour (&bs);
    if (!isPermutation (tour, n)) okHull = 0;
    boxFree (&bs);
  }
  report ("small random sets (300, 4-12 cities, some duplicates): ring checks pass", okCheck, "");
  report ("small random sets: every ring of 4+ cities is an optimal loop of its cities", okRing, okRing ? "" : msg);
  report ("small random sets: M1-M5 and H0 give permutations, length = rings + merge cost", okMerge, okMerge ? "" : msg2);
  report ("small random sets: hull order alone is a permutation", okHull, "");
}

/* ---------------------------------------------------------------- 3. grids and lines */

static void testGrid (int side, const char *want) {
  int x, y, c = 1;
  BoxSet bs; BoxStats st; RingStats r;
  char got[220], name[80];
  newInstance (side * side);
  for (x = 0; x < side; ++x)
    for (y = 0; y < side; ++y) { xx[c] = x * 1000; yy[c] = y * 1000; c++; }
  ringBuild (&bs, &st, &r);
  ringSizes (&bs, got);
  sprintf (name, "%dx%d grid: collinear boundary points stay in their ring", side, side);
  report (name, !strcmp (got, want) && ringCheck (&bs, 1) == 0, got);
  boxFree (&bs);
}

static void testLine (void) {
  int c, n = 12, lo = 1, hi = 1;
  long long want;
  BoxSet bs; BoxStats st; RingStats r;
  char msg[120];
  newInstance (n);                             /* 12 cities on 11 spots of one line */
  for (c = 1; c <= n; ++c) { xx[c] = (c * 37) % 11 * 3000; yy[c] = (c * 37) % 11 * 4000; }
  for (c = 1; c <= n; ++c) { if (xx[c] < xx[lo]) lo = c; if (xx[c] > xx[hi]) hi = c; }
  want = 2 * (long long) distance (lo, hi);    /* out to the far end and back */
  ringBuild (&bs, &st, &r);
  sprintf (msg, "%d ring(s), loop %lld, 2 x span %lld", bs.nb, cycleLength (bs.loops, n), want);
  report ("points on one line (with duplicates): one ring, out and back", bs.nb == 1
	  && isPermutation (bs.loops, n) && cycleLength (bs.loops, n) == want, msg);
  boxFree (&bs);
}

static void testTiny (void) {
  int n, ok = 1, m;
  long long added;
  BoxSet bs; BoxStats st; RingStats r;
  for (n = 1; n <= 3; ++n) {
    makeInstance (n, UNIFORM);
    if (n == 3) { xx[3] = (xx[1] + xx[2]) / 2; yy[3] = (yy[1] + yy[2]) / 2; xx[1] = 2 * xx[3] - xx[2]; yy[1] = 2 * yy[3] - yy[2]; }
    if (n > 1) buildNeighbors (n - 1);
    ringBuild (&bs, &st, &r);
    if (bs.nb != 1 || ringCheck (&bs, 1)) ok = 0;
    for (m = 0; m <= MERGE_H0 && n > 1; ++m) {
      mergeLoopsAdded (m, &bs, &st, &added);
      if (!isPermutation (tour, n)) ok = 0;
    }
    boxFree (&bs);
  }
  report ("1, 2 and 3 cities (3 on a line): one ring, valid tours", ok, "");
}

/* ---------------------------------------------------------------- 4. larger maps */

static void testLarge (int n, int type, const char *label) {
  int m, okMerge = 1, f;
  long long added;
  BoxSet bs; BoxStats st; RingStats r;
  char name[160], msg[160];
  makeInstance (n, type);
  buildNeighbors (16);
  ringBuild (&bs, &st, &r);
  f = ringCheck (&bs, 1);
  sprintf (msg, "%d rings, %.2f x N^(2/3), largest %d", r.rings, r.rings / pow (n, 2.0 / 3.0), r.largest);
  sprintf (name, "%s %d: rings convex, later rings strictly inside, every city once", label, n);
  report (name, f == 0, msg);
  for (m = 0; m <= MERGE_H0; ++m) {
    mergeLoopsAdded (m, &bs, &st, &added);
    if (!isPermutation (tour, n) || cycleLength (tour, n) != st.loopsSum + added) okMerge = 0;
  }
  hullOrderTour (&bs);
  if (!isPermutation (tour, n)) okMerge = 0;
  sprintf (name, "%s %d: every merge and hull order give valid tours of the right length", label, n);
  report (name, okMerge, "");
  boxFree (&bs);
}

int main (void) {
  testConvexOptimal();
  testSmallRandom();
  testGrid (10, "36,28,20,12,4");
  testGrid (9, "32,24,16,8,1");
  testLine();
  testTiny();
  testLarge (5000, UNIFORM, "uniform");
  testLarge (5000, CLUSTERED, "clustered");
  testLarge (3000, DUPLICATES, "duplicate-heavy");
  printf ("%s: %d failure(s)\n", failures ? "FAILED" : "ALL PASS", failures);
  return failures ? 1 : 0;
}
