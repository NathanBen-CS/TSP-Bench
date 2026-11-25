#include "testlib3d.h"
#include "neighbors.h"
#include "curves3d.h"
#include "boxes.h"

/*********************************************************************

3D unit tests: reader, neighbour lists, curves.  (Boxes, merges and cleanup
are in boxtest3d.c.)

	unit3d             run every test, print PASS/FAIL per test,
	                   exit code 0 only if all pass

 *********************************************************************/

/* ---------------------------------------------------------------- 1. reader */

static void testReader (void) {

  const char *file = "unit3d_tmp.tsp";
  FILE *fp = fopen (file, "wb");
  int ok;
  char msg[200];

  /* header in star-file style, real and negative coordinates, cities out of
     order, CRLF line ends, a blank line, no EOF line */
  fprintf (fp, "NAME: odd\r\nTYPE: TSP\r\nCOMMENT : NAME inside a comment\r\nCOMMENT\r\n"
	   "DIMENSION: 4\r\nEDGE_WEIGHT_TYPE: EUC_3D\r\nNODE_COORD_SECTION\r\n"
	   "2 1.5 -2.25 1e2\r\n1 0.000000 0 0\r\n\r\n4 3 4 12\r\n3 0.3 0.4 0\r\n");
  fclose (fp);
  readInstance ((char *) file);
  remove (file);
  ok = ncities == 4 && xx[2] == 1.5 && yy[2] == -2.25 && zz[2] == 100
    && xx[3] == 0.3 && zz[4] == 12 && isEuclidean && notIntegral && !strcmp (instanceName, "odd");
  report ("reader: star-style header, real coordinates, any city order, no EOF line", ok, "");

  /* (3,4,12) is 13 from the origin; (0.3,0.4,0) is 0.5, which rounds up */
  sprintf (msg, "d(1,4)=%d (13), d(1,3)=%d (0.5 -> 1), d(3,4)=%d (12.53 -> 13)",
	   distance (1,4), distance (1,3), distance (3,4));
  report ("reader: nearest-integer 3D distance", distance (1,4) == 13 && distance (1,3) == 1
	  && distance (3,4) == 13 && distance (4,3) == 13 && distance (2,2) == 0, msg);
}

/* ---------------------------------------------------------------- 2. neighbours */

static double trueDist2 (int a, int b) {
  double dx = xx[a]-xx[b], dy = yy[a]-yy[b], dz = zz[a]-zz[b];
  return dx*dx + dy*dy + dz*dz;
}

static int dcmp (const void *a, const void *b) {
  double x = *(const double *) a, y = *(const double *) b;
  return x < y ? -1 : x > y;
}

/* Is list[0..k-1] a correct k-nearest list of city c among set[0..m-1]?
   Ties may come in any order, so the sorted distances are compared; every
   entry must be in the set, not c itself, and listed once. */
static int listOk (int c, const int *list, int k, const int *set, int m, const char *inSet) {
  int i, j, n = 0, ok = 1;
  double *d = (double *) malloc (m * sizeof(double));
  for (i = 0; i < m; ++i) if (set[i] != c) d[n++] = trueDist2 (c, set[i]);
  qsort (d, n, sizeof(double), dcmp);
  for (i = 0; i < k; ++i) {
    if (list[i] < 1 || list[i] > ncities || list[i] == c || !inSet[list[i]]) { ok = 0; break; }
    if (trueDist2 (c, list[i]) != d[i]) ok = 0;
    for (j = 0; j < i; ++j) if (list[j] == list[i]) ok = 0;
  }
  free (d);
  return ok;
}

static void testNeighbors (void) {

  int types[] = { UNIFORM, CLUSTERED, DUPLICATES }, ti, ki, ks[] = { 1, 8, 16, 32 };
  int n = 300, c, i, m, bad, *set, *out;
  char *inSet, name[160], msg[120];

  for (ti = 0; ti < 5; ++ti)
    for (ki = 0; ki < 4; ++ki) {
      if (ti < 3) makeInstance (n, types[ti]);
      else makeInstance (n, UNIFORM);
      if (ti == 3) for (c = 1; c <= n; ++c) zz[c] = 5;                /* flat: one plane */
      if (ti == 4) for (c = 1; c <= n; ++c) if (c % 10) {             /* tiny core, far rim */
	xx[c] = 500000 + c % 7; yy[c] = 500000 + c % 11; zz[c] = 500000 + c % 5; }
      set = (int *) malloc (n * sizeof(int));
      inSet = (char *) calloc (n + 1, 1);
      for (c = 1; c <= n; ++c) { set[c-1] = c; inSet[c] = 1; }
      buildNeighbors (ks[ki]);
      for (bad = 0, c = 1; c <= n; ++c) if (!listOk (c, neighbors (c), nbK, set, n, inSet)) bad++;

      /* nearestAmong on every third city */
      memset (inSet, 0, n + 1);
      for (m = 0, c = 1; c <= n; c += 3) { set[m++] = c; inSet[c] = 1; }
      out = (int *) malloc (m * ks[ki] * sizeof(int));
      nearestAmong (set, m, ks[ki], out);
      for (i = 0; i < m; ++i) if (!listOk (set[i], out + i * ks[ki], ks[ki], set, m, inSet)) bad++;
      sprintf (name, "neighbours = brute force, 300 cities (%s), K=%d, and a subset of 100",
	       ti < 3 ? typeName[types[ti]] : ti == 3 ? "one plane" : "dense core + far points", ks[ki]);
      sprintf (msg, "%d wrong lists", bad);
      report (name, bad == 0, msg);
      free (set); free (inSet); free (out); free (nbList);
    }

  makeInstance (5, UNIFORM);                        /* fewer cities than K */
  buildNeighbors (16);
  report ("neighbours: K is cut to N-1 on tiny instances", nbK == 4, "");
  free (nbList);

  makeInstance (40, UNIFORM);                       /* 40 cities on one spot */
  for (c = 1; c <= 40; ++c) xx[c] = yy[c] = zz[c] = 7;
  buildNeighbors (16);
  for (bad = 0, c = 1; c <= 40; ++c)
    for (i = 0; i < nbK; ++i) if (neighbors (c)[i] == c) bad++;
  report ("neighbours: all cities on one spot", bad == 0, "");
  free (nbList);
}

/* ---------------------------------------------------------------- 3. curves */

static void testCubeCurves (void) {

  int curves[] = { CURVE_HILBERT3, CURVE_MORTON3, CURVE_PEANO3, CURVE_SNAKE3, CURVE_MOORE3 };
  int ci, c, L, levels, i, ok, closes;
  unsigned long long side, n, k, touching, sub, cells, lo, hi, badCells, *minKey, *maxKey, *cnt;
  unsigned x, y, z, *cx, *cy, *cz;
  char *seen, name[160], msg[160];

  for (ci = 0; ci < 5; ++ci)
    for (levels = 1; levels <= (curves[ci] == CURVE_PEANO3 ? 3 : 5); ++levels) {
      c = curves[ci];
      side = curveGridSide3 (c, levels);
      n = side * side * side;
      cx = (unsigned *) malloc (n * sizeof *cx);
      cy = (unsigned *) malloc (n * sizeof *cy);
      cz = (unsigned *) malloc (n * sizeof *cz);
      seen = (char *) calloc (n, 1);

      ok = 1;                                          /* every cell its own key */
      for (x = 0; x < side; ++x)
	for (y = 0; y < side; ++y)
	  for (z = 0; z < side; ++z) {
	    k = curveGridKey3 (c, x, y, z, levels);
	    if (k >= n || seen[k]) ok = 0;
	    else { seen[k] = 1; cx[k] = x; cy[k] = y; cz[k] = z; }
	  }
      sprintf (name, "%s: bijection on a %llu^3 grid", curveInfo[c].name, side);
      report (name, ok, "");
      if (!ok) { free (cx); free (cy); free (cz); free (seen); continue; }

      touching = 0;                                    /* steps to a face neighbour */
      for (k = 1; k < n; ++k)
	if (abs ((int) cx[k] - (int) cx[k-1]) + abs ((int) cy[k] - (int) cy[k-1])
	    + abs ((int) cz[k] - (int) cz[k-1]) == 1) touching++;
      sprintf (name, "%s: consecutive cells share a face (%llu^3)", curveInfo[c].name, side);
      sprintf (msg, "%llu of %llu steps (%.1f%%)", touching, n - 1, 100.0 * touching / (n - 1));
      if (c == CURVE_MORTON3)                          /* expected to jump */
	printf ("NOTE  %s - %s; expected for Morton (negative control)\n", name, msg);
      else report (name, touching == n - 1, msg);

      closes = abs ((int) cx[n-1] - (int) cx[0]) + abs ((int) cy[n-1] - (int) cy[0])
	+ abs ((int) cz[n-1] - (int) cz[0]) == 1;
      if (curveInfo[c].closed) {
	sprintf (name, "%s: closed, last cell shares a face with the first (%llu^3)", curveInfo[c].name, side);
	report (name, closes, "");
      }

      /* at every level the cells are contiguous, aligned key ranges (the
	 property the box split relies on).  Snake3's base-8 ranges are not
	 cubes, so for it the check is on the key ranges only, which any
	 bijection passes; it is skipped. */
      if (c != CURVE_SNAKE3) {
	badCells = 0;
	for (L = 1; L < levels; ++L) {
	  sub = curveGridSide3 (c, levels - L);        /* side of a cell at level L */
	  cells = n / (sub * sub * sub);
	  minKey = (unsigned long long *) malloc (cells * sizeof *minKey);
	  maxKey = (unsigned long long *) calloc (cells, sizeof *maxKey);
	  cnt = (unsigned long long *) calloc (cells, sizeof *cnt);
	  for (k = 0; k < cells; ++k) minKey[k] = ~0ULL;
	  for (k = 0; k < n; ++k) {
	    i = (int) ((cx[k] / sub) * (side / sub) * (side / sub) + (cy[k] / sub) * (side / sub) + cz[k] / sub);
	    if (k < minKey[i]) minKey[i] = k;
	    if (k > maxKey[i]) maxKey[i] = k;
	    cnt[i]++;
	  }
	  for (k = 0; k < cells; ++k) {
	    lo = minKey[k];
	    hi = maxKey[k];
	    if (hi - lo + 1 != cnt[k] || lo % cnt[k] != 0) badCells++;
	  }
	  free (minKey); free (maxKey); free (cnt);
	}
	if (levels > 1) {
	  sprintf (name, "%s: cells at every level are contiguous key ranges (%llu^3)", curveInfo[c].name, side);
	  sprintf (msg, "%llu bad cells", badCells);
	  report (name, badCells == 0, msg);
	}
      }
      free (cx); free (cy); free (cz); free (seen);
    }

  /* Morton3 against the definition; Hilbert3's first steps */
  ok = 1;
  for (x = 0; x < 16; ++x) for (y = 0; y < 16; ++y) for (z = 0; z < 16; ++z) {
    k = 0;
    for (i = 3; i >= 0; --i) k = 8*k + 4*((z>>i)&1) + 2*((y>>i)&1) + ((x>>i)&1);
    if (curveGridKey3 (CURVE_MORTON3, x, y, z, 4) != k) ok = 0;
  }
  report ("morton3: bits of z, y, x interleaved (z most significant)", ok, "");
  sprintf (msg, "(0,0,0) -> %llu, (7,0,0) -> %llu", curveGridKey3 (CURVE_HILBERT3,0,0,0,3),
	   curveGridKey3 (CURVE_HILBERT3,7,0,0,3));
  report ("hilbert3: starts at (0,0,0), ends at (n-1,0,0)", curveGridKey3 (CURVE_HILBERT3,0,0,0,3) == 0
	  && curveGridKey3 (CURVE_HILBERT3,7,0,0,3) == 511, msg);
}

/* Full resolution: keys sorted, order a permutation, key = the key of the
   city's own cell, and the largest key fits the curve's digits. */
static void testCurveSort (void) {

  int c, type, n = 5000, k, ok, *order;
  unsigned long long *keys, top;
  double span;
  char name[160], msg[120];

  for (c = 0; c < NCURVES3; ++c)
    for (type = 0; type < NTYPES; ++type) {
      makeInstance (n, type);
      order = (int *) malloc (n * sizeof(int));
      keys = (unsigned long long *) malloc (n * sizeof *keys);
      curveSort (c, order, keys);
      span = pow (curveInfo[c].base, curveInfo[c].digits);
      ok = isPermutation (order, n);
      top = 0;
      curveFrame();
      for (k = 0; k < n; ++k) {
	if (k && keys[k] < keys[k-1]) ok = 0;
	if (keys[k] != curveKeyAt (c, xx[order[k]], yy[order[k]], zz[order[k]])) ok = 0;
	if (keys[k] > top) top = keys[k];
      }
      if ((double) top >= span) ok = 0;
      sprintf (name, "%s: keys rise along the sorted order and fit the digits (%s)", curveInfo[c].name, typeName[type]);
      sprintf (msg, "largest key %.3f of the key range", top / span);
      report (name, ok, msg);
      free (order); free (keys);
    }
}

/* ---------------------------------------------------------------- 4. tetrahedra */

static double tetVolume (double t[4][3]) {
  double a[3], b[3], c[3];
  int d;
  for (d = 0; d < 3; ++d) { a[d] = t[1][d]-t[0][d]; b[d] = t[2][d]-t[0][d]; c[d] = t[3][d]-t[0][d]; }
  return fabs (a[0]*(b[1]*c[2]-b[2]*c[1]) - a[1]*(b[0]*c[2]-b[2]*c[0]) + a[2]*(b[0]*c[1]-b[1]*c[0])) / 6;
}

static int sharedCorners (double a[4][3], double b[4][3]) {
  int i, j, s = 0;
  for (i = 0; i < 4; ++i)
    for (j = 0; j < 4; ++j)
      if (a[i][0] == b[j][0] && a[i][1] == b[j][1] && a[i][2] == b[j][2]) s++;
  return s;
}

static void testTetra (void) {

  int levels, ok, d, shared[5], s;
  unsigned long long k, n;
  double cur[4][3], prev[4][3], first[4][3], mid[3], vol;
  char name[160], msg[200];

  for (levels = 3; levels <= 15; levels += 3) {
    n = 6ULL << (levels - 3);
    ok = 1;
    vol = 0;
    memset (shared, 0, sizeof shared);
    for (k = 0; k < n; ++k) {
      curveTetra (k, levels, cur);
      for (d = 0; d < 3; ++d) mid[d] = (cur[0][d] + cur[1][d] + cur[2][d] + cur[3][d]) / 4;
      if (curveTetraKey (mid, levels) != k) ok = 0;
      vol += tetVolume (cur);
      if (fabs (tetVolume (cur) * n - 1) > 1e-9) ok = 0;
      if (k == 0) memcpy (first, cur, sizeof cur);
      else shared[sharedCorners (prev, cur)]++;
      memcpy (prev, cur, sizeof cur);
    }
    sprintf (name, "sierpinski3: %llu equal tetrahedra fill the cube, each centre maps to its own key", n);
    sprintf (msg, "total volume %.6f", vol);
    report (name, ok && fabs (vol - 1) < 1e-9, msg);
    sprintf (name, "sierpinski3: consecutive tetrahedra touch (%llu cells)", n);
    sprintf (msg, "share a face %d, an edge %d, a corner %d, nothing %d", shared[3], shared[2], shared[1], shared[0]);
    report (name, shared[0] == 0, msg);
    s = sharedCorners (prev, first);
    sprintf (name, "sierpinski3: closed, last tetrahedron touches the first (%llu cells)", n);
    sprintf (msg, "%d shared corners", s);
    report (name, s > 0, msg);
  }
}

int main (void) {
  testReader();
  testNeighbors();
  testCubeCurves();
  testCurveSort();
  testTetra();
  return finish();
}
