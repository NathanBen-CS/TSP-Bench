#include <stdlib.h>
#include <time.h>

#include "Instance.h"
#include "Utilities.h"
#include "tour.h"
#include "curves.h"
#include "merges.h"
#include "boxes.h"

/*********************************************************************

Box start tour, ported from the prototype proto.c and adapted to the bench:

  * distances come from the bench's distance() (DIMACS integer rounding),
    so every length here matches codes/length.exe;
  * a box that reaches the depth cap with more than t cities (this only
    happens when many cities share one grid cell, i.e. duplicate or
    near-identical points) is NOT given to the DP, which would overflow;
    those cities keep their curve order instead;
  * any curve from curves.c and any merge method from merges.c.

This is a heuristic in the spirit of Karp (1977): no worst-case guarantee.

 *********************************************************************/

static double secondsSince (clock_t start) {
  return (double)(clock() - start) / CLOCKS_PER_SEC;
}

/* ================================================================ 1. boxes

A space-filling curve visits its cells in order: in base b, the first
1/b of the key range is the first child cell, the next 1/b the second,
and so on, at every level.  So every cell is simply a contiguous run of
the key-sorted city list, and splitting a cell into its b children just
means cutting its run at b-1 key values.  No tree is ever built.        */

static const unsigned long long *sKeys;   /* sorted curve keys             */
static int  *sStart, *sLevel;             /* output: box starts / depths   */
static int   sCount, sDigits, sBase, sT;
static unsigned long long sPow[65];       /* sPow[k] = base^k              */

static void splitRange (int lo, int hi, int level, unsigned long long base) {

  unsigned long long q, end;
  int c, cur, j;

  if (hi == lo) return;                          /* empty cell: no box      */
  if (hi - lo <= sT || level == sDigits) {       /* small enough, or depth  */
    if (sLevel) sLevel[sCount] = level;          /*   cap reached: a box    */
    sStart[sCount++] = lo;
    return;
  }

  q = sPow[sDigits - level - 1];                 /* keys in each child cell */
  cur = lo;
  for (c = 0; c < sBase; ++c) {                  /* the children, in curve order */
    end = base + (c + 1) * q;                    /* first key of the next child  */
    j = cur;
    while (j < hi && sKeys[j] < end) j++;        /* cities of child c: cur..j-1  */
    splitRange (cur, j, level + 1, base + c * q);
    cur = j;
  }
}

int boxSplitBase (const unsigned long long *keys, int m, int base, int digits,
		  int t, int *start, int *level) {
  int k;
  sKeys = keys;
  sStart = start;
  sLevel = level;
  sCount = 0;
  sDigits = digits;
  sBase = base;
  sT = t;
  sPow[0] = 1;
  for (k = 1; k <= digits; ++k) sPow[k] = sPow[k - 1] * base;
  splitRange (0, m, 0, 0);      /* the root covers the whole key range */
  start[sCount] = m;            /* end marker: box j = start[j] .. start[j+1]-1 */
  return sCount;
}

int boxSplit (const unsigned long long *keys, int m, int bits, int t, int *start) {
  return boxSplitBase (keys, m, 4, bits, t, start, 0);
}

/* ================================================================ 2. exact loop per box

Held-Karp dynamic programming.  Fix city 0 of the box as the start.
best[S][j] = length of the shortest path that starts at city 0, visits
exactly the cities in set S (a bitmask over cities 1..s-1), and ends at
city j (which is in S).  Then

	best[{j}][j]  = d(0, j)
	best[S][j]    = min over i in S-{j} of  best[S-{j}][i] + d(i, j)
	loop length   = min over j of  best[all][j] + d(j, 0)

That is 2^(s-1) * (s-1) table entries, each taking O(s) work, so a box of
10 cities costs about 512 * 9 * 9 = 41,000 steps.  Exact, unlike everything
else in this folder.                                                       */

static long long best[1 << (BOX_MAXT - 1)][BOX_MAXT - 1];
static unsigned char parent[1 << (BOX_MAXT - 1)][BOX_MAXT - 1];

long long heldKarpLoop (const int *cities, int s, int *loop) {

  int d[BOX_MAXT][BOX_MAXT];
  int m = s - 1, full, mask, prev, i, j, bj, cur, c, k;
  long long v, b, total;
  int back[BOX_MAXT];

  if (s <= 3) {                       /* 1-3 cities: every order is the same loop */
    total = 0;
    for (i = 0; i < s; ++i) {
      loop[i] = cities[i];
      if (s > 1) total += distance (cities[i], cities[(i + 1) % s]);
    }
    return total;
  }

  for (i = 0; i < s; ++i)             /* distances inside the box, once */
    for (j = 0; j < s; ++j)
      d[i][j] = distance (cities[i], cities[j]);

  for (mask = 1; mask < (1 << m); ++mask)        /* every set S, small to big */
    for (j = 0; j < m; ++j) {                    /* every end city j in S     */
      if (!(mask & (1 << j))) continue;
      prev = mask & ~(1 << j);                   /* S without j               */
      if (prev == 0) {                           /* S = {j}: straight from 0  */
	best[mask][j] = d[0][j + 1];
	parent[mask][j] = 255;                   /* 255 = "came from city 0"  */
	continue;
      }
      b = -1;
      for (i = 0; i < m; ++i)                    /* best city i to come from  */
	if (prev & (1 << i)) {
	  v = best[prev][i] + d[i + 1][j + 1];
	  if (b < 0 || v < b) { b = v; parent[mask][j] = (unsigned char) i; }
	}
      best[mask][j] = b;
    }

  full = (1 << m) - 1;                           /* all cities visited:       */
  b = -1;
  bj = 0;
  for (j = 0; j < m; ++j) {                      /* close the loop back to 0  */
    v = best[full][j] + d[j + 1][0];
    if (b < 0 || v < b) { b = v; bj = j; }
  }

  /* walk the parent pointers backwards to recover the order */
  c = 0;
  cur = bj;
  mask = full;
  while (1) {
    back[c++] = cur + 1;                         /* box-local city number */
    prev = parent[mask][cur];
    if (prev == 255) break;
    mask &= ~(1 << cur);
    cur = prev;
  }
  loop[0] = cities[0];                           /* start at city 0 ...   */
  for (k = 0; k < c; ++k) loop[1 + k] = cities[back[c - 1 - k]];   /* ... forwards */
  return b;
}

/* ================================================================ 3. build the boxes */

void boxBuild (int curve, int t, BoxSet *bs, BoxStats *st) {

  unsigned long long *keys;
  clock_t clk;

  /* ---- sort the cities along the curve ---- */
  clk = clock();
  bs->order = (int *) malloc (ncities * sizeof(int));
  keys = (unsigned long long *) malloc (ncities * sizeof(unsigned long long));
  if (!bs->order || !keys) fatal ("Out of memory in boxBuild");
  curveSort (curve, bs->order, keys);
  st->sortTime = secondsSince (clk);

  /* ---- boxes = curve cells = runs of the sorted list ---- */
  clk = clock();
  bs->start = (int *) malloc ((ncities + 1) * sizeof(int));
  if (!bs->start) fatal ("Out of memory in boxBuild");
  bs->nb = boxSplitBase (keys, ncities, curveInfo[curve].base, curveInfo[curve].digits,
			 t, bs->start, 0);
  st->splitTime = secondsSince (clk);
  free (keys);

  boxSolve (t, bs, st);                            /* exact loop in every box */
}

void boxSolve (int t, BoxSet *bs, BoxStats *st) {

  int j, lo, s, k;
  clock_t clk = clock();

  st->leaves = bs->nb;
  bs->loops = (int *) malloc (ncities * sizeof(int));   /* box j's loop sits in */
  if (!bs->loops) fatal ("Out of memory in boxBuild");  /*   loops[start[j]..]  */
  st->loopsSum = 0;
  st->capped = 0;
  st->largest = 0;
  for (j = 0; j < bs->nb; ++j) {
    lo = bs->start[j];
    s = bs->start[j + 1] - lo;
    if (s > st->largest) st->largest = s;
    if (s <= t) {                                  /* normal box: exact loop  */
      st->loopsSum += heldKarpLoop (bs->order + lo, s, bs->loops + lo);
    }
    else {                                         /* depth-capped box: all its */
      st->capped++;                                /*   cities share one tiny   */
      for (k = 0; k < s; ++k) {                    /*   cell, keep curve order  */
	bs->loops[lo + k] = bs->order[lo + k];
	st->loopsSum += distance (bs->order[lo + k], bs->order[lo + (k + 1) % s]);
      }
    }
  }
  st->dpTime = secondsSince (clk);
}

void boxFree (BoxSet *bs) {
  free (bs->order);
  free (bs->start);
  free (bs->loops);
}

/* ================================================================ 4. the whole tour */

void boxTourWith (int curve, int merge, int t, BoxStats *st) {

  BoxSet bs;
  clock_t clk;

  boxBuild (curve, t, &bs, st);
  clk = clock();
  mergeLoops (merge, &bs, st);                   /* writes tour[] and pos[] */
  st->glueTime = secondsSince (clk);
  boxFree (&bs);
}

void boxTour (int t, BoxStats *st) {
  boxTourWith (CURVE_HILBERT, MERGE_M1, t, st);
}
