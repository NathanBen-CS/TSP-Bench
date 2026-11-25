#include <math.h>
#include <stdio.h>
#include <stdlib.h>
#include <time.h>

#include "Instance.h"
#include "Utilities.h"
#include "tour.h"
#include "boxes.h"
#include "rings2d.h"

/*********************************************************************

Convex layers by repeated Andrew monotone chain.

  1. Sort the cities once by (x, y), then city number.  Cities with the
     same coordinates form one "point"; the hull is computed on points, and
     a point's cities are written next to each other in its ring.
  2. One pass per ring over the points still alive, in sorted order:
       lower chain left to right, upper chain right to left, popping only
       on a strict right turn (cross < 0), so points ON the boundary
       (collinear) are kept, in boundary order.
     lower + upper without their shared end points = the ring, counter-
     clockwise, starting at the lowest of the leftmost points.
  3. Drop the ring's points from the alive list and repeat.

Special cases: 1 or 2 points left = one ring of those points; all points
left on one line = one ring in sorted order (out along the line and back,
which is optimal for collinear points).

Cost: O(points alive) per ring, so O(N * rings) in total plus the sort.
On uniform data the number of rings grows like N^(2/3) (Dalal 2004), so
this is about N^(5/3): measured, not assumed (ring_time in the output).

Coordinates are doubles from the instance; the cross products are exact
for integer coordinates below about 2^26, which covers every map here.

 *********************************************************************/

static double secondsSince (clock_t start) {
  return (double)(clock() - start) / CLOCKS_PER_SEC;
}

static int byXY (const void *p, const void *q) {
  int a = *(const int *) p, b = *(const int *) q;
  if (xx[a] != xx[b]) return xx[a] < xx[b] ? -1 : 1;
  if (yy[a] != yy[b]) return yy[a] < yy[b] ? -1 : 1;
  return a - b;
}

static double *gx, *gy;            /* coordinates of each distinct point */

static double cross (int a, int b, int c) {     /* > 0: a -> b -> c turns left */
  return (gx[b] - gx[a]) * (gy[c] - gy[a]) - (gy[b] - gy[a]) * (gx[c] - gx[a]);
}

void ringBuild (BoxSet *bs, BoxStats *st, RingStats *rs) {

  int n = ncities, i, k, g, ng, m, h, hu, nb = 0, out = 0, s, lo;
  int *city, *gStart, *alive, *lower, *upper;
  char *used;
  clock_t clk;

  /* ---- 1. presort once by (x, y); group equal coordinates into points ---- */
  clk = clock();
  city = (int *) malloc (n * sizeof(int));
  gStart = (int *) malloc ((n + 1) * sizeof(int));
  gx = (double *) malloc (n * sizeof(double));
  gy = (double *) malloc (n * sizeof(double));
  if (!city || !gStart || !gx || !gy) fatal ("Out of memory in ringBuild");
  for (i = 0; i < n; ++i) city[i] = i + 1;
  qsort (city, n, sizeof(int), byXY);
  ng = 0;
  for (i = 0; i < n; ++i)
    if (i == 0 || xx[city[i]] != xx[city[i - 1]] || yy[city[i]] != yy[city[i - 1]]) {
      gStart[ng] = i;
      gx[ng] = xx[city[i]];
      gy[ng] = yy[city[i]];
      ng++;
    }
  gStart[ng] = n;
  st->sortTime = secondsSince (clk);

  /* ---- 2. peel ---- */
  clk = clock();
  alive = (int *) malloc (ng * sizeof(int));
  lower = (int *) malloc (ng * sizeof(int));
  upper = (int *) malloc (ng * sizeof(int));
  used  = (char *) calloc (ng, 1);
  bs->order = (int *) malloc (n * sizeof(int));
  bs->start = (int *) malloc ((n + 1) * sizeof(int));
  if (!alive || !lower || !upper || !used || !bs->order || !bs->start)
    fatal ("Out of memory in ringBuild");
  for (i = 0; i < ng; ++i) alive[i] = i;
  m = ng;

#define EMIT(G) do { int e_; used[G] = 1;                                \
    for (e_ = gStart[G]; e_ < gStart[(G) + 1]; ++e_) bs->order[out++] = city[e_]; } while (0)

  while (m > 0) {
    bs->start[nb++] = out;
    if (m <= 2) {                                  /* 1 or 2 points: one ring */
      for (i = 0; i < m; ++i) EMIT (alive[i]);
      break;
    }
    h = 0;                                         /* lower chain, left to right */
    for (i = 0; i < m; ++i) {
      g = alive[i];
      while (h >= 2 && cross (lower[h - 2], lower[h - 1], g) < 0) h--;
      lower[h++] = g;
    }
    if (h == m) {                                  /* every point on the lower    */
      for (i = 1; i < m - 1; ++i)                  /*   chain: are they all on    */
	if (cross (alive[0], alive[m - 1], alive[i]) != 0) break;   /* one line? */
      if (i >= m - 1) {                            /* yes: sorted order, the last ring */
	for (i = 0; i < m; ++i) EMIT (alive[i]);
	break;
      }
    }
    hu = 0;                                        /* upper chain, right to left */
    for (i = m - 1; i >= 0; --i) {
      g = alive[i];
      while (hu >= 2 && cross (upper[hu - 2], upper[hu - 1], g) < 0) hu--;
      upper[hu++] = g;
    }
    for (i = 0; i < h; ++i) EMIT (lower[i]);       /* counter-clockwise ring */
    for (i = 1; i < hu - 1; ++i) EMIT (upper[i]);
    k = 0;                                         /* keep the points still alive */
    for (i = 0; i < m; ++i)
      if (!used[alive[i]]) alive[k++] = alive[i];
    m = k;
  }
#undef EMIT
  bs->start[nb] = n;
  bs->nb = nb;
  if (out != n) fatal ("ringBuild: lost cities");
  st->splitTime = secondsSince (clk);

  /* ---- 3. each ring's loop is its order; stats ---- */
  bs->loops = (int *) malloc (n * sizeof(int));
  if (!bs->loops) fatal ("Out of memory in ringBuild");
  st->leaves = nb;
  st->capped = 0;
  st->largest = 0;
  st->loopsSum = 0;
  st->dpTime = 0;
  rs->rings = nb;
  rs->small = 0;
  rs->edges = 0;
  for (k = 0; k < nb; ++k) {
    lo = bs->start[k];
    s = bs->start[k + 1] - lo;
    if (s > st->largest) st->largest = s;
    if (s <= 2) rs->small++;
    if (s >= 2) rs->edges += s;
    for (i = 0; i < s; ++i) {
      bs->loops[lo + i] = bs->order[lo + i];
      if (s >= 2) st->loopsSum += distance (bs->order[lo + i], bs->order[lo + (i + 1) % s]);
    }
  }
  rs->largest = st->largest;
  rs->sortTime = st->sortTime;
  rs->peelTime = st->splitTime;

  free (city); free (gStart); free (gx); free (gy);
  free (alive); free (lower); free (upper); free (used);
}

/* ================================================================ hull order alone */

void hullOrderTour (const BoxSet *bs) {

  int j, lo, s, i, r, prev = -1, tn = 0;
  long long d, best;

  for (j = 0; j < bs->nb; ++j) {
    lo = bs->start[j];
    s = bs->start[j + 1] - lo;
    r = 0;                                         /* entry: nearest to the last */
    if (prev >= 0) {                               /*   city of the ring before  */
      best = -1;
      for (i = 0; i < s; ++i) {
	d = distance (prev, bs->loops[lo + i]);
	if (best < 0 || d < best) { best = d; r = i; }
      }
    }
    for (i = 0; i < s; ++i) tour[tn++] = bs->loops[lo + (r + i) % s];   /* once round */
    prev = tour[tn - 1];
  }
  tourSetPositions();
}

/* ================================================================ debug checks */

static double crossC (int a, int b, int c) {       /* same as cross(), on cities */
  return (xx[b] - xx[a]) * (yy[c] - yy[a]) - (yy[b] - yy[a]) * (xx[c] - xx[a]);
}

int ringCheck (const BoxSet *bs, int containment) {

  int n = ncities, j, i, lo, s, nv, a, b, c, p, q, fails = 0, collinear;
  int *v = (int *) malloc (n * sizeof(int));
  char *seen = (char *) calloc (n + 1, 1);
  double turn, dx1, dy1, dx2, dy2;

  for (i = 0; i < n; ++i) {                        /* every city exactly once */
    c = bs->order[i];
    if (c < 1 || c > n || seen[c]) { printf ("ring check: city %d twice or bad\n", c); fails++; }
    else seen[c] = 1;
    if (bs->loops[i] != c) { printf ("ring check: loop differs from order at %d\n", i); fails++; }
  }

  for (j = 0; j < bs->nb; ++j) {
    lo = bs->start[j];
    s = bs->start[j + 1] - lo;
    nv = 0;                    /* runs of equal coordinates -> one point each; a  */
    for (i = 0; i < s; ++i) {  /* point split into two runs (duplicates apart)   */
      c = bs->loops[lo + i];   /* fails the one-full-turn test below             */
      if (nv == 0 || xx[c] != xx[v[nv - 1]] || yy[c] != yy[v[nv - 1]]) v[nv++] = c;
    }
    while (nv > 1 && xx[v[nv - 1]] == xx[v[0]] && yy[v[nv - 1]] == yy[v[0]]) nv--;
    if (nv < 3) {
      if (j < bs->nb - 1) { printf ("ring check: ring %d has %d points but is not the last\n", j, nv); fails++; }
      continue;
    }
    collinear = 1;
    for (i = 2; i < nv; ++i) if (crossC (v[0], v[1], v[i]) != 0) collinear = 0;
    if (collinear) {
      if (j < bs->nb - 1) { printf ("ring check: ring %d is on one line but is not the last\n", j); fails++; }
      continue;
    }
    turn = 0;                                      /* no right turn; one full turn */
    for (i = 0; i < nv; ++i) {
      a = v[(i + nv - 1) % nv]; b = v[i]; c = v[(i + 1) % nv];
      if (crossC (a, b, c) < 0) { printf ("ring check: ring %d turns right at city %d\n", j, b); fails++; break; }
      dx1 = xx[b] - xx[a]; dy1 = yy[b] - yy[a];
      dx2 = xx[c] - xx[b]; dy2 = yy[c] - yy[b];
      turn += atan2 (dx1 * dy2 - dy1 * dx2, dx1 * dx2 + dy1 * dy2);
    }
    if (fabs (turn - 6.283185307179586) > 1e-6) {
      printf ("ring check: ring %d turns %.6f, not one full turn\n", j, turn); fails++;
    }
    if (containment)                               /* later cities strictly inside */
      for (p = bs->start[j + 1]; p < n; ++p) {
	c = bs->order[p];
	for (q = 0; q < nv; ++q)
	  if (crossC (v[q], v[(q + 1) % nv], c) <= 0) {
	    printf ("ring check: city %d of a later ring is not strictly inside ring %d\n", c, j);
	    fails++;
	    p = n;
	    break;
	  }
      }
  }
  free (v); free (seen);
  return fails;
}
