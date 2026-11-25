#include <stdlib.h>
#include <time.h>

#include "instance3d.h"
#include "Utilities.h"
#include "curves3d.h"
#include "kdboxes3d.h"

/* See kdboxes3d.h. */

static double secondsSince (clock_t start) {
  return (double)(clock() - start) / CLOCKS_PER_SEC;
}

static int *list;                /* the cities, cut in place             */
static int *cutStart, nCuts;     /* piece j = list[cutStart[j] ..]       */
static int  kdT;

static double coord (int c, int d) { return d == 0 ? xx[c] : d == 1 ? yy[c] : zz[c]; }

static int sortAxis;             /* qsort comparator: one axis, ties by city */
static int byAxis (const void *a, const void *b) {
  int p = *(const int *) a, q = *(const int *) b;
  double x = coord (p, sortAxis), y = coord (q, sortAxis);
  if (x != y) return x < y ? -1 : 1;
  return p - q;
}

static void cut (int first, int last) {

  double lo[3], hi[3], v;
  int i, d, widest, n = last - first, pieces, mid;

  if (n <= kdT) { cutStart[nCuts++] = first; return; }
  for (d = 0; d < 3; ++d) lo[d] = hi[d] = coord (list[first], d);
  for (i = first + 1; i < last; ++i)
    for (d = 0; d < 3; ++d) {
      v = coord (list[i], d);
      if (v < lo[d]) lo[d] = v;
      if (v > hi[d]) hi[d] = v;
    }
  widest = 0;
  for (d = 1; d < 3; ++d) if (hi[d] - lo[d] > hi[widest] - lo[widest]) widest = d;

  sortAxis = widest;
  qsort (list + first, n, sizeof(int), byAxis);
  pieces = (n + kdT - 1) / kdT;                  /* pieces this run will give; */
  mid = first + (pieces / 2) * kdT;              /*   the left half is full    */
  cut (first, mid);
  cut (mid, last);
}

typedef struct { unsigned long long key; int piece; } Keyed;

static int byKey (const void *a, const void *b) {
  const Keyed *p = (const Keyed *) a, *q = (const Keyed *) b;
  if (p->key != q->key) return p->key < q->key ? -1 : 1;
  return p->piece - q->piece;
}

void boxBuildKd (int curve, int t, BoxSet *bs, BoxStats *st) {

  Keyed *kd;
  int i, j, k, n, lo, s;
  double c[3];
  clock_t clk = clock();

  /* ---- 1. cut the cities into pieces of t ---- */
  kdT = t;
  list = (int *) malloc (ncities * sizeof(int));
  cutStart = (int *) malloc ((ncities + 1) * sizeof(int));
  if (!list || !cutStart) fatal ("Out of memory in boxBuildKd");
  for (i = 0; i < ncities; ++i) list[i] = i + 1;
  nCuts = 0;
  cut (0, ncities);
  cutStart[nCuts] = ncities;
  st->splitTime = secondsSince (clk);

  /* ---- 2. order the pieces by the curve position of their centres ---- */
  clk = clock();
  curveFrame();
  kd = (Keyed *) malloc (nCuts * sizeof(Keyed));
  for (j = 0; j < nCuts; ++j) {
    c[0] = c[1] = c[2] = 0;
    s = cutStart[j + 1] - cutStart[j];
    for (i = cutStart[j]; i < cutStart[j + 1]; ++i)
      for (k = 0; k < 3; ++k) c[k] += coord (list[i], k);
    kd[j].key = curveKeyAt (curve, c[0] / s, c[1] / s, c[2] / s);
    kd[j].piece = j;
  }
  qsort (kd, nCuts, sizeof(Keyed), byKey);

  bs->nb = nCuts;
  bs->order = (int *) malloc (ncities * sizeof(int));
  bs->start = (int *) malloc ((ncities + 1) * sizeof(int));
  if (!bs->order || !bs->start) fatal ("Out of memory in boxBuildKd");
  n = 0;
  for (j = 0; j < nCuts; ++j) {
    bs->start[j] = n;
    lo = cutStart[kd[j].piece];
    s = cutStart[kd[j].piece + 1] - lo;
    for (k = 0; k < s; ++k) bs->order[n++] = list[lo + k];
  }
  bs->start[nCuts] = n;
  free (kd);
  free (list);
  free (cutStart);
  st->sortTime = secondsSince (clk);

  boxSolve (t, bs, st);                          /* ---- 3. exact loops ---- */
}
