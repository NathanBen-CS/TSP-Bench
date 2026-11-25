#include <stdio.h>
#include <stdlib.h>

#include "Instance.h"
#include "Utilities.h"
#include "tour.h"
#include "neighbors.h"
#include "diag.h"

/* Box diagnostics, any dimension (see diag.h).  Not part of any solver run
   that is timed: boxdiag and boxtsp3d --diag call it. */

static int *boxOfCity (const BoxSet *bs) {
  int i, j, *boxOf = (int *) malloc ((ncities + 1) * sizeof(int));
  if (!boxOf) fatal ("Out of memory in diagnostics");
  for (j = 0; j < bs->nb; ++j)
    for (i = bs->start[j]; i < bs->start[j + 1]; ++i) boxOf[bs->order[i]] = j;
  return boxOf;
}

double outsideShare (const BoxSet *bs, double *nearestOutside) {
  int c, k, *boxOf = boxOfCity (bs);
  long long out = 0, first = 0;
  for (c = 1; c <= ncities; ++c)
    for (k = 0; k < nbK; ++k)
      if (boxOf[neighbors (c)[k]] != boxOf[c]) {
	out++;
	if (k == 0) first++;
      }
  free (boxOf);
  *nearestOutside = (double) first / ncities;
  return (double) out / ((double) ncities * nbK);
}

long long crossingEdges (const BoxSet *bs, long long *length) {
  int i, u, v, *boxOf = boxOfCity (bs);
  long long n = 0;
  *length = 0;
  for (i = 0; i < ncities; ++i) {
    u = tour[i];
    v = tour[(i + 1) % ncities];
    if (boxOf[u] != boxOf[v]) { n++; *length += distance (u, v); }
  }
  free (boxOf);
  return n;
}

void printDiagnostics (const BoxSet *bs, const BoxStats *st) {
  long long len, n;
  double nearest, out = outsideShare (bs, &nearest);
  printf ("diag_boxes=%d\ndiag_avg_box=%.3f\ndiag_largest_box=%d\ndiag_capped_boxes=%d\n",
	  bs->nb, (double) ncities / bs->nb, st->largest, st->capped);
  printf ("diag_k=%d\ndiag_outside_k=%.4f\ndiag_outside_nearest=%.4f\n", nbK, out, nearest);
  if (tour) {
    n = crossingEdges (bs, &len);
    printf ("diag_cross_edges=%lld\ndiag_cross_share=%.4f\ndiag_cross_length_share=%.4f\n",
	    n, (double) n / ncities, len / tourLength());
  }
}
