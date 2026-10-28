#include <stdlib.h>

#include "Instance.h"
#include "Utilities.h"
#include "tour.h"
#include "hilbert.h"

/*********************************************************************

Space-filling curve tour (Bartholdi & Platzman idea, using the Hilbert
curve rather than their Sierpinski curve).

A Hilbert curve is a single continuous line that passes through every cell
of a 2^ORDER x 2^ORDER grid, and it has a useful property: cells that are
close together along the curve are close together on the map.  So if we
  1. scale every city onto that grid,
  2. ask "how far along the curve is this city's cell?" (its index), and
  3. sort the cities by that index,
the sorted order is a reasonable tour.  It is ~25-40% above optimal, but it
needs no distance calculations at all and runs in O(N log N), so it works
even on the 10-million-city instances.

 *********************************************************************/

#define ORDER 16                         /* grid is 65536 x 65536 cells */

typedef struct {
  unsigned long long key;                /* position along the curve */
  int city;
} CurvePoint;

/*******************************************************************/

/* Position of grid cell (x, y) along the Hilbert curve on an n x n grid
   (n a power of 2).  Standard algorithm: look at one bit of x and y at a
   time, from the most significant down.  Each bit pair says which of the
   four quadrants the cell is in; quadrants are visited in the order
   lower-left, upper-left, upper-right, lower-right, so that quadrant adds
   0, 1, 2 or 3 times (quadrant size) to the index.  Then the coordinates are
   rotated/flipped so the next, smaller level sees its quadrant the right way
   up. */

unsigned long long hilbertIndex (unsigned x, unsigned y, unsigned n) {

  unsigned long long d = 0;
  unsigned s, rx, ry, t;

  for (s = n/2; s > 0; s /= 2) {
    rx = (x & s) > 0;                    /* right half at this level? */
    ry = (y & s) > 0;                    /* upper half at this level? */
    d += (unsigned long long) s * s * ((3 * rx) ^ ry);   /* quadrant offset */
    if (ry == 0) {                       /* lower quadrants are rotated: */
      if (rx == 1) {                     /*   lower-right is also flipped */
	x = n - 1 - x;
	y = n - 1 - y;
      }
      t = x; x = y; y = t;               /*   swap x and y (transpose) */
    }
  }
  return d;
}

/*******************************************************************/

static int byKey (const void *a, const void *b) {     /* qsort comparator */

  const CurvePoint *p = (const CurvePoint *) a;
  const CurvePoint *q = (const CurvePoint *) b;

  if (p->key != q->key) return (p->key < q->key) ? -1 : 1;   /* curve order */
  return p->city - q->city;              /* same cell: lower city number first */
}

/*******************************************************************/

/* Sort all cities along the Hilbert curve on a 2^order x 2^order grid.
   cities[i] = the i-th city along the curve; keys[i] = its curve position
   (keys may be NULL if not needed).  Used by hilbertTour (order 16) and by
   the box solver (order 20, see boxes.c). */

void hilbertSort (int order, int *cities, unsigned long long *keys) {

  double minx, maxx, miny, maxy, side, scale;
  unsigned n = 1u << order, gx, gy;
  CurvePoint *pts;
  int c, i;

  if (!isEuclidean) fatal ("Space-filling curve needs city coordinates");

  /* 1. bounding box of all cities */
  minx = maxx = xx[1];
  miny = maxy = yy[1];
  for (c=2; c<=ncities; ++c) {
    if (xx[c] < minx) minx = xx[c];
    if (xx[c] > maxx) maxx = xx[c];
    if (yy[c] < miny) miny = yy[c];
    if (yy[c] > maxy) maxy = yy[c];
  }

  /* 2. one scale for both axes, so the map isn't stretched */
  side = (maxx - minx > maxy - miny) ? maxx - minx : maxy - miny;
  scale = (side > 0) ? (n - 1) / side : 0;

  /* 3. each city's grid cell -> its position along the curve */
  pts = (CurvePoint *) malloc (ncities * sizeof(CurvePoint));
  if (pts == 0) fatal ("Out of memory in hilbertSort");
  for (c=1; c<=ncities; ++c) {
    gx = (unsigned) ((xx[c] - minx) * scale);
    gy = (unsigned) ((yy[c] - miny) * scale);
    pts[c-1].key  = hilbertIndex (gx, gy, n);
    pts[c-1].city = c;
  }

  /* 4. sort by curve position */
  qsort (pts, ncities, sizeof(CurvePoint), byKey);
  for (i=0; i<ncities; ++i) {
    cities[i] = pts[i].city;
    if (keys) keys[i] = pts[i].key;
  }

  free (pts);
}

/*******************************************************************/

void hilbertTour (void) {

  hilbertSort (ORDER, tour, 0);          /* the curve order is the tour */
  tourSetPositions();
}
