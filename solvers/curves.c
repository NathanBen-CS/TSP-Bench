#include <stdlib.h>
#include <string.h>

#include "Instance.h"
#include "Utilities.h"
#include "hilbert.h"
#include "curves.h"

/*********************************************************************

Space-filling curves for the box solver (2D).  See curves.h.

Square-cell curves work on an integer grid (side x side cells) and give
each cell a key; triangle curves work on the unit square directly by
cutting triangles in half.  Either way, all points of one curve cell
share a key prefix, which is what the box splitting needs.

Excluded from the plan's list (with reasons, see docs/box_tsp_2d.md):
Beta-Omega (no verified xy -> index formula found), H-index (optional,
not built), Gosper and Dragon (cannot tile a square / not one-to-one).

 *********************************************************************/

const CurveInfo curveInfo[NCURVES] = {
  /* name          base digits closed triangles */
  { "hilbert",     4,   20,    0,     0 },
  { "moore",       4,   20,    1,     0 },
  { "sierpinski",  2,   41,    1,     1 },   /* 1 bit: which half, then 40 */
  { "knopp",       2,   42,    0,     1 },   /* one triangle twice the square */
  { "peano",       9,   13,    0,     0 },   /* 3^13 = 1,594,323 per side */
  { "morton",      4,   20,    0,     0 },
  { "strip",       4,   20,    0,     0 },
};

int curveByName (const char *name) {
  int c;
  for (c = 0; c < NCURVES; ++c)
    if (!strcmp (name, curveInfo[c].name)) return c;
  return -1;
}

/* ================================================================ square-cell curves */

unsigned long long curveGridSide (int curve, int levels) {
  unsigned long long side = 1;
  int i;
  for (i = 0; i < levels; ++i)
    side *= (curve == CURVE_PEANO) ? 3 : 2;      /* Peano cells split 3 ways per side */
  return side;
}

/* Moore curve: a closed loop made of four Hilbert curves, one per quadrant.
   The standard Hilbert curve (hilbert.c) runs along the bottom of its
   square, from (0,0) to (n-1,0).  The Moore curve visits the quadrants
   lower-left, upper-left, upper-right, lower-right; the two left ones run
   upward along their right edge, the two right ones downward along their
   left edge, so it starts just left of the bottom centre and ends just
   right of it: the two ends are neighbours, so the curve is closed.     */

static unsigned long long mooreKey (unsigned x, unsigned y, int levels) {

  unsigned n = 1u << (levels - 1);               /* side of one quadrant     */
  int right = x >= n, top = y >= n, q;
  unsigned lx = right ? x - n : x;               /* position inside quadrant */
  unsigned ly = top ? y - n : y;
  unsigned hx, hy;

  q = !right ? (top ? 1 : 0) : (top ? 2 : 3);    /* quadrant's place in order */
  if (!right) { hx = ly;         hy = n - 1 - lx; }   /* undo a 90 degree turn  */
  else        { hx = n - 1 - ly; hy = lx;         }   /*   (left / right side)  */
  return (unsigned long long) q * n * n + hilbertIndex (hx, hy, n);
}

/* Morton (Z-order): interleave the bits of y and x.  Cells are visited in
   a Z inside every square, so the curve jumps: the negative control.    */

static unsigned long long mortonKey (unsigned x, unsigned y, int levels) {
  unsigned long long key = 0;
  int i;
  for (i = levels - 1; i >= 0; --i)
    key = key * 4 + ((((y >> i) & 1u) << 1) | ((x >> i) & 1u));
  return key;
}

/* Strip ("snake rows"): row by row, alternating left-to-right and
   right-to-left.  Trivial, and boxes become long thin strips.           */

static unsigned long long stripKey (unsigned x, unsigned y, int levels) {
  unsigned long long side = 1ULL << levels;
  return (unsigned long long) y * side + ((y & 1u) ? side - 1 - x : x);
}

/* Peano: each square is cut into 3 x 3.  The basic pattern walks the
   columns as a snake: column 0 upward, column 1 downward, column 2
   upward.  Inside cell (dx, dy) the smaller copy is mirrored left-right
   when dy is odd and up-down when dx is odd, so each copy starts where
   the previous one ended.  The mirrors pile up level by level (fx, fy). */

static unsigned long long peanoKey (unsigned x, unsigned y, int levels) {

  unsigned long long key = 0, p = 1;
  unsigned dx, dy, oy;
  int l, fx = 0, fy = 0;

  for (l = 1; l < levels; ++l) p *= 3;           /* size of a top-level cell */
  for (l = 0; l < levels; ++l, p /= 3) {
    dx = (unsigned) ((x / p) % 3);               /* which of the 3 x 3 cells */
    dy = (unsigned) ((y / p) % 3);
    if (fx) dx = 2 - dx;                         /* apply the mirrors so far */
    if (fy) dy = 2 - dy;
    oy = (dx % 2 == 0) ? dy : 2 - dy;            /* columns 0, 2 up; 1 down  */
    key = key * 9 + 3 * dx + oy;                 /* cell's place in the walk */
    if (dx % 2 == 1) fy ^= 1;                    /* mirror the smaller copy  */
    if (dy % 2 == 1) fx ^= 1;
  }
  return key;
}

unsigned long long curveGridKey (int curve, unsigned x, unsigned y, int levels) {
  switch (curve) {
    case CURVE_HILBERT: return hilbertIndex (x, y, 1u << levels);
    case CURVE_MOORE:   return mooreKey (x, y, levels);
    case CURVE_MORTON:  return mortonKey (x, y, levels);
    case CURVE_STRIP:   return stripKey (x, y, levels);
    case CURVE_PEANO:   return peanoKey (x, y, levels);
    default:            return 0;                /* triangle curves: no grid */
  }
}

/* ================================================================ triangle curves

A right isosceles triangle is stored as P, R, Q: the curve enters at P,
leaves at Q (the ends of the long side) and R is the right-angle corner.
Cutting from R to the middle M of PQ gives two smaller triangles of the
same shape:

	first half:  enters at P, leaves at R, right angle at M
	second half: enters at R, leaves at Q, right angle at M

so the curve runs through the first half, then the second, and each half
is cut the same way again.  One key bit per cut: 0 = first half.          */

static void halve (double P[2], double R[2], double Q[2], int bit) {
  double M[2] = { (P[0] + Q[0]) / 2, (P[1] + Q[1]) / 2 };
  if (bit == 0) { Q[0] = R[0]; Q[1] = R[1]; }     /* first half:  P -> R */
  else          { P[0] = R[0]; P[1] = R[1]; }     /* second half: R -> Q */
  R[0] = M[0]; R[1] = M[1];
}

static unsigned long long triDescend (double u, double v, double P[2], double R[2],
				      double Q[2], int levels, unsigned long long key) {
  int l, bit;
  double M[2], side, sideP;
  for (l = 0; l < levels; ++l) {
    M[0] = (P[0] + Q[0]) / 2;
    M[1] = (P[1] + Q[1]) / 2;
    /* which side of the cut R-M is the point on?  P's side = first half */
    side  = (M[0] - R[0]) * (v - R[1])    - (M[1] - R[1]) * (u - R[0]);
    sideP = (M[0] - R[0]) * (P[1] - R[1]) - (M[1] - R[1]) * (P[0] - R[0]);
    bit = (side == 0 || (side > 0) == (sideP > 0)) ? 0 : 1;   /* on the cut: first */
    halve (P, R, Q, bit);
    key = key * 2 + bit;
  }
  return key;
}

/* Sierpinski (Platzman & Bartholdi): the square is cut along its diagonal
   into two triangles.  The curve runs through the lower-right one from
   (0,0) to (1,1), then the upper-left one from (1,1) back to (0,0), so it
   ends where it started: closed.                                        */

static void sierpinskiStart (int half, double P[2], double R[2], double Q[2]) {
  if (half == 0) { P[0] = 0; P[1] = 0; R[0] = 1; R[1] = 0; Q[0] = 1; Q[1] = 1; }
  else           { P[0] = 1; P[1] = 1; R[0] = 0; R[1] = 1; Q[0] = 0; Q[1] = 0; }
}

/* Knopp: one triangle, right angle at (0,0), long side from (0,2) to
   (2,0), so it covers the whole unit square (and twice its area).  The
   curve's two ends are far apart: open.                                 */

static void knoppStart (double P[2], double R[2], double Q[2]) {
  P[0] = 0; P[1] = 2; R[0] = 0; R[1] = 0; Q[0] = 2; Q[1] = 0;
}

void curveTriangle (int curve, unsigned long long key, int levels, double tri[3][2]) {
  double P[2], R[2], Q[2];
  int l, first = levels - 1;
  if (curve == CURVE_SIERPINSKI) {
    sierpinskiStart ((int) ((key >> first) & 1), P, R, Q);   /* top bit: which half */
    first--;
  }
  else knoppStart (P, R, Q);
  for (l = first; l >= 0; --l) halve (P, R, Q, (int) ((key >> l) & 1));
  tri[0][0] = P[0]; tri[0][1] = P[1];
  tri[1][0] = R[0]; tri[1][1] = R[1];
  tri[2][0] = Q[0]; tri[2][1] = Q[1];
}

/* ================================================================ keys and sorting */

unsigned long long curveKey (int curve, double u, double v) {

  double P[2], R[2], Q[2];
  unsigned long long side;
  int half;

  if (curve == CURVE_SIERPINSKI) {
    half = (v <= u) ? 0 : 1;                     /* below the diagonal: first */
    sierpinskiStart (half, P, R, Q);
    return triDescend (u, v, P, R, Q, curveInfo[curve].digits - 1, (unsigned long long) half);
  }
  if (curve == CURVE_KNOPP) {
    knoppStart (P, R, Q);
    return triDescend (u, v, P, R, Q, curveInfo[curve].digits, 0);
  }
  side = curveGridSide (curve, curveInfo[curve].digits);
  return curveGridKey (curve, (unsigned) (u * (side - 1)), (unsigned) (v * (side - 1)),
		       curveInfo[curve].digits);
}

typedef struct { unsigned long long key; int city; } KeyedCity;

static int byKey (const void *a, const void *b) {               /* qsort comparator */
  const KeyedCity *p = (const KeyedCity *) a, *q = (const KeyedCity *) b;
  if (p->key != q->key) return (p->key < q->key) ? -1 : 1;     /* curve order      */
  return p->city - q->city;                                    /* ties: city number */
}

void curveSort (int curve, int *cities, unsigned long long *keys) {

  double minx, maxx, miny, maxy, extent, scale;
  unsigned long long side;
  KeyedCity *kc;
  int c, i;

  if (!isEuclidean) fatal ("Space-filling curves need city coordinates");

  /* 1. bounding box, and one scale for both axes (no stretching) */
  minx = maxx = xx[1];
  miny = maxy = yy[1];
  for (c = 2; c <= ncities; ++c) {
    if (xx[c] < minx) minx = xx[c];
    if (xx[c] > maxx) maxx = xx[c];
    if (yy[c] < miny) miny = yy[c];
    if (yy[c] > maxy) maxy = yy[c];
  }
  extent = (maxx - minx > maxy - miny) ? maxx - minx : maxy - miny;

  /* 2. every city's key */
  kc = (KeyedCity *) malloc (ncities * sizeof(KeyedCity));
  if (!kc) fatal ("Out of memory in curveSort");
  if (curveInfo[curve].triangles) {              /* triangle curves: unit square */
    for (c = 1; c <= ncities; ++c) {
      double u = extent > 0 ? (xx[c] - minx) / extent : 0;
      double v = extent > 0 ? (yy[c] - miny) / extent : 0;
      kc[c-1].key = curveKey (curve, u, v);
      kc[c-1].city = c;
    }
  }
  else {                                         /* grid curves: same scaling as */
    side = curveGridSide (curve, curveInfo[curve].digits);   /* hilbertSort       */
    scale = (extent > 0) ? (side - 1) / extent : 0;
    for (c = 1; c <= ncities; ++c) {
      unsigned gx = (unsigned) ((xx[c] - minx) * scale);
      unsigned gy = (unsigned) ((yy[c] - miny) * scale);
      kc[c-1].key = curveGridKey (curve, gx, gy, curveInfo[curve].digits);
      kc[c-1].city = c;
    }
  }

  /* 3. sort along the curve */
  qsort (kc, ncities, sizeof(KeyedCity), byKey);
  for (i = 0; i < ncities; ++i) {
    cities[i] = kc[i].city;
    if (keys) keys[i] = kc[i].key;
  }
  free (kc);
}
