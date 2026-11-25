#include <stdlib.h>
#include <string.h>

#include "instance3d.h"
#include "Utilities.h"
#include "curves3d.h"

/*********************************************************************

Space-filling curves for the box solver (3D).  See curves.h, curves3d.h.

Cube-cell curves work on an integer grid (side x side x side cells) and
give each cell a key; all points of one curve cell share a key prefix,
which is what the box splitting (boxSplitBase) needs.

Keys fit in 64 bits: base 8 with 20 digits is 60 bits, base 27 with 13
digits is 3^39 < 2^62, base 2 with 60 digits is 60 bits.

 *********************************************************************/

const int nCurves = NCURVES3;

const CurveInfo curveInfo[] = {
  /* name          base digits closed tetrahedra */
  { "hilbert3",     8,   20,    0,     0 },
  { "morton3",      8,   20,    0,     0 },
  { "peano3",      27,   13,    0,     0 },   /* 3^13 = 1,594,323 per side */
  { "snake3",       8,   20,    0,     0 },
  { "moore3",       8,   20,    1,     0 },
  { "sierpinski3",  2,   60,    1,     1 },   /* 3 digits pick one of 6 tetrahedra */
};

int curveByName (const char *name) {
  int c;
  for (c = 0; c < NCURVES3; ++c)
    if (!strcmp (name, curveInfo[c].name)) return c;
  return -1;
}

unsigned long long curveGridSide3 (int curve, int levels) {
  unsigned long long side = 1;
  int i;
  for (i = 0; i < levels; ++i) side *= (curve == CURVE_PEANO3) ? 3 : 2;
  return side;
}

/* ================================================================ cube-cell curves */

/* Interleave `levels` bits of a, b, c into base-8 digits, a most significant. */
static unsigned long long interleave (unsigned a, unsigned b, unsigned c, int levels) {
  unsigned long long key = 0;
  int i;
  for (i = levels - 1; i >= 0; --i)
    key = key * 8 + ((((a >> i) & 1u) << 2) | (((b >> i) & 1u) << 1) | ((c >> i) & 1u));
  return key;
}

/* Hilbert3: J. Skilling, "Programming the Hilbert curve", AIP Conf. Proc.
   707 (2004), AxestoTranspose for 3 dimensions.  X[0] holds the most
   significant bit of every digit.  Starts at (0,0,0) and ends at
   (n-1,0,0); consecutive cells always share a face.                      */

static unsigned long long hilbert3Key (unsigned x, unsigned y, unsigned z, int levels) {

  unsigned X[3] = { x, y, z }, M = 1u << (levels - 1), P, Q, t;
  int i;

  for (Q = M; Q > 1; Q >>= 1) {                  /* undo the excess work   */
    P = Q - 1;
    for (i = 0; i < 3; ++i)
      if (X[i] & Q) X[0] ^= P;                   /*   invert low bits      */
      else {                                     /*   or exchange low bits */
	t = (X[0] ^ X[i]) & P;
	X[0] ^= t;
	X[i] ^= t;
      }
  }
  for (i = 1; i < 3; ++i) X[i] ^= X[i-1];        /* Gray encode            */
  t = 0;
  for (Q = M; Q > 1; Q >>= 1)
    if (X[2] & Q) t ^= Q - 1;
  for (i = 0; i < 3; ++i) X[i] ^= t;
  return interleave (X[0], X[1], X[2], levels);
}

/* Moore3: a closed curve made of eight Hilbert3 curves, one per octant (the
   2D Moore curve is four Hilbert curves).  The octants are visited in the
   cycle below; neighbours in the cycle differ in one axis (x y x z x y x z).
   Each octant's Hilbert3 curve runs parallel to the x axis, along the edge
   of the octant that lies on the central line y = z = middle, so the end
   of one octant touches the start of the next, and the last the first.    */

static const unsigned octant[8][3] = {
  {0,0,0}, {1,0,0}, {1,1,0}, {0,1,0}, {0,1,1}, {1,1,1}, {1,0,1}, {0,0,1}
};

static unsigned long long moore3Key (unsigned x, unsigned y, unsigned z, int levels) {

  unsigned n = 1u << (levels - 1);               /* side of one octant       */
  unsigned ox = x >= n, oy = y >= n, oz = z >= n;
  unsigned lx = x - ox * n, ly = y - oy * n, lz = z - oz * n;
  int q;

  for (q = 0; q < 8; ++q)                        /* octant's place in the cycle */
    if (octant[q][0] == ox && octant[q][1] == oy && octant[q][2] == oz) break;

  /* start corner: y and z at the central line; x at the outer face for even
     places, at the central plane for odd ones.  Mirror so it becomes 0.     */
  if ((q % 2 == 0) ? ox : !ox) lx = n - 1 - lx;
  if (!oy) ly = n - 1 - ly;
  if (!oz) lz = n - 1 - lz;
  if (levels == 1) return q;
  return (unsigned long long) q * n * n * n + hilbert3Key (lx, ly, lz, levels - 1);
}

/* Peano3: every cube is cut into 3 x 3 x 3, walked as a snake: x slowest,
   z fastest.  Peano's rule: an axis' digit is mirrored when the digits of
   the other two axes so far add up to an odd number, so consecutive cells
   always share a face.                                                   */

static unsigned long long peano3Key (unsigned x, unsigned y, unsigned z, int levels) {

  unsigned long long key = 0;
  unsigned p = 1, kx, ky, kz, sx = 0, sy = 0, sz = 0;   /* digit sums so far */
  int l;

  for (l = 1; l < levels; ++l) p *= 3;
  for (l = 0; l < levels; ++l, p /= 3) {
    kx = (x / p) % 3;  if ((sy + sz) & 1) kx = 2 - kx;  sx += kx;
    ky = (y / p) % 3;  if ((sx + sz) & 1) ky = 2 - ky;  sy += ky;
    kz = (z / p) % 3;  if ((sx + sy) & 1) kz = 2 - kz;  sz += kz;
    key = key * 27 + 9 * kx + 3 * ky + kz;
  }
  return key;
}

/* Snake3 (boustrophedon): row by row along x, the rows of a layer along y,
   the layers along z, every direction alternating.  The 3D counterpart of
   the 2D "strip" curve, and like it at full grid resolution: a layer holds
   about one city, so the order is close to "sort by z".                   */

static unsigned long long snake3Key (unsigned x, unsigned y, unsigned z, int levels) {
  unsigned long long side = 1ULL << levels, row;
  row = z * side + ((z & 1u) ? side - 1 - y : y);
  return row * side + ((row & 1u) ? side - 1 - x : x);
}

unsigned long long curveGridKey3 (int curve, unsigned x, unsigned y, unsigned z, int levels) {
  switch (curve) {
    case CURVE_HILBERT3: return hilbert3Key (x, y, z, levels);
    case CURVE_MORTON3:  return interleave (z, y, x, levels);
    case CURVE_PEANO3:   return peano3Key (x, y, z, levels);
    case CURVE_SNAKE3:   return snake3Key (x, y, z, levels);
    case CURVE_MOORE3:   return moore3Key (x, y, z, levels);
    default:             return 0;               /* tetrahedra: no grid */
  }
}

/* ================================================================ tetrahedra (Sierpinski3)

The 2D Sierpinski curve cuts a triangle in half again and again and walks the
halves so that the path never jumps.  This is the same idea in 3D.

The unit cube is cut into the 6 Kuhn tetrahedra: one per order of the axes
(a,b,c), with corners 0, e_a, e_a+e_b, (1,1,1).  They are visited in a cycle
in which neighbours share a face, alternately from (0,0,0) to (1,1,1) and
back, so the curve ends where it started: closed.

Each tetrahedron is cut in half at the midpoint m of one edge (Maubach's
bisection): a cell is an ordered list of corners v0 v1 v2 v3 with a type
k = 3, 2, 1; the cut edge is v0 - vk.

	half 0: v0 .. v(k-1), m, v(k+1) ..     type k-1 (after 1 comes 3)
	half 1: v1 .. vk,     m, v(k+1) ..

The curve enters a cell at one corner and leaves at another.  For each
(type, entry, exit) the table tetRule says which half comes first and at
which corner of the shared face the path crosses over (4 = the midpoint m),
so consecutive cells always share at least that corner.  Any such table gives
a continuous curve.  This one was picked by hill climbing: 6 random
tables improved one entry at a time for the shortest plain curve order through
8,000 uniform random points (60,000 other points gave the same value).  It
is our own choice, not one of the curves in Haverkort's papers.             */

static const int kuhn[6][3] = {
  {0,1,2}, {1,0,2}, {1,2,0}, {2,1,0}, {2,0,1}, {0,2,1}
};

static const signed char tetRule[4][4][4][2] = {   /* [type][entry][exit] = {first half, corner} */
  {{{0}}},
  { {{0,0},{0,3},{0,3},{0,2}}, {{1,3},{0,0},{1,3},{1,2}}, {{1,3},{0,3},{0,0},{0,4}}, {{1,4},{0,4},{1,4},{0,0}} },
  { {{0,0},{0,4},{0,1},{0,1}}, {{1,4},{0,0},{0,3},{0,4}}, {{1,1},{1,4},{0,0},{1,1}}, {{1,1},{1,4},{0,4},{0,0}} },
  { {{0,0},{0,2},{0,4},{0,1}}, {{1,4},{0,0},{1,4},{0,4}}, {{1,4},{0,4},{0,0},{0,1}}, {{1,1},{1,4},{1,1},{0,0}} }
};

typedef struct { double v[4][3]; int type, in, out; } Tet;

static void tetRoot (int r, Tet *T) {
  int i, d;
  for (d = 0; d < 3; ++d) T->v[0][d] = 0;
  for (i = 0; i < 3; ++i) {
    for (d = 0; d < 3; ++d) T->v[i+1][d] = T->v[i][d];
    T->v[i+1][kuhn[r][i]] = 1;
  }
  T->type = 3;
  T->in  = (r & 1) ? 3 : 0;                      /* every second one runs back */
  T->out = (r & 1) ? 0 : 3;
}

/* Where parent corner i (4 = midpoint) sits in the corner list of a half. */
static int tetCorner (int k, int half, int i) {
  if (i == 4) return k;
  return (i > k || half == 0) ? i : i - 1;
}

/* Cut T and keep half 0 (with v0) or half 1 (with v_k).  Returns the key bit:
   0 if the curve visits that half first, 1 if second. */
static int tetCut (Tet *T, int half) {
  int d, i, k = T->type, bit;
  const signed char *rule = tetRule[k][T->in][T->out];
  Tet N;
  bit = (half != rule[0]);
  for (i = 0; i < 4; ++i)
    for (d = 0; d < 3; ++d)
      N.v[i][d] = (half == 1 && i < k) ? T->v[i+1][d] : T->v[i][d];
  for (d = 0; d < 3; ++d) N.v[k][d] = (T->v[0][d] + T->v[k][d]) / 2;
  N.type = (k == 1) ? 3 : k - 1;
  N.in  = tetCorner (k, half, bit ? rule[1] : T->in);
  N.out = tetCorner (k, half, bit ? T->out : rule[1]);
  *T = N;
  return bit;
}

/* Is p on v0's side of the cutting plane (the plane through m and the two
   corners that are not on the cut edge)?  Points on the plane count as v0's. */
static int tetSide (const Tet *T, const double p[3]) {
  double m[3], a[3], b[3], n[3], sp, s0;
  const double *u[2];
  int d, i, c = 0, k = T->type;
  for (d = 0; d < 3; ++d) m[d] = (T->v[0][d] + T->v[k][d]) / 2;
  for (i = 1; i < 4; ++i) if (i != k) u[c++] = T->v[i];
  for (d = 0; d < 3; ++d) { a[d] = u[0][d] - m[d]; b[d] = u[1][d] - m[d]; }
  n[0] = a[1]*b[2] - a[2]*b[1];
  n[1] = a[2]*b[0] - a[0]*b[2];
  n[2] = a[0]*b[1] - a[1]*b[0];
  sp = s0 = 0;
  for (d = 0; d < 3; ++d) {
    sp += n[d] * (p[d] - m[d]);
    s0 += n[d] * (T->v[0][d] - m[d]);
  }
  return sp == 0 || (sp > 0) == (s0 > 0);
}

unsigned long long curveTetraKey (const double p[3], int levels) {
  Tet T;
  unsigned long long key;
  int r, l, half;
  for (r = 0; r < 6; ++r)                        /* which Kuhn tetrahedron */
    if (p[kuhn[r][0]] >= p[kuhn[r][1]] && p[kuhn[r][1]] >= p[kuhn[r][2]]) break;
  tetRoot (r, &T);
  key = r;
  for (l = 3; l < levels; ++l) {
    half = tetSide (&T, p) ? 0 : 1;
    key = key * 2 + tetCut (&T, half);
  }
  return key;
}

int curveTetra (unsigned long long key, int levels, double tet[4][3]) {
  Tet T;
  int l, first, r = (int) (key >> (levels - 3));
  if (r > 5) return 0;
  tetRoot (r, &T);
  for (l = levels - 4; l >= 0; --l) {
    first = tetRule[T.type][T.in][T.out][0];
    tetCut (&T, ((key >> l) & 1) ? 1 - first : first);
  }
  memcpy (tet, T.v, sizeof T.v);
  return 1;
}

/* ================================================================ keys and sorting */

static double minx, miny, minz, extent;          /* the frame: cities' bounding box */

void curveFrame (void) {
  double maxx, maxy, maxz;
  int c;
  if (!isEuclidean) fatal ("Space-filling curves need city coordinates");
  minx = maxx = xx[1];
  miny = maxy = yy[1];
  minz = maxz = zz[1];
  for (c = 2; c <= ncities; ++c) {
    if (xx[c] < minx) minx = xx[c];
    if (xx[c] > maxx) maxx = xx[c];
    if (yy[c] < miny) miny = yy[c];
    if (yy[c] > maxy) maxy = yy[c];
    if (zz[c] < minz) minz = zz[c];
    if (zz[c] > maxz) maxz = zz[c];
  }
  extent = maxx - minx;                          /* one scale for all axes */
  if (maxy - miny > extent) extent = maxy - miny;
  if (maxz - minz > extent) extent = maxz - minz;
}

unsigned long long curveKeyAt (int curve, double x, double y, double z) {
  int digits = curveInfo[curve].digits;
  double scale, p[3];
  if (curveInfo[curve].triangles) {
    p[0] = extent > 0 ? (x - minx) / extent : 0;
    p[1] = extent > 0 ? (y - miny) / extent : 0;
    p[2] = extent > 0 ? (z - minz) / extent : 0;
    return curveTetraKey (p, digits);
  }
  scale = extent > 0 ? (curveGridSide3 (curve, digits) - 1) / extent : 0;
  return curveGridKey3 (curve, (unsigned) ((x - minx) * scale), (unsigned) ((y - miny) * scale),
			(unsigned) ((z - minz) * scale), digits);
}

typedef struct { unsigned long long key; int city; } KeyedCity;

static int byKey (const void *a, const void *b) {               /* qsort comparator */
  const KeyedCity *p = (const KeyedCity *) a, *q = (const KeyedCity *) b;
  if (p->key != q->key) return (p->key < q->key) ? -1 : 1;     /* curve order      */
  return p->city - q->city;                                    /* ties: city number */
}

void curveSort (int curve, int *cities, unsigned long long *keys) {

  KeyedCity *kc;
  int c, i;

  curveFrame();
  kc = (KeyedCity *) malloc (ncities * sizeof(KeyedCity));
  if (!kc) fatal ("Out of memory in curveSort");
  for (c = 1; c <= ncities; ++c) {
    kc[c-1].key = curveKeyAt (curve, xx[c], yy[c], zz[c]);
    kc[c-1].city = c;
  }
  qsort (kc, ncities, sizeof(KeyedCity), byKey);
  for (i = 0; i < ncities; ++i) {
    cities[i] = kc[i].city;
    if (keys) keys[i] = kc[i].key;
  }
  free (kc);
}
