#include <stdlib.h>

#include "instance3d.h"
#include "Utilities.h"
#include "neighbors.h"

/*********************************************************************

3D nearest-neighbour lists (same API as d2/neighbors2d.c).

2D uses a uniform grid.  A uniform 3D grid is a bad fit for the star
instances: hyg119614 has 109,399 of its 119,614 stars within 1% of the
width of the bounding box, so nearly all of them would share a few cells.
So 3D uses a kd-tree instead, which adapts to the density:

  * Build: cut the cities at the median of the widest axis, again and
    again, until a piece holds at most LEAF cities.  Every cell keeps its
    tight bounding box.
  * Search: go down to the query's own leaf first, then back up, and skip
    any cell whose bounding box is no closer than the k-th best so far.

The cities are copied in tree order, and the queries run in tree order,
so nearby searches touch nearby memory (same trick as the 2D grid).

Same conventions as 2D: lists are sorted by true (unrounded) distance, the
city itself is left out, ties keep whichever was found first, and the
result is deterministic.  With more than k exact duplicates of a city, its
list is filled by duplicates, as in 2D.

 *********************************************************************/

#define LEAF 8

int  nbK;
int *nbList;

/* ---------------------------------------------------------------- the tree */

typedef struct {
  double lo[3], hi[3];           /* bounding box of the cell's cities  */
  int first, last;               /* its cities: items first .. last-1  */
  int left, right;               /* children, or -1, -1 for a leaf     */
} Node;

typedef struct { double x[3]; int idx; } Item;     /* idx = index in cities[] */

static Node *nodes;
static int   nNodes;
static Item *items;

/* Ties are broken by idx, so the order along an axis is total and the tree
   is deterministic. */
static int before (const Item *p, const Item *q, int d) {
  if (p->x[d] != q->x[d]) return p->x[d] < q->x[d];
  return p->idx < q->idx;
}

static void swapItems (int i, int j) { Item t = items[i]; items[i] = items[j]; items[j] = t; }

/* Partially sort items[a..b-1] so that items[m] is in its sorted place along
   axis d (quickselect with a median-of-three pivot). */
static void nthElement (int a, int b, int m, int d) {
  int i, j, mid;
  Item pivot;
  while (b - a > 1) {
    mid = a + (b - a) / 2;
    if (before (&items[mid], &items[a], d))   swapItems (mid, a);
    if (before (&items[b-1], &items[a], d))   swapItems (b-1, a);
    if (before (&items[b-1], &items[mid], d)) swapItems (b-1, mid);
    pivot = items[mid];
    i = a;
    j = b - 1;
    while (i <= j) {
      while (before (&items[i], &pivot, d)) i++;
      while (before (&pivot, &items[j], d)) j--;
      if (i <= j) { swapItems (i, j); i++; j--; }
    }
    if (m <= j) b = j + 1;                          /* keep the side holding m */
    else if (m >= i) a = i;
    else return;
  }
}

static int build (int first, int last) {

  int me = nNodes++, d, widest, i, mid;
  Node *n = &nodes[me];

  for (d = 0; d < 3; ++d) n->lo[d] = n->hi[d] = items[first].x[d];
  for (i = first + 1; i < last; ++i)
    for (d = 0; d < 3; ++d) {
      if (items[i].x[d] < n->lo[d]) n->lo[d] = items[i].x[d];
      if (items[i].x[d] > n->hi[d]) n->hi[d] = items[i].x[d];
    }
  n->first = first;
  n->last = last;
  n->left = n->right = -1;
  if (last - first <= LEAF) return me;

  widest = 0;
  for (d = 1; d < 3; ++d)
    if (n->hi[d] - n->lo[d] > n->hi[widest] - n->lo[widest]) widest = d;
  mid = first + (last - first) / 2;
  nthElement (first, last, mid, widest);
  n->left = build (first, mid);
  n->right = build (mid, last);
  return me;
}

/* ---------------------------------------------------------------- one query */

static const double *qx;         /* query point                           */
static int     qIdx, qK, qCount; /* its index, k, how many found so far   */
static int    *qIds;             /* best so far: indices in cities[] ...  */
static double *qDist;            /* ... and squared distances, ascending  */

static double boxDist (const Node *n) {             /* squared distance to a cell */
  double s = 0, t;
  int d;
  for (d = 0; d < 3; ++d) {
    t = (qx[d] < n->lo[d]) ? n->lo[d] - qx[d] : (qx[d] > n->hi[d]) ? qx[d] - n->hi[d] : 0;
    s += t * t;
  }
  return s;
}

static void search (int node) {

  const Node *n = &nodes[node];
  int i, j, nearChild, farChild;
  double dx, dy, dz, d, dNear, dFar;

  if (n->left < 0) {                                /* leaf: try every city */
    for (i = n->first; i < n->last; ++i) {
      if (items[i].idx == qIdx) continue;
      dx = items[i].x[0] - qx[0];
      dy = items[i].x[1] - qx[1];
      dz = items[i].x[2] - qx[2];
      d = dx*dx + dy*dy + dz*dz;
      j = qCount;
      if (j == qK) {                                /* list full            */
	if (d >= qDist[qK-1]) continue;
	j = qK - 1;
      }
      else qCount++;
      while (j > 0 && qDist[j-1] > d) { qDist[j] = qDist[j-1]; qIds[j] = qIds[j-1]; j--; }
      qDist[j] = d;
      qIds[j] = items[i].idx;
    }
    return;
  }

  nearChild = n->left;                              /* nearer child first   */
  farChild = n->right;
  dNear = boxDist (&nodes[nearChild]);
  dFar = boxDist (&nodes[farChild]);
  if (dFar < dNear) {
    j = nearChild; nearChild = farChild; farChild = j;
    d = dNear; dNear = dFar; dFar = d;
  }
  if (qCount < qK || dNear < qDist[qK-1]) search (nearChild);
  if (qCount < qK || dFar  < qDist[qK-1]) search (farChild);
}

/* ---------------------------------------------------------------- API */

/* For each of the m cities in cities[], its k nearest cities among cities[]
   itself, nearest first, in out[i*k .. i*k+k-1].  Requires k <= m-1. */

void nearestAmong (const int *cities, int m, int k, int *out) {

  int i, j, c;

  if (!isEuclidean) fatal ("nearestAmong needs city coordinates");
  items = (Item *) malloc ((size_t) m * sizeof(Item));
  nodes = (Node *) malloc (((size_t) m / (LEAF / 2) * 2 + 4) * sizeof(Node));
  qIds  = (int *) malloc (k * sizeof(int));
  qDist = (double *) malloc (k * sizeof(double));
  if (!items || !nodes || !qIds || !qDist) fatal ("Out of memory building kd-tree");

  for (i = 0; i < m; ++i) {
    c = cities[i];
    items[i].x[0] = xx[c];
    items[i].x[1] = yy[c];
    items[i].x[2] = zz[c];
    items[i].idx = i;
  }
  nNodes = 0;
  build (0, m);

  qK = k;
  for (j = 0; j < m; ++j) {                         /* queries in tree order */
    qx = items[j].x;
    qIdx = items[j].idx;
    qCount = 0;
    search (0);
    for (i = 0; i < k; ++i) out[(long long) qIdx * k + i] = cities[qIds[i]];
  }

  free (items);
  free (nodes);
  free (qIds);
  free (qDist);
}

void buildNeighbors (int k) {

  int *all, c;

  nbK = (k < ncities - 1) ? k : ncities - 1;        /* can't have more than N-1 */
  nbList = (int *) malloc ((long long)(ncities + 1) * nbK * sizeof(int));
  all = (int *) malloc (ncities * sizeof(int));
  if (!nbList || !all) fatal ("Out of memory allocating neighbour lists");
  for (c = 1; c <= ncities; ++c) all[c-1] = c;
  nearestAmong (all, ncities, nbK, nbList + nbK);   /* city c's list starts at */
  free (all);                                       /*   nbList[c*nbK]         */
}
