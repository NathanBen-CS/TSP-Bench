#include <stdlib.h>

#include "Instance.h"
#include "Utilities.h"
#include "tour.h"
#include "neighbors.h"
#include "greedystart.h"

/*********************************************************************

Greedy (multi-fragment) start tour, built from candidate edges.

The idea is the same as the DIMACS greedy code (codes/gmain.c):
  look at edges from shortest to longest, and keep an edge whenever
    - neither end city already has 2 tour edges, and
    - it doesn't close a loop (its ends are in different pieces).
That grows many path "fragments" that merge into longer ones.

Instead of searching for true nearest neighbours like codes/gmain.c, we
work in rounds over short candidate lists:

  round 1   candidates = each city's nearest neighbours (neighbor lists)
  round 2+  candidates = each fragment END's nearest other fragment ends
            (geometric instances; recomputed each round with the grid)
  finally   any fragments still left are chained nearest-end-first

Each round merges most remaining fragments, so only a handful reach the
final chaining step, which is a linear scan per fragment.

Why bother: 2-opt is very sensitive to its start tour.  From a greedy start
it ends ~5% above the Held-Karp bound on random uniform instances; from a
space-filling-curve start it ends far worse (~18%), because those tours
contain long "seams" that 2-opt with neighbour lists can't repair.

 *********************************************************************/

#define JOIN_K 8                /* candidate partners per end in join rounds */

typedef struct {
  int len;                      /* edge length          */
  int u, v;                     /* end cities, u < v    */
} Edge;

static int *degree;             /* number of tour edges at each city (0..2) */
static int *adj;                /* adj[2c], adj[2c+1]: c's tour partners    */
static int *parent;             /* union-find: which fragment is c in?      */
static int  fragments;          /* how many fragments exist right now       */

/*******************************************************************/

static int find (int c) {               /* root of c's fragment */
  while (parent[c] != c) {
    parent[c] = parent[parent[c]];      /* path halving: shortcut as we go */
    c = parent[c];
  }
  return c;
}

static int byLength (const void *a, const void *b) {   /* qsort comparator */

  const Edge *p = (const Edge *) a;
  const Edge *q = (const Edge *) b;

  if (p->len != q->len) return (p->len < q->len) ? -1 : 1;  /* shortest first */
  if (p->u   != q->u)   return p->u - q->u;                  /* then fixed ties */
  return p->v - q->v;
}

static void setEdge (Edge *e, int a, int b) {   /* store edge (a,b) with u < v */
  e->len = distance (a, b);
  e->u = (a < b) ? a : b;
  e->v = (a < b) ? b : a;
}

/*******************************************************************/

/* Sort the candidate edges and keep each one that is still legal.
   Returns how many were added to the tour. */

static int addGreedily (Edge *edges, int nEdges) {

  int i, u, v, added = 0;

  qsort (edges, nEdges, sizeof(Edge), byLength);   /* shortest first */

  for (i = 0; i < nEdges; ++i) {
    u = edges[i].u;
    v = edges[i].v;
    if (degree[u] == 2 || degree[v] == 2) continue;  /* a city would get 3 edges */
    if (find (u) == find (v)) continue;              /* would close a loop early */
    adj[2*u + degree[u]++] = v;                      /* record the edge at u ... */
    adj[2*v + degree[v]++] = u;                      /* ... and at v             */
    parent[find (u)] = find (v);                     /* the two pieces merge     */
    fragments--;
    added++;
  }
  return added;
}

/*******************************************************************/

/* Round 1: each city's k nearest neighbours as candidate edges. */

static void neighborRound (int k) {

  int c, i, j, n, nEdges = 0, dup;
  int *nb;
  Edge *edges = (Edge *) malloc ((long long) ncities * k * sizeof(Edge));

  if (edges == 0) fatal ("Out of memory in greedyTour");

  for (c = 1; c <= ncities; ++c) {
    for (i = 0; i < k; ++i) {
      n = neighbors(c)[i];
      if (n < c) {                              /* (n,c) was already added    */
	nb = neighbors (n);                     /*   from n's list if c is in */
	for (dup = 0, j = 0; j < k; ++j)        /*   n's first k neighbours   */
	  if (nb[j] == c) { dup = 1; break; }
	if (dup) continue;
      }
      setEdge (&edges[nEdges++], c, n);
    }
  }

  addGreedily (edges, nEdges);
  free (edges);
}

/*******************************************************************/

/* Rounds 2+: join fragment ends to their nearest other fragment ends.
   Repeats until one fragment is left or a round adds nothing. */

static void joinRounds (void) {

  int *ends, *cand;
  int c, i, j, m, k, nEdges;
  Edge *edges;

  ends = (int *) malloc (ncities * sizeof(int));
  while (fragments > 1) {

    m = 0;                                       /* collect current ends: */
    for (c = 1; c <= ncities; ++c)               /*   cities with < 2 edges */
      if (degree[c] < 2) ends[m++] = c;
    if (m < 3) break;

    k = (JOIN_K < m - 1) ? JOIN_K : m - 1;
    cand = (int *) malloc ((long long) m * k * sizeof(int));
    edges = (Edge *) malloc ((long long) m * k * sizeof(Edge));
    if (cand == 0 || edges == 0) fatal ("Out of memory in greedyTour");

    nearestAmong (ends, m, k, cand);             /* each end's k nearest ends */

    nEdges = 0;
    for (i = 0; i < m; ++i)                      /* every (end, near end) pair */
      for (j = 0; j < k; ++j)                    /*   is a candidate edge;     */
	setEdge (&edges[nEdges++], ends[i], cand[i*k + j]);  /* duplicates are */
                                                 /* rejected by the loop check */
    i = addGreedily (edges, nEdges);
    free (cand);
    free (edges);
    if (i == 0) break;                           /* no progress: stop rounds */
  }
  free (ends);
}

/*******************************************************************/

/* Final step: walk the remaining fragments into tour[], always continuing
   with the fragment whose end is closest to where we are. */

static void chainFragments (void) {

  int *otherEnd, *ends, *endIndex;
  int c, i, j, nEnds, e, f, prev, cur, next, best, bestd, d;

  /* each fragment's two ends: walk from one end to the other
     (a city with 0 edges is a fragment on its own, both of its ends) */
  otherEnd = (int *) malloc ((ncities + 1) * sizeof(int));
  for (c = 1; c <= ncities; ++c) {
    if (degree[c] == 2) continue;                /* middle city: not an end */
    if (degree[c] == 0) { otherEnd[c] = c; continue; }
    prev = c;
    cur = adj[2*c];
    while (degree[cur] == 2) {                   /* follow the path ...     */
      next = (adj[2*cur] == prev) ? adj[2*cur + 1] : adj[2*cur];
      prev = cur;
      cur = next;
    }
    otherEnd[c] = cur;                           /* ... to the far end      */
  }

  /* to-do list of ends; endIndex[c] = c's slot in ends[] */
  ends = (int *) malloc (ncities * sizeof(int));
  endIndex = (int *) malloc ((ncities + 1) * sizeof(int));
  nEnds = 0;
  for (c = 1; c <= ncities; ++c)
    if (degree[c] < 2) { endIndex[c] = nEnds; ends[nEnds++] = c; }

  j = 0;                                         /* next free tour position */
  e = ends[0];                                   /* start at the first end  */
  while (1) {
    f = otherEnd[e];

    /* take both ends of this fragment off the to-do list (swap-remove) */
    i = endIndex[e];  ends[i] = ends[--nEnds];  endIndex[ends[i]] = i;
    if (f != e) {
      i = endIndex[f];  ends[i] = ends[--nEnds];  endIndex[ends[i]] = i;
    }

    /* copy the fragment's cities into the tour, walking from e to f */
    prev = 0;
    cur = e;
    while (1) {
      tour[j++] = cur;
      if (cur == f) break;
      next = (adj[2*cur] != prev) ? adj[2*cur] : adj[2*cur + 1];  /* the partner */
      prev = cur;                                       /* we didn't come from */
      cur = next;
    }

    if (nEnds == 0) break;                       /* every fragment placed */

    best = ends[0];                              /* closest remaining end to f */
    bestd = distance (f, best);
    for (i = 1; i < nEnds; ++i) {
      d = distance (f, ends[i]);
      if (d < bestd) { bestd = d; best = ends[i]; }
    }
    e = best;                                    /* continue from there */
  }

  if (j != ncities) fatal ("greedyTour: internal error, tour incomplete");

  free (otherEnd);
  free (ends);
  free (endIndex);
}

/*******************************************************************/

void greedyTour (int candidatesPerCity) {

  int c, k = (candidatesPerCity < nbK) ? candidatesPerCity : nbK;

  degree = (int *) calloc (ncities + 1, sizeof(int));
  adj    = (int *) calloc (2 * (ncities + 1), sizeof(int));
  parent = (int *) malloc ((ncities + 1) * sizeof(int));
  if (!degree || !adj || !parent) fatal ("Out of memory in greedyTour");
  for (c = 1; c <= ncities; ++c) parent[c] = c;  /* every city on its own */
  fragments = ncities;

  neighborRound (k);                   /* 1. nearest-neighbour candidates  */
  if (isEuclidean) joinRounds();       /* 2. join fragment ends, in rounds */
  chainFragments();                    /* 3. chain whatever is left        */
  tourSetPositions();

  free (degree);
  free (adj);
  free (parent);
}
