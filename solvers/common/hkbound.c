#include <math.h>
#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <time.h>

#include "Instance.h"
#include "Utilities.h"
#include "tour.h"
#include "neighbors.h"
#include "greedystart.h"
#include "localsearch.h"

/*********************************************************************

Our own Held-Karp 1-tree lower bound (any dimension).

	hkbound instance_file [--period P] [--step S] [--k K] [--rounds R]

Give every city v a number pi[v] and use the edge costs c(u,v) + pi[u] +
pi[v].  Every tour gets longer by exactly 2 * sum(pi), so

	W(pi) = (cheapest 1-tree under these costs) - 2 * sum(pi)

is a lower bound on the shortest tour, for ANY pi.  Subgradient ascent
(Held & Karp 1971) looks for a good pi; the best value is the Held-Karp
bound, and what this program reports is at most that.

A minimum 1-tree here = a minimum spanning tree plus the second-cheapest
edge at one of its leaves, the leaf where that edge is longest (leaf v's
tree edge is its cheapest edge, so this is the cheapest 1-tree rooted at v).

Two steps:
  1. ascent on a sparse graph (K nearest neighbours + a greedy TO tour):
     fast, but a sparse spanning tree can miss cheaper edges, so its value
     is NOT a proven bound;
  2. the final pi is evaluated on the COMPLETE graph (Prim, O(N^2)).  That
     value is a proven bound, and is the one reported as hk_bound.
Edges the complete tree uses that the sparse graph lacks are added, and
both steps repeat (--rounds) so the two values meet.

Distances are the bench's rounded integers; pi are doubles, and the result
is printed rounded down to one decimal (the error of the sums is far
smaller than that).

 *********************************************************************/

typedef struct { int u, v; } E;

static E     *edges;             /* sparse graph                              */
static int    nE, capE;
static int   *adjStart, *adj;    /* adjacency: adj[adjStart[v]..] = neighbours */
static int   *adjLen;            /*   and the edge lengths                    */
static double *pi, *bestPi;
static int   *deg, *oldG;

static int edgeCmp (const void *a, const void *b) {
  const E *p = (const E *) a, *q = (const E *) b;
  if (p->u != q->u) return p->u < q->u ? -1 : 1;
  return p->v < q->v ? -1 : p->v > q->v;
}

static void addEdge (int a, int b) {
  if (nE == capE) {
    capE = capE ? capE * 2 : 1024;
    edges = (E *) realloc (edges, capE * sizeof(E));
    if (!edges) fatal ("Out of memory");
  }
  edges[nE].u = a < b ? a : b;
  edges[nE].v = a < b ? b : a;
  nE++;
}

static void buildAdjacency (void) {
  int i, j, n = 0, *fill;
  qsort (edges, nE, sizeof(E), edgeCmp);
  for (i = 0; i < nE; ++i)                       /* drop duplicates */
    if (i == 0 || edges[i].u != edges[i-1].u || edges[i].v != edges[i-1].v) edges[n++] = edges[i];
  nE = n;
  free (adjStart); free (adj); free (adjLen);
  adjStart = (int *) calloc (ncities + 2, sizeof(int));
  adj    = (int *) malloc (2 * (size_t) nE * sizeof(int));
  adjLen = (int *) malloc (2 * (size_t) nE * sizeof(int));
  fill = (int *) malloc ((ncities + 2) * sizeof(int));
  if (!adjStart || !adj || !adjLen || !fill) fatal ("Out of memory");
  for (i = 0; i < nE; ++i) { adjStart[edges[i].u + 1]++; adjStart[edges[i].v + 1]++; }
  for (i = 1; i <= ncities + 1; ++i) adjStart[i] += adjStart[i-1];
  memcpy (fill, adjStart, (ncities + 2) * sizeof(int));
  for (i = 0; i < nE; ++i) {
    j = distance (edges[i].u, edges[i].v);
    adj[fill[edges[i].u]] = edges[i].v;  adjLen[fill[edges[i].u]++] = j;
    adj[fill[edges[i].v]] = edges[i].u;  adjLen[fill[edges[i].v]++] = j;
  }
  free (fill);
}

/* ---------------------------------------------------------------- sparse 1-tree */

static double *key;              /* Prim: cheapest link to the tree           */
static int    *par;              /*       and where it leads                  */
static int    *heap, *hpos, hn;  /* binary heap of cities by key              */

static void siftUp (int i) {
  int c = heap[i], p;
  while (i > 0 && key[heap[p = (i - 1) / 2]] > key[c]) { heap[i] = heap[p]; hpos[heap[i]] = i; i = p; }
  heap[i] = c;
  hpos[c] = i;
}

static int popMin (void) {
  int top = heap[0], c = heap[--hn], i = 0, ch;
  while ((ch = 2 * i + 1) < hn) {
    if (ch + 1 < hn && key[heap[ch + 1]] < key[heap[ch]]) ch++;
    if (key[heap[ch]] >= key[c]) break;
    heap[i] = heap[ch];
    hpos[heap[i]] = i;
    i = ch;
  }
  if (hn > 0) { heap[i] = c; hpos[c] = i; }
  hpos[top] = -2;
  return top;
}

/* Value of the sparse 1-tree under pi; fills deg[].  *leaf = its special city. */
static double sparseTree (int *leaf) {
  int c, v, a, n = 0;
  double total = 0, w, second, bestSecond = -1;
  int bestV = 0, bestW = 0, secondW;

  for (c = 1; c <= ncities; ++c) { hpos[c] = -1; deg[c] = 0; }
  hn = 0;
  key[1] = 0; par[1] = 0; heap[0] = 1; hpos[1] = 0; hn = 1;
  while (hn > 0) {
    c = popMin();
    n++;
    if (par[c]) { total += key[c]; deg[c]++; deg[par[c]]++; }
    for (a = adjStart[c]; a < adjStart[c + 1]; ++a) {
      v = adj[a];
      if (hpos[v] == -2) continue;               /* already in the tree */
      w = adjLen[a] + pi[c] + pi[v];
      if (hpos[v] == -1) { key[v] = w; par[v] = c; heap[hn] = v; siftUp (hn++); }
      else if (w < key[v]) { key[v] = w; par[v] = c; siftUp (hpos[v]); }
    }
  }
  if (n < ncities) fatal ("sparse graph is not connected");

  for (c = 1; c <= ncities; ++c) {               /* second-cheapest edge at leaves */
    if (deg[c] != 1) continue;
    double m1 = 1e300, m2 = 1e300;
    secondW = 0;
    int firstW = 0;
    for (a = adjStart[c]; a < adjStart[c + 1]; ++a) {
      w = adjLen[a] + pi[c] + pi[adj[a]];
      if (w < m1) { m2 = m1; secondW = firstW; m1 = w; firstW = adj[a]; }
      else if (w < m2) { m2 = w; secondW = adj[a]; }
    }
    second = m2;
    if (second < 1e299 && (bestV == 0 || second > bestSecond)) { bestSecond = second; bestV = c; bestW = secondW; }
  }
  if (bestV) { total += bestSecond; deg[bestV]++; deg[bestW]++; }
  *leaf = bestV;
  for (c = 1; c <= ncities; ++c) total -= 2 * pi[c];
  return total;
}

/* ---------------------------------------------------------------- exact 1-tree */

#define TRY 64                   /* leaves whose exact second edge is computed */

/* Value of the minimum 1-tree on the complete graph under pi: a proven lower
   bound.  Adds its edges to the sparse graph; *missing = how many were new. */
static double exactTree (int *missing) {

  int *rest, m, i, j, c, v, best, *d = (int *) calloc (ncities + 1, sizeof(int));
  double *k = (double *) malloc ((ncities + 1) * sizeof(double));
  int *p = (int *) malloc ((ncities + 1) * sizeof(int)), a, found;
  double total = 0, w, bw, m1, m2, secondBest = -1e300;
  int tryV[TRY], nTry = 0;
  double tryS[TRY];

  rest = (int *) malloc (ncities * sizeof(int));
  if (!d || !k || !p || !rest) fatal ("Out of memory");
  m = 0;
  for (c = 2; c <= ncities; ++c) { rest[m++] = c; k[c] = 1e300; p[c] = 1; }
  c = 1;
  *missing = 0;
  while (m > 0) {                                /* Prim: c just joined the tree */
    best = 0;
    bw = 1e301;
    for (i = 0; i < m; ++i) {
      v = rest[i];
      w = distance (c, v) + pi[c] + pi[v];
      if (w < k[v]) { k[v] = w; p[v] = c; }
      if (k[v] < bw) { bw = k[v]; best = i; }
    }
    c = rest[best];
    rest[best] = rest[--m];
    total += k[c];
    d[c]++;
    d[p[c]]++;
    found = 0;                                   /* in the sparse graph? */
    for (a = adjStart[c]; a < adjStart[c + 1]; ++a) if (adj[a] == p[c]) { found = 1; break; }
    if (!found) { addEdge (c, p[c]); (*missing)++; }
  }

  /* the leaves with the longest sparse second edge: compute theirs exactly
     (any leaf gives a valid bound, so a short list is enough) */
  for (v = 1; v <= ncities; ++v) {
    if (d[v] != 1) continue;
    m1 = m2 = 1e300;
    for (a = adjStart[v]; a < adjStart[v + 1]; ++a) {
      w = adjLen[a] + pi[v] + pi[adj[a]];
      if (w < m1) { m2 = m1; m1 = w; } else if (w < m2) m2 = w;
    }
    for (i = nTry; i > 0 && tryS[i-1] < m2; --i)
      if (i < TRY) { tryS[i] = tryS[i-1]; tryV[i] = tryV[i-1]; }
    if (i < TRY) { tryS[i] = m2; tryV[i] = v; if (nTry < TRY) nTry++; }
  }
  for (i = 0; i < nTry; ++i) {
    v = tryV[i];
    m1 = m2 = 1e300;
    for (j = 1; j <= ncities; ++j) {
      if (j == v) continue;
      w = distance (v, j) + pi[v] + pi[j];
      if (w < m1) { m2 = m1; m1 = w; } else if (w < m2) m2 = w;
    }
    if (m2 > secondBest) secondBest = m2;
  }
  total += secondBest;
  for (c = 1; c <= ncities; ++c) total -= 2 * pi[c];
  free (rest); free (d); free (k); free (p);
  return total;
}

/* ---------------------------------------------------------------- main */

int main (int argc, char **argv) {

  int i, k = 10, period0 = 0, rounds = 8, r, c, leaf, missing = 1, x, p, period, initial, total;
  long long m2, mo;
  double ub, w, best, exact, bestExact = 0, t, t0, avg, step0 = 0.5;
  clock_t clk = clock();

  if (argc < 2) fatal ("Usage: hkbound instance_file [--period P] [--step S] [--k K] [--rounds R]");
  for (i = 2; i + 1 < argc; i += 2) {
    if      (!strcmp (argv[i], "--period")) period0 = atoi (argv[i+1]);
    else if (!strcmp (argv[i], "--step"))   step0 = atof (argv[i+1]);
    else if (!strcmp (argv[i], "--k"))      k = atoi (argv[i+1]);
    else if (!strcmp (argv[i], "--rounds")) rounds = atoi (argv[i+1]);
    else fatal ("unknown option");
  }
  readInstance (argv[1]);
  tourAlloc();
  buildNeighbors (16);
  greedyTour (10);
  localSearch (1, &m2, &mo);
  ub = tourLength();
  avg = ub / ncities;                            /* mean tour edge */
  if (period0 == 0) period0 = ncities / 2 < 100 ? 100 : ncities / 2 > 2000 ? 2000 : ncities / 2;
  if (k > nbK) k = nbK;
  for (c = 1; c <= ncities; ++c)
    for (i = 0; i < k; ++i) addEdge (c, neighbors (c)[i]);
  for (i = 0; i < ncities; ++i) addEdge (tour[i], tour[(i + 1) % ncities]);
  buildAdjacency();

  pi     = (double *) calloc (ncities + 1, sizeof(double));
  bestPi = (double *) calloc (ncities + 1, sizeof(double));
  key    = (double *) malloc ((ncities + 1) * sizeof(double));
  par    = (int *) calloc (ncities + 1, sizeof(int));
  heap   = (int *) malloc ((ncities + 1) * sizeof(int));
  hpos   = (int *) malloc ((ncities + 1) * sizeof(int));
  deg    = (int *) malloc ((ncities + 1) * sizeof(int));
  oldG   = (int *) calloc (ncities + 1, sizeof(int));
  if (!pi || !bestPi || !key || !par || !heap || !hpos || !deg || !oldG) fatal ("Out of memory");
  printf ("n=%d\ntour=%.0f\n", ncities, ub);

  for (r = 0; r < rounds; ++r) {
    best = -1e300;
    memcpy (pi, bestPi, (ncities + 1) * sizeof(double));
    memset (oldG, 0, (ncities + 1) * sizeof(int));
    total = 0;
    t0 = avg * step0;
    initial = 1;
    /* step sizes as in LKH (Helsgaun 2000): double the step while that keeps
       helping, then halve the step and the length of the phase in turn     */
    for (period = period0, t = t0; period > 0 && t > avg * 1e-6; period /= 2, t /= 2) {
      for (p = 1; p <= period; ++p, ++total) {
	w = sparseTree (&leaf);
	if (w > best) {
	  best = w;
	  memcpy (bestPi, pi, (ncities + 1) * sizeof(double));
	  if (initial && total > 0) t *= 2;
	  if (p == period && period * 2 <= period0) period *= 2;
	}
	else if (initial && p > period / 2) {      /* stopped helping: end of phase 1 */
	  initial = 0;
	  p = 0;
	  t = 3 * t / 4;
	}
	for (c = 1; c <= ncities; ++c) {           /* step, with some momentum */
	  x = deg[c] - 2;
	  pi[c] += t * (0.7 * x + 0.3 * oldG[c]);
	  oldG[c] = x;
	}
      }
    }
    memcpy (pi, bestPi, (ncities + 1) * sizeof(double));
    exact = exactTree (&missing);
    if (exact > bestExact) bestExact = exact;
    fprintf (stderr, "round %d: %d iterations, sparse %.1f, exact %.1f, %d missing edges, %.0f s\n",
	     r + 1, total, best, exact, missing, (double)(clock() - clk) / CLOCKS_PER_SEC);
    printf ("sparse_value_%d=%.1f\nexact_value_%d=%.1f\nmissing_edges_%d=%d\n",
	    r + 1, best, r + 1, exact, r + 1, missing);
    if (missing == 0 && r > 0) break;
    buildAdjacency();
  }
  printf ("hk_bound=%.1f\ntour_gap_pct=%.3f\ntime=%.1f\n", floor (bestExact * 10) / 10,
	  (ub / bestExact - 1) * 100, (double)(clock() - clk) / CLOCKS_PER_SEC);
  return 0;
}
