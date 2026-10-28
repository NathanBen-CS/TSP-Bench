#include <stdlib.h>
#include <string.h>

#include "Instance.h"
#include "Utilities.h"
#include "tour.h"
#include "neighbors.h"
#include "boxes.h"
#include "merges.h"

/*********************************************************************

Merging box loops into one tour.  See merges.h for the five methods.

The basic move (used by M3, M4, M5) joins two loops into one:

	remove edge {x,y} from loop A and edge {p,q} from loop B, then add
	either {x,q} + {y,p}   (option A)   or   {x,p} + {y,q}   (option B)

	cost = new edges - removed edges
	     = d(x,e1) + d(y,e2) - d(x,y) - d(p,q),   {e1,e2} = {q,p} or {p,q}

Loops are stored as cycles with next/prev links (nxt[], prv[]) so a merge
just relinks four cities.  A loop that runs the wrong way for the chosen
option is reversed first: always the smaller of the two, so reversals cost
O(N log N) in total.  A union-find structure says which loop a city is in.

 *********************************************************************/

const char *mergeName[NMERGES] = { "M1", "M2", "M3", "M4", "M5" };

int mergeByName (const char *name) {
  int m;
  for (m = 0; m < NMERGES; ++m)
    if (!strcmp (name, mergeName[m])) return m;
  return -1;
}

int mergeNeedsNeighbors (int merge) {
  return merge == MERGE_M3 || merge == MERGE_M4 || merge == MERGE_M5;
}

#define D(a,b) ((long long) distance ((a), (b)))
#define BIG 0x3fffffffffffffffLL

/* ================================================================ M1: snake glue

Walk the boxes in curve order.  Each loop is opened at one edge (u, v):
enter at u from the previous box, go round, leave from v.  Each box picks
the u and direction minimising
	d(prev, u) - d(u, v) + d(v, first city of the next box)
on its own (the first box is opened at a fixed spot).  This is the
original glue from boxes.c, unchanged.                                  */

static long long mergeM1 (const BoxSet *bs) {

  int j, lo, s, r, dir, bu, bdir, u, v, prev, nextRef, tn, k;
  const int *g;
  long long cost, bestCost, added = 0;

  tn = 0;
  prev = -1;                                       /* no previous box yet */
  for (j = 0; j < bs->nb; ++j) {
    lo = bs->start[j];
    s = bs->start[j + 1] - lo;
    g = bs->loops + lo;                            /* this box's loop     */
    nextRef = (j + 1 < bs->nb) ? bs->order[bs->start[j + 1]] : -1;
    bu = 0;
    bdir = 1;
    if (s > 1 && prev >= 0) {
      bestCost = 0;
      for (r = 0; r < s; ++r)                      /* every entry city u  */
	for (dir = -1; dir <= 1; dir += 2) {       /* and both directions */
	  u = g[r];
	  v = g[((r - dir) % s + s) % s];          /* exit = u's other loop neighbour */
	  cost = D(prev, u) - D(u, v)
	       + (nextRef >= 0 ? D(v, nextRef) : D(v, tour[0]));
	  if ((r == 0 && dir == -1) || cost < bestCost) {   /* first option, or */
	    bestCost = cost;                                /*   a cheaper one  */
	    bu = r;
	    bdir = dir;
	  }
	}
    }
    if (s > 1) added -= D(g[bu], g[((bu - bdir) % s + s) % s]);  /* edge we cut */
    if (prev >= 0) added += D(prev, g[bu]);        /* edge into this box  */
    for (k = 0; k < s; ++k)                        /* round the loop from u  */
      tour[tn++] = g[((bu + bdir * k) % s + s) % s];
    prev = tour[tn - 1];
  }
  added += D(prev, tour[0]);                       /* close the tour */
  return added;
}

/* ================================================================ M2: snake, cuts by DP

Same chain of boxes as M1, but every box's cut is chosen jointly.  A
"choice" for a box is (entry u, direction): it fixes the exit v and the
cut edge (u, v).  The snake's extra length is

	sum over boxes of  -d(u_j, v_j)  +  sum over joins of  d(v_j, u_j+1)
	+ the closing join d(v_last, u_0)

Dynamic programming along the chain: f(choice of box j) = best value of
boxes 0..j ending with that choice.  The closing join ties the last box
to the first, so the DP runs once for every choice of box 0 and keeps the
best.  That makes M2 the best possible snake for these loops and order.  */

static long long mergeM2 (const BoxSet *bs) {

  int nb = bs->nb, j, lo, s, r, dir, c0, ch, pc, k, tn, nc, best0 = 0, last = 0;
  int *cStart, *cEntry, *cExit, *cR, *cDir, *back;
  long long *cCut, *f, val, bestVal = BIG, total, bv;
  const int *g;

  /* ---- 1. list every box's choices ---- */
  cStart = (int *) malloc ((nb + 1) * sizeof(int));
  nc = 0;
  for (j = 0; j < nb; ++j) {
    cStart[j] = nc;
    s = bs->start[j + 1] - bs->start[j];
    nc += (s == 1) ? 1 : 2 * s;
  }
  cStart[nb] = nc;
  cEntry = (int *) malloc (nc * sizeof(int));
  cExit  = (int *) malloc (nc * sizeof(int));
  cR     = (int *) malloc (nc * sizeof(int));
  cDir   = (int *) malloc (nc * sizeof(int));
  cCut   = (long long *) malloc (nc * sizeof(long long));
  f      = (long long *) malloc (nc * sizeof(long long));
  back   = (int *) malloc (nc * sizeof(int));
  if (!cEntry || !cExit || !cR || !cDir || !cCut || !f || !back) fatal ("Out of memory in M2");

  for (j = 0; j < nb; ++j) {
    lo = bs->start[j];
    s = bs->start[j + 1] - lo;
    g = bs->loops + lo;
    ch = cStart[j];
    if (s == 1) {                                  /* one city: enter = leave */
      cEntry[ch] = cExit[ch] = g[0]; cR[ch] = 0; cDir[ch] = 1; cCut[ch] = 0;
      continue;
    }
    for (r = 0; r < s; ++r)
      for (dir = -1; dir <= 1; dir += 2, ++ch) {
	cR[ch] = r;
	cDir[ch] = dir;
	cEntry[ch] = g[r];
	cExit[ch] = g[((r - dir) % s + s) % s];
	cCut[ch] = D(cEntry[ch], cExit[ch]);
      }
  }

  /* ---- 2. DP for every choice of box 0; remember the best ---- */
  for (c0 = cStart[0]; c0 < cStart[1]; ++c0) {
    if (nb == 1) { val = -cCut[c0] + D(cExit[c0], cEntry[c0]); last = c0; }
    else {
      for (ch = cStart[1]; ch < cStart[2]; ++ch)   /* box 1 follows the fixed box 0 */
	f[ch] = -cCut[c0] + D(cExit[c0], cEntry[ch]) - cCut[ch];
      for (j = 2; j < nb; ++j)
	for (ch = cStart[j]; ch < cStart[j + 1]; ++ch) {
	  bv = BIG;
	  for (pc = cStart[j - 1]; pc < cStart[j]; ++pc) {
	    long long v2 = f[pc] + D(cExit[pc], cEntry[ch]);
	    if (v2 < bv) bv = v2;
	  }
	  f[ch] = bv - cCut[ch];
	}
      val = BIG;
      for (ch = cStart[nb - 1]; ch < cStart[nb]; ++ch) {  /* close back to box 0 */
	long long v2 = f[ch] + D(cExit[ch], cEntry[c0]);
	if (v2 < val) val = v2;
      }
    }
    if (val < bestVal) { bestVal = val; best0 = c0; }
  }

  /* ---- 3. rerun with the best box-0 choice, keeping back pointers ---- */
  if (nb > 1) {
    for (ch = cStart[1]; ch < cStart[2]; ++ch) {
      f[ch] = -cCut[best0] + D(cExit[best0], cEntry[ch]) - cCut[ch];
      back[ch] = best0;
    }
    for (j = 2; j < nb; ++j)
      for (ch = cStart[j]; ch < cStart[j + 1]; ++ch) {
	bv = BIG;
	for (pc = cStart[j - 1]; pc < cStart[j]; ++pc) {
	  long long v2 = f[pc] + D(cExit[pc], cEntry[ch]);
	  if (v2 < bv) { bv = v2; back[ch] = pc; }
	}
	f[ch] = bv - cCut[ch];
      }
    total = BIG;
    for (ch = cStart[nb - 1]; ch < cStart[nb]; ++ch) {
      long long v2 = f[ch] + D(cExit[ch], cEntry[best0]);
      if (v2 < total) { total = v2; last = ch; }
    }
  }
  else last = best0;

  /* ---- 4. walk the choices back from the last box, then write the tour ---- */
  {
    int *pick = (int *) malloc (nb * sizeof(int));
    ch = last;
    for (j = nb - 1; j >= 1; --j) { pick[j] = ch; ch = back[ch]; }
    pick[0] = best0;
    tn = 0;
    for (j = 0; j < nb; ++j) {
      lo = bs->start[j];
      s = bs->start[j + 1] - lo;
      g = bs->loops + lo;
      ch = pick[j];
      for (k = 0; k < s; ++k)                      /* from u round to v */
	tour[tn++] = g[((cR[ch] + cDir[ch] * k) % s + s) % s];
    }
    free (pick);
  }

  free (cStart); free (cEntry); free (cExit); free (cR); free (cDir);
  free (cCut); free (f); free (back);
  return bestVal;
}

/* ================================================================ loops as cycles (M3-M5) */

static int *nxt, *prv;            /* cycle links                          */
static int *uf, *ufSize;          /* union-find over cities; loop sizes   */
static int  loopsLeft;

static int find (int c) {
  while (uf[c] != c) { uf[c] = uf[uf[c]]; c = uf[c]; }   /* path halving */
  return c;
}

static void reverseCycle (int c) {                 /* flip a loop's direction */
  int x = c, t;
  do { t = nxt[x]; nxt[x] = prv[x]; prv[x] = t; x = t; } while (x != c);
}

static void initCycles (const BoxSet *bs) {
  int j, lo, s, k, a, b;
  nxt = (int *) malloc ((ncities + 1) * sizeof(int));
  prv = (int *) malloc ((ncities + 1) * sizeof(int));
  uf  = (int *) malloc ((ncities + 1) * sizeof(int));
  ufSize = (int *) malloc ((ncities + 1) * sizeof(int));
  if (!nxt || !prv || !uf || !ufSize) fatal ("Out of memory in merge");
  for (j = 0; j < bs->nb; ++j) {                   /* each box loop: one cycle */
    lo = bs->start[j];
    s = bs->start[j + 1] - lo;
    for (k = 0; k < s; ++k) {
      a = bs->loops[lo + k];
      b = bs->loops[lo + (k + 1) % s];
      nxt[a] = b;
      prv[b] = a;
      uf[a] = bs->loops[lo];                       /* first city is the root */
    }
    ufSize[bs->loops[lo]] = s;
  }
  loopsLeft = bs->nb;
}

static void freeCycles (void) {
  free (nxt); free (prv); free (uf); free (ufSize);
}

/* Best way to merge using edge {x,y} of one loop and {p,q} of another.
   Returns the cost; *e1, *e2 receive the partners of x and y.          */
static long long mergeCost (int x, int y, int p, int q, int *e1, int *e2) {
  long long base = D(x, y) + D(p, q);
  long long a = D(x, q) + D(y, p) - base;          /* option A: {x,q} {y,p} */
  long long b = D(x, p) + D(y, q) - base;          /* option B: {x,p} {y,q} */
  if (a <= b) { *e1 = q; *e2 = p; return a; }
  *e1 = p; *e2 = q; return b;
}

/* Apply a merge: remove {x,y} and {e1,e2}, add {x,e1} and {y,e2}. */
static void mergeCycles (int x, int y, int e1, int e2) {
  int t, rx, re;
  if (nxt[x] != y) { t = x; x = y; y = t; t = e1; e1 = e2; e2 = t; }  /* x -> y   */
  if (nxt[e2] != e1) {                             /* other loop runs e1 -> e2: */
    rx = find (x); re = find (e1);                 /*   reverse the smaller one */
    if (ufSize[re] <= ufSize[rx]) reverseCycle (e1);
    else {
      reverseCycle (x);                            /* now y -> x: swap names,   */
      t = x; x = y; y = t; t = e1; e1 = e2; e2 = t;/*   same new edges          */
    }
  }
  /* now x -> y in one loop and e2 -> e1 in the other: x -> e1 ... e2 -> y */
  nxt[x] = e1;  prv[e1] = x;
  nxt[e2] = y;  prv[y] = e2;
  rx = find (x); re = find (e1);                   /* the two loops are one now */
  if (ufSize[rx] < ufSize[re]) { t = rx; rx = re; re = t; }
  uf[re] = rx;
  ufSize[rx] += ufSize[re];
  loopsLeft--;
}

/* Best merge between the loop of box A's cities and the loop of box B's
   cities, trying every current edge at those cities.  Returns the cost. */
static long long mergeBoxes (const BoxSet *bs, int boxA, int boxB) {
  int ia, ib, ea, eb, x, y, p, q, e1, e2;
  int bx = 0, by = 0, b1 = 0, b2 = 0;
  long long c, bestC = BIG;
  for (ia = bs->start[boxA]; ia < bs->start[boxA + 1]; ++ia)
    for (ea = 0; ea < 2; ++ea) {
      x = bs->loops[ia];
      y = ea ? prv[x] : nxt[x];                    /* both edges at x */
      for (ib = bs->start[boxB]; ib < bs->start[boxB + 1]; ++ib)
	for (eb = 0; eb < 2; ++eb) {
	  p = bs->loops[ib];
	  q = eb ? prv[p] : nxt[p];
	  c = mergeCost (x, y, p, q, &e1, &e2);
	  if (c < bestC) { bestC = c; bx = x; by = y; b1 = e1; b2 = e2; }
	}
    }
  mergeCycles (bx, by, b1, b2);
  return bestC;
}

static void cyclesToTour (int startCity) {         /* read the one loop off */
  int i, c = startCity;
  for (i = 0; i < ncities; ++i) { tour[i] = c; c = nxt[c]; }
}

/* ================================================================ M3: best pair in curve order */

static long long mergeM3 (const BoxSet *bs) {

  int j, lo, s, i, k, c, nc, ea, x, y, ib, p, q, e1, e2;
  int bx = 0, by = 0, b1 = 0, b2 = 0, stamp = 0;
  int *inR, *mark, *cand;
  long long cost, bestC, added = 0;

  initCycles (bs);
  inR  = (int *) calloc (ncities + 1, sizeof(int));
  mark = (int *) calloc (ncities + 1, sizeof(int));
  cand = (int *) malloc ((size_t) (BOX_MAXT * nbK + ncities) * sizeof(int));
  for (i = bs->start[0]; i < bs->start[1]; ++i) inR[bs->loops[i]] = 1;  /* R = box 0 */

  for (j = 1; j < bs->nb; ++j) {
    lo = bs->start[j];
    s = bs->start[j + 1] - lo;

    /* running-loop cities near the new box: its cities' neighbours in R */
    stamp++;
    nc = 0;
    for (i = lo; i < lo + s; ++i)
      for (k = 0; k < nbK; ++k) {
	c = neighbors (bs->loops[i])[k];
	if (inR[c] && mark[c] != stamp) { mark[c] = stamp; cand[nc++] = c; }
      }
    if (nc == 0)                                   /* none near: use the previous box */
      for (i = bs->start[j - 1]; i < lo; ++i) cand[nc++] = bs->loops[i];

    /* every candidate R edge x every edge of the new loop, both options */
    bestC = BIG;
    for (i = 0; i < nc; ++i)
      for (ea = 0; ea < 2; ++ea) {
	x = cand[i];
	y = ea ? prv[x] : nxt[x];
	for (ib = lo; ib < lo + s; ++ib) {
	  p = bs->loops[ib];
	  q = nxt[p];                              /* each loop edge once */
	  cost = mergeCost (x, y, p, q, &e1, &e2);
	  if (cost < bestC) { bestC = cost; bx = x; by = y; b1 = e1; b2 = e2; }
	}
      }
    mergeCycles (bx, by, b1, b2);
    added += bestC;
    for (i = lo; i < lo + s; ++i) inR[bs->loops[i]] = 1;
  }

  cyclesToTour (bs->loops[0]);
  free (inR); free (mark); free (cand);
  freeCycles();
  return added;
}

/* ================================================================ M4: greedy global cheapest */

typedef struct { long long cost; int a, c, x, y, e1, e2; } Cand;

static Cand *heap;
static int heapN, heapCap;

static void heapPush (Cand v) {
  int i, p;
  if (heapN == heapCap) {
    heapCap = heapCap ? heapCap * 2 : 1024;
    heap = (Cand *) realloc (heap, heapCap * sizeof(Cand));
    if (!heap) fatal ("Out of memory in M4");
  }
  i = heapN++;
  while (i > 0 && heap[p = (i - 1) / 2].cost > v.cost) { heap[i] = heap[p]; i = p; }
  heap[i] = v;
}

static Cand heapPop (void) {
  Cand top = heap[0], last = heap[--heapN];
  int i = 0, ch;
  while ((ch = 2 * i + 1) < heapN) {
    if (ch + 1 < heapN && heap[ch + 1].cost < heap[ch].cost) ch++;
    if (heap[ch].cost >= last.cost) break;
    heap[i] = heap[ch];
    i = ch;
  }
  heap[i] = last;
  return top;
}

/* Cheapest merge using the edges at city a and the edges at city c. */
static int bestAt (int a, int c, Cand *out) {
  int ea, ec, x, y, p, q, e1, e2;
  long long cost;
  if (find (a) == find (c)) return 0;
  out->cost = BIG;
  for (ea = 0; ea < 2; ++ea)
    for (ec = 0; ec < 2; ++ec) {
      x = a; y = ea ? prv[a] : nxt[a];
      p = c; q = ec ? prv[c] : nxt[c];
      cost = mergeCost (x, y, p, q, &e1, &e2);
      if (cost < out->cost) {
	out->cost = cost; out->a = a; out->c = c;
	out->x = x; out->y = y; out->e1 = e1; out->e2 = e2;
      }
    }
  return 1;
}

static void pushCity (int a) {                     /* candidates of city a */
  int k;
  Cand cd;
  for (k = 0; k < nbK; ++k)
    if (bestAt (a, neighbors (a)[k], &cd)) heapPush (cd);
}

static int adjacent (int a, int b) { return nxt[a] == b || nxt[b] == a; }

static long long mergeM4 (const BoxSet *bs) {

  int a, j;
  long long added = 0;
  Cand cd, fresh;

  initCycles (bs);
  heap = 0; heapN = heapCap = 0;
  for (a = 1; a <= ncities; ++a) pushCity (a);

  while (loopsLeft > 1 && heapN > 0) {
    cd = heapPop();                                /* cheapest candidate      */
    if (find (cd.a) == find (cd.c)) continue;      /* already in one loop     */
    if (!adjacent (cd.x, cd.y) || !adjacent (cd.e1, cd.e2)) {
      if (bestAt (cd.a, cd.c, &fresh)) heapPush (fresh);   /* edges changed: */
      continue;                                    /*   re-price and requeue  */
    }
    mergeCycles (cd.x, cd.y, cd.e1, cd.e2);
    added += cd.cost;
    pushCity (cd.x); pushCity (cd.y);              /* new edges, new chances  */
    pushCity (cd.e1); pushCity (cd.e2);
  }

  /* fallback (rare): loops whose neighbours are all in their own loop are
     joined in curve order, box j to box j-1, best pair over their edges */
  for (j = 1; j < bs->nb && loopsLeft > 1; ++j)
    if (find (bs->loops[bs->start[j]]) != find (bs->loops[bs->start[j - 1]]))
      added += mergeBoxes (bs, j - 1, j);

  cyclesToTour (bs->loops[0]);
  free (heap); heap = 0;
  freeCycles();
  return added;
}

/* ================================================================ M5: loop MST */

typedef struct { long long d; int ba, bb; } Link;

static int linkCmp (const void *p, const void *q) {
  const Link *a = (const Link *) p, *b = (const Link *) q;
  if (a->d != b->d) return a->d < b->d ? -1 : 1;
  if (a->ba != b->ba) return a->ba - b->ba;
  return a->bb - b->bb;
}

static int *bf;                                    /* union-find over boxes */
static int bfind (int b) { while (bf[b] != b) { bf[b] = bf[bf[b]]; b = bf[b]; } return b; }

static long long mergeM5 (const BoxSet *bs) {

  int nb = bs->nb, j, i, k, a, c, ra, rb, nl = 0, ia, ib;
  int *boxOf;
  Link *links;
  long long added = 0, dmin, dd;

  /* ---- 1. candidate links between boxes: neighbour pairs, plus every
	     pair of boxes next to each other on the curve (keeps it connected) */
  boxOf = (int *) malloc ((ncities + 1) * sizeof(int));
  for (j = 0; j < nb; ++j)
    for (i = bs->start[j]; i < bs->start[j + 1]; ++i) boxOf[bs->loops[i]] = j;
  links = (Link *) malloc (((size_t) ncities * nbK + nb) * sizeof(Link));
  for (a = 1; a <= ncities; ++a)
    for (k = 0; k < nbK; ++k) {
      c = neighbors (a)[k];
      if (boxOf[a] != boxOf[c]) {
	links[nl].d = D(a, c); links[nl].ba = boxOf[a]; links[nl].bb = boxOf[c]; nl++;
      }
    }
  for (j = 1; j < nb; ++j) {                       /* closest pair of curve neighbours */
    dmin = BIG;
    for (ia = bs->start[j - 1]; ia < bs->start[j]; ++ia)
      for (ib = bs->start[j]; ib < bs->start[j + 1]; ++ib)
	if ((dd = D(bs->loops[ia], bs->loops[ib])) < dmin) dmin = dd;
    links[nl].d = dmin; links[nl].ba = j - 1; links[nl].bb = j; nl++;
  }

  /* ---- 2. Kruskal: shortest links first, keep those joining new groups;
	     merge the two loops at once with the best edge pair at those boxes */
  qsort (links, nl, sizeof(Link), linkCmp);
  bf = (int *) malloc (nb * sizeof(int));
  for (j = 0; j < nb; ++j) bf[j] = j;
  initCycles (bs);
  for (i = 0; i < nl && loopsLeft > 1; ++i) {
    ra = bfind (links[i].ba);
    rb = bfind (links[i].bb);
    if (ra == rb) continue;                        /* already connected */
    bf[ra] = rb;
    added += mergeBoxes (bs, links[i].ba, links[i].bb);
  }

  cyclesToTour (bs->loops[0]);
  freeCycles();
  free (bf); free (links); free (boxOf);
  return added;
}

/* ================================================================ entry point */

void mergeLoopsAdded (int merge, const BoxSet *bs, BoxStats *st, long long *added) {

  int i, j, *boxOf;
  long long a = 0, len = 0, g = 0;

  switch (merge) {
    case MERGE_M1: a = mergeM1 (bs); break;
    case MERGE_M2: a = mergeM2 (bs); break;
    case MERGE_M3: a = mergeM3 (bs); break;
    case MERGE_M4: a = mergeM4 (bs); break;
    case MERGE_M5: a = mergeM5 (bs); break;
    default: fatal ("unknown merge method");
  }
  tourSetPositions();
  if (added) *added = a;

  /* measured independently of the method: tour edges between boxes, and
     how much longer the tour is than the box loops added up */
  boxOf = (int *) malloc ((ncities + 1) * sizeof(int));
  for (j = 0; j < bs->nb; ++j)
    for (i = bs->start[j]; i < bs->start[j + 1]; ++i) boxOf[bs->loops[i]] = j;
  for (i = 0; i < ncities; ++i) {
    int u = tour[i], v = tour[(i + 1) % ncities];
    len += D(u, v);
    if (boxOf[u] != boxOf[v]) g += D(u, v);
  }
  st->glueG = g;
  st->joinCost = len - st->loopsSum;
  free (boxOf);
}

void mergeLoops (int merge, const BoxSet *bs, BoxStats *st) {
  mergeLoopsAdded (merge, bs, st, 0);
}
