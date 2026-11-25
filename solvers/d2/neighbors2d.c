#include <math.h>
#include <stdlib.h>

#include "Instance.h"
#include "Utilities.h"
#include "neighbors.h"

/*********************************************************************

Builds the nbK-nearest-neighbour list of every city.

Geometric instances (EUC_2D, CEIL_2D) use a uniform grid: the bounding box
is cut into square cells holding ~2 cities each, and each city searches
outward ring by ring around its own cell until no unseen city can be closer
than its current k-th best.  Roughly O(N * k) for evenly spread cities.
The same grid search also works on any subset of the cities
(nearestAmong), which greedystart.c uses to join leftover fragments.

Matrix instances have no coordinates, so each city scans its whole row of
the distance matrix: O(N^2), fine for the matrix sizes in the testbed.

 *********************************************************************/

int  nbK;
int *nbList;

/*******************************************************************/

/* Offer candidate "id" at distance "d" to a sorted best-k list
   (ids[0..*count-1], dist ascending).  Keeps only the k closest. */

static void offer (int *ids, double *dist, int *count, int k, int id, double d) {

  int j = *count;

  if (j == k) {                       /* list full ...                       */
    if (d >= dist[k-1]) return;       /* ... and not better than the worst   */
    j = k - 1;                        /* overwrite the worst slot            */
  }
  else (*count)++;                    /* list not full: grow it              */

  while (j > 0 && dist[j-1] > d) {    /* insertion sort: shift bigger ones up */
    dist[j] = dist[j-1];
    ids[j]  = ids[j-1];
    j--;
  }
  dist[j] = d;                        /* drop the new candidate into place   */
  ids[j]  = id;
}

/*******************************************************************/

static void matrixNeighbors (void) {

  int c, j, count;
  double *dist = (double *) malloc (nbK * sizeof(double));

  for (c=1; c<=ncities; ++c) {               /* for every city c ...          */
    count = 0;
    for (j=1; j<=ncities; ++j)               /* ... look at every other city  */
      if (j != c)
	offer (neighbors(c), dist, &count, nbK, j, distance (c, j));
  }
  free (dist);
}

/*******************************************************************/

/* For each of the m cities in cities[], find its k nearest cities among
   cities[] itself (geometric instances only) and store them, nearest first,
   in out[i*k .. i*k+k-1] for cities[i].  Requires k <= m-1. */

void nearestAmong (const int *cities, int m, int k, int *out) {

  double minx, maxx, miny, maxy, side, cell, d, dx, dy;
  int g, gx, gy, ncells, c, i, j, t, r, cx, cy, ix, iy, count, p, q;
  int *cellOf, *cellStart, *cellMember;
  double *dist;

  if (!isEuclidean) fatal ("nearestAmong needs city coordinates");

  /* 1. bounding box of these cities */
  minx = maxx = xx[cities[0]];
  miny = maxy = yy[cities[0]];
  for (i=1; i<m; ++i) {
    c = cities[i];
    if (xx[c] < minx) minx = xx[c];
    if (xx[c] > maxx) maxx = xx[c];
    if (yy[c] < miny) miny = yy[c];
    if (yy[c] > maxy) maxy = yy[c];
  }

  /* 2. square cells, about 2 cities per cell on average */
  side = (maxx - minx > maxy - miny) ? maxx - minx : maxy - miny;
  if (side <= 0) side = 1;                       /* all cities on one spot */
  g = (int) ceil (sqrt (m / 2.0));               /* cells along the longer side */
  if (g < 1) g = 1;
  cell = side / g;                               /* width of one cell */
  gx = (int) ((maxx - minx) / cell) + 1;         /* grid columns */
  gy = (int) ((maxy - miny) / cell) + 1;         /* grid rows    */
  ncells = gx * gy;

  /* 3. bucket the cities by cell (counting sort):
        cellMember[cellStart[q] .. cellStart[q+1]-1] = indices i in cell q */
  cellOf     = (int *) malloc (m * sizeof(int));
  cellStart  = (int *) calloc (ncells + 1, sizeof(int));
  cellMember = (int *) malloc (m * sizeof(int));
  if (!cellOf || !cellStart || !cellMember) fatal ("Out of memory building grid");

  for (i=0; i<m; ++i) {
    c = cities[i];
    cx = (int) ((xx[c] - minx) / cell);  if (cx >= gx) cx = gx - 1;
    cy = (int) ((yy[c] - miny) / cell);  if (cy >= gy) cy = gy - 1;
    cellOf[i] = cy * gx + cx;                    /* row-major cell number */
    cellStart[cellOf[i] + 1]++;                  /* count cities per cell */
  }
  for (q=0; q<ncells; ++q) cellStart[q+1] += cellStart[q];   /* running totals */
  for (i=0; i<m; ++i) cellMember[cellStart[cellOf[i]]++] = i;  /* place them  */
  for (q=ncells; q>0; --q) cellStart[q] = cellStart[q-1];    /* undo the shift */
  cellStart[0] = 0;

  /* 4. for each city, search rings of cells around its own cell.
        Cities are visited cell by cell (cellMember order), so consecutive
        searches look at the same nearby cells: much friendlier to the CPU
        cache than visiting cities in number order on million-city inputs. */
  dist = (double *) malloc (k * sizeof(double));
  for (t=0; t<m; ++t) {
    i = cellMember[t];                           /* next city, in cell order */
    c = cities[i];
    cx = cellOf[i] % gx;                         /* this city's column ... */
    cy = cellOf[i] / gx;                         /* ... and row            */
    count = 0;

    for (r=0; ; ++r) {                           /* ring r = cells exactly r away */
      for (iy = cy - r; iy <= cy + r; ++iy) {
	if (iy < 0 || iy >= gy) continue;        /* off the grid */
	for (ix = cx - r; ix <= cx + r; ++ix) {
	  if (ix < 0 || ix >= gx) continue;
	  if (iy != cy - r && iy != cy + r       /* inside rows: only the two */
	      && ix != cx - r && ix != cx + r)   /* edge columns are on ring r */
	    continue;
	  q = iy * gx + ix;
	  for (j = cellStart[q]; j < cellStart[q+1]; ++j) {   /* cities in cell */
	    p = cities[cellMember[j]];
	    if (p == c) continue;
	    dx = xx[p] - xx[c];
	    dy = yy[p] - yy[c];
	    d = sqrt (dx*dx + dy*dy);
	    offer (out + (long long) i * k, dist, &count, k, p, d);
	  }
	}
      }
      /* Every city not yet seen is in ring r+1 or beyond, so at least
	 r*cell away.  Stop once our k-th best is already that close,
	 or once the rings cover the whole grid. */
      if (count == k && dist[k-1] <= r * cell) break;
      if (r > gx && r > gy) break;
    }
  }

  free (dist);
  free (cellOf);
  free (cellStart);
  free (cellMember);
}

/*******************************************************************/

void buildNeighbors (int k) {

  int *all, c;

  nbK = (k < ncities - 1) ? k : ncities - 1;     /* can't have more than N-1 */
  nbList = (int *) malloc ((long long)(ncities + 1) * nbK * sizeof(int));
  if (nbList == 0) fatal ("Out of memory allocating neighbour lists");

  if (!isEuclidean) {
    matrixNeighbors();
    return;
  }

  all = (int *) malloc (ncities * sizeof(int));  /* every city: 1..ncities */
  for (c = 1; c <= ncities; ++c) all[c-1] = c;
  nearestAmong (all, ncities, nbK, nbList + nbK); /* city c's list starts at */
  free (all);                                     /*   nbList[c*nbK]         */
}
