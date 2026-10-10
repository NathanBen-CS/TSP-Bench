#ifndef NEIGHBORS_H
#define NEIGHBORS_H

/*********************************************************************

Nearest-neighbour candidate lists.

	neighbors(c)[0 .. nbK-1]   the nbK cities closest to city c,
	                           nearest first

Local search only tries moves that connect a city to one of its nbK
nearest neighbours.  Good tours almost never use longer edges, so this
loses little quality and turns an O(N^2) search into O(N * nbK).

 *********************************************************************/

extern int  nbK;          /* neighbours per city */
extern int *nbList;       /* nbList[c*nbK + k]   */

static inline int *neighbors (int c) { return nbList + (long long)c * nbK; }

void buildNeighbors (int k);

/* k nearest among a subset: for each cities[i] (i < m), its k nearest cities
   within cities[] go to out[i*k ..].  Geometric instances only; k <= m-1. */
void nearestAmong (const int *cities, int m, int k, int *out);

#endif
