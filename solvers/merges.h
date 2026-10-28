#ifndef MERGES_H
#define MERGES_H

#include "boxes.h"

/*********************************************************************

Ways to merge the box loops (each box's exact tour) into one tour.

	M1  snake glue: cut each loop once, chain the lines in curve order,
	    each cut chosen on its own with a one-box look-ahead (original)
	M2  snake with chosen cuts: same chain, but every loop's cut and
	    direction chosen together by dynamic programming (exact for the
	    snake shape, including the closing edge)
	M3  best edge pair in curve order: merge box 1 with box 2 over all
	    edge pairs, then the running loop with box 3, and so on; only
	    running-loop edges near the new box are tried
	M4  greedy global cheapest: candidate merges from neighbour lists,
	    always apply the cheapest one anywhere, until one loop is left
	M5  loop MST: minimum spanning tree over the loops by their closest
	    cities, then best-pair merges along the tree edges

M3-M5 need the neighbour lists (buildNeighbors) to exist.

 *********************************************************************/

enum { MERGE_M1, MERGE_M2, MERGE_M3, MERGE_M4, MERGE_M5, NMERGES };

extern const char *mergeName[NMERGES];          /* "M1" .. "M5" */

int mergeByName (const char *name);             /* -1 if unknown */
int mergeNeedsNeighbors (int merge);

/* Merge the loops of bs into tour[] / pos[].  Fills st->glueG (length of
   tour edges between different boxes), st->joinCost (tour length minus
   the sum of box loops) and *added (the extra length the method itself
   computed while merging; must equal st->joinCost).                     */
void mergeLoopsAdded (int merge, const BoxSet *bs, BoxStats *st, long long *added);
void mergeLoops (int merge, const BoxSet *bs, BoxStats *st);

#endif
