#ifndef RINGS2D_H
#define RINGS2D_H

#include "boxes.h"

/*********************************************************************

Hull rings ("onion peeling", convex layers) as a box set, 2D only.

  ring 0 = convex hull of all cities; ring 1 = hull of the rest; ...

Each ring's loop is its hull boundary in counter-clockwise order, with
collinear boundary points kept in boundary order and duplicate points kept
next to each other.  For points in convex position that order is an
optimal tour of the ring, so no Held-Karp is needed.  The rings fill a
BoxSet (boxes.h) outermost first, so the merges in merges.c join them
exactly as they join box loops.

See rings2d.c for the algorithm.  No worst-case guarantee is claimed.

 *********************************************************************/

typedef struct {
  int       rings;        /* number of rings (convex layers)                 */
  int       largest;      /* most cities in any ring                         */
  int       small;        /* rings of 1 or 2 cities                          */
  long long edges;        /* ring edges: s for a ring of s >= 2, 0 for s = 1 */
  double    sortTime;     /* presort by (x, y)                               */
  double    peelTime;     /* all the monotone-chain passes                   */
} RingStats;

/* Peel every city into rings: fills bs (order = loops = ring orders,
   outermost first; start[]), and st->leaves, largest, loopsSum, the sort
   and split times (dpTime = 0, capped = 0).                              */
void ringBuild (BoxSet *bs, BoxStats *st, RingStats *rs);

/* Hull order alone (the reference, no merge optimisation): ring 0 from
   its first city round counter-clockwise, then each ring entered at its
   city nearest to the last city of the ring before, once round.  Writes
   tour[] and pos[].                                                       */
void hullOrderTour (const BoxSet *bs);

/* Debug checks; prints each failure and returns how many there were.
     - every city in exactly one ring
     - every ring convex and counter-clockwise: no right turn, and its
       turning angles add up to one full turn
     - with containment != 0 (O(N * ring size) per ring): every city of a
       later ring lies strictly inside the ring, so no boundary point was
       left for a later ring                                              */
int ringCheck (const BoxSet *bs, int containment);

#endif
