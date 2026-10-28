#ifndef BOXES_H
#define BOXES_H

/*********************************************************************

Box start tour (Karp-style partitioning heuristic; no worst-case guarantee):

  cities -> sort along a space-filling curve (curves.c)
         -> boxes = cells of the curve with at most t cities each
         -> exact shortest loop inside each box (Held-Karp DP)
         -> merge the box loops into one tour (merges.c, methods M1-M5)

See boxes.c for the details.

 *********************************************************************/

#define BOX_BITS 20      /* Hilbert grid 2^20 x 2^20; also the quadtree depth cap */
#define BOX_MAXT 12      /* largest box the DP is allowed to solve exactly        */

typedef struct {
  int       leaves;      /* number of boxes                                     */
  int       capped;      /* boxes over t cities (depth cap: duplicate points)   */
  int       largest;     /* most cities in any box                              */
  long long loopsSum;    /* sum of the exact box-loop lengths                   */
  long long glueG;       /* total length of tour edges that join two boxes      */
  long long joinCost;    /* start tour length - loopsSum (net cost of merging)  */
  double    sortTime, splitTime, dpTime, glueTime;    /* seconds per stage    */
} BoxStats;

/* The boxes and their loops, before merging. */
typedef struct {
  int  nb;               /* number of boxes                                     */
  int *order;            /* cities in curve order                               */
  int *start;            /* box j = positions start[j] .. start[j+1]-1          */
  int *loops;            /* box j's exact loop, in loops[start[j] ..]           */
} BoxSet;

/* Split the sorted keys[0..m-1] into curve cells of at most t cities.
   A cell at depth L covers base^(digits-L) consecutive keys.  Box j covers
   sorted positions start[j] .. start[j+1]-1 (start[] needs m+1 entries).
   level[] (may be NULL) receives each box's depth.  Returns the box count. */
int boxSplitBase (const unsigned long long *keys, int m, int base, int digits,
		  int t, int *start, int *level);

/* Hilbert / quadtree version: base 4, `bits` digits. */
int boxSplit (const unsigned long long *keys, int m, int bits, int t, int *start);

/* Exact shortest closed loop through cities[0..s-1] (s <= BOX_MAXT).
   Writes the visiting order (city numbers) to loop[] and returns its length. */
long long heldKarpLoop (const int *cities, int s, int *loop);

/* Sort along `curve`, cut into boxes, solve every box exactly. */
void boxBuild (int curve, int t, BoxSet *bs, BoxStats *st);
void boxFree (BoxSet *bs);

/* Whole box tour with the original method (Hilbert curve, merge M1) into
   tour[] / pos[] (see tour.h). */
void boxTour (int t, BoxStats *st);

/* Any curve and merge method (merge ids in merges.h). */
void boxTourWith (int curve, int merge, int t, BoxStats *st);

#endif
